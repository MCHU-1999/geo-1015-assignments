# -- my_code_hw02.py
# -- geo1015.2024.hw02
# -- Ming-Chieh Hu
# -- 6186416


import numpy as np
import startinpy
from tqdm import tqdm
import math


def distance_3d(pt1: np.ndarray, pt2: np.ndarray):
    vec = pt2 - pt1
    return math.sqrt(sum([x**2 for x in vec]))

def distance_2d(pt1: np.ndarray, pt2: np.ndarray):
    vec = pt2[:2] - pt1[:2]
    return math.sqrt(sum([x**2 for x in vec]))

def cross_2d(x, y):
    return x[0] * y[1] - x[1] * y[0]

def normal_vec_2d(pt1: np.ndarray, pt2: np.ndarray):
    vec = pt2[:2] - pt1[:2]
    return np.array([-vec[1], vec[0]])

def sarea_2d(pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray):
    n1 = pt2[:2] - pt1[:2]
    n2 = pt3[:2] - pt1[:2]
    return cross_2d(n1, n2) / 2.0

def circumcenter_2d(pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray):

    def line_intersection(mid1: np.ndarray, n1: np.ndarray, mid2: np.ndarray, n2: np.ndarray):
        A = np.column_stack((n1, -n2))
        b = mid2 - mid1
        
        try:    # Solve Ax = b
            t = np.linalg.solve(A, b)[0]
            return mid1 + t * n1
        except np.linalg.LinAlgError:
            # Handles parallel lines (this shouldn't happen)
            return None

    # Get mid points and normal vectors of two sides
    mid_12 = (pt1[:2] + pt2[:2]) / 2.0
    mid_13 = (pt1[:2] + pt3[:2]) / 2.0
    n1 = normal_vec_2d(pt1, pt2)
    n2 = normal_vec_2d(pt1, pt3)

    return line_intersection(mid_12, n1, mid_13, n2)

def invalid_error(error):
    """If error lower than 10^-5 I'd say it's neglectable.
    """
    return abs(error) > 0.00001

def interpolate_linear(dt, x, y):
    """Function that interpolates at location (x,y) in a DT with the linear in TIN interpolation.
    Parameters:
        dt: the startinpy DT
        x:  coordinate x to interpolate
        y:  coordinate y to interpolate
    Returns:
        z: the estimated value for z
        - np.nan if outside the convex hull (impossible to interpolate)
        (NaN: https://numpy.org/doc/stable/reference/constants.html#numpy.nan)
    """

    if not dt.is_inside_convex_hull([x, y]):
        return np.nan

    triangle = dt.locate([x, y])
    pt1 = dt.get_point(triangle[0])
    pt2 = dt.get_point(triangle[1])
    pt3 = dt.get_point(triangle[2])
    P = np.array([x, y])

    # Duplication Check
    for pt in [pt1, pt2, pt3]:
        if pt[0] == x and pt[1] == y:
            return pt[2]
    
    # Solving the equation with 3 coefficients from sarea_2d()
    tri_area = sarea_2d(pt1, pt2, pt3)
    w1 = sarea_2d(P, pt2, pt3) / tri_area
    w2 = sarea_2d(pt1, P, pt3) / tri_area
    w3 = sarea_2d(pt1, pt2, P) / tri_area
    result_pt = w1*pt1 + w2*pt2 + w3*pt3

    # gt = dt.interpolate({"method": "TIN"}, [[x, y]])[0]
    # if invalid_error(gt - result_pt[2]):
    #     print("expect:", gt)
    #     print("result:", result_pt[2])
    #     raise Exception("abs(error) > 10^-5, you may have a wrong answer")

    return result_pt[2]


def interpolate_laplace(dt, x, y):
    """Function that interpolates at location (x,y) in a DT with the Laplace interpolation.
    Parameters:
        dt: the startinpy DT
        x:  coordinate x to interpolate
        y:  coordinate y to interpolate
    Returns:
        z: the estimated value for z
        - np.nan if outside the convex hull (impossible to interpolate)
        (NaN: https://numpy.org/doc/stable/reference/constants.html#numpy.nan)
    """

    if not dt.is_inside_convex_hull([x, y]):
        return np.nan

    triangle = dt.locate([x, y])
    pt1 = dt.get_point(triangle[0])
    pt2 = dt.get_point(triangle[1])
    pt3 = dt.get_point(triangle[2])
    P = np.array([x, y])

    for pt in [pt1, pt2, pt3]:
        # Duplication Check
        if pt[0] == x and pt[1] == y:
            return pt[2]
        
        # Tolerance Check
        dis = distance_2d(P, pt)
        if dt.snap_tolerance > dis:
            dt.snap_tolerance = dis / 2.0
            print(f"dt.snap_tolerance updated: {dis / 2.0}")
        
    p_i, inserted, _ = dt.insert_one_pt([x, y, 0])
    if not inserted:
        raise Exception("Point insertion failed.")
    
    triangles = dt.incident_triangles_to_vertex(p_i)
    voronoi_vertices, dt_vertices, weights = [], [], []

    for triangle in triangles:
        pt1 = dt.get_point(triangle[0])
        pt2 = dt.get_point(triangle[1])
        pt3 = dt.get_point(triangle[2])
        voronoi_vertices.append(circumcenter_2d(pt1, pt2, pt3))
        dt_vertices.append(pt2)
    
    for i in range(len(voronoi_vertices)):
        te = distance_2d(P, dt_vertices[i])
        ve = distance_2d(voronoi_vertices[i], voronoi_vertices[i-1])
        weights.append(ve/te)
    
    dt_vertices_z = [vertex[2] for vertex in dt_vertices]
    z = sum(np.multiply(weights, dt_vertices_z))/ sum(weights)
    dt.remove(p_i)

    # gt = dt.interpolate({"method": "Laplace"}, [[x, y]])[0]
    # if invalid_error(gt - z):
    #     print("expect:", gt)
    #     print("result:", z)
    #     raise Exception("abs(error) > 10^-5, you may have a wrong answer")

    return z




def primary_tin(pts: np.ndarray, resolution: int = 20):
    """Generate the rudimentary initial TIN needed for ground filtering.
    Returns:
        dt: the rudimentary DT (lowest points in each grid cell)
    """

    def get_filter_grid(pts: np.ndarray, resolution: int = 20):
        """Get filter grid from points and resolution.
        Returns:
            tuple: (xmin, ymin), grid
            1. starting coordinates (x_min, y_min)
            2. grid: the grid containing all the lowest points, np.array with shape(width, height, 3).
        """
        x_max = np.max(pts[:,0])
        x_min = np.min(pts[:,0])
        y_max = np.max(pts[:,1])
        y_min = np.min(pts[:,1])

        width = int((x_max - x_min) // resolution) + 1
        height = int((y_max - y_min) // resolution) + 1

        default = np.full(3, np.inf)
        grid = np.full((height, width, 3), default)
        return (x_min, y_min), grid


    (x_min, y_min), grid = get_filter_grid(pts, resolution)
    
    for pt in pts:
        w_pos = int((pt[0] - x_min) // resolution)
        h_pos = int((pt[1] - y_min) // resolution)
        
        if pt[2] < grid[h_pos][w_pos][2]:
            grid[h_pos][w_pos] = pt

    sampled_pts = grid.reshape(-1, 3)
    dt = startinpy.DT()
    dt.insert(sampled_pts)

    return dt

def angle_3pt(pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray):
    """Calculate angle of ∠123 (pt2 is the pivot).
    Returns:
        float: angle in radians (between 0 and pi).
    """
    v1 = pt1 - pt2
    v2 = pt3 - pt2
    magnitudes = distance_3d(pt1, pt2) * distance_3d(pt3, pt2)
    if magnitudes == 0:
        return 0
    
    cosine = np.dot(v1, v2) / magnitudes
    cosine = min(cosine, 1) if cosine > 1 else max(cosine, -1)
    
    return math.acos(cosine)


def distance_angle_p(pt: np.ndarray, pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray):
    """Calculate distance and the max angle from a point to a plane (formed by 3 points).
    Returns:
        tuple: (dis, ang)
        1. dis: perpendicular distance from point pt to triangle (pt1, pt2, pt3)
        2. max angle from a point (pt) to 3 points (pt1, pt2, pt3).
    """
    
    vec1 = pt2 - pt1
    vec2 = pt3 - pt1
    A = np.column_stack((vec1, vec2))
    b = pt - pt1

    try:    # Find the least square solution of Ax = b
        w1, w2 = np.linalg.lstsq(A, b)[0]
    except np.linalg.LinAlgError:
        # Handles parallel lines (this shouldn't happen)
        return None
    
    lstsq_pt = pt1 + w1 * vec1 + w2 * vec2

    dis = distance_3d(pt, lstsq_pt)
    
    ang = max(
        angle_3pt(pt, pt1, lstsq_pt),
        angle_3pt(pt, pt2, lstsq_pt),
        angle_3pt(pt, pt3, lstsq_pt)
    )

    return dis, ang
    

def gftin(pts: np.ndarray, resolution: int, max_dist: float, max_angle: float):
    """Function that performs ground filtering using TIN refinement and returns the DT.
    Parameters:
        pts:        the Nx3 numpy array of the PC
        resolution: resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
        max_dist:   distance_2d threshold used in the ground filtering algorithm,
        max_angle:  angle threshold used in the ground filtering algorithm in degrees,
    Returns:
        dt: the startinpy DT of the ground
    """
    dt = primary_tin(pts, resolution)
    pts_copy = np.copy(pts)
    loop_count = 1

    while True:

        failed = np.full(len(pts_copy), False)
        to_delete = []
        skip_count = 0
        print(f"densifying... loop {loop_count}", end="\r")

        for i, pt in enumerate(pts_copy):
            try:
                triangle = dt.locate(pt[:2])
            except Exception as e:
                # print("Point not in convex hull, skipped.")
                failed[i] = True
                skip_count += 1
                continue

            pt1 = dt.get_point(triangle[0])
            pt2 = dt.get_point(triangle[1])
            pt3 = dt.get_point(triangle[2])
            dis, ang = distance_angle_p(pt, pt1, pt2, pt3)

            if dis <= max_dist and ang <= max_angle:
                dt.insert_one_pt(pt)
                to_delete.append(i)
            else:
                failed[i] = True
        
        pts_copy = np.delete(pts_copy, to_delete, axis=0)
        loop_count += 1
        
        if failed.all():
            break
    
    print(f"Done - None of the remaining points pass the ground test.")
    print(f"{skip_count} points out of convex hull skipped.")
    return dt
