# -- my_code_hw02.py
# -- geo1015.2024.hw02
# -- Ming-Chieh Hu
# -- 6186416


import time

import numpy as np
import startinpy
from tqdm import tqdm
import math

def invalid_error(error):
    """If error lower than 10^-5 I'd say it's neglectable.
    """
    return abs(error) > 0.00001

def cross_2d(x, y):
    return x[0] * y[1] - x[1] * y[0]

def distance_2d(pt1: np.ndarray, pt2: np.ndarray):
    vec = pt2[:2] - pt1[:2]
    return math.sqrt(sum([x**2 for x in vec]))

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

def interpolate_linear(dt, x, y):
    """Function that interpolates at location (x,y) in a DT with the linear in TIN interpolation.
    
    Inputs:
        dt: the startinpy DT
        x:  coordinate x to interpolate
        y:  coordinate y to interpolate
    
    Output:
        - the estimated value for z
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
    w0 = sarea_2d(P, pt2, pt3) / tri_area
    w1 = sarea_2d(pt1, P, pt3) / tri_area
    w2 = sarea_2d(pt1, pt2, P) / tri_area
    result_pt = w0*pt1 + w1*pt2 + w2*pt3

    # gt = dt.interpolate({"method": "TIN"}, [[x, y]])[0]
    # if invalid_error(gt - result_pt[2]):
    #     print("expect:", gt)
    #     print("result:", result_pt[2])
    #     raise Exception("abs(error) > 10^-5, you may have a wrong answer")

    return result_pt[2]


def interpolate_laplace(dt, x, y):
    """Function that interpolates at location (x,y) in a DT with the Laplace interpolation.

    Inputs:
        dt: the startinpy DT
        x:  coordinate x to interpolate
        y:  coordinate y to interpolate

    Output:
        - the estimated value for z
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


def gftin(pts, resolution, max_dist, max_angle):
    """Function that performs ground filtering using TIN refinement and returns the DT.

    Inputs:
        pts:        the Nx3 numpy array of the PC
        resolution: resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
        max_dist:   distance_2d threshold used in the ground filtering algorithm,
        max_angle:  angle threshold used in the ground filtering algorithm in degrees,

    Output:
        the startinpy DT of the ground
    """
    # -- generate 100 points randomly in the plane
    rng = np.random.default_rng(seed=42)
    pts = rng.random((100, 3))
    pts = pts * 100
    dt = startinpy.DT()
    dt.insert(pts, insertionstrategy="AsIs")
    # -- showcase for tqdm package to see progress of a loop
    for i in tqdm(range(500)):
        time.sleep(0.01)
    return dt
