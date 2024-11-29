# -- geo1015.2024.hw01
# -- Ming-Chieh Hu
# -- 6186416

import argparse
import json
import sys
import math
import numpy as np
import numpy.ma as ma
import rasterio
import startinpy
from typing import Iterable


# ===================================================================================================
# == Class definition (isoline segment)
# ===================================================================================================

class LineSegment:
    def __init__(self, tail: np.ndarray, head: np.ndarray):
        """Construct a 2D line segment with seperated height value.
        Takes two 3D points (ndarray) to define the head and the tail (point).
        """
        self.tail = tail[:2]
        self.head = head[:2]
        self.vector = head - tail
        self.height = head[2]

    def __repr__(self):
        return f"LineSegment{self.tail, self.head, self.height}"

    def __str__(self):
        return f"LineSegment{self.tail, self.head, self.height}"
    
    def __eq__(self, value):
        if isinstance(value, LineSegment):
            return np.array_equal(self.tail, value.tail) and np.array_equal(self.head, value.head)
        elif isinstance(value, Iterable):
            return len(value) == 2 and np.array_equal(self.tail, value[0]) and np.array_equal(self.head, value[1])
        else:
            return False

    def __lt__(self, value):
        if isinstance(value, LineSegment):
            if self.__eq__(value):
                return False
            elif self.tail[0] < value.tail[0]:
                return True
            else:
                return self.head.sum() < value.head.sum()
        else:
            raise Exception("Comparing <LineSegment object> is only possible with another <LineSegment object>.")
            
    def __gt__(self, value):
        if isinstance(value, LineSegment):
            if self.__eq__(value):
                return False
            else:
                return not self.__lt__(value)
        else:
            raise Exception("Comparing <LineSegment object> is only possible with another <LineSegment object>.")
    
    def reverse(self):
        """Reverse the direction of this line segment (vector) directly.
        """
        self.tail, self.head = (self.head, self.tail)
        self.vector = self.vector * (-1)

    def twin(self):
        """Get the twin (reverse) of this line segment (vector).
        """
        return LineSegment(self.head, self.tail)
        
    def normal_vec_2d(self):
        """Produces a 2D unit normal vector (x, y), rotate counter clockwise for 90 degree.
        """
        norm_vec = np.array([(-1)*self.vector[1], self.vector[0]])

        return norm_vec
    
    def lightness(self):
        """Return the lightness value in the range [0, 100] as QGIS expects.
        """
        l = abs(((self.azimuth() - 225) % 360) - 180)
        return (l / 180) * 100
    
    def azimuth(self):
        """Return the azimuth angle (North=0; East=90; South=180; West=270).
        """
        # Notice here using (x, y) instead of (y, x) to align North=0
        angle_rad = math.atan2(self.vector[0], self.vector[1])
        angle_deg = math.degrees(angle_rad)
        azimuth = (angle_deg + 360) % 360
        
        return azimuth
    
    def feature_obj(self):
        """Create a feature dict object.
        """
        return {
            "type": "Feature",
            "geometry": {
                "type": "LineString",
                "coordinates": [self.tail.tolist(), self.head.tolist()]
            },
            "properties": {
                "height": float(self.height),
                "azimuth": self.azimuth(),
                "lightness": self.lightness()
            }
        }


# ===================================================================================================
# == Function definition
# ===================================================================================================

def normal_vec_triangle(pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray):
    """Calculate the normal vector of a 3D triangle given three points, this vector will always facing upward.
    Returns:
        numpy.ndarray: A unit normal vector to the triangle as a NumPy array.
    """
    edge1 = pt2 - pt1
    edge2 = pt3 - pt1
    
    norm_vec = np.cross(edge1, edge2)
    if norm_vec[2] < 0:
        norm_vec = norm_vec * (-1)

    return norm_vec


def random_sample_points(raster: rasterio.io.DatasetReader, thinning: float):
    """Random sample (pixels * thinning) points.
    Returns:
        numpy.ndarray: Points with coordinates [x, y, z] stacked in a NumPy array.
    """
    # Read raster data as numpy.ma.MaskedArray (can be divided into data & mask).
    dem_with_mask = raster.read(1, masked=True)
    dem = ma.getdata(dem_with_mask)
    mask = ma.getmaskarray(dem_with_mask)

    # Reshape, invert, convert to float
    mask = ~mask.reshape(-1)
    float_mask = mask.astype(float)

    # Count numbers
    valid_pixel_num = np.count_nonzero(float_mask)
    pixel_num = raster.width * raster.height
    sample_num = math.ceil(pixel_num * thinning)

    if valid_pixel_num < sample_num:
        raise Exception("Number of pixels with data less than desired sample pixels (valid < pixel number * thinning).")
    else:
        prob_mask = float_mask / valid_pixel_num

    rng = np.random.default_rng()
    sampled_indexes = rng.choice(pixel_num, size=sample_num, replace=False, p=prob_mask)

    print("Pixel count:", pixel_num)
    print("No-data count:", pixel_num - valid_pixel_num)
    print("Sample count:", sampled_indexes)

    sampled_col = [index % raster.width for index in sampled_indexes]
    sampled_row = [index // raster.width for index in sampled_indexes]

    zz = []
    for row, col in zip(sampled_row, sampled_col):
        zz.append(dem[row][col])

    xx, yy = raster.xy(sampled_row, sampled_col, None)

    return np.vstack((xx, yy, zz)).transpose()


def is_intersected(pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray, target_z: float):
    """Determine whether a target z value is in a triangle (takes 3 points as input).
    Returns:
        bool: True if intersected, else False.
    """
    if min([pt1[2], pt2[2], pt3[2]]) > target_z:
        return False
    elif max([pt1[2], pt2[2], pt3[2]]) < target_z:
        return False
    else:
        return True


def interpolate_xy_with_z(pt1: np.ndarray, pt2: np.ndarray, target_z: float):
    """Interpolate x, y coordinates with target z value in an edge of triangle (takes 2 end points as input).
    Returns:
        numpy.ndarray: Points with target z value in an edge of triangle.
    """
    if np.array_equal(pt1, pt2):
        print("Error: pt1 cannot be the same as pt2")
        return np.array([])
    if not ((pt1[2] >= target_z and pt2[2] <= target_z) or (pt1[2] <= target_z and pt2[2] >= target_z)):
        return np.array([])
    else:
        scale_ratio = (target_z - pt1[2]) / (pt2[2] - pt1[2])
        return (pt2 - pt1) * scale_ratio + pt1


def extract_isoline_from_triangle(pt1: np.ndarray, pt2: np.ndarray, pt3: np.ndarray, target_z: float):
    """Extract line segments with target z value in a triangle (takes 3 vertices as input).
    Returns:
        numpy.ndarray: LineSegment object with target z value in a triangle.
    """
    # duplicates_flag = True
    if pt1[2] == pt2[2] and pt1[2] == pt3[2] and pt1[2] == target_z:
        pt_array = np.vstack((pt1, pt2, pt3))
    elif pt1[2] == pt2[2] and pt1[2] == target_z:
        pt_array = np.vstack((pt1, pt2))
    elif pt2[2] == pt3[2] and pt2[2] == target_z:
        pt_array = np.vstack((pt2, pt3))
    elif pt1[2] == pt3[2] and pt3[2] == target_z:
        pt_array = np.vstack((pt3, pt1))
    else:
        # duplicates_flag = False
        result_12 = interpolate_xy_with_z(pt1, pt2, target_z)
        result_23 = interpolate_xy_with_z(pt2, pt3, target_z)
        result_31 = interpolate_xy_with_z(pt3, pt1, target_z)

        pt_array = np.concatenate((result_12, result_23, result_31)).reshape((-1, 3))
        pt_array = np.unique(pt_array, axis=0)

    # Map different conditions with number of points
    match len(pt_array):
        case 1:
            # drop the data with only 1 point
            return np.array([])

        case 2:
            # intersects with 2 points (common case)
            line_seg = np.array([])
            seg = LineSegment(pt_array[0], pt_array[1])
            line_seg = np.append(line_seg, seg)

            # re-orient
            normal_tri = normal_vec_triangle(pt1, pt2, pt3)
            cp = np.cross(seg.vector, normal_tri)
            if normal_tri[0] * cp[0] >= 0 and normal_tri[1] * cp[1] >= 0:
                pass
            else:
                seg.reverse()
                        
            return line_seg

        case 3:
            # intersects with 3 points (whole triangle)
            return np.array([])

        case _:
            # default pattern (error handling)
            print(f"Error: Length of 'pt_array' should be exact 1, 2 or 3, but got length of {len(pt_array)}")
            print("pt_array:", pt_array)
            sys.exit()


def extract_isolines_from_dt(dt: startinpy.DT, z_range: range):
    all_segments = {}

    for target_z in z_range:
        print("target z:", target_z, end="\r")
        line_seg_arr = np.array([])
        for triangle in dt.triangles:
            pt1 = dt.get_point(triangle[0])
            pt2 = dt.get_point(triangle[1])
            pt3 = dt.get_point(triangle[2])

            if is_intersected(pt1, pt2, pt3, target_z):
                line_seg = extract_isoline_from_triangle(pt1, pt2, pt3, target_z)
                line_seg_arr = np.append(line_seg_arr, line_seg, axis=0)
        
        all_segments[target_z] = np.unique(line_seg_arr)
    print("Extracted contours at:", list(z_range), end="\r")
    # print(all_segments)

    return all_segments


def create_geojson(line_segments: dict[any, np.ndarray], crs: str, output_file="output.geojson"):
    """Create a GeoJSON file with the given features.
    """
    # Base GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "features": [],
        "crs": {
            "type": "name",
            "properties":
            {
                "name": f"{crs}"
            }
        }
    }
    for line_seg_arr in line_segments.values():
        for line_seg in line_seg_arr:
            geojson["features"].append(line_seg.feature_obj())
    
    # print(json.dumps(geojson, indent=2))
    with open(output_file, "w") as f:
        json.dump(geojson, f, indent=2, separators=(",", ": "))


def maximum_tolerance(raster: rasterio.io.DatasetReader):
    """To find the maximum snapping tolerance needed to split all the points.
    Returns:
        float: maximum snapping tolerance
    """
    default_tolerance = 0.001
    raster_width = raster.width
    xy_width = raster.bounds.right - raster.bounds.left
    needed_tolerance = xy_width / raster_width

    while default_tolerance > needed_tolerance:
        default_tolerance *= 0.1
    
    return default_tolerance


# ===================================================================================================
# == Main function definition
# ===================================================================================================

def main():
    parser = argparse.ArgumentParser(description="My GEO1015.2024 hw01")
    parser.add_argument("inputfile", type=str, help="GeoTIFF")
    parser.add_argument(
        "thinning", type=float, help="Thinning factor (between 0 and 1)"
    )
    parser.add_argument(
        "range", type=str, help="a Python range for the contours, eg: (0, 1000, 100)"
    )

    args = parser.parse_args()

    # Validate thinning factor
    if not 0 <= args.thinning <= 1:
        parser.error("Thinning factor must be between 0 and 1")
    # Validate the range
    try:
        tmp = list(map(int, args.range.strip("() ").split(",")))
    except:
        parser.error("range invalid")
    myrange = range(tmp[0], tmp[1], tmp[2])

    # -- load in memory the input GeoTIFF
    try:
        d = rasterio.open(args.inputfile)
    except Exception as e:
        print(e)
        sys.exit()

    print("\nRaster file metedata:")
    print("name:", d.name)
    print("crs:", d.crs.to_string())
    print("size:", d.shape)
    print("bounds:", d.bounds)
    print("no_data:", d.nodata)
    CRS = d.crs.to_string()

    # Sampling points
    print("\nSampling points...")
    points = random_sample_points(d, args.thinning)
    tolerance = maximum_tolerance(d)
    d.close()

    # Building delaunay triangulation surface
    print("\nBuilding DT...")
    print("Snapping tolerance:", tolerance)
    dt = startinpy.DT()
    dt.snap_tolerance = tolerance
    dt.insert(points, insertionstrategy="BBox")
    print("Number of vertices:", dt.number_of_vertices())
    print("Number of triangles:", dt.number_of_triangles())
    
    # Extract isoline
    print("\nExtracting isoline...")
    isolines = extract_isolines_from_dt(dt, myrange)

    # Export files
    print("\nExporting files...")
    create_geojson(isolines, CRS, "mycontours.geojson")
    dt.write_ply("mydt.ply")


if __name__ == "__main__":
    main()
    