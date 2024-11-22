# -- geo1015.2024.hw01
# -- [YOUR NAME]
# -- [YOUR STUDENT NUMBER]

import argparse
import json
import sys

import rasterio
import startinpy


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
    print("Extracting the contours at:", list(myrange))

    # -- load in memory the input GeoTIFF
    try:
        # -- this gives you a Rasterio dataset
        # -- https://rasterio.readthedocs.io/en/latest/quickstart.html
        d = rasterio.open(args.inputfile)
    except Exception as e:
        print(e)
        sys.exit()

    print("name:", d.name)
    print("crs:", d.crs)
    print("size:", d.shape)
    print("no_data:", d.nodata)

    # -- create the TIN
    # -- full docs: https://startinpy.readthedocs.io
    dt = startinpy.DT()
    dt.insert_one_pt([2, 3, 4])
    dt.insert_one_pt([12, 3, 42])
    dt.insert_one_pt([6, 1, 4])

    # -- visit the triangles
    trs = dt.triangles
    one_triangle = trs[0]
    first_vertex = one_triangle[0]
    print("z-coordinate of first vertex: ", dt.points[first_vertex][2])

    # -- create a dummy GeoJSON file
    mygeojson = create_dummy_geojson()

    # -- write the triangulation to a PLY file
    dt.write_ply("mydt.ply")
    print("File 'mydt.ply' created.")
    # -- write the contours to a GeoJSON file
    with open("mycontours.geojson", "w") as file:
        file.write(json.dumps(mygeojson, indent=2))
        print("File 'mycontour.geojson' created.")


def create_dummy_geojson():
    mygeojson = {}
    mygeojson["type"] = "FeatureCollection"
    mygeojson["features"] = []
    f = {}
    f["type"] = "Feature"
    f["geometry"] = {"type": "LineString", "coordinates": [[0, 0], [10, 10]]}
    f["properties"] = {"height": 100, "azimuth": 45, "lightness": 0.75}
    mygeojson["features"].append(f)
    return mygeojson


if __name__ == "__main__":
    main()
