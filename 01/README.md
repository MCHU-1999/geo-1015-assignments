## Assignment 01: Tanaka contours

### Aim of the Assignment
The goal is to extract isolines from a triangulated TIN, structure them, and create Tanaka contours using QGIS.

The Tanaka contours method, named after cartographer Kitiro Tanaka, shades isolines for better terrain relief visualization. Illuminated lines face the light source (North-West by convention). For this, you need to:
1. Read a gridded DTM (GeoTIFF format).
2. Randomly sample cells.
3. Create Delaunay triangulation (DT) of samples.
4. Output a GeoJSON of contours and DT (PLY format).

---

## The input/output of your Python program

### Starting Code
Help/starting code is in `/hw/01/` folder in the [TUDelft GitLab repository](https://gitlab.tudelft.nl/3d/geo1015.2024/) (NetID login required).

### Required Libraries
Install:
```bash
pip install -U startinpy
pip install -U rasterio
pip install -U numpy
```

## Marking

| Criterion                              | Points |
|----------------------------------------|--------|
| Followed all rules; compiles/runs without modifications | 2.0    |
| Isolines extracted (no consistent orientation)         | 2.0    |
| Complex cases handled (e.g., horizontal triangles)     | 1.0    |
| Orientation of isolines (counter-clockwise)            | 2.0    |
| Azimuth calculated                                     | 1.0    |
| Lightness calculated                                   | 1.0    |
| Valid output GeoJSON file                              | 1.0    |
