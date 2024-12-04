## Assignment 02: Constructing DTMs from AHN4

### Aim of the Assignment
The goal is to create three DTMs from an input AHN(4) LiDAR dataset using different methods. The tasks include:
1. Reading an AHN4 LAZ file.
2. Outputting a gridded DTM interpolated with linear in TIN interpolation of ground-classified points.
3. Outputting a gridded DTM interpolated with Laplace interpolation of ground-classified points.
4. Extracting ground points using the GFTIN algorithm, then interpolating with linear in TIN.

**The linear in TIN interpolation, Laplace interpolation, and GFTIN algorithm must be fully developed by you.**

---

## The input/output of your Python program

### Starting Code
Code is available in the `/hw/02/` folder of the [TUDelft GitLab repository](https://gitlab.tudelft.nl/3d/geo1015.2024/) (NetID login required).

### Running the Program
Run the program with:
```bash
python geo1015_hw02.py myfile.laz
```

Installing dependencies with:
```python
pip install -r requirements.txt
```

### Optional Parameters
The program accepts these optional parameters:
1. `--input_thinning`: Thinning of the input LAZ file (default=100).
2. `--resolution`: GFTIN resolution for the ground grid (default=20).
3. `--max_dist`: GFTIN max distance parameter (default=0.5).
4. `--max_angle`: GFTIN max angle parameter (default=20.0).

### Output Files
The program outputs the following in the same folder:
1. `dtm_c_linear.tiff`: A 1mX1m GeoTIFF created using AHN4 ground classification (`ground==2`) with linear in TIN interpolation.
2. `dtm_c_laplace.tiff`: A 1mX1m GeoTIFF created using AHN4 ground classification (`ground==2`) with Laplace interpolation.
3. `dtm_gftin.tiff`: A 1mX1m GeoTIFF created using GFTIN to extract ground points and linear in TIN interpolation.

---

## Libraries Allowed
No external libraries may be used except those in `requirements.txt`. Standard Python libraries are permitted (e.g., `math`, `json`, `sys`).

---

## What You Need to Do
1. **Do not modify the main file (`geo1015_hw02.py`)**; the main function will be tested independently.
2. Complete these functions in `my_code_hw02.py`:
   - `interpolate_linear(dt, x, y)`
   - `interpolate_laplace(dt, x, y)`
   - `gftin(pts, resolution, max_dist, max_angle)`
3. Adhere to the function APIs/signatures provided.

---

## Good to Know
- All code must go in `my_code_hw02.py`. Additional helper functions can be added to this file.
- Program runtime will not affect marks.
- Pre-built functions that solve the problem cannot be used.
- Do not use `startinpy.DT.interpolate()`; interpolation steps must be manually implemented.

---

## Marking

| Criterion                              | Points |
|----------------------------------------|--------|
| Followed all rules; compiles/runs without modifications | 2.0    |
| Linear in TIN interpolation implementation           | 1.0    |
| Laplace interpolation implementation                 | 3.0    |
| GFTIN implementation                                 | 4.0    |

---

## What to Submit and How to Submit It

- Submit **one file**: `my_code_hw02.py`.
- Include your name and student number at the top of the file.
- **No report necessary.**

### Submission Options
1. Upload the file to this [Dropbox page](https://www.dropbox.com/request/kjWVsSfgADMQY8m8EeoH).
2. Alternatively, create a private GitHub repository, add the instructor as a collaborator (`@hugoledoux`), and upload a simple text file with your name, student number, and the repository URL. The latest version of `my_code_hw02.py` will be used for marking.

_Last updated: 2024-11-29 07:37_
