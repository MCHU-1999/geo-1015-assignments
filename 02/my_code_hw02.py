# -- my_code_hw02.py
# -- geo1015.2024.hw02
# -- [YOUR NAME]
# -- [YOUR STUDENT NUMBER]


import time

import numpy as np
import startinpy
import tqdm


def interpolate_linear(dt, x, y):
    """
    !!! TO BE COMPLETED !!!
    !!! You are free to subdivide the functionality of this function into several functions !!!

    Function that interpolates at location (x,y) in a DT with the linear in TIN interpolation.

    Inputs:
      dt: the startinpy DT
      x:  coordinate x to interpolate
      y:  coordinate y to interpolate

    Output:
      - the estimated value for z
      - np.nan if outside the convex hull (impossible to interpolate)
        (NaN: https://numpy.org/doc/stable/reference/constants.html#numpy.nan)
    """
    # -- !!! dt.interpolate() is illegal to use for this assignment
    # -- !!! you need to write your own code, and you can use the functions in startinpy
    return dt.interpolate({"method": "TIN"}, [[x, y]])
    # -- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def interpolate_laplace(dt, x, y):
    """
    !!! TO BE COMPLETED !!!
    !!! You are free to subdivide the functionality of this function into several functions !!!

    Function that interpolates at location (x,y) in a DT with the Laplace interpolation.

    Inputs:
      dt: the startinpy DT
      x:  coordinate x to interpolate
      y:  coordinate y to interpolate

    Output:
      - the estimated value for z
      - np.nan if outside the convex hull (impossible to interpolate)
        (NaN: https://numpy.org/doc/stable/reference/constants.html#numpy.nan)
    """
    # -- !!! dt.interpolate() is illegal to use for this assignment
    # -- !!! you need to write your own code, and you can use the functions in startinpy
    return dt.interpolate({"method": "Laplace"}, [[x, y]])
    # -- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def gftin(pts, resolution, max_dist, max_angle):
    """
    !!! TO BE COMPLETED !!!
    !!! the code written below is just dummy code, replace it !!!
    !!! You are free to subdivide the functionality of this function into several functions !!!

    Function that performs ground filtering using TIN refinement and returns the DT.

    Inputs:
      pts:        the Nx3 numpy array of the PC
      resolution: resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
      max_dist:   distance threshold used in the ground filtering algorithm,
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
    for i in tqdm.tqdm(range(500)):
        time.sleep(0.01)
    return dt
