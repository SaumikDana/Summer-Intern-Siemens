## Overview

This Python code simulates the temperature distribution in a domain during an additive manufacturing process, where a laser beam is used to melt a material and create a microslice layer by layer. The code uses the finite element method to solve the heat equation with convection, conduction, and heat generation terms. The toolpath information is read in from an Excel file, and the temperature is computed for each segment of the toolpath.
Requirements

The code requires the following Python packages to be installed:

    dolfin
    matplotlib
    xlrd

## Usage

The toolpath information should be provided in an Excel file named data1.xlsx. The file should have the following columns:

    Column 1: Time (in milliseconds)
    Column 2: x-coordinate of the toolpath
    Column 3: y-coordinate of the toolpath
    Column 4: Layer number (not used in this code)

To run the code, execute the following command in the terminal:

python3 additive_manufacturing.py

The output will be a sequence of contour plots showing the temperature distribution and the toolpath overlay for each segment of the toolpath.
Outputs

The code generates contour plots showing the temperature distribution and the toolpath overlay for each segment of the toolpath. The plots are saved in PNG format in the current working directory.

## References

    FEniCS Project. (n.d.). FEniCS Documentation. Retrieved from https://fenicsproject.org/documentation/
    Triangulation class — Matplotlib 3.5.0 documentation. (n.d.). Retrieved from https://matplotlib.org/stable/api/tri_api.html
    xlrd documentation. (n.d.). Retrieved from https://xlrd.readthedocs.io/en/latest/
