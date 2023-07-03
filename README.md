## *Overview*

This Python code simulates the temperature distribution in a domain during an additive manufacturing process, where a laser beam is used to melt a material and create a microslice layer by layer. The code uses the finite element method to solve the heat equation with convection, conduction, and heat generation terms. The toolpath information is read in from an Excel file, and the temperature is computed for each segment of the toolpath.
Requirements

The code requires the following Python packages to be installed:

- dolfin
- matplotlib
- xlrd

## *FEniCS/dolfin*

The FEniCS Project is a collection of free and open-source software components with the common goal to enable automated solution of differential equations. The components provide scientific computing tools for working with computational meshes, finite-element variational formulations of ordinary and partial differential equations, and numerical linear algebra.

The FEniCS Project is designed as an umbrella project for a collection of interoperable components. The core components include:

- UFL (Unified Form Language): A domain-specific language embedded in Python for specifying finite element discretizations of differential equations in terms of finite element variational forms.
- FIAT (Finite element Automatic Tabulator): The finite element backend of FEniCS, a Python module for generation of arbitrary order finite element basis functions on simplices.
- FFC (FEniCS Form Compiler): A compiler for finite element variational forms taking UFL code as input and generating UFC output.
- UFC (Unified Form-assembly Code): A C++ interface consisting of low-level functions for evaluating and assembling finite element variational forms.
- Instant: A Python module for inlining C and C++ code in Python.
- DOLFIN: A C++/Python library providing data structures and algorithms for finite element meshes, automated finite element assembly, and numerical linear algebra.

DOLFIN, the computational high-performance C++ backend of FEniCS, functions as the main problem-solving environment (in both C++ and Python) and user interface. Its functionality integrates the other FEniCS components and handles communication with external libraries such as PETSc, Trilinos and Eigen for numerical linear algebra, ParMETIS and SCOTCH for mesh partitioning, and MPI and OpenMP for distributed computing.

The FEniCS Project was initiated in 2003 as a research collaboration between the University of Chicago and Chalmers University of Technology. Several institutions have been involved in the development of the project. Since 2019, a refactoring of code is in progress.

## *Usage*

The toolpath information should be provided in an Excel file named data1.xlsx. The file should have the following columns:

    Column 1: Time (in milliseconds)
    Column 2: x-coordinate of the toolpath
    Column 3: y-coordinate of the toolpath
    Column 4: Layer number (not used in this code)

To run the code, execute the following command in the terminal:

    python3 benchmark.py

For the benchmark model, And

    python3 homogenized.py

For the reduced order model

The output will be a sequence of contour plots showing the temperature distribution and the toolpath overlay for each segment of the toolpath.

## *Outputs*

The code generates contour plots showing the temperature distribution and the toolpath overlay for each segment of the toolpath. The plots are saved in PNG format in the current working directory.

## *References*

- FEniCS Project. (n.d.). FEniCS Documentation. Retrieved from https://fenicsproject.org/documentation/

- Triangulation class — Matplotlib 3.5.0 documentation. (n.d.). Retrieved from https://matplotlib.org/stable/api/tri_api.html

- xlrd documentation. (n.d.). Retrieved from https://xlrd.readthedocs.io/en/latest/

