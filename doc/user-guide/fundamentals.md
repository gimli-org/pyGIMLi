---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  language: python
  name: python3
---

# Fundamentals

## Software design
In applied geophysics, a lot of research efforts are directed towards the integration of different physical models and combined inversion approaches to estimate multi-physical subsurface parameters. Such approaches often require coupling of different software packages and file-based data exchange. The idea of pyGIMLi is to present a very flexible framework for geophysical modelling and inversion, which allows standard and customized modelling and inversion workflows to be realized in a reproducible manner.

The software is written in **Python** on top of a **C++ core library**, which allows a combination of flexible scripting and numerical efficiency. pyGIMLi uses selected **external dependencies** for quality constrained **mesh generation** in 2D (Triangle) and 3D (Tetgen) and **visualization** in 2D (Matplotlib) and 3D (Paraview) for example. For solving linear systems we use the open-source collection **SuiteSparse** [Chen et al., 2009], which contains multi-frontal direct and iterative solvers as well as reordering algorithms.

pyGIMLi is organized in three different abstraction levels:

**Application level**

In the application level, ready-to-use method managers and frameworks are provided. Method managers (`pygimli.manager`) hold all relevant functionality related to a geophysical method. A method manager can be initialized with a data set and used to analyze and visualize this data set, create a corresponding mesh, and carry out an inversion. Various method managers are available in `pygimli.physics`. Frameworks (`pygimli.frameworks`) are generalized abstractions of standard and advanced inversions tasks such as time-lapse or joint inversion for example. Since frameworks communicate through a unified interface, they are method independent.

**Modelling level**

In the modelling level, users can set up customized forward operators that map discretized parameter distributions to a data vector. Once defined, it is straightforward to set up a corresponding inversion workflow or combine the forward operator with existing ones.

**Equation level**

The underlying equation level allows to directly access the finite element (`pygimli.solver.solveFiniteElements()`) and finite volume (`pygimli.solver.solveFiniteVolume()`) solvers to solve various partial differential equations on unstructured meshes, i.e. to approach various physical problems with possibly complex 2D and 3D geometries.


## Module/method overview

| method  | Forward | Inverse  |  Dimension  | 
| :---------------- | :------: |:--------:|:--------:|
| {py:class}`em <pygimli.physics.em>`     | **✓** | **✓** | *1D* |
| {py:class}`ert <pygimli.physics.ert>`   | **✓** | **✓** | *2D, 3D, 4D* |
| {py:class}`gravimetry <pygimli.physics.gravimetry>`    | **✓** | **✓** | *2D, 3D* |
| {py:class}`magnetics <pygimli.physics.gravimetry>`    | **✓** | **✓** | *2D, 3D* |
| {py:class}`petro <pygimli.physics.petro>`     | **✓** | **✓** | *dimensionless* |
| {py:class}`seismics <pygimli.physics.seismics>` | **✓** | **-** | "1D" |
| {py:class}`SIP <pygimli.physics.SIP>` | **✓** | **✓** | *1D, 2D* |
| {py:class}`sNMR <pygimli.physics.sNMR>` | **✓** | **✓** | *1D* |
| {py:class}`traveltime <pygimli.physics.traveltime>` | **✓** | **✓** | *2D, (3D)* |
| {py:class}`ves <pygimli.physics.ves>` | **✓** | **✓** | *1D* |
:::


## Basic pyGIMLi classes

| pyGIMLi class              | Description |
| :---------------- | :------: |
| {py:class}` Matrix <pygimli.Matrix>` | All elements are stored column-wise, i.e. all rows *A[i]* are of type `pg.Vector`. This matrix is used for storing dense data (like ERT Jacobians) and doing simple algebra.   |
| {py:class}` RVector <pygimli.RVector>` |  One dimensional array aka Vector of limited size to store data, like ERT sensitivities.   |
| {py:class}`SparseMatrix <pygimli.SparseMatrix>` |  Used for numerical approximation of partial differential equations like finite-element or finite volume. Not typically used unless efficiency is of importance. It exists also complex-valued as pg.matrix.CSparseMatrix  |
| {py:class}` BlockMatrix <pygimli.BlockMatrix>`     |   Arbitrary matrices are combined into a logical matrix. This is of importance for inversion frameworks, e.g., concatenated Jacobian matrices during joint inversions.   |
| {py:class}`DataContainer <pygimli.DataContainer>`     |   Data container storing the individual data values as well as any description how they were obtained, e.g. the geometry of source and receivers.   |
:::

## Viewer interface !
