(sec-glossary)=

# Glossary

:::{glossary}
GCC

: The GNU Compiler Collection (GCC). Default compiler for most linux
  systems. The Windows port is {term}`MinGW` See: <http://gcc.gnu.org/>

FEniCS

: The FEniCSx computing platform. <https://fenicsproject.org/>

Gmsh

: Gmsh: a three-dimensional finite element mesh generator with built-in
  pre- and post-processing facilities
  {cite}`GeuzaineRemacle2009`
  <http://gmsh.info/>.
  See: {py:func}`pygimli.meshtools.mesh.readGmsh`

BERT

: Boundless electrical resistivity tomography <http://www.resistivity.net/>

IPython

: An improved {term}`Python` shell that integrates nicely with
  {term}`Matplotlib`. See <http://ipython.org/>.

Matplotlib

: Matplotlib {term}`Python` package displays publication quality results.
  It displays both 1D X-Y type plots and 2D contour plots for structured
  and unstructured data. It works on all common platforms and produces
  publication quality hard copies. <http://matplotlib.org>

MinGW

: MinGW, a contraction of "Minimalist GNU for Windows", is a minimalist
  development environment for native Microsoft Windows applications. See:
  <http://www.mingw.org/>

MSYS

: MSYS, a contraction of "Minimal SYStem", is a Bourne Shell command line
  interpreter system. Offered as an alternative to Microsoft's cmd.exe,
  this provides a general purpose command line environment, which is
  particularly suited to use with MinGW, for porting of many Open Source
  applications to the MS-Windows platform. See: <http://www.mingw.org/>

NumPy

: NumPy is the fundamental package for scientific computing in Python. It
  is a Python library that provides a multidimensional array object and
  various functionalities related to it. <https://docs.scipy.org/doc/numpy/>

Paraview

: Is an open-source, multi-platform data analysis and visualization
  application. See: <http://paraview.org/>

Python

: The programming language that {term}`pyGIMLi` (and your scripts) are
  written in. See: <https://www.python.org/>

Pylab

: Meta-package importing several core packages such as {term}`NumPy`,
  {term}`SciPy`, {term}`Matplotlib`, etc. into a single namespace. This is
  usually not recommended due to possible name conflicts but provides a
  quick way to get MATLAB-like functionality.

PyVista

: 3D visualization tool based on VTK: <https://www.pyvista.org>

SciPy

: Scientific Computing Tools for Python - Open-source library with many
  numerical routines but the term is often used as a synonym for the
  scientific python community, several conferences, and the *SciPy
  Stack*, i.e. a set of core packages. <https://scipy.org/about.html>

Sphinx

: The tools used to generate the {term}`GIMLi` documentation. See:
  <http://sphinx-doc.org>

STL

: Unstructured triangulated surface file format native to the "stereolithography"
  CAD software created by 3D Systems. <https://en.wikipedia.org/wiki/STL_%28file_format%29>

SuiteSparse

: SuiteSparse is a single archive that contains packages for solving
  large sparse problems using Sparse Cholesky factorization.
  <http://faculty.cse.tamu.edu/davis/suitesparse.html>

Triangle

: A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.
  {cite}`Shewchuk96b`
  <http://www.cs.cmu.edu/~quake/triangle.html>
  See: {py:func}`pygimli.meshtools.mesh.readTriangle`

Tetgen

: A Quality Tetrahedral Mesh Generator and a 3D Delaunay Triangulator.
  {cite}`Si2015`
  <http://tetgen.org/>
  See: {py:func}`pygimli.meshtools.mesh.readTetgen`
:::