---
file_format: mystnb
kernelspec:
  name: python3
---

(first-steps)=
# First steps

The modelling as well as the inversion part of {term}`pyGIMLi` often requires a
spatial discretization for the domain of interest, the so called
{gimliapi}`GIMLI::Mesh`.
This tutorial shows some basic aspects of handling a mesh.

First, the library needs to be imported.
To avoid name clashes with other libraries we suggest to `import pygimli` and
alias it to the simple abbreviation `pg`: CR

```{code-cell} ipython3
import pygimli as pg
```

Every part of the c++ namespace {gimliapi}`GIMLI` is bound to python and can
be used with the leading `pg.`

For instance get the current version for `pyGIMLi` with:

```{code-cell} ipython3
print(pg.__version__)
```

Now that we know the name space {gimliapi}`GIMLI`, we can create a first mesh.
A mesh is represented by a collection of nodes, cells and boundaries,
i.e., geometrical entities.

:::{note}
A regularly spaced mesh consisting of rectangles or hexahedrons is
usually called a grid. However, a grid is just a special variant of a mesh
so GIMLi treats it the same. The only difference is how they are created.
:::

GIMLi provides a collection of tools for mesh import, export and generation.
A simple grid generation is built-in but we also provide wrappers for
unstructured mesh generations, e.g., {term}`Triangle`, {term}`Tetgen` and
{term}`Gmsh`. To create a 2d grid you need to give two arrays/lists of sample points
in x and y direction, in that order, or just numbers.

```{code-cell} ipython3
grid = pg.createGrid(x=[-1.0, 0.0, 1.0, 4.0], y=[-1.0, 0.0, 1.0, 4.0])
```

The returned object `grid` is an instance of {gimliapi}`GIMLI::Mesh` and provides
various methods for modification and io-operations. General information about the
grid can be printed using the simple print() function.

```{code-cell} ipython3
print(grid)
```

Or you can access them manually using different methods:

```{code-cell} ipython3
print('Mesh: Nodes:', grid.nodeCount(),
      'Cells:', grid.cellCount(),
      'Boundaries:', grid.boundaryCount())
```

You can iterate through all cells of the general type {gimliapi}`GIMLI::Cell`
that also provides a lot of methods. Here we list the number of nodes and the
node ids per cell:

```{code-cell} ipython3
for cell in grid.cells():
    print("Cell", cell.id(), "has", cell.nodeCount(),
          "nodes. Node IDs:", [n.id() for n in cell.nodes()],
          "Position:", cell.center())

print(type(grid.cell(0)))
```

Similarly, you can iterate through all

```{code-cell} ipython3
for node in grid.nodes()[:3]:
    print("Node", node.id(), "at", node.pos())
```

or boundaries

```{code-cell} ipython3
for bnd in grid.boundaries()[:3]:
    print("Boundary", bnd.id(),
          "at center", bnd.center(),
          "directing at", bnd.norm())

print(type(grid.boundary(0)))
```

Note that the node numbers can be more easily accessed by `pg.x(grid)` etc.

To generate the input arrays `x` and `y`, you can use the
built-in {gimliapi}`GIMLI::Vector` (pre-defined with values that are type double as
`pg.Vector`), standard python lists or {term}`numpy` arrays,
which are widely compatible with {term}`pyGIMLi` vectors.

```{code-cell} ipython3
import numpy as np

grid = pg.createGrid(x=np.linspace(-1.0, 1.0, 10),
                     y=1.0 - np.logspace(np.log10(1.0), np.log10(2.0), 10))
```

This new `grid` contains

```{code-cell} ipython3
print(grid.cellCount())
```

rectangles of type {gimliapi}`GIMLI::Quadrangle` derived from the
base type {gimliapi}`GIMLI::Cell`, edges of type {gimliapi}`GIMLI::Edge`,
which are boundaries of the general type {gimliapi}`GIMLI::Boundary`, counting to

```{code-cell} ipython3
print(grid.boundaryCount())
```

The mesh can be saved and loaded in our binary mesh format `.bms`.
Or exported into `.vtk` format for 2D or 3D visualization using
{term}`Paraview`.

However, we recommend visualizing 2-dimensional content using Python scripts
that provide better exports to graphics files (e.g., png, pdf, svg).
{term}`pygimli` provides basic post-processing routines using the {term}`matplotlib` visualization package.
The main visualization call is {py:mod}`pygimli.viewer.show`, or shortly `pg.show()`
which is sufficient for most meshes, fields, models and streamline views.

```{code-cell} ipython3
pg.show(grid)
```

```{code-cell} ipython3
pg.wait()
```
