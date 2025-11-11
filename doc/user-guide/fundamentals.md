---
file_format: mystnb
kernelspec:
  name: python3
---

# Fundamentals

## Software design

In applied geophysics, a lot of research efforts are directed towards the
integration of different physical models and combined inversion approaches to
estimate multi-physical subsurface parameters. Such approaches often require
coupling of different software packages and file-based data exchange. The idea
of pyGIMLi is to present a very flexible framework for geophysical modelling and
inversion, which allows standard and customized modelling and inversion
workflows to be realized in a reproducible manner.

The software is written in Python on top of a C++ core library, which allows a
combination of flexible scripting and numerical efficiency. pyGIMLi uses
selected external dependencies for quality constrained mesh generation in 2D
**Triangle** and 3D **Tetgen** and visualization in 2D
(**Matplotlib**) and 3D (**Paraview**) for example. For solving linear
systems we use the open-source collection **SuiteSparse**, which contains multi-frontal direct and iterative
solvers as well as reordering algorithms.

(fig-gimliblock)=
:::{figure} ../_static/pg_design.png
:align: center
:class: wrap-fig
:::

pyGIMLi is organized in three different abstraction levels:

::::{tab-set}
:::{tab-item} **Application level**

In the application level, ready-to-use method managers and frameworks are provided. Method managers (`pygimli.manager`) hold all relevant functionality related to a geophysical method. A method manager can be initialized with a data set and used to analyze and visualize this data set, create a corresponding mesh, and carry out an inversion. Various method managers are available in `pygimli.physics`. Frameworks (`pygimli.frameworks`) are generalized abstractions of standard and advanced inversions tasks such as time-lapse or joint inversion for example. Since frameworks communicate through a unified interface, they are method independent.

:::

:::{tab-item} **Modelling level**

In the modelling level, users can set up customized forward operators that map discretized parameter distributions to a data vector. Once defined, it is straightforward to set up a corresponding inversion workflow or combine the forward operator with existing ones.

:::

:::{tab-item} **Equation level**

The underlying equation level allows to directly access the finite element (`pygimli.solver.solveFiniteElements()`) and finite volume (`pygimli.solver.solveFiniteVolume()`) solvers to solve various partial differential equations on unstructured meshes, i.e. to approach various physical problems with possibly complex 2D and 3D geometries.

:::
::::

## Module/method overview

```{eval-rst}
.. autosummary::
    :nosignatures:

    pygimli.physics.em
    pygimli.physics.ert
    pygimli.physics.gravimetry
    pygimli.physics.petro
    pygimli.physics.seismics
    pygimli.physics.SIP
    pygimli.physics.sNMR
    pygimli.physics.traveltime
    pygimli.physics.ves
```

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

## Basic pyGIMLi classes

```{eval-rst}
.. autosummary::
    :nosignatures:

    pygimli.Matrix
    pygimli.Vector
    pygimli.SparseMatrix
    pygimli.BlockMatrix
    pygimli.DataContainer
```

```{note}

Please add docstrings for abovementioned classes and then use autosummary function for user guide!
```

| pyGIMLi class              | Description |
| :---------------- | :------: |
| {py:class}`Matrix <pygimli.Matrix>` | All elements are stored column-wise, i.e. all rows *A[i]* are of type `pg.Vector`. This matrix is used for storing dense data (like ERT Jacobians) and doing simple algebra.   |
| {py:class}`Vector <pygimli.Vector>` |  One dimensional array aka Vector of limited size to store data, like ERT sensitivities.   |
| {py:class}`SparseMatrix <pygimli.SparseMatrix>` |  Used for numerical approximation of partial differential equations like finite-element or finite volume. Not typically used unless efficiency is of importance. It exists also complex-valued as pg.matrix.CSparseMatrix  |
| {py:class}`BlockMatrix <pygimli.BlockMatrix>`     |   Arbitrary matrices are combined into a logical matrix. This is of importance for inversion frameworks, e.g., concatenated Jacobian matrices during joint inversions.   |
| {py:class}`DataContainer <pygimli.DataContainer>`     |   Data container storing the individual data values as well as any description how they were obtained, e.g. the geometry of source and receivers.   |

## Viewer interface

In pygimli we provide some basic post-processing routines using the matplotlib visualization framework. The main visualization call is `pygimli.viewer.show()` which is sufficient for most meshes, fields, models and streamline views. It forwards the object to be plotted to a known visualization function. `pygimli.viewer.showMesh()` is the typical call for the visualization of 2D data.

However, the underlying `show` functions only provide an input instance and are not directly responsible for plotting the data. Depending on the data type, `pg.viewer.showMesh()` then again forwards to the following instances of `pg.viewer.mpl`:

| **Data type**              | **Draw function** |
| :---------------- | :------: |
| 2D meshes | {py:class}`drawMesh <pygimli.viewer.mpl.drawModel>` |
| Value per mesh cell - *model patch* | {py:class}`drawModel <pygimli.viewer.mpl.drawModel>` |
| Value per mesh node - *scalar field* | {py:class}`drawField <pygimli.viewer.mpl.drawField>` |
| Iterable of type `[float,float]` - *vector field* | {py:class}`drawStreams <pygimli.viewer.mpl.drawStreams>` |
| {py:class}`PosVector <pg.PosVector>` - vector field | {py:class}`drawStreams <pygimli.viewer.mpl.drawStreams>` |
| {py:class}`Pos <pg.core.Pos>` - sensor position vector | {py:class}`drawSensors <pygimli.viewer.mpl.drawSensors>` |

An empty show call creates an empty ax window:

```{code-cell} ipython3
:tags: [hide-input]

import pygimli as pg
import pygimli.meshtools as mt
world = mt.createWorld(start=[-10, 0], end=[10, -10],
                       layers=[-3, -7], worldMarker=True)
mesh = mt.createMesh(world, quality=32, area=0.2, smooth=[1, 10])

res = [[1,10],[2,100],[3,999]]

cell_markers = mesh.cellMarkers()
resistivity = [res[marker-1][1] for marker in cell_markers]
nodeDta = pg.meshtools.cellDataToNodeData(mesh, resistivity)

pg.show()
```

The functions `pygimli.viewer.show()` and `pygimli.viewer.showMesh()` return the matplotlib axis object as well as the colorbar object. This allows for further customization of the plot. The following example plots the mesh with cell and boundary markers as "minimal working example". The orientation of the colorscale can either be parallel to the x- (`orientation="horizontal"`) or y-axis (`orientation="vertical"`). Further customization of axes or the colorbar can be done by accessing the returned objects directly:

```{code-cell} ipython3
:tags: [hide-input]

ax,cbar = pg.viewer.show(mesh, markers=True, showMesh=True, orientation="vertical")
ax.set_title("Mesh with cell and boundary markers", fontweight="bold")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
```

If you want to plot data on top of a mesh, you can use the `data` argument. The following example shows electrical resistivities as cell data:

```{code-cell} ipython3
:tags: [hide-input]

ax, cbar = pg.viewer.show(mesh, data=res, logScale=True, cMin=1, cMax=1000, orientation="vertical", label="Resistivity [$\Omega$m]")
ax.set_title("Electrical resistivities as cell data", fontweight="bold")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
```

If you want to plot node data, this can be done with the same argument `data`. The following example shows the same electrical resistivities as node data:

```{code-cell} ipython3
:tags: [hide-input]

ax, cbar = pg.viewer.show(mesh, data=nodeDta, logScale=True, cMin=1, cMax=1000, orientation = "vertical", label="Resistivity [$\Omega$m]")
ax.set_title("Electrical resistivities as node data", fontweight="bold")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
```

:::{admonition} Convert cell to node data
:class: info
By using the function `pg.meshtools.cellDataToNodeData()` you can convert cell data to node data. Please refer to the documentation for further details.
:::
