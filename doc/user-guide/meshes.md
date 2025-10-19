---
file_format: mystnb
kernelspec:
  name: python3
---

+++ {"raw_mimetype": ""}

# Meshes

This part of the user guide covers mesh-related topics, starting with a [general introduction](#basics-of-the-mesh-class) to the mesh class. Moreover, this section introduces general operations to [create](#mesh-creation) or [import](#mesh-import) meshes. Moreover, the general aspects of [visualization](#mesh-visualization) are covered within this section.


## Basics of the mesh class


We start off by looking at the general anatomy of a pyGIMLi mesh. It is represented by a collection of nodes, cells and boundaries, i.e., geometrical entities:

```{code-cell}
:tags: [hide-input]

import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.meshtools as mt
import numpy as np

from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke


def circle(x, y, text=None, radius=0.15, c="blue"):
    circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=1,
                    edgecolor='black', facecolor=(0, 0, 0, .0125),
                    path_effects=[withStroke(linewidth=5, foreground='w')])
    ax.add_artist(circle)
    ax.plot(x, y, color=c, marker=".")
    ax.text(x, y - (radius + 0.05), text, backgroundcolor="white",
            ha='center', va='top', weight='bold', color=c)

m = pg.createGrid(4,4)

# Create matplotlib figure and set size
fig, ax = plt.subplots(figsize=(8,8), dpi=90)
ax.set_title("The anatomy of a pyGIMLi mesh", fontweight="bold")

# Visualize mesh with the generic pg.show command
pg.show(m, ax=ax)

# Number all cells
for cell in m.cells():
    node = cell.allNodes()[2]
    ax.text(node.x() - 0.5, node.y() - .05, "m.cell(%d)" % cell.id(),
            fontsize=11, ha="center", va="top", color="grey")

# Mark horizontal and vertical extent
circle(m.xmin(), m.xmax(), "m.xmin(),\nm.ymax() ", c="grey")
circle(m.xmin(), m.ymin(), "m.xmin(),\nm.ymin() ", c="grey")
circle(m.xmax(), m.ymin(), "m.xmax(),\nm.ymin() ", c="grey")
circle(m.xmax(), m.ymax(), "m.xmax(),\nm.ymax() ", c="grey")

# Mark center of a cell
cid = 2
circle(m.cell(cid).center().x(), m.cell(cid).center().y(),
       "Mesh cell center")

# Mark node
circle(m.node(1).x(), m.node(6).y(), "Mesh node position", c="green")


# Mark boundary center
bid = 15
boundary_center = m.boundary(bid).center()
circle(boundary_center.x(), boundary_center.y(),
       "Mesh boundary center", c="red")

# Create boundary around cell
bid = [14,17,21,20]
bounds = []
for i in bid:
    b = m.boundaries()[i] # equivalent to mesh.boundary(17)
    bounds.append(b)
for bound in bounds:
    n1 = bound.allNodes()[0]
    n2 = bound.allNodes()[1]
    ax.plot([n1.x(), n2.x()], [n1.y(), n2.y()], "b-", lw=3, zorder=10)
ax.annotate('Mesh cell', xy=(b.center().x(), b.center().y()),
            xycoords='data', xytext=(8, 25), textcoords='offset points',
            weight='bold', color="blue",
            arrowprops=dict(arrowstyle='->',
                            connectionstyle="arc",
                            color="blue"))

ax.text(3.0, -0.5, "Made with matplotlib & pyGIMLi",
        fontsize=10, ha="right", color='.5')
fig.tight_layout(pad=5.7)
pg.wait()
```

The mesh class holds either geometric definitions (piece-wise linear complex - **PLC**) or discretizations of the subsurface. It contains <span style="color:green">**nodes**</span>, <span style="color:red">**boundaries**</span> (edges in 2D, faces in 3D) and <span style="color:blue">**cells**</span> (triangles, quadrangles in 2D, hexahedra or tetrahedra in 3D) with associated markers and arbitrary data for nodes or cells.


:::: {admonition} Working with meshes
:class: tip

:::{table} **General commands for accessing different entities of a mesh**
:widths: auto
:align: center

| Function              | Usage |
| :---------------- | :------: |
| `mesh.cells()`          |   Allows to access mesh cells   |
| `mesh.cell(id:int)`     | Allows to access a single cell (using the cell ID) |
| `mesh.cell(id:int).center()`  | Allows to point to the center of a specific mesh cell |
| `mesh.cell(id:int).center().x()/.y()` | Allows to extract the x / y coordinate of a cell center |
| `mesh.boundaries()`  |  Accesses mesh boundaries   |
| `mesh.boundary(id:int)` | Points to a specific boundary of the mesh (using boundary IDs)
| `mesh.boundary(id:int).center()` | Allows to point to the center of a specific mesh boundary |
| `mesh.boundary(id:int).center().x/.y` | Allows to extract x / y coordinate of a boundary
| `mesh.nodes()` |  Accesses mesh nodes   |
| `mesh.node(id:int)` | Refers to one specific node within the mesh |
| `mesh.node(id:int).x()/.y()` | Allows to extract x / y coordinate of a node
:::

::::


It is common practice to classify meshes into two main types - structured and unstructured.

**Structured** meshes are meshes with implicit connectivit, exhibiting a  well-known pattern in which the cells are arranged. As the cells are in a particular order, the topology of such mesh is regular. Such meshes enable easy identification of neighboring cells and points, because of their formation and structure. Often structured meshes have orthogonal quadrilateral (2D) or hexahedral (3D) elements.

:::{admonition} Structured meshes
:class: info

A regularly spaced mesh consisting of rectangles or hexahedrons is usually called a grid. However, a grid is just a special variant of a mesh so GIMLi treats it the same. The only difference is how they are created.
:::

**Unstructured** meshes, as the name suggests, are more general and can randomly form any geometry shape. Unlike structured meshes, the connectivity pattern is not fixed hence unstructured meshes do not follow a uniform pattern. However, unstructured meshes are more flexible and thus allow for more complex applications.

```{code-cell}
---
jupyter:
  source_hidden: true
---
from pygimli.meshtools import polytools as plc

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,8), dpi=90)
ax1.set_title("Structured mesh", fontweight="bold")
ax2.set_title("Unstructured mesh", fontweight="bold")

# Create geometry
world = plc.createWorld(start=[-10, 0], end=[10, -10], marker=1,
                        worldMarker=False)
c1 = plc.createCircle(pos=[0.0, -5.0], radius=3.0, area=.3)

# Structured mesh
xmin, xmax = -10, 10.
zmin, zmax = -10., 0.

xreg = np.linspace(xmin, xmax, 20)
zreg = np.linspace(zmin, zmax, 10)
m_reg = pg.meshtools.createGrid(xreg,zreg,marker=2)

pg.show(m_reg, ax=ax1)
# Unstructured mesh
m = pg.meshtools.createMesh(world, quality=34, area=1)
pg.show(m, ax=ax2)
fig.tight_layout()
```


## Mesh creation

### Creating a regular mesh / grid in 2D

To create a regular 2D grid, pyGIMLi offers a variety of tools that help with the task. To create a regular grid, we first of all have to create the extent of the mesh in x and z direction. For this example, we create a mesh of _20 x 10 m_ with a regular cell size of _1 x 1 m_. After defining the extent, we can simply call the pg.meshtools.createGrid() and get an overview of the number of nodes, cells and boundaries:

```{code-cell}
xmin, xmax = -10, 10.
zmin, zmax = -10., 0.

xreg = np.linspace(xmin, xmax, 20)
zreg = np.linspace(zmin, zmax, 10)
m_reg = pg.meshtools.createGrid(xreg,zreg,marker=2)
m_reg
```

To show the mesh, we can simply call the pg.show() function:

```{code-cell}
pg.show(m_reg)
```

:::: {admonition} Regular grids in 2D
:class: tip

:::{table}
:widths: 200px
:align: center

| Function              | Usage |
| :---------------- | :------: |
| {py:class}`createGrid <pygimli.meshtools.createWorld>`  |   Creates a regular grid mesh from the provided x- and z-coordinates   |
:::

::::


### Creating an irregular mesh with pyGIMLi

After covering the basics of regular grids, we want to dive into the world of irregular meshes.

However, we first of all have to create a **geometry** that is used as underlying susurface model for the mesh creation.

```{code-cell}
# dimensions of the world
left = -30
right = 30
depth = 25

world = mt.createWorld(start=[left, 0],
                       end=[right, -depth],
                       layers=[-5], worldMarker=False)
print(world)
pg.show(world)
```

:::{admonition} Region markers
:class: info

A mesh can have different regions, which are defined by **region markers** for each cell. Region markers can be used to assign properties for forward modelling as well as to control the inversion behavior.

For more information, see the later section on **Markers**
:::


We are using the mt.createWorld() function to create a world based on the gíven x- & z-coordinates. The following table lists all handy functions that can be utilized when creating a geometry in pyGIMLi:

:::: {admonition} PLC creation in pyGIMLi
:class: tip
:::{table} **General commands for geometry creations**
:widths: 200px
:align: center

| Function              | Usage |
| :---------------- | :------: |
| {py:class}`createWorld <pygimli.meshtools.createWorld>`        |   Creates a PLC out of a given geometry   |
| {py:class}`createCircle <pygimli.meshtools.createCircle>`        |   Creates a circular PLC   |
| {py:class}`createCube <pygimli.meshtools.createCube>`          |   Creates a cubic PLC   |
| {py:class}`createCylinder <pygimli.meshtools.createCylinder>`       |   Creates a cylindric PLC   |
| {py:class}`createLine <pygimli.meshtools.createLine>`    |   Creates a line polygon  |
| {py:class}`createPolygon <pygimli.meshtools.createPolygon>` |   Creates a polygon from a list of vertices   |
| {py:class}`createRectangle <pygimli.meshtools.createRectangle>`      |   Creates a rectangular PLC   |
:::
::::

pyGIMLi has different ways to create meshes. mt.createMesh creates a mesh using Triangle, a two-dimensional constrained Delaunay mesh generator.

The additional input parameters control the maximum triangle area and the mesh smoothness. The quality factor prescribes the minimum angles allowed in the final mesh. This can improve numerical accuracy, however, fine meshes lead to increased computational costs. Notice that we are now using showMesh which is convenient for 2D mesh visualization.

```{code-cell}
from pygimli.viewer import showMesh

mesh = mt.createMesh(world,
                     area=2.0,
                     quality=33,
                     smooth=[2, 4] # [0:no smoothing or 1:node center or 2:weighted node center, # of iter]
                    )
showMesh(mesh, markers=True, showMesh=True);
```


## Utilizing markers

A mesh can have different regions, which are defined by **region markers** for each cell. Region markers can be used to assign properties for forward modelling as well as to control the inversion behavior. Moreover, a pyGIMLi mesh holds specific markers for all boundaries (outer mesh boundaries as well as interfaces between layers and boundaries around structures, see figure above). Those boundary markers are used to define specific properties or conditions on the boundaries of a mesh. They are essential for setting up boundary conditions in simulations, such as specifying fixed values or fluxes. Properly assigning boundary markers ensures accurate representation of physical boundaries and interactions in the model.

When constructing complex geometries out of basic geometric shapes (e.g., circle, rectangle, …) we need to be careful with both the region and boundary markers. In this user guide, we explain the differences as well as the usage of both types of markers.

### Region markers

Let's take a closer look at the region markers of the two-layer case presented above. Every region marker of a specific mesh region has a position, which it is assigned to.  By using ``, we can visualize the position in the PLC:

```{code-cell}

fig, ax = plt.subplots(1,1)
pg.show(world, ax=ax)
for nr, marker in enumerate(world.regionMarkers()):
    print('Position marker number {}:'.format(nr + 1), marker.x(), marker.y(),
          marker.z())
    ax.text(marker.x()+2,marker.y()+2, f"Region marker position {nr+1}",
            fontweight="bold",fontsize=11, ha="left", va="top", color="black")
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')

```

By default, pyGIMLi assigns unique markers to every region of the PLC, starting from 1. A marker value of 0 is assigned to a region if no region marker is found - but even if this happens, we could still manually add a marker to the specific region. For now, we want to add a polygon to the already existing 2-layer PLC. This is done using the `createPolygon()` function. Within this function, we have two adjusting variables for the region markers - *marker* and *markerPosition*. The first variable specifies the marker number and the second the position of the marker. Note, that the **position** of a region marker has to lie **within its region** to work! In this case, we specify the newly defined region to have the marker 3. The marker position is set to *[-10, -10]*.

```{code-cell}
poly = mt.createPolygon([[left, -20], [right,-15],[right,-5],[left,-5]], isClosed=True, marker=3, markerPosition=[-10, -10])
geometry = world + poly
fig, ax = plt.subplots(1,1)
pg.show(geometry, ax=ax)
for nr, marker in enumerate(geometry.regionMarkers()):
    print('Position marker number {}:'.format(nr + 1), marker.x(), marker.y(),
          marker.z())
    ax.text(marker.x()+2,marker.y()+2, f"Region marker position {nr+1}",
            fontweight="bold",fontsize=11, ha="left", va="top", color="black")
    ax.scatter(marker.x(), marker.y(), s=(4 - nr) * 20, color='k')
```

Of course it is also possible to assign the same region marker to multiple parts of your PLC. This gets especially handy as soon as you start assigning physical parameters to your meshed geometry for modelling and/or inversion purposes. For more detailed information on this matter, please take a look at the [modelling section](modelling.md) of the user guide.

### Boundary markers
If we talk about markers in a PLC, one often forgets about the boundary markers. If you create a geometry in pyGIMLi, no matter whether it is a circle, rectangle or layered world, you need to deal with boundaries of different kinds. pyGIMLi automatically assigns numbers to the boundaries and, doing so, follows a specific pattern. In order to change the default boundary markers, please refer to the [modelling section](modelling.md) of the user guide.

If we assume a simple rectangular world without additional layers, we have two different options on how to enumerate the boundaries:

1) Set `createWorld(worldMarker=False)`, which leads to enumerated boundaries (*1 = left, 2=right, 3=bottom, 4=top, 5=front and 6=back*).

2) Set `createWorld(worldMarker=True)`, which leads to the top holding the boundary marker -1 and other outer boundaries -2.


```{code-cell}
start = [left,0]
end = [right, -depth]
world1 = mt.createWorld(start=start,
                       end=end, worldMarker=False)

world2 = mt.createWorld(start=start,
                       end=end, worldMarker=True)

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,12), dpi=90)
ax1.set_title("worldMarker = False", fontweight="bold")
ax2.set_title("worldMarker = True", fontweight="bold")

pg.show(world1, markers=True, ax=ax1)
pg.show(world2, markers=True, ax=ax2)
```

Let's assume that we want to add an additional feature to our PLC, such as a rectangle. Doing so, the circle is automatically assigned a boundary marker - in this case 1. However, assuming the case where we set `worldMarker=False`, this marker is already assigned to the left boundary of the PLC (which could cause problems in any modelling application later on). So we need to adapt the boundary marker of the circle, which is simply done as follows:

```{code-cell}

c1 = plc.createRectangle(start=[-10,-5], end=[10,-10], marker=2)
geom1 = world1+c1
# %%
for i in range(geom1.boundaryCount()):
    if i >= 4:
        geom1.boundary(i).setMarker(5)

pg.show(geom1, markers=True)
```

As you can see, the boundaries around the rectangle were all changed to 5, as defined in the for-loop above. `PLC.boundaryCount()` gives us the total count of boundaries in the PLC, while `PLC.boundary().setMarker()` allows us to change a specific boundary marker.

To see applications where the indices of the PLC / mesh boundary are utilized, please refer to the [modelling section](modelling.md) of the user guide.

## Mesh import

### Import options for meshes in pyGIMLi

A broad variety of functions to import and convert different mesh types into a GIMLI mesh object exists within pyGIMLi. The following functions are the most commonly used ones:

:::: {admonition} PLC creation in pyGIMLi
:class: tip
:::{table}
:widths: 200px
:align: center

| Function              | Usage |
| :---------------- | :------: |
| {py:class}`FenicsHDF5Mesh <pygimli.meshtools.readFenicsHDF5Mesh>`        |   Read FEniCS mesh from .h5 format and return a GIMLI mesh object  |
| {py:class}`Gmsh <pygimli.meshtools.readGmsh>`        |   Read GMSH ASCII file and return a GIMLI mesh object   |
| {py:class}`HDF5Mesh <pygimli.meshtools.readHDF5Mesh>`          |   Load a mesh from HDF5 file format   |
| {py:class}`Hydrus2dMesh <pygimli.meshtools.readHydrus2dMesh>`       |   Import mesh from Hydrus 2D   |
| {py:class}`Hydrus3dMesh <pygimli.meshtools.readHydrus3dMesh>`    |   Import mesh from Hydrus 3D  |
| {py:class}`MeshIO <pygimli.meshtools.readMeshIO>` |   Read generic meshIO mesh   |
| {py:class}`STL <pygimli.meshtools.readSTL>`      |   Read STL surface mesh and converts to GIMLI mesh object   |
| {py:class}`Tetgen <pygimli.meshtools.readTetgen>` |   Read and convert a mesh from the basic Tetgen output   |
| {py:class}`Triangle <pygimli.meshtools.readTriangle>`      |   Read Triangle mesh   |
:::
::::

### Example: mesh generation using Gmsh

When the scientific task requires a complex finite-element discretization (i.e. incorporation of structural information, usage of a complete electrode model (CEM), etc.), external meshing tools with visualization capabilities may be the option of choice for some users. In general, the bindings provided by pygimli allow to interface any external mesh generation software. For examples on how to use external meshing software, refer to the [examples section](../examples/1_meshing/) on the pyGIMLi website.


## Mesh modification


pyGIMLi provides a variety of operators to modify your mesh. The following table gives an overview of the most important functions:

:::{admonition} Mesh modifications
:class: tip
:::{table}
:widths: 200px
:align: center

| Function              | Usage |
| :---------------- | :------: |
| {py:class}`merge2Meshes <pygimli.meshtools.merge2Meshes>`        |   Merges two meshes   |
| {py:class}`mergeMeshes <pygimli.meshtools.mergeMeshes>`        |   Merges two or more meshes   |
| `mesh.translate()`          |   Allows to translate a mesh   |
| `mesh.scale()`       |   Scale a mesh with provided factors   |
| `mesh.rotate()`    |   Rotate a provided mesh  |
:::


### Merging meshes

In some cases, the modelling domain may require different degrees of flexibility in separate mesh regions. In the following, we demonstrate this for a two-dimensional mesh consisting of a region with regularly spaced quadrilaterals and a region with unstructured triangles.

pyGIMLi offers the possibility to merge two meshes by calling {py:class}`merge2Meshes <pygmli.meshtools.merge2Meshes()>`.

```{code-cell}
xmin, xmax = -30, 30.
zmin, zmax = -50, -25.

xreg = np.linspace(xmin, xmax, 30)
zreg = np.linspace(zmin, zmax, 25)

mesh2merge = mt.createGrid(xreg, zreg, marker=3)
mergedMesh = mt.merge2Meshes(mesh,mesh2merge)
pg.show(mergedMesh, markers=True, showMesh=True)
```

The merged meshes appear as a singular hybrid mesh now, so we can append a triangle boundary as for non-hybrid meshes:

```{code-cell}
mesh_append = mt.appendTriangleBoundary(mergedMesh, xbound=50., ybound=50., quality=31, smooth=True,
                                 marker=4, isSubSurface=True, addNodes=10)

ax, cb = pg.show(mesh_append)
```

When merging more than two meshes, the function {py:class}`mergeMeshes() <pygmli.meshtools.mergeMeshes()>` can be utilized instead of {py:class}`merge2Meshes <pygmli.meshtools.merge2Meshes()>`.


### Translating meshes

To perform simple changes of the meshs x- and z-coordinates, we can make use of the {py:class}`translate <pygmli.meshtools.mesh().translate()>` function. The following lines of code move the mesh 500 m in x and 25 m in z direction:

```{code-cell}
translated_mesh = pg.Mesh(mesh)
translated_mesh.translate([500, 25])
pg.show(translated_mesh)
```

### Scaling meshes

Apart from moving the mesh along its axes, pyGIMLi also provides a tool to scale the mesh along specified axes - in this example, we scale along the z-axis with a factor of 2:

```{code-cell}
scaled_mesh = pg.Mesh(mesh)
scaled_mesh.scale([1, 2])
pg.show(scaled_mesh)
```

### Rotating meshes
Another valuable mesh modification tool is the mesh.rotate function. By providing a rotation angle and the rotation plane, we can adjust the orientation angle of our mesh:

```{code-cell}
import numpy as np
rotated_mesh = pg.Mesh(mesh)
rotated_mesh.rotate([0, 0, np.deg2rad(-20)])
pg.show(rotated_mesh)
```

## Mesh export


Suppose we want to continue working on our GIMLi mesh object in a different meshing tool - pyGIMLi provides a variety of export functions to transfer your GIMLi mesh into a different format:

:::: {admonition} Mesh export functions
:class: tip
:::{table}
:widths: 200px
:align: center

| Function              | Usage |
| :---------------- | :------: |
| {py:class}` FenicsHDF5Mesh <pygimli.meshtools.exportFenicsHDF5Mesh>`        |   Exports GIMLi mesh in HDF5 format suitable for Fenics   |
| {py:class}`HDF5Mesh <pygimli.meshtools.exportHDF5Mesh>`        |   Writes given GIMLI::Mesh in a hdf5 format file   |
| {py:class}` PLC <pygimli.meshtools.exportPLC>`     |   	Export a piece-wise linear complex (PLC) to a .poly file (2D or 3D)   |
| {py:class}`STL <pygimli.meshtools.exportSTL>`     |   Write STL surface mesh   |
| {py:class}`STL <pygimli.meshtools.mesh().exportVTK>`     |   Save mesh alongside properties as vtk   |
:::
::::

```{code-cell}
# Save rotated mesh from above as vtk file

rotated_mesh.exportVTK('rotated_mesh.vtk')
```
