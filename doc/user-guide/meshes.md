---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"editable": true, "slideshow": {"slide_type": ""}, "raw_mimetype": ""}

# User guide - mesh section



This part of the user guide covers mesh-related topics, starting with a [general introduction](#basics-of-the-mesh-class) to the mesh class. Moreover, this section introduces general operations to [create](#mesh-creation) or [import](#mesh-import) meshes. Moreover, the general aspects of [visualization](#mesh-visualization) are covered within this section.


## Basics of the mesh class


We start off by looking at the general anatomy of a pyGIMLi mesh. It is represented by a collection of nodes, cells and boundaries, i.e., geometrical entities:

```{code-cell} ipython3
---
editable: true
jupyter:
  source_hidden: true
slideshow:
  slide_type: ''
---
import matplotlib.pyplot as plt
import pygimli as pg

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


:::{admonition} Working with meshes
:class: tip
:::{table} **General commands for accessing different entities of a mesh**
:widths: auto
:align: center

| Command              | Useage |
| :---------------- | :------: |
| mesh.cells()          |   Allows to access mesh cells   |
| mesh.cell()          | Allows to access a single cell (using the cell ID) |
| mesh.cell().center()  | Allows to point to the center of a specific mesh cell |
| mesh.cell().center().x()/.y() | Allows to extract the x / y coordinate of a cell center |
| mesh.boundaries()  |  Accesses mesh boundaries   |
| mesh.boundary() | Points to a specific boundary of the mesh (using boundary IDs)
| mesh.boundary().center() | Allows to point to the center of a specific mesh boundary |
| mesh.boundary().center().x/.y | Allows to extract x / y coordinate of a boundary
| mesh.nodes() |  Accesses mesh nodes   |
| mesh.node() | Refers to one specific node within the mesh |
| mesh.node().x()/.y() | Allows to extract x / y coordinate of a node
:::

+++

It is common practice to classify meshes into two main types - structured and unstructured.

**Structured** meshes are meshes with implicit connectivit, exhibiting a  well-known pattern in which the cells are arranged. As the cells are in a particular order, the topology of such mesh is regular. Such meshes enable easy identification of neighboring cells and points, because of their formation and structure. Often structured meshes have orthogonal quadrilateral (2D) or hexahedral (3D) elements.

:::{admonition} Structured meshes
:class: info

A regularly spaced mesh consisting of rectangles or hexahedrons is usually called a grid. However, a grid is just a special variant of a mesh so GIMLi treats it the same. The only difference is how they are created.
:::

**Unstructured** meshes, as the name suggests, are more general and can randomly form any geometry shape. Unlike structured meshes, the connectivity pattern is not fixed hence unstructured meshes do not follow a uniform pattern. However, unstructured meshes are more flexible and thus allow for more complex applications.

```{code-cell} ipython3
---
jupyter:
  source_hidden: true
---
from pygimli.meshtools import polytools as plc
import numpy as np

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

A mesh can have different regions, which are defined by **region markers** for each cell. Region markers can be used to assign properties for forward modelling as well as to control the inversion behavior.

### Add chapter on markers ###

+++

## Mesh creation

+++

### Creating a regular mesh / grid in 2D

To create a regular 2D grid, pyGIMLi offers a variety of tools that help with the task. To create a regular grid, we first of all have to create the extent of the mesh in x and z direction. For this example, we create a mesh of _20 x 10 m_ with a regular cell size of _1 x 1 m_. After defining the extent, we can simply call the pg.meshtools.createGrid() and get an overview of the number of nodes, cells and boundaries:

```{code-cell} ipython3
xmin, xmax = -10, 10.
zmin, zmax = -10., 0.

xreg = np.linspace(xmin, xmax, 20)
zreg = np.linspace(zmin, zmax, 10)
m_reg = pg.meshtools.createGrid(xreg,zreg,marker=2)
m_reg
```

To show the mesh, we can simply call the pg.show() function:

```{code-cell} ipython3
pg.show(m_reg)
```

:::{admonition} Regular grids in 2D
:class: tip
:::{table} **General commands for 2D regular grids**
:widths: 200px
:align: center

| Command              | Useage |
| :---------------- | :------: |
| mesh.cells()          |   Allows to access mesh cells   |
| mesh.cell()          | Allows to access a single cell (using the cell ID) |
:::

+++

### Creating an irregular mesh with pyGIMLi

After covering the basics of regular meshes

+++

## Mesh import

+++

## Mesh visualization

+++

## Mesh modification

+++

## Mesh export

```{code-cell} ipython3
---
editable: true
slideshow:
  slide_type: ''
---

```
