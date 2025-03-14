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

(first-modelling)=
# First modelling

This first modeling tutorial solves the stationary heat equation in 2D. The example is
taken from the pyGIMLi paper (<https://cg17.pygimli.org>).

```{code-cell} ipython3
import pygimli as pg
import pygimli.meshtools as mt
```

Create geometry definition for the modelling domain.

```{code-cell} ipython3
world = mt.createWorld(start=[-20, 0], end=[20, -16], layers=[-2, -8],
                       worldMarker=False)
# Create a heterogeneous block
block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
                           marker=4,  boundaryMarker=10, area=0.1)
# Merge geometrical entities
geom = world + block
ax, cb = pg.show(geom, markers=True)
```

Create a mesh from based on the geometry definition.
When calling the {func}`pg.meshtools.createMesh` function, a quality parameter
can be forwarded to Triangle, which prescribes the minimum angle allowed in
the final mesh. 

```{code-cell} ipython3
mesh = mt.createMesh(geom, quality=33, area=0.2, smooth=[1, 10])
ax, _ = pg.show(mesh)
```

Call {py:func}`pygimli.solver.solveFiniteElements` to solve the heat
diffusion equation $\nabla\cdot(a\nabla T)=0$ with $T(bottom)=0$
(boundary marker 4) and $T(top)=1$ (boundary marker 8), where $a$
is the thermal diffusivity and $T$ is the temperature distribution.
We assign thermal diffusivities to the four # regions using their marker
numbers in a dictionary (a) and the fixed temperatures at the boundaries
using Dirichlet boundary conditions with the respective markers in another
dictionary (bc)

```{code-cell} ipython3
T = pg.solver.solveFiniteElements(mesh,
                                  a={1: 1.0, 2: 2.0, 3: 3.0, 4:0.1},
                                  bc={'Dirichlet': {8: 1.0, 4: 0.0}}, verbose=True)
ax, _ = pg.show(mesh, data=T, label='Temperature $T$',
                cMap="hot_r", nCols=8, contourLines=False)

ax, _ = pg.show(geom, ax=ax, fillRegion=False)
```