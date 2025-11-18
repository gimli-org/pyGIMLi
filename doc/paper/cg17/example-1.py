#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Minimal example of using pygimli to simulate the steady heat equation.
"""
import pygimli as pg
import pygimli.meshtools as mt

# Create geometry definition for the modelling domain
world = mt.createWorld(start=[-20, 0], end=[20, -16], layers=[-2, -8],
                       worldMarker=False)
# Create a heterogeneous block
block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0], marker=4,
                           boundaryMarker=10, area=0.1)
# Merge geometrical entities
geom = world + block
pg.show(geom, boundaryMarker=True, savefig='geometry.pdf')

# Create a mesh from the geometry definition
mesh = mt.createMesh(geom, quality=33, area=0.2, smooth=[1, 10])
pg.show(mesh, savefig='mesh.pdf')

# $\diverg(a\grad T)=0$ with $T(bottom)=1$, $T(top)=0$
T = pg.solver.solve(mesh,
                    a=[[1, 1.0], [2, 2.0], [3, 3.0], [4, 0.1]],
                    bc={'Dirichlet': {8: 1.0, 4: 0.0}},
                    verbose=True)

ax, _ = pg.show(mesh, data=T, label='Temperature $T$', cMap="hot_r",
                nLevs=11, showBoundary=True, savefig='T_field.pdf')
