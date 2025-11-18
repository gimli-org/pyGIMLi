#!/usr/bin/env python
# encoding: utf-8

r"""
Meshing the Omega aka. BERT logo
================================

This is a fun example creating a logo for the BERT software. It illustrates the
possibility to hand over matplotlib path objects to the TriangleWrapper."""

import matplotlib as mpl
import matplotlib.pyplot as plt

import pygimli as pg

###############################################################################
# We start by generating a matplotlib path respresenting the :math:`\Omega`
# character.

logo_path = mpl.textpath.TextPath((0, 0), r'$\Omega$', size=5)
patch = mpl.patches.PathPatch(logo_path)

###############################################################################
# The vertices of the path are defined as mesh nodes and connected with edges.

poly = pg.Mesh(2)

for n in patch.get_verts() * 10:
    poly.createNodeWithCheck(n)

for i in range(poly.nodeCount() - 1):
    poly.createEdge(poly.node(i), poly.node(i + 1))

edge = poly.createEdge(poly.node(poly.nodeCount() - 1), poly.node(0))

###############################################################################
# We create mesh from the polygone and set the x values as the
# data for a color transition.

mesh = pg.meshtools.createMesh(poly, area=5)

###############################################################################
# Last, we create a BERT caption, visualize the mesh and fine-tune the figure.

fig, ax = plt.subplots(figsize=(4, 3))
ax.axis('off')
offset = -10
t = ax.text(mesh.xmin() + (mesh.xmax()-mesh.xmin())/2, offset, 'BERT',
            horizontalalignment='center', size=40, fontweight='bold')
pg.show(mesh, pg.x(mesh.cellCenters()), ax=ax, cMap='Spectral_r',
        logScale=False, showLater=True, showMesh=True, colorBar=False)
ax.set_ylim(offset, mesh.ymax())
