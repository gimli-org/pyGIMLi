#!/usr/bin/env python
# encoding: utf-8
"""
3D visualization with PyVista
==============================

pyGIMLi uses `PyVista <https://pyvista.org>`_ for interactive 3D
visualization of meshes and inversion results.  This tutorial collects
common tasks and recipes so that new users can quickly get productive
with 3D plots without having to hunt through scattered examples.

Topics covered:

* Creating and displaying a 3D mesh
* Mapping scalar data onto cells
* Changing styles (surface / wireframe)
* Clipping, thresholding and slicing
* Controlling the camera
* Exporting a static image to a file
"""

# sphinx_gallery_thumbnail_number = 1
import os
import tempfile
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.viewer import pv

# %%%
# Creating a simple 3D mesh
# -------------------------
# We start with a minimal 3D example: a box-shaped world with a smaller
# resistive cube inside it.  The geometry is described as a Piecewise
# Linear Complex (PLC) and meshed with tetrahedral cells.
#

world = mt.createCube(size=[20, 20, 20], pos=[0, 0, -10], marker=1)
anomaly = mt.createCube(size=[6, 6, 4], pos=[0, 0, -8], marker=2)
geom = world + anomaly
mesh = mt.createMesh(geom, quality=1.25, area=4)
print(mesh)

# %%%
# Assigning a resistivity model
# -----------------------------
# We assign synthetic resistivity values: 100 Ω·m for the background and
# 500 Ω·m for the resistive anomaly, using the cell markers.
#

rhomap = {1: 100., 2: 500.}
res = np.array([rhomap[m] for m in mesh.cellMarkers()])
mesh["res"] = res

# %%%
# Basic 3D surface plot
# ---------------------
# :py:func:`pg.show` automatically uses PyVista when given a 3D mesh.
# We pass ``hold=True`` so that we can add more elements before calling
# ``pl.show()``.
#

pl, _ = pg.show(mesh, "res", label="Resistivity (Ω·m)", cMap="Spectral_r",
                style="surface", hold=True, cMin=50, cMax=600)
pl.camera_position = "xz"
pl.camera.azimuth = -30
pl.camera.elevation = 20
pl.camera.zoom(1.1)
_ = pl.show()

# %%%
# Wireframe overlay
# -----------------
# Setting ``showMesh=True`` draws cell edges on top of the surface, which
# is helpful for inspecting mesh quality.
#

pl, _ = pg.show(mesh, "res", cMap="Spectral_r", style="surface",
                showMesh=True, hold=True)
pl.camera_position = "xz"
pl.camera.elevation = 15
_ = pl.show()

# %%%
# Thresholding
# ------------
# The ``filter`` keyword accepts any PyVista filter name together with its
# arguments as a dictionary.  Here we threshold to only show cells with a
# resistivity above 200 Ω·m, revealing the anomaly while hiding the
# background.
#

pl, _ = pg.show(mesh, "res", cMap="Spectral_r", hold=True, style="surface",
                filter={"threshold": dict(value=200, scalars="res")})
pl.camera_position = "xz"
pl.camera.azimuth = -30
pl.camera.elevation = 20
_ = pl.show()

# %%%
# Slicing through the volume
# --------------------------
# A vertical slice reveals the interior of the model.  We can add several
# elements to the same plotter for a composite view: a semi-transparent
# outer surface and an opaque cross-section through the anomaly.
#

pl, _ = pg.show(mesh, "res", cMap="Spectral_r", hold=True, style="surface",
                alpha=0.15)
pv.drawMesh(pl, mesh, label="res", cMap="Spectral_r",
            filter={"slice": dict(normal=[0, 1, 0], origin=[0, 0, -8])})
pl.camera_position = "xz"
pl.camera.azimuth = -20
pl.camera.elevation = 20
pl.camera.zoom(1.1)
_ = pl.show()

# %%%
# Clipping the mesh
# -----------------
# A clip operation cuts the mesh along a plane and exposes the interior.
# For unstructured grids, ``crinkle=True`` is applied automatically so
# that full cells are shown rather than interpolated cut faces.
#

pl, _ = pg.show(mesh, "res", cMap="Spectral_r", hold=True, style="surface",
                filter={"clip": dict(normal=[1, 1, 0], origin=[0, 0, -8])})
pl.camera_position = "xz"
pl.camera.azimuth = -30
pl.camera.elevation = 25
_ = pl.show()

# %%%
# Saving a screenshot to a file
# -----------------------------
# Use :py:meth:`pyvista.Plotter.screenshot` to write a PNG image without
# displaying an interactive window.  This is particularly useful in
# automated scripts or when running on a headless server.
#

pl, _ = pg.show(mesh, "res", cMap="Spectral_r", hold=True, style="surface",
                filter={"threshold": dict(value=200, scalars="res")})
pl.camera_position = "xz"
pl.camera.azimuth = -30
pl.camera.elevation = 20

outfile = os.path.join(tempfile.gettempdir(), "pygimli_3d_result.png")
pl.screenshot(outfile)
print("Screenshot saved to:", outfile)
_ = pl.show()

pg.wait()
