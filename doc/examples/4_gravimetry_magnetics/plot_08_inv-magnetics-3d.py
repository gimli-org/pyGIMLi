#!/usr/bin/env python
# -*- coding: utf-8 -*-
# %%
r"""
3D magnetics modelling and inversion
====================================

Based on the synthetic model of Li & Oldenburg (1996), we demonstrate 3D
inversion of magnetic data. The forward operator bases on the formula given by
Holstein et al. (2007).

In the following, we will build the model, create synthetic data, and do
inversion using a depth-weighting function as outlined in the paper.

"""
# sphinx_gallery_thumbnail_number = 2
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.viewer import pv
from pygimli.physics.gravimetry import MagneticsModelling

# %%
# Synthetic model and data generation
# -----------------------------------
# The grid is 1x1x0.5 km with a spacing of 50 m.
#

dx = 50
x = np.arange(0., 1001, dx)
y = np.arange(0., 1001, dx)
z = np.arange(-500., .1, dx)
grid = pg.createGrid(x, y, z)
print(grid)

# %%%
# We create a 3D matrix that is later filled as vector into the grid. The
# model consists of zeros and patches of 5x6 cells per depth slice that
# are shifted by one cell for subsequent cells.
#

v = np.zeros((len(z)-1, len(y)-1, len(x)-1))
for i in range(7):
    v[1+i, 11-i:16-i, 7:13] = 0.05

grid["synth"] = v.ravel()

# %%%
# We show the model making use of the pyvista package that can be called by
# :py:func:`pygimli.show`. The mesh itself is shown as a wireframe, the anomaly
# is plotted as a surface plot using a threshold filter. After the first call
# with ``hold=True``, the plotter is used to draw any subsequent plots that can
# also be slices or clips. Moreover, the camera position is set so that the
# vertical axis is going downwards (x is Northing and y is Easting as common in
# magnetics).
#

pl, _ = pg.show(grid, style="wireframe", hold=True)
pv.drawMesh(pl, grid, label="synth", style="surface", cMap="Spectral_r",
            filter={"threshold": dict(value=0.05, scalars="synth")})
pl.camera_position = "yz"
pl.camera.roll = 90
pl.camera.azimuth = 180 - 15
pl.camera.elevation = 10
pl.camera.zoom(1.2)
_ = pl.show()

# %%%
# For the computation of the total field, we define the global magnetic
# field using the IGRF (total field, inclination and declination) settings
# given in the paper. Any global field can also be retrieved by the
# ``pyIGRF`` module.
#

F, I, D = 50000, 75, 25  # total field in nT
H = F * np.cos(np.deg2rad(I))
X = H * np.cos(np.deg2rad(D))
Y = H * np.sin(np.deg2rad(D))
Z = F * np.sin(np.deg2rad(I))
igrf = [D, I, H, X, Y, Z, F]

# Alternatively one could use pyIGRF at a specific position
# import pyIGRF
# igrf = pyIGRF.igrf_value(lat=50.59465, lon=12.64139)

py, px = np.meshgrid(x, y)
px = px.ravel()
py = py.ravel()
points = np.column_stack((px, py, np.ones_like(px)*20))

# The forward operator
cmp = ["TFA"]  # ["Bx", "By", "Bz"]
fop = MagneticsModelling(mesh=grid, points=points, cmp=cmp, igrf=igrf)
model = pg.Vector(grid.cellCount(), 1.0)
data = fop.response(grid["synth"])

# %%%
# Just like in the paper, the data are contaminated with an error model
# consisting of relative and absolute noise of 2% and 1 nT, respectively.
#

absError = np.abs(data) * 0.02 + 1
data += np.random.randn(len(data)) * absError

# %%%
# Depth weighting
# ---------------
#
# In the paper of Li & Oldenburg (1996), they propose a depth weighting of the
# constraints with the formula
#
# .. math::
#
#
#    w_z = \frac{1}{(z+z_0)^{3/2}}
#

# depth weighting
bz = np.array([abs(b.center().z()) for b in grid.boundaries() if not b.outside()])
z0 = 25
wz = 100 / (bz+z0)**1.5
print(min(wz), max(wz))

# %%%
# Inversion
# ---------
#
# The inversion is rather straightforward using the standard inversion
# framework :py:class:`pygimli.Inversion`.
#

inv = pg.Inversion(fop=fop, verbose=True)  # , debug=True)
inv.setRegularization(limits=[0, 0.07])  # to limit values
inv.setConstraintWeights(wz)
invmodel = inv.run(data, absoluteError=absError, lam=200, startModel=1e-3, verbose=True)
grid["inv"] = invmodel

# %%%
# Visualization
# -------------
#
# For showing the model, we again use the threshold filter with a value of
# 0.02. For comparison with the synthetic model, we plot the latter as a
# wireframe.
#

pl, _ = pg.show(grid, label="synth", style="wireframe", hold=True,
                filter={"threshold": dict(value=0.025, scalars="synth")})
pv.drawMesh(pl, grid, label="inv", style="surface", cMap="Spectral_r",
            filter={"threshold": dict(value=0.02, scalars="inv")})
pv.drawMesh(pl, grid, label="inv", style="surface", cMap="Spectral_r",
            filter={"slice": dict(normal=[-1, 0, 0], origin=[500, 600, 250])})
pl.camera_position = "yz"
pl.camera.roll = 90
pl.camera.azimuth = 180 - 15
pl.camera.elevation = 10
pl.camera.zoom(1.2)
_ = pl.show()

# %%%
# The model can nicely outline the top part of the anomalous body, but not
# its depth extent.
#
# We compare the data and model response by means of scatter plots:
#

fig, ax = plt.subplots(ncols=2, figsize=(12, 5), sharex=True, sharey=True)
vals = data
mm = np.max(np.abs(vals))
ax[0].scatter(px, py, c=vals, cmap="bwr", vmin=-mm, vmax=mm)
im = ax[1].scatter(px, py, c=inv.response, cmap="bwr", vmin=-mm, vmax=mm)
cb = plt.colorbar(im, ax=ax[1])
cb.set_label("B (nT)")

# %%%
# Alternatively, we can also plot the error-weighted misfit.
#

misfit = (inv.response - data) / absError
im = plt.scatter(py, px, c=misfit, cmap="bwr", vmin=-3, vmax=3)
cb = plt.colorbar(im)
cb.set_label("misfit / error")

# %%%
# References
# ----------
#
# -  Li, Y. & Oldenburg, D. (1996): 3-D inversion of magnetic data.
#    Geophysics 61(2), 394-408.
# -  Holstein, H., Sherratt, E.M., Reid, A.B.  (2007): Gravimagnetic field
#    tensor gradiometry formulas for uniform polyhedra, SEG Ext. Abstr.
#
