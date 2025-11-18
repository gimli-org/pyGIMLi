#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Computing resolution properties
===============================
This example demonstrates how block and smooth TEM modelling is done and how to
access resolution matrices and kernels.
"""
# sphinx_gallery_thumbnail_number = 2

# %%%
# We import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import em
from pygimli.viewer.mpl import drawModel1D

# %%%
# Next, we define a time vector from 100µs to 10ms
# and assume a loop of 100x100m.
# We use the TDEMBlockModelling operator to generate
# synthetic data for a three-layer case and plot it.
# Note that the output is apparent resistivity.
# We could now setup a Marquard-type block inversion
# as done in the VES example.
#
t = np.logspace(-4, -2, 20)
fopkw = dict(times=t, txArea=100**2)
synth = [50, 100, 300, 30, 300]
fopSynth = em.TDEMBlockModelling(nLayers=(len(synth)+1)//2, **fopkw)
data = fopSynth(synth)
plt.loglog(data, t, "x-")
plt.grid()

# %%%
# Instead, we want to do a smooth (Occam) style
# inversion where the layer thicknesses are known.
# For this, we use the TDEMSmoothModelling operator
# and logarithmic transforms for both model and data.
#
thk = np.logspace(0.7, 1.7, 15)
print(sum(thk))
fop = em.TDEMSmoothModelling(thk=thk, **fopkw)
inv = pg.Inversion(fop=fop)
inv.dataTrans = pg.trans.TransLog()
inv.modelTrans = pg.trans.TransLog()
model = inv.run(data, relativeError=0.03, verbose=True)

# %%%
# The inverted model is plotted along with the
# synthetic model, it resembles a smoothed image of it.
#
fig, ax = plt.subplots()
drawModel1D(ax, model=synth, label="synthetic")
drawModel1D(ax, plot="semilogx", thickness=thk, values=model, label="inverted")
ax.grid(which='minor')
_ = ax.legend()

# %%%
# We are now interested in the resolution properties of the inverse problem.
# We call the function resolution matrix and obtain both the model resolution
# and the data resolution matrices. For only one, there are the functions
# `modelResolutionMatrix` and `dataResolutionMatrix`.
from pygimli.frameworks.resolution import resolutionMatrix
RM, RD = resolutionMatrix(inv, returnRD=True)
z = np.hstack([np.cumsum(thk)-thk/2, sum(thk)])

# %%%
# We display the model resolution matrix using imshow.
#
fig, ax = plt.subplots()
im = ax.imshow(RM, vmin=-0.3, vmax=0.3, cmap="RdBu_r")
ticks = np.arange(0, len(z), 2)
labels = [f"{zi:.0f}" for zi in z[ticks]]
ax.set_xticks(ticks)
ax.set_xticklabels(labels)
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
cb = plt.colorbar(im)

# %%%
# Except the last layer, the highest values occur on the
# main diagonal with maximum values at depths of about 100-200m.
# We can plot this formal resolution
#
fig, ax = plt.subplots()
drawModel1D(ax, plot="plot", thickness=thk, values=np.diag(RM), label="resolution")

# %%%
# A single resolution kernel can be obtained by extracting a column of the
# matrix. For large-scale inverse problems, the model size might be too big to
# compute the whole matrix. One can then compute single resolution kernels by
# solving an inverse problem using `modelResolutionKernel`.
#
from pygimli.frameworks.resolution import modelResolutionKernel
fig, ax = plt.subplots()
nr = 10
drawModel1D(ax, plot="plot", thickness=thk, values=RM[:, nr], label="column")
knl = modelResolutionKernel(inv, nr)
drawModel1D(ax, plot="plot", thickness=thk, values=knl, label="kernel", ls="--")
ax.legend()
ax.hlines(z[nr], *ax.get_xlim())

# %%%
# We also have a look at the data resolution matrix telling us about the
# inter-dependency and information content of the individual data.
#
im = plt.imshow(RD, vmin=-0.3, vmax=0.3, cmap="RdBu_r")
plt.colorbar(im)

# %%%
# As expected, the neighboring time gates are highly correlated.
# The most important gates are the lowest and highest one as they define the
# first and last layers and cannot be replaced by others.
#