#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. _ex:koenigsee:

Field data inversion ("Koenigsee")
==================================

This minimalistic example shows how to use the Refraction Manager to invert
a field data set. Here, we consider the Koenigsee data set, which represents
classical refraction seismics data set with slightly heterogeneous overburden
and some high-velocity bedrock. The data file can be found in the `pyGIMLi
example data repository
<https://github.com/gimli-org/example-data/blob/master/traveltime/koenigsee.sgt>`_.
"""
# sphinx_gallery_thumbnail_number = 4

# We import pyGIMLi and the traveltime module.

import matplotlib.pyplot as plt
import pygimli as pg
import pygimli.physics.traveltime as tt

###############################################################################
# The helper function `pg.getExampleData` downloads the data set to a temporary
# location and loads it. Printing the data reveals that there are 714 data
# points using 63 sensors (shots and geophones) with the data columns s (shot),
# g (geophone), and t (traveltime). By default, there is also a validity flag.

data = pg.getExampleData("traveltime/koenigsee.sgt", verbose=True)
print(data)

###############################################################################
# Let's have a look at the data in the form of traveltime curves.

fig, ax = plt.subplots()
lines = tt.drawFirstPicks(ax, data)

###############################################################################
# We initialize the refraction manager.

mgr = tt.TravelTimeManager(data)

# Alternatively, one can plot a matrix plot of apparent velocities which is the
# more general function also making sense for crosshole data.

ax, cbar = mgr.showData()

###############################################################################
# Finally, we call the `invert` method and plot the result.The mesh is created
# based on the sensor positions on-the-fly.

mgr.invert(secNodes=3, paraMaxCellSize=5.0,
           zWeight=0.2, vTop=500, vBottom=5000, verbose=1)

###############################################################################
# Look at the fit between measured (crosses) and modelled (lines) traveltimes.

mgr.showFit(firstPicks=True)

###############################################################################
# You can plot only the model and customize with a bunch of keywords

ax, cbar = mgr.showResult(logScale=False, cMin=500, cMax=3000, cMap="plasma_r",
                          coverage=mgr.standardizedCoverage())
rays = mgr.drawRayPaths(ax=ax, color="k", lw=0.3, alpha=0.5)

# mgr.coverage() yields the ray coverage in m and standardizedCoverage as 0/1

###############################################################################
# You can play around with the gradient starting model (`vTop` and `vBottom`
# arguments) and the regularization strength `lam` and customize the mesh.
