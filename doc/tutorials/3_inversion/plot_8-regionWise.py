#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Region-wise regularization
==========================

In this tutorial we like to demonstrate how to control the regularization of
subsurface regions individually by using an ERT field case. The data is a 2d
profile that was measured in 2005 on the bottom of a lake. The water body is
of course influencing the fields and needs to be treated accordingly.

We first import pygimli and the modules for ERT and mesh building.
"""

# sphinx_gallery_thumbnail_number = 12
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import ert
import pygimli.meshtools as mt

# %%%
# Data and geometry
# ~~~~~~~~~~~~~~~~~
#

# %%%
# The data was measured across a shallow lake with the most electrodes being on
# the bottom of a lake. We used cables with 2m spaced takeouts.
#

data = pg.getExampleData("ert/lake.ohm")
print(data)

# %%%
# The data consists of 658 data with current and voltage using 48 electrodes.
# We first have a look at the electrode positions measured by a stick.
#

plt.plot(pg.x(data), pg.z(data), "o-")
plt.grid()

# %%%
# On both sides, two electrodes are on shore, but the others are on the bottom
# of a shallow lake with a maximum depth of 2.5m. We first compute a geometric factor using
# the analytical formula, thereby treating the electrodes as subsurface electrodes.
#

data["k"] = ert.createGeometricFactors(data, numerical=False, dim=2)
data["rhoa"] = data["u"] / data["i"] * data["k"]
ax, cb = data.show()

# %%%
# We combined Wenner-Schlumberger (top) and Wenner-beta (bottom) data. The
# lowest resistivities correspond with the water resistivity of 22.5
# :math:`\Omega`m.
#
# The contained errors are measured standard devitations and should not be used
# for inversion. Instead, we estimate new errors using 2% and 100microVolts.
#

data["err"] = ert.estimateError(data, relativeError=0.02, absoluteUError=1e-4)
print(max(data["err"]))
# pg.show(data, data["err"]*100, label="error (%)");

# %%%
# Building a mesh with the water body
# -----------------------------------
#

# We create a piece-wise linear complex (PLC) as for a case with topography
plc = mt.createParaMeshPLC(data, paraDepth=20, boundary=1)
ax, _ = pg.show(plc, markers=True)
for i, n in enumerate(plc.nodes()[:12]):
    ax.text(n.x(), n.y(), str(i))
    print(i, n.x(), n.y())

# %%%
# So node number 10 is the left one at the shore
#

for i in range(95, plc.nodeCount()):
    print(i, plc.node(i).x(), plc.node(i).y())

# %%
# and 100 the first on the other side. We connect nodes 10 and 100 by an edge
#

plc.createEdge(plc.node(10), plc.node(100), marker=-1)
plc.addRegionMarker([50, -0.1], marker=3)
ax, _ = pg.show(plc, markers=True)

# %%%
# As the lake bottom is not a surface boundary (-1) anymore, but an inside
# boundary, we set its marker >0 by iterating through all boundaries.
#

mesh = mt.createMesh(plc, quality=34.2)
for b in mesh.boundaries():
    if b.marker() == -1 and not b.outside():
        b.setMarker(2)

print(mesh)
ax, _ = pg.show(mesh, markers=True, showMesh=True)

# %%%
# Inversion with the ERT manager
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

mgr = ert.ERTManager(data, verbose=True)
mgr.setMesh(mesh)  # use this mesh for all subsequent runs
mgr.invert()
# mgr.invert(mesh=mesh) would only temporally use the mesh

# %%%
# The fit is obviously not perfect. So we have a look at data and model
# response.
#

ax = mgr.showFit()

# %%%
# Both look very similar, but let us look at the misfit function in detail.
#

mgr.showMisfit(errorWeighted=True)

# %%%
# There is still systematics in the misfit. Ideally it should be a random
# distribution of Gaussian noise.
#

kw = dict(cMin=20, cMax=300, logScale=True, cMap="Spectral_r")
ax, cb = mgr.showResult(**kw)

# %%%
# Apparently, the two regions are already decoupled from each other which
# makes sense. Let us look in detail at the water cells by extracting the
# water body.
# Note. The manager class performs a model value permutation to fit
# the parametric mesh cell. So if you want to relate model values to the input
# mesh, you need to use the unpermutated model values directly from the
# inversion framework instance: `mgr.fw.model``
#

water = mesh.createSubMesh(mesh.cells(mesh.cellMarkers() == 3))
resWater = mgr.fw.model[len(mgr.model)-water.cellCount():]
ax, cb = pg.show(water, resWater)


# %%%
# Apparently, all values are below the expected 22.5 :math:`\Omega`m
# and some are implausibly low. Therefore we should try to limit them.
# Moreover, the subsurface structures do not look very “layered”, which is
# why we make the smoothness anisotropic.
#

mgr.inv.setRegularization(zWeight=0.1)
mgr.invert()
# mgr.invert(zWeight=0.1)  # only temporarily
ax, cb = mgr.showResult(**kw)

# %%%
# Region-specific regularization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# We first want to limit the resistivity of the water between some
# plausible bounds around our measurements.
#

mgr.inv.setRegularization(3, limits=[20, 25], trans="log")
mgr.invert()
ax, cb = mgr.showResult(**kw)

# %%%
# As a result of the log-log transform, we have a homogeneous body but
# below the lake bottom values below 20, maybe due to clay content or
# maybe as compensation of limiting the water resistivity too strong. We
# could limit the subsurface, too.
#

mgr.inv.setRegularization(2, limits=[20, 2000], trans="log")
mgr.invert()
ax, cb = mgr.showResult(**kw)

# %%%
# Apparently, this makes it harder to fit the data accurately. So maybe an
# increased clay content can be responsible for resistivity below
# 20\ :math:`\Omega`m in the mud.
#

# %%%
# Model reduction
# ~~~~~~~~~~~~~~~
#

# %%%
# Another option is to treat the water body as a homogeneous body with
# only one unknown in the inversion.
#

mgr.inv.setRegularization(limits=[0, 0], trans="log")
mgr.inv.setRegularization(3, single=True)
mgr.invert()
ax, cb = mgr.showResult(**kw)

# %%%
# The last value represents the value for the lake, close to our measurement.
# This value can, however, also be set beforehand.

mgr.inv.setRegularization(3, fix=22.5)
mgr.invert()
ax, cb = mgr.showResult(**kw)

# %%%
# We see that the lake does not appear anymore as it is not a part of the
# inversion mesh ``mgr.paraDomain`` anymore.
#

# %%%
# Instead of the standard smoothness we use geostatistical regularization.
#
mgr.inv.setRegularization(2, correlationLengths=[30, 2])
mgr.invert()
ax, cb = mgr.showResult(**kw)

# %%%
# Region coupling
# ~~~~~~~~~~~~~~~
#

# %%%
# In case (does not make sense here) the two regions should be coupled to
# each other, you can set so-called inter-region constraints.
#

mgr = ert.ERTManager(data, verbose=True)
mgr.setMesh(mesh)
print("Number of regions: ", mgr.fop.regionManager().regionCount())
mgr.inv.setRegularization(cType=1, zWeight=0.2)
mgr.fop.setInterRegionCoupling(2, 3, 1.0)  # normal coupling
mgr.invert()
ax, cb = mgr.showResult(**kw)

# %%%
# The general image is of course similar, but the structures are mirrored
# around the lake bottom. Moreover the resistivity in the lake is far too
# high. Note that all of the obtained images are equivalent with respect
# to data and errors.
#

# %%%
# Take-away messages
# ------------------
#
# -  always have a look at the data fit and get hands on data errors
# -  a lot of different models are able to fit the data, particularly in
#    the full space
# -  regions can be very specifically controlled
# -  constrain or fix whenever possibly (and reliable)
# -  sometimes geostatistic constraints outperform classical smoothness
#    but sometimes not
# -  play with regularization and keep looking at data fit
#
