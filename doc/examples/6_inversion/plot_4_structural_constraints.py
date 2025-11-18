#!/usr/bin/env python
# encoding: utf-8
"""
Inversion with structural constraints
=====================================

Any inversion suffers from ambiguity and missing resolution at depth.
Structural information can provide valuable information about geological
boundaries in the subsurface. Such information can come from wave methods
like seismics (Tanner et al. 2020), ground-penetrating radar (Doetsch et
al. 2012, Jiang et al. 2020) or even piece-wise from boreholes (Wunderlich
et al. 2018). This can be done in both 2D (Jiang et al. 2020) using lines
or in 3D using facets (Doetsch et al. 2012).

We demonstrate on a very simple example from a bedrock detection how
a structural interface from a seismic refraction can lead to models that
are far easier to interpret.
"""

# We import the numpy and matplotlib
import numpy as np
# Next we import pyGIMLi and the ERT and mesh tools modules
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert

# %%%
# First we load the data and the interface
#
# sphinx_gallery_thumbnail_number = 4

data = pg.getExampleData("ert/struct.dat")
xz = pg.getExampleData("ert/struct.txt")
ert.show(data)

# %%%
# We estimate an error and do a normal inversion using the ERT manager.
# As we expect more or less layered structures, we choose a vertical
# weighting factor for the regularization of 0.5.
#

data.estimateError()
mgr = ert.Manager(data)
mgr.invert(paraDepth=150, zWeight=0.5, verbose=True)

# %%%
# We plot the inversion result using a defined color scale.
# Additionally, we plot the measured line on top of it and see that generally
# there is a change from conductive overburden to resistive bedrock.
#

ax, cb = mgr.showResult(cMin=50, cMax=1500)
ax.plot(xz[:, 0], xz[:, 1], "k-")

# %%%
# We now want to include this line into the inversion mesh.
# Therefore we generate the geometry of the inversion, a polygone from the line
# and merge both objects before we send it to the mesh generation.
#

plc = mt.createParaMeshPLC(data, paraDepth=150, boundary=1)
line = mt.createPolygon(xz, marker=1)
plc += line
mesh = mt.createMesh(plc, quality=34.3)
pg.show(mesh, markers=True, showMesh=True)

# %%%
# We see the outer boundaries have markers of -1 (homogeneous Neumann, no-flow)
# and -2 (mixed) important to define the boundary conditions of the forward
# problem. The interface between the inversion domain (yellow) and the outer
# region (blue) and the added interface are >0 and do not affect the forward.
# We now use the generated mesh in a new ERT manager by setting it.
#

mgr = ert.Manager(data)
mgr.setMesh(mesh)
mgr.invert(verbose=True)
ax, cb = mgr.showResult(cMin=50, cMax=1500)

# %%%
# As a result, we see sharp resistivity jumps across the line in the center
# of the profile, i.e. the resistivity structure follows the given hints.
# At other positions like at the beginning of the profile, it does not show
# such strong contrasts indicating that the bedrock interface is probably more
# shallow and there is likely a weathered zone.
#

# %%%
# References
# Doetsch, J., Linde, N., Pessognelli, M., Green, A.G. & Günther, T. (2012):
# Constraining 3-D electrical resistance tomography with GPR data for improved
# aquifer characterization. Journal of Applied Geophysics 78, 68-76,
# doi:10.1016/j.jappgeo.2011.04.008.
# Wunderlich, T., Fischer, P., Wilken, D., Hadler, H., Erkul, E., Mecking, R.,
# Günther, T., Heinzelmann, M., Vött, A. & Rabbel, W. (2018): Constraining
# Electric Resistivity Tomography by Direct Push Electric Conductivity
# logs and vibracores: An exemplary study of the Fiume Morto silted riverbed
# (Ostia Antica, Western Italy). Geophysics 83(3), B87-B103, doi:10.1190/geo2016-0660.1.
# Tanner, D.C., Buness, H., Igel, J., Günther, T., Gabriel, G., Skiba, P.,
# Plenefisch, T., Gestermann, N. & Walter, T. (2020): Fault Detection. in:
# Tanner, C.D. & Brandes, C. (Eds.): Understanding Faults, 380p., Elsevier,
# 81-146, doi:10.1016/B978-0-12-815985-9.00003-5.
# Jiang, C., Igel, J., Dlugosch, R., Müller-Petke, M., Günther, T., Helms, J.,
# Lang, J. & Winsemann (2020): Magnetic resonance tomography constrained by
# ground-penetrating radar for improved hydrogeophysical characterisation,
# Geophysics 85(6), JM13-JM26, doi:10.1190/geo2020-0052.1.
# Norooz, R., Olsson, P.-I., Dahlin, T., Günther, T. & Bernstone, C. (2021): A
# geoelectrical pre-study of Älvkarleby test embankment dam: 3D forward
# modelling and effects of structural constraints on the 3D inversion model of
# zoned embankment dams. J. Appl. Geophys. 191, 104355, doi:10.1016/j.jappgeo.2021.104355.
#