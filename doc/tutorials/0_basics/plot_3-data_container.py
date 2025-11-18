#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
The DataContainer class
=======================
Data are often organized in a data container storing the individual data values
as well as any description how they were obtained, e.g. the geometry of source
and receivers.

So first a data container holds vectors like in a dictionary, however, all of
them need to have the same length defined by the .size() method.
Assume we want to store Vertical Electrical Sounding (VES) data.
"""
# sphinx_gallery_thumbnail_number = 2

# We start off with the typical imports
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import VESManager

# %%%
# We define logarithmically equidistant AB/2 spacings
#

ab2 = np.logspace(0, 3, 11)
print(ab2)

# %%%
# We create an empty data container
ves = pg.DataContainer()
print(ves)

# %%%
# We feed it into the data container just like in a dictionary.
#

ves["ab2"] = ab2
ves["mn2"] = ab2 / 3
print(ves)

# %%%
# We now want to do a VES simulation and use the VES Manager for this task.
#

mgr = VESManager()
model = [10, 10, 100, 10, 1000]
ves["rhoa"] = mgr.simulate(model, ab2=ves["ab2"], mn2=ves["mn2"])
print(ves)

# %%%
# We can plot the sounding curve by assessing its fields
#

fig, ax = plt.subplots()
ax.loglog(ves["rhoa"], ves["ab2"], "x-");
ax.set_ylim(ax.get_ylim()[::-1])
ax.grid(True)
ax.set_xlabel("Apparent resistivity (Ohmm)")
ax.set_ylabel("AB/2 (m)");

# %%%
# A data container can be saved to disk
#

ves.save("ves.data")
print(open("ves.data").read())

# %%%
# The data are (along with a valid flat) in the second section.
# We can add arbitrary entries to the data container but define what to save.
#

ves["flag"] = pg.Vector(ves["rhoa"] > 100) + 1
print(ves)
ves.save("ves.data", "ab2 mn2 rhoa")
print(open("ves.data").read())

# %%%
# We can mask or unmask the data with a boolean vector.

ves.markValid(ves["ab2"] > 2)
ves.save("ves.data", "ab2 rhoa")  # note that only valid data are saved!
print(ves)

# %%%
# Data containers with indexed data
# ---------------------------------
#
# Assume we have data associate with a transmitter, receivers and a property U.
# The transmitter (Tx) and receiver (Rx) positions are stored separately and we
# refer them with an Index (integer). Therefore we define these fields index.
#

data = pg.DataContainer()
data.registerSensorIndex("Tx")
data.registerSensorIndex("Rx")
print(data)

# %%%
# Create a list of 10 sensors, 2m spacing
#

for x in np.arange(10):
    data.createSensor([x*2, 0])

print(data)

# %%%
# We want to use all of them (and two more!) as receivers and a constant
# transmitter of number 2.
#

data["Rx"] = np.arange(12)
# data["Tx"] = np.arange(9) # does not work as size matters!
data["Tx"] = pg.Vector(data.size(), 2)
print(data)
data.save("TxRx.data")
print(open("TxRx.data").read())

# %%%
# Again, we can mark the data validity.
#

data.markValid(data["Rx"] >= 0)
print(data["valid"])
print(data["Rx"])

# %%%
# or check the data validity automatically.
#

data.checkDataValidity()
print(data["valid"])
data.removeInvalid()
print(data)
# data.save("TxRx.data");

# %%%
# Suppose we want to compute the horizontal offset between Tx and Rx.
# We first retrieve the x position and use Tx and Rx as indices.
#

sx = pg.x(data)
data["dist"] = np.abs(sx[data["Rx"]] - sx[data["Tx"]])
print(data["dist"])

# %%%
# It might be useful to only use data where transmitter is not receiver.
#

data.markInvalid(data["Rx"] == data["Tx"])
print(data)
# data.save("TxRx.data");

# %%%
# They are still there but can be removed.
#

data.removeInvalid()
print(data)

# %%%
# At any stage we can create a new sensor
#

data.createSensor(data.sensors()[-1])
print(data)  # no change

# %%%
# , however, not at a position where already a sensor is
#

data.createSensor(data.sensors()[-1]+0.1)
print(data)
# data.save("TxRx.data")

# %%%
# Any DataContainer (indexed or not) can be visualized as matrix plot
#

pg.viewer.mpl.showDataContainerAsMatrix(data, "Rx", "Tx", "dist");

# %%%
# Instead of marking and filtering one can remove directly
#

print(data["dist"])
data.remove(data["dist"] > 11)
print(data["dist"])
print(data)

# %%%
# Similar to the nodes of a mesh, the sensor positions can be changed.
#

data.scale([2, 1])
data.translate([10, 0])
data.save("TxRx.data")

# %%%
# Suppose a receiver has not been used
#

data["Rx"][5] = data["Rx"][4]
data.removeUnusedSensors()
print(data)

# %%%
# or any measurement with it (as Rx or Tx) is corrupted
#

data.removeSensorIdx(2)
print(data)

# %%%
# There are specialized data containers with predefined indices like
# pg.DataContainerERT having indices for a, b, m and b electrodes.
# One can also add alias translators like C1, C2, P1, P2, so that
# dataERT["P1"] will return dataERT["m"]
#
