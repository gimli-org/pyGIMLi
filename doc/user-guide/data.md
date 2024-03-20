---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# User guide - Data in pyGIMLi

+++

## What is a `DataContainer`?

+++ {"editable": true, "slideshow": {"slide_type": ""}}

Data are often organized in a data container storing the individual data values as well as any description how they were obtained, e.g. the geometry of source and receivers. DataContainer is essentially a class that stores all of this information. It is different for every method and can be then filtered accordingly. It stores two primary things: sensor information and data measurements.


A data container holds vectors like in a dictionary, however, all of them need to have the same length defined by the .size() method. Let's first go over how you can create a DataContainer from scratch. Assume we want to create and store Vertical Electrical Sounding (VES) data.

```{code-cell} ipython3
---
editable: true
slideshow:
  slide_type: ''
tags: [hide-cell]
---
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import VESManager
```

+++ {"editable": true, "slideshow": {"slide_type": ""}}

We define logarithmically equidistant AB/2 spacings. Which is the midpoint of current electrode spacing for each measurement. 

```{code-cell} ipython3
ab2 = np.logspace(0, 3, 11)
print(ab2)
```

We create an empty data container. In this case not using any method yet.

```{code-cell} ipython3
ves_data = pg.DataContainer()
print(ves_data)
```

We feed it into the data container just like in a dictionary.

```{code-cell} ipython3
ves_data["ab2"] = ab2
ves_data["mn2"] = ab2 / 3
print(ves_data)
```

One can also use `showInfos()` to see the content of the data container with more wording. 

```{code-cell} ipython3
ves_data.showInfos()
```

As you can see from the print out there is no sensor information. In the next subsection we will explain how to add sensor information to a data container.

+++ {"editable": true, "slideshow": {"slide_type": ""}}

:::{admonition} Note
:class: tip

Data containers can also be initialized from different method managers. These have the custom names for sensors and data types of each method. For example `pygimli.physics.ert.DataContainer` already has ['a', 'b', 'm', 'n'] entries. One can also add alias translators like C1, C2, P1, P2, so that dataERT[“P1”] will return dataERT[“m”]
:::

+++

## Creating Sensors in DataContainer

+++ {"editable": true, "slideshow": {"slide_type": ""}}

Assume we have data associate with a transmitter, receivers and a property U. The transmitter (Tx) and receiver (Rx) positions are stored separately and we refer them with an Index (integer). Therefore we define these fields index. 

```{code-cell} ipython3
---
editable: true
slideshow:
  slide_type: ''
---
data = pg.DataContainer()
data.registerSensorIndex("Tx")
data.registerSensorIndex("Rx")
```

Then we create a list of 10 sensors with a 2m spacing. We can create sensors at any moment as long as it is not in the same position of an existing sensor.

```{code-cell} ipython3
for x in np.arange(10):
    data.createSensor([x*2, 0])

print(data)
```

We want to use all of them (and two more!) as receivers and a constant transmitter of number 2.

```{code-cell} ipython3
data["Rx"] = np.arange(12)
# data["Tx"] = np.arange(9) # does not work as size matters!
data["Tx"] = pg.Vector(data.size(), 2)
print(data.sensors())
```

:::{admonition} Note
:class: warning

The positions under the sensor indexes must be of the same size.
:::

+++

You can check the validity of the measurements using a given condition. We can mask or unmask the data with a boolean vector. For example, below we would like to mark valid all receivers that are larger or equal to 0. 

```{code-cell} ipython3
data.markValid(data["Rx"] >= 0)
print(data["valid"])
print(len(data["Rx"]))
```

That adds a 'valid' entry to the data container that contains 1 and 0 values. You can also check the data validity by using `checkDataValidity`. It automatically removes values that are out of the sensor index bounds and reports a validity check. In this case it will remove the two additional values that were added. 

```{code-cell} ipython3
data.checkDataValidity()
```

At any point we can add sensors in the data container. 

```{code-cell} ipython3
data.showInfos()
```

## File export

+++

This data container can also be saved on local disk using the method `.save()` usually in a .dat format. It can then be read using `open('filename').read()`. This is also a useful way to see the data file in Python. 

```{code-cell} ipython3
data.save("data.dat")
print(open("data.dat").read())
```

## File format import

+++

Now we will go over the case if you have your own data and want to first import it using pyGIMLi and assign it to a data container. You can manually do this by importing data via Python (data must be assigned as Numpy arrays) and assign the values to the different keys in the data container. 

However, pyGIMLi makes it easier for you using the different physics managers. For example, the `pygimli.physics.ert` has the `load()` function which supports most data formats. `pygimli.physics.em` has a `readusffile` function that reads data from single universal sounding file. Below is a table of the current loading utilities for every method

+++ {"editable": true, "slideshow": {"slide_type": ""}, "tags": ["hide-cell"]}

:::{admonition} Table of available import functions
:class: tip

:::{table} Commands to import data types depending on the physics manager
:widths: auto
:align: center

| physics manager | available loading functions |
| --- | --- |
| em | {py:func}`importMaxminData <pygimli.physics.em.importMaxminData>`, {py:func}`readusffile <pygimli.physics.em.readusffile>`, {py:func}`readHEMData <pygimli.physics.em.FDEM.readHEMData>`, {py:func}`importEmsysAsciiData <pygimli.physics.em.FDEM.importEmsysAsciiData>` |
| ert |{py:func}`load <pygimli.physics.ert.load>` |
| SIP | {py:func}`load <pygimli.physics.SIP.load>` |
| traveltime | {py:func}`load <pygimli.physics.traveltime.load>`|
| ves |  {py:func}`loadData <pygimli.physics.traveltime.loadData>` |
:::

:::

+++

## Visualization

+++

## Processing

```{code-cell} ipython3

```
