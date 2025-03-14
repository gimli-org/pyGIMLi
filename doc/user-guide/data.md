---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  language: python
  name: python3
---

# Data


## What is a pyGIMLi DataContainer?


Data are often organized in a data container storing the individual data values as well as any description how they were obtained, e.g. the geometry of source and receivers. DataContainer is essentially a class that stores all of this information. It is different for every method and can be then filtered accordingly. It stores two primary things: sensor information and data measurements.


A data container holds vectors like in a dictionary, however, all of them need to have the same length defined by the .size() method. Let's first go over how you can create a DataContainer from scratch. Assume we want to create and store Vertical Electrical Sounding (VES) data.

```{code-cell} ipython3
:tags: [hide-cell]

import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import VESManager
```

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

One can also use `.showInfos()` to see the content of the data container with more wording.

```{code-cell} ipython3
ves_data.showInfos()
```

As you can see from the print out there is no sensor information. In the next subsection we will explain how to add sensor information to a data container.

+++

:::{admonition} Note
:class: tip

Data containers can also be initialized from different method managers. These have the custom names for sensors and data types of each method. For example {py:class}`pygimli.physics.ert.DataContainer` already has ['a', 'b', 'm', 'n'] entries. One can also add alias translators like C1, C2, P1, P2, so that dataERT[“P1”] will return dataERT[“m”]
:::


## Creating Sensors in DataContainer


Assume we have data associate with a transmitter, receivers and a property U. The transmitter (Tx) and receiver (Rx) positions are stored separately and we refer them with an Index (integer). Therefore we define these fields index.

```{code-cell} ipython3
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
print(data)
```

:::{admonition} Note
:class: warning

The positions under the sensor indexes must be of the same size.
:::

If the sensor positions are given by another file (for example a GPS file), you can transform this to a NumPy array and set the sensor positions using [`.setSensorPositions()`](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#a3c58caec5486d5390a2b6c2d8056724f) method of the DataContainer. 

## File export

This data container can also be saved on local disk using the method `.save()` usually in a .dat format. It can then be read using `open('filename').read()`. This is also a useful way to see the data file in Python.

```{code-cell} ipython3
data.save("data.data")
print(open("data.data").read())
```

## File format import

Now we will go over the case if you have your own data and want to first import it using pyGIMLi and assign it to a data container. You can manually do this by importing data via Python (data must be assigned as Numpy arrays) and assign the values to the different keys in the data container. 

pyGIMLi also uses {py:func}`pygimli.load` that loads meshes and data files. It should handle most data types since it detects the headings and file extensions to get a good guess on how to load the data. 

Most methods also have the `load` function to load common data types used for the method. Such as, {py:func}`pygimli.physics.ert.load`. Method specific load functions assign the sensors if specified in the file. For a more extensive list of data imports please refer to [pybert importer package](http://resistivity.net/bert/_api/pybert.importer.html#module-pybert.importer).


## Processing

To start processing the data for inversion, you can filter out and analyze the data container by applying different methods available to all types of data containers. This is done to the data container and in my cases the changes happen in place, so it is recommended to view the data in between the steps to observe what changed. 

You can check the validity of the measurements using a given condition. We can mask or unmask the data with a boolean vector. For example, below we would like to mark valid all receivers that are larger or equal to 0.

```{code-cell} ipython3
data.markValid(data["Rx"] >= 0)
print(data["valid"])
print(len(data["Rx"]))
```

That adds a 'valid' entry to the data container that contains 1 and 0 values. You can also check the data validity by using `.checkDataValidity()`. It automatically removes values that are 0 in the valid field and writes the `invalid.data` file to your local directory. In this case it will remove the two additional values that were marked invalid.

```{code-cell} ipython3
data.checkDataValidity()
```

You can pass more information about the data set into the data container. For example, here we calculate the distance between transmitter and receiver.

```{code-cell} ipython3
sx = pg.x(data)
data["dist"] = np.abs(sx[data["Rx"]] - sx[data["Tx"]])
print(data["dist"])
```

+++ {"tags": ["hide-cell"]}

You can also do some pre-processing using the validity option again. For example, here we would like to mark as **invalid** where the receiver is the same as the transmitter.

```{code-cell} ipython3
data.markInvalid(data["Rx"] == data["Tx"])
data.checkDataValidity()
```

then we can remove the invalid data and see the information of the remaining data. 

```{code-cell} ipython3
data.removeInvalid()
data.showInfos()
```


Below there is a table with the most useful methods, for a full list of methods of data container, please refer to [DataContainer class reference](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html)

:::{admonition} Table of useful methods for DataContainer
:class: tip

:::{table} Methods to add, remove and sort information for a data container. 
data represents a DataContainer
:widths: auto
:align: center

| Method | Description |
| :--- | :---: |
| [data.remove()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#a3d07be3931623a5ebef2bf14d06d4f50) | Remove data from index vector. Remove all data that are covered by idx. Sensors are preserved. (Inplace - DataContainer is overwritten inplace) |
| [data.add()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#af550aeba4f21ba26fd04c9bd3a3800ac) | Add data to this DataContainer and snap new sensor positions by tolerance snap. Data fields from this data are preserved.  |
| [data.clear()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#ad5aa50883c00a3989aed97401592f41b) | Clear the container, remove all sensor locations and data.|
| [data.findSensorIndex()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#adeff1f8755b09ce80bb92b890fb0d837) | Translate a RVector into a valid IndexArray for the corresponding sensors. |
| [data.sortSensorsIndex()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#a3d07be3931623a5ebef2bf14d06d4f50) | Sort all data regarding there sensor indices and sensorIdxNames. Return the resulting permuation index array. |
| [data.sortSensorsX()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#a8d917f4e6049799bda190dfd42efef6b) | Sort all sensors regarding their increasing coordinates. Set inc flag to False to sort respective coordinate in decreasing direction.|
| [data.registerSensorIndex()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#a955c7c33ff8118ff9c3f5a7a78b75283) | Mark the data field entry as sensor index. |
| [data.removeUnusedSensor()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#a017137847f635e56a0eb5f84cbc58f5d) |  Remove all unused sensors from this DataContainer and recount data sensor index entries. |
| [data.removeSensorIdx()](https://www.pygimli.org/gimliapi/classGIMLI_1_1DataContainer.html#ab0207d2be4491338818a6c67d1ed78e3) |  Remove all data that contains the sensor and the sensor itself. |
:::



## Visualization

You can visualize the data in many ways depending on the physics manager. To simply view the data as a matrix you can use `pg.viewer.mpl.showDataContainerAsMatrix`. This visualizes a matrix of receivers and transmitters pairs with the associated data to plot : 'dist'. 

```{code-cell} ipython3ç
pg.viewer.mpl.showDataContainerAsMatrix(data, "Rx", "Tx", 'dist')
```


```{code-cell} ipython3
data.markValid(data["Rx"] >= 0)
print(data["valid"])
print(len(data["Rx"]))
```
There are more formal ways of plotting different data containers depending on the method. As seen in [Fundamentals](fundamentals.md) here we can focus on showing the data container. Most of the method managers use `.show()` that is specific to their method, however, you can always use the main function {py:func}`pg.show`, and it will detect the data type and plot it accordingly. For more information on visualizing data please refer to [data visualization](visualization.md).


