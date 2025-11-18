# -*- coding: utf-8 -*-
"""Extensions to the core DataContainer class[es]."""
import numpy as np
from . logger import critical, verbose
from .core import (RVector, RVector3, DataContainer, DataContainerERT)
from .core import (yVari, zVari, swapXY, swapYZ, x, y, z)


def __DataContainer_str(self):
    return "Data: Sensors: " + str(self.sensorCount()) + " data: " + \
        str(self.size()) + ", nonzero entries: " + \
        str([d for d in self.dataMap().keys() if self.isSensorIndex(d) or
             self.haveData(d)])


DataContainer.__repr__ = __DataContainer_str
DataContainer.__str__ = __DataContainer_str


def __DataContainer_setSensors(self, sensors):
    """Set Sensor positions.

    Set all sensor positions.
    This is just syntactic sugar for setSensorPositions.

    Parameters
    ----------
    sensors: iterable
        Iterable that can be converted into a pg.Pos.

    Tests
    -----
    >>> import pygimli as pg
    >>> d = pg.DataContainerERT()
    >>> d.setSensors(pg.utils.grange(0.0, 3, n=4))
    >>> assert d.sensorCount() == 4
    """
    for i, s in enumerate(sensors):
        nS = s
        if isinstance(s, float) or isinstance(s, int):
            nS = RVector3(s, 0.0)

        if i > self.sensorCount():
            self.createSensor(nS)
        else:
            self.setSensorPosition(i, nS)


DataContainer.setSensors = __DataContainer_setSensors


def __DataContainer_copy(self):
    return type(self)(self)


DataContainer.copy = __DataContainer_copy


def __DC_setVal(self, key, val):
    """Set datacontainer values for specific token: data[token] = x."""
    if isinstance(val, (float, int)):
        val = RVector(self.size(), val)

    if len(val) > self.size():
        verbose("DataContainer resized to:", len(val))
        self.resize(len(val))

    self.set(key, val)


DataContainer.__setitem__ = __DC_setVal


def __DC_getVal(self, key):
    if self.isSensorIndex(key):
        return np.array(self(key), dtype=int)
    # return self(key).array() // d['a'][2] = 0.0, would be impossible
    return self(key)


DataContainer.__getitem__ = __DC_getVal


def __DataContainer_ensure2D(self):
    sen = self.sensors()
    if ((zVari(sen) or max(abs(z(sen))) > 0) and
            (not yVari(sen) and max(abs(y(sen))) < 1e-8)):
        swapYZ(sen)
        self.setSensorPositions(sen)


DataContainer.ensure2D = __DataContainer_ensure2D


def __DataContainer_swapXY(self):
    sen = self.sensors()
    swapYZ(sen)
    self.setSensorPositions(sen)


DataContainer.swapXY = __DataContainer_swapXY


def __DataContainerERT_addFourPointData(self, *args,
                                        indexAsSensors=False, **kwargs):
    """Add a new data point to the end of the dataContainer.

    Add a new 4 point measurement to the end of the dataContainer and increase
    the data size by one. The index of the new data point is returned.

    Parameters
    ----------
    *args: [int]
        At least four index values for A, B, M and N.
    indexAsSensors: bool [False]
        Indices A, B, M and N are interpreted as sensor position in [m, 0, 0].
    **kwargs: dict
        Named values for the data configuration.

    Returns
    -------
    ret: int
        Index of this new data point.

    Examples
    --------
    >>> import pygimli as pg
    >>> d = pg.DataContainerERT()
    >>> d.setSensors(pg.utils.grange(0, 3, n=4))
    >>> d.addFourPointData(0,1,2,3)
    0
    >>> d.addFourPointData([3,2,1,0], rhoa=1.0)
    1
    >>> print(d)
    Data: Sensors: 4 data: 2, nonzero entries: ['a', 'b', 'm', 'n', 'rhoa', 'valid']
    >>> print(d('rhoa'))
    2 [0.0, 1.0]
    """
    try:
        if len(args) == 1:
            a, b, m, n = args[0][:]
        else:
            [a, b, m, n] = args

        if indexAsSensors:
            a = self.createSensor([float(a), 0.0, 0.0])
            b = self.createSensor([float(b), 0.0, 0.0])
            m = self.createSensor([float(m), 0.0, 0.0])
            n = self.createSensor([float(n), 0.0, 0.0])
        idx = self.createFourPointData(self.size(), a, b, m, n)

    except Exception as e:
        print(e)
        print("args:", args, len(args))
        critical("Can't interpret arguments:", *args)

    for k, v in kwargs.items():
        if not self.haveData(k):
            self.add(k)
        self.ref(k)[idx] = v
    return idx


DataContainerERT.addFourPointData = __DataContainerERT_addFourPointData

def __DataContainer_show(self, *args, **kwargs):
    """Use data.show(**) instead of pg.show(data, *) syntactic sugar."""
    import pygimli as pg
    return pg.show(self, *args, **kwargs)

DataContainer.show = __DataContainer_show


def __DataContainer_getIndices(self, **kwargs):
    """Return indices for all data keys equalling values."""
    good = np.ones(self.size(), dtype=bool)
    for k, v in kwargs.items():
        good = np.bitwise_and(good, self[k] == v)

    return np.nonzero(good)[0]

DataContainer.getIndices = __DataContainer_getIndices


def __DataContainer_removeData(self, **kwargs):
    self.markInvalid(self.getIndices(**kwargs))
    self.removeInvalid()

DataContainer.removeData = __DataContainer_removeData


def __DataContainer_subset(self, **kwargs):
    """Return a subset for which all kwarg conditions hold.

    Parameters
    ----------
    data : DataContainer
        pyGIMLi data container or derived class
    kwargs : dict
        dictionary forwarded to getIndices marking validity
    x/y/z : float
        positions are extended to all sensor indices

    Returns
    -------
    out : DataContainer
        filtered data container (of the same class)
    """
    new = self.copy()
    new["valid"] = 0
    remSen = False
    for xyz in "xyz":
        xx = kwargs.pop(xyz, None)
        if xx is not None:
            remSen = True
            ex = eval(xyz)(new.sensorPositions())
            for key in new.dataMap().keys():
                if new.isSensorIndex(key):
                    new[xyz+key] = ex[new[key]]
                    kwargs[xyz+key] = xx

    new.markValid(new.getIndices(**kwargs))
    new.removeInvalid()
    if remSen:
        new.removeUnusedSensors()

    return new

DataContainer.subset = __DataContainer_subset

def __DataContainer_markSensorInvalid(self, idx):
    """Mark all measurements using a specific sensor invalid."""
    for tok in self.dataMap().keys():
        if self.isSensorIndex(tok):
            self.markInvalid(self[tok] == idx)

DataContainer.markSensorInvalid = __DataContainer_markSensorInvalid

def __DataContainer_removeSensorData(self, idx):
    """Remove all measurements using a specific sensor."""
    self.markSensorInvalid(idx)
    self.removeInvalid()  # TODO4Nico: do it correctly
    ## It will also remove any previously invalid measurements
    ## So you better use only remove or markValid

DataContainer.removeSensorData = __DataContainer_removeSensorData

def __DataContainer_removeSensor(self, idx):
    """Remove a specific sensor."""
    self.removeSensorData(idx)
    self.removeUnusedSensors()
    ## Same argument as above, there might be other unused ones

DataContainer.removeSensor = __DataContainer_removeSensor
