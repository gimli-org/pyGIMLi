#!/usr/bin/env python
"""Tools for traveltime/refraction tomography."""

import numpy as np
import pygimli as pg


def createCrossholeData(sensors=None, x=None, z=None):
    """
    Create crosshole scheme assuming two boreholes with equal sensor numbers.

    Parameters
    ----------
    sensors : array (Nx2)
        Array with position of sensors. Alternatively, specify x and z
    x : iterable [M]
        horizontal positions of boreholes
    z : iterable [M]
        vertical positions of sensors in each borehole

    Returns
    -------
    scheme : DataContainer
        Data container with `sensors` predefined sensor indices 's' and 'g'
        for shot and receiver numbers.
    """
    from itertools import product
    if sensors is None:
        if x is None or z is None:
            raise ValueError("Both x and z must be given")

        sensors = np.array([(xi, zi) for zi in z for xi in x])

    if len(sensors) % 2 > 0:
        pg.error("createCrossholeData is only defined for an equal number of"
                 " sensors in two boreholes.")
    n = len(sensors) // 2
    numbers = np.arange(n)
    rays = np.array(list(product(numbers, numbers + n)))

    # Empty container
    scheme = pg.DataContainer()

    # Add sensors
    for sen in sensors:
        scheme.createSensor(sen)

    # Add measurements
    scheme.resize(len(rays))
    scheme["s"] = rays[:, 0]
    scheme["g"] = rays[:, 1]
    scheme["valid"] = np.ones(len(rays))
    scheme.registerSensorIndex("s")
    scheme.registerSensorIndex("g")
    return scheme


def shotReceiverDistances(data, full=True):
    """Distance vector between shot and for each 's' and 'g' in data.

    Parameters
    ----------
    data : pg.DataContainerERT

    full : bool [True]
        Get distances between shot and receiver position when full is True or
        only form x coordinate if full is False

    Returns
    -------
    dists :  array
        Array of distances

    """
    if full:
        pos = data.sensors()
        s, g = data.id("s"), data.id("g")
        off = [pos[s[i]].distance(pos[g[i]]) for i in range(data.size())]
        return np.absolute(off)
    else:
        px = pg.x(data)
        gx = np.array([px[g] for g in data.id("g")])
        sx = np.array([px[s] for s in data.id("s")])
        return np.absolute(gx - sx)


def createRAData(sensors, shotDistance=1):
    """Create a refraction data container.

    Default data container for shot and geophon at every sensor position.
    Predefined sensor indices's 's' as shot position and 'g' as
    geophon position.

    Parameters
    ----------
    sensors: ndarray | R3Vector
        Geophon and shot positions (same)
    shotDistances: int [1]
        Distance between shot indices.

    Returns
    -------
    data : DataContainer
        Data container with predefined sensor indices 's' and 'g'
        for shot and receiver numbers.
    """
    data = pg.DataContainer()
    data.registerSensorIndex("s")
    data.registerSensorIndex("g")

    if isinstance(sensors, np.ndarray):
        if len(sensors.shape) == 1:
            for x in sensors:
                data.createSensor([x, 0.0, 0.0])
        else:
            data.setSensorPositions(sensors)

    else:
        data.setSensorPositions(sensors)

    S, G = [], []
    for s in range(0, data.sensorCount(), shotDistance):
        for g in range(data.sensorCount()):
            if s is not g:
                S.append(s)
                G.append(g)

    data.resize(len(S))
    data.set("s", S)
    data.set("g", G)
    data.set("valid", np.abs(np.sign(data("g") - data("s"))))

    return data


def createGradientModel2D(data, mesh, vTop, vBot, flat=False):
    """Create 2D velocity gradient model.

    Creates a smooth, linear, starting model that takes the slope
    of the topography into account. This is done by fitting a straight line
    and using the distance to that as the depth value.

    Parameters
    ----------
    data: pygimli DataContainer
        The topography list is in here.
    mesh: pygimli.Mesh
        The parametric mesh used for the inversion
    vTop: float
        The velocity at the surface of the mesh
    vBot: float
        The velocity at the bottom of the mesh

    Returns
    -------
    model: pygimli Vector, length M
        A numpy array with slowness values that can be used to start
        the inversion.
    """
    xVals = pg.x(data)
    yVals = pg.y(data)
    if mesh.dim() == 3 or abs(min(yVals)) < 1e-8 and abs(max(yVals)) < 1e-8:
        yVals = pg.z(data)

    x = pg.x(mesh.cellCenters())
    z = pg.y(mesh.cellCenters())
    if mesh.dim() == 3 or abs(min(z)) < 1e-8 and abs(max(z)) < 1e-8:
        z = pg.z(mesh.cellCenters())

    pos = np.column_stack((x, z))

    if flat:
        p = np.polyfit(xVals, yVals, deg=1)  # slope-intercept form
        n = np.asarray([-p[0], 1.0])  # normal vector
        nLen = np.sqrt(np.dot(n, n))
        d = np.array([np.abs(np.dot(pos[i, :], n) - p[1]) / nLen for i
                      in range(pos.shape[0])])
    else:
        d = np.interp(x, xVals, yVals) - z

    return 1.0 / np.interp(d, [min(d), max(d)], [vTop, vBot])
