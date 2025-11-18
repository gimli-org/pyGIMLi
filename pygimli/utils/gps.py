# -*- coding: utf-8 -*-
"""Utility methods to read GPS data and convert them via pyproj."""

import sys
import os

from math import floor

import numpy as np

import pygimli as pg


def findUTMZone(lon, lat):
    """Find utm zone for lon and lat values.

    lon -180 -- -174 -> 1 ... 174 -- 180 -> 60
    lat < 0 hemisphere = S, > 0 hemisphere = N

    Parameters
    ----------
    lon : float
    lat : float

    Returns
    -------
    str :
        zone + hemisphere
    """
    zone = int((int(lon) + 180) / 6) + 1

    if lat > 0:
        return str(zone) + 'N'
    return str(zone) + 'S'


def getUTMProjection(zone, ellps='WGS84'):
    """Create and return the current coordinate projection.

    This is a proxy for pyproj.

    Parameters
    ----------
    utmZone : str
        Zone for for UTM

    ellipsoid : str
        ellipsoid based on ['wgs84']

    Returns
    -------
    Pyproj Projection

    """
    pyproj = pg.optImport('pyproj', 'Coordinate transformations.')

    return pyproj.Proj(proj='utm', zone=zone, ellps=ellps)

def getProjection(name, ref=None, **kwargs):
    """Syntactic sugar to get some default Projections."""
    pyproj = pg.optImport('pyproj', 'Coordinate transformations.')

    if name == 'utm':
        return getUTMProjection(**kwargs)
    elif name == 'RD83':
        return pyproj.Proj(init="epsg:4745")
    elif name == 'gk2':
        return pyproj.Proj(init="epsg:31466")
    elif name == 'gk3':
        return pyproj.Proj(init="epsg:31467")
    elif name == 'gk4':
        return pyproj.Proj(init="epsg:31468")
    elif name == 'gk5':
        return pyproj.Proj(init="epsg:31469")
    elif name == 'Soldner':
        return pyproj.Proj(init="epsg:3068")

def _getXMLData(ele, name, default):
    ret = default
    if ele.getElementsByTagName(name):
        ret = ele.getElementsByTagName(name)[0].childNodes[0].data

    return ret

def _extractWPTS(wpts):
    """Handler for Waypoints in gpx xml-dom."""
    w = dict()

    for wpt in wpts:
        if wpt.hasAttribute('lat'):
            lat = float(wpt.getAttribute('lat'))
        else:
            continue
        if wpt.hasAttribute('lon'):
            lon = float(wpt.getAttribute('lon'))
        else:
            continue

        name = _getXMLData(wpt, 'name', 'WP ' + str(len(w.keys())))

        w[name] = {'lat': lat,
                   'lon': lon,
                   'time': _getXMLData(wpt, 'time', 0),
                   'desc': _getXMLData(wpt, 'desc', None),
                   }
    return w


def readGPX(fileName):
    """Extract GPS Waypoint from GPS Exchange Format (GPX).

    Currently only simple waypoint extraction is supported.

    <gpx version="1.0" creator="creator">

    <metadata>
    <name>Name</name>
    </metadata>
    <wpt lat="51." lon="11.">
    <name>optional</name>
    <time>optional</time>
    <description>optional</description>
    </wpt>

    </gpx>

    """
    from xml.dom.minidom import parse

    dom = parse(fileName)

    name = fileName

    meta = dom.getElementsByTagName("metadata")
    if len(meta) > 0:
        name = _getXMLData(meta[0], 'name', fileName)

    wpts = dom.getElementsByTagName("wpt")

    return _extractWPTS(wpts), name


def readSimpleLatLon(filename, verbose=False):
    """Read Lat Lon coordinates from file.

    Try converting automatic.
    To be sure, provide format without "d" to ensure floating point format:

    lon lat

    or

    marker lat lon

    Returns
    -------
    list: []
        lon lat name time
    """
    def conv_(deg):
        """Convert degree into floating vpoint."""
        ret = 0.0
        if 'd' in deg:
            # 10d???
            d = deg.split('d')
            ret = float(d[0])

            if "'" in d[1]:
                # 10d23'2323.44''
                m = d[1].split("'")
                if len(m) > 1:
                    ret += float(m[0]) / 60.
                    ret += float(m[1]) / 3600.
            else:
                # 10d23.2323
                ret += float(d[1]) / 60.
        else:
            # 10.4432323
            ret = float(deg)

        return ret
    # def conv_(...):

    w = []

    with open(filename, 'r') as fi:
        content = fi.readlines()

    for line in content:
        if line[0] == '#':
            continue

        vals = line.split()

        if len(vals) == 2:  # lon lat
            w.append((conv_(vals[1]), conv_(vals[0]), '', 'time'))
        if len(vals) == 3:
            # marker lat lon
            w.append((conv_(vals[2]), conv_(vals[1]), vals[0], 'time'))

        if verbose:
            print(w[-1])

    return w


def GK2toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 2 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=2)


def GK3toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 3 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=3)


def GK4toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 4 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=4)

def GK5toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 5 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=5)


def GKtoUTM(ea, no=None, zone=32, gk=None, gkzone=None):
    """Transform any Gauss-Krueger to UTM autodetect GK zone from offset."""
    if gk is None and gkzone is None:
        if no is None:
            rr = ea[0][0]
        else:
            if isinstance(ea, list) or isinstance(ea, tuple):
                rr = ea[0]
            else:
                rr = ea

        gkzone = int(floor(rr * 1e-6))

        if gkzone <= 0 or gkzone >= 5:
            print("cannot detect valid GK zone")

    pyproj = pg.optImport('pyproj', 'coordinate transformations')
    if pyproj is None:
        return None

    gk = pyproj.Proj(init="epsg:"+str(31464+gkzone))
    wgs84 = pyproj.Proj(init="epsg:4326")  # pure ellipsoid to double transform

    utm = getUTMProjection(zone=zone, ellps='WGS84')  # UTM

    if no is None:  # two-column matrix
        lon, lat = pyproj.transform(gk, wgs84, ea[0], ea[1])
    else:
        lon, lat = pyproj.transform(gk, wgs84, ea, no)

    return utm(lon, lat)


def convddmm(num):
    """Convert numeric position into degree and minute."""
    dd = np.floor(num / 100.)
    r1 = num - dd * 100.
    return dd + r1 / 60.


def readGeoRefTIF(file_name):
    """Read geo-referenced TIFF file and return image and bbox.

    plt.imshow(im, ext = bbox.ravel()), bbox might need transform.
    """
    try:
        import gdal
        from osgeo.gdalconst import GA_ReadOnly
    except ImportError:
        sys.stderr.write("no module osgeo\n")

    dataset = gdal.Open(file_name, GA_ReadOnly)
    geotr = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    im = np.flipud(mpimg.imread(file_name))

    tifx, tify, dx = geotr[0], geotr[3], geotr[1]

    bbox = [[tifx, tifx + im.shape[1] * dx],
            [tify - im.shape[0] * dx, tify]]

    return im, bbox, projection

def readCoordsFromKML(xmlfile, zone=32, ellps="WGS84"):
    """Read coordinates from KML file and return UTM coordinates.

    Parameters
    ----------
    xmlfile : str
        XML or KML file
    zone : int [32]
        UTM zone

    Returns
    -------
    pos : np.array (Nx3)
        matrix of x, y, z positions
    """
    import utm
    import xml.etree.ElementTree as ET

    tree = ET.parse(xmlfile)
    root = tree.getroot()
    X, Y, Z = [], [], []
    for line in root.iter("*"):
        if line.tag.find("coordinates") >= 0:
            try:
                lin = line.text.replace("\n", "").replace("\t", "")
            except AttributeError:
                lin = root[0][4][2][1].text.replace("\n", "").replace("\t", "")
            for col in lin.split(" "):
                if col.find(",") > 0:
                    vals = np.array(col.split(","), dtype=float)
                    if len(vals) > 2:
                        X.append(vals[0])
                        Y.append(vals[1])
                        Z.append(vals[2])

    return np.vstack((*utm.from_latlon(X, Y), Z))
