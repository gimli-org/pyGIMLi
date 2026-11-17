#!/usr/bin/env python
"""Test the OSM getter helper (currently without active test code)."""

# write a correct test!
import unittest

# from pygimli.io import readGPX, readSimpleLatLon

# import pyproj


class TestOSMGetter(unittest.TestCase):
    """Placeholder for OSM/map getter tests."""

    def test_OSM(self):
        """Keep a placeholder test that does nothing."""
        # zone = 32
        # vendor = 'OSM'
        # pnts = readGPX( 'rieseltag.gpx' )

        # zone = 31
        # vendor = 'GM'
        # pnts = readSimpleLatLon('bleicherode.dat', True)

        ##zone = 29
        ##pnts = readGPX( 'gps.gpx' )

        # proj = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')

        ## Gauss-Krueger projection, 4th meridian strip (GK 4),
        ## ellipsoid = "Bessel"
        # proj = pyproj.Proj(proj='tmerc', init='epsg:31468', ellps='WGS84')
        ##proj = pyproj.Proj( proj = 'tmerc', init = 'epsg:31468',
        ##                   ellps = 'bessel' )

        # for p in pnts:
        # x, y = proj(p[0], p[1])
        # print(x, y)
        # plt.plot(x, y, '.', markersize=18)
        # plt.text(x, y, p[2])

        # axes = plt.gca()

        # underlayMap(
        # axes,
        # proj,
        # vendor=vendor,
        # zoom=10,
        # pixelLimit=[
        # 1024,
        # 1024],
        # verbose=True,
        # fitMap=True)

        # axes.grid()

        # plt.show()
        pass


if __name__ == "__main__":
    unittest.main()
