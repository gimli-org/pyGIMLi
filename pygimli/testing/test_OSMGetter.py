#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

from pygimli.viewer.mpl.overlayimage import underlayMap
#from pygimli.io import readGPX, readSimpleLatLon

import matplotlib.pyplot as plt
#import pyproj

class TestOSMGetter(unittest.TestCase):

    def test_OSM(self):
        #zone = 32
        #vendor = 'OSM'
        #pnts = readGPX( 'rieseltag.gpx' )

        #zone = 31
        #vendor = 'GM'
        #pnts = readSimpleLatLon('bleicherode.dat', True)

        ##zone = 29
        ##pnts = readGPX( 'gps.gpx' )

        #proj = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')

        ## gauss krüger projektion 4. Meridianstreifen (GK 4), Ellipsoid = "Bessel"
        #proj = pyproj.Proj(proj='tmerc', init='epsg:31468', ellps='WGS84')
        ##proj = pyproj.Proj( proj = 'tmerc', init = 'epsg:31468', ellps = 'bessel' )

        #for p in pnts:
            #x, y = proj(p[0], p[1])
            #print(x, y)
            #plt.plot(x, y, '.', markersize=18)
            #plt.text(x, y, p[2])

        #axes = plt.gca()

        #underlayMap(
            #axes,
            #proj,
            #vendor=vendor,
            #zoom=10,
            #pixelLimit=[
                #1024,
                #1024],
            #verbose=True,
            #fitMap=True)

        #axes.grid()

        #plt.show()
        pass

if __name__ == '__main__':

    unittest.main()

