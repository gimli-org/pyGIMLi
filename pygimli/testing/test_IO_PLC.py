#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import unittest

import pygimli as pg
import numpy as np


class TestPLCIO(unittest.TestCase):

    def test_io_triangle(self):
        """
        """
        # create tempfile in most secure manner, only accessible by this process
        # id no execution allowed at all, will be deleted as soon as this
        # process stops
        try:
            import tempfile as tmp
        except ImportError:
            return

        _, name2D = tmp.mkstemp(suffix='.poly')

        # 2D, creating test trinangle poly
        m = pg.Mesh(2)
        nodes = []
        nodes.append(m.createNode([-1, -2]))
        nodes.append(m.createNode([-1, 2]))
        nodes.append(m.createNode([1, 2]))
        nodes.append(m.createNode([1, -2]))
        for i in range(4):
            m.createEdge(nodes[i], nodes[(i + 1) % 4])

        nodes.append(m.createNode([-1, -3]))
        nodes.append(m.createNode([1, -3]))

        m.createEdge(nodes[0], nodes[4])
        m.createEdge(nodes[1], nodes[5])
        m.createEdge(nodes[4], nodes[5])

        m.addRegionMarker([0., 0.], -3, area=-1)
        m.addRegionMarker([0., -2.5], 1, area=42.42)
        pg.meshtools.exportPLC(m, name2D)

        # 2D, test triangle mesh for consistincy
        poly = pg.meshtools.readPLC(name2D)

        np.testing.assert_allclose(poly.regionMarkers()[0].array(),
                                   np.array([0., -2.5, 0.]))
        np.testing.assert_allclose(poly.holeMarkers()[0].array(),
                                   np.array([0., 0., 0.]))
        np.testing.assert_allclose(np.sort(poly.positions().array()),
                                   np.sort(m.positions().array()))
        np.testing.assert_equal(poly.regionMarkers()[0].area(), 42.42)

        try:
            os.remove(name2D)
        except:
            print("can't remove:", name2D)


    def test_io_tetgen(self):
        """
        """
        try:
            import tempfile as tmp
        except ImportError:
            return

        _, name3D = tmp.mkstemp(suffix='.poly')
        #print(name3D)

        # 3D, creating test tetgen poly
        minx = -2.0
        maxx = np.pi
        miny = 0.0
        maxy = 23
        minz = -3
        maxz = 1e-4
        cube = pg.Mesh(3, isGeometry=True)

        ob0 = cube.createNode(minx, miny, minz)
        ob1 = cube.createNode(maxx, miny, minz)
        ob2 = cube.createNode(maxx, maxy, minz)
        ob3 = cube.createNode(minx, maxy, minz)
        ob4 = cube.createNode(minx, miny, maxz)
        ob5 = cube.createNode(maxx, miny, maxz)
        ob6 = cube.createNode(maxx, maxy, maxz)
        ob7 = cube.createNode(minx, maxy, maxz)

        cube.createQuadrangleFace(ob0, ob1, ob2, ob3)
        cube.createQuadrangleFace(ob1, ob2, ob6, ob5)
        cube.createQuadrangleFace(ob0, ob1, ob5, ob4)
        cube.createQuadrangleFace(ob0, ob3, ob7, ob4)
        cube.createQuadrangleFace(ob2, ob3, ob7, ob6)
        b = cube.createBoundary([n.id() for n in [ob4, ob5, ob6, ob7]])
        b.addHoleMarker([0.0, 0.0, 1.0])

        cube.addRegionMarker([0., 0., -2.0], -3, area=-1)
        cube.addRegionMarker([-1.99, 0.001, 1e-6], 1, area=42.42)
        cube.addRegionMarker([-0.99, 1, 1e-7], 1, area=1)
        cube.addRegionMarker([0.99, 11, 1e-8], 1, area=2)
        cube.addRegionMarker([1.99, 22, 1e-9], 1, area=1245535642455)

        cube.exportPLC(name3D)
        # print(len(cube.boundaries()[-1].holeMarker()))
        # cube.exportPLC('tmp')

        # 3D, test tetgen plc for validity
        poly = pg.meshtools.readPLC(name3D)

        np.testing.assert_allclose(poly.regionMarkers()[0].array(),
                                   np.array([-1.99, 0.001, 1e-6]))
        np.testing.assert_allclose(poly.holeMarkers()[0].array(),
                                   cube.holeMarkers()[0].array())
        np.testing.assert_allclose(np.sort(poly.positions().array()),
                                   np.sort(cube.positions().array()))
        np.testing.assert_equal(poly.regionMarkers()[0].area(), 42.42)
        
        np.testing.assert_equal(poly.boundaries()[-1].holeMarkers()[0], [0.0, 0.0, 1.0])
    
        try:
            os.remove(name3D)
        except:
            print("can't remove:", name3D)


    def test_io_STL(self):
        try:
            import tempfile as tmp
        except ImportError:
            return

        str = r"""solid name
    facet normal 1.000000e+000 0.000000e+000 0.000000e+000 
        outer loop 
            vertex   6.100000e+002 -1.046773e+002 -2.818708e+003 
            vertex   6.100000e+002 -1.249950e+002 -2.814064e+003 
            vertex   6.100000e+002 -2.507565e+000 -2.930000e+003 
        endloop
    endfacet
endsolid
solid name
    facet normal -2.568735e-007 1.396183e-002 -9.999025e-001
        outer loop
            vertex   1.350000e+002 3.839764e+001 -2.929429e+003
            vertex   6.100000e+002 7.930283e+001 -2.928858e+003
            vertex   6.100000e+002 -2.507565e+000 -2.930000e+003
        endloop
    endfacet
endsolid"""
        _, fileName = tmp.mkstemp(suffix='.stl')
        fi = open(fileName, 'w')
        fi.write(str)
        fi.close()

        mesh = pg.load(fileName, verbose=True)

        np.testing.assert_equal(mesh.cellCount(), 0)
        np.testing.assert_equal(mesh.nodeCount(), 5)
        np.testing.assert_equal(mesh.boundaryCount(), 2)

        np.testing.assert_equal(np.array(pg.unique(pg.sort(mesh.boundaryMarkers()))),
                               [0, 1])
        
        try:
            os.remove(fileName)
        except:
            print("can't remove:", fileName)

            
if __name__ == '__main__':
    unittest.main()
