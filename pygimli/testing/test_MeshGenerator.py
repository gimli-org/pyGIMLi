#!/usr/bin/env python
"""Test cases for mesh generation."""
import sys

import unittest
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt


_show_ = False

class TestMeshGenerator(unittest.TestCase):
    """Test cases for mesh generation."""

    def testMeshAccess(self):
        """Test mesh access."""
        x = [.5, 1, 2, 3, 42]
        y = [.5, 1, 2, 3]
        z = [.5, 1, 3]

        mesh = pg.createGrid(x, y, z)
        self.assertEqual(len(mesh.nodes()), 60)
        self.assertEqual(len(mesh.cells()), 24)
        self.assertEqual(len(mesh.boundaries()), 98)


    def testMeshGenRValueProblem(self):
        """Test RValue input for mesh generation."""
        dx = 100
        x = np.arange(-1000, 1001, dx)
        z = np.arange(-800, 1, dx)
        grid = pg.meshtools.createGrid(x=x, y=x, z=z)
        self.assertEqual(len(grid.nodes()), 3969)
        self.assertEqual(len(grid.cells()), 3200)
        self.assertEqual(len(grid.boundaries()), 10320)


    def testTriangle(self):
        """Test interface to triangle mesh generation."""
        plc = pg.meshtools.createRectangle()
        mesh = pg.meshtools.createMesh(plc)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 2)
        self.assertEqual(mesh.boundaryCount(), 5)


    def _testPLC(self):
        """Test PLC merging."""
        w = mt.createCube(marker=1)
        for i, b in enumerate(w.boundaries()):
            b.setMarker(i + 1)

        c = mt.createCube(marker=2)
        for i, b in enumerate(c.boundaries()):
            b.setMarker(6 + i + 1)

        c.translate([c.xmax()-w.xmin(), 0.0])

        return mt.mergePLC3D([w, c])


    def _testTetMesh(self, mesh, plc, tgVersion='1.5'):
        """Test tet mesh."""
        if _show_:
            pg.show(mesh, showMesh=True, markers=True)

        print(f"Testing tetgen version {tgVersion} mesh generation.")
        if "1.6" in tgVersion:
            self.assertEqual(mesh.nodeCount(), 107)
            self.assertEqual(mesh.cellCount(), 351)
            self.assertEqual(mesh.boundaryCount(), 783)
        elif "1.5.1" in tgVersion:
            self.assertEqual(mesh.nodeCount(), 583)
            self.assertEqual(mesh.cellCount(), 2154)
            self.assertEqual(mesh.boundaryCount(), 4716)
        else:
            self.assertEqual(mesh.nodeCount(), 567)
            self.assertEqual(mesh.cellCount(), 2069)
            self.assertEqual(mesh.boundaryCount(), 4541)

        self.assertEqual(sorted(pg.unique(mesh.boundaryMarkers())),
                         sorted([0, *pg.unique(plc.boundaryMarkers())]))
        self.assertEqual(sorted(pg.unique(mesh.cellMarkers())),
                         sorted([1, 2]))


        for marker in pg.unique(pg.sort(plc.boundaryMarkers())):
            b1 = plc.boundaries(plc.boundaryMarkers() == marker)[0]
            b2 = mesh.boundaries(mesh.boundaryMarkers() == marker)[0]
            if b2.outside():
                np.testing.assert_array_equal(b1.norm(), b2.norm())


    def testTetgen(self, version=None):
        """Test interface to tetgen mesh generation."""
        tetgen = f"tetgen-{version}" if version else "tetgen"
        plc = self._testPLC()
        try:
            mesh = pg.meshtools.createMesh(plc,
                                           syscall=True, verbose=False,
                                           area=0.01, quality=1.12,
                                           tetgen=tetgen,
                                           )
        except RuntimeError:
            self.skipTest(f"{tetgen} binary not found in PATH")

        if version is None:
            import subprocess
            version = subprocess.getoutput(f"{tetgen} -h | grep Version").split()[-1]
            if version == '1.5':
                version = '1.5.1' ## default for testing server

        self._testTetMesh(mesh, plc, tgVersion=version)


    def testTetgen150(self):
        """Test interface to tetgen mesh generation with version 1.5."""
        self.testTetgen(version="1.5.0")


    def testTetgen151(self):
        """Test interface to tetgen mesh generation with version 1.5.1."""
        self.testTetgen(version="1.5.1")


    def testTetgen160(self):
        """Test interface to tetgen mesh generation with version 1.6.0."""
        self.testTetgen(version="1.6.0")


    def testPyTetgen(self):
        """Test interface to tetgen mesh generation with syscall=False."""
        plc = self._testPLC()
        try:
            mesh = pg.meshtools.createMesh(plc,
                                           syscall=False, verbose=False,
                                           area=0.01, quality=1.12)

        except ImportError as e:
            print(e)
            self.skipTest("tetgen python wrapper not installed")
        except OSError as e:
            if "This file was not able to be automatically read by pyvista."\
                in str(e):
                self.skipTest("tetgen wrapper probably to old.")
        except Exception as e:
            print(e)

        self._testTetMesh(mesh, plc, tgVersion="1.6")


    def testTriangleMAC(self):
        """There seems to be an issue on Mac for higher quality settings.

        which produce different meshes on mac vs. Linux/Windows.
        Needs observation and/or clarification.
        """
        plc = mt.createCircle(nSegments=24)
        l = mt.createLine(start=[0, -1], end=[0, -0.1], boundaryMarker=2)
        mesh = mt.createMesh([plc, l], area=0.1, quality=30)
        #print(mesh)
        # On Linux and Mac: Mesh: Nodes: 43 Cells: 60 Boundaries: 102
        self.assertEqual(mesh.nodeCount(), 43)

        mesh = mt.createMesh([plc, l], area=0.1, quality=32)
        #print(mesh)
        if sys.platform == "darwin":
            # On Mac: Mesh: Nodes: 46 Cells: 66 Boundaries: 111
            self.assertEqual(mesh.nodeCount(), 46)
        else:
            # On Linux Mesh: Nodes: 43 Cells: 60 Boundaries: 102
            # (same as for quality=30)
            self.assertEqual(mesh.nodeCount(), 43)


    def testCreateGrid(self):
        """Test createGrid function."""
        mesh = pg.createGrid(3)
        self.assertEqual(mesh.xMax(), 2.0)
        mesh = pg.createGrid(3, 3)
        self.assertEqual(mesh.yMax(), 2.0)
        mesh = pg.createGrid(3, 3, 3)
        self.assertEqual(mesh.zMax(), 2.0)
        # mesh = pg.meshtools.createMesh1D(10, 1)
        # print(mesh)
        # self.assertEqual(mesh.cellCount(), 10.0)

        # test also RValue conversion to correct arrays signed arrays
        x = np.arange(-10, 11)
        self.assertEqual(pg.createGrid(x,x,x).bb(),
                         [pg.Pos(-10.0, -10.0, -10.0),
                          pg.Pos(10.0, 10.0, 10.0)])


    def testCreateMesh1D(self):
        """Test createMesh1D function."""
        mesh = pg.meshtools.createMesh1D(10, 1)
        self.assertEqual(mesh.cellCount(), 10.0)
        self.assertEqual(mesh.xMax(), 10.0)

        mesh = pg.meshtools.createMesh1D(nCells=10)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh1D(nCells=5, nProperties=2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh1D(5, 2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh1D(10)
        self.assertEqual(mesh.cellCount(), 10.0)


    def testCreateMesh1DBlock(self):
        """Test createMesh1DBlock function."""
        mesh = pg.meshtools.createMesh1DBlock(nLayers=5)
        self.assertEqual(mesh.cellCount(), 9.0)

        mesh = pg.meshtools.createMesh1DBlock(5)
        self.assertEqual(mesh.cellCount(), 9.0)

        mesh = pg.meshtools.createMesh1DBlock(5, 1)
        self.assertEqual(mesh.cellCount(), 9.0)

        mesh = pg.meshtools.createMesh1DBlock(nLayers=4, nProperties=2)
        self.assertEqual(mesh.cellCount(), 11.0)

        mesh = pg.meshtools.createMesh1DBlock(4, 2)
        self.assertEqual(mesh.cellCount(), 11.0)


    def testCreateMesh2D(self):
        """Test createMesh2D function."""
        mesh = pg.meshtools.createMesh2D(xDim=5, yDim=2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh2D(5, 2)
        self.assertEqual(mesh.cellCount(), 10.0)

        mesh = pg.meshtools.createMesh2D(np.linspace(0, 1, 6),
                                         np.linspace(0, 1, 3))
        self.assertEqual(mesh.cellCount(), 10.0)

        phi = 0.025*2*np.pi
        m = pg.meshtools.createGrid([1, 2], degree=phi/(2*np.pi)*360, h=0)
        #m.show(markers=True, showMesh=True)
        for b in m.boundaries():
            if b.marker() != 0:
                if b.leftCell() is None:
                    pg._r(b)

        self.assertEqual(pg.meshtools.checkMeshConsistency(m), True)

        mesh = pg.meshtools.createMesh(pg.meshtools.createPolygon([[0.0, 0.0],
                                                                [10, 0.0],
                                                                [10, 3.0],
                                                                [0.0, 3.0]],
                                                                isClosed=True,
                                                 boundaryMarker=[1,2,3,4]),
                          marker=1)

        print(mesh.boundaryMarkers())
        mesh.show(markers=True)
        pg.wait()
        #mesh.show(markers=True, showMesh=False)




    def testCreateMesh3D(self):
        """Test createMesh3D function."""
        mesh = pg.meshtools.createMesh3D(xDim=5, yDim=3, zDim=2)
        self.assertEqual(mesh.cellCount(), 30.0)


    def testCreateMesh3DExtrude(self):
        """Test createMesh3D function with extrude option."""
        m3 = pg.meshtools.createMesh3D(xDim=1, yDim=1, zDim=1)
        self.assertEqual(pg.meshtools.checkMeshConsistency(m3), True)

        m2 = pg.meshtools.createMesh2D(xDim=2, yDim=2)
        self.assertEqual(pg.meshtools.checkMeshConsistency(m2), True)

        m3 = pg.meshtools.createMesh3D(m2, [0, 1])
        self.assertEqual(pg.meshtools.checkMeshConsistency(m3), True)
        # # print(m3)
        # print(m3)

        # for b in m3.boundaries():
        #     if b.marker() != 0:

        #         print(b.outside(), b.rightCell(), b.leftCell())
        #         if b.leftCell() is None:
        #             pg._r(b)
        #         #     #self.assertEqual(b.outside(), True)
        #         #     print(b)
        #         #     print(b.marker(), b.outside(), b.norm())
        #         #self.assertEqual(b.outside(), True)


    def testCreateMesh3DFrustums(self):
        slice = 0.025

        mesh = pg.meshtools.createFrustums([1.0, 1.2],
                                           phi=slice*2*np.pi, h=0)

        mesh.show(markers=True, showMesh=True)
        pg.meshtools.checkMeshConsistency(mesh)

        for b in mesh.boundaries():
            if b.leftCell() is None:
                print(b)


    def testCreatePartMesh(self):
        """Test createPartMesh function."""
        mesh = pg.meshtools.createMesh1D(np.linspace(0, 1, 10))
        self.assertEqual(mesh.cellCount(), 9)

        mesh2 = mesh.createMeshByCellIdx(
            pg.find(pg.x(mesh.cellCenters()) < 0.5))
        self.assertEqual(mesh2.cellCount(), 4)
        self.assertEqual(mesh2.cellCenters()[-1][0] < 0.5, True)


    def testMeshCreatePolyList(self):
        """Test createMesh with PolyList input."""
        pos = [[0, 0], [1, 0], [1, -1], [0, -1]]
        poly = pg.meshtools.createPolygon(pos, isClosed=0)
        mesh = pg.meshtools.createMesh(poly, quality=20, area=0.001)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 0)
        poly = pg.meshtools.createPolygon(pos, isClosed=1)
        mesh = pg.meshtools.createMesh(poly, quality=0, area=0.)
        self.assertEqual(mesh.nodeCount(), 4)
        self.assertEqual(mesh.cellCount(), 2)
        self.assertEqual(mesh.boundaryCount(), 5)


    def testMeshCreateSecNodes(self):
        """Test createSecondaryNodes function."""
        x = [0, 1, 2, 3, 42]
        y = [0, 1, 2, 3]
        z = [0, 1, 3]

        mesh = pg.createGrid(x, y, z)
        mesh.createSecondaryNodes(n=1)
        self.assertEqual(mesh.secondaryNodeCount(), mesh.boundaryCount() + \
                                                    (len(x)-1)*len(y)*len(z) + \
                                                    (len(y)-1)*len(z)*len(x) + \
                                                    (len(z)-1)*len(x)*len(y))

    def testMeshStr(self):
        """Test __str__ method of Mesh class."""
        mesh = pg.createGrid(2,2,2)
        #print(mesh.node(0))


    def testMeshDataAccess(self):
        """Test data access in Mesh class."""
        mesh = pg.Mesh()
        a = pg.Vector(10, 1.0)
        b = [pg.Vector(10, 1.0)]*3
        c = np.array(b).T

        mesh['a'] = a
        mesh['b'] = b
        mesh['v'] = c
        mesh['vs'] = [c, c, c]

        # pg.core.setDeepDebug(True)
        # pg.core.setDeepDebug(False)

        np.testing.assert_array_equal(mesh['a'], a)
        np.testing.assert_array_equal(mesh['b'], b)
        np.testing.assert_array_equal(mesh['v'], c)
        np.testing.assert_array_equal(mesh['vs'], [c, c, c])

        #mesh['c'] = pg.PosList(10, [1.0, 0., 0.0])


    def test_findCell(self):
        """Fix for fail """
        np.random.seed(1337)
        x = np.linspace(-2, 2, 5)
        m1 = pg.createGrid(x, x)

        m2 = pg.Mesh(m1)
        dx = m2.h()[0]
        for n in m2.nodes():
            if not n.onBoundary():
                n.setPos(n.pos() + [np.random.randn()*dx/8, np.random.randn()*dx/8, .0])

        pos = [2.0, -1.875]
        np.testing.assert_equal(m2.findCell(pos) != None, True)
        #pg._g(m2.findCell(pos) != None)

        # ax =  pg.show(m2)[0]
        # ax.plot(2.0, -1.875, 'o')

        # for c in m2.cells():
        #     if c.shape().touch(pos):
        #         print(c.shape().isInside(pos))
        #         print(c.shape().isInside(pos, True))


    def testAppendTetrahedron(self):
        """Test appendTetrahedron function."""
        grid = pg.meshtools.createGrid(5,5,5)
        mesh = pg.meshtools.appendBoundary(grid, xbound=5, ybound=5, zbound=5,
                                 isSubSurface=False)
        ax, _ = pg.show(mesh, mesh.cellMarkers(), showMesh=True,
                        filter={'clip':{}})


    def testMeshBMS(self):
        """Test mesh saving and loading in BMS format."""
        # text bms version v3 which stores geometry flag
        mesh = pg.Mesh(2, isGeometry=True)

        import tempfile as tmp
        _, fn = tmp.mkstemp()

        mesh.save(fn)
        mesh2 = pg.load(fn+'.bms', verbose=True)

        self.assertEqual(mesh.isGeometry(), mesh2.isGeometry())


    def test_meshBMSStream(self):
        """ Test mesh streaming for serialization."""
        mesh = pg.meshtools.createGrid(2)
        mesh['test node data'] = pg.Vector(mesh.nodeCount(), 33.0)
        mesh['test cell data'] = pg.Vector(mesh.cellCount(), 42.0)
        import tempfile as tmp
        _, fn = tmp.mkstemp()

        def compare(m1, m2):
            self.assertEqual(m1.dimension(), m2.dimension())
            self.assertEqual(m1.isGeometry(), m2.isGeometry())
            self.assertEqual(m1['test cell data'], m2['test cell data'])
            self.assertEqual(hash(m1), hash(m2))

        mesh.save(fn)
        mesh2 = pg.load(fn + '.bms', verbose=True)
        compare(mesh, mesh2)

        ms = mesh.serialize()
        mesh3 = pg.Mesh()
        mesh3.deserialize(ms)
        compare(mesh, mesh3)

        import pickle
        p = pickle.dumps(mesh)
        mesh4 = pickle.loads(p)
        compare(mesh, mesh4)


    def testVTKDataRead(self):
        """Test data reading from VTK files."""
        grid = pg.createGrid(np.arange(4), np.arange(3), np.arange(2))
        cM = np.arange(grid.cellCount())
        grid.setCellMarkers(cM)

        #grid['n-field'] = np.ones((4, grid.nodeCount()))
        grid['c-field'] = np.ones((4, grid.cellCount()))

        import tempfile as tmp
        _, fn = tmp.mkstemp(suffix='.vtk')


        fn = "l.vtk"
        grid.exportVTK(fn)
        mesh = pg.load(fn)

        np.testing.assert_array_equal(mesh.cellMarkers(), cM)
        np.testing.assert_array_equal(mesh['Marker'], cM)

        mesh = pg.meshtools.readMeshIO(fn)
        np.testing.assert_array_equal(mesh['Marker'], cM)

        fn = pg.getExampleFile('meshes/test_tetgen_dataCol.vtk')
        mesh = pg.load(fn)
        np.testing.assert_array_equal(mesh.cellMarkers(), cM)
        np.testing.assert_array_equal(mesh['Marker'], cM)

        mesh = pg.meshtools.readMeshIO(fn)
        np.testing.assert_array_equal(mesh['Marker'], cM)

        # pg._g('pg import vtk')
        # print(mesh)
        # print(mesh["Marker"])
        # print(mesh.cellMarkers())

        # pg._g('pg import tetgen vtk')
        # mesh = pg.load("grid1.vtk")
        # print(mesh)
        # print(mesh["Marker"])
        # print(mesh.cellMarkers())

        # pg._g('meshio import vtk')
        # print(mesh)
        # print(mesh["Marker"])
        # #print(mesh.cellMarkers())

        # pg._g('meshio import tetgen vtk')
        # mesh = pg.meshtools.readMeshIO("grid1.vtk")
        # print(mesh)
        # print(mesh["Marker"])


    def testVTKExportVTU(self):
        """Test to fix export bug."""
        mesh = pg.createGrid(4,4,4)
        mesh.exportBoundaryVTU("bounds.vtu")


    def testSimpleMeshExport(self):
        """Test simple mesh export."""
        mesh = pg.createGrid(3, 3)
        verts = mesh.positions()
        cellIds = [c.ids() for c in mesh.cells()]

        mesh2 = pg.Mesh(2)
        mesh2.createNodes(verts)
        mesh2.createCells(cellIds)

        np.testing.assert_array_equal(mesh2.nodeCount(), mesh.nodeCount())
        np.testing.assert_array_equal(mesh2.cellCount(), mesh.cellCount())


    def testRefine(self):
        """Test surface mesh refinement."""
        p = pg.meshtools.createCube()
        # test quad refine
        m = p.createH2()
        self.assertEqual(m.cellCount(), 0)
        self.assertEqual(m.nodeCount(), 26)
        self.assertEqual(m.boundaryCount(), 24)

        #mesh = pg.meshtools.createMesh(m, quality=34, area=1)

        # test tri refine
        p = pg.Mesh(3, isGeometry=True)
        # Create a simple tetrahedron geometry
        p.createNode([0, 0, 0])
        p.createNode([1, 0, 0])
        p.createNode([0, 1, 0])
        p.createNode([0, 0, 1])
        p.createBoundary([0, 1, 2])
        p.createBoundary([0, 1, 3])
        p.createBoundary([0, 2, 3])
        p.createBoundary([1, 2, 3])
        m = p.createH2()
        #print(m)
        self.assertEqual(m.nodeCount(), 10)
        self.assertEqual(m.cellCount(), 0)
        self.assertEqual(m.boundaryCount(), 16)
        #mesh = pg.meshtools.createMesh(m, quality=34, area=1)
        #self.assertEqual(m.cellCount(), 0)
        #self.assertEqual(mesh.nodeCount(), 27)
        #self.assertEqual(m.boundaryCount(), 24)

        #pg.show(m, showMesh=True, markers=True)


    def testSphere(self):
        """Test sphere generation."""
        s1 = pg.meshtools.createSphere(var='uvsphere', pos=[0,0,0])
        s2 = pg.meshtools.createSphere(var='qsphere', pos=[1,0,0],
                                       refine=3, triFaces=False)
        s3 = pg.meshtools.createSphere(var='qsphere', pos=[2,0,0],
                                       refine=3, triFaces=True)
        s4 = pg.meshtools.createSphere(var='icosphere', pos=[3,0,0],
                                       refine=2)

        pg.show([s1, s2, s3, s4], showMesh=True, markers=True)


if __name__ == '__main__':

    import sys
    if 'show' in sys.argv:
        sys.argv.remove('show')
        _show_ = True

    unittest.main()
