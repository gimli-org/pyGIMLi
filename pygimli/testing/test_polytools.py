"""Tests for pygimli.meshtools.polytools
"""
import numpy as np
import unittest

import pygimli as pg
import pygimli.meshtools as mt


class TestMisc(unittest.TestCase):

    def test_appendTriangleBoundary(self):
        geom = mt.createWorld(start=[-10, 0], end=[10, -10], layers=[-5, -10])
        mesh = mt.createMesh(geom, area=1)

        mesh2 = mt.appendTriangleBoundary(mesh, marker=0)

        # test if boundary markers are preserved
        np.testing.assert_array_equal(
            pg.unique(pg.sort(mesh2.boundaryMarkers())),
            [-2, -1, 0, 2, 7, 8])
        # pg.show(mesh, markers=True)
        # pg.show(mesh2, markers=True)


class TestCreateRectangle(unittest.TestCase):
    def test_region_marker_position_basics(self):
        rect1 = mt.createRectangle(
            start=[0.0, 0.0],
            end=[2.0, -1.0],
            isClosed=True,
            marker=1,
        )
        # by default the region marker should be located at
        # minPos + (maxPos - minPos) * 0.2)
        assert rect1.regionMarkers()[0].x() == 0.4
        assert rect1.regionMarkers()[0].y() == -0.8
        assert rect1.regionMarkers()[0].marker() == 1

    def test_region_marker_position_translation_scale(self):
        rect1 = mt.createRectangle(
            pos=[1.0, 2.0],
            size=[2.0, 4],
            isClosed=True,
            marker=20,
        )
        # by default the region marker should be located at
        # minPos + (maxPos - minPos) * 0.2)
        assert rect1.regionMarkers()[0].x() == 0.4
        assert rect1.regionMarkers()[0].y() == 0.8
        assert rect1.regionMarkers()[0].marker() == 20

    def test_region_marker_position_two_ways_v1(self):
        rect1 = mt.createRectangle(
            start=[-2.0, -1.0],
            end=[2.0, 1.0],
            isClosed=True,
            marker=1,
        )

        # scaled version
        rect2 = mt.createRectangle(
            pos=[0.0, 0.0],
            size=[4.0, 2],
            isClosed=True,
            marker=2,
        )
        assert rect1.regionMarkers()[0].x() == -2 + 0.2*4
        assert rect1.regionMarkers()[0].y() == -1 + 0.2*2

        assert rect1.regionMarkers()[0].marker() == 1
        assert rect2.regionMarkers()[0].marker() == 2

        assert rect1.regionMarkers()[0] == rect2.regionMarkers()[0]

    def test_region_marker_position_two_ways_v2(self):
        rect1 = mt.createRectangle(
            start=[0.0, -1.5],
            end=[4.0, -4.0],
            isClosed=True,
            marker=1,
        )

        # scaled version
        rect2 = mt.createRectangle(
            pos=[2.0, -2.75],
            size=[4.0, 2.5],
            isClosed=True,
            marker=2,
        )
        assert rect1.regionMarkers()[0].x() == 0 + 4 * 0.2
        assert rect1.regionMarkers()[0].y() == -4.0 + 2.5 * 0.2

        assert rect1.regionMarkers()[0].marker() == 1
        assert rect2.regionMarkers()[0].marker() == 2

        assert rect1.regionMarkers()[0] == rect2.regionMarkers()[0]

    def test_region_markerposition_start_end(self):
        rect1 = mt.createRectangle(
            start=[0.0, -1.5],
            end=[4.0, -4.0],
            isClosed=True,
            marker=10,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarkers()[0].x() == 0
        assert rect1.regionMarkers()[0].y() == 0

        assert rect1.regionMarkers()[0].marker() == 10

    def test_region_markerposition_pos_size(self):
        # scaled version
        rect1 = mt.createRectangle(
            pos=[2.0, -2.75],
            size=[4.0, 2.5],
            isClosed=True,
            marker=5,
            markerPosition=[0.0, 0.0],
        )
        assert rect1.regionMarkers()[0].x() == 0
        assert rect1.regionMarkers()[0].y() == 0
        assert rect1.regionMarkers()[0].marker() == 5


class TestCreateCircle(unittest.TestCase):
    def test_default_create(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=1.0,
            marker=1
        )
        assert abs(circle.regionMarkers()[0].dist(circle.node(0).pos()) -
                   0.001) < 1e-8
        assert circle.regionMarkers()[0].marker() == 1

    def test_default_create_with_scaling(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=3.0,
            marker=6
        )
        assert abs(circle.regionMarkers()[0].dist(circle.node(0).pos()) -
                   0.001) < 1e-8
        assert circle.regionMarkers()[0].marker() == 6

    def test_create_with_custom_markerPosition(self):
        circle = mt.createCircle(
            pos=[0.0, 0.0],
            radius=3.0,
            marker=9,
            markerPosition=[2.0, -1.0],
        )
        assert circle.regionMarkers()[0].x() == 2.0
        assert circle.regionMarkers()[0].y() == -1
        assert circle.regionMarkers()[0].marker() == 9


class TestCreatePolygon(unittest.TestCase):
    def test_default_create(self):
        polygon = mt.createPolygon(
            [[-1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1]],
            isClosed=True,
            marker=1,
        )

        assert abs(polygon.regionMarkers()[0].dist(polygon.node(0).pos()) -
                   0.001) < 1e-8
        assert polygon.regionMarkers()[0].marker() == 1

    def test_default_create_different_center(self):
        polygon = mt.createPolygon(
            [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
            isClosed=True,
            marker=2,
        )
        assert abs(polygon.regionMarkers()[0].dist(polygon.node(0).pos()) -
                   0.001) < 1e-8
        assert polygon.regionMarkers()[0].marker() == 2

    def test_create_with_custom_markerPosition(self):
        polygon = mt.createPolygon(
            [[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1]],
            isClosed=True,
            marker=7,
            markerPosition=[0.1, -0.1],
        )

        assert polygon.regionMarkers()[0].x() == 0.1
        assert polygon.regionMarkers()[0].y() == -0.1
        assert polygon.regionMarkers()[0].marker() == 7


class Test3DMerge(unittest.TestCase):

    @pg.skipOnDefaultTest
    def test_PosTouchFace(self):
        m = pg.Mesh(dim=3, isGeometry=True)
        D = 0.1
        H1 = 2.
        H = 0.1

        m.createNode([-D/2, +D/2, -(H1-H)/2 - H1/2])
        m.createNode([+D/2, +D/2, -(H1-H)/2 - H1/2])
        m.createNode([+D/2, +D/2, +(H1-H)/2 - H1/2])
        m.createNode([-D/2, +D/2, +(H1-H)/2 - H1/2])
        m.createNode([-D/2, -D/2, -(H1-H)/2 - H1/2])
        m.createNode([+D/2, -D/2, -(H1-H)/2 - H1/2])
        m.createNode([+D/2, -D/2, +(H1-H)/2 - H1/2])
        m.createNode([-D/2, -D/2, +(H1-H)/2 - H1/2])
        m.createBoundary([0, 1, 2, 3])

        n = [-D/6, D/2, -H/2]
        m.createNode(n)
        self.assertEqual(m.boundary(0).nodeCount(), 5)

    def test_cubeBasics(self):
        plc = mt.createCube()
        for i, b in enumerate(plc.boundaries()):
            b.setMarker(i+1)

        mesh = mt.createMesh(plc)

        for marker in pg.unique(pg.sort(plc.boundaryMarkers())):
            b1 = plc.boundaries(plc.boundaryMarkers() == marker)[0]
            b2 = mesh.boundaries(mesh.boundaryMarkers() == marker)[0]

            np.testing.assert_array_equal(b1.norm(), b2.norm())

    def test_cube_cube_same(self):
        c1 = mt.createCube()
        c2 = mt.createCube()
        m = mt.mergePLC3D([c1, c2])
        self.assertEqual(c1.nodeCount(), m.nodeCount())
        self.assertEqual(c1.boundaryCount(), m.boundaryCount())

    def test_cube_cube_equalface(self):
        w = mt.createCube(marker=1)
        c = mt.createCube(marker=2)
        c.translate([c.xmax()-w.xmin(), 0.0])

        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+4)
        self.assertEqual(w.boundaryCount(), 6+5)

        c = mt.createCube(marker=3)
        c.translate([0.0, w.ymax()-c.ymin(), 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+4+4)
        self.assertEqual(w.boundaryCount(), 6+5+5)

        c = mt.createCube(marker=4)
        c.translate([0.0, 0.0, c.zmax()-w.zmin()])
        w = mt.mergePLC3D([c, w])
        self.assertEqual(w.nodeCount(), 8+4+4+4)
        self.assertEqual(w.boundaryCount(), 6+5+5+5)

        c = mt.createCube(marker=5)
        c.translate([0.0, w.ymax()-c.ymin(), c.zmax()-w.zmin()])
        w = mt.mergePLC3D([c, w])
        self.assertEqual(w.nodeCount(), 8+4+4+4+6)
        self.assertEqual(w.boundaryCount(), 6+5+5+5+6)

        c = mt.createCube(marker=6)
        c.translate([0.0, c.ymax()-w.ymin(), c.zmax()-w.zmin()])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+4+4+4+6+0)
        self.assertEqual(w.boundaryCount(), 6+5+5+5+6+3)

        # w.exportPLC('t.poly')
        # pg.show(mt.createMesh(w))

    def test_cube_cube_coplanar_touchface(self):
        w = mt.createCube(marker=1)
        w.scale([2.0, 2.0, 2.0])

        c = mt.createCube(marker=2)
        c.translate([1.5, 0.0, 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8)
        self.assertEqual(w.boundaryCount(), 6+5)

        c = mt.createCube(marker=3)
        c.translate([-1.5, 0.0, 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8+8)
        self.assertEqual(w.boundaryCount(), 6+5+5)

        c = mt.createCube(marker=4)
        c.translate([0.0, 1.5, 0.0])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8+8+8)
        self.assertEqual(w.boundaryCount(), 6+5+5+5)

        c = mt.createCube(marker=5)
        c.translate([0.0, 0.0, -1.5])
        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8+8+8+8)
        self.assertEqual(w.boundaryCount(), 6+5+5+5+5)

        # pg.show(w)
        # w.exportPLC('t.poly')
        # pg.show(mt.createMesh(w))

    def test_smallcube_in_bigcube(self):
        """
        A small cube in a bigger one, creating two subfaces.
        author: @frodo4fingers
        """
        w = mt.createCube(marker=1)
        c = mt.createCube(size=[0.5, 1.0, 1.0], marker=2)

        w = mt.mergePLC3D([w, c])
        self.assertEqual(w.nodeCount(), 8+8)
        self.assertEqual(w.boundaryCount(), 8)

        # will not work until edge intersection is working
        # d = mt.createCube(size=[0.8, 1.0, 1.0],
        #                   pos=[0.1, 0.0, 1.0],
        #                   marker=3)
        # w = mt.mergePLC3D([w, d])
        # self.assertEqual(w.nodeCount(), 8+8)
        # self.assertEqual(w.boundaryCount(), 8)

        # print(w)
        # pg.show(w)
        # pg.show(mt.createMesh(w))

    def test_cyl_on_cyl(self):
        # merge only works if smaller face merged into larger on contact plane
        segs = 12
        c1 = mt.createCylinder(radius=2, marker=1,
                               nSegments=segs, boundaryMarker=1)
        c2 = mt.createCylinder(radius=1, marker=2,
                               nSegments=segs, boundaryMarker=2)
        c1.translate([0, 0, 0.5])
        c2.translate([0, 0, -0.5])

        w = mt.mergePLC3D([c1, c2])

        self.assertEqual(w.nodeCount(), segs*2 * 2)
        self.assertEqual(w.boundaryCount(), segs*2 + 3)

        # w.exportBoundaryVTU('w')
        # pg.show(w)
        m = mt.createMesh(w)
        pg.show(m, m.cellMarkers())
        # pg.show(mt.createMesh(w))

    def test_face_in_face(self):
        """Test subface with different marker constructed with hole marker."""
        w = mt.createCube(marker=1, boundaryMarker=1)
        w.translate([0.5, 0.5, 0.5])
        b = w.boundary(2)

        pad = mt.createFacet(mt.createCircle(radius=0.2, nSegments=4,
                                             isHole=True))
        b2 = pad.boundary(0)

        # rotate to match target norm and pos
        rot = pg.core.getRotation(b2.norm(), b.norm())
        pad.transform(rot)
        pad.translate(b.center())

        # create a boundary with new marker as hole
        w.copyBoundary(b2)

        # fill the hole with new boundary
        w.createBoundary(w.nodes([w.createNode(n.pos()).id()
                                  for n in b2.nodes()]), marker=2)

        # print(w.boundaryMarkers())
        # pg.show(w, markers=True, showMesh=True)
        # w.exportPLC('pad.poly')

        mesh = mt.createMesh(w)
        # print(mesh)
        # mesh.exportBoundaryVTU('b.vtu')
        pg.show(mesh, markers=True, showMesh=True)

        # mesh.exportBoundaryVTU('b.vtu')
        self.assertEqual(w.boundaryCount(), 7)

        np.testing.assert_array_equal(
            pg.unique(pg.sort(mesh.boundaryMarkers())),
            [0, 1, 2])

    def test_cube_cube_halfside(self):
        """Add half size cube on another cube"""
        D = 1
        H = 1
        W = 3
        c1 = mt.createCube([D, D, H/2], pos=[0.0, 0.0, H/4])
        c1L = pg.Mesh(c1)
        c1L.addRegionMarker([0.0, 0.0, 0.0], marker=2)
        # c1R = pg.Mesh(c1L)
        # c1R.translate([0.0, W-D, 0.0])

        c2 = mt.createCube([D, W-2*D, H], pos=[0, -(W)/2+D/2-0.0, 0], marker=2)
        # plc = mt.merge([c2, c1L, c1R])
        plc = mt.merge([c2, c1L])
        # plc.exportPLC('cubecut')
        pg.show(plc)

        self.assertEqual(plc.nodeCount(), 14)
        self.assertEqual(plc.boundaryCount(), 11)

        m = mt.createMesh(plc)
        # m = mt.createMesh(plc, tetgen='tetgen-1.4.3')
        pg.show(m, m.cellMarkers())

    @pg.skipOnDefaultTest
    def test_cube_cube_thirdside(self):
        """Add half size cube on another cube"""
        D = 0.1
        H1 = 2
        H = 0.1

        c1 = mt.createCube([D, D, H1-H], marker=1)
        c2 = mt.createCube([D/3, D, H], marker=2)
        c1.translate([0.0, 0.0, -H1/2])

        # c1 = mt.createCube([D, D, D], pos=[0.0, 0.0, 0.0], marker=1)
        # c2 = mt.createCube([D, D/4, D], pos=[0, 0.0, D], marker=2)

        plc = mt.merge([c1, c2])
        plc.exportPLC('ccq')
        # pg.show(plc)

        self.assertEqual(plc.nodeCount(), 16)
        self.assertEqual(plc.boundaryCount(), 11)
        self.assertEqual(plc.boundary(2).nodeCount(), 8)
        self.assertEqual(plc.boundary(4).nodeCount(), 6)
        self.assertEqual(plc.boundary(5).nodeCount(), 6)

        m = mt.createMesh(plc)
        # m = mt.createMesh(plc, tetgen='tetgen-1.4.3')
        pg.show(m, m.cellMarkers())

    def test_cube_cut_cube(self):
        """Test cube cut from cube on two neighbour faces."""
        c1 = mt.createCube(marker=1)
        c2 = mt.createCube(marker=2)
        c2.scale([0.5, 0.5, 0.5])
        c2.translate([0.25, -0.25, 0.25])

        plc = mt.mergePLC3D([c1, c2])

        self.assertEqual(plc.nodeCount(), 15)
        self.assertEqual(plc.boundaryCount(), 9)

        c3 = mt.createCube(isHole=True, marker=3)
        c3.scale([0.3, 0.3, 0.3])
        c3.translate([0.35, 0.35, 0.35])
        plc = mt.mergePLC3D([plc, c3])

        c4 = mt.createCube(isHole=False, marker=4)
        c4.scale([0.2, 0.2, 1.0])
        c4.translate([-0.4, 0.0, 0.0])
        plc = mt.mergePLC3D([plc, c4])

        c5 = mt.createCube(isHole=False, marker=5)
        c5.scale([0.8, 0.2, 0.2])
        c5.translate([0.1, 0.0, -0.4])
        plc = mt.mergePLC3D([c1, c4, c5])

        # plc.exportPLC('cubecut')
        # pg.show(plc)

        # m = mt.createMesh(plc, tetgen='tetgen-1.4.3')
        m = mt.createMesh(plc)

        pg.show(m, m.cellMarkers())

        # m.exportVTK('cubecut')

    @pg.skipOnDefaultTest
    def test_cube_mult_cut1(self):
        c = mt.createCube(marker=1)
        c.translate([0.5, 0.5, 0.5])
        c1 = mt.createCube(size=[0.1, 1, 0.1], marker=2)
        c1.translate([0.5 - 0.05, 0, 0])
        c1.translate([0.5, 0.5, 0.5])

        c3 = mt.createCube(size=[1-2*0.1, 0.1, 0.1], marker=3)
        c3.translate([0, -0.5 + 0.05, 0])
        c3.translate([0.5, 0.5, 0.5])

        c = mt.merge([c, c1, c3])

        self.assertEqual(c.nodeCount(), 22)
        self.assertEqual(c.boundaryCount(), 13)
        c.exportPLC('c')

        # pg.show(c, showMesh=True)
        m = mt.createMesh(c)
        self.assertEqual(m.nodeCount() > 0, True)
        # pg.show(m, m.cellMarkers())

    @pg.skipOnDefaultTest
    def test_cube_mult_cut2(self):
        c = mt.createCube(marker=1)
        c.translate([0.5, 0.5, 0.5])
        c1 = mt.createCube(size=[0.1, 1, 0.1], marker=2)
        c1.translate([0.5 - 0.05, 0, 0])
        c1.translate([0.5, 0.5, 0.5])

        c2 = mt.createCube(size=[0.1, 1, 0.1], marker=2)
        c2.translate([-0.5 + 0.05, 0, 0])
        c2.translate([0.5, 0.5, 0.5])

        c3 = mt.createCube(size=[1-2*0.1, 0.1, 0.1], marker=3)
        c3.translate([0, -0.5 + 0.05, 0])
        c3.translate([0.5, 0.5, 0.5])

        c = mt.merge([c, c1, c2, c3])

        c.exportPLC('c')
        self.assertEqual(c.nodeCount(), 28)
        self.assertEqual(c.boundaryCount(), 15)

        # pg.show(c, showMesh=True)
        m = mt.createMesh(c)
        pg.show(m, m.cellMarkers())

        c.exportPLC('c')

        # pg.show(c, showMesh=True)
        m = mt.createMesh(c)
        self.assertEqual(m.nodeCount() > 0, True)

    def test_face_inCube(self):
        # plc = mt.createCube()

        # face = mt.createSurface(mt.createGrid(x=np.linspace(-0.5, 0.5, 2),
        #                                       y=np.linspace(-0.5, 0.5, 2)))

        # m = mt.mergePLC3D([plc, face])
        # mt.exportPLC(m, 'tmp.poly')
        # mesh = mt.createMesh(m, verbose=True)
        # pg.show(mesh)
        #
        # incomplete!!! Fix me, tetgen segfaults
        pass

    def test_appendTetrahedron(self):
        grid = mt.createGrid(5, 5, 5)
        # pg.show(grid)
        mesh = mt.appendBoundary(grid, xbound=5, ybound=5, zbound=5,
                                 isSubSurface=False)
        # ax, _ = pg.show(mesh, mesh.cellMarkers(), hold=True, opacity=0.5)

        # try:
        #     fallback mpl axes have no show
        #     ax.show()
        # except:
        #     pass

    @pg.skipOnDefaultTest
    def test_patch_on_facet(self):
        world = mt.createWorld(start=[-1, -1, -1], end=[1, 1, 0])

        pad = mt.createFacet(mt.createCircle(pos=[0, 0, 0], radius=0.5, isHole=True))
        plc = (world + pad)
        plc.createBoundary(plc.nodes([plc.createNode(n.pos()).id() for n in pad.nodes() ]), marker=2)

        # instead of the above  -- IMPLEMENTME
        #pad = mt.createFacet(mt.createCircle(pos=[0, 0, 0], radius=0.5, marker=2))
        #plc = world + pad

        mesh = mt.createMesh(plc)
        np.testing.assert_array_equal(np.unique(mesh.boundaryMarkers()),
                                     [-2, -1,  0, 2])
        #pg.show(mesh, markers=True, showMesh=True)

        plc = mt.createWorld(start=[-1, -1, -1], end=[2, 1, 0])
        for i in range(4):
            pad = mt.createFacet(mt.createCircle(pos=[i*0.5, 0, 0], radius=0.1, isHole=True))
            plc += pad
            plc.createBoundary([plc.createNode(n.pos()).id() for n in pad.nodes()], marker=i+1)

        mesh = mt.createMesh(plc)

        #pg.show(mesh, showMesh=True, markers=True, backend="trame")
        ## colorbar wrong - FIXME
        np.testing.assert_array_equal(np.unique(mesh.boundaryMarkers()),
                                     [-2, -1,  0,  1,  2,  3,  4])



if __name__ == '__main__':
    unittest.main()
