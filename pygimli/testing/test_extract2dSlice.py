"""Regression tests for slicing a 3D mesh into a 2D triangle mesh.

The slice is cut by VTK, whose output contains degenerate triangles (zero
area) where the plane runs through grid nodes. Since the shape functions are
computed when a cell is created, such a triangle made the whole slice fail
with an exception from the least-squares fit in the shape function cache. In
addition the points were put into the 2D mesh with their 3D coordinates, so
every plane that is not parallel to the x/y plane collapsed to a line.

Both is covered here: every orientation has to give a usable 2D mesh whose
cells have a non-zero area and whose data fields survive. The node data is
checked against a field that follows directly from the position (distance to
the origin), and the cell data against the cell index of the source grid.

Set ``GIMLI_TEST_SHOW=1`` to have a look at the slices, e.g.::

    MPLBACKEND=TkAgg GIMLI_TEST_SHOW=1 pytest test_extract2dSlice.py
"""

import os
import unittest

import numpy as np

import pygimli as pg
from pygimli.meshtools import extract2dSlice

DX = 50
SHOW = os.environ.get("GIMLI_TEST_SHOW", "0") == "1"


class TestExtract2dSlice(unittest.TestCase):
    """Slices of a 3D grid in horizontal, vertical and diagonal direction."""

    def setUp(self):
        """Create a 1x1x0.5 km grid with node and cell data fields."""
        x = np.arange(0.0, 1001, DX)
        z = np.arange(0.0, 501, DX)
        self.grid = pg.createGrid(x, x, z)
        pos = np.asarray(self.grid.positions())
        self.grid["dist"] = np.linalg.norm(pos, axis=1)
        self.grid["px"] = pos[:, 0]
        self.grid["py"] = pos[:, 1]
        self.grid["pz"] = pos[:, 2]
        self.grid["idx"] = np.arange(self.grid.cellCount(), dtype=float)

    def extract(self, **kwargs):
        """Slice the grid and optionally show the result."""
        mesh = extract2dSlice(self.grid, **kwargs)
        if SHOW:
            ax, _ = pg.show(mesh, "dist", showMesh=True)
            ax.set_title(str(kwargs))
            pg.wait()
        return mesh

    def assertValidSlice(self, mesh, minCells=1):
        """Check that the slice is a 2D mesh with area, node and cell data."""
        self.assertEqual(mesh.dim(), 2)
        self.assertGreaterEqual(mesh.cellCount(), minCells)
        self.assertIn("dist", mesh.dataKeys())
        self.assertIn("idx", mesh.dataKeys())
        self.assertEqual(len(np.asarray(mesh["dist"])), mesh.nodeCount())
        self.assertEqual(len(np.asarray(mesh["idx"])), mesh.cellCount())
        for cell in mesh.cells():
            self.assertGreater(cell.size(), 0.0)
        # The nodes of the dropped zero-area triangles must be removed as
        # well, a mesh with nodes that belong to no cell is not usable.
        used = {node.id() for cell in mesh.cells() for node in cell.nodes()}
        self.assertEqual(len(used), mesh.nodeCount())

    def assertNodeValues(self, mesh, atol):
        """Node data is the distance to the origin of the sliced position."""
        points = np.column_stack(
            [np.asarray(mesh[key]) for key in ("px", "py", "pz")]
        )
        np.testing.assert_allclose(
            np.asarray(mesh["dist"]),
            np.linalg.norm(points, axis=1),
            atol=atol,
        )

    def test_all_orientations(self):
        """Horizontal, vertical and diagonal planes all give a 2D mesh."""
        cases = {
            "horizontal": dict(origin=[500, 500, 200], normal=[0, 0, 1]),
            "y": dict(origin=[500, 500, 0], normal=[0, 1, 0]),
            "x": dict(origin=[500, 500, 0], normal=[1, 0, 0]),
            "diagonal": dict(origin=[500, 500, 0], normal=[1, 1, 0]),
        }
        for name, kwargs in cases.items():
            with self.subTest(orientation=name):
                self.assertValidSlice(self.extract(**kwargs))

    def test_axis_keywords(self):
        """The axis keywords are shortcuts for the axis-parallel planes."""
        for kwargs in (dict(x=500), dict(y=500), dict(z=200)):
            with self.subTest(**kwargs):
                self.assertValidSlice(self.extract(**kwargs))

    def test_node_values_on_axis_planes(self):
        """Middle and border slices carry the exact distance values.

        The planes run through grid nodes, so the interpolation is exact and
        the values can be compared to the position without tolerance.
        """
        cases = {
            "middle-x": dict(origin=[500, 500, 0], normal=[1, 0, 0]),
            "middle-y": dict(origin=[500, 500, 0], normal=[0, 1, 0]),
            "middle-z": dict(origin=[500, 500, 200], normal=[0, 0, 1]),
            "border-x": dict(origin=[DX, 0, 0], normal=[1, 0, 0]),
            "border-y": dict(origin=[0, 950, 0], normal=[0, 1, 0]),
            "border-z": dict(origin=[0, 0, DX], normal=[0, 0, 1]),
        }
        for name, kwargs in cases.items():
            with self.subTest(slice=name):
                mesh = self.extract(**kwargs)
                self.assertValidSlice(mesh)
                self.assertNodeValues(mesh, atol=1e-9)

    def test_node_values_on_tilted_planes(self):
        """Tilted slices keep the position dependent node values.

        The cut points are no grid nodes here, so the values are interpolated
        and only accurate up to the interpolation error of the grid.
        """
        cases = {
            "diagonal": dict(origin=[500, 500, 0], normal=[1, 1, 0]),
            "angle": dict(origin=[500, 500, 0], normal=[0, 1, 0], angle=30),
            "dip": dict(origin=[500, 500, 0], normal=[0, 1, 0], dip=20),
            "angle-dip": dict(
                origin=[500, 500, 0], normal=[0, 1, 0], angle=30, dip=20
            ),
        }
        for name, kwargs in cases.items():
            with self.subTest(slice=name):
                mesh = self.extract(**kwargs)
                self.assertValidSlice(mesh)
                self.assertNodeValues(mesh, atol=DX)

    def test_cell_values_are_kept(self):
        """The cell data of the slice comes from the sliced cells."""
        known = np.asarray(self.grid["idx"])
        cases = {
            "horizontal": dict(origin=[500, 500, 200], normal=[0, 0, 1]),
            "vertical": dict(origin=[500, 500, 0], normal=[0, 1, 0]),
            "diagonal": dict(origin=[500, 500, 0], normal=[1, 1, 0]),
            "angle": dict(origin=[500, 500, 0], normal=[0, 1, 0], angle=30),
        }
        for name, kwargs in cases.items():
            with self.subTest(slice=name):
                mesh = self.extract(**kwargs)
                self.assertTrue(np.all(np.isin(np.asarray(mesh["idx"]), known)))

    def test_horizontal_keeps_map_view(self):
        """A horizontal slice keeps the x/y coordinates of the grid."""
        mesh = extract2dSlice(self.grid, origin=[0, 0, 200], normal=[0, 0, 1])
        self.assertAlmostEqual(min(pg.x(mesh.positions())), 0.0)
        self.assertAlmostEqual(max(pg.x(mesh.positions())), 1000.0)
        self.assertAlmostEqual(min(pg.y(mesh.positions())), 0.0)
        self.assertAlmostEqual(max(pg.y(mesh.positions())), 1000.0)

    def test_vertical_second_axis_is_depth(self):
        """A vertical slice becomes (distance along slice, z)."""
        mesh = extract2dSlice(self.grid, origin=[500, 500, 0], normal=[0, 1, 0])
        self.assertAlmostEqual(min(pg.x(mesh.positions())), -500.0)
        self.assertAlmostEqual(max(pg.x(mesh.positions())), 500.0)
        self.assertAlmostEqual(min(pg.y(mesh.positions())), 0.0)
        self.assertAlmostEqual(max(pg.y(mesh.positions())), 500.0)

    def test_plane_on_boundary_is_empty(self):
        """A plane on the domain boundary gives an empty mesh, no error."""
        mesh = extract2dSlice(self.grid, origin=[500, 500, 0], normal=[0, 0, 1])
        self.assertEqual(mesh.cellCount(), 0)
