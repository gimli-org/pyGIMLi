"""Regression test for interpolating to positions outside the mesh.

Positions outside the source mesh must give the fill value (0.0 by default).
Falling back to the nearest cell (introduced with ``findNearestCell`` in
v2.0beta2) extrapolates with the shape functions instead. For the
hydrogeophysical example this produced non-positive resistivities in the
coarse boundary cells far outside the model domain, which made the ERT
forward operator log a critical error and the doc build fail.
"""

import unittest

import numpy as np

import pygimli as pg


class TestInterpolateOutside(unittest.TestCase):
    """Positions outside the mesh have to give the fill value."""

    def test_fill_value_outside_mesh(self):
        """Inside gives the data value, outside gives 0.0."""
        grid = pg.createGrid(x=[0.0, 1.0], y=[0.0, 1.0])
        u = pg.Vector(grid.nodeCount(), 1.0)

        pnts = [[0.5, 0.5], [-1.0, 0.5], [2.0, 0.5], [0.5, -1.0], [0.5, 2.0]]
        np.testing.assert_allclose(
            pg.interpolate(grid, u, pnts), [1.0, 0.0, 0.0, 0.0, 0.0]
        )

    def test_fill_value_outside_mesh_cells(self):
        """Cell data outside the mesh must not be extrapolated."""
        grid = pg.createGrid(x=[0.0, 1.0, 2.0], y=[0.0, 1.0])
        u = pg.Vector(grid.cellCount(), 100.0)
        u[0] = 1.0

        pnts = [[-5.0, 4.0], [20.0, -3.0]]
        uI = np.asarray(pg.interpolate(grid, u, pnts))
        np.testing.assert_allclose(uI, 0.0)


if __name__ == "__main__":
    unittest.main()
