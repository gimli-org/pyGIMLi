# -*- coding: utf-8 -*-
"""
Test for interpolation matrix.
"""
import numpy as np

import unittest

import pygimli as pg


class TestInterpolate(unittest.TestCase):

    def test_Interpolate(self):
        grid = pg.createGrid(x=[0.0, 1.0], y=[0.0, 1.0])
        u = pg.Vector(grid.nodeCount(), 1.)

        # test with pg.interpolate
        queryPos = [0.2, 0.2]
        uI = pg.interpolate(srcMesh=grid,
                            inVec=u, destPos=[queryPos,queryPos])

        np.testing.assert_allclose(uI[0], 1.)

        # test manual interpolation
        c = grid.findCell(queryPos)
        uI = c.pot(queryPos, u)
        np.testing.assert_allclose(uI, 1.)

        # test with manual interpolationMatrix generation
        I = pg.matrix.SparseMapMatrix(1, grid.nodeCount())
        cI = c.N(c.shape().rst(queryPos))
        for i in range(c.nodeCount()):
            I.addVal(0, c.node(i).id(), cI[i])

        uI = I.mult(u)
        np.testing.assert_allclose(uI[0], 1)

        # test with automatic interpolationMatrix generation
        I = grid.interpolationMatrix([[0.0, 0.0],
                                      [1.0, 0.0],
                                      [1.0, 1.0],
                                      [0.0, 1.0]])
        uI = I * u
        np.testing.assert_allclose(uI, u)


        # api test https://github.com/gimli-org/gimli/issues/131
        x = np.linspace(grid.xmin(), grid.xmax(), 11)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid), x), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x=x), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x, x*0.), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x=x, y=x*0), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x, x*0, x*0), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x=x, y=x*0, z=x*0), x)
        x = pg.Vector(x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x=x), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x, x*0.), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x=x, y=x*0), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x, x*0, x*0), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), x=x, y=x*0, z=x*0), x)

        pnts = [[0.  ,0. ],
                [0.1 ,0. ],
                [0.2 ,0. ],
                [0.3 ,0. ],
                [0.4 ,0. ],
                [0.5 ,0. ],
                [0.6 ,0. ],
                [0.7 ,0. ],
                [0.8 ,0. ],
                [0.9 ,0. ],
                [1.  ,0. ]]
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), pnts), x)
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), np.array(pnts)), x)
        
        pnts = np.linspace((min(x), 0.0), (max(x), 0.0), len(x))
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), pnts), x)

        pnts = np.linspace((min(x), 0.0, 0.0), (max(x), 0.0, 0.0), len(x))
        np.testing.assert_allclose(pg.interpolate(grid, pg.x(grid.positions()), pnts), x)


if __name__ == '__main__':
    #pg.setDebug(1)
    unittest.main()
