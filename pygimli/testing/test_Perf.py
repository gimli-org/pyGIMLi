#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest
import timeit

import pygimli as pg
import numpy as np

setup = """\
import pygimli as pg
import numpy as np

N=10001

np1 = np.linspace(1.1, 1.1, N)
np2 = np.linspace(2.1, 2.1, N)

pg1 = pg.Vector(np1)
pg2 = pg.Vector(np2)

st = range(N)

"""

class TestPerf(unittest.TestCase):

    def test_Performance(self):
        """
        """
        #pg.setDebug(True)

        sw = pg.Stopwatch(True)
        #print(timeit.repeat('r = grid.cellSizes() * np1', setup=setup, number=1000))
        #print(timeit.repeat('r = c * np1', setup=setup, number=1000))

        print(("np(np)", timeit.repeat('np.array(np1)', setup=setup, number=5000)))
        print(("pg(pg)", timeit.repeat('pg.Vector(pg1)', setup=setup, number=5000)))
        print(("pg(np)", timeit.repeat('pg.Vector(np1)', setup=setup, number=5000)))
        print(("np(pg)", timeit.repeat('np.array(pg1)', setup=setup, number=5000)))
        print(("np * np", timeit.repeat('np3 = np1 * np2', setup=setup, number=5000)))
        print(("pg * pg", timeit.repeat('pg3 = pg1 * pg2', setup=setup, number=5000)))
        print(("pg * np", timeit.repeat('pg1 * np1', setup=setup, number=5000)))
        print(("np * pg", timeit.repeat('np1 * pg1', setup=setup, number=5000)))
        print(("sum(np)", timeit.repeat('sum(np1)', setup=setup, number=300)))
        print(("sum(pg)", timeit.repeat('sum(pg1)', setup=setup, number=300)))
        print(("pg.sum(pg)", timeit.repeat('pg.sum(pg1)', setup=setup, number=300)))
        print(("pg.sum(st)", timeit.repeat('pg.sum(st)', setup=setup, number=300)))

        print(("s", sw.duration(True)))

        N = 10001
        np1 = np.linspace(1.1, 1.2, N)
        np2 = np.linspace(2.1, 2.1, N)

        pg1 = pg.Vector(np1)
        pg2 = pg.Vector(np2)

        # print(sw.duration(True))
        print((sum(np1 * np1)))
        print((sum(pg1 * pg1)))
        print((sum(np1 * pg1)))
        print((sum(pg1 * np1)))


class TestMT(unittest.TestCase):

    def test_WayMatrix(self):
        """ MT-OMP Test
        run with env OMP_NUM_THREADS=4 GIMLI_NUM_THREADS=3 PG_USE_OMP=1

        OMP THREADS if USE_OMP=1 else DISTR_CALC with GIMLI_NUM_THREADS

        USE_OMP=0 and GIMLI_NUM_THREADS > 1 segfaults on py312
        """
        from pygimli.physics import traveltime as tt

        x, y = 25, 25
        mesh = pg.createGrid(x, y)
        data = tt.createRAData([(0, 0)] + [(x, i) for i in range(y)],
                               shotDistance=y+1)
        data["t"] = 1.0
        mgr = tt.Manager(data, verbose=False)
        mgr.applyMesh(mesh, secNodes=10)
        pg.tic()
        mgr.fop.createJacobian(np.full(mgr.fop.parameterCount, 1))
        pg.toc()


    def test_ERTJacobian(self):
        """ MT-OMP Test
        run with env OMP_NUM_THREADS=4 GIMLI_NUM_THREADS=3 PG_USE_OMP=1

        OMP THREADS if USE_OMP=1 else DISTR_CALC with GIMLI_NUM_THREADS

        USE_OMP=0 and GIMLI_NUM_THREADS > 1 segfaults on py312
        """
        from pygimli.physics import ert

        dat = pg.getExampleFile('ert/gallery.dat', load=True, verbose=True)
        dat['k'] = ert.createGeometricFactors(dat)
        mesh = pg.meshtools.createParaMesh(dat.sensors(),
                                           quality=30.0,
                                           paraDX=0.3,
                                           paraMaxCellSize=0.5, paraDepth=8)

        mgr = ert.ERTManager(sr=True, useBert=True, verbose=False, debug=False)
        mgr.setMesh(mesh)
        mgr.setData(dat)
        pg.tic()
        mgr.fop.createJacobian(np.full(mgr.fop.parameterCount, 1))
        pg.toc()


if __name__ == '__main__':

    unittest.main()

