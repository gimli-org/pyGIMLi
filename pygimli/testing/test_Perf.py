#!/usr/bin/env python
"""Test the performance of core vector and timing operations."""

# write a correct test!
import unittest
import timeit
import time

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
    """Test the performance of conversions, arithmetic and timing."""

    def test_Performance(self):
        """Print timings for vector conversions and arithmetic."""
        # pg.setDebug(True)

        sw = pg.Stopwatch(True)
        # print(timeit.repeat('r = grid.cellSizes() * np1', setup=setup,
        #                    number=1000))
        # print(timeit.repeat('r = c * np1', setup=setup, number=1000))

        print(
            ("np(np)", timeit.repeat("np.array(np1)", setup=setup, number=5000))
        )
        print(
            (
                "pg(pg)",
                timeit.repeat("pg.Vector(pg1)", setup=setup, number=5000),
            )
        )
        print(
            (
                "pg(np)",
                timeit.repeat("pg.Vector(np1)", setup=setup, number=5000),
            )
        )
        print(
            ("np(pg)", timeit.repeat("np.array(pg1)", setup=setup, number=5000))
        )
        print(
            (
                "np * np",
                timeit.repeat("np3 = np1 * np2", setup=setup, number=5000),
            )
        )
        print(
            (
                "pg * pg",
                timeit.repeat("pg3 = pg1 * pg2", setup=setup, number=5000),
            )
        )
        print(("pg * np", timeit.repeat("pg1 * np1", setup=setup, number=5000)))
        print(("np * pg", timeit.repeat("np1 * pg1", setup=setup, number=5000)))
        print(("sum(np)", timeit.repeat("sum(np1)", setup=setup, number=300)))
        print(("sum(pg)", timeit.repeat("sum(pg1)", setup=setup, number=300)))
        print(
            (
                "pg.sum(pg)",
                timeit.repeat("pg.sum(pg1)", setup=setup, number=300),
            )
        )
        print(
            ("pg.sum(st)", timeit.repeat("pg.sum(st)", setup=setup, number=300))
        )

        print(("s", sw.duration(True)))

        N = 10001
        np1 = np.linspace(1.1, 1.2, N)

        pg1 = pg.Vector(np1)

        # print(sw.duration(True))
        print(sum(np1 * np1))
        print(sum(pg1 * pg1))
        print(sum(np1 * pg1))
        print(sum(pg1 * np1))

    def test_tictoc1(self):
        """Test core tic toc."""
        pg.tic()
        time.sleep(0.1)
        print(pg.dur())

        time.sleep(0.1)
        pg.toc()
        print(pg.timings())

    def test_tictoc2(self):
        """Test core TicToc with global trace."""
        with pg.tictoc("foo1"):
            pass
        print(pg.timings())
        with pg.tictoc("foo2"):
            pass
        print(pg.timings())
        # return

        def bar1():
            with pg.tictoc("bar1"):
                time.sleep(0.01)

        def bar2():
            with pg.tictoc("bar2"):
                time.sleep(0.01)
                bar1()

        def foo():
            with pg.tictoc("foo"):
                for _i in range(4):
                    bar1()
                bar2()

        with pg.tictoc("foo"):
            with pg.tictoc("foobar1"):
                for _i in range(2):
                    time.sleep(0.01)  # uncovered
                    foo()

            with pg.tictoc("foobar2"):
                for _i in range(2):
                    time.sleep(0.01)  # uncovered
                    foo()
        print(pg.timings("foo"))

        with pg.tictoc("foo1"):
            foo()

        # print(pg.timings())
        print(pg.timings("foo"))
        print(pg.timings("foo/foobar1"))
        print(pg.timings("foo/foobar2"))
        print(pg.timings())


class TestMT(unittest.TestCase):
    """Test the traveltime and ERT Jacobian computation."""

    def test_WayMatrix(self):
        """MT-OMP Test.

        Run with env OMP_NUM_THREADS=4 GIMLI_NUM_THREADS=3 PG_USE_OMP=1

        OMP THREADS if USE_OMP=1 else DISTR_CALC with GIMLI_NUM_THREADS

        USE_OMP=0 and GIMLI_NUM_THREADS > 1 segfaults on py312
        """
        from pygimli.physics import traveltime as tt

        x, y = 25, 25
        mesh = pg.createGrid(x, y)
        data = tt.createRAData(
            [(0, 0)] + [(x, i) for i in range(y)], shotDistance=y + 1
        )
        data["t"] = 1.0
        mgr = tt.Manager(data, verbose=False)
        mgr.applyMesh(mesh, secNodes=10)
        pg.tic()
        mgr.fop.createJacobian(np.full(mgr.fop.parameterCount, 1))
        pg.toc()

    def test_ERTJacobian(self):
        """MT-OMP Test.

        Run with env OMP_NUM_THREADS=4 GIMLI_NUM_THREADS=3 PG_USE_OMP=1

        OMP THREADS if USE_OMP=1 else DISTR_CALC with GIMLI_NUM_THREADS

        USE_OMP=0 and GIMLI_NUM_THREADS > 1 segfaults on py312
        """
        from pygimli.physics import ert

        dat = pg.getExampleFile("ert/gallery.dat", load=True, verbose=True)
        dat["k"] = ert.createGeometricFactors(dat)
        mesh = pg.meshtools.createParaMesh(
            dat.sensors(),
            quality=30.0,
            paraDX=0.3,
            paraMaxCellSize=0.5,
            paraDepth=8,
        )

        mgr = ert.ERTManager(sr=True, useBert=True, verbose=False, debug=False)
        mgr.setMesh(mesh)
        mgr.setData(dat)
        pg.tic()
        mgr.fop.createJacobian(np.full(mgr.fop.parameterCount, 1))
        pg.toc()


if __name__ == "__main__":
    unittest.main()
