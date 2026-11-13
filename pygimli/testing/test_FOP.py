#!/usr/bin/env python
"""Basic tests."""

import unittest

import pygimli as pg
import pygimli.meshtools as mt

run_solve = True


class ModellingMT(pg.core.ModellingBase):
    """Modelling class used by the FOP tests."""

    def __init__(self, nPars, verbose):
        """Set up the model with `nPars` parameter values."""
        pg.core.ModellingBase.__init__(self, verbose)
        self.regionManager().setParameterCount(nPars)

    def dummySolve(self, info):
        """Solve a dummy problem and print its info."""
        world = mt.createWorld(
            start=[-10, 0], end=[10, -10], marker=1, worldMarkers=False
        )
        c1 = mt.createCircle(pos=[0.0, -5.0], radius=3.0, area=0.1, marker=2)
        mesh = pg.meshtools.createMesh([world, c1], quality=34.3)

        u = pg.solver.solveFiniteElements(
            mesh, a=[[1, 100], [2, 1]], bc={"Dirichlet": [[4, 1.0], [2, 0.0]]}
        )
        print(info, mesh, u)

    def response(self, par):
        """Return the model response for `par`."""
        if run_solve:
            self.dummySolve("ST")
        return par * 1.0

    def response_mt(self, par, i=0):
        """Return the MT response of component `i` for `par`."""
        if run_solve:
            self.dummySolve("MT" + str(i))
            # pg.solver.solve(pg.createGrid(range(100), range(100)), a=1, f=1,
            # bc={'Dirichlet': 0})

        # print(i)
        return par * 2.0


class TestFOP(unittest.TestCase):
    """Tests for the forward operator base classes."""

    def test_FOP(self):
        """Test FOP."""
        # ab2 = pg.Vector(2, 1.0)
        # ab2[1] = 2.0
        # mn2 = pg.Vector(2, 3.0)
        # mn2[1] = 4.0

        # nlay = 2
        # model = pg.Vector(3, 10.)

        # F = pg.core.DC1dModelling(nlay, ab2, mn2)

        # print(F.response(model))

        pass

    def test_multiResponseMT(self):
        """Test FOP response - mt."""
        ####### temporary deactivated -- test me
        # m = pg.Vector(nPars, 1)
        # fop = ModellingMT(nPars, verbose=True)

        # ms = np.array([m*2, m*3, m*4, m*5])
        # print(ms)
        # fop.setMultiThreadJacobian(1)
        # print('#')

        # ds1 = np.zeros((len(ms), len(ms[0])))
        # #print('#')
        # fop.responses(ms, ds1)

        # print('#')

        # fop.setMultiThreadJacobian(nPars)
        # ds2 = np.zeros((len(ms), len(ms[0])))
        # fop.responses(ms, ds2)

        # np.testing.assert_array_equal(ds1, ds2)
        ####### temporary deactivated -- test me

    def test_MT(self):
        """Test the MT response."""
        ####### temporary deactivated -- test me
        # m = pg.Vector(nPars, 1)

        # fop = ModellingMT(nPars, verbose=False)
        # fop.setMultiThreadJacobian(1)
        # fop.setThreadCount(1)
        # fop.createJacobian(m)
        # J1 = pg.Matrix(fop.jacobian())

        # fop.setMultiThreadJacobian(nPars)
        # fop.setThreadCount(1)
        # fop.createJacobian(m)
        # J2 = fop.jacobian()

        # np.testing.assert_array_equal(J1 * 2.0, J2)
        #######  temporary deactivated  -- test me


if __name__ == "__main__":
    fop = TestFOP()
    fop.test_MT()
    # fop.test_multiResponseMT()

    # unittest.main()
