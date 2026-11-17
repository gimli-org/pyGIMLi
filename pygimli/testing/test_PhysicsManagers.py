#!/usr/bin/env python
"""Test the physics method managers (ERT, traveltime, VMD, VES)."""

# write a correct test!
import unittest
import numpy as np

import pygimli as pg
from pygimli.physics import ert
from pygimli.physics import VESManager
from pygimli.physics.em import VMDTimeDomainModelling

# pg.setTestingMode(True)
np.random.seed(1337)


class TestManagers(unittest.TestCase):
    """Test the ERT, traveltime, VMD and VES managers."""

    def test_ERT(self, showProgress=False):
        """Test ERT inversion with and without skin effect correction."""
        dat = pg.getExampleFile("ert/gallery.dat", load=True, verbose=True)
        dat["k"] = ert.createGeometricFactors(dat)

        mesh = pg.meshtools.createParaMesh(
            dat.sensors(),
            quality=30.0,
            paraDX=0.3,
            paraMaxCellSize=0.5,
            paraDepth=8,
        )

        # with SR
        mgr = ert.ERTManager(sr=True, useBert=True, verbose=False, debug=False)
        mgr.invert(dat, mesh=mesh, maxIter=20, lam=10)
        np.testing.assert_approx_equal(mgr.inv.chi2(), 1.040, significant=3)

        # without SR
        mgr = ert.ERTManager(sr=False, useBert=True, verbose=False)
        mgr.invert(dat, mesh=mesh, maxIter=20, lam=10)
        np.testing.assert_approx_equal(mgr.inv.chi2(), 1.033, significant=3)

    def test_TT(self, showProgress=False):
        """Keep a placeholder test for the traveltime manager."""
        pass

    @pg.skipOnDefaultTest
    def test_VMD(self, showProgress=False):
        """Test the VMD time domain modelling and inversion."""
        t = np.logspace(-5.5, -2.2, 20)
        verbose = False
        synthModel = np.array([25.0, 5.0, 100.0, 150.0, 1.0, 10.0, 4.0])
        nLayers = (len(synthModel) + 1) // 2
        fop = VMDTimeDomainModelling(
            times=t,
            txArea=10000.0,
            rxArea=10000.0,
            nLayers=nLayers,
            verbose=verbose,
        )
        # [thick[3], res[4]] nLay=4
        vmdMgr = pg.frameworks.MethodManager1d(fop)

        ra = vmdMgr.simulate(synthModel)

        err = abs(np.log(t) / 2) * 0.01
        ra *= 1.0 + pg.randn(len(ra), seed=1337) * err

        vmdMgr.invert(
            ra,
            err,
            nLayers=nLayers,
            layerLimits=[2, 500],
            maxIter=50,
            showProgress=showProgress,
            verbose=verbose,
        )

        if showProgress is True:
            fop.drawModel(ax=vmdMgr.inv.axs[0], model=synthModel, label="Synth")
        np.testing.assert_array_less(vmdMgr.fw.chi2(), 1.5)
        print("VMD chi2:", vmdMgr.fw.chi2())

    def test_VES(self, showProgress=False):
        """Test VES inversion for a layered halfspace."""
        thicks = [2.0, 10.0]
        res = [100.0, 5.0, 30]  # Ohm m
        phi = [0.0, 20.0, 0.0]  # neg mrad

        # model fails
        # thicks = [2., 6., 10.]
        # res = [100., 500., 20., 800.]
        # phi = [0., 20., 50., 0]

        synthModel = pg.cat(thicks, res)
        ab2 = np.logspace(np.log10(1.5), np.log10(100.0), 25)

        mgr = VESManager(verbose=False, debug=False)

        if showProgress:
            mgr.verbose = True
            fig, axs = pg.plt.subplots(2, 4, figsize=(12, 7))
            mgr.inv.axs = [axs[0][0], axs[1][0]]

        # Test -- basic
        ra, err = mgr.simulate(
            synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01, seed=0
        )
        mgr.exportData("synth.ves", ra, err)

        mgr.invert(
            ra,
            err,
            nLayers=4,
            lam=100,
            layerLimits=False,
            showProgress=showProgress,
        )
        if showProgress is True:
            mgr.fop.drawModel(ax=axs[0][0], model=synthModel, label="Synth")
        np.testing.assert_array_less(mgr.fw.chi2(), 1)

        # Test -- reinit with new parameter count
        if showProgress is True:
            mgr.inv.axs = [axs[0][1], axs[1][1]]
        mgr.invert(
            ra, err, nLayers=5, layerLimits=False, showProgress=showProgress
        )
        if showProgress is True:
            mgr.fop.drawModel(ax=axs[0][1], model=synthModel, label="Synth")
        # axs[0][1].legend()
        np.testing.assert_array_less(mgr.inv.inv.chi2(), 1)

        # Test -- reinit with new data basis
        ab2_2 = np.logspace(np.log10(1.5), np.log10(50.0), 10)
        ra, err = mgr.simulate(
            synthModel, ab2=ab2_2, mn2=1.0, noiseLevel=0.01, seed=0
        )

        if showProgress is True:
            mgr.inv.axs = [axs[0][2], axs[1][2]]

        mgr.invert(
            ra,
            err,
            nLayers=4,
            ab2=ab2_2,
            mn2=1.0,
            layerLimits=False,
            showProgress=showProgress,
        )
        if showProgress is True:
            mgr.fop.drawModel(ax=axs[0][2], model=synthModel, label="Synth")
        # axs[0][2].legend()
        # np.testing.assert_approx_equal(mgr.inv.inv.chi2(), 0.524220118768226,
        #                                significant=3)

        # Test -- reinit with complex resistivies
        mgr.complex = True
        synthModel = pg.cat(synthModel, np.array(phi) * 1e-3)

        ra, err = mgr.simulate(
            synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01, seed=1337
        )

        if showProgress is True:
            mgr.inv.axs = [axs[0][3], axs[1][3]]
        mgr.invert(
            ra, err, layerLimits=False, showProgress=showProgress, maxIter=50
        )

        np.testing.assert_array_less(mgr.inv.inv.chi2(), 1.09)

    def test_ERT_Simulate(self):
        """Test ERT simulation for a model with reversed K-factors."""
        import numpy as np
        import pygimli as pg
        import pygimli.meshtools as mt
        from pygimli.physics import ert

        world = mt.createWorld(
            start=[-50, 0], end=[50, -50], layers=[-1, -5], worldMarkers=True
        )

        scheme = pg.DataContainerERT()
        elec_positions = [
            [-2, 0, 0],
            [-1, 0, 0],
            [+0, 0, 0],
            [+1, 0, 0],
            [+2, 0, 0],
        ]
        for x, y, z in elec_positions:
            scheme.createSensor((x, y, z))

        # Define number of measurements
        scheme.resize(3)

        # first two measurements have reversed geometric factor!
        scheme.set("a", [0, 0, 3])
        scheme.set("b", [1, 1, 2])
        scheme.set("m", [3, 2, 1])
        scheme.set("n", [2, 3, 0])

        for pos in scheme.sensorPositions():
            world.createNode(pos)
            # world.createNode(pos + pg.RVector3(0, -0.1))

        mesh = mt.createMesh(world, quality=34)

        rhomap = [
            [1, 99.595 + 8.987j],
            [2, 99.595 + 8.987j],
            [3, 59.595 + 8.987j],
        ]

        # pg.show(mesh, markers=True, showMesh=True)
        data = ert.simulate(mesh, res=rhomap, scheme=scheme, verbose=True)
        rhoa = data["rhoa"]
        phia = data["phia"]

        # make sure all computed responses are equal, especially the first
        # two, which only differ in the sign of their geometrical factor
        np.testing.assert_allclose(rhoa, rhoa[0])
        np.testing.assert_allclose(phia, phia[0])

        # make sure the phases are positive (for positive input)
        self.assertTrue(np.all(phia > 0))

        # make sure rhoa is also positive
        self.assertTrue(np.all(rhoa > 0))


if __name__ == "__main__":
    unittest.main()
