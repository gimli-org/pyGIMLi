#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

    def test_ERT(self, showProgress=False):
        dat = pg.getExampleFile('ert/gallery.dat', load=True, verbose=True)
        dat['k'] = ert.createGeometricFactors(dat)
        mesh = pg.meshtools.createParaMesh(
            dat.sensors(), quality=30.0,
            paraDX=0.3, paraMaxCellSize=0.5, paraDepth=8)
        print(mesh)
        # with SR
        mgr = ert.ERTManager(sr=True, useBert=True, verbose=False, debug=False)
        mod = mgr.invert(dat, mesh=mesh, maxIter=20, lam=10)
        np.testing.assert_approx_equal(mgr.inv.chi2(), 1.040, significant=3)

        # without SR
        mgr = ert.ERTManager(sr=False, useBert=True, verbose=False)
        mod = mgr.invert(dat, mesh=mesh, maxIter=20, lam=10)
        np.testing.assert_approx_equal(mgr.inv.chi2(), 1.033, significant=3)

    def test_TT(self, showProgress=False):
        pass

    @pg.skipOnDefaultTest
    def test_VMD(self, showProgress=False):
        t = np.logspace(-5.5, -2.2, 20)
        verbose = False
        synthModel = np.array([25., 5., 100., 150., 1., 10., 4.])
        nLayers = (len(synthModel)+1) // 2
        fop = VMDTimeDomainModelling(times=t, txArea=10000.0, rxArea=10000.0,
                                     nLayers=nLayers, verbose=verbose)
        # [thick[3], res[4]] nLay=4

        vmdMgr = pg.frameworks.MethodManager1d(fop)

        ra = vmdMgr.simulate(synthModel)

        err = abs(np.log(t)/2) * 0.01
        ra *= 1. + pg.randn(len(ra), seed=1337) * err

        model = vmdMgr.invert(ra, err, nLayers=nLayers, layerLimits=[2, 500],
                              maxIter=50,
                              showProgress=showProgress, verbose=verbose)

        if showProgress is True:
            fop.drawModel(ax=vmdMgr.inv.axs[0],
                          model=synthModel, label='Synth')
        np.testing.assert_array_less(vmdMgr.fw.chi2(), 1.5)

    def test_VES(self, showProgress=False):
        """
        """
        thicks = [2., 10.]
        res = [100., 5., 30]  # Ohm m
        phi = [0., 20., 0.]  # neg mrad

        # model fails
        # thicks = [2., 6., 10.]
        # res = [100., 500., 20., 800.]
        # phi = [0., 20., 50., 0]

        synthModel = pg.cat(thicks, res)
        ab2 = np.logspace(np.log10(1.5), np.log10(100.), 25)

        mgr = VESManager(verbose=False, debug=False)

        if showProgress:
            mgr.verbose = True
            fig, axs = pg.plt.subplots(2, 4, figsize=(12, 7))
            mgr.inv.axs = [axs[0][0], axs[1][0]]

        # Test -- basic
        ra, err = mgr.simulate(synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01,
                               seed=0)
        mgr.exportData('synth.ves', ra, err)

        mgr.invert(ra, err, nLayers=4, lam=100, layerLimits=False,
                   showProgress=showProgress)
        if showProgress is True:
            mgr.fop.drawModel(ax=axs[0][0], model=synthModel, label='Synth')
        np.testing.assert_array_less(mgr.fw.chi2(), 1)

        # Test -- reinit with new parameter count
        if showProgress is True:
            mgr.inv.axs = [axs[0][1], axs[1][1]]
        mgr.invert(ra, err, nLayers=5, layerLimits=False,
                   showProgress=showProgress)
        if showProgress is True:
            mgr.fop.drawModel(ax=axs[0][1], model=synthModel, label='Synth')
        # axs[0][1].legend()
        np.testing.assert_array_less(mgr.inv.inv.chi2(), 1)

        # Test -- reinit with new data basis
        ab2_2 = np.logspace(np.log10(1.5), np.log10(50.), 10)
        ra, err = mgr.simulate(synthModel, ab2=ab2_2, mn2=1.0, noiseLevel=0.01,
                               seed=0)

        if showProgress is True:
            mgr.inv.axs = [axs[0][2], axs[1][2]]

        mgr.invert(ra, err, nLayers=4, ab2=ab2_2, mn2=1.0, layerLimits=False,
                   showProgress=showProgress)
        if showProgress is True:
            mgr.fop.drawModel(ax=axs[0][2], model=synthModel, label='Synth')
        # axs[0][2].legend()
        # np.testing.assert_approx_equal(mgr.inv.inv.chi2(), 0.524220118768226,
        #                                significant=3)

        # Test -- reinit with complex resistivies
        mgr.complex = True
        synthModel = pg.cat(synthModel, np.array(phi)*1e-3)

        ra, err = mgr.simulate(synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01,
                               seed=1337)

        if showProgress is True:
            mgr.inv.axs = [axs[0][3], axs[1][3]]
        mgr.invert(ra, err, layerLimits=False, showProgress=showProgress,
                   maxIter=50)

        np.testing.assert_array_less(mgr.inv.inv.chi2(), 1.09)


if __name__ == '__main__':
    unittest.main()
