#!/usr/bin/env python
import unittest

import numpy as np
import pygimli as pg

from pygimli.physics import TravelTimeManager

class TestTT(unittest.TestCase):

    def setUp(self):
        # Dummy data container
        self.data = pg.DataContainer()
        self.data.createSensor([0.0, 0.0])
        self.data.createSensor([1.0, 2.0])
        self.data.resize(1)
        self.data.set("s", [1])
        self.data.set("g", [2])
        self.data.registerSensorIndex("s")
        self.data.registerSensorIndex("g")

        # Without secondary nodes
        self.mesh = pg.createGrid([0,1,2],[0,1,2])

        # Slowness
        self.slo = [1,2,1,4]
        
        self.mgr = TravelTimeManager()

    def test_withoutSecNodes(self):
        fop = self.mgr.fop
        fop.setData(self.data)
        fop.setMesh(self.mesh, 
                    ignoreRegionManager=True)
        t = fop.response(self.slo)
        np.testing.assert_allclose(t, 1 + np.sqrt(2))

        pg.show(self.mesh)
        data = self.mgr.simulate(slowness=self.slo, scheme=self.data, 
                                 mesh=self.mesh, secNodes=0)
        np.testing.assert_allclose(data['t'], 1 + np.sqrt(2))
        

    def test_withSecNodes(self):
        fop = self.mgr.fop
        fop.setData(self.data)
        fop.setMesh(self.mesh.createMeshWithSecondaryNodes(n=3), 
                    ignoreRegionManager=True)

        t = fop.response(self.slo)
        np.testing.assert_allclose(t, np.sqrt(5)) # only works for odd secNodes
        
        data = self.mgr.simulate(slowness=self.slo, scheme=self.data, 
                                 mesh=self.mesh, secNodes=3)
        np.testing.assert_allclose(data['t'], np.sqrt(5)) # only works for odd secNodes


    def test_Jacobian(self):
        fop = self.mgr.fop
        fop.setData(self.data)
        fop.setMesh(self.mesh.createMeshWithSecondaryNodes(n=5), 
                    ignoreRegionManager=True)

        fop.createJacobian(self.slo)
        J = fop.jacobian()
        np.testing.assert_allclose(J * self.slo, np.sqrt(5))

if __name__ == '__main__':

    # fop  = TestTT()
    # fop.test_MT()
    
    unittest.main()