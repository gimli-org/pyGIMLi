#!/usr/bin/env python
"""
"""
import unittest

import numpy as np
import pygimli as pg
from pygimli.frameworks import MeshModelling

class DummyMod(MeshModelling):
    def __init__(self, mesh, verbose=True):
        super(DummyMod, self).__init__()
        self.meshlist = []
        for i in range(2):
            for cell in mesh.cells():
                cell.setMarker(i + 1)
            self.meshlist.append(pg.Mesh(mesh))
            self.regionManager().setMesh(self.meshlist[i])
            self.regionManager().addRegion(i + 1, self.meshlist[i])
            self.regionManager().region(i + 1).setConstraintType(1)

        # self.setMesh(self.meshlist[0], ignoreRegionManager=True)


class TestRM(unittest.TestCase):

    def test_Constraints(self):
        """ Test FOP """

        grid = pg.createGrid(x=[0., 1., 2.], y=[0., 1., 2.])
        fop = DummyMod(grid)

        fop.createConstraints()
        # pg.solver.showSparseMatrix(fop.constraints(), full=True)

        fop.regionManager().setConstraintType(2)
        fop.createConstraints()
        # pg.solver.showSparseMatrix(fop.constraints(), full=True)

    def test_zweight(self):
        # what's this test for?
        mesh = pg.meshtools.createGrid(x=np.arange(0, 4), y=np.arange(0, 3))

        # fop = pg.Modelling()

        # fop.setMesh(mesh)
        # fop.createConstraints()

        # rm = fop.regionManager()
        # rm.addRegion(1, mesh, 0)
        # rm.addRegion(2, mesh, 0)

        # rm.setZWeight(0.1)

        # # check distribution of zWeight
        # self.assertTrue(np.isclose(rm.region(0).zWeight(), 0.1))
        # self.assertTrue(np.isclose(rm.region(1).zWeight(), 0.1))
        # self.assertTrue(np.isclose(rm.region(2).zWeight(), 0.1))

        # w0 = rm.region(0).constraintWeights()
        # w1 = rm.region(1).constraintWeights()
        # w2 = rm.region(2).constraintWeights()

        # print(w0.array())
        # print(w1.array())
        # print(w2.array())

        # # check actual constraint weight values
        # self.assertTrue(np.isclose(np.min(w0), 0.1))

        # # check distribution of zWeight
        # self.assertTrue(np.allclose(w0, w1))
        # self.assertTrue(np.allclose(w0, w2))

    def test_zweight_2meshes(self):

        # marker = 0
        mesh = pg.meshtools.createGrid(x=np.arange(0, 4), y=np.arange(0, 3))
        mesh.setCellMarkers(np.zeros(mesh.cellCount(), dtype=int))

        # marker = 1
        mesh2 = pg.Mesh(mesh)
        mesh2.setCellMarkers(np.ones(mesh.cellCount(), dtype=int))

        together = pg.meshtools.merge2Meshes(
            mesh,
            mesh2.translate(pg.Vector([0., 0., 1.])))

        fop = pg.Modelling()

        fop.setMesh(together)
        fop.createConstraints()

        rm = fop.regionManager()

        rm.setZWeight(0.1)

        # check distribution of zWeight
        self.assertTrue(np.isclose(rm.region(0).zWeight(), 0.1))
        self.assertTrue(np.isclose(rm.region(1).zWeight(), 0.1))
        # self.assertTrue(np.isclose(rm.region(2).zWeight(), 0.1))

        w0 = rm.region(0).constraintWeights()
        w1 = rm.region(1).constraintWeights()

        # print(w0, w1)
        # w2 = rm.region(2).constraintWeights()

        # check actual constraint weight values
        self.assertTrue(np.isclose(np.min(w0), 0.1))

        # check distribution of zWeight
        self.assertTrue(np.allclose(w0, w1))
        # self.assertTrue(np.allclose(w0, w2))


if __name__ == '__main__':

    unittest.main()
