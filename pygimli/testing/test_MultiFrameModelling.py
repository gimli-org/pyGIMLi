#!/usr/bin/env python
"""
Regression test for a broken 4D (timelapse) inversion.

In MultiFrameModelling the per-frame Jacobians are still empty when
prepareJacobian() builds the block matrix. As recalcMatrixSize() was only
called there, the block matrix kept the row/column offsets of the last
frame instead of its full size, i.e. for nFrames=4 it reported 1395x2124
(3 frames) instead of 1860x2832 (4 frames). The solver then worked with
"J.rows != nData 1395 / 1860", the line search degenerated and the
inversion diverged into NaNs.
"""

import unittest

import numpy as np

import pygimli.meshtools as mt
from pygimli.frameworks import MultiFrameModelling
from pygimli.physics import ert


class TestMultiFrameModelling(unittest.TestCase):
    """Test the block Jacobian of the full-frame forward operator."""

    def setUp(self):
        """Create a small mesh, one data set and two frames."""
        self.nFrames = 2
        mesh = mt.createMesh(
            mt.createWorld(start=[-6, 0], end=[12, -6], area=1.0), quality=34
        )
        self.mesh = mesh
        scheme = ert.createData(
            elecs=9, spacing=1, schemeName="dd", maxSeparation=4
        )
        self.data = ert.simulate(
            mesh=mesh,
            res=100.0,
            scheme=scheme,
            noiseLevel=0.0,
            noiseAbs=0.0,
            seed=42,
            verbose=False,
        )
        self.fop = MultiFrameModelling(ert.ERTModelling)
        self.fop.setData([self.data] * self.nFrames)
        self.fop.setMesh(mesh)
        # like TimelapseERT.fullInversion: calling mesh() creates the forward
        # mesh and hands it over to the frame operators
        self.fop.mesh()

    def test_JacobianCoversAllFrames(self):
        """Jacobian must have the size of all frames, not of all but one."""
        dataVec = np.concatenate([self.data["rhoa"]] * self.nFrames)
        model = self.fop.createStartModel(dataVec)
        self.fop.createJacobian(model)

        self.assertEqual(model.size(), self.nFrames * self.mesh.cellCount())
        self.assertEqual(
            self.fop.jacobian().rows(), self.nFrames * self.data.size()
        )
        self.assertEqual(self.fop.jacobian().cols(), model.size())


if __name__ == "__main__":
    unittest.main()
