"""Gravimetry modelling operator."""
import numpy as np
import pygimli as pg
from .kernel import SolveGravMagHolstein


class GravityModelling(pg.frameworks.MeshModelling):
    """Magnetics modelling operator using Holstein (2007)."""

    def __init__(self, mesh, points, cmp=None):
        """Set up forward operator.

        Parameters
        ----------
        mesh : pygimli:mesh
            tetrahedral or hexahedral mesh
        points : list|array of (x, y, z)
            measuring points
        cmp : list of str
            component of: gx, gy, gz
        """
        # check if components do not contain g!
        super().__init__(mesh=mesh)
        self.createRefinedForwardMesh(refine=False, pRefine=False)
        self.mesh_ = mesh
        self.sensorPositions = points
        if cmp is None:
            cmp = ["gz"]

        self.components = cmp
        # self.footprint = foot
        self.kernel = None
        self.J = pg.matrix.BlockMatrix()
        self.createKernel()

    def createKernel(self):
        """Create computational kernel."""
        self.kernel = SolveGravMagHolstein(self.mesh_,
                                           pnts=self.sensorPositions,
                                           cmp=self.components)
        self.J = pg.matrix.BlockMatrix()
        self.Ki = []
        self.Ji = []
        for iC in range(self.kernel.shape[1]):
            self.Ki.append(np.squeeze(self.kernel[:, iC, :]))
            self.Ji.append(pg.matrix.NumpyMatrix(self.Ki[-1]))
            self.J.addMatrix(self.Ji[-1], iC*self.kernel.shape[0], 0)

        self.J.recalcMatrixSize()
        self.setJacobian(self.J)

    def setMesh(self, mesh, ignoreRegionManager=False):
        """Set the mesh."""
        super().setMesh(mesh, ignoreRegionManager=False)
        self.createKernel(mesh)

    def response(self, model):
        """Compute forward response."""
        return self.J.dot(model)

    def createJacobian(self, model):
        """Do nothing as this is a linear problem."""
