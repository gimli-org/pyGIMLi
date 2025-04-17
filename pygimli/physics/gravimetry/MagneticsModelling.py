"""Magnetics forward operator."""
import numpy as np
import pygimli as pg
from .kernel import SolveGravMagHolstein


class MagneticsModelling(pg.frameworks.MeshModelling):
    """Magnetics modelling operator using Holstein (2007)."""

    def __init__(self, mesh=None, points=None, cmp=["TFA"], igrf=[50, 13]):
        """Setup forward operator.

        Parameters
        ----------
        mesh : pygimli:mesh
            tetrahedral or hexahedral mesh
        points : list|array of (x, y, z)
            measuring points
        cmp : list of str
            component of: gx, gy, gz, TFA, Bx, By, Bz, Bxy, Bxz, Byy, Byz, Bzz
        igrf : list|array of size 3 or 7
            international geomagnetic reference field, either
            [D, I, H, X, Y, Z, F] - declination, inclination, horizontal field,
                                    X/Y/Z components, total field OR
            [X, Y, Z] - X/Y/Z components OR
            [lat, lon] - latitude, longitude (automatic by pyIGRF)
        """
        # check if components do not contain g!
        super().__init__()
        self._refineH2 = False
        # self.createRefinedForwardMesh(refine=False, pRefine=False)
        self.mesh_ = mesh
        self.sensorPositions = points

        self.components = cmp
        self.igrf = None
        if hasattr(igrf, "__iter__"):
            if len(igrf) == 2: # lat lon
                import pyIGRF
                self.igrf = pyIGRF.igrf_value(*igrf)
            else: # 3 (x,y,z) or 7 (T,H,X,Y,Z,D,I)
                self.igrf = igrf

        self.kernel = None
        self.J = pg.matrix.BlockMatrix()
        if self.mesh_ is not None:
            self.setMesh(self.mesh_)

    def computeKernel(self):
        """Compute the kernel."""
        points = np.column_stack([self.sensorPositions[:, 1],
                                  self.sensorPositions[:, 0],
                                  -np.abs(self.sensorPositions[:, 2])])
        self.kernel = SolveGravMagHolstein(self.mesh().NED(),
                                           pnts=points, igrf=self.igrf,
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

    # better move the latter to
    # self.createKernel

    # def setMesh(self, mesh):
        # self.createKernel(mesh)

    def response(self, model):
        """Compute forward response."""
        if self.kernel is None:
            self.computeKernel()

        return self.J.dot(model)

    def createJacobian(self, model):
        """Do nothing as this is a linear problem."""
        pass
