"""Magnetics forward operator."""
import numpy as np
import pygimli as pg
from .kernel import SolveGravMagHolstein


class MagneticsModelling(pg.frameworks.MeshModelling):
    """Magnetics modelling operator using Holstein (2007)."""

    def __init__(self, mesh=None, points=None, cmp=None, igrf=None,
                 verbose=False):
        """Set up forward operator.

        Parameters
        ----------
        mesh: pygimli:mesh
            Tetrahedral or hexahedral mesh.

        points: list|array of (x, y, z)
            Measuring points.

        cmp: [str,]
            Component of: gx, gy, gz, TFA, Bx, By, Bz, Bxy, Bxz, Byy, Byz, Bzz

        igrf: list|array of size 3 or 7
            International geomagnetic reference field.
            either:

            * [D, I, H, X, Y, Z, F] - declination,
                inclination, horizontal field, X/Y/Z components, total field OR
            * [X, Y, Z] - X/Y/Z
                components OR
            * [lat, lon] - latitude,
                longitude (automatic by pyIGRF)
        """
        # check if components do not contain g!
        super().__init__()
        if cmp is None:
            cmp = ["TFA"]
        if igrf is None:
            igrf = [0, 0, 50000]  # vertical inducing field

        self._refineH2 = False
        # self.createRefinedForwardMesh(refine=False, pRefine=False)
        self.mesh_ = mesh
        self.sensorPositions = points

        self.components = cmp
        self.igrf = None

        if hasattr(igrf, "__iter__"):
            if len(igrf) == 2: # lat lon
                pyIGRF = pg.optImport('pyIGRF', requiredFor="use igrf support"
                                      f" for {self.__class__.__name__}. "
                            "Please install pyIGRF with: pip install pyIGRF")

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

    def response(self, model):
        """Compute forward response.

        Arguments
        ---------
        model: array-like
            Model parameters.
        """
        if self.kernel is None:
            self.computeKernel()

        return self.J.dot(model)

    def createJacobian(self, model):
        """Do nothing as this is a linear problem.

        Abstract method to create the Jacobian matrix.
        Need to be implemented in derived classes.

        Arguments
        ---------
        model: array-like
            Model parameters.
        """
        pass


class RemanentMagneticsModelling(MagneticsModelling):
    """Remanent magnetics modelling operator for arbitrary magnetization."""

    def __init__(self, mesh, points, cmp=None, igrf=None):
        self.mesh_ = mesh
        self.mesh_["marker"] = 0
        if cmp is None:
            cmp = ["Bx", "By", "Bz"]
        if igrf is None:
            igrf = [0, 0, 50000]

        super().__init__(mesh=self.mesh_, points=points, igrf=igrf, cmp=cmp)
        self.magX = MagneticsModelling(
            self.mesh_, points, igrf=[1, 0, 0], cmp=cmp)
        self.magX.computeKernel()
        self.magY = MagneticsModelling(
            self.mesh_, points, igrf=[0, 1, 0], cmp=cmp)
        self.magY.computeKernel()
        self.magZ = MagneticsModelling(
            self.mesh_, points, igrf=[0, 0, 1], cmp=cmp)
        self.magZ.computeKernel()
        self.m1 = pg.Mesh(self._baseMesh)
        self.m2 = pg.Mesh(self._baseMesh)
        self.regionManager().addRegion(1, self.m1, 0)
        self.regionManager().addRegion(2, self.m2, 0)
        self.fak = 4e-7 * np.pi * 1e9  # H to B and T to nT
        self.JJ = pg.matrix.hstack([self.magX.jacobian(),
                                    self.magY.jacobian(),
                                    self.magZ.jacobian()])
        self.JJ.recalcMatrixSize()
        self.J = pg.matrix.ScaledMatrix(self.JJ, self.fak)
        self.setJacobian(self.J)

    def createJacobian(self, model):
        """Do nothing."""
        pass

    def response(self, model):
        """Add together all three responses."""
        modelXYZ = np.reshape(model, [3, -1])
        return (self.magX.response(modelXYZ[0]) +
                self.magY.response(modelXYZ[1]) +
                self.magZ.response(modelXYZ[2])) * self.fak
