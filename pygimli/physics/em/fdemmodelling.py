#!/usr/bin/env python
"""Frequency Domain Electromagnetics (FDEM) functions and class."""

from pathlib import Path
import numpy as np

import pygimli as pg
from .hemmodelling import HEMmodelling

class FDEM2dFOPold(pg.core.ModellingBase):
    """Old variant of 2D FOP (to be deleted)."""

    def __init__(self, data, nlay=2, verbose=False):
        """Initialize with data and number of layers."""
        pg.core.ModellingBase.__init__(self, verbose)
        self.nlay = nlay
        self.FOP1d = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        self.mesh_ = pg.meshtools.createMesh1D(self.nx, 2 * nlay - 1)
        self.setMesh(self.mesh_)

    def response(self, model):
        """Yield forward model response."""
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        resp = pg.Vector(0)
        for modi in modA:
            resp = pg.cat(resp, self.FOP1d.response(modi))

        return resp


class FDEM2dFOP(pg.core.ModellingBase):
    """FDEM 2d-LCI modelling class based on BlockMatrices."""

    def __init__(self, data, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers."""
        super().__init__(verbose)
        self.nlay = nlay
        self.header = {}
        self.pos, self.z, self.topo = None, None, None
        self.FOP = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        npar = 2 * nlay - 1
        self.mesh1d = pg.meshtools.createMesh1D(self.nx, npar)
        self.mesh_ = pg.meshtools.createMesh1D(self.nx, 2 * nlay - 1)
        self.setMesh(self.mesh_)

        # self.J = NDMatrix(self.nx, self.nf*2, npar)
        self.J = pg.matrix.BlockMatrix()
        self.FOP1d = []
        for i in range(self.nx):
            self.FOP1d.append(pg.core.FDEM1dModelling(
                nlay, data.freq(), data.coilSpacing, -data.height))
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, self.nf * 2 * i, npar * i)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())

    def response(self, model):
        """Cut together forward responses of all soundings."""
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        resp = pg.Vector(0)
        for modi in modA:
            resp = pg.cat(resp, self.FOP.response(modi))

        return resp

    def createJacobian(self, model):
        """Create Jacobian matrix by creating individual Jacobians."""
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class HEM1dWithElevation(pg.core.ModellingBase):
    """Airborne FDEM modelling including variable bird height."""

    def __init__(self, frequencies, coilspacing, nlay=2, verbose=False):
        """Set up class by frequencies and geometries."""
        pg.core.ModellingBase.__init__(self, verbose)
        self.nlay_ = nlay  # real layers (actually one more!)
        self.FOP_ = pg.core.FDEM1dModelling(nlay + 1, frequencies,
                                            coilspacing, self.height)
        self.mesh_ = pg.meshtools.createMesh1D(nlay, 2)
        # thicknesses & res
        self.mesh_.cell(0).setMarker(2)
        self.setMesh(self.mesh_)

    def response(self, model):
        """Return forward response for a given model."""
        thk = model(0, self.nlay)  # thicknesses including height
        res = model(self.nlay - 1, self.nlay * 2)
        res[0] = 10000.
        return self.FOP_.response(pg.cat(thk, res))
