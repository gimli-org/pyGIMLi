#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pygimli as pg
from .modelling import MeshModelling


class MultiFrameModelling(MeshModelling):
    """Full frame (multiple fop parallel) forward modelling."""

    def __init__(self, modellingOperator, scalef=1.0, **ini):
        """Init class and jacobian matrix."""
        super().__init__()
        self.ini = ini
        self.modellingOperator = modellingOperator
        self.jac = pg.matrix.BlockMatrix()
        self.scalef = scalef
        self.fops = []
        self.nf = 0
        self.nm = 0


    def setData(self, data, modellingOperator=None, **ini):
        """Distribute the data containers amongst the fops."""
        modellingOperator = modellingOperator or self.modellingOperator
        ini = self.ini or ini
        self.fops = []
        for datai in data:
            fopi = modellingOperator(**ini)
            fopi.setData(datai)
            # fopi.setData(pg.Vector(datai))
            self.fops.append(fopi)


    def setMeshPost(self, mesh):
        """Set mesh to all forward operators."""
        for fop in self.fops:
            fop.setMesh(mesh, ignoreRegionManager=True)

        self.prepareJacobian()



    def setDefaultBackground(self):
        """Set the default background behaviour."""
        super().setDefaultBackground()
        # regionIds = self.regionManager().regionIdxs()
        # pg.info("Found {} regions.".format(len(regionIds)))
        # if len(regionIds) > 1:
        #     bk = pg.sort(regionIds)[0]
        #     pg.info("Region with smallest marker ({0}) "
        #             "set to background".format(bk))
        #     self.setRegionProperties(bk, background=True)


    @property
    def parameterCount(self):
        """Return number of parameters."""
        pc = self.regionManager().parameterCount()
        if pc == 0:
            pc = self.fops[0].parameterCount

        return pc * len(self.fops)

    def prepareJacobian(self):
        """Build up Jacobian block matrix (once the sizes are known)."""
        self.jac.clear()
        self.nm = self.regionManager().parameterCount()
        if self.nm == 0:
            self.nm = self.fops[0].parameterCount

        self.nf = len(self.fops)
        print(self.nm, "model cells")
        nd = 0
        for i, fop in enumerate(self.fops):
            self.jac.addMatrix(fop.jacobian(), nd, i*self.nm)
            nd += fop.data.size()

        self.jac.recalcMatrixSize()
        self.setJacobian(self.jac)

    def createConstraints(self, C=None):
        """Create constraint matrix (special type for this)."""
        if C is not None:
            self.C1 = C
        elif isinstance(super().createConstraints(), pg.SparseMapMatrix):
            self.C1 = pg.SparseMapMatrix(self.constraintsRef())
            # make a copy because it will be overwritten
        else:
            self.C1 = self.constraints()

        self.C = pg.matrix.FrameConstraintMatrix(self.C1,
                                                 len(self.fops),
                                                 self.scalef)
        self.setConstraints(self.C)
        # cw = self.regionManager().constraintWeights()
        # self.regionManager().setConstraintsWeights(np.tile(cw, self.nf))
        # switch off automagic inside core.inversion which checks for
        # local modeltransform of the regionManager
        self.regionManager().setLocalTransFlag(False)

    def response(self, model):
        """Forward response."""
        mod = np.reshape(model, [len(self.fops), -1])
        return np.concatenate([fop.response(mo) for fop, mo in
                               zip(self.fops, mod)])

    def createJacobian(self, model):
        """Create Jacobian matrix."""
        mod = np.reshape(model, [len(self.fops), -1])
        for i, fop in enumerate(self.fops):
            fop.createJacobian(mod[i])

    def createDefaultStartModel(self):  # , dataVals):
        """Create standard starting model."""
        return pg.Vector(self.nm*self.nf, 10.0)  # look up in fop

    def createStartModel(self, dataVals):
        """Create a starting model from mean data values."""
        # return np.concatenate([f.createStartModel() for f in self.fops])
        self.nm = self.regionManager().parameterCount()
        return pg.Vector(self.nm*len(self.fops), np.median(dataVals))
