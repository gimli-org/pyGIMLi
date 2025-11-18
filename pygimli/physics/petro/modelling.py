#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Special meta forward operator for modelling with petrophysical relations."""

import pygimli as pg
from pygimli.frameworks import MethodManager


class PetroModelling(pg.Modelling):
    """
    Combine petrophysical relation m(p) with modelling class f(p).

    Combine petrophysical relation m(p) with modelling class f(p) to invert
    for m (or any inversion transformation) instead of p.
    """

    def __init__(self, fop, trans, mesh=None, verbose=False):
        """Save forward class and transformation, create Jacobian matrix."""
        pg.warn('do not use')
        super().__init__(verbose=verbose)
        self.fop = fop
        self.trans = trans  # class defining m(p)
        # self.setData(self.fop.data())

        if mesh is not None:
            self.setMesh(mesh)

    def setData(self, data):
        """TODO."""
        pg.core.ModellingBase.setData(self, data)
        self.fop.setData(data)

    def setMesh(self, mesh):
        """TODO."""
        if mesh is None and self.fop.mesh() is None:
            raise BaseException("Please provide a mesh for "
                                "this forward operator")

        if mesh is not None:
            self.fop.setMesh(mesh)
            self.fop.createRefinedForwardMesh(refine=False)

        # self.setMesh(f.mesh(), ignoreRegionManager=True) # not really nessary
        self.setRegionManager(self.fop.regionManagerRef())

        self.nModel = self.regionManager().parameterCount()

        self.jac = pg.MultRightMatrix(self.fop.jacobian())
        self.setJacobian(self.jac)

    def response(self, model):
        """Use inverse transformation to get p(m) and compute response."""
        tModel = self.trans(model)
        ret = self.fop.response(tModel)
        return ret

    def createJacobian(self, model):
        """Fill the individual jacobian matrices."""
        self.fop.createJacobian(self.trans(model))
        self.jac.r = self.trans.deriv(model)  # set inner derivative


class PetroJointModelling(pg.Modelling):
    """Cumulative (joint) forward operator for petrophysical inversions."""

    def __init__(self, f=None, p=None, mesh=None, verbose=True):
        """Constructor."""
        pg.warn('do not use')
        super().__init__(verbose=verbose)

        self.fops = None
        self.jac = None
        self.jacI = None
        self.mesh = None

        if f is not None and p is not None:
            self.setFopsAndTrans(f, p)

    def setFopsAndTrans(self, fops, trans):
        """TODO."""
        self.fops = [PetroModelling(fi, pi, self.mesh)
                     for fi, pi in zip(fops, trans)]

    def setMesh(self, mesh):
        """TODO."""
        self.mesh = mesh
        for fi in self.fops:
            fi.setMesh(mesh)

        self.setRegionManager(self.fops[0].regionManagerRef())
        self.initJacobian()

    def setData(self, data):
        """TODO."""
        for i, fi in enumerate(self.fops):
            fi.setData(data[i])

        self.initJacobian()

    def initJacobian(self):
        """TODO."""
        self.jac = pg.matrix.BlockMatrix()
        nData = 0
        for fi in self.fops:
            self.jac.addMatrix(fi.jacobian(), nData, 0)
            nData += fi.data().size()  # update total vector length
        self.setJacobian(self.jac)

    def response(self, model):
        """Create concatenated response for fop stack with model."""
        resp = []
        for f in self.fops:
            resp.extend(f.response(model))
        return resp

    def createJacobian(self, model):
        """Creating individual Jacobian matrices."""
        self.initJacobian()
        for f in self.fops:
            f.createJacobian(model)


class JointPetroInversion(MethodManager):  # bad name: no inversion framework!
    """TODO."""

    def __init__(self, managers, trans, verbose=False, debug=False, **kwargs):
        """TODO."""
        pg.warn('do not use')
        MethodManager.__init__(self, verbose=verbose, debug=debug, **kwargs)

        self.managers = managers
        self.trans = trans
        self.fops = []
        self.dataVals = pg.Vector(0)
        self.dataErrs = pg.Vector(0)
        self.mod = pg.Vector(0)  # resulting model
        self.data = None

        self.tD = pg.trans.TransCumulative()
        self.tM = managers[0].tM

        for mgr in self.managers:
            fop = mgr.createFOP(verbose)
            fop.setVerbose(verbose=verbose)
            self.fops.append(fop)

        self.fop.setFopsAndTrans(self.fops, self.trans)

    @staticmethod
    def createFOP(verbose=False):
        """Create forward operator."""
        fop = PetroJointModelling(verbose)

        return fop

    def createInv(self, fop, verbose=True, doSave=False):
        """TODO."""
        inv = pg.Inversion(verbose, doSave)
        inv.setForwardOperator(fop)

        return inv

    def model(self):
        return self.mod

    def setData(self, data):
        """TODO."""
        if isinstance(data, list):
            if len(data) == len(self.managers):
                self.tD.clear()
                self.dataVals.clear()
                self.dataErrs.clear()

                self.fop.setData(data)

                for i, mgr in enumerate(self.managers):
                    t = mgr.tD
                    self.tD.add(t, data[i].size())

                    self.dataVals = pg.cat(self.dataVals,
                                           data[i](mgr.dataToken()))

                    if mgr.errIsAbsolute:
                        self.dataErrs = pg.cat(
                            self.dataErrs,
                            data[i]('err') / data[i](mgr.dataToken()))
                    else:
                        self.dataErrs = pg.cat(self.dataErrs, data[i]('err'))

                self.data = data

                self.inv.setTransData(self.tD)
                self.inv.setTransModel(self.tM)
            else:
                raise BaseException("To few datacontainer given")

    def setMesh(self, mesh):
        """TODO."""
        self.fop.setMesh(mesh)

    def invert(self, data=None, mesh=None, lam=20, limits=None, **kwargs):
        """TODO."""

        if 'verbose' in kwargs:
            self.setVerbose(kwargs.pop('verbose'))

        self.setData(data)
        self.setMesh(mesh)

        nModel = self.fop.regionManager().parameterCount()
        startModel = None

        if limits is not None:
            if hasattr(self.tM, 'setLowerBound'):
                if self.verbose:
                    print('Lower limit set to', limits[0])
                self.tM.setLowerBound(limits[0])
            if hasattr(self.tM, 'setUpperBound'):
                if self.verbose:
                    print('Upper limit set to', limits[1])
                self.tM.setUpperBound(limits[1])

            startModel = pg.Vector(nModel, (limits[1]-limits[0])/2.0)
        else:
            for i in range(len(self.managers)):
                startModel += pg.Vector(nModel, pg.math.median(
                    self.trans[i].inv(
                        self.managers[i].createApparentData(self.data[i]))))
            startModel /= len(self.managers)

        self.inv.setModel(startModel)
        self.fop.regionManager().setZWeight(1.0)

        self.inv.setData(self.dataVals)
        self.inv.setRelativeError(self.dataErrs)
        self.inv.setLambda(lam)

        self.mod = self.inv.run()
        return self.mod

    def showModel(self, **showkwargs):
        """TODO."""
        if len(showkwargs):
            pg.show(self.fop.regionManager().paraDomain(),
                    self.mod, **showkwargs)


class PetroInversion(JointPetroInversion):
    """TODO."""

    def __init__(self, manager, trans, **kwargs):
        """TODO."""
        JointPetroInversion.__init__(self, [manager], [trans], **kwargs)

    def invert(self, data, **kwargs):
        """TODO."""
        return JointPetroInversion.invert(self, [data], **kwargs)


if __name__ == "__main__":
    pass
