#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""DOCUMENTME.

@author: Guenther.T
"""

import numpy as np

import pygimli as pg
from pygimli.utils import gmat2numpy


def iterateBounds(inv, dchi2=0.5, maxiter=100, change=1.02):
    """Find parameter bounds by iterating model parameter.

    Find parameter bounds by iterating model parameter until error
    bound is reached

    Parameters
    ----------
    inv :
        gimli inversion object
    dchi2 :
        allowed variation of chi^2 values [0.5]
    maxiter :
        maximum iteration number for parameter iteration [100]
    change:
        changing factor of parameters [1.02, i.e. 2%]
    """
    f = inv.forwardOperator()

    model = inv.model()
    resp = inv.response()

    nd, nm = len(resp), len(model)
    modelU = np.zeros(nm)
    modelL = np.zeros(nm)
    maxchi2 = inv.chi2() + dchi2

    for im in range(nm):
        model1 = pg.Vector(model)
        chi2 = .0
        it = 0

        while (chi2 < maxchi2) & (it < maxiter):
            it += 1
            model1[im] *= change
            resp1 = f(model1)
            chi2 = inv.getPhiD(resp1) / nd

        modelU[im] = model1[im]

        model2 = pg.Vector(model)
        chi2 = 0.0
        it = 0

        while (chi2 < maxchi2) & (it < maxiter):
            it += 1
            model2[im] /= change
            resp2 = f(model2)
            chi2 = inv.getPhiD(resp2) / nd

        modelL[im] = model2[im]

    return modelL, modelU

def scaledJacobianMatrix(inv):
    """Return error-weighted transformation-scaled Jacobian.

    Parameters
    ----------
    inv : pg.Inversion (pygimli.framework.Inversion)

    Returns
    -------
    DJ : numpy full matrix
    """
    J = inv.fop.jacobian()  # sensitivity matrix
    d = 1. / inv.dataTrans.error(inv.response, inv.errorVals)
    left = np.reshape(inv.dataTrans.deriv(inv.response) / d, [-1, 1])
    right = np.reshape(1 / inv.modelTrans.deriv(inv.model), [1, -1])
    if isinstance(J, pg.Matrix):  # e.g. ERT
        return left * pg.utils.gmat2numpy(J) * right
    elif isinstance(J, pg.SparseMapMatrix):  # e.g. Traveltime
        return left * pg.utils.sparseMat2Numpy.sparseMatrix2Dense(J) * right
    else:
        raise TypeError("Matrix type cannot be converted")



def modelResolutionMatrix(inv):
    """Formal model resolution matrix (MRM) from inversion.

    Parameters
    ----------
    inv : pg.Inversion (pygimli.framework.Inversion)

    Returns
    -------
    MR : pg.Matrix (pg.matrix.core.RMatrix dense matrix)
    """

    DJ = scaledJacobianMatrix(inv)
    JTJ = DJ.T.dot(DJ)


def modelCovariance(inv):
    """Formal model covariance matrix (MCM) from inversion.

    Parameters
    ----------
    inv : pygimli inversion object

    Returns
    -------
    var  : variances (inverse square roots of MCM matrix)
    MCMs : scaled MCM (such that diagonals are 1.0)

    Examples
    --------
    >>> # import pygimli as pg
    >>> # import matplotlib.pyplot as plt
    >>> # from matplotlib.cm import bwr
    >>> # INV = pg.Inversion(data, f)
    >>> # par = INV.run()
    >>> # var, MCM = modCovar(INV)
    >>> # i = plt.imshow(MCM, interpolation='nearest',
    >>> #                 cmap=bwr, vmin=-1, vmax=1)
    >>> # plt.colorbar(i)
    """

    DJ = scaledJacobianMatrix(inv)
    JTJ = DJ.T.dot(DJ)
    try:
        MCM = np.linalg.inv(JTJ)   # model covariance matrix

        varVG = np.sqrt(np.diag(MCM))  # standard deviations from main diagonal
        di = (1.0 / varVG)  # variances as column vector

        # scaled model covariance (=correlation) matrix
        MCMs = di.reshape(len(di), 1) * MCM * di
        return varVG, MCMs

    except BaseException as e:
        print(e)
        import traceback
        import sys

        traceback.print_exc(file=sys.stdout)
        return np.zeros(len(inv.model()),), 0

def modCovarCoreInv(inv):
    """Formal model covariance matrix (MCM) from inversion.

    var, MCMs = modCovar(inv)

    Parameters
    ----------
    inv : pygimli inversion object

    Returns
    -------
    var  : variances (inverse square roots of MCM matrix)
    MCMs : scaled MCM (such that diagonals are 1.0)

    Examples
    --------
    >>> # import pygimli as pg
    >>> # import matplotlib.pyplot as plt
    >>> # from matplotlib.cm import bwr
    >>> # INV = pg.Inversion(data, f)
    >>> # par = INV.run()
    >>> # var, MCM = modCovar(INV)
    >>> # i = plt.imshow(MCM, interpolation='nearest',
    >>> #                 cmap=bwr, vmin=-1, vmax=1)
    >>> # plt.colorbar(i)
    """
    td = np.asarray(inv.transData().deriv(inv.response()))
    tm = np.asarray(inv.transModel().deriv(inv.model()))

    J = td.reshape(len(td), 1) * \
        gmat2numpy(inv.forwardOperator().jacobian()) * (1. / tm)
    d = 1. / np.asarray(inv.transData().error(inv.response(), inv.error()))

    DJ = d.reshape(len(d), 1) * J
    JTJ = DJ.T.dot(DJ)
    try:
        MCM = np.linalg.inv(JTJ)   # model covariance matrix

        varVG = np.sqrt(np.diag(MCM))  # standard deviations from main diagonal
        di = (1.0 / varVG)  # variances as column vector

        # scaled model covariance (=correlation) matrix
        MCMs = di.reshape(len(di), 1) * MCM * di
        return varVG, MCMs

    except BaseException as e:
        print(e)
        import traceback
        import sys

        traceback.print_exc(file=sys.stdout)
        return np.zeros(len(inv.model()),), 0
