#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Methods to calculate the model resolution.
"""

import numpy as np

import pygimli as pg


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
    d = inv.dataTrans.error(inv.response, inv.errorVals)
    left = np.reshape(inv.dataTrans.deriv(inv.response) / d, [-1, 1])
    right = np.reshape(1 / inv.modelTrans.deriv(inv.model), [1, -1])
    if isinstance(J, pg.Matrix):  # e.g. ERT
        return left * pg.utils.gmat2numpy(J) * right
    elif isinstance(J, pg.SparseMapMatrix):  # e.g. Traveltime
        return left * pg.utils.sparseMat2Numpy.sparseMatrix2Dense(J) * right
    else:
        raise TypeError("Matrix type cannot be converted")

def resolutionMatrix(inv, returnRD=False):
    r"""Formal model (and data) resolution matrix (MCM) from inversion.

    .. math::

        \mathbf{R_M} = (\mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J} + \alpha
        \mathbf{C}^T\mathbf{C})^{-1}
        \mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J}

        \mathbf{R_D} = \mathbf{D}\mathbf{J}
        (\mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J} + \alpha
        \mathbf{C}^T\mathbf{C})^{-1}
        \mathbf{J}^T\mathbf{D}^T

    Parameters
    ----------
    inv : pg.Inversion
        pygimli inversion instance after inversion
    returnRD : bool [false]
        also return data resolution (information) matrix

    Returns
    -------
    RM : numpy.array
        model resolution matrix
    RD : numpy.array
        data resolution matrix
    """
    DJ = scaledJacobianMatrix(inv)
    C = pg.utils.sparseMat2Numpy.sparseMatrix2Dense(inv.fop.constraints())
    cw = inv.fop.regionManager().constraintWeights()
    CC = np.reshape(cw ,[-1, 1]) * C
    JTJ = DJ.T @ DJ
    JI = np.linalg.inv(JTJ + CC.T @ CC * inv.lam)
    RM = JI @ JTJ
    if returnRD:
        RD = DJ @ JI @ DJ.T
        return RM, RD
    else:
        return RM

def modelResolutionMatrix(inv):
    r"""Formal model resolution matrix (MCM) from inversion.

    .. math::

        \mathbf{R_m} = (\mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J} + \alpha
        \mathbf{C}^T\mathbf{C})^{-1}
        \mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J}

    Parameters
    ----------
    inv : pg.Inversion
        pygimli inversion instance after inversion

    Returns
    -------
    RM : numpy.array
        model resolution matrix
    """
    return resolutionMatrix(inv)

def modelResolutionKernel(inv, nr=0, maxiter=50):
    """Compute single resolution kernel by solving an inverse problem.

    Parameters
    ----------
    inv : pg.Inversion
        inversion instance
    nr : int
        parameter/cell number
    maxiter : int
        maximum iterations for LSQR solver

    Returns
    -------
    reskernel : np.array
        resolution
    """
    from pygimli.solver.leastsquares import lsqr
    td = inv.dataTrans  # data transformation (e.g. lin/log/symlog)
    tm = inv.modelTrans  # model transformation (typically log or logLU)
    C = inv.fop.constraints()  # (sparse) regularization matrix
    left = td.deriv(inv.response) / inv.errorVals
    right = 1 / tm.deriv(inv.model)
    DS = pg.matrix.MultLeftRightMatrix(inv.fop.jacobian(), left, right)
    JC = pg.BlockMatrix()
    JC.addMatrix(DS, 0, 0)
    JC.addMatrix(C, DS.rows(), 0, np.sqrt(inv.lam))
    JC.recalcMatrixSize()
    if isinstance(nr, int):
        invec = pg.cat(pg.math.matrix.matrixColumn(DS, nr),
                       pg.Vector(C.rows()))
        return lsqr(JC, invec, maxiter=50)

def modelResolutionRadius(inv, nr=None, RM=None):
    """Compute resolution radius from model resolution matrix diagonal.

    According to Friedel (2003), it is defined as the radius of a circle (2D) or
    sphere (3D) having a resolution of 1, i.e. Sr = Se / RM(i, i) where Sr and
    are the sizes (area or volume) of the circle/sphere and the element.

    Parameters
    ----------
    nr : int|None [None]
        compute only resolution radius for a single cell (otherwise all cells)
    RM : numpy.matrix [None]
        already existing resolution matrix, otherwise compute it
    """
    pd = inv.paraDomain
    cs = pd.cellSizes()
    if nr is not None:  # only a single one
        if RM is not None:
            rm = RM[:, nr]
        else:
            rm = modelResolutionKernel(inv, nr=nr)[nr]

        cs = cs[nr]
    else:
        if RM is None:
            RM = modelResolutionMatrix(inv)

        rm = np.diag(RM)

    if pd.dim() == 2:
        return np.abs(cs/rm/np.pi)**0.5
    elif pd.dim() == 3:
        return np.abs(cs/rm*3/4/np.pi)**(1/3)
    else: # 1D
        return cs/rm

def computeR(J, C, alpha=0.5):
    r"""Return diagional of model resolution matrix.

    Calculates the formal model resolution matrix deterministically following:

    .. math::

        \mathbf{R_m} = (\mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J} + \alpha
        \mathbf{C}^T\mathbf{C})^{-1}
        \mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J}

    .. note::

        The current implementation assumes that :math:`\mathbf{D}` is the
        identitiy matrix, i.e. equal data weights.

    Parameters
    ----------
    J : array
        Jacobian matrix.
    C : array
        Constraint matrix.
    alpha : float
        Regularization strength :math:`\alpha`.
    """
    lin = pg.optImport("scipy.linalg", "Calculate model resolution matrices.")
    if lin:
        import scipy.sparse as sp
        from scipy.sparse import coo_matrix, csc_matrix
        from scipy.sparse.linalg import aslinearoperator, LinearOperator, lsqr

    JTJ = J.T.dot(J)
    CM_inv = C.T.dot(C)
#    Jsharp = lin.solve(JTJ + alpha * CM_inv, J.T)
#    R = np.diag(Jsharp.dot(J))
    RM = lin.solve(JTJ + alpha * CM_inv, JTJ)
    R = np.diag(RM)
    return R

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

#
# # self-made imports
# import pygimli as g
# from bert_tools import *
#
# #general imports
# import numpy as np
# from scipy.io import loadmat
# import scipy.linalg as lin
# import scipy.sparse as sp
# from scipy.sparse import coo_matrix, csc_matrix
# from scipy.sparse.linalg import aslinearoperator, LinearOperator, lsqr
# import matplotlib.pyplot as plt
#
# def load_constraint_matrix(fname):
#     """ Load constraint matrix in sparse format """
#     global mesh
#     i, j, data = np.loadtxt(fname, unpack=True, dtype=int)
#     C = coo_matrix((data, (i,j)), shape=(mesh.boundaryCount(), mesh.cellCount()), dtype=int)
#     return C
#
# def Rmult(x):
#     """ Performs matrix-vector multiplication y = R * x """
#     global J
#     global C
#     global linop
#     r, s = C.shape
#     Jx = J.dot(x)
#     Jxpad = np.hstack((Jx, np.zeros(r)))
#     y = lsqr(linop, Jxpad, atol=1.0e-4, iter_lim=5000, show=False)
#     return y[0]
#
# def diagestrand(dim):
#     """
#         Diagonal estimation based on the algorithm of Bekas et al. (2007)
#     """
#     global Rmult
#
#     chi = 0.01  # stopping criterion (based on convergence tests)
#     d = np.zeros(dim)
#     t = d.copy()
#     q = d.copy()
#     norm = 999
#     k = 0
#
#     while norm > chi:
#         # draw random vector from normal distribution
#         v = np.random.randn(dim)
#         t += Rmult(v) * v
#         q += v * v
#         d_k = t / q
#
#         # evalute difference norm to check convergence
#         norm = np.linalg.norm(d - d_k)
#         d = d_k
#         k +=1
#
#     print "Stopped stochastic estimation after %d iterations with chi = %s." % (k, chi)
#     return d
#
# def compute_res(J, C, alpha=0.5):
#     """
#     Determinisctic computation of model resolution """
#
#     print "Starting deterministic estimation of resolution diagonal..."
#
#     JTJ = J.T.dot(J)
#     CM_inv = C.T.dot(C)
#     Jsharp = lin.solve(JTJ + alpha * CM_inv, J.T)
#     R = np.diag(Jsharp.dot(J))
#     return R
#
# def compare_res(J, C, alpha=0.5):
#     """ Compare determinstic to stochastic resolution diagonal estimate """
#
#     R_det = compute_res(J, C, alpha)
#     R_est = diagestrand(J.shape[1])
#
#     fig = plt.figure(figsize=(12,5))
#     ax1 = fig.add_subplot(121)
#     ax1.plot(R_det, R_est, 'k.', alpha=0.5)
#     ax1.plot((0,1), 'k--')
#     plt.xlabel('Deterministic estimate of resolution diagonal')
#     plt.ylabel('Stochastic estimate of resolution diagonal')
#     plt.axis('scaled')
#     plt.xlim(0, R_det.max()+0.1)
#     plt.ylim(0, R_det.max()+0.1)
#
#     ax2 = fig.add_subplot(122)
#     ax2.hist(R_det - R_est, 20, color='gray')
#     plt.xlabel('Deterministic - stochastic estimate of resolution diagonal')
#     plt.xlim(-0.15, 0.15)
#     plt.ylabel('Count')
#     plt.draw()
#     plt.show()
#
#     return fig
#
# if __name__ == '__main__':
#     alpha = 0.5
#     mesh = g.Mesh('meshParaDomain.bms')
#     J = loadsens('sens.bmat')
#     C = load_constraint_matrix('constraint.matrix')
#
#     p, q = J.shape
#     r, s = C.shape
#     linop = aslinearoperator(sp.vstack((J, C), 'csc'))
#
#     fig = compare_res(J,C)
#
#
#
# #diagR = diagestrand(q, 256)
# #
# #ntimes = 3
#
# #diags=np.zeros((ntimes, q))
#
# #for k in range(ntimes):
# ## use 256 random vectors for each estimate
# #numcols=256;
# #print 'iteration', k
#
# ## estimate the diagonal of r.
# #diagr=diagestrand(q, numcols);
#
# ## store the estimate.
# #diags[k,:] = diagr
#
# ## take median value.
# #diagrest = np.median(diags, 0)
