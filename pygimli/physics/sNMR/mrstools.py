# -*- coding: utf-8 -*-
"""Tools for magnetic resonance sounding (MRS)."""

import numpy as np

#import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.viewer.mpl import draw1dmodel
from pygimli.utils import rndig


class MRS1dBlockQTModelling(pg.core.ModellingBase):
    """
    MRS1dBlockQTModelling - pygimli modelling class for block-mono QT inversion
    f=MRS1dBlockQTModelling(lay, kr, ki, zvec, t, verbose = False )
    """

    def __init__(self, nlay, kr, ki, zvec, t, verbose=False):
        """Initialize modelling class.

        Parameters
        ----------
        nlay : int [3]
            number of layers
        kr/ki : RMatrix
            real/imaginary kernel functions
        zvec : array
            vertical discretization of the kernel
        t : array
            time discretization
        """
        mesh = pg.meshtools.createMesh1DBlock(nlay, 2)  # thk, wc, T2*
        pg.core.ModellingBase.__init__(self, mesh, verbose)
        self.kr_ = kr
        self.ki_ = ki
        self.zv_ = np.array(zvec)
        self.nl_ = nlay
        self.nq_ = len(kr)
        self.t_ = np.array(t)
        self.nt_ = len(t)

    def response(self, par):
        """Compute response function."""
        nl = self.nl_
        thk = par(0, nl - 1)
        wc = par(nl - 1, 2 * nl - 1)
        t2 = par(2 * nl - 1, 3 * nl - 1)
        zthk = np.cumsum(thk)
        zv = self.zv_
        lzv = len(zv)
        izvec = np.zeros(nl + 1, np.int32)
        rzvec = np.zeros(nl + 1)
        for i in range(nl - 1):
            ii = (zv < zthk[i]).argmin()
            izvec[i + 1] = ii
            if ii <= len(zv):
                rzvec[i + 1] = (zthk[i] - zv[ii - 1]) / (zv[ii] - zv[ii - 1])
        izvec[-1] = lzv - 1
        A = np.zeros((self.nq_, self.nt_), dtype=complex)

        for i in range(nl):
            wcvec = np.zeros(lzv - 1)
            wcvec[izvec[i]:izvec[i + 1]] = wc[i]
            if izvec[i + 1] < lzv:
                wcvec[izvec[i + 1] - 1] = wc[i] * rzvec[i + 1]
            if izvec[i] > 0:
                wcvec[izvec[i] - 1] = wc[i] * (1 - rzvec[i])
            amps = np.dot(self.kr_, wcvec) + np.dot(self.ki_, wcvec) * 1j

            for ii in range(len(A)):
                A += np.exp(-self.t_ / t2[i]) * amps[ii]

        return np.abs(A).ravel()  # formerly pg as vector


def loadmrsproject(mydir):
    """
    load mrs project from given directory (zkernel.ve)

    (datacube.dat, KR/KI.bmat, zkernel.vec)
    """
    if mydir is None:
        mydir = '.'
    if mydir[-1] != '/':
        mydir = mydir + '/'
    # load files from directory
    zvec = pg.Vector(mydir + 'zkernel.vec')
    KR = pg.Matrix(mydir + 'KR.bmat')
    KI = pg.Matrix(mydir + 'KI.bmat')
    A = pg.Matrix()
    pg.loadMatrixCol(A, mydir + 'datacube.dat')
    t = np.array(A[0])
    # data
    datvec = pg.Vector()
    for i in range(1, len(A)):
        datvec = pg.cat(datvec, A[i])
    # print len(A), len(t)
    return KR, KI, zvec, t, datvec


def qtblockmodelling(mydir, nlay,
                     startvec=None, lowerbound=None, upperbound=None):
    """Loads data from dir, creates forward operator (old style)

    Parameters
    ----------
    mydir : string
        directory to load the files (kernel, data) from
    nlay : int
        number of layers
    startvec : array
        starting vector
    lowerbound : tuple/list of 3 values
        lower bound for thickness, water content and relaxation time
    upperbound : tuple/list of 3 values
        upper bound for thickness, water content and relaxation time

    Returns
    -------
    t : array
        time vector
    datvec : array
        data cube in a vector
    f : pygimli modelling class
        forward operator

    Examples
    --------

    t,datvec,f = qtblockmodelling(mydir,nlay,startvec=(10,0.3,0.2),
                                  lowerbound=(0.1,0,0.02),
                                  upperbound(100.,0.45,,0.5))
    """
    KR, KI, zvec, t, datvec = loadmrsproject(mydir)
    if startvec is None:
        startvec = [10., 0.3, 0.2]
    if lowerbound is None:
        lowerbound = [0.1, 0.0, 0.02]
    if upperbound is None:
        upperbound = [100., 0.45, 0.5]
    # modelling operator
    f = MRS1dBlockQTModelling(nlay, KR, KI, zvec, t)  # A[0])
    f.region(0).setStartValue(startvec[0])
    f.region(1).setStartValue(startvec[1])
    f.region(2).setStartValue(startvec[2])
    # Model transformations
    f.transTH = pg.trans.TransLogLU(lowerbound[0], upperbound[0])
    f.transWC = pg.trans.TransLogLU(lowerbound[1], upperbound[1])
    f.transT2 = pg.trans.TransLogLU(lowerbound[2], upperbound[2])
    f.region(0).setTransModel(f.transTH)
    f.region(1).setTransModel(f.transWC)
    f.region(2).setTransModel(f.transT2)
    return t, datvec, f


def showqtresultfit(thk, wc, t2, datvec, resp, t,
                    islog=True, clim=None, nu=3, nv=2):
    """Show mrs qt result and data fit.

    Parameters
    ----------
    thk : array of length nlay-1
        thickness vector
    wc : array of length nlay
        water content vector
    t2 : array of length nlay
        relaxation time vector
    datvec : array of length len(t)*len(Q)
        data vector
    t : array
        time vector
    islog : bool [True]
        use logarithms for colour scale
    clim : tuple of 2 floats
        color scale of data cube (in nV)
    nu/nv : int
        number of rows/columns for subplot

    Returns
    -------
    ax : mpl.Axes object

    Examples
    --------
    showqtresultfit(thk,wc,t2,datvec,resp,t,islog=True,clim=None,nu=3,nv=2)"""
    if clim is None:
        cma = max(datvec)
        cmi = min(datvec)
        if islog:
            cma = np.log10(cma)
            cmi = cma - 1.5
        clim = (cmi, cma)

    nt = len(t)
    nq = len(datvec) / nt
    si = (nq, nt)

    fig = plt.figure(1)
    ax1 = fig.add_subplot(nu, nv, 1)

    draw1dmodel(wc, thk, islog=False, xlab=r'$\theta$')
    ax3 = fig.add_subplot(nu, nv, 3)
    draw1dmodel(t2, thk, xlab='T2* in ms')
    ax3.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5])
    ax3.set_xticklabels(('0.02', '0.05', '0.1', '0.2', '0.5'))
    ax2 = fig.add_subplot(nu, nv, 2)
    if islog:
        plt.imshow(np.log10(np.array(datvec).reshape(si)),
                   interpolation='nearest', aspect='auto')
    else:
        plt.imshow(np.array(datvec).reshape(si),
                   interpolation='nearest', aspect='auto')
    plt.clim(clim)
    ax4 = fig.add_subplot(nu, nv, 4)
    if islog:
        plt.imshow(np.log10(resp.reshape(si)),
                   interpolation='nearest', aspect='auto')
    else:
        plt.imshow(resp.reshape(si),
                   interpolation='nearest', aspect='auto')
    misfit = np.array(datvec - resp)
    plt.clim(clim)
    ax5 = fig.add_subplot(nu, nv, 5)
    plt.hist(misfit, bins=30)
    plt.axis('tight')
    plt.grid(which='both')
    plt.text(plt.xlim()[0], np.mean(plt.ylim()),
             ' std=%g nV' % rndig(np.std(misfit), 3))
    ax6 = fig.add_subplot(nu, nv, 6)
    plt.imshow(misfit.reshape(si), interpolation='nearest', aspect='auto')
    ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    return ax
