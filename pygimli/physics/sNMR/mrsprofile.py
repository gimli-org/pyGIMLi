#!/usr/bin/env python
"""Classes for inverting profile data with magnetic resonance soundings (MRS).

The preferred LCI type of MRS inversion was published in:
Costabel, S., Günther, T., Dlugosch, R. & Müller-Petke, M. (2016):
Torus-nuclear magnetic resonance: Quasi-continuous airborne magnetic resonance
profiling by using a helium-filled balloon. Geophysics 81(4), W119-W129,
doi:10.1190/geo2015-0467.1.
"""

import sys
import os
from glob import glob
import numpy as np

import pygimli as pg
from pygimli.viewer.mpl import showStitchedModels
from . mrs import MRS
from . modelling import MRS1dBlockQTModelling


class MultiFOP(pg.core.ModellingBase):  # classical joint FOP => frameworks
    def __init__(self, mrsAll, nlay=3):
        mesh = pg.meshtools.createMesh1DBlock(nlay, 2)  # thk, wc, T2*
        pg.core.ModellingBase.__init__(self, mesh)
        self.fAll = []
        for mrs in mrsAll:
            mrs.createFOP(nlay)
            self.fAll.append(mrs.f)

    def response(self, model):
        """compute response by concatenating individual responses."""
        resp = pg.Vector()
        for f in self.fAll:
            resp = pg.cat(resp, f.response(model))
        return resp


class JointMRSModelling(MRS1dBlockQTModelling):
    """MRS Laterally constrained modelling based on BlockMatrices."""

    def __init__(self, mrs, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers"""
        super(JointMRSModelling, self).__init__(nlay, mrs[0].K, mrs[0].z,
                                                mrs[0].t, verbose)
        self.mrs = mrs
        for mrsi in mrs:
            mrsi.createFOP(nlay)

    def response(self, model):
        """response function (all responses together)"""
        response = pg.Vector()
        for mrs in self.mrs:
            response = pg.cat(response, mrs.f(model))


class MRSLCI(pg.core.ModellingBase):
    """MRS Laterally constrained modelling based on BlockMatrices."""

    def __init__(self, profile, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers"""
        super(MRSLCI, self).__init__(verbose)
        self.nlay = nlay
        self.nx = len(profile)
        self.np = 3 * nlay - 1
#        self.mesh2d = pg.meshtools.createMesh2D(npar, self.nx)
#        self.mesh2d = pg.meshtools.createMesh2D(self.nx, self.np)
        self.mesh2d = pg.meshtools.createMesh2D(range(self.np+1), range(self.nx+1))
        self.mesh2d.rotate(pg.RVector3(0, 0, -np.pi/2))
        self.setMesh(self.mesh2d)
        self.J = pg.matrix.BlockMatrix()
        self.FOP1d = []
        ipos = 0
        for i, mrs in enumerate(profile):
            mrs.createFOP(nlay, verbose=False)
            self.FOP1d.append(mrs.f)
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, ipos, self.np * i)
            ipos += len(mrs.data)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())
        self.setJacobian(self.J)

    def response(self, model):
        """cut-together forward responses of all soundings"""
        modA = np.asarray(model).reshape(self.nx, self.np)
#        modA = np.reshape(model, (self.nx, self.np))
        resp = pg.Vector(0)
        for i, modi in enumerate(modA):
            resp = pg.cat(resp, self.FOP1d[i].response(modi))

        return resp

    def createJacobian(self, model):
        """ fill individual blocks of Block-Jacobian matrix """
        modA = np.asarray(model).reshape(self.nx, self.np)
#        modA = np.reshape(model, (self.nx, self.np))
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class MRSprofile():
    """Manager class for several MRS data along a profile for joint inversion.

    Attributes
    ----------
    mrs: list of MRS objects (single soundings)
    x: list of positions for the soundings


    Methods
    -------
    load - load mrs files from a directory

    set X - set x vector

    showData - show MRS data

    independentBlock1dInversion - perform independent 1D block inversion

    block1dInversion - 1D block inversion of all data sets together

    blockLCInversion - 1D block laterally constrained inversion of all data

    printFits - print total misfit (chi^2, rms) and individual values

    showModel - show LCI model
    """

    def __init__(self, filename=None, x=None, dx=1, x0=0, **kwargs):
        """Initialize profile object by mrs objects and optional positions.

        Arguments
        ---------
        filename : list of str | str
            list of files OR filenames(with wildcard '*') OR directory to load
        x : iterable
            position vector of individual soundings
        x0 : float [0]
            starting position
        dx : position [1]
            position increment
        """
        self.mrs = []
        self.nData = 0
        self.figs = {}
        self.totalChi2 = None
        self.totalRMS = None
        self.WMOD = None
        self.TMOD = None
        self.RMSvec = None
        self.Chi2vec = None
        if '*' in filename:  # a filename with asterisks
            files = glob(filename)
        elif os.path.isdir(filename):  # a directory with all files to take
            files = glob(filename+'/*.mrsi')
            if len(files) == 0:
                files = glob(filename+'/*.mrsd')
        elif hasattr(filename, '__iter__'):  # already a list
            files = filename
        else:
            if filename is not None:
                print('Do not know what to do with filename')
            return
        self.load(files, **kwargs)
        if x is not None:
            self.setX(x)
        else:
            self.x = np.arange(len(self.mrs)) * dx + x0

    def __repr__(self):
        return "MRS profile with "+str(len(self.mrs))+" soundings"

    def setX(self, x=None, x0=0, dx=1):
        """define positions for soundings and sort accordingly"""
        if x is None:
            x = np.arange(len(self.mrs)) * dx + x0
            print(x)
        ind = np.argsort(x)
        self.mrs = self.mrs(ind)
        self.x = np.sort(x)

    def load(self, filenames, **kwargs):
        """ load mrs files in a list of (single) MRS handlers
            filename can be a list of mrsi files or a directory to search
            Additional parameters: usereal, mint, maxt (see MRS.load)
        """
        self.mrs = [MRS(filename, **kwargs) for filename in filenames]

    def loadKernel(self, kernelfile):
        """ load one kernel file for all soundings """
        self.mrs[0].loadKernel(kernelfile)
        for mrsi in self.mrs:
            mrsi.z = self.mrs[0].z
            mrsi.K = self.mrs[0].K

    def showData(self, figsize=(15, 10), nc=0, nr=0, clim=None):
        """show all data cubes in subplots"""
        from math import sqrt, ceil
        nsond = len(self.mrs)
        if nc == 0:
            nc = ceil(sqrt(nsond*3))
        if nr == 0:
            nr = ceil(nsond/nc)
        fig, ax = plt.subplots(nrows=int(nr), ncols=int(nc),
                               figsize=figsize)
        for i, mrs in enumerate(self.mrs):
            mrs.showCube(ax=ax.flat[i], vec=mrs.data*1e9, islog=False,
                         clim=clim)
        self.figs['data'] = fig
        return fig, ax

    def showInitialValues(self):
        """ show initial values of whole profile """
        IVI = np.zeros((len(self.mrs[0].q), len(self.mrs)))
        for i, mrs in enumerate(self.mrs):
            IVI[:, i] = mrs.dcube[:, 0]

        fig, ax = plt.subplots(figsize=(15, 10))
        im = ax.imshow(IVI*1e9, interpolation='nearest')
        plt.colorbar(im, ax=ax, orientation='horizontal')
        self.figs['ivi'] = fig
        return fig, ax

    def independentBlock1dInversion(self, nlay=2, lam=100, startModel=None):
        """Independent inversion of all soundings.

        Parameters
        ----------
        nlay : int [2]
            number of layers
        lam : float
            regularization parameter
        startModel : array/vector
            starting model (see MRS.run parameters)
        """
        self.WMOD, self.TMOD = [], []
        self.RMSvec, self.Chi2vec, self.nData = [], [], []
        for i, mrs in enumerate(self.mrs):
            if hasattr(startModel[0], '__iter__'):
                startvec = startModel[i]
            else:
                startvec = startModel
            if startvec is not None:
                nlay = (len(startvec)-1) // 2
            mrs.run(nlay=nlay, startvec=startvec, lam=lam, verbose=False)
            mrs.INV.echoStatus()
            thk, wc, t2 = mrs.result()
            self.WMOD.append(np.hstack((thk, wc)))
            self.TMOD.append(np.hstack((thk, t2)))
            self.RMSvec.append(mrs.INV.absrms()*1e9)  # in nV
            self.Chi2vec.append(mrs.INV.chi2())
            self.nData.append(len(mrs.data))

        self.totalChi2 = sum(np.array(self.Chi2vec)*np.array(self.nData)) / \
            sum(self.nData)
        self.totalRMS = np.sqrt(sum(np.array(self.RMSvec)**2 *
                                    np.array(self.nData)) / sum(self.nData))

    def block1dInversionOld(self, nlay=2, startModel=None, verbose=True,
                            uncertainty=False, **kwargs):
        """Invert all data together by one 1D model (variant 1 - all equal)."""
        self.mrsall = MRS()
        self.mrsall.z = self.mrs[0].z
        self.mrsall.t = self.mrs[0].t
        self.mrsall.data, self.mrsall.error, self.mrsall.q = [], [], []
        self.mrsall.K = np.zeros((0, self.mrs[0].K.shape[1]))
        for mrs in self.mrs:
            self.mrsall.data = np.hstack((self.mrsall.data, mrs.data))
            self.mrsall.error = np.hstack((self.mrsall.error,
                                           np.real(mrs.error)))
            self.mrsall.q = np.hstack((self.mrsall.q, mrs.q))
            self.mrsall.K = np.vstack((self.mrsall.K, mrs.K))

        self.mrsall.run(nlay, stvec=startModel, verbose=verbose, **kwargs)
        if verbose:
            self.mrsall.showResult()
        return self.mrsall.model

    def block1dInversion(self, nlay=2, lam=100., show=False, verbose=True,
                         uncertainty=False):
        """Invert all data together by a 1D model (more general solution)."""
        data, error = pg.Vector(), pg.Vector()
        for mrs in self.mrs:
            data = pg.cat(data, mrs.data)
            error = pg.cat(error, np.real(mrs.error))

#        f = JointMRSModelling(self.mrs, nlay)
        f = MultiFOP(self.mrs, nlay)
        mrsobj = self.mrs[0]
        for i in range(3):
            f.region(i).setParameters(mrsobj.startval[i], mrsobj.lowerBound[i],
                                      mrsobj.upperBound[i])

        INV = pg.Inversion(data, f, verbose)
        INV.setLambda(lam)
        INV.setMarquardtScheme(0.8)
#        INV.stopAtChi1(False)  # should be already in MarquardtScheme
        INV.setDeltaPhiAbortPercent(0.5)
        INV.setAbsoluteError(error)
        model = INV.run()
        m0 = self.mrs[0]
        m0.model = np.asarray(model)
        if uncertainty:
            from pygimli.utils import iterateBounds
            m0.modelL, m0.modelU = iterateBounds(
                INV, dchi2=INV.chi2() / 2, change=1.2)
        if show:
            self.show1dModel()
        # %% fill up 2D model (for display only)
        self.WMOD, self.TMOD = [], []
        thk = model[0:nlay-1]
        wc = model[nlay-1:2*nlay-1]
        t2 = model[2*nlay-1:3*nlay-1]
        for i in range(len(self.mrs)):
            self.WMOD.append(np.hstack((thk, wc)))
            self.TMOD.append(np.hstack((thk, t2)))

        return model

    def show1dModel(self):
        """Show 1D model (e.g. of joint block inversion)."""
        self.figs['1dmodel'], ax = self.mrs[0].showResult()

    def blockLCInversion(self, nlay=2, startModel=None, **kwargs):
        """Laterally constrained (piece-wise 1D) block inversion."""
        data, error, self.nData = pg.Vector(), pg.Vector(), []
        for mrs in self.mrs:
            data = pg.cat(data, mrs.data)
            error = pg.cat(error, mrs.error)
            self.nData.append(len(mrs.data))

        fop = MRSLCI(self.mrs, nlay=nlay)
        fop.region(0).setZWeight(kwargs.pop('zWeight', 0))
        fop.region(0).setConstraintType(kwargs.pop('cType', 1))
        transData, transMod = pg.trans.Trans(), pg.trans.TransLog()  # LU(1., 500.)
        if startModel is None:
            startModel = self.block1dInversion(nlay, verbose=False)
        model = kwargs.pop('startvec', np.tile(startModel, len(self.mrs)))
        INV = pg.Inversion(data, fop, transData, transMod, True, False)
        INV.setModel(model)
        INV.setReferenceModel(model)
        INV.setAbsoluteError(error)
        INV.setLambda(kwargs.pop('lam', 100))
        INV.setMaxIter(kwargs.pop('maxIter', 20))
#        INV.stopAtChi1(False)
        INV.setLambdaFactor(0.9)
        INV.setDeltaPhiAbortPercent(0.1)
        model = INV.run()
        self.WMOD, self.TMOD = [], []
        for par in np.reshape(model, (len(self.mrs), 3*nlay-1)):
            thk = par[0:nlay-1]
            self.WMOD.append(np.hstack((thk, par[nlay-1:2*nlay-1])))
            self.TMOD.append(np.hstack((thk, par[2*nlay-1:3*nlay-1])))

        ind = np.hstack((0, np.cumsum(self.nData)))
        resp = INV.response()
        misfit = data - resp
        emisfit = misfit / error
        misfit *= 1e9
        self.totalChi2 = INV.chi2()
        self.totalRMS = INV.absrms()*1e9
        self.RMSvec, self.Chi2vec = [], []
        for i in range(len(self.mrs)):
            self.RMSvec.append(np.sqrt(np.mean(misfit[ind[i]:ind[i+1]]**2)))
            self.Chi2vec.append(np.mean(emisfit[ind[i]:ind[i+1]]**2))

    def printFits(self):
        """Show single fits and total fit."""
        np.set_printoptions(precision=2)
        print("Single RMS [nV]:")
        print(np.array(self.RMSvec))
        print("Single Chi^2:")
        print(np.array(self.Chi2vec))
        print('Total RMS/Chi^2 value:')
        print(np.round(self.totalRMS, 2), np.round(self.totalChi2, 2))

    def showWC(self, wlim=(0, 0.5), ax=None, cmap='Spectral',
               title=r'$\theta$ (-)'):
        """Show water content distribution as stitched model section."""
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))
            self.figs['wc'] = fig

        showStitchedModels(self.WMOD, ax=ax, x=self.x, islog=False, cmap=cmap,
                           title=title, cmin=wlim[0], cmax=wlim[1])
        return fig, ax

    def showT2(self, tlim=(0.05, 0.5), ax=None, cmap='Spectral',
               title=r'$T_2^*$ (s)'):
        """Show relaxation time distribution as stitched model section."""
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))
            self.figs['t2'] = fig

        showStitchedModels(self.TMOD, ax=ax, x=self.x, islog=True, cmap=cmap,
                           cmin=tlim[0], cmax=tlim[1], title=title)
        return fig, ax

    def showFits(self, ax=None):
        """Show chi-square and rms fits of individual soundings."""
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))
            self.figs['fits'] = fig
        axb = ax.twinx()
        axb.set_ylabel('rms (nV)')
        ax.plot(self.x, self.Chi2vec, 'rx', label=r'$\chi^2$')
        axb.plot(self.x, self.RMSvec, 'bx', label='rms [nV]')
        ax.legend(numpoints=1, loc=2)
        axb.legend(numpoints=1, loc=1)

    def showModel(self, showFit=0, cmap='Spectral', figsize=(13, 12),
                  wlim=(0, 0.5), tlim=(0.05, 0.5)):
        """Show 2d model as stitched 1d models along with fit."""
        fig, ax = plt.subplots(nrows=2+showFit, figsize=figsize, sharex=True)
        self.showWC(wlim, ax=ax[-2], cmap=cmap)
        self.showT2(tlim, ax=ax[-1], cmap=cmap)
        xl = ax[-1].get_xlim()
        ax[0].set_xlabel('x (m)')
        ax[0].xaxis.set_label_position('top')
        ax[0].xaxis.tick_top()
        ax[0].set_ylabel(r'$\chi^2$ (-)')
        for axi in ax[-2:]:
            axi.set_ylabel('z (m)')
        if showFit > 0:
            self.showFits(ax[0])
            ax[0].set_xlim(xl)

        return fig, ax

    def saveFigs(self, basename="out", extension="pdf"):
        """Save all figures to (pdf) files."""
        for key in self.figs:
            self.figs[key].savefig(basename+"-"+key+"."+extension,
                                   bbox_inches='tight')

if __name__ == "__main__":
    name = sys.argv[-1]
    numlay = 3
    prof = MRSprofile(name)  # directory or list of names
    print(prof)
    figure, axes = prof.showData()  # subplots with data cubes
    mymodel = prof.block1dInversion(nlay=numlay)  # all in one model
    prof.independentBlock1dInversion(nlay=2)  # each separately
    prof.printFits()
    figure, axes = prof.showModel(showFit=1)
    figure.savefig(name+'-Ind-N'+str(numlay)+'.pdf', bbox_inches='tight')
    prof.blockLCInversion(nlay=numlay, lam=100, startModel=mymodel)
    prof.printFits()
    figure, axes = prof.showModel(showFit=1)
    figure.savefig(name+'-LCI-N'+str(numlay)+'.pdf', bbox_inches='tight')
