#!/usr/bin/env python
"""Spectral induced polarization (SIP) spectrum class and modules."""

import sys
import codecs

from math import log10, exp, pi
import numpy as np

import pygimli as pg
from pygimli.utils import isComplex, squeezeComplex, toComplex, KramersKronig

from .importData import readTXTSpectrum, readFuchs3File, readRadicSIPFuchs

from .plotting import drawAmplitudeSpectrum, drawPhaseSpectrum, showSpectrum
from .models import DebyePhi, DebyeComplex, relaxationTerm, DoubleColeCole
from .tools import fitCCEMPhi, fitCCC
from .tools import fitCCCC, fitCCPhi, fit2CCPhi

from pygimli.frameworks import MethodManager
from pygimli.frameworks import ParameterModelling


class SpectrumModelling(ParameterModelling):
    """Modelling framework with an array of freqencies as data space.

    Attributes
    ----------
    params: dict
    function: callable
    complex: bool

    """

    def __init__(self, funct=None, **kwargs):
        """Initialize.

        Parameters
        ----------
        func : function
            modelling function
        complex : bool
            complex function
        frequencies : iterable
            frequency vector
        """
        self._complex = kwargs.pop("complex", False)
        super().__init__(funct=funct, **kwargs)
        self.defaultModelTrans = 'log'
        self._freqs = kwargs.pop("frequencies", None)

    @property
    def complex(self):
        """Return if spectrum is complex."""
        return self._complex

    @complex.setter
    def complex(self, c):
        self._complex = c

    @property
    def freqs(self):
        """Return frequency vector."""
        if self._freqs is None:
            pg.critical("No frequencies defined.")
        return self._freqs

    @freqs.setter
    def freqs(self, f):
        self._freqs = f
        self.dataSpace = self._freqs

    def response(self, params):
        """Model response.

        Parameters
        ----------
        params : iterable
            model

        Returns
        -------
        array_like
            model response vector

        """
        ret = super().response(params)
        if self.complex:
            return squeezeComplex(ret)
        return ret

    def drawData(self, ax, data, err=None, **kwargs):
        """Draw data."""
        if self.complex:
            Z = toComplex(data)
            showSpectrum(self.freqs, np.abs(Z), -np.angle(Z)*1000,
                         axs=ax, **kwargs)
        else:
            ax.semilogx(self.freqs, data)
            ax.legend()


class SpectrumManager(MethodManager):
    """Manager to work with spectra data."""

    def __init__(self, fop=None, **kwargs):
        """Set up spectrum manager.

        Parameters
        ----------
        fop : pg.Modelling operator, optional
            Forward operator. The default is None.
        **kwargs : TYPE
            passed to SpectrumManager.

        """
        self._funct = fop
        super().__init__(**kwargs)

    def setFunct(self, fop, **kwargs):
        """Set forward modelling function."""
        self._funct = fop
        self.reinitForwardOperator(**kwargs)

    def createForwardOperator(self, **kwargs):
        """Create a Forward operator."""
        if isinstance(self._funct, SpectrumModelling):
            return self._funct

        fop = SpectrumModelling(self._funct, **kwargs)
        return fop

    def createInversionFramework(self, **kwargs):
        """Create inversion framework."""
        return pg.frameworks.MarquardtInversion(**kwargs)

    def simulate(self):
        """Make a simulation."""
        pass

    def setData(self, freqs=None, amp=None, phi=None, eAmp=0.03, ePhi=0.001):
        """Set data for chosen sip model.

        Parameters
        ----------
        freqs: iterable
            Array-like frequencies.
        amp: iterable
            Array-like amplitudes to work with.
        phi: iterable
            Array-like phase angles to work with.
        eAmp: float|iterable
            Relative error for amplitudes.
        ePhi: float|iterable
            Absolute error for phase angles.
        """
        self.fop.freqs = freqs
        if phi is not None:
            self.fw.dataVals = self._ensureData(toComplex(amp, phi))
            if isinstance(ePhi, float):
                ePhi = abs(ePhi/phi)
        else:
            self.fw.dataVals = self._ensureData(amp)

        if isinstance(eAmp, float):
            eAmp = np.ones(len(freqs)) * eAmp

        err = np.asarray([*eAmp, *ePhi]) if phi is not None else eAmp

        self.fw.errorVals = self._ensureError(err, self.fw.dataVals)

    def _ensureData(self, data):
        """Check data validity."""
        if isinstance(data, pg.DataContainer):
            pg.critical("Implement me")

        if data is None:
            data = self.fw.dataVals

        vals = data
        if isComplex(data):
            self.fop.complex = True
            vals = squeezeComplex(data)

        if abs(min(vals)) < 1e-12:
            print(min(vals), max(vals))
            pg.critical("There are zero data values.")

        return vals

    def _ensureError(self, err, dataVals=None):
        """Check data validity."""
        if isinstance(err, pg.DataContainer):
            pg.critical("Implement me")

        if err is None:
            err = self.fw.errorVals

        vals = err
        if vals is None:
            return self._ensureError(0.01, dataVals)

        if abs(min(vals)) < 1e-12:
            print(min(vals), max(vals))
            pg.critical("There are zero data values.")

        return vals

    def invert(self, data=None, f=None, **kwargs):
        """Invert the spectrum."""
        if f is not None:
            self.fop.freqs = f

        limits = kwargs.pop('limits', {})

        for k, v in limits.items():
            self.fop.setRegionProperties(k, limits=v)

            if 'startmodel' not in kwargs:
                sm = (v[1] + v[0]) / 2
                if v[0] > 0:
                    sm = np.exp(np.log(v[0]) +
                                (np.log(v[1]) - np.log(v[0])) / 2.)

                self.fop.setRegionProperties(k, startModel=sm)

        return super().invert(data, **kwargs)

    def showResult(self):
        """Show resulting data."""
        ax = None
        if self.fop.complex:
            fig, ax = pg.plt.subplots(nrows=2, ncols=1)
        else:
            fig, ax = pg.plt.subplots(nrows=1, ncols=1)

        self.fop.drawModel(ax, self.fw.model)
        self.fop.drawData(ax, self.fw.dataVals, label='data')
        self.fop.drawData(ax, self.fw.response, label='response')

        return ax


class SIPSpectrum:
    """SIP spectrum data analysis."""

    def __init__(self, filename=None, unify=False, onlydown=True,
                 f=None, amp=None, phi=None, k=1, sort=True,
                 basename='new', **kwargs):
        """Init SIP class with either filename to read or data vectors.

        Examples
        --------
        >>> #sip = SIPSpectrum('sipexample.txt', unify=True) # unique f values
        >>> #sip = SIPSpectrum(f=f, amp=R, phi=phase, basename='new')
        """
        self._verbose = False
        self.basename = basename
        self.fig = {}
        self.k = k
        self.f = None  # mandatory frequencies
        self.amp = None  # mandatory amplitudes
        self.phi = None  # mandatory phases

        self.epsilon0 = 8.854e-12

        if filename is not None:
            self.loadData(filename, **kwargs)
        else:
            if f is not None:
                self.f = np.asarray(f)
            if amp is not None:
                self.amp = np.asarray(amp)
            if phi is not None:
                self.phi = np.asarray(phi)

        if unify and self.amp is not None:
            self.unifyData(onlydown)
        if sort and self.amp is not None:
            self.sortData()

        self.ampOrg = None
        self.phiOrg = None
        self.phiCC = None  # better make a struct of m, amp, phi (tau)
        self.ampCC = None
        self.mCC = None
        self.phiDD = None
        self.ampDD = None
        self.mDD = None
        self.tau = None

    def __repr__(self):
        """Human readable string representation of the class."""
        out = self.__class__.__name__ + " object"
        if self.f is not None and hasattr(self.f, '__iter__'):
            out += "\nnf=" + str(len(self.f)) + " min/max="
            out += str(min(self.f)) + "/" + str(max(self.f))
        return out

    def loadData(self, filename, **kwargs):
        """Import spectral data.

        Import Data and try to assume the file format.
        """
        verbose = kwargs.pop('verbose', self._verbose)

        with codecs.open(filename, 'r', encoding='iso-8859-15',
                         errors='replace') as f:
            firstLine = f.readline()

        fnLow = filename.lower()
        self.basename = filename[:-4]
        if 'SIP Fuchs III' in firstLine:
            if verbose:
                pg.info("Reading SIP Fuchs III file")
            self.f, self.amp, self.phi, self.header = readFuchs3File(
                filename, verbose=verbose, **kwargs)
            self.phi *= -np.pi/180.
        elif 'SIP-Quad' in firstLine:
            if verbose:
                pg.info("Reading SIP Quad file")
            self.f, self.amp, self.phi, self.header = readFuchs3File(
                filename, nfr=9, namp=10, nphi=11, nk=7,
                verbose=verbose, **kwargs)
            self.phi *= -np.pi/180.
        elif 'SIP-Fuchs' in firstLine:
            if verbose:
                pg.info("Reading SIP Fuchs file")
            self.f, self.amp, self.phi, drhoa, dphi = readRadicSIPFuchs(
                filename, verbose=verbose, quad='SIP-Quad' in firstLine,
                **kwargs)
            self.phi *= -np.pi/180.
        elif fnLow.endswith('.txt') or fnLow.endswith('.csv'):
            self.f, self.amp, self.phi = readTXTSpectrum(filename, **kwargs)
        else:
            try:
                out = np.genfromtxt(filename, names=True)
                self.f = out["FreqHz"]
                self.amp = out["AppResOhmm"]
                self.phi = -out["Phasedeg"] * np.pi / 180
            except BaseException as e:
                raise ValueError("Don't know how to read data.") from e

        self.amp *= self.k
        return self.f, self.amp, self.phi

    def unifyData(self, onlydown=False):
        """Unify data (only one value per frequency) by mean or selection."""
        if self.f is None:
            return
        fu = np.unique(self.f)
        if len(fu) < len(self.f) or onlydown:
            if onlydown:
                nonzero = np.nonzero(np.diff(self.f) > 0)[0]
                if len(nonzero) > 0:
                    wende = min(nonzero)
                    if wende > 0:
                        self.f = self.f[wende::-1]
                        self.amp = self.amp[wende::-1]
                        self.phi = self.phi[wende::-1]
            else:
                amp = np.zeros(fu.shape)
                phi = np.zeros(fu.shape)
                for i in range(len(fu)):
                    ind = np.nonzero(self.f == fu[i])[0]
                    amp[i] = np.mean(self.amp[ind])
                    phi[i] = np.mean(self.phi[ind])

                self.f = fu
                self.amp = amp
                self.phi = phi

    def sortData(self):
        """Sort data along increasing frequency (e.g. useful for KK)."""
        if self.f is None:
            return

        ind = np.argsort(self.f)
        self.amp = self.amp[ind]
        self.phi = self.phi[ind]
        self.f = self.f[ind]

    def cutF(self, fcut=1e99, down=False):
        """Cut (delete) frequencies above a certain value fcut."""
        ind = self.f >= fcut if down else self.f <= fcut

        self.amp = self.amp[ind]
        self.phi = self.phi[ind]
        if np.any(self.phiOrg):
            self.phiOrg = self.phiOrg[ind]
        if np.any(self.ampOrg):
            self.ampOrg = self.ampOrg[ind]

        self.f = self.f[ind]
#        self.amp = self.amp[self.f <= fcut]
#        self.phi = self.phi[self.f <= fcut]
#        if np.any(self.phiOrg):
#            self.phiOrg = self.phiOrg[self.f <= fcut]
#        # finally cut f
#        self.f = self.f[self.f <= fcut]

    def omega(self):
        """Angular frequency."""
        return self.f * 2 * np.pi

    def realimag(self, cond=False):
        """Real and imaginary part of resistivity/conductivity (cond=True)."""
        amp = 1.0 / self.amp if cond else self.amp
        return amp * np.cos(self.phi), amp * np.sin(self.phi)

    def zNorm(self):
        """Normalize real (difference) and imag. z :cite:`NordsiekWel2008`."""
        re, im = self.realimag()
        R0 = max(self.amp)
        zNormRe = 1. - re / R0
        zNormIm = im / R0
        return zNormRe, zNormIm

    def showPhase(self, ax=None, **kwargs):
        """Plot phase spectrum (-phi over frequency)."""
        if ax is None:
            fig, ax = pg.plt.subplots()
            self.fig['phase'] = fig

        drawPhaseSpectrum(ax, self.f, self.phi*1000, **kwargs)
        return ax

    def showData(self, reim=False, znorm=False, cond=False, nrows=2, ax=None,
                 **kwargs):
        """Show amplitude and phase spectrum in two subplots.

        Parameters
        ----------
        reim : bool
            show real/imaginary part instead of amplitude/phase

        znorm : bool (true forces reim)
            use normalized real/imag parts after Nordsiek&Weller (2008)

        nrows - use nrows subplots (default=2)

        Returns
        -------
        fig, ax : mpl.figure, mpl.axes array
        """
        if reim or znorm or cond:
            addstr = ''
            if znorm:
                re, im = self.zNorm()
                addstr = ' (norm)'
            else:
                re, im = self.realimag(cond=cond)

            fig, ax = showSpectrum(self.f, re, im, ylog=cond, nrows=nrows,
                                   axs=ax, **kwargs)
            self.fig['data'] = fig
            ax[0].set_ylabel('real part'+addstr)
            ax[1].set_ylabel('imaginary part'+addstr)
        else:
            fig, ax = showSpectrum(self.f, self.amp, self.phi*1000,
                                   axs=ax, **kwargs)
            self.fig['data'] = fig

        ax[0].set_title(kwargs.pop("title", self.basename))
        return fig, ax

    def getKK(self, use0=False):
        """Compute Kramers-Kronig impedance values (re->im and im->re)."""
        re, im = self.realimag()
        if False:
            ind = np.argsort(self.f)
            reKK, imKK = KramersKronig(self.f[ind], re[ind], im[ind],
                                       usezero=use0)
        else:
            fsort, ind = np.unique(self.f, return_index=True)

        reKK, imKK = KramersKronig(fsort, re[ind], im[ind], usezero=use0)
        re[ind] = reKK  # sort back
        im[ind] = imKK
        return re, im

    def getPhiKK(self, use0=False):
        """Compute phase from Kramers-Kronig quantities."""
        _, imKK = self.getKK(use0)
        re, _ = self.realimag()
        return np.arctan2(imKK, re)

    def showDataKK(self, use0=False):
        """Show data as real/imag subplots along with Kramers-Kronig curves."""
        fig, ax = self.showData(reim=True)
        reKK, imKK = self.getKK(use0)
        ax[0].semilogx(self.f, reKK, label='KK')
        ax[1].semilogx(self.f, imKK, label='KK')
        for a in ax:
            a.set_yscale('linear')
            a.legend()

        self.fig['dataKK'] = fig
        return fig, ax

    def checkCRKK(self, useEps=False, use0=False, ax=None):
        """Check coupling removal (CR) by Kramers-Kronig (KK) relation."""
        if ax is None:
            fig, ax = pg.plt.subplots()
            self.fig['dataCRKK'] = fig

        ax.semilogx(self.f, self.phi*1000, "+-", label='org')
        ax.semilogx(self.f, self.getPhiKK(use0)*1000, "x-", label='orgKK')
        if useEps:
            self.removeEpsilonEffect()
        else:
            self.fitCCEM()

        ax.semilogx(self.f, self.phi*1000, "+--", label='corr')
        ax.semilogx(self.f, self.getPhiKK(use0)*1000, "+--", label='corrKK')
        ax.grid(True)
        ax.legend(loc='best')

    def showPolarPlot(self, cond=False):
        """Show data in a polar plot (imaginary vs. real parts)."""
        re, im = self.realimag(cond=cond)
        fig, ax = pg.plt.subplots()
        self.fig['polar'] = fig
        ax.plot(re, im, 'b.')
        ax.set_aspect(1)
        ax.grid(True)
        for i in range(0, len(re), 5):
            fi = self.f[i]
            mul = 10**np.floor(np.log10(fi) - 2)  # 3 counting digits
            ax.text(re[i], im[i], str(int(fi / mul) * mul))

        return fig, ax

    def epsilonR(self):
        """Calculate relative permittivity from imaginary conductivity."""
        _, sigmaI = self.realimag(cond=True)

        return sigmaI / (self.f * 2 * pi * self.epsilon0)

    def determineEpsilon(self, mode=0, sigmaR=None, sigmaI=None):
        """Retrieve frequency-independent epsilon for f->Inf.

        Parameters
        ----------
        mode : int
            Operation mode:
                =0 - extrapolate using two highest frequencies (default)
                <0 - take last -n frequencies
                >0 - take n-th frequency
        sigmaR/sigmaI : float
            real and imaginary conductivity (if not given take data)

        Returns
        -------
        er : float
            relative permittivity (epsilon) value (dimensionless)
        """
        if sigmaR is None or sigmaI is None:
            sigmaR, sigmaI = self.realimag(cond=True)

        epsr = sigmaI / self.omega() / self.epsilon0
        nmax = np.argmax(self.f)
        if mode == 0:
            f1 = self.f * 1
            a = f1[nmax] * 1.
            f1[nmax] = 0
            nmax1 = np.argmax(f1)
            a = a / f1[nmax1]
            er = (a*epsr[nmax] - epsr[nmax1]) / (a-1)
        elif mode < 0:
            er = np.median(epsr[mode:])
        else:
            er = epsr[nmax]

        return er

    def removeEpsilonEffect(self, er=None, mode=0):
        """Remove effect of (constant high-frequency) epsilon from sigma.

        Parameters
        ----------
        er : float
            relative epsilon to correct for (else automatically determined)
        mode : int
            automatic epsilon determination mode (see determineEpsilon)

        Returns
        -------
        er : float
            determined permittivity (see determineEpsilon)
        """
        sigR, sigI = self.realimag(cond=True)
        if er is None:  #
            er = self.determineEpsilon(mode=mode, sigmaR=sigR, sigmaI=sigI)
            print("detected epsilon of ", er)

        sigI -= er * self.omega() * self.epsilon0
        self.phiOrg = self.phi
        self.phi = np.arctan(sigI/sigR)
        self.ampOrg = self.amp
        self.amp = 1. / np.sqrt(sigR**2 + sigI**2)
        return er

    def fitCCPhi(self, ePhi=0.001, lam=1000., mpar=(0, 0, 1),
                 taupar=(0, 1e-5, 100), cpar=(0.3, 0, 1), verbose=False):
        """Fit a Cole-Cole term (to phase only).

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        mpar, taupar, cpar : list[3]
            inversion parameters (starting value, lower bound, upper bound)
            for Cole-Cole parameters (m, tau, c) and EM relaxation time (em)

        """
        if taupar[0] == 0:
            taupar = (1.0 / self.f[np.argmax(self.phi)] / 2.0 / pi,
                      taupar[1], taupar[2])
#            taupar[0] = 1.0 / self.f[np.argmax(self.phi)] / 2.0 / pi
            print("taupar", taupar)
        if mpar[0] == 0:
            mpar = (1. - min(self.amp)/max(self.amp), mpar[1], mpar[2])
#            mpar[0] = 1. - min(self.amp)/max(self.amp)
            print("mpar", mpar)
        self.mCC, self.phiCC, self.chi2 = fitCCPhi(
            self.f, self.phi, ePhi, lam, mpar=mpar, taupar=taupar, cpar=cpar)

    def fit2CCPhi(self, ePhi=0.001, lam=1000., mpar=(0, 0, 1),
                  taupar1=(0, 1e-5, 1), taupar2=(0, 1e-1, 1000),
                  cpar=(0.5, 0, 1), verbose=False):
        """Fit two Cole-Cole terms (to phase only).

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        mpar : list[3]
            starting value, lower bound, upper bound for chargeability
        taupar1 / taupar2 : list[3]
            starting value, lower bound, upper bound for 2 time constants
        cpar1 / cpar2 : list[3]
            starting value, lower bound, upper bound for 2 relaxation exponents
        """
        if taupar1[0] == 0:
            taupar1 = (np.sqrt(taupar1[1]*taupar1[2]), taupar1[1], taupar1[2])
            print("taupar1", taupar1)
        if taupar2[0] == 0:
            taupar2 = (np.sqrt(taupar2[1]*taupar2[2]), taupar2[1], taupar2[2])
            print("taupar2", taupar2)
#            taupar1[0] = 1.0 / self.f[np.argmax(self.phi)] / 2.0 / pi
        if mpar[0] == 0:
            mpar = (1. - min(self.amp)/max(self.amp), mpar[1], mpar[2])  # *2
            print("mpar", mpar)
        self.mCC, self.phiCC = fit2CCPhi(self.f, self.phi, ePhi, lam,
                                         mpar1=mpar, mpar2=mpar,
                                         cpar1=cpar, cpar2=cpar,
                                         taupar1=taupar1, taupar2=taupar2)

    def fitCCEM(self, ePhi=0.001, lam=1000., remove=True,
                mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100),
                cpar=(0.25, 0, 1), empar=(1e-7, 1e-9, 1e-5), verbose=False):
        """Fit a Cole-Cole term with additional EM term to phase.

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        remove: bool
            remove EM term from data
        mpar, taupar, cpar, empar : list[3]
            inversion parameters (starting value, lower bound, upper bound)
            for Cole-Cole parameters (m, tau, c) and EM relaxation time (em)
        """
        self.mCC, self.phiCC = fitCCEMPhi(self.f, self.phi, ePhi, lam, mpar,
                                          taupar, cpar, empar, verbose=verbose)
        # correct EM term from data
        if remove:
            self.phiOrg = self.phi
            self.phi = self.phi + \
                np.angle(relaxationTerm(self.f, self.mCC[3]))

    def fitColeCole(self, useCond=False, **kwargs):
        """Fit a Cole-Cole model to the phase data.

        Parameters
        ----------
        useCond : bool
            use conductivity form of Cole-Cole model instead of resistivity
        error : float [0.01]
            absolute phase error
        lam : float [1000]
            initial regularization parameter
        mpar : tuple/list (0, 0, 1)
            inversion parameters for chargeability: start, lower, upper bound
        taupar : tuple/list (1e-2, 1e-5, 100)
            inversion parameters for time constant: start, lower, upper bound
        cpar : tuple/list (0.25, 0, 1)
            inversion parameters for Cole exponent: start, lower, upper bound
        """
        if useCond:  # use conductivity formulation instead of resistivity
            self.mCC, self.ampCC, self.phiCC, self.chi2 = fitCCCC(
                self.f, self.amp, self.phi, **kwargs)
            self.mCC[0] = 1. / self.mCC[0]
        else:
            self.mCC, self.ampCC, self.phiCC, self.chi2 = fitCCC(
                self.f, self.amp, self.phi, **kwargs)

    def fitDoubleColeCole(self, ePhi=0.001, eAmp=0.01, lam=1000., robust=False,
                          verbose=True, useRho=True, useMult=False, aphi=True,
                          mpar1=(0.2, 0, 1), mpar2=(0.2, 0, 1), tauRho=True,
                          taupar1=(1e-2, 1e-5, 100), taupar2=(1e-4, 1e-5, 100),
                          cpar1=(0.5, 0, 1), cpar2=(0.5, 0, 1)):
        """Fit double Cole-Cole term to complex resistivity or phase.

        Parameters
        ----------
        ePhi : float [0.001]
            absolute error of phase angle in rad
        eAmp : float [0.01 = 1%]
            absolute error of phase angle
        lam : float
            regularization parameter
        robust : bool [False]
            use robustData (L1 norm on data side)
        useRho : bool [True]
            Cole-Cole defined for impedance/resistivity, otherwise conductivity
        useMult : bool [False]
            the two terms are combined as product (otherwise sum)
        tauRho : bool [False]
            in case of useRho=False the time constant is defined like for rho
        mpar1/2, taupar1/2, cpar1/2 : list[3]
            inversion parameters (starting value, lower bound, upper bound)
            for Cole-Cole parameters (m, tau, c)

        """
        f2CC = DoubleColeCole(self.f, rho=useRho, aphi=aphi, tauRho=False,
                              mult=useMult)
        if useRho:
            rhoStart = min(self.amp)
            f2CC.region(0).setParameters(rhoStart, 0., rhoStart*10)
        else:
            sigStart = 1./max(self.amp)
            f2CC.region(0).setParameters(sigStart, 0., sigStart*10)

        f2CC.region(1).setParameters(*mpar1)    # m (start,lower,upper)
        f2CC.region(2).setParameters(*taupar1)  # tau
        f2CC.region(3).setParameters(*cpar1)   # c
        f2CC.region(4).setParameters(*mpar2)    # m (start,lower,upper)
        f2CC.region(5).setParameters(*taupar2)  # tau
        f2CC.region(6).setParameters(*cpar2)   # c
        if aphi:
            amp = self.amp if useRho else 1./self.amp
            data = np.hstack((amp, self.phi))
            error = np.hstack((np.ones_like(amp)*eAmp,
                               ePhi / np.abs(self.phi)))
        else:
            re, im = self.realimag(not useRho)
            data = np.hstack((re, im))
            error = np.ones(len(self.f)*2) * eAmp

        ICC = pg.core.Inversion(data, f2CC, False)  # set up inversion class
        ICC.setRelativeError(error)  # 1 mrad
        ICC.setLambda(lam)  # start with large damping and cool later
        ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
        ICC.setRobustData(robust)
        ICC.setDeltaPhiAbortPercent(1)
    #    ICC.setMaxIter(0)
        self.mCC = ICC.run()  # run inversion
        self.chi2 = ICC.chi2()
        if verbose:
            ICC.echoStatus()

        one, two = np.split(ICC.response(), 2)
        if aphi:
            self.ampCC = one
            self.phiCC = two

    def fitDebyeModel(self, ePhi=0.001, lam=1e3, lamFactor=0.8,
                      tau=None, mint=None, maxt=None, nt=None, useComplex=True,
                      showFit=False, verbose=False, **kwargs):
        """Fit a (smooth) continuous Debye model (Debye decomposition).

        Parameters
        ----------
        ePhi : float
            absolute error of phase angle
        lam : float
            regularization parameter
        lamFactor : float
            regularization factor for subsequent iterations
        mint/maxt : float
            minimum/maximum tau values to use (else automatically from f)
        nt : int
            number of tau values (default number of frequencies * 2)
        new : bool
            new implementation (experimental)
        showFit : bool
            show fit
        cType : int
            constraint type (1/2=smoothness 1st/2nd order, 0=minimum norm)
        phi : iterable
            use phi instead of self.phi
        """
        nf = len(self.f)
        if tau is None:
            if mint is None:
                mint = .1 / max(self.f)
            if maxt is None:
                maxt = .5 / min(self.f)
            if nt is None:
                nt = nf*2
            self.tau = np.logspace(log10(mint), log10(maxt), nt)
        else:
            self.tau = tau
        # discretize tau, setup DD and perform DD inversion
        startModel = pg.Vector(len(self.tau), 0.01)
        new = kwargs.pop("new", useComplex)  # renamed kwargs
        if new:
            reNorm, imNorm = self.zNorm()
            fDD = DebyeComplex(self.f, self.tau)
            Znorm = pg.cat(reNorm, imNorm)
            IDD = pg.Inversion(fop=fDD)
            absErr = max(Znorm)*0.003+ePhi
            self.mDD = IDD.run(Znorm, absoluteError=absErr,
                               startModel=startModel,
                               lam=lam, lambdaFactor=lamFactor,
                               **kwargs)
            IDD.echoStatus()
        else:
            fDD = DebyePhi(self.f, self.tau)
            IDD = pg.Inversion(fop=fDD)
            phi = kwargs.pop("phi", self.phi)
            self.mDD = IDD.run(phi, absoluteError=ePhi, startModel=startModel,
                               lam=lam, lambdaFactor=lamFactor, **kwargs)

        self.invDD = IDD
        if new:
            if verbose:
                pg.info("ARMS=", IDD.absrms(), "RRMS=", IDD.absrms()/max(Znorm)*100)

            resp = np.array(IDD.response)
            respRe = resp[:nf]
            respIm = resp[nf:]
            respC = ((1 - respRe) + respIm * 1j) * max(self.amp)
            self.phiDD = np.angle(respC)
            self.ampDD = np.abs(respC)
            if showFit:
                fig, ax = self.showData(znorm=True, nrows=3)
                self.fig['DebyeFit'] = fig
                ax[0].plot(self.f, respRe, 'r-')
                ax[1].plot(self.f, respIm, 'r-')
                ax[2].semilogx(self.tau, self.mDD, 'r-')
                ax[2].set_xlim(max(self.tau), min(self.tau))
                ax[2].set_ylim(0., max(self.mDD))
                ax[2].grid(True)
                ax[2].set_xlabel(r'$\tau$ (s)')
                ax[2].set_ylabel('$m$ (-)')

        else:
            self.phiDD = IDD.response
            if showFit:
                fig, ax = self.showData(nrows=3)
                self.fig['DebyeSpectrum'] = fig
                ax[2].semilogx(self.tau, self.mDD, 'r-')

    def totalChargeability(self):
        """Total chargeability (sum) from Debye curve."""
        return sum(self.mDD)

    def logMeanTau(self):
        """Mean logarithmic relaxation time (50% cumulative log curve)."""
        return exp(np.sum(np.log(self.tau) * self.mDD) / sum(self.mDD))

    def createDecay(self, t=None):
        """Create decay from Debye decomposition."""
        v = np.zeros(len(t))
        for tau, m in zip(self.tau, self.mDD, strict=False):
            v += np.exp(-t/tau) * m

        return v

    def showAll(self, save=False, ax=None):
        """Plot spectrum, Cole-Cole fit and Debye distribution."""
        if np.any(self.mCC):  # generate title strings
            mCC = self.mCC
            rstr = r'$\rho$={:.4f} '
            cstr = r'CC: m={:.3f} $\tau$={:.1e}s c={:.2f} '
            if len(mCC) == 7:  # double Cole-Cole
                tstr = rstr + cstr + cstr
            elif len(mCC) == 6:  # double Cole-Cole only Phi
                tstr = cstr + cstr
            elif len(mCC) == 4:
                tstr = rstr + cstr
            else:
                tstr = cstr

            tCC = tstr.format(*mCC)

        if ax is None:
            fig, ax = pg.plt.subplots(nrows=2+(self.mDD is not None),
                                      figsize=(12, 12))
        else:
            fig = ax[0].figure

        self.fig['all'] = fig
        fig.subplots_adjust(hspace=0.25)
        # amplitude
        drawAmplitudeSpectrum(ax[0], self.f, self.amp,
                              label='data', ylog=0)
        if np.any(self.ampDD):
            ax[0].plot(self.f, self.ampDD, 'm-', label='DD response')
        if np.any(self.ampCC):
            ax[0].semilogx(self.f, self.ampCC, 'r-', label='CC model')
        ax[0].legend(loc='best')
        # phase
        if np.any(self.ampOrg):
            ax[0].semilogx(self.f, self.ampOrg, 'cx-', label='org. data')
            ax[0].legend(loc='best')
        if np.any(self.phiOrg):
            ax[1].semilogx(self.f, self.phiOrg * 1e3, 'cx-', label='org. data')

        ax[1].semilogx(self.f, self.phi * 1e3, 'b+-', label='data')
        if np.any(self.phiCC):
            ax[1].semilogx(self.f, self.phiCC * 1e3, 'r-', label='CC model')
        if np.any(self.phiDD):
            ax[1].semilogx(self.f, self.phiDD * 1e3, 'm-', label='DD model')

        ax[1].grid(True)
        ax[1].legend(loc='best')
        ax[1].set_xlabel('f (Hz)')
        ax[1].set_ylabel('-phi (mrad)')
        if np.any(self.mCC):
            ax[1].set_title(tCC, loc='left')
        if np.any(self.mDD):
            mtot = self.totalChargeability()
            lmtau = self.logMeanTau()
            tDD = rf'DD: m={mtot:.3f} $\tau$={lmtau:.1e}s'
            ax[2].semilogx(self.tau, self.mDD * 1e3)
            ax[2].set_xlim(ax[2].get_xlim()[::-1])
            ax[2].grid(True)
            ax[2].set_xlabel(r'$\tau$ (s)')
            ax[2].set_ylabel('m (mV/V)')
            ax[2].set_title(tDD, loc='left')

        if save:
            savename = save if isinstance(save, str) else self.basename + '.pdf'

            fig.savefig(savename, bbox_inches='tight')

        pg.plt.show(block=False)
        return fig, ax

    def saveFigures(self, name=None, ext='pdf'):
        """Save all existing figures to files using file basename."""
        if name is None:
            name = self.basename
        if name is None or not any(name):
            name = 'out'
        for key in self.fig:
            self.fig[key].savefig(name+'-'+key+'.'+ext, bbox_inches='tight')


def run_SIPSPectrum(myfile):
    """Run typical SIP spectrum workflow for a given file."""
    sip = SIPSpectrum(myfile)
    # sip.showData(znorm=True)
    if True:  # Pelton
        sip.fitCCEM()
    else:
        sip.removeEpsilonEffect()
        sip.fitColeCole(useCond=False)

    sip.fitDebyeModel()  # , showFit=True)
    # create titles and plot data, fit and model
    sip.showAll(save=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No filename given, falling back to test case")
        run_SIPSPectrum('sipexample.txt')
    else:
        run_SIPSPectrum(sys.argv[1])
