#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Different spectral induced polarization (SIP) tools."""

import numpy as np
import pygimli as pg
from .models import ColeColeComplex, ColeColeComplexSigma, PeltonPhiEM
from .models import ColeColeAbs, ColeColePhi, DoubleColeColePhi


def fitCCEMPhi(f, phi, ePhi=0.001, lam=1000.,
               mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100),
               cpar=(0.25, 0, 1), empar=(1e-7, 1e-9, 1e-5), verbose=True):
    """Fit a Cole-Cole term with additional EM term to phase."""
    fCCEM = PeltonPhiEM(f)
    fCCEM.region(0).setParameters(*mpar)    # m (start,lower,upper)
    fCCEM.region(1).setParameters(*taupar)  # tau
    fCCEM.region(2).setParameters(*cpar)   # c
    fCCEM.region(3).setParameters(*empar)   # tau-EM
    ICC = pg.core.Inversion(phi, fCCEM, verbose, verbose)  # set up inversion
    ICC.setAbsoluteError(ePhi)  # 1 mrad
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = ICC.run()  # run inversion
    if verbose:
        ICC.echoStatus()
    return model, np.asarray(ICC.response())


def fitCCPhi(f, phi, ePhi=0.001, lam=1000., verbose=True, robust=False,
             mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100), cpar=(0.25, 0, 1)):
    """Fit a Cole-Cole term with additional EM term to phase."""
    fCCEM = ColeColePhi(f)
    fCCEM.region(0).setParameters(*mpar)    # m (start,lower,upper)
    fCCEM.region(1).setParameters(*taupar)  # tau
    fCCEM.region(2).setParameters(*cpar)   # c
    ICC = pg.core.Inversion(phi, fCCEM, False)  # set up inversion class
    ICC.setAbsoluteError(ePhi)  # 1 mrad
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    ICC.setRobustData(robust)
#    ICC.setMaxIter(0)
    model = ICC.run()  # run inversion
    if verbose:
        ICC.echoStatus()
    return model, np.asarray(ICC.response()), ICC.chi2()
#    return np.append(model, ICC.chi2()), np.asarray(ICC.response())


def fit2CCPhi(f, phi, ePhi=0.001, lam=1000., verbose=True, robust=False,
              mpar1=(0.2, 0, 1), taupar1=(1e-2, 1e-5, 100), cpar1=(0.2, 0, 1),
              mpar2=(0.2, 0, 1), taupar2=(1e-2, 1e-5, 100), cpar2=(0.2, 0, 1)):
    """Fit a double Cole-Cole term to phase."""
    f2CC = DoubleColeColePhi(f)
    f2CC.region(0).setParameters(*mpar1)    # m (start,lower,upper)
    f2CC.region(1).setParameters(*taupar1)  # tau
    f2CC.region(2).setParameters(*cpar1)   # c
    f2CC.region(3).setParameters(*mpar2)    # m (start,lower,upper)
    f2CC.region(4).setParameters(*taupar2)  # tau
    f2CC.region(5).setParameters(*cpar2)   # c
    ICC = pg.core.Inversion(phi, f2CC, False)  # set up inversion class
    ICC.setAbsoluteError(ePhi)  # 1 mrad
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    ICC.setRobustData(robust)
    ICC.setDeltaPhiAbortPercent(1)
#    ICC.setMaxIter(0)
    model = ICC.run()  # run inversion
    if verbose:
        ICC.echoStatus()
    return np.asarray(model), np.asarray(ICC.response())


def fitCCAbs(f, amp, error=0.01, lam=1000., mstart=None,
             taupar=(1e-2, 1e-5, 100), cpar=(0.5, 0, 1)):
    """Fit amplitude spectrum by Cole-Cole model."""
    fCC = ColeColeAbs(f)
    tLog = pg.trans.TransLog()
    fCC.region(0).setStartValue(max(amp))
    if mstart is None:  # compute from amplitude decay
        mstart = 1. - min(amp) / max(amp)
    fCC.region(1).setParameters(mstart, 0, 1)    # m (start,lower,upper)
    fCC.region(2).setParameters(*taupar)  # tau
    fCC.region(3).setParameters(*cpar)   # c
    ICC = pg.core.Inversion(amp, fCC, tLog, tLog, False)  # set up inversion
    ICC.setRelativeError(error)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = np.asarray(ICC.run())  # run inversion
    ICC.echoStatus()
    response = np.asarray(ICC.response())
    return model, response


def fitCCC(f, amp, phi, eRho=0.01, ePhi=0.001, lam=1000., mstart=None,
           taupar=(1e-2, 1e-5, 100), cpar=(0.5, 0, 1), verbose=False):
    """Fit complex spectrum by Cole-Cole model."""
    fCC = ColeColeComplex(f)
    tLog = pg.trans.TransLog()
    fCC.region(0).setStartModel(max(amp))
    if mstart is None:  # compute from amplitude decay
        mstart = 1. - min(amp) / max(amp)
    fCC.region(1).setParameters(mstart, 0, 1)    # m (start,lower,upper)
    fCC.region(2).setParameters(*taupar)  # tau
    fCC.region(3).setParameters(*cpar)   # c
    data = pg.cat(amp, phi)
    ICC = pg.core.Inversion(data, fCC, verbose, verbose)  # set up inversion
    ICC.setTransModel(tLog)
    error = pg.cat(eRho*amp, pg.Vector(len(f), ePhi))
    ICC.setAbsoluteError(error)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = ICC.run()  # run inversion
    ICC.echoStatus()
    response = ICC.response()
    return model, response[:len(f)], response[len(f):], ICC.chi2()


def fitCCCC(f, amp, phi, error=0.01, lam=10., taupar=[1e-2, 1e-5, 100],
            cpar=[0.25, 0, 1], mpar=[0, 0, 1]):
    """Fit complex spectrum by Cole-Cole model based on sigma."""
    fCC = ColeColeComplexSigma(f)
    tLog = pg.trans.TransLog()
    fCC.region(0).setStartValue(1./max(amp))
    if mpar[0] == 0:
        mpar[0] = 1. - min(amp)/max(amp)
    fCC.region(1).setParameters(*mpar)    # m (start,lower,upper)
    fCC.region(2).setParameters(*taupar)  # tau
    fCC.region(3).setParameters(*cpar)   # c
    data = pg.cat(1./amp * np.cos(phi), 1./amp * np.sin(phi))
    ICC = pg.core.Inversion(data, fCC, False)  # set up inversion class
    ICC.setTransModel(tLog)
    ICC.setAbsoluteError(data*error+max(data)*0.0001)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = np.asarray(ICC.run())  # run inversion
    ICC.echoStatus()
    response = np.asarray(ICC.response())
    rRe, rIm = response[:len(f)], response[len(f):]
    rAmp = 1./np.sqrt(rRe**2+rIm**2)
    return model, rAmp, np.arctan(rIm/rRe), ICC.chi2()


if __name__ == "__main__":
    pass
