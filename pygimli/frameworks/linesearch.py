# -*- coding: utf-8 -*-
"""pyGIMLi - Linesearch.

Linesearch procedures used by various inversion frameworks.
"""
import numpy as np


def tauVector(taumin=0.01, taumax=1, logScale=False, n=21):
    """Generate a vector of tau values."""
    if logScale:
        return np.logspace(np.log10(taumin), np.log10(taumax), n)
    else:
        return np.linspace(taumin, taumax, n)

def lineSearchArmijo(inv, dM, tau=100, c=1e-4, tauFactor=0.8):
    """Line search using Armijos rule, typical for Steepest Descent.

    Start with a "largish" value alpha and decrease it by a factor until the
    Armijo rul is fulfilled (c is a small constant between 0 and 1):

    Phi(m+tau*m) <= Phi(m) + c alpha v^T dM
    
    where n is the negative normed gradient, i.e. -dM/norm(dM)
    """
    v = -dM / np.linalg.norm(dM)
    Elim = inv.phi() + sum(v*dM) * c * tau
    while Enew > Elim:
        newModel = inv.modelTrans.update(inv.model, v*tau)
        newResponse = inv.fop.response(newModel)
        Enew = inv.phi(newModel, newResponse)
        tau *= tauFactor

    return tau/tauFactor, newResponse

def lineSearchExact(inv, dM, taus=None, show=False, **kwargs):
    """Line search by exact forward response.

    Parameters
    ----------
    inv : pg.Inversion
        pygimli Inversion (or any derived class) instance
    taus : array
        array containing the tau values to test, alternatively:
    taumin : float
        minimum value
    taumax : float
        maximum value
    logScale : bool
        use logarithmic scaling, otherwise linear
    show : bool
        show curve
    """
    if taus is None:
        taus = tauVector(**kwargs)

    phis = np.zeros_like(taus)
    for i, tau in enumerate(taus):
        newModel = inv.modelTrans.update(inv.model, dM*tau)
        newResponse = inv.fop.response(newModel)
        phis[i] = inv.phi(newModel, newResponse)

    if show:
        import matplotlib.pyplot as plt
        if kwargs.get("logScale", False):
            plt.semilogx(taus, phis)
        else:
            plt.plot(taus, phis)

    return taus[np.argmin(phis)], newResponse

def lineSearchInter(inv, dM, taus=None, show=False, **kwargs):
    """Optimizes line search parameter by linear response interpolation.

    Parameters
    ----------
    inv : pg.Inversion
        pygimli Inversion (or any derived class) instance
    taus : array
        array containing the tau values to test, alternatively:
    taumin : float
        minimum value
    taumax : float
        maximum value
    logScale : bool
        use logarithmic scaling, otherwise linear
    show : bool
        show curve
    """
    if taus is None:
        taus = tauVector(**kwargs)

    phis = np.zeros_like(taus)
    dT = inv.dataTrans
    oldResponse = dT(inv.response)
    fullModel = inv.modelTrans.update(inv.model, dM)
    fullResponse = dT(inv.fop.response(fullModel))
    for i, tau in enumerate(taus):
        newModel = inv.modelTrans.update(inv.model, dM*tau)
        newResponse = dT.inv(oldResponse + (fullResponse - oldResponse) * tau)
        phis[i] = inv.phi(newModel, newResponse)

    if show:
        import matplotlib.pyplot as plt
        if kwargs.get("logScale", False):
            plt.semilogx(taus, phis)
        else:
            plt.plot(taus, phis)

    return taus[np.argmin(phis)], newResponse


def lineSearchInterOld(inv, dM, nTau=100, maxTau=1.0):
    """Optimizes line search parameter by linear response interpolation."""
    tD = inv.dataTrans
    tM = inv.modelTrans
    model = inv.model
    response = inv.response
    modelLS = tM.update(model, dM * maxTau)
    responseLS = inv.fop.response(modelLS)
    taus = np.arange(1, nTau+1) / nTau * maxTau
    phi = np.ones_like(taus) * inv.phi()
    phi[-1] = inv.phi(modelLS, responseLS)
    t0 = tD.fwd(response)
    t1 = tD.fwd(responseLS)
    for i, tau in enumerate(taus):
        modelI = tM.update(model, dM*tau)
        responseI = tD.inv(t1*tau+t0*(1.0-tau))
        phi[i] = inv.phi(modelI, responseI)

    return taus[np.argmin(phi)], responseLS

def lineSearchQuad(inv, dm, tautest=0.3, tau1=1, show=False, **kwargs):
    """Optimize line search by fitting parabola by Phi(tau) curve."""
    y0 = inv.phi()
    x1 = tau1
    fullModel = np.exp(dm*x1)*inv.model
    fullResponse = inv.fop.response(fullModel)
    xt = tautest
    testModel = np.exp(dm*xt)*inv.model
    testResponse = inv.fop.response(testModel)
    y1 = inv.phi(fullModel, fullResponse)
    yt = inv.phi(testModel, testResponse)
    rt = (yt-y0) / xt
    r1 = (y1-y0) / x1
    a = (rt - r1) / (xt - x1)
    b = (-rt*x1 + r1*xt) / (xt - x1)
    xopt = -b/a/2
    # xopt = (rt*x1 - r1*xt) / (rt - r1) / 2
    if show:
        taus = np.arange(0, 1.001, 0.01)
        import matplotlib.pyplot as plt
        ax = plt.subplots()[1]
        ax.plot(taus, taus**2*a+taus*b+y0, label="parabola")
        ax.plot(0, y0, "*", label="start")
        ax.plot(xt, yt, "*", label="test")
        ax.plot(x1, y1, "*", label="full")
        ax.plot(xopt, a*xopt**2+b*xopt+y0, "*", label="min")
        ax.grid()
        ax.legend()

    return xopt, None

def lineSearch(inv, dm, method='auto', **kwargs):
    """Carry out line search.

    Optimize step length s such that
    m + s*dm
    is minimized.

    Parameter
    ---------
    inv : pg.Inversion
        Inversion instance
    dm : iterable
        model update direction
    method : str ['auto']
        Method to be used:
        'exact' : function evaluation for every step
        'interp' : linear interpolation of response
        'quad' : fitting a parabola through 3 points
        'auto': first try 'inter', then 'quad', else 0.1
    taus : array [None]
        array containing the tau values to test, alternatively:
    taumin : float [0.01]
        minimum value
    taumax : float [1]
        maximum value
    logScale : bool [False]
        use logarithmic scaling, otherwise linear
    show : bool [False]
        show line search curve
    """
    if method.lower().startswith("exact"):
        return lineSearchExact(inv, dm, **kwargs)
    elif method.lower().startswith("int"):
        return lineSearchInter(inv, dm, **kwargs)
    elif method.lower().startswith("quad"):
        return lineSearchQuad(inv, dm, **kwargs)
    else:
        tau, response = lineSearchInter(inv, dm, **kwargs)
        if tau > 0.01 and tau <= 1:
            return tau, response

        tau, response = lineSearchQuad(inv, dm, **kwargs)
        if tau > 0.01 and tau <= 1:
            return tau, response
        else:
            return 0.1, None

if __name__ == "__main__":
    pass