#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pygimli base functions to handle complex arrays"""

from math import pi
from pygimli.core.base import isComplex

import numpy as np
import pygimli as pg


isComplex = pg.isComplex

def toComplex(amp, phi=None):
    """Convert real values into complex (z = a + ib) valued array.

    If no phases phi are given, assuming z = amp[0:N] + i amp[N:2N].

    If phi is given in (rad) complex values are generated:
    z = amp*(cos(phi) + i sin(phi))

    Parameters
    ----------
    amp: iterable (float)
        Amplitudes or real unsqueezed real valued array.
    phi: iterable (float)
        Phases in neg radiant
    Returns
    -------
    z: ndarray(dtype=np.complex)
        Complex values
    """
    if phi is not None:
        return np.array(amp) * (np.cos(phi) + 1j *np.sin(phi))
    N = len(amp) // 2
    return np.array(amp[0:N]) + 1j * np.array(amp[N:])
    #return np.array(pg.math.toComplex(amp[0:N], amps[N:]))

def toPolar(z):
    """Convert complex (z = a + ib) values array into amplitude and phase in radiant

    If z is real valued we assume its squeezed.

    Parameters
    ----------
    z: iterable (floats, complex)
        If z contains of floats and squeezedComplex is assumed [real, imag]

    Returns
    -------
    amp, phi: ndarray
        Amplitude amp and phase angle phi in radiant.

    """
    if isComplex(z):
        return np.abs(z), np.angle(z)
    else:
        return toPolar(toComplex(z))

def squeezeComplex(z, polar=False, conj=False):
    """Squeeze complex valued array into [real, imag] or [amp, phase(rad)]"""
    if isinstance(z, (pg.matrix.CSparseMapMatrix,
                      pg.matrix.CSparseMatrix,
                      pg.matrix.CMatrix)):
        return toRealMatrix(z, conj=conj)

    if isComplex(z):
        vals = np.array(z)
        if conj:
            vals = np.conj(vals)

        if polar is True:
            vals = pg.cat(*toPolar(z))
        else:
            vals = pg.cat(vals.real, vals.imag)
        return vals
    return z

def toRealMatrix(C, conj=False):
    """Convert complex valued matrix into a real valued Blockmatrix

    Parameters
    ----------
    C: CMatrix
        Complex valued matrix
    conj: bool [False]
        Fill the matrix as complex conjugated matrix

    Returns
    -------
    R : pg.matrix.BlockMatrix()

    """
    R = pg.matrix.BlockMatrix()
    Cr = pg.math.real(A=C)
    Ci = pg.math.imag(A=C)

    rId = R.addMatrix(Cr)
    iId = R.addMatrix(Ci)
    # we store the mats in R to keep the GC happy after leaving the scope

    R.addMatrixEntry(rId, 0,         0,         scale=1.0)
    R.addMatrixEntry(rId, Cr.rows(), Cr.cols(), scale=1.0)
    if conj == True:
        pg.warn('Squeeze conjugate complex matrix.')
        R.addMatrixEntry(iId, 0,         Cr.cols(), scale=1.0)
        R.addMatrixEntry(iId, Cr.rows(), 0,         scale=-1.0)
    else:
        R.addMatrixEntry(iId, 0,         Cr.cols(), scale=-1.0)
        R.addMatrixEntry(iId, Cr.rows(), 0,         scale=1.0)
    return R

def KramersKronig(f, re, im, usezero=False):
    """Return real/imaginary parts retrieved by Kramers-Kronig relations.

    formulas including singularity removal according to Boukamp (1993)
    """
    from scipy.integrate import simpson

    x = 2. * pi * f  # omega
    im2 = np.zeros(im.shape)
    re2 = np.zeros(im.shape)
    drdx = np.diff(re) / np.diff(x)  # d Re/d omega
    didx = np.diff(im) / np.diff(x)  # d Im/d omega
    dRedx = np.hstack((drdx[0], (drdx[:-1] + drdx[1:]) / 2, drdx[-1]))
    dImdx = np.hstack((didx[0], (didx[:-1] + didx[1:]) / 2, didx[-1]))
    for num, w in enumerate(x):
        x2w2 = x**2 - w**2
        x2w2[num] = 1e-12

        fun1 = (re - re[num]) / x2w2
        fun1[num] = dRedx[num] / 2 / w
        im2[num] = -2./pi * w * simpson(fun1, x=x)

        if usezero:
            fun2 = (im * w / x - im[num]) / x2w2
            re2[num] = 2./pi * w * simpson(fun2, x=x)  + re[0]
        else:
            fun2 = (im * x - im[num] * w) / x2w2
            fun2[num] = (im[num] / w + dImdx[num]) / 2
            re2[num] = 2./pi * simpson(fun2, x=x) + re[-1]

    return re2, im2

