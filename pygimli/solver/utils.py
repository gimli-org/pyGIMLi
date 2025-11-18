#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions for the FE/FV solver."""

import numpy as np
import pygimli as pg


def createAnisotropyMatrix(lon, trans, theta):
    """Create anisotropy matrix with desired properties.

    Anistropy tensor from longitudinal value lon,
    transverse value trans and the angle theta of the symmetry axis relative to the vertical after  cite:WieseGreZho2015
    https://www.researchgate.net/publication/249866312_Explicit_expressions_for_the_Frechet_derivatives_in_3D_anisotropic_resistivity_inversion

    TODO
    ----
        * 3D
    """
    C = np.zeros((2,2))
    C[0,0] = lon * np.cos(theta)**2 + trans * np.sin(theta)**2
    C[0,1] = 0.5 * (-lon + trans) * np.sin(theta) * np.cos(theta)
    C[1,0] = 0.5 * (-lon + trans) * np.sin(theta) * np.cos(theta)
    # Check what is correct .. the papers are divergend
    # C[0,1] = 0.5 * (-valL + valT) * np.sin(2*theta) * np.cos(theta)
    # C[1,0] = 0.5 * (-valL + valT) * np.sin(2*theta) * np.cos(theta)
    C[1,1] = lon * np.sin(theta)**2 + trans * np.cos(theta)**2
    return C

@pg.renamed(createAnisotropyMatrix)
def anisotropyMatrix(*args, **kwrags):
    return createAnisotropyMatrix(*args, **kwrags)

class ConstitutiveMatrix(np.ndarray):
    def __new__(cls, input_array, voigtNotation=False):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.voigtNotation = voigtNotation
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.voigtNotation = getattr(obj, 'voigtNotation', False)


def createConstitutiveMatrix(lam=None, mu=None, E=None, nu=None, dim=2,
                             voigtNotation=False):
    """Create constitutive matrix for 2 or 3D isotropic media.

    Either give lam and mu or E and nu.

    TODO
    ----
        * dim == 1
        * Tests
        * Examples
        * comparision Voigts/Kelvin notation

    Parameters
    ----------
    lam: float [None]
        1. Lame' constant
    mu: float [None]
        2. Lame' constant = G = Schubmodul
    E: float [None]
        Young's Modulus
    nu: float [None]
        Poisson's ratio
    voigtNotation: bool [False]
        Return in Voigt's notation instead of Kelvin's notation [default].

    Returns
    -------
    C: mat
        Either 3x3 or 6x6 matrix depending on the dimension
    """
    if voigtNotation is True:
        # pg._r('C-Voigt')
        a = 1 # Voigts notation
    else:
        # pg._y('C-Kelvin')
        a = 2 # Kelvins' notation

    if lam is not None and mu is not None:
        nu = (lam) / (2*(lam + mu))
        E = mu * (3*lam + 2*mu) / (lam + mu)

    if E is not None and nu is not None:
        if nu == 0.5 or nu >= 1.0:
            pg.critical('nu should be greater or smaller than 0.5 and < 1')
        lam = (E * nu) / ((1 + nu) * (1-2*nu))
        mu  = E / (2*(1 + nu))

    if dim == 1:
        pg.critical('not yet implemented')
    elif dim == 2:
        ## for pure 2d plane stress
        C = ConstitutiveMatrix(np.zeros((3, 3)),
                               voigtNotation=voigtNotation)
        C[0, 0] = 1
        C[1, 1] = 1
        C[0, 1] = nu
        C[1, 0] = nu
        C[2, 2] = (1-nu)/2 * a
        C *= E/(1-nu**2)
        return C

    # C = np.zeros((6, 6))
    # C[0, 0] = 1-nu
    # C[1, 1] = 1-nu
    # C[2, 2] = 1-nu
    # C[0, 1] = nu
    # C[1, 2] = nu
    # C[0, 2] = nu
    # C[1, 0] = nu
    # C[2, 1] = nu
    # C[2, 0] = nu
    # C[3, 3] = (1-2*nu)/2 * a
    # C[4, 4] = (1-2*nu)/2 * a
    # C[5, 5] = (1-2*nu)/2 * a
    # C *= E/((1+nu)*(1-2*nu))
    # print('c1', C)

    C = ConstitutiveMatrix(np.zeros((6, 6)),
                            voigtNotation=voigtNotation)
    C[0][0:3] = lam
    C[1][0:3] = lam
    C[2][0:3] = lam
    C[0][0] += 2. * mu
    C[1][1] += 2. * mu
    C[2][2] += 2. * mu

    C[3][3] = mu * a
    C[4][4] = mu * a
    C[5][5] = mu * a

    #print('c2', C)

    return C

@pg.renamed(createConstitutiveMatrix)
def constitutiveMatrix(*args, **kwrags):
    return createConstitutiveMatrix(*args, **kwrags)
