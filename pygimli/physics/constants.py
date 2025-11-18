#!/usr/bin/env python
# -*- coding: utf-8 -*-
from math import pi, sqrt


class Constants(object):
    """Container class for some constants repeatedly used in geophysics."""
    # magnetic constant, vacuum permeability
    mu0 = 4.0 * pi * 1e-7  # [(kg * m) / (A^2 * s^2)]

    # electric constant, vacuum permittivity
    e0 = 8.85418781762e-12  # [(A^2 * s^4)/(kg m^3)]

    # speed of light in vacuum
    c0 = 1. / sqrt(mu0 * e0)
    # speed of light in vacuum
    c = c0

    # gravimetric constant
    G = 6.6742e-11  # [m^3/(kg s^2)]
    GmGal = G / 1e-5  # mGal

    Darcy = 9.86923e-13  # [m^2]
    # unit (median) acceleration on earth
    g = 9.80665  # [m/s^2] ex. 9.798
    # also a g function of latitude (and altitude)?

    # gyromagnetic ratio of protons in water
    gammaP = 2.67515255e8  # T/s
    # Larmor frequency of water
    fLarmor = gammaP / 2 / pi  # in Hz
    fLarmorMhZ = fLarmor * 1e-6
    # nepers to dB multiplication factor
    nepers2dB = 8.686
