#!/usr/bin/env python
"""Physical constants repeatedly used in geophysics.

The :class:`Constants` container collects the vacuum permeability and
permittivity, the speed of light, the gravimetric constant and the
gyromagnetic ratio of protons in water.
"""

from math import pi, sqrt


class Constants:
    """Container class for some constants repeatedly used in geophysics."""

    # magnetic constant, vacuum permeability
    mu0 = 4.0 * pi * 1e-7  # [(kg * m) / (A^2 * s^2)]

    # electric constant, vacuum permittivity
    e0 = 8.85418781762e-12  # [(A^2 * s^4)/(kg m^3)]

    # speed of light in vacuum
    c0 = 1.0 / sqrt(mu0 * e0)
    # speed of light in vacuum
    c = c0

    # gravimetric constant
    G = 6.6742e-11  # [m^3/(kg s^2)]
    GmGal = G / 1e-5  # mGal

    Darcy = 9.86923e-13  # [m^2]

    # unit (median) acceleration on earth
    g = 9.80665  # [m/s^2] ex. 9.798
    # also a function of latitude (and altitude)?

    # gyromagnetic ratio of protons in water in T/s
    gammaP = 2.67515255e8  # noqa: N815 -- public API name

    # Larmor frequency of water
    fLarmor = gammaP / 2 / pi  # noqa: N815 -- public API name, in Hz
    fLarmorMhZ = fLarmor * 1e-6  # noqa: N815 -- public API name

    # nepers to dB multiplication factor
    nepers2dB = 8.686  # noqa: N815 -- public API name
