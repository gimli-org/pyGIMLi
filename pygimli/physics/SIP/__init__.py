#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Spectral induced polarization (SIP) measurements and fittings.

* SIPSpectrum - main class for working with data

Modelling classes

* ColeColePhi
* DoubleColeColePhi
* ColeColeAbs
* ColeColeComplex
* ColeColeComplexSigma
* PeltonPhiEM
* DebyePhi
* DebyeComplex

Models:

* modelColeColeRho
* modelColeColeRhoDouble
* modelColeColeSigma
* modelColeColeSigmaDouble
* modelColeColeEpsilon
* modelColeDavidson

Functions:

* tauRhoToTauSigma
* showSpectrum
* drawPhaseSpectrum
* drawAmplitudeSpectrum
"""

from .importData import load
from .sipspectrum import SpectrumManager, SpectrumModelling, SIPSpectrum

from .models import (
    modelColeColeRho,
    modelColeColeRhoDouble,
    modelColeColeSigma,
    modelColeColeSigmaDouble,
    modelColeColeEpsilon,
    modelColeDavidson,

    ColeColePhi,
    DoubleColeColePhi,
    ColeColeAbs,
    ColeColeComplex,
    ColeColeComplexSigma,
    PeltonPhiEM,
    DebyePhi,
    DebyeComplex,

    tauRhoToTauSigma
)

from .plotting import showSpectrum, drawPhaseSpectrum, drawAmplitudeSpectrum

#__all__ = [name for name in dir() if '_' not in name]
#__all__ = [SIPSpectrum, ColeColeRho, ColeColePhi, DoubleColeColePhi]
