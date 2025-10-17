#!/usr/bin/env python
"""Module containing submodules for various geophysical methods."""

from .constants import Constants
from .complexSpectrum import ComplexSpectrum
from .ert import ERTManager, ERTModelling
from .ves import VESManager
from .em import VMDTimeDomainModelling
from .traveltime import TravelTimeManager

Refraction = TravelTimeManager  # Backward compatibility to pg 1.0

from .em import FDEM, TDEM
from .SIP import SIPSpectrum, SpectrumManager
from .sNMR import MRS

constants = Constants

# from . gravimetry import Gravimetry
# from . seismics import *

# __all__ = ["FDEM", "TDEM", "MRS", "SIPSpectrum", "Refraction"]
