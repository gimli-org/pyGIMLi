#!/usr/bin/env python
"""Frequency-domain (FD) or time-domain (TD) semi-analytical 1D solutions."""

from .vmd import VMDTimeDomainModelling

from .fdem import FDEM
from .tdem import TDEM, rhoafromB, rhoafromU
from .tdem import VMDTimeDomainModelling, TDEMSmoothModelling

TDEMBlockModelling = VMDTimeDomainModelling  # better name
TDEMOccamModelling = TDEMSmoothModelling  # alias

from .hemmodelling import HEMmodelling
from .io import readusffile, importMaxminData
from .tools import cmapDAERO, xfplot, FDEMsystems
