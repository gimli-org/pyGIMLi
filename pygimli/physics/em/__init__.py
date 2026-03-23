"""
Electromagnetics (EM) in 1D (for 3D we refer to custEM(x))

Frequency-domain (FD) or time-domain (TD) semi-analytical 1D solutions.

Manager classes:

* FDEM - Frequency-domain (two-loop) EM
* TDEM - Time-domain (in-loop) EM

Modelling operators:

* HEMmodelling - airborne EM modelling (with elevation)
* MT1dBlockModelling - magnetotelluric (MT) block model
* MT1dSmoothModelling - magnetotelluric (MT) smooth model
* TDEMBlockModelling - time-domain block model
* TDEMSmoothModelling - time-domain block model

Utility functions:

* readusffile, importMaxminData - import files
* rhoafromB, rhoafromU - conversion to apparent resictivity
"""

from .vmd import VMDTimeDomainModelling

from .fdem import FDEM
from .tdem import TDEM, rhoafromB, rhoafromU
from .tdem import VMDTimeDomainModelling, TDEMSmoothModelling
from .mt1dmodelling import MT1dBlockModelling, MT1dSmoothModelling

MT1dModelling = MT1dBlockModelling  # default
TDEMBlockModelling = VMDTimeDomainModelling  # better name
TDEMOccamModelling = TDEMSmoothModelling  # alias

from .hemmodelling import HEMmodelling
from .io import readusffile, importMaxminData
from .tools import cmapDAERO, xfplot, FDEMsystems
