#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unified and method independent inversion frameworks."""


from .modelling import (Modelling, Block1DModelling, MeshModelling,
                        JointModelling, PriorModelling, LinearModelling,
                        PetroModelling, LCModelling, ParameterModelling)

from .inversion import (Inversion, MarquardtInversion,
                        Block1DInversion,
                        LCInversion)

from .lsqrinversion import LSQRInversion  # circular import, why?

from .methodManager import (fit, MethodManager, MethodManager1d,
                            ParameterInversionManager,
                            MeshMethodManager,
                            PetroInversionManager,
                            JointPetroInversionManager)

from .timelapse import MultiFrameModelling

from .linesearch import lineSearch

from .resolution import resolutionMatrix

from .harmfit import HarmFunctor, harmfit, harmfitNative

__all__ = ['HarmFunctor', 'harmfitNative', 'harmfit']
