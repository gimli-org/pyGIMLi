# -*- coding: utf-8 -*-
"""General physics independent solver interface."""

from .utils import (anisotropyMatrix, constitutiveMatrix,
                    createAnisotropyMatrix, createConstitutiveMatrix)

from .green import greenDiffusion1D

from .solver import *
from .solver import cellValues
from .solverFiniteVolume import solveFiniteVolume


__all__ = []

class WorkSpace:
    pass

