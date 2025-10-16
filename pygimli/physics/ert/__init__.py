# -*- coding: utf-8 -*-
"""Electrical Resistivity Tomography (ERT).

Direct-Current (DC) Resistivity and Induced Polarisation (IP)

This package contains tools, modelling operators, and managers for
Electrical Resistivity Tomography (ERT) & Induced polarization (IP)

Main entry functions or classes:
* simulate - synthetic (real or complex-valued) modelling
* createData - generate data sets for synthetic modelling
* ERTModelling - Modelling operator
* ERTManager - data inversion and modelling for real resistivity
* ERTIPManager - extension to IP (either frequency or time domain)
* TimelapseERT - processing and inversion of timelapse ERT data
* CrossholeERT - timelapse ERT in crosshole environments
"""

import pygimli as pg
from .ert import (simulate, estimateError,
                  createGeometricFactors, createInversionMesh)
from .ertManager import ERTManager
from .ertIPManager import ERTIPManager
from .ertModelling import ERTModelling, ERTModellingReference
from .ertScheme import createData
from .processing import (uniqueERTIndex, generateDataFromUniqueIndex,
                         reciprocalIndices, fitReciprocalErrorModel,
                         reciprocalProcessing, combineMultipleData,
                         removeDuplicates)
from .timelapse import TimelapseERT
from .crosshole import CrossholeERT
from pygimli.physics.ves import VESManager  # backward compatibility
from pygimli.physics.ves.vesModelling import VESModelling
# , VESCModelling, VESRhoModelling

from .visualization import showERTData, drawERTData, generateDataPDF
from .importData import load


@pg.renamed(createData, '1.3')  # 20210302
def createERTData(*args, **kwargs):
    pass

showData = showERTData
show = showERTData  # better create a function that can also handle mgr
pg.core.DataContainerERT.show.__doc__ = showERTData.__doc__
geometricFactor = pg.core.geometricFactors  # even more to be removed
geometricFactors = geometricFactor  # to be removed

# def __DataContainerERT_createGeometricFactors(self, *args, **kwargs):
#     self['k'] = createGeometricFactors(self, *args, **kwargs)
#     return self['k']

# __DataContainerERT_createGeometricFactors.__doc__ = createGeometricFactors.__doc__
# pg.core.DataContainerERT.createGeometricFactors = __DataContainerERT_createGeometricFactors

# Module prototypes
DataContainer = pg.core.DataContainerERT
Manager = ERTManager
Modelling = ERTModelling
coverageERT = pg.core.coverageDCtrans
