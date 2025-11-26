#!/usr/bin/env python
"""
First-arrival traveltime

e.g. refraction or crosshole seismics and GPR

Classes:

* TravelTimeManager
* RefractionNLayer
* DataContainer[TT]

Main entry functions:

* load - load or import from various formats
* createRAData - create refraction data scheme
* simulate - synthetic model computation
* show - show first arrivals als curves or image
"""

from pygimli.core import Dijkstra
from .importData import load
from .tt import simulate, DataContainerTT, show
from .plotting import drawFirstPicks, drawTravelTimeData, drawVA, showVA
from .utils import (createGradientModel2D, createRAData, shotReceiverDistances,
                    createCrossholeData)
#from .refraction import Refraction, Tomography # will be removed(201909)
from .refraction1d import RefractionNLayer, RefractionNLayerFix1stLayer
from .TravelTimeManager import TravelTimeDijkstraModelling, TravelTimeManager

Manager = TravelTimeManager
DataContainer = DataContainerTT
TravelTimeModelling = TravelTimeDijkstraModelling

__all__ = [
    'drawTravelTimeData',
    'drawVA',
    'showVA',
    'simulate',
    'show',
    'load',
    'Dijkstra',
    'drawFirstPicks',
    'createRAData',
    'createGradientModel2D',
    'createCrossholeData',
    'RefractionNLayer',
    'RefractionNLayerFix1stLayer',
    'shotReceiverDistances',
    'TravelTimeManager',
    'TravelTimeDijkstraModelling'
]
