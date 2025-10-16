#!/usr/bin/env python
"""Refraction seismics or first arrival traveltime calculations."""


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
    'drawFirstPicks',
    'createRAData',
    'createGradientModel2D',
    'RefractionNLayer',
    'RefractionNLayerFix1stLayer',
    'shotReceiverDistances',
    'TravelTimeManager',
    'TravelTimeDijkstraModelling'
]
