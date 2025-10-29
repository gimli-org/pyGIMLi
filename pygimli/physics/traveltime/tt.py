#!/usr/bin/env python
"""General convenience functions for traveltime module."""

# general purpose
import numpy as np
import pygimli as pg
from pygimli.viewer.mpl import createColorBar
from .TravelTimeManager import TravelTimeManager
from .plotting import drawVA, drawFirstPicks
Manager = TravelTimeManager  # convenience alias


def simulate(mesh, scheme, slowness=None, **kwargs):
    """Simulate traveltime data."""
    mgr = TravelTimeManager()
    return mgr.simulate(mesh=mesh, scheme=scheme, slowness=slowness, **kwargs)


simulate.__doc__ = TravelTimeManager.simulate.__doc__


class DataContainerTT(pg.DataContainer):
    """Data Container for traveltime."""

    def __init__(self, data=None, **kwargs):
        """Initialize empty data container, load or copy existing."""
        if isinstance(data, pg.DataContainer):
            super().__init__(data, **kwargs)
            self.registerSensorIndex("s")
            self.registerSensorIndex("g")
            self.setSensorIndexOnFileFromOne(True)
        else:
            super().__init__(**kwargs)
            self.registerSensorIndex("s")
            self.registerSensorIndex("g")
            if isinstance(data, str):
                self.load(data)


def show(data, **kwargs):
    """Show data."""
    ax, _ = pg.show(ax=kwargs.pop("ax", None))
    va = kwargs.pop("va", None)
    if va is None:  # check if refraction
        va = len(np.unique(pg.x(data))) < data.sensorCount()

    if va:
        gci = drawVA(ax, data=data, **kwargs)
        cBar = createColorBar(gci, **kwargs)
    else:
        drawFirstPicks(ax, data, tt=kwargs.pop("t", None), **kwargs)
        # drawTravelTimeData(ax, data, **kwargs)
        cBar = None

    return ax, cBar
