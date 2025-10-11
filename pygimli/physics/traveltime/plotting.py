#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plotting functions for traveltime."""
import numpy as np

import pygimli as pg
from pygimli.viewer.mpl import createColorBar
from .utils import shotReceiverDistances


def drawTravelTimeData(ax, data, t=None):
    """Draw first arrival traveltime data into mpl ax a.

    data of type pg.DataContainer must contain sensorIdx 's' and 'g'
    and thus being numbered internally [0..n)
    """
    x = pg.x(data.sensorPositions())
    # z = pg.z(data.sensorPositions())

    shots = pg.unique(pg.sort(data('s')))
    geoph = pg.unique(pg.sort(data('g')))

    startOffsetIDX = 0

    if min(min(shots), min(geoph)) == 1:
        startOffsetIDX = 1

    tShow = data('t')
    if t is not None:
        tShow = t

    ax.set_xlim([min(x), max(x)])
    ax.set_ylim([max(tShow), -0.002])
    ax.figure.show()  # a draw function should never trigger a figure show

    for shot in shots:
        gIdx = pg.find(data('s') == shot)
        sensorIdx = [int(i__ - startOffsetIDX) for i__ in data('g')[gIdx]]
        ax.plot(x[sensorIdx], tShow[gIdx], 'x-')

    yPixel = ax.transData.inverted().transform_point((1, 1))[1] - \
        ax.transData.inverted().transform_point((0, 0))[1]
    xPixel = ax.transData.inverted().transform_point((1, 1))[0] - \
        ax.transData.inverted().transform_point((0, 0))[0]

    # draw shot points
    ax.plot(x[[int(i__ - startOffsetIDX) for i__ in shots]],
            np.zeros(len(shots)) + 8. * yPixel, 'gv', markersize=8)

    # draw geophone points
    ax.plot(x[[int(i__ - startOffsetIDX) for i__ in geoph]],
            np.zeros(len(geoph)) + 3. * yPixel, 'r^', markersize=8)

    ax.grid()
    ax.set_ylim([max(tShow), + 16. * yPixel])
    ax.set_xlim([min(x) - 5. * xPixel, max(x) + 5. * xPixel])

    ax.set_xlabel('x-Coordinate [m]')
    ax.set_ylabel('Traveltime [ms]')


def plotFirstPicks(ax, data, tt=None, plotva=False, marker='x-'):
    """Naming convention. drawFOO(ax, ... )."""
    pg.deprecated("use drawFirstPicks")
    return drawFirstPicks(ax=ax, data=data, tt=tt, plotva=plotva,
                          marker=marker)


def drawFirstPicks(ax, data, tt=None, plotva=False, **kwargs):
    """Plot first arrivals as lines.

    Parameters
    ----------
    ax : matplotlib.axes
        axis to draw the lines in
    data : :gimliapi:`GIMLI::DataContainer`
        data containing shots ("s"), geophones ("g") and traveltimes ("t").
        (`:py:class:pygimli.physics.traveltime.DataContainerTT.`)
    Return
    ------
    gci : list
        list of plotting items (matplotlib lines)
    """
    px = pg.x(data)
    gx = np.array([px[int(g)] for g in data["g"]])
    sx = np.array([px[int(s)] for s in data["s"]])
    if tt is None:
        tt = data["t"]

    if plotva:
        tt = np.absolute(gx - sx) / tt

    uns = np.unique(sx)

    cols = pg.plt.cm.tab10(np.arange(10))
    kwargs.setdefault('marker', 'x')
    kwargs.setdefault('markersize', 8)
    kwargs.setdefault('linestyle', '-')
    plotSource = kwargs.pop('plotSource', True)
    GCI = []
    for i, si in enumerate(uns):
        ti = tt[sx == si]
        gi = gx[sx == si]
        ii = gi.argsort()
        GCI.append(ax.plot(gi[ii], ti[ii], color=cols[i % 10], **kwargs))
        if plotSource:
            ax.plot(si, 0., 's', color=cols[i % 10])

    ax.grid(True)
    if plotva:
        ax.set_ylabel("Apparent velocity (m/s)")
    else:
        ax.set_ylabel("Traveltime (s)")

    ax.set_xlabel("x (m)")
    ax.invert_yaxis()
    return ax


# better be renamed to showData and optionaly show first pick curves
def showVA(data, usePos=True, ax=None, **kwargs):
    """Show apparent velocity as image plot.

    Parameters
    ----------
    data : pg.DataContainer()
        Datacontainer with 's' and 'g' Sensorindieces and 't' traveltimes.
    """
    ax, _ = pg.show(ax=ax)
    gci = drawVA(ax, data=data, usePos=usePos, **kwargs)

    cBar = createColorBar(gci, **kwargs)

    return ax, cBar


def drawVA(ax, data, vals=None, usePos=True, pseudosection=False, **kwargs):
    """Draw apparent velocities as matrix into an axis.

    Parameters
    ----------
    ax : mpl.Axes

    data : pg.DataContainer()
        Datacontainer with 's' and 'g' Sensorindieces and 't' traveltimes.

    usePos: bool [True]
        Use sensor positions for axes tick labels

    pseudosection : bool [False]
        Show in pseudosection style.

    vals : iterable
        Traveltimes, if None data need to contain 't' values.
    """
    if isinstance(vals, str):
        vals = data(vals)

    if vals is None:
        vals = data('t')

    px = pg.x(data)
    gx = np.asarray([px[g] for g in data.id("g")])
    sx = np.asarray([px[s] for s in data.id("s")])

    offset = shotReceiverDistances(data, full=True)

    if min(vals) < 1e-10:
        print(vals)
        pg.error('zero traveltimes found.')
    va = offset / vals

    if pseudosection:
        midpoint = (gx + sx) / 2
        gci = pg.viewer.mpl.dataview.drawVecMatrix(ax, midpoint, offset, va,
                                                   queeze=True,
                                                   label=pg.unit('as'))
    else:
        gci = pg.viewer.mpl.dataview.drawVecMatrix(ax, data["g"], data["s"],
                                                   va, squeeze=True,
                                                   label=pg.unit('as'))

    # A = np.ones((data.sensorCount(), data.sensorCount())) * np.nan
    # for i in range(data.size()):
    #     A[int(data('s')[i]), int(data('g')[i])] = va[i]
    # gci = ax.imshow(A, interpolation='nearest')
    # ax.grid(True)

    if usePos:
        nt = np.maximum(data.sensorCount() // 50, 10)
        xt = np.arange(0, data.sensorCount(), nt)
        ax.set_xticks(xt)
        ax.set_xticklabels([str(int(px[xti])) for xti in xt])
        ax.set_yticks(xt)
        ax.set_yticklabels([str(int(px[xti])) for xti in xt])

    return gci
