#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""View ERT data."""

from math import pi
import numpy as np
from numpy import ma

import pygimli as pg
from pygimli.viewer.mpl.dataview import showValMapPatches
from pygimli.viewer.mpl import showDataContainerAsMatrix


def generateDataPDF(data, filename="data.pdf", **kwargs):
    """Generate a multi-page pdf showing all data properties."""
    if isinstance(data, str):
        filename = data.replace('.txt', '-data.pdf')
        data = pg.load(data)

    from matplotlib.backends.backend_pdf import PdfPages
    logToks = ["Uout(V)", "u", "i", "r", "rhoa"]
    with PdfPages(filename) as pdf:
        fig = pg.plt.figure()
        for tok in data.tokenList().split():
            if data.haveData(tok):
                vals = data[tok]
                logScale = min(vals) > 0 and tok in logToks
                ax = fig.add_subplot()
                pg.show(data, vals, ax=ax, label=tok, logScale=logScale, **kwargs)
                fig.savefig(pdf, format='pdf')
                fig.clf()

def showERTData(data, vals=None, **kwargs):
    """Plot ERT data as pseudosection matrix (position over separation).

    Creates figure, axis and draw a pseudosection.

    Parameters
    ----------
    data : :gimliapi:`BERT::DataContainerERT`

    **kwargs :

        * vals : array[nData] | str
            Values to be plotted. Default is data['rhoa'].
            Can be array or string whose data field is extracted.
        * axes : matplotlib.axes
            Axes to plot into. By default (None), a new figure with
            a single Axes is created.
        * x and y : str | list(str)
            forces using matrix plot (drawDataContainerAsMatrix)
            x, y define the electrode number on x and y axis
            can be strings ("a", "m", "mid", "sep") or lists of them ["a", "m"]
        * style : str
            predefined styles for choosing x and y arguments (x/y overrides)
            - "ab-mn" (default):  any combination of current/potential electrodes
            - "a-m" : only a and m electrode (for unique dipole spacings like DD)
            - "a-mn" : a and combination of mn electrode (PD with different MN)
            - "ab-m" : a and combination of mn electrode
            - "sepa-m" : current dipole length with a and m (multi-gradient)
            - "a-sepm" : a and potential dipole length with m
        * switchxy : bool
            exchange x and y axes before plotting

    Returns
    -------
    ax : matplotlib.axes
        axis containing the plots
    cb : matplotlib.colorbar
        colorbar instance
    """
    # remove ax keyword globally (problems with colorbar)
    ax = kwargs.pop('ax', None)
    if ax is None:
        fig = pg.plt.figure()
        ax = None
        axTopo = None
        if 'showTopo' in kwargs:  # can be removed?
            ax = fig.add_subplot(1, 1, 1)
#            axs = fig.subplots(2, 1, sharex=True)
#            # Remove horizontal space between axes
#            fig.subplots_adjust(hspace=0)
#            ax = axs[1]
#            axTopo = axs[0]
        else:
            ax = fig.add_subplot(1, 1, 1)

    pg.checkAndFixLocaleDecimal_point(verbose=False)

    if vals is None:
        vals = 'rhoa'

    if isinstance(vals, str):
        if data.haveData(vals):
            kwargs.setdefault('label', pg.utils.unit(vals))
            vals = data(vals)
        else:
            pg.critical('field not in data container: ', vals)

    kwargs.setdefault('cMap', pg.utils.cMap('rhoa'))  # better vals?
    kwargs.setdefault('label', pg.utils.unit('rhoa'))
    kwargs.setdefault('logScale', min(vals) > 0.0)

    sty = kwargs.pop("style", None)
    if isinstance(sty, str):
        sty = sty.lower()
    if sty == "a-m":
        kwargs.setdefault("y", "a")
        kwargs.setdefault("x", "m")
    elif sty == "a-mn":
        kwargs.setdefault("y", "a")
        kwargs.setdefault("x", ["m", "n"])
    elif sty == "ab-m":
        kwargs.setdefault("y", ["a", "b"])
        kwargs.setdefault("x", "m")
    elif sty == "sepa-m":
        data["ab"] = np.abs(data["b"] - data["a"])
        kwargs.setdefault("y", ["ab", "a"])
        kwargs.setdefault("x", "m")
    elif sty == "a-sepm":
        data["mn"] = np.abs(data["n"] - data["m"])
        kwargs.setdefault("x", "a")
        kwargs.setdefault("y", ["mn", "m"])
    elif sty is not None and sty != 0:
        kwargs.setdefault("y", ["a", "b"])
        kwargs.setdefault("x", ["m", "n"])

    if "x" in kwargs and "y" in kwargs:
        if kwargs["y"] == "mid":
            kwargs["y"] = (data["a"] + data["b"]) / 2
        elif kwargs["y"] == "sep":
            kwargs["y"] = np.abs(data["a"] - data["b"])

        if kwargs["x"] == "mid":
            kwargs["x"] = (data["m"] + data["n"]) / 2
        elif kwargs["x"] == "sep":
            kwargs["x"] = np.abs(data["m"] - data["n"])

        if kwargs.pop("switchxy", False):
            kwargs["x"], kwargs["y"] = kwargs["y"], kwargs["x"]

        ax, cbar = showDataContainerAsMatrix(data, v=vals, ax=ax, **kwargs)
    else:
        equidistant = kwargs.pop("equidistant", False)
        if not equidistant:
            try:
                ax, cbar = drawERTData(ax, data, vals=vals, **kwargs)
            except Exception:
                pg.warning('Something gone wrong while drawing data. '
                           'Try fallback with equidistant electrodes.')
                equidistant = True

        if equidistant:
            d = pg.DataContainerERT(data)
            sc = data.sensorCount()
            d.setSensors(list(zip(range(sc), np.zeros(sc))))
            ax, cbar = drawERTData(ax, d, vals=vals, **kwargs)

    # TODO here cbar handling like pg.show

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])

    if 'showTopo' in kwargs:
        # if axTopo is not None:
        print(ax.get_position())
        axTopo = pg.plt.axes([ax.get_position().x0,
                              ax.get_position().y0,
                              ax.get_position().x0+0.2,
                              ax.get_position().y0+0.2])

        x = pg.x(data)
        x *= (ax.get_xlim()[1] - ax.get_xlim()[0]) / (max(x)-min(x))
        x += ax.get_xlim()[0]
        axTopo.plot(x, pg.z(data), '-o', markersize=4)
        axTopo.set_ylim(min(pg.z(data)), max(pg.z(data)))
        axTopo.set_aspect(1)

    pg.viewer.mpl.updateAxes(ax)
    return ax, cbar


def drawERTData(ax, data, vals=None, **kwargs):
    """Plot ERT data as pseudosection matrix (position over separation).

    Parameters
    ----------
    data : DataContainerERT
        data container with sensorPositions and a/b/m/n fields
    vals : iterable of data.size() [data['rhoa']]
        vector containing the vals to show
    ax : mpl.axis
        axis to plot, if not given a new figure is created
    cMin/cMax : float
        minimum/maximum color vals
    logScale : bool
        logarithmic colour scale [min(A)>0]
    label : string
        colorbar label

    **kwargs:
        * dx : float
            x-width of individual rectangles
        * ind : integer iterable or IVector
            indices to limit display
        * circular : bool
            Plot in polar coordinates when plotting via patchValMap
    
    Returns
    -------
    ax:
        The used Axes
    cbar:
        The used Colorbar or None
    """
    if vals is None:
        vals = data['rhoa']

    valid = data.get("valid").array().astype("bool")
    vals = ma.array(vals, mask=~valid)

    ind = kwargs.pop('ind', None)
    sw = kwargs.pop("switch", False)
    if ind is not None:
        vals = vals[ind]
        mid, sep = midconfERT(data, ind, switch=sw)
    else:
        mid, sep = midconfERT(data, circular=kwargs.get('circular', False),
                              switch=sw)

    cbar = None
    # dx = kwargs.pop('dx', np.median(np.diff(np.unique(mid))))
    dx = kwargs.pop('dx', np.median(np.diff(pg.x(data))))
    ax, cbar, ymap = showValMapPatches(vals, xVec=mid, yVec=sep,
                                       dx=dx, ax=ax, **kwargs)

    if kwargs.get('circular', False):
        sM = np.mean(data.sensors(), axis=0)
        a = np.array([np.arctan2(s[1]-sM[1], s[0]-sM[0])
                      for s in data.sensors()])
        p = list(range(len(a)))
        # p.append(0)
        ax.plot(np.cos(a)[p], np.sin(a)[p], 'o', color='black')

        for i in range(len(a)):
            ax.text(1.15 * np.cos(a[i]),
                    1.15 * np.sin(a[i]), str(i+1),
                    horizontalalignment='center',
                    verticalalignment='center')
        ax.set_axis_off()
        ax.set_aspect(1)

    else:
        ytl = generateConfStr(np.sort([int(k) for k in ymap]), switch=sw)
        # if only DD1/WE1 in WB/SL data rename to WB/SL
        if 'DD1' in ytl and 'WB2' in ytl and 'DD2'not in ytl:
            ytl[ytl.index('DD1')] = 'WB1'
        if 'WA1' in ytl and 'SL2' in ytl and 'WA2'not in ytl:
            ytl[ytl.index('WA1')] = 'SL1'

        yt = ax.get_yticks()
        yt = np.unique(yt.clip(0, len(ytl)-1))
#        if yt[0] == yt[1]:
#            yt = yt[1:]
        dyt = np.diff(yt)
        if len(dyt) > 1 and dyt[-1] < dyt[-2]:
            yt = yt[:-1]

        ax.set_yticks(yt)
        ax.set_yticklabels([ytl[int(yti)] for yti in yt])
    return ax, cbar


def midconfERT(data, ind=None, rnum=1, circular=False, switch=False):
    """Return the midpoint and configuration key for ERT data.

    Return the midpoint and configuration key for ERT data.

    Parameters
    ----------
    data : DataContainerERT
        data container with sensorPositions and a/b/m/n fields

    ind : []
        Documentme

    rnum : []
        Documentme

    circular : bool
        Return midpoint in degree (rad) instead if meter.

    Returns
    -------
    mid : np.array of float
        representative midpoint (middle of MN, AM depending on array)
    conf : np.array of float
        configuration/array key consisting of
        1) array type (Wenner-alpha/beta, Schlumberger, PP, PD, DD, MG)
            00000: pole-pole
            10000: pole-dipole or dipole-pole
            30000: Wenner-alpha
            40000: Schlumberger or Gradient
            50000: dipole-dipole or Wenner-beta
        2) potential dipole length (in electrode spacings)
            .XX..: dipole length
        3) separation factor (current dipole length or (di)pole separation)
            ...XX: pole/dipole separation (PP,PD,DD,GR) or separation
    """
#    xe = np.hstack((pg.x(data.sensorPositions()), np.nan))  # not used anymore
    x0 = data.sensorPosition(0).x()
    xe = pg.x(data.sensorPositions()) - x0
    ux = pg.unique(xe)
    mI, mO, mT = 1, 100, 10000
    if switch:
        mI, mO = mO, mI

    if len(ux) * 2 > data.sensorCount() and not circular:  # 2D with topography case
        dx = np.array(pg.utils.diff(pg.utils.cumDist(data.sensorPositions())))
        dxM = pg.mean(dx)
        if min(pg.y(data)) != max(pg.y(data)) or \
           min(pg.z(data)) != max(pg.z(data)):
            # Topography case
            if (max(abs(dx-dxM)) < dxM*0.9):
                # if the maximum spacing < meanSpacing/2 we assume equidistant
                # spacing and no missing electrodes
                dx = np.ones(len(dx)) * dxM
            else:
                # topography with probably missing electrodes
                dx = np.floor(dx/np.round(dxM)) * dxM

        if max(dx) < 0.5:
            pg.debug("Detecting small distances, using mm accuracy")
            rnum = 3
        xe = np.hstack((0., np.cumsum(np.round(dx, rnum)), np.nan))

        de = np.median(np.diff(xe[:-1])).round(rnum)
        ne = np.round(xe/de)
    else:  # 3D (without topo) case => take positions directly
        de = np.median(np.diff(ux)).round(1)
        ne = np.array(xe/de, dtype=int)

    # a, b, m, n = data['a'], data['b'], data['m'], data['n']
    # check if xe[a]/a is better suited (has similar size)
    if circular:
        # for circle geometry
        center = np.mean(data.sensorPositions(), axis=0)
        r = data.sensors()[0].distance(center)
        s0 = data.sensors()[0]-center
        s1 = data.sensors()[1]-center
        p0 = np.arctan2(s0[1], s0[0])
        p1 = np.arctan2(s1[1], s1[0])
        if p1 > p0:
            # rotate left
            x = np.cos(np.linspace(0, 2*pi, data.sensorCount()+1)+p0)[:-1] * r
            y = np.sin(np.linspace(0, 2*pi, data.sensorCount()+1)+p0)[:-1] * r
        else:
            x = np.cos(np.linspace(2*pi, 0, data.sensorCount()+1)+p0)[:-1] * r
            y = np.sin(np.linspace(2*pi, 0, data.sensorCount()+1)+p0)[:-1] * r

        a = np.array([np.arctan2(y[i], x[i]) for i in data['a']])
        b = np.array([np.arctan2(y[i], x[i]) for i in data['b']])
        m = np.array([np.arctan2(y[i], x[i]) for i in data['m']])
        n = np.array([np.arctan2(y[i], x[i]) for i in data['n']])

        a = np.unwrap(a) % (np.pi*2)
        b = np.unwrap(b) % (np.pi*2)
        m = np.unwrap(m) % (np.pi*2)
        n = np.unwrap(n) % (np.pi*2)

    else:
        a = np.array([ne[int(i)] for i in data['a']])
        b = np.array([ne[int(i)] for i in data['b']])
        m = np.array([ne[int(i)] for i in data['m']])
        n = np.array([ne[int(i)] for i in data['n']])

    if ind is not None:
        a = a[ind]
        b = b[ind]
        m = m[ind]
        n = n[ind]

    anan = np.isnan(a)
    a[anan] = b[anan]
    b[anan] = np.nan
    ab, am, an = np.abs(a-b), np.abs(a-m), np.abs(a-n)
    bm, bn, mn = np.abs(b-m), np.abs(b-n), np.abs(m-n)

    if circular:
        for v in [ab, mn, bm, an]:
            v[v > pi] = 2*pi - v[v > pi]

    # 2-point (default) 00000
    sep = np.abs(a-m)  # * mI # does not make sense here
    mid = (a+m) / 2

    # 3-point (PD, DP) (now only b==-1 or n==-<1, check also for a and m)
    imn = np.isfinite(n)*np.isnan(b)
    mid[imn] = (m[imn]+n[imn]) / 2
    sep[imn] = np.minimum(am[imn], an[imn]) * mI + mT + mO * (mn[imn]-1) + \
        (np.sign(a[imn]-m[imn])/2+0.5) * mT
    iab = np.isfinite(b)*np.isnan(n)
    mid[iab] = (a[iab]+b[iab]) / 2  # better 20000 or -10000?
    sep[iab] = np.minimum(am[iab], bm[iab]) * mI + mT + mO * (ab[iab]-1) + \
        (np.sign(a[iab]-n[iab])/2+0.5) * mT
    #  + 10000*(a-m)

    # 4-point alpha: 30000 (WE) or 4000 (SL)
    iabmn = np.isfinite(a) & np.isfinite(b) & np.isfinite(m) & np.isfinite(n)
    ialfa = np.copy(iabmn)
    ialfa[iabmn] = (ab[iabmn] >= mn[iabmn]+2)  # old
    mnmid = (m[iabmn] + n[iabmn]) / 2
    ialfa[iabmn] = np.sign((a[iabmn]-mnmid)*(b[iabmn]-mnmid)) < 0

    mid[ialfa] = (m[ialfa] + n[ialfa]) / 2
    spac = np.minimum(bn[ialfa], bm[ialfa])
    abmn3 = np.round((3*mn[ialfa]-ab[ialfa])*mT)/mT
    sep[ialfa] = spac * mI + (mn[ialfa]-1) * mO * (abmn3 != 0) + \
        3*mT + (abmn3 < 0)*mT
    # gradient

    # 4-point beta
    ibeta = np.copy(iabmn)
    ibeta[iabmn] = (bm[iabmn] >= mn[iabmn]) & (~ialfa[iabmn])

    if circular:
        # print(ab[ibeta])
        ibeta = np.copy(iabmn)

        def _averageAngle(vs):
            sumsin = 0
            sumcos = 0

            for v in vs:
                sumsin += np.sin(v)
                sumcos += np.cos(v)

            return np.arctan2(sumsin, sumcos)

        abC = _averageAngle([a[ibeta], b[ibeta]])
        mnC = _averageAngle([m[ibeta], n[ibeta]])

        mid[ibeta] = _averageAngle([abC, mnC])

        # special case when dipoles are completely opposite
        iOpp = abs(abs((mnC - abC)) - np.pi) < 1e-3
        mid[iOpp] = _averageAngle([b[iOpp], m[iOpp]])

        minAb = min(ab[ibeta])
        sep[ibeta] = 5 * mT + (np.round(ab[ibeta]/minAb)) * mO + \
            np.round(np.minimum(np.minimum(am[ibeta], an[ibeta]),
                                np.minimum(bm[ibeta], bn[ibeta])) / minAb) * mI
    else:
        mid[ibeta] = (a[ibeta] + b[ibeta] + m[ibeta] + n[ibeta]) / 4

        sep[ibeta] = 5 * mT + (ab[ibeta]-1) * mO + np.minimum(
            np.minimum(am[ibeta], an[ibeta]),
            np.minimum(bm[ibeta], bn[ibeta])) * mI

    # %% 4-point gamma
    # multiply with electrode distance and add first position
    if not circular:
        mid *= de
        mid += x0

    return mid, sep


def generateConfStr(yy, switch=False):
    """Generate configuration string to characterize array."""
    mI, mO, mT = 1, 100, 10000
    types = ['PP', 'PD', 'DP', 'WA', 'SL', 'DD']  # base types
    typ = np.round(yy//mT)
    if switch:
        dip = yy % mO  # source-receiver distance
        spac = np.round(yy//mO) % mO  # MN dipole length
    else:
        spac = yy % mO  # source-receiver distance
        dip = np.round(yy//mO) % mO  # MN dipole length

    # check if SL is actually GR (multi-gradient)

    # check if DD-n-n should be renamed
    rendd = (np.mean(spac / (dip+1)) < 2.1)
    keys = []
    for s, d, t in zip(spac, dip, typ):
        key = types[t]
        if d > 0:
            if rendd and d+1 == s and t == 5:
                key = 'WB'
            else:
                key = key + str(d+1) + '-'
        key = key + "{:2d}".format(s)  # str(s)
        keys.append(key)

    return keys
