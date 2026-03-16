#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plotting functions for magnetic resonance data."""

import numpy as np
from pygimli.viewer.mpl import drawModel1D

# def drawModel1D(ax, thickness, values, plotfunction='plot',
#                 xlabel='', *args, **kwargs):
#  """Draw 1d block model into axis ax defined by values and thickness vectors
#     using plotfunction."""
#
#     nLayers = len(thickness) + 1
#     px = np.zeros(nLayers * 2)
#     pz = np.zeros(nLayers * 2)
#     z1 = np.cumsum(thickness)
#
#     for i in range(nLayers):
#         px[2 * i] = values[i]
#         px[2 * i + 1] = values[i]
#
#         if i == nLayers - 1:
#             pz[2 * i + 1] = z1[i - 1] * 1.2
#         else:
#             pz[2 * i + 1] = z1[i]
#             pz[2 * i + 2] = z1[i]
#
#     if plotfunction == 'loglog' or plotfunction == 'semilogy':
#         pz[0] = thickness[0] * 0.8
#
#     try:
#         plot = getattr(ax, plotfunction)
#         plot(px, pz, *args, **kwargs)
#     except Exception as e:
#         print(e)
#
#     ax.set_ylabel('Depth (m)')
#     ax.set_xlabel(xlabel)
#     ax.set_ylim(pz[-1], pz[0])
#     ax.grid(True)


def showErrorBars(ax, thk, val, thkL, thkU, valL, valU, *args, **kwargs):
    """Plot wc and t2 models with error bars."""
    zb = np.cumsum(thk)
    zm = np.hstack((zb - thk / 2, zb[-1] * 1.2))  # zb[-1]+thk[-1]/2))
    valm = (val[:-1] + val[1:]) / 2
    xerr = [val - valL, valU - val]
    yerr = [thk - thkL, thkU - thk]
    ax.errorbar(val, zm, fmt='.', xerr=xerr, ecolor='r', **kwargs)
    ax.errorbar(valm, zb, fmt='.', yerr=yerr, ecolor='g', **kwargs)
    ax.set_ylim(bottom=zm[-1] * 1.02, top=0)


def showWC(ax, thk, wc, wmin=0., wmax=0.45, maxdep=0., dw=0.05, **kwargs):
    """Show water content function nicely."""
    drawModel1D(ax, thk, wc, xlabel=r'$\theta$')
    ax.set_xlim(0., 0.45)
    if maxdep > 0.:
        ax.set_ylim(maxdep, 0.)
    wt = np.arange(wmin, wmax, dw)
    ax.set_xticks(wt)
    ax.set_xticklabels([str(wi) for wi in wt])


def showT2(ax, thk, t2, maxdep=0., **kwargs):
    """Show T2 function nicely."""
    drawModel1D(ax, thk, t2*1e3, xlabel=r'$T_2^*$ [ms]',
                plot='semilogx')
    tmin = min(20, min(t2) * 0.9e3)
    tmax = max(500, max(t2) * 1.1e3)
    ax.set_xlim(tmin, tmax)
    if maxdep > 0.:
        ax.set_ylim(maxdep, 0.)
    xt = [20, 50, 100, 200, 500]
    ax.set_xticks(xt)
    ax.set_xticklabels([str(ai) for ai in xt])

if __name__ == "__main__":
    pass
