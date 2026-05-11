#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create output for example-2.

Use it with the synthetic dataset:

python example-2-visualisation.py synth/

or with the final inversion result:

python example-2-visualisation.py permModel_h-3/
"""

import sys

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.viewer.mpl import saveAxes

from matplotlib.offsetbox import AnchoredText

days = [0, 2, 4, 6, 8, 10]


def add_inner_title(ax, title, loc=2, color="k", **kwargs):
    at = AnchoredText(title, loc=loc,
                      pad=0., borderpad=0.5, prop={"color": color},
                      frameon=False, **kwargs)
    ax.add_artist(at)
    return at


def savefig(mesh, geom, data=None, label='', out=None, **kwargs):
    """Little shortcut to plot mesh with associated geometry."""

    pg.viewer.mpl.hold(1)
    ax = kwargs.pop('ax', None)
    print(ax)

    if data is not None:
        ax, _ = pg.show(mesh, data, hold=1, pad=0.55, label=label,
                        ax=ax, **kwargs)
    if out:
        ax, _ = pg.show(geom, ax=ax, fillRegion=False, **kwargs)
        saveAxes(ax, out, adjust=True)
    else:
        ax, _ = pg.show(geom, fillRegion=False, ax=ax, **kwargs)

    if kwargs.pop('adjust', True):
        pg.viewer.mpl.adjustWorldAxes(ax)
    return ax


def showSynthData(synthPath):
    geom = pg.load(synthPath + 'synthGeom.bms')
    mesh = pg.load(synthPath + 'synth.bms')
    k = np.load(synthPath + 'synthK.npy')
    vel = np.load(synthPath + 'synthVel.npy')
    sat = np.load(synthPath + 'synthSat.npy')

    scheme = pg.DataContainer(synthPath + 'synth.shm', 'a b m n')
#    rhoaR = np.load(synthPath + 'synthRhoaRatio.npy')
#    rhoa = np.load(synthPath + 'synthRhoa.npy')

    row = 3
    col = 2

    # START model perm + input
    savefig(mesh, geom, k, 'Hydraulic conductivity $K$ in m$/$s',
            out='hydrConductModel',
            cMin=1e-5, cMax=1e-2, nLevs=4, cmap='viridis')

    # START velocity
    axVel, _ = pg.show(mesh, np.sqrt(vel[0]**2 + vel[1]**2),
                       logScale=0, colorBar=1, pad=0.55,
                       label='Velocity $|v|$ in m$/$s', hold=1)

    meshC = pg.meshtools.createMesh(geom, quality=33, area=0.5, smooth=[1, 10])
    pg.show(mesh, data=vel, ax=axVel, coarseMesh=meshC, color='black',
            linewidth=0.5, dropTol=1e-6)
    pg.show(geom, ax=axVel, fillRegion=False)
    saveAxes(axVel, 'hydrVelocity', adjust=True)

    # START Saturation
    axs = plt.subplots(row, col, sharex=True, sharey=True,
                       figsize=(10.*0.65, 7.25*0.65))[1].flatten()

    satScale = 0.001
    for i, a in enumerate(axs):
        savefig(mesh, geom,
                sat[i*len(sat)/(len(axs))+1] * satScale,  # /mesh.cellSizes(),
                label=None,
                out=None,
                cMin=0, cMax=2.5,
                ax=a,
                adjust=True)
        pg.viewer.mpl.drawSensors(a, scheme.sensorPositions(),
                                 diam=0.15, color='green')

        add_inner_title(a, "t = %d days" % days[i], 3, color="w")

        if i < (row-1)*col:
            a.set_xlabel('')
        if i % col:
            a.set_ylabel('')
        a.set_ylim([-16, 0])

    pg.viewer.mpl.saveFigure(axs[0].figure, "hydrSaturation")
    pg.viewer.mpl.createColorBarOnly(cMin=0, cMax=2.5, logScale=False,
                                    cMap='Spectral_r',
                                    nLevs=5,
                                    label=r'Concentration $c$ in g$/$l',
                                    orientation='horizontal',
                                    savefig='hydrSaturationCbar')

    # END Saturation
    pg.wait()


def showModel(outPath):

    # geom = pg.load('synth/' + 'synthGeom.bms')

    # pd = pg.load(outPath + '/paraDomain.bms')
    paraMesh = pg.load(outPath + '/paraMesh.bms')
    model = pg.load(outPath + "/model.vector")

    fopMesh = pg.load(outPath + "fopMesh.bms")
    fopModel = pg.load(outPath + "fopModel.vector")

    allModel = np.zeros(len(model)+2)
    allModel[2:] = model

#    allModel[0] = fopModel[fopMesh.findCellByMarker(
#        pg.MARKER_FIXEDVALUE_REGION - 0)[0].id()]
    allModel[1] = fopModel[fopMesh.findCellByMarker(
        pg.MARKER_FIXEDVALUE_REGION - 1)[0].id()]

    ax = savefig(paraMesh, None, allModel[paraMesh.cellMarkers()],
                 'Hydraulic conductivity $K$ in m$/$s',
                 cMin=1e-5, cMax=1e-2,
                 nLevs=4, cMap='viridis')
    pg.wait()
    pg.show(fopMesh, ax=ax, linewidth=0.2)

    paraMesh.createNeighborInfos()
    bs = []
    for b in paraMesh.boundaries():
        if b.leftCell() and b.rightCell():
            if b.leftCell().marker() > 1 or b.rightCell().marker() > 1:
                bs.append(b)
    pg.viewer.mpl.drawSelectedMeshBoundaries(
            ax, bs, color=(0.0, 0.0, 0.0, 1.0), linewidth=1.0)

    pg.viewer.mpl.saveFigure(ax.figure, 'hydrInversionModel')

    print('Done!')
    pg.wait()

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Give result path!")
        sys.exit()

    outPath = sys.argv[1]
    print(outPath)

    if 'synth' in outPath:
        showSynthData(outPath)
        # pg.wait()
    else:
        showModel(outPath)
