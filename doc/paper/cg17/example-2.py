#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example for fully coupled hydrogeophysical Inversion."""

import os

import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import ert
from pygimli.physics.petro import resistivityArchie


class WorkSpace(object):
    """Empty class to store some data."""
    pass


def solveDarcy(mesh, k=None, p0=1, verbose=False):
    """Darcy flow."""
    if verbose:
        print("Solve Darcy equation ...")

    uDir = {2: p0,  # left aquifer
            3: p0,  # left bedrock
            # [4, 0],  # bottom (paper)
            5: 0,  # right bedrock
            6: 0,  # right aquiver
            7: 0,  # right top
            }

    p = pg.solver.solve(mesh, a=k, bc={'Dirichlet': uDir}, verbose=True)
    vel = -pg.solver.grad(mesh, p) * np.asarray([k, k, k]).T
    mvel = mt.cellDataToNodeData(mesh, vel)
    return mesh, mvel, p, k, np.asarray([pg.x(vel), pg.y(vel)])


def solveAdvection(mesh, vel, times, diffusion, verbose=False):
    """Solve Diffusion/Advection equation"""
    if verbose:
        print("Solve for concentration movement on", len(times),
              "time steps ...")

    S = pg.Vector(mesh.cellCount(), 0.0)
    injectPos = [-19.1, -4.6]
    sourceCell = mesh.findCell(injectPos)
    S[sourceCell.id()] = 1.0/sourceCell.size()

    t = times[:int(len(times)/2)]
    t = np.linspace(t[0], t[-1], len(t))

    c1 = pg.solver.solveFiniteVolume(mesh, a=diffusion, f=S, vel=vel,
                                     times=t,  # uB=[1, 0],
                                     bc={'Dirichlet': {1: 0}},
                                     scheme='PS', verbose=0)

    c2 = pg.solver.solveFiniteVolume(mesh, a=diffusion, f=0., vel=vel,
                                     times=t,  # uB=[1, 0],
                                     bc={'Dirichlet': {1: 0}},
                                     u0=c1[-1],
                                     scheme='PS', verbose=0)
    c = np.vstack((c1, c2))

    return c[::]


def solveERT(mesh, concentration, verbose=0):
    """Simulate resistivity distribution for given nonsteady concentration."""
    if verbose:
        print("Solve for ERT ...")

    ertScheme = ert.createData(pg.utils.grange(-20, 20, dx=1.0),
                               schemeName='dd')

    meshERT = mt.createParaMesh(ertScheme, quality=33, paraMaxCellSize=0.2,
                                boundaryMaxCellSize=50, smooth=[1, 2])

    scale = 0.001
    concentration *= scale  # mg/m²

    # apply saturation model to simulate unsaturated topsoil
    sat = np.zeros(mesh.cellCount())
    for c in mesh.cells():
        if c.center()[1] < -8:
            sat[c.id()] = 1.
        elif c.center()[1] < -2:
            sat[c.id()] = 1.
        else:
            sat[c.id()] = .5

    cWater = 1./100.
    conductivity = concentration * 0.1 + cWater

    rArchie = resistivityArchie(rFluid=1. / conductivity,
                                porosity=0.3, sat=sat, m=1.3,
                                mesh=mesh, meshI=meshERT, fill=1)

    # apply background resistivity model
    rho0 = np.zeros(meshERT.cellCount())
    for c in meshERT.cells():
        if c.center()[1] < -8:
            rho0[c.id()] = 150.
        elif c.center()[1] < -2:
            rho0[c.id()] = 500.
        else:
            rho0[c.id()] = 1000.

    resis = pg.Matrix(rArchie)

    for i, rbI in enumerate(rArchie):
        resis[i] = 1. / ((1./rbI) + 1./rho0)

    ertScheme['k'] = ert.createGeometricFactors(ertScheme)

    errPerc = 0.01
    errVolt = 1e-5

    rhoa = ert.simulate(mesh=meshERT, scheme=ertScheme, res=resis,
                        verbose=0, returnArray=True)

    voltage = rhoa / ertScheme('k')
    err = np.abs(errVolt / voltage) + errPerc

    dRhoa = rhoa[1:] / rhoa[0]
    dErr = err[1:]

    return meshERT, ertScheme, resis, rhoa, dRhoa, dErr


class HydroGeophysicalModelling(pg.Modelling):
    """Forward Operator for fully coupled hydrogeophysical inversion."""

    def __init__(self, verbose=False, **kwargs):
        """Constructor."""
        pg.core.ModellingBase.__init__(self, verbose=verbose)
        self.init(kwargs.pop('mesh', None),
                  kwargs.pop('tMax', 50000),
                  kwargs.pop('satSteps', 50),
                  kwargs.pop('ertSteps', 5))
        self.iter = 0

    def init(self, mesh, tMax, satSteps, ertSteps):
        """Initialize some settings."""

        if mesh is not None:
            self.parMesh = pg.Mesh(mesh)
            self.setMesh(mesh)
            self.createRefinedForwardMesh(refine=False, pRefine=False)

        self.tMax = tMax
        self.satSteps = satSteps
        self.ertSteps = ertSteps
        self.timesAdvection = np.linspace(1, tMax, satSteps)
        self.timesERT = pg.core.IndexArray(np.floor(
            np.linspace(0, len(self.timesAdvection)-1, self.ertSteps)))

        self._J = pg.Matrix()
        self.setJacobian(self._J)
        self.ws = WorkSpace()

    def response(self, par):
        """Return response.

        Simulate resistivity data for a given hydraulic conductivity.
        """
        return self.response_mt(par)

    def response_mt(self, par, i=0):
        """Return response (multi threaded)."""

        verbose = 1
        if i == 0:
            ws = self.ws
        else:
            ws = WorkSpace()

        mesh = pg.Mesh(self.mesh())

        k = self.createMappedModel(par)

        ws.mesh, ws.vel, ws.p, ws.k, ws.velC = solveDarcy(mesh, k=k, p0=0.75,
                                                          verbose=verbose)

        ws.sat = solveAdvection(ws.mesh, ws.vel, self.timesAdvection,
                                diffusion=pg.abs(ws.velC.T) * 1e-2,
                                verbose=verbose)

        ws.meshERT, ws.scheme, ws.resis, ws.rhoa, ws.rhoaR, ws.derr = \
            solveERT(ws.mesh, ws.sat[self.timesERT], verbose=verbose)

        return ws.rhoaR.flatten()


def simulateSynth(model, tMax=5000, satSteps=150, ertSteps=10, area=0.1,
                  synthPath='synth/'):
    """Create synthetic example."""

    if not os.path.exists('synth/'):
        os.mkdir(synthPath)

    world = mt.createWorld(start=[-20, 0], end=[20, -16], layers=[-2, -8],
                           worldMarker=False)
    for i, b in enumerate(world.boundaries()):
        b.setMarker(i + 1)

    block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0], marker=4,
                               boundaryMarker=11, area=area)
    geom = world + block
    geom.save(synthPath + 'synthGeom')
    # pg.show(geom, boundaryMarker=1)

    paraMesh = pg.meshtools.createMesh(geom, quality=32, area=area,
                                       smooth=[1, 10])

    # translate 1 2 3 4 - > 0 1 2 3
    mapMarker = np.array([0, 0, 1, 2, 3], 'float')
    paraMesh.setCellMarkers(mapMarker[np.array(paraMesh.cellMarkers())])
    paraMesh.save(synthPath + 'synth.bms')

    fop = HydroGeophysicalModelling(mesh=paraMesh, tMax=tMax,
                                    satSteps=satSteps,
                                    ertSteps=ertSteps,
                                    verbose=1)

    # openblas have some problems with to high thread count ..
    # we need to dig into
    print("ThreadCount:", pg.core.threadCount())
    pg.setThreadCount(4)

    print('##### Simulate synthetic data ' + '#'*50)
    pg.tic()
    rhoaR = fop.response(pg.Vector(model)[paraMesh.cellMarkers()])
    pg.toc()
    print('#'*100)

    # add some noise here
    rand = pg.Vector(len(rhoaR))
    pg.math.randn(rand)

    rhoaR *= (1.0 + rand * fop.ws.derr.flatten())
    fop.ws.rhoaR = rhoaR.reshape(fop.ws.derr.shape)

    # fop.ws.mesh.save(synthPath + 'synth.bms')
    np.save(synthPath + 'synthK', fop.ws.k)
    np.save(synthPath + 'synthVel', fop.ws.vel)
    np.save(synthPath + 'synthSat', fop.ws.sat)

    fop.ws.scheme.save(synthPath + 'synth.shm', 'a b m n')
    np.save(synthPath + 'synthRhoaRatio', fop.ws.rhoaR)
    np.save(synthPath + 'synthRhoa', fop.ws.rhoa)
    np.save(synthPath + 'synthErr', fop.ws.derr)


def createFopWithParaDomain(paraRefine=0, ncpu=6):
    """Create Forward operator and synthetic data."""
    tMax = 345600 * 3
    satSteps = 800 * 2
    ertSteps = 10

    synthPath = 'synth/'

    synthArea = 0.1
    paraArea = 0.1

    simulateSynth(model=[1e-8, 5e-3, 1e-4, 8e-4],
                  tMax=tMax,
                  satSteps=satSteps,
                  ertSteps=ertSteps,
                  area=synthArea,
                  synthPath=synthPath)

    fop = HydroGeophysicalModelling(mesh=None, tMax=tMax,
                                    satSteps=satSteps,
                                    ertSteps=ertSteps,
                                    verbose=1)

    rhoaR = np.load(synthPath + 'synthRhoaRatio.npy')
    err = np.load(synthPath + 'synthErr.npy')

    fop.setVerbose(True)
    fop.setMultiThreadJacobian(ncpu)

    paraMesh = pg.createGrid(x=[-20, -10, 0, 10, 20], y=[-16, -8, -2, 0])
    # top Boundary Marker
    for i, b in enumerate(paraMesh.findBoundaryByMarker(4)):
        b.setMarker(8)
    # right Boundary Marker
    for i, b in enumerate(paraMesh.findBoundaryByMarker(2)):
        b.setMarker(7 - i)
    # bottom Boundary Marker
    for i, b in enumerate(paraMesh.findBoundaryByMarker(3)):
        b.setMarker(4)
    # left boundary Marker
    for i, b in enumerate(paraMesh.findBoundaryByMarker(1)):
        b.setMarker(3-i)

#    ax, _ = pg.show(paraMesh)
#    pg.viewer.mpl.drawMeshBoundaries(ax=ax, mesh=paraMesh)
#    pg.wait()

    paraMesh.cell(0).setMarker(0)  # bottom
    paraMesh.cell(1).setMarker(0)  # bottom
    paraMesh.cell(2).setMarker(0)  # bottom
    paraMesh.cell(3).setMarker(0)  # bottom
    paraMesh.cell(4).setMarker(2)  # center
    paraMesh.cell(5).setMarker(2)  # center
    paraMesh.cell(6).setMarker(2)  # center
    paraMesh.cell(7).setMarker(2)  # center
    paraMesh.cell(8).setMarker(1)  # top
    paraMesh.cell(9).setMarker(1)  # top
    paraMesh.cell(10).setMarker(1)  # top
    paraMesh.cell(11).setMarker(1)  # top

    if paraRefine:
        for i in range(paraRefine):
            paraMesh = paraMesh.createH2()

        paraCount = 2
        for c in paraMesh.cells():
            if c.marker() > 1:
                c.setMarker(paraCount)
                paraCount += 1

        fop.setMesh(paraMesh)

        # for i in range(fop.regionManager().regionCount()):
        #     fop.regionManager().region(i).setSingle(1)
    else:
        fop.setMesh(paraMesh)
        fop.createRefinedForwardMesh(refine=False, pRefine=False)

    fop.setRegionProperties("*", single=True)
    fopMesh = pg.meshtools.createMesh(fop.regionManager().paraDomain(),
                                      area=paraArea, smooth=[1, 10])

    fop.setMesh(fopMesh)  # , ignoreRegionManager=False)
    fop.setVerbose(True)

    return fop, rhoaR, err, paraMesh


if __name__ == '__main__':

    paraRefine = 2
    ncpu = 6

    fop, rhoaR, err, paraMesh = createFopWithParaDomain(paraRefine=paraRefine,
                                                        ncpu=ncpu)
    fop.setRegionProperties(1, fix=1e-8)
    fop.setRegionProperties(0, fix=1e-4)
    # fop.regionManager().region(1).setFixValue(1e-8)  # top layer
    # fop.regionManager().region(0).setFixValue(1e-4)  # bedrock
    # fop.regionManager().region(0).setSingle(1)  # bedrock

    # Reflect the fix value setting here!!!!
    fop.createRefinedForwardMesh(refine=False, pRefine=False)

    # Connect all regions
    for i in range(2, fop.regionManager().regionCount()):
        for j in range(i+1, fop.regionManager().regionCount()):
            fop.regionManager().setInterRegionConstraint(i, j, 1.0)

    startModel = pg.Vector(fop.regionManager().parameterCount(), 1e-3)

    fop.setStartModel(startModel)

    inv = pg.Inversion(rhoaR.flatten(), fop, verbose=1, dosave=0)

    tD = pg.trans.TransLog()
    tM = pg.trans.TransLogLU(1e-9, 1e-2)
    inv.setTransData(tD)
    inv.setTransModel(tM)

    inv.setRelativeError(err.flatten())
    inv.setMaxIter(50)
    inv.setLineSearch(True)
    inv.setLambda(1000)

    outPath = "permModel_h-" + str(paraRefine)
    if not os.path.exists(outPath):
        os.mkdir(outPath)

    paraMesh.save(outPath + '/paraMesh')
    fop.mesh().save(outPath + "/fopMesh")
    fop.regionManager().paraDomain().save(outPath + '/paraDomain')

    coeff = inv.start()
    coeff.save(outPath + "/model.vector")

    allModel = fop.createMappedModel(coeff)
    allModel.save(outPath + "/fopModel.vector")
