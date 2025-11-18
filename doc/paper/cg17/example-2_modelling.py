#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Extended minimal pygimli example to simulate Darcy velocity,
    mass transport and time-lapse ERT measurements
"""
import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert
import pygimli.physics.petro as petro


# Create geometry definition for the modelling domain
world = mt.createWorld(start=[-20, 0], end=[20, -16], layers=[-2, -8],
                       worldMarker=False)
# Create a heterogeneous block
block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
                           marker=4, boundaryMarker=10, area=0.1)
# Merge geometrical entities
geom = world + block
# pg.show(geom, boundaryMarker=True, savefig='geometry.pdf')

# Create a mesh from the geometry definition
mesh = mt.createMesh(geom, quality=32, area=0.2, smooth=[1, 10])
# pg.show(mesh, savefig='mesh.pdf')

print('Solve Darcy equation ... ')
# Map regions to hydraulic conductivity in $m/s$
kMap = [[1, 1e-8], [2, 5e-3], [3, 1e-4], [4, 8e-4]]
# Map conductivity value per region to each cell in the given mesh
K = pg.solver.parseMapToCellArray(kMap, mesh)
# Dirichlet conditions for hydraulic potential
# pBound = [[[1, 2, 3], 0.75], [[5, 6, 7], 0.0]]
pBound = {1: 0.75, 2: 0.75, 3: 0.75, 5: 0, 6: 0, 7: 0}
# Solve for hydraulic potential
p = pg.solver.solveFiniteElements(mesh, a=K, bc={'Dirichlet': pBound})
# Solve velocity as gradient of hydraulic potential
vel = -pg.solver.grad(mesh, p) * np.asarray([K, K, K]).T

ax, _ = pg.show(mesh, data=K, label='Hydraulic conductivity $K$ in m$/$s',
                cMin=1e-5, cMax=1e-2, nLevs=4, cMap='viridis')
ax, _ = pg.show(mesh, data=pg.abs(vel), logScale=0,
                label='Velocity $v$ in m$/$s')
ax, _ = pg.show(mesh, data=vel, ax=ax, color='black', linewidth=0.5,
                dropTol=1e-6)

print('Solve Advection-diffusion equation ...')
S = pg.Vector(mesh.cellCount(), 0.0)
# Fill injection source vector for a fixed injection position
sourceCell = mesh.findCell([-19.1, -4.6])
S[sourceCell.id()] = 1.0 / sourceCell.size() / 1000  # g/(l s)
# Choose 800 time steps for 6 days in seconds
# t = pg.utils.grange(0, 6 * 24 * 3600, n=800)
t = np.arange(0, 6*24*3600+1, 600)
# Create dispersitivity, depending on the absolute velocity
dispersion = pg.abs(vel) * 1e-2
# Solve for injection time, but we need velocities on cell nodes
vel = mt.cellDataToNodeData(mesh, vel)
c1 = pg.solver.solveFiniteVolume(mesh, a=dispersion, f=S, vel=vel,
                                 times=t,
                                 bc={'Dirichlet': {1: 0}},
                                 scheme='PS', verbose=0)
# Solve without injection starting with last result
c2 = pg.solver.solveFiniteVolume(mesh, a=dispersion, f=0, vel=vel,
                                 u0=c1[-1], times=t,
                                 bc={'Dirichlet': {1: 0}},
                                 scheme='PS', verbose=0)
# Stack results together
c = np.vstack((c1, c2))
# Select 10 time frame to simulate ERT data
# timesERT = np.array(np.linspace(0, len(c)-1, 10), dtype=int)
timesERT = np.arange(0, 11, 2) * 6 * 24
# %%
axs = pg.plt.subplots(3, 2, sharex=True, sharey=True)[1].flatten()
for i in range(6):
    pg.show(mesh, c[timesERT[i]], cMin=0, cMax=2.5, ax=axs[i],
            label='Concentration $c$ in g$/$l', orientation="vertical")
# %%
print('Solve ERT modelling ...')

# Create survey measurement scheme
ertScheme = ert.createData(pg.utils.grange(-20, 20, dx=1.0),
                           schemeName='dd')
# Create suitable mesh for ert forward calculation
meshERT = mt.createParaMesh(ertScheme, quality=33, paraMaxCellSize=0.2,
                            boundaryMaxCellSize=50, smooth=[1, 2])
# %%Create conductivity of fluid for salt concentration $c$
# sigmaFluid = c[timesERT] / 11.23 * 0.001 + 0.01
sigmaFluid = c[timesERT] * 0.1 + 0.01
# Calculate bulk resistivity based on Archie's Law
resBulk = petro.resistivityArchie(rFluid=1./sigmaFluid,
                                  porosity=0.3, m=1.3,
                                  mesh=mesh, meshI=meshERT, fill=1)
# %% apply background resistivity model
rho0 = np.zeros(meshERT.cellCount()) + 1000.
for cell in meshERT.cells():
    if cell.center()[1] < -8:
        rho0[cell.id()] = 150.
    elif cell.center()[1] < -2:
        rho0[cell.id()] = 500.

resis = pg.Matrix(resBulk)
for i, rbI in enumerate(resBulk):
    resis[i] = 1. / ((1./rbI) + 1./rho0)
# %%Run  simulation for  the apparent resistivities
rhoa = ert.simulate(meshERT, scheme=ertScheme, res=resis, verbose=0,
                    returnArray=True)
# %%Solve the electrical forward problem using the ERT method manager
axs = pg.plt.subplots(3, 2, sharex=True, sharey=True)[1].flatten()
for i in range(6):
    ert.showData(ertScheme, rhoa[i], ax=axs[i], cMin=20, cMax=500,
                 colorBar=1, cMap="Spectral_r", orientation="vertical")

# just hold figure windows open if run outside from spyder, ipython or similar
pg.wait()
