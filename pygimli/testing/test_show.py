#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np

import pygimli as pg
from pygimli import plt

import pygimli.meshtools as mt


def testShowVariants():
    # Create geometry definition for the modelling domain
    world = mt.createWorld(start=[-10, 0], end=[10, -16],
                           layers=[-8], worldMarker=False)

    # Create a heterogeneous block
    block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
                               marker=4, boundaryMarker=10, area=0.1)

    circ = mt.createCircle(pos=[0, -11], radius=2, boundaryMarker=11,
                           isHole=True)

    # Merge geometrical entities
    geom = world + block + circ
    mesh = mt.createMesh(geom)

    fig, axs = plt.subplots(3, 5)

    pg.show(geom, ax=axs[0][0])
    axs[0][0].set_title('plc, (default)')
    pg.show(geom, fillRegion=False, ax=axs[0][1])
    axs[0][1].set_title('plc, fillRegion=False')
    pg.show(geom, showBoundary=False, ax=axs[0][2])
    axs[0][2].set_title('plc, showBoundary=False')
    pg.show(geom, markers=True, ax=axs[0][3])
    axs[0][3].set_title('plc, markers=True')
    pg.show(mesh, ax=axs[0][4], showBoundary=False)
    axs[0][4].set_title('mesh, showBoundary=False')

    pg.show(mesh, ax=axs[1][0])
    axs[1][0].set_title('mesh, (default)')
    pg.show(mesh, mesh.cellMarkers(), label='Cell markers', ax=axs[1][1])
    axs[1][1].set_title('mesh, cells, (default)')
    pg.show(mesh, markers=True, ax=axs[1][2])
    axs[1][2].set_title('mesh, cells, markers=True')
    pg.show(mesh, mesh.cellMarkers(), label='Cell markers', showMesh=True,
            ax=axs[1][3])
    axs[1][3].set_title('mesh, cells, showMesh=True')
    pg.show(mesh, mesh.cellMarkers(), label='Cell markers', showBoundary=False,
            ax=axs[1][4])
    axs[1][4].set_title('mesh, cells, showBoundary=False')

    pg.show(mesh, pg.x(mesh), label='Nodes (x)', ax=axs[2][0])
    axs[2][0].set_title('mesh, nodes, (default)')
    pg.show(mesh, pg.x(mesh), label='Nodes (x)', showMesh=True, ax=axs[2][1])
    axs[2][1].set_title('mesh, nodes, showMesh=True')
    pg.show(mesh, pg.x(mesh), label='Nodes (x)', showBoundary=True,
            ax=axs[2][2])
    axs[2][2].set_title('mesh, nodes, showBoundary=True')
    pg.show(mesh, pg.y(mesh.cellCenters()), label='Cell center (y)',
            tri=True, shading='flat', ax=axs[2][3])
    axs[2][3].set_title('mesh, cells, tri=True, shading=flat')
    pg.show(mesh, pg.y(mesh.cellCenters()), label='Cell center (y)',
            tri=True, shading='gouraud', ax=axs[2][4])
    axs[2][4].set_title('mesh, cells, tri=True, shading=gouraud')
    # pg.show(mesh, mesh.cellMarker(), label(markers), axs[1][1])
    axs[2][4].figure.tight_layout()
    fig.tight_layout()


def testColorbar():

    grid = pg.createGrid(x=np.linspace(10., 110., 11)-5,
                         y=np.linspace(0., 20, 2))

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(10, 6))
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()),
                       label='log x', ax=axs[0][0], showMesh=True,
                       cMap='Paired', logScale=True)
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())+5.,
                       label='log x', ax=axs[1][0], showMesh=True,
                       cMap='Paired', logScale=True)

    pg.viewer.mpl.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))

    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()), logScale=True,
                       ax=axs[2][0], showMesh=True, cMap='Paired')
    ###########################################################################
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()),
                       label='lin x, cMap=Paired',
                       ax=axs[0][1], logScale=False, cMap='Paired')
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())+5,
                       label='lin x, cMap=Paired (same like above)',
                       ax=axs[1][1], logScale=False, cMap='Paired')

    ##!!// thisfails and changes cmap
    pg.viewer.mpl.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))

    ###########################################################################
    grid = pg.createGrid(x=np.linspace(10., 110., 11)-25,
                         y=np.linspace(0., 20, 2))
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter()),
                       label='log with neg. x', ax=axs[0][2], showMesh=True,
                       cMap='Paired', logScale=True)
    ax, cbar = pg.show(grid, data=pg.x(grid.cellCenter())-45,
                       label='log with neg. x', ax=axs[1][2], showMesh=True,
                       cMap='Paired', logScale=True)
    pg.viewer.mpl.setMappableData(cbar.mappable, pg.x(grid.cellCenter()))
    ax.figure.tight_layout()

    ax, cbar = pg.show(grid, data=pg.x(grid.positions())-5,
                       label='white at 0 in the middle', ax=axs[2][2], showMesh=True,
                       cMap='RdBu', #cMin= -80, cMax= 80,
                       nLevs=6, logScale=False)

    return
    pg.show(grid, pg.x(grid.cellCenter()), tri=True, shading='gouraud',
            cMap='Spectral_r', logScale=False, cMin=0.01, cMax=10,
            levels=[10, 55, 100],
            orientation="vertical",
            colorBar=True)

    pg.show(grid, pg.x(grid.cellCenter()), tri=True, shading='gouraud',
            cMap='Spectral_r', logScale=True, cMin=0.01, cMax=10,
            levels=[10, 55, 100],
            orientation="vertical",
            colorBar=True)

    pg.show(grid, pg.x(grid.cellCenter()), tri=True, shading='gouraud',
            cMap='Spectral_r', logScale=False, cMin=0.01, cMax=10,
            levels=[10, 55, 100],
            orientation="horizontal",
            colorBar=True)


def testColRange():
    n = 5
    mesh = pg.createGrid(n, n)

    fig, ax = plt.subplots(2, 3, figsize=(8, 8))

    ax[0][0].set_title("pyGIMLi (CellData)")
    data = np.arange(mesh.cellCount())
    pg.show(mesh, data, ax=ax[0][0], label='default (nLevs=5, nCols=256)')
    pg.show(mesh, data, ax=ax[1][0], nCols=4, nLevs=5,
            label="nLevs=5 nCols=4, \n(aka. poor man's contour)")

    ax[0][1].set_title("pyGIMLi (CellData) shading='gouraud'")
    data = np.arange(mesh.cellCount())
    pg.show(mesh, data, ax=ax[0][1], shading='gouraud',
            label='default (nLevs=5, nCols=256)')
    pg.show(mesh, data, ax=ax[1][1], shading='gouraud', nCols=4, nLevs=5,
            label="nLevs=5 nCols=4, \n(aka. poor man's contour)")

    ax[0][2].set_title("pyGIMLi (NodeData)")
    data = np.arange(mesh.nodeCount())
    pg.show(mesh, data, ax=ax[0][2],
            label='default (nLevs=5, nCols=4)')
    pg.show(mesh, data, ax=ax[1][2], nCols=12, nLevs=4,
            label='nLevs=4, nCols=12')

    for a in ax.flatten():
        a.yaxis.set_visible(False)
        a.xaxis.set_visible(False)

    fig.tight_layout()


def testCBarLevels():
    """Test colorbar levels with pygimli.show.

    Expectations:

    axs[0, 0]: show regions with plc
        Show needs to deliver the regions with Set3 colormap. Each tick on the
        colobar should be in the middle of the related color section.

    axs[1, 0]: show regions with mesh
        really the same thing as on axs[0, 0] but with mesh.

    ax[0, 1]: show mesh with cell data
        if nLevs is given i would expect that the colormap then is levelled.
        currently that is not the fact. but at least its the full range. labels
        need to be at begin/end of each color section.

    ax[1, 1]: show mesh with node data
        the colorbar range misses parts of its full range. labels need to be
        at begin/end of each color section.
    """
    # create a geometry
    world = mt.createWorld(
        start=[-10, 0], end=[10, -16], layers=[-8], worldMarker=False)

    block = mt.createRectangle(
        start=[-6, -3.5], end=[6, -6.0], marker=4, boundaryMarker=10, area=0.1)

    circ = mt.createCircle(pos=[0, -11], marker=5, radius=2, boundaryMarker=11)

    poly = mt.mergePLC([world, block, circ])
    mesh = mt.createMesh(poly, quality=34)

    # create random data
    rhomap = [[1, 10], [2, 4], [4, 20], [5, 8]]
    # map data to cell/node count
    cell_data = pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)
    node_data = mt.cellDataToNodeData(mesh, cell_data)

    # plot everything
    fig, axs = pg.plt.subplots(2, 2, figsize=(20, 10))

    pg.show(poly, ax=axs[0, 0])

    pg.show(mesh, ax=axs[1, 0], markers=True)

    pg.show(mesh, cell_data, ax=axs[0, 1], colorBar=True, nLevs=7)

    pg.show(mesh, node_data, ax=axs[1, 1], colorBar=True, nLevs=7)


def testColorBarFalse():
    from pygimli.physics import ert
    data = ert.createData(6, "dd")

    rho = np.arange(data.size())
    fig, ax = pg.plt.subplots(ncols=3, nrows=2)

    mesh = pg.createGrid(3, 2)
    pg.show(mesh, rho, ax=ax[0][0], colorBar=True, cMap="plasma")
    pg.show(mesh, rho, ax=ax[0][1], colorBar=False, cMap="plasma")
    pg.show(mesh, rho, ax=ax[0][2], colorBar=False, cMap="plasma",
            cMin=1, cMax=2)

    ert.show(data, rho, ax=ax[1][0], colorBar=True, cMap="plasma")
    ert.show(data, rho, ax=ax[1][1], colorBar=False, cMap="plasma")
    ert.show(data, rho, ax=ax[1][2], colorBar=False, cMap="plasma",
             cMin=1, cMax=2)


def testShowPV():
    """Test pyvista backend for pygimli.show.

    import pygimli as pg
    from pygimli.testing.test_show import testShowPV
    import pyvista as pv
    print(pv.__version__)
    testShowPV()
    """
    m1 = mt.createCube()
    m1.setBoundaryMarkers(range(m1.boundaryCount()))

    m2 = mt.createCube(boundaryMarker=8)
    m2.scale(0.5)

    pg.rc['view3D'] = 'pyvista'
    print('Show Boundary:', m1)
    pg.show(m1+m2, bc='#DDDDFF', alpha=0.5)
    m1 = mt.createMesh(m1+m2, area=0.1)

    print('Show Cells:', m1)
    m1.setCellMarkers(range(m1.cellCount()))
    pg.show(m1, showMesh=True)

    print('Show Markers with vtk filters:', m1)
    pg.show(m1, markers=True, filter={'clip': {'origin': (0, 0, 0)}})

    print('Show Field (x)')
    pg.show(m1, data=pg.x(m1), label='x', cMap='Spectral_r')

    print('Show threshold and slice')
    m2 = mt.createGrid(4, 4, 4)
    m2["m"] = pg.utils.grange(1, m2.cellCount(), 1)
    ax, _ = pg.show(m2, style="wireframe", hold=True)
    pg.viewer.pv.drawMesh(ax, m2, label="m", filter={
        "threshold": dict(value=14, scalars="m", invert=True)})
    pg.viewer.pv.drawMesh(ax, m2, label="m", filter={
        "slice": dict(origin=[1, 1, 1], normal=[1, 1, 0])})
    ax.show()

    print('Show TetP2 Field (x)')
    m2 = mt.refineHex2Tet(m2)
    m2 = m2.createP2()
    pg.show(m2, data=pg.z(m2), showMesh=True, label='x', cMap='Spectral_r')

    print('Show Streams (x)')
    u = pg.x(m1)
    vel = -pg.solver.grad(m1, u)

    ax, _ = pg.show(m1, data=u, label='x', alpha=0.2, hold=True)
    pg.viewer.pv.drawStreamLines(ax, m1, vel, radius=0.01,
                                 source_radius=0.5,
                                 source_center=[-.5, -0., -0.])
    pg.viewer.pv.drawSlice(ax, m1, data=u, label='x', normal=[0, 1, 0])
    ax.show()

def testShowPV2():
    print('Show HexP2 Mesh (x)')
    h1 = mt.createGrid(2, 2, 2)
    h1 = h1.createP2()
    h1.show(pg.x(h1))



def testPVBackends():
    """Test different available backends for pygimli.show.

    import pygimli as pg
    from pygimli.testing.test_show import testPVBackends
    import pyvista as pv
    testPVBackends()
    """
    m1 = mt.createCube()


    # pg.rc['view3D'] = 'fallback'
    # pg.show(m1)

    pg.rc['view3D'] = 'pyvista'
    print('Default')
    pg.show(m1)
    print('Trame-client')
    pg.show(m1, backend='client')
    # print('ipyvtklink')
    # pg.show(m1, backend='ipyvtklink')


def testCoverage():
    grid = pg.createGrid(10, 10)
    cov = pg.y(grid.cellCenters()).array()
    cov[-9:] = 0  # remove first row
    data = pg.Vector(grid.cellCount(), 1.0)
    pg.show(grid, data, coverage=cov)


def testAnimations():
    mesh = pg.createGrid(20,20)
    data = np.random.randn(10, mesh.cellCount())
    pg.show(mesh, data, cMin=0, cMax=1)

    # data = np.random.randn(10, mesh.nodeCount())
    # pg.show(mesh, data)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        locals()[sys.argv[1]]()
    else:
        pass
        #testShowVariants()
        testColorbar()
        #testColorBarFalse()
        #testShowPV()
        #testCBarLevels()
        #testColRange()
        # testCoverage()
