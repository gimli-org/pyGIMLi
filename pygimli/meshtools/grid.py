#!/usr/bin/env python3
"""General grid generation and maintenance."""

import numpy as np
import pygimli as pg


def createGrid(x=None, y=None, z=None, **kwargs):
    """Create grid style mesh.

    Generate simple grid with defined node positions for each dimension.
    The resulting grid depends on the amount of given coordinate arguments
    and consists out of edges (1D - x), quads (2D- x and y), or
    hexahedrons(3D- x, y, and z). For grids with world boundary markers
    (-1 surface and -2 subsurface) the y or z array need to be in increasing
    order.

    Parameters
    ----------
    kwargs:
        x: array
            x-coordinates for all Nodes (1D, 2D, 3D)
        y: array
            y-coordinates for all Nodes (2D, 3D)
        z: array
            z-coordinates for all Nodes (3D)
        marker: int = 0
            Marker for resulting cells.
        worldBoundaryMarker : bool = False
            Boundaries are enumerated with world marker, i.e., Top = -1
            All remaining = -2.
            Default markers: left=1, right=2, top=3, bottom=4, front=5, back=6

    Returns
    -------
    :gimliapi:`GIMLI::Mesh`
        Either 1D, 2D or 3D mesh depending the input.


    Examples
    --------
    >>> import pygimli as pg
    >>> mesh = pg.meshtools.createGrid(x=[0, 1, 1.5, 2], y=[-1, -.5, -.25, 0],
    ...                                marker=2)
    >>> print(mesh)
    Mesh: Nodes: 16 Cells: 9 Boundaries: 24
    >>> fig, axs = pg.plt.subplots(1, 2)
    >>> _ = pg.show(mesh, markers=True, showMesh=True, ax=axs[0])
    >>> mesh = pg.meshtools.createGrid(x=[0, 1, 1.5, 2], y=[-1, -.5, -0.25, 0],
    ...                                worldBoundaryMarker=True, marker=2)
    >>> print(mesh)
    Mesh: Nodes: 16 Cells: 9 Boundaries: 24
    >>> _ = pg.show(mesh, markers=True, showBoundaries=True,
    ...             showMesh=True, ax=axs[1])
    """
    if 'degree' in kwargs:
        return createGridPieShaped(x, **kwargs)

    if x is not None:
        if isinstance(x, int):
            x = list(range(x))
        kwargs['x'] = x
    if y is not None:
        if isinstance(y, int):
            y = list(range(y))
        kwargs['y'] = y
    if z is not None:
        if isinstance(z, int):
            z = list(range(z))
        kwargs['z'] = z

    return pg.core.createGrid(**kwargs)


def createGridPieShaped(x, degree=10.0, h=2, marker=0):
    """Create a 2D pie shaped grid (segment from annulus or cirlce).

    TODO
    ----
        * degree: > 90 .. 360

    Arguments
    ---------
    x: array
        x-coordinates for all Nodes (2D). If you need it 3D, you can apply
        :py:mod:`pygimli.meshtools.extrudeMesh` on it.

    degree: float [None]
        Create a pie shaped grid for a value between 0 and 90.
        Creates an optional inner boundary (marker=2) for a annulus with
        x[0] > 0. Outer boundary marker is 1. Optional h refinement.
        Center node is the first for circle segment.

    h: int [2]
        H-Refinement for degree option.

    marker: int = 0
        Marker for resulting cells.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    Examples
    --------
    >>> import pygimli as pg
    >>> mesh = pg.meshtools.createGridPieShaped(x=[0, 1, 3], degree=45, h=3)
    >>> print(mesh)
    Mesh: Nodes: 117 Cells: 128 Boundaries: 244
    >>> _ = pg.show(mesh)
    >>> mesh = pg.meshtools.createGridPieShaped(x=[1, 2, 3], degree=45, h=3)
    >>> print(mesh)
    Mesh: Nodes: 153 Cells: 128 Boundaries: 280
    >>> _ = pg.show(mesh)
    """
    mesh = pg.Mesh(dim=2)

    for i in range(0, len(x)):
        mesh.createNodeWithCheck([x[i], 0.0])

        mesh.createNodeWithCheck([x[i]*np.cos(degree*np.pi/180),
                                  x[i]*np.sin(degree*np.pi/180)])

    if abs(x[0]) < 1e-6:
        mesh.createCell([0, 1, 2])
        for i in range(0, (len(x)-2)*2-1, 2):
            c = mesh.createCell([i+1, i+3, i+4, i+2])
    else:
        for i in range(0, len(x)*2-2, 2):
            c = mesh.createCell([i, i+2, i+3, i+1])
        mesh.createBoundary([0, 1], marker=1)

    mesh.createBoundary([mesh.nodeCount()-2, mesh.nodeCount()-1], marker=2)

    for i in range(h):
        mesh = mesh.createH2()
    mesh.createNeighbourInfos()

    for b in mesh.boundaries():
        if b.outside() and b.marker() == 0:
            if b.norm()[1] == 0.0:
                b.setMarker(4)  # bottom
            else:
                b.setMarker(3)

    meshR = pg.Mesh(mesh)

    # move all nodes on the inner boundary to rw
    for b in mesh.boundaries():
        line = pg.Line(b.node(0).pos(), b.node(1).pos())

        rSoll = line.intersect([0.0, 0.0], [1.0, 0.0])[0]
        if rSoll > 1e-4:
            for n in b.nodes():
                scale = rSoll/n.pos().abs()
                if scale > 1:
                    meshR.node(n.id()).setPos(pg.Line([0.0, 0.0],
                                                      n.pos()).at(scale))

    if marker != 0:
        for c in meshR.cells():
            c.setMarker(marker)

    return meshR


def appendBoundary(mesh, **kwargs):
    """Append Boundary to a given mesh.

    Syntactic sugar for :py:mod:`pygimli.meshtools.appendTriangleBoundary`
    and :py:mod:`pygimli.meshtools.appendTetrahedronBoundary`.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        "2d or 3d Mesh to which the boundary will be appended.

    Keyword Args
    ------------
    **kwargs forwarded to :py:mod:`pygimli.meshtools.appendTriangleBoundary`
    or :py:mod:`pygimli.meshtools.appendTetrahedronBoundary`.

    Returns
    -------
    :gimliapi:`GIMLI::Mesh`
        A new 2D or 3D mesh containing the original mesh and a boundary around.
    """
    if mesh.dim() == 2:
        return appendTriangleBoundary(mesh, **kwargs)
    elif mesh.dim() == 3:
        return appendTetrahedronBoundary(mesh, **kwargs)

    pg.critical("Don't know how to append boundary to: ", mesh)


def appendTriangleBoundary(mesh, xbound=10, ybound=10, marker=1,
                           isSubSurface=True,
                           **kwargs):
    """Add a triangle mesh boundary to a given mesh.

    Returns a new mesh that contains a triangulated box around a given mesh
    suitable for geo-simulation (surface boundary with marker == -1  at the top
    and marker == -2 in the inner subsurface). The old boundary marker from the
    input mesh are preserved, except for marker == -2 which will be changed
    to +2 as we assume marker == -2 is the world marker for outer boundaries in
    the subsurface.

    Note, this all will only work stable if the mesh generator (triangle)
    preserve all input boundaries. However this will lead to bad quality meshes
    for the boundary region so its a good idea to play with the addNodes keyword
    argument to manually refine the newly created outer boundaries.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh to which the triangle boundary should be appended.
    xbound: float, optional
        Absolute horizontal prolongation distance. Need to be greater 2.
    ybound: float, optional
        Absolute vertical prolongation distance.
    marker: int, optional
        Marker of new cells.
    isSubSurface: boolean [True]
        Apply boundary conditions suitable for geo-simulation and prolongate
        mesh to the surface if necessary.

    Keyword Args
    ------------
    ** kargs forwarded to pg.createMesh

    quality : float, optional
        Triangle quality.
    area: float, optional
        Triangle max size within the boundary.
    smooth : boolean, optional
        Apply mesh smoothing.

    addNodes : int[5], iterable
        Add additional nodes on the outer boundaries. Or for each boundary
        if given 5 values (isSubsurface=True) or 4 for isSubsurface=False

    Returns
    -------
    :gimliapi:`GIMLI::Mesh`
        A new 2D mesh containing the original mesh and a boundary arround.

    See Also
    --------
    :py:mod:`pygimli.meshtools.appendBoundary`
    :py:mod:`pygimli.meshtools.appendTetrahedronBoundary`

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawMesh, drawModel
    >>> import pygimli.meshtools as mt
    >>> inner = pg.createGrid(range(5), range(5), marker=1)
    >>> fig, axs = plt.subplots(2,3)
    >>> ax, _ = pg.show(inner, markers=True, showBoundaries=True, showMesh=True, ax=axs[0][0])
    >>> m = mt.appendTriangleBoundary(inner, xbound=3, ybound=3, marker=2, addNodes=0, isSubSurface=False)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[0][1])
    >>> m = mt.appendTriangleBoundary(inner, xbound=4, ybound=1, marker=2, addNodes=5, isSubSurface=False)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[0][2])
    >>> m = mt.appendTriangleBoundary(inner, xbound=4, ybound=4, marker=2, addNodes=0, isSubSurface=True)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[1][0])
    >>> m = mt.appendTriangleBoundary(inner, xbound=4, ybound=4, marker=2, addNodes=5, isSubSurface=True)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[1][1])
    >>> surf = mt.createPolygon([[0, 0],[5, 3], [10, 1]], boundaryMarker=-1, addNodes=5, interpolate='spline')
    >>> m = mt.appendTriangleBoundary(surf, xbound=4, ybound=4, marker=2, addNodes=5, isSubSurface=True)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[1][2])
    """
    poly = pg.Mesh(isGeometry=True)

    if isSubSurface:

        bs = mesh.findBoundaryByMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

        if len(bs) == 0:
            for b in mesh.boundaries():
                if b.outside() and b.norm()[1] == 1.0:
                    bs.append(b)

        paths = mesh.findPaths(bs)

        if len(paths) > 0:
            startPoint = mesh.node(paths[0][0]).pos()
            endPoint = mesh.node(paths[0][-1]).pos()

            if startPoint[0] > endPoint[0]:
                startPoint = endPoint
                endPoint = mesh.node(paths[0][0]).pos()
        else:
            pg.critical("Can't identify upper part of the mesh to be moved to"
                        "the surface. Maybe define them with Marker==-1")

        addNodes = kwargs.pop('addNodes', 5)

        boundPoly = [pg.Pos(startPoint)]

        if isinstance(addNodes, (float, int)) and addNodes > 0:
            addNodes = np.full(5, addNodes)

        if hasattr(addNodes, '__len__') and len(addNodes) == 5:
            boundPoly.extend([boundPoly[-1] - pg.Pos(x, 0)
                              for x in pg.utils.grange(1, xbound,
                                                       n=addNodes[0]+1,
                                                       log=True)])
            boundPoly.extend([boundPoly[-1] - pg.Pos(0, y)
                              for y in np.linspace(
                0, mesh.ymax() - mesh.ymin() + ybound, addNodes[1] + 1)[1:]])
            boundPoly.extend(
                [boundPoly[-1] + pg.Pos(x, 0)
                 for x in np.linspace(0, (endPoint-startPoint)[0] + 2*xbound,
                                      addNodes[2]+1)[1:]])
            boundPoly.extend(
                [boundPoly[-1] + pg.Pos(0, y)
                 for y in np.linspace(0, endPoint[1]-boundPoly[-1][1],
                                      addNodes[3]+1)[1:]])
            boundPoly.extend(
                [boundPoly[-1] - pg.Pos(xbound-x, 0)
                 for x in pg.utils.grange(1, xbound, n=addNodes[4]+1,
                                          log=True)[::-1][1:]])
        else:
            boundPoly.append(boundPoly[-1] - pg.Pos(xbound, 0))
            boundPoly.append(boundPoly[-1] - pg.Pos(0, mesh.ymax() -
                                                    mesh.ymin() + ybound))
            boundPoly.append(boundPoly[-1] + pg.Pos((endPoint-startPoint)[0] +
                                                    2*xbound, 0))
            boundPoly.append(pg.Pos(endPoint) + pg.Pos(xbound, 0))

        boundPoly.append(pg.Pos(endPoint))

        poly = pg.meshtools.createPolygon(boundPoly, isClosed=False)

        poly.addRegionMarker(pg.Pos([poly.xmin(), poly.ymin()]) +
                             [xbound/100, ybound/100],
                             marker=marker)

        if mesh.cellCount() > 0:
            poly.addHoleMarker(pg.Pos([mesh.xmin(), mesh.ymin()]) +
                               [0.001, 0.001])

    else:  # no isSubSurface

        boundPoly = [[mesh.xmin() - xbound, mesh.ymin() - ybound],
                     [mesh.xmin() - xbound, mesh.ymax() + ybound],
                     [mesh.xmax() + xbound, mesh.ymax() + ybound],
                     [mesh.xmax() + xbound, mesh.ymin() - ybound]]

        poly = pg.meshtools.createPolygon(boundPoly, isClosed=True,
                                          marker=marker,
                                          addNodes=kwargs.pop('addNodes', 5))

    for b in mesh.boundaries():
        if b.outside() or b.marker() == -1:
            poly.copyBoundary(b)

    preserveSwitch = 'Y'
    #pg.show(poly, fillRegion=False)
    #pg.show(poly, boundaryMarkers=True, showNodes=True)
    #pg.wait()

    mesh2 = pg.meshtools.createMesh(poly, preserveBoundary=preserveSwitch,
                                    **kwargs)

    # pg.show(mesh2, boundaryMarkers=True, showNodes=True)

    # start extracting all cells with marker from mesh2 and all orginal cells
    mesh3 = pg.Mesh(2)

    for c in mesh2.cells():
        if c.marker() == marker:
            mesh3.copyCell(c)

    # ! map copies the cell not the reference, this should not happen
    # **TODO check 20210305
    # map(lambda cell: mesh2.copyCell(cell), mesh2.cells())
    for c in mesh.cells():
        mesh3.copyCell(c)

    # we need to delete the old boundary markers or the new neighbour infos
    # will fail for old outside boundaries
    mesh3.setBoundaryMarkers(np.zeros(mesh3.boundaryCount()))
    mesh3.createNeighborInfos(force=True)

    for b in mesh.boundaries():
        if b.marker() != 0:
            b2 = mesh3.copyBoundary(b)

            # some automagic: original mesh contains bmarker == -2 which means
            # mixed condition, this special marker will be switched to 2
            if b.marker() == -2:
                b2.setMarker(2)

    for b in mesh3.boundaries():
        if b.outside() and b.marker() > -1:
            if b.norm().x() != 0 or b.norm().y() == -1.0 or not isSubSurface:
                b.setMarker(pg.core.MARKER_BOUND_MIXED)
            else:
                b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

    return mesh3


def appendBoundaryGrid(grid, xbound=None, ybound=None, zbound=None,
                       marker=1, isSubSurface=True, **kwargs):
    """Return a copy of grid surrounded by a boundary grid.

    Note, the input grid needs to be a 2d or 3d grid with quad/hex cells.

    TODO
    ----
        * preserve inner boundaries
        * add subsurface setting

    Args
    ----
    grid: :gimliapi:`GIMLI::Mesh`
        2D or 3D Mesh that must contain structured quads or hex cells
    xbound: iterable of type float [None]
        Needed for 2D or 3D grid prolongation and will be added on the left
        side in opposit order and on the right side in normal order.
    ybound: iterable of type float [None]
        Needed for 2D or 3D grid prolongation and will be added (2D bottom,
        3D front) in opposit order and (2D top, 3D back) in normal order.
    zbound: iterable of type float [None]
        Needed for 3D grid prolongation and will be added the bottom side in
        opposite order on the top side in normal order.
    marker: int [1]
        Cellmarker for the cells in the boundary region
    isSubSurface : boolean, optional
        Apply boundary conditions suitable for geo-simulaion and prolongate
        mesh to the surface if necessary, e.i., no boundary on top of the grid.

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> grid = mt.createGrid(5,5)
    ...
    >>> g1 = mt.appendBoundaryGrid(grid,
    ...                            xbound=[1, 3, 6],
    ...                            ybound=[1, 3, 6],
    ...                            marker=2,
    ...                            isSubSurface=False)
    >>> ax,_ = pg.show(g1, markers=True, showMesh=True)
    >>> grid = mt.createGrid(5,5,5)
    ...
    >>> g2 = mt.appendBoundaryGrid(grid,
    ...                            xbound=[1, 3, 6],
    ...                            ybound=[1, 3, 6],
    ...                            zbound=[1, 3, 6],
    ...                            marker=2,
    ...                            isSubSurface=False)
    >>> ax, _ = pg.show(g2, g2.cellMarkers(), showMesh=True,
    ...                 filter={'clip':{}});
    """
    if isSubSurface:
        pg.critical('Implement me')

    def _concat(v, vBound):
        if (not pg.isArray(vBound)):
            pg.critical("please give bound array")

        v = np.append(-np.array(vBound)[::-1] + v[0], v)
        v = np.append(v, v[-1] + np.array(vBound))
        return v

    x = None
    y = None
    z = None

    if grid.dim() > 1:
        if grid.dim() == 2:
            if any([c.nodeCount() != 4 for c in grid.cells()]):
                pg.critical("Grid have other cells than quads. "
                            "Can't refine it with a grid")

        x = pg.utils.unique(pg.x(grid))
        y = pg.utils.unique(pg.y(grid))
        x = _concat(x, xbound)
        y = _concat(y, ybound)

        if grid.dim() == 3:
            if any([c.nodeCount() != 8 for c in grid.cells()]):
                pg.critical("Grid have other cells than hex's. "
                            "Can't refine it with a grid")

            z = pg.utils.unique(pg.z(grid))
            z = _concat(z, zbound)

    mesh = pg.meshtools.createGrid(x=x, y=y, z=z, marker=marker)

    mesh.setCellMarkers(pg.interpolate(grid,
                                       grid.cellMarkers(),
                                       mesh.cellCenters(),
                                       fallback=marker))
    return mesh


def appendTetrahedronBoundary(mesh, xbound=10, ybound=10, zbound=10,
                              marker=1, isSubSurface=True, **kwargs):
    """Return a copy of mesh surrounded by a tetrahedron mesh as boundary.

    Returns a new mesh that contains a tetrahedron mesh box around a given mesh
    suitable for geo-simulation (surface boundary with marker = -1  at top and
    marker = -2 in the inner subsurface). The old boundary marker from mesh
    will be preserved, except for marker == -2 which will be switched to 2 as
    we assume -2 is the world marker for outer boundaries in the subsurface.

    Note
    ----
    This method will only work stable if the mesh generator (Tetgen) preserves
    all input boundaries. This will lead to bad quality meshes for the boundary
    region so its a good idea to play with the addNodes keword argument to
    manually refine the newly created outer boundaries.

    If the input mesh consists of hexahedrons a small inconsistency will arise
    because a quad boundary element will be split by 2 triangle boundaries from
    the boundary tetrahedrons. The effect of this hanging edges are unclear,
    also createNeighbourInfos may fail. We need to implement/test pyramid cells
    to handle this.

    TODO
    ----
        * set correct boundary conditions
        * isSubSurface
        * pyramid cells as connecting cells
        * need for preserve Boundary check
        * preserve Boundary support
        * addNodes support

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        3D Mesh to which the tetrahedron boundary should be appended.
    xbound: float [10]
        Horizontal prolongation distance in meter at x-direction.
        Need to be >= 0.
    ybound: float [10]
        Horizonal prolongation distance in meter at y-direction.
        Need to be greater 0.
    zbound: float [10]
        Vertical prolongation distance in meter at z-direction (>0).
    marker: int, optional
        Marker of new cells.
    addNodes: float, optional
        Triangle quality.
    isSubSurface : boolean, optional
        Apply boundary conditions suitable for geo-simulaion and prolongate
        mesh to the surface if necessary.
    verbose : boolean, optional
        Be verbose.

    Returns
    -------
    :gimliapi:`GIMLI::Mesh`
        A new 3D mesh containing the original mesh and a boundary arround.

    See Also
    --------
    :py:mod:`pygimli.meshtools.appendBoundary`,
    :py:mod:`pygimli.meshtools.appendTriangleBoundary`

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> grid = mt.createGrid(5,5,5)
    ...
    >>> mesh = mt.appendBoundary(grid, xbound=5, ybound=5, zbound=5,
    ...                          isSubSurface=False)
    >>> ax, _ = pg.show(mesh, mesh.cellMarkers(), showMesh=True,
    ...                 filter={'clip':{}})
    """
    if isSubSurface:
        pg.critical('Implement me')

    meshBoundary = pg.Mesh(3, isGeometry=True)
    for b in mesh.boundaries():
        if b.outside() or b.marker() == -1:
            meshBoundary.copyBoundary(b)

    bb = meshBoundary.bb()
    meshBoundary.addHoleMarker(bb[0] + (bb[1]-bb[0])/1000.)

    if not any([xbound > 0, ybound > 0, zbound > 0]):
        pg.critical('all boundaries need to be greater 0.')

    startPos = bb[0] - [xbound, ybound, zbound]
    endPos = bb[1] + [xbound, ybound, zbound]

    boundaryBox = pg.meshtools.createCube(start=startPos, end=endPos)

    boundMesh = pg.meshtools.createMesh(boundaryBox + meshBoundary)
    boundMesh.setCellMarkers(np.ones(boundMesh.cellCount()) * marker)

    outMesh = pg.Mesh(mesh)

    for c in boundMesh.cells():
        outMesh.copyCell(c)

    return outMesh


if __name__ == "__main__":
    pass
