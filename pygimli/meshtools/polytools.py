# -*- coding: utf-8 -*-
"""Tools to create or manage PLC.

Please note there is currently no collision or intersection check at all.

Volunteers welcome to help creating, adapting or interfacing a basic
geometry system. A lot of things are needed:

    * 2D
    * 3D
    * More geometric primitives
    * Boolean operations (union, intersection, difference)
    * Collision recognizing
    * Cubic spline interpolation for polygons (partly done)
    * GUI .. interactive creation
    *
"""

import math
import os
from os import system

import functools
import numpy as np
import pygimli as pg


def _polyCreateDefaultEdges(poly, boundaryMarker=1, isClosed=True, **kwargs):
    """INTERNAL."""
    nEdges = poly.nodeCount() - 1 + isClosed
    bm = None
    if hasattr(boundaryMarker, '__len__'):
        if len(boundaryMarker) == nEdges:
            bm = boundaryMarker
        else:
            raise Exception("marker length != nEdges", len(boundaryMarker),
                            nEdges)
    else:
        bm = [boundaryMarker] * nEdges

    for i in range(poly.nodeCount() - 1):
        poly.createEdge(poly.node(i), poly.node(i + 1), bm[i])

    if isClosed:
        poly.createEdge(poly.node(poly.nodeCount() - 1), poly.node(0), bm[-1])


def setPolyRegionMarker(poly, marker=1, area=0.0, markerPosition=None,
                        isHole=False, **kwargs):
    """Internal to set region markers to single elementary geometry.

    Internal to set region markers.
    If no absolute marker position is given.
    The region marker is placed 1mm beside the first node in direction to the
    geometry center.

    Parameters
    ----------
    poly : :gimliapi:`GIMLI::Mesh`
        The Polygon that will get a marker
    marker : int[1]
        The region marker, every resulting mesh cell will get this marker.
    area : float[0]
        The region max cell size, every resulting mesh cell will get a cell
        size lower than area in m² or m³ for 3D, respectively.
    markerPosition : pg.Pos
        Absolute marker position if you don't want the marker in the center of
        the geometry.
    isHole : bool [False]
        Marks the geometry as a hole and will be cut in any merge mesh.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs
    """
    pos = None
    if markerPosition is not None:
        pos = markerPosition
    else:
        center = pg.center(poly.positions())
        p0 = poly.node(0).pos()
        # region marker near node 0; 1mm in direction to the center
        # should be safer than the center itself
        pos = p0 + (center-p0).norm() * 0.001

    if isHole is True:
        poly.addHoleMarker(pos)
    else:
        poly.addRegionMarker(pos, marker=marker, area=area)


def createRectangle(start=None, end=None, pos=None, size=None, **kwargs):
    """Create rectangle polygon.

    Create rectangle with start position and a given size.
    Give either start and end OR pos and size.

    Parameters
    ----------
    start : [x, y]
        Left upper corner. Default [-0.5, 0.5]
    end : [x, y]
        Right lower corner. Default [0.5, -0.5]
    pos : [x, y]
        Center position. The rectangle will be moved.
    size : [x, y]
        Factors for x and y by which the rectangle, defined by **start** and
        **width**, are scaled.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs

        marker : int [1]
            Marker for the resulting triangle cells after mesh generation
        markerPosition : floats [x, y] [pos + (end - start) * 0.2]
            Absolute position of the marker (works for both regions and
            holes).
        area : float [0]
            Maximum cell size for resulting triangles after mesh generation
        isHole : bool [False]
            The polygon will become a hole instead of a triangulation
        boundaryMarker : int [1]
            Marker for the resulting boundary edges
        leftDirection : bool [True]
            TODO Rotational direction
        pnts: [[x, y],]
            Return squared rectangle of origin-aligned boundingbox for pnts.
        minBB: False
            Return squared rectangle of non-origin-aligned minimum bounding box
            for pnts.
        minBBOffset: [1.0, 1.0]
            Offset for minimal boundingbox in x and y direction in relative
            extent .. whatever that means for non-aligned boxes.

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> # no need to import matplotlib, pygimli show does.
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> r1 = mt.createRectangle(pos=[1, -1], size=[4.0, 4.0],
    ...                      marker=1, area=0.1, markerPosition=[0, -2])
    >>> r2 = mt.createRectangle(start=[0.5, -0.5], end=[2, -2],
    ...                      marker=2, area=1.1)
    >>> pnts3 = [[-0.5, 0], [0, 0.2], [-0.2, 0.5], [-0.3, 0.25]]
    >>> r3 = mt.createRectangle(pnts=pnts3, marker=3, area=0.2)
    >>> pnts4 = [[1.5, 0], [2.0, 0.2], [1.8, 0.5], [1.7, 0.25]]
    >>> r4 = mt.createRectangle(pnts=pnts4, marker=4, minBB=True)
    >>> pnts5 = [[-0.5, -1], [0, -0.8], [-0.2, -0.5], [-0.3, -0.75]]
    >>> r5 = mt.createRectangle(pnts=pnts5, marker=5,
    ...                         minBB=True, minBBOffset=[1.2, 1.2])
    >>> ax, _ = pg.show(mt.mergePLC([r1, r2, r3, r4, r5]))
    >>> pg.viewer.mpl.drawSensors(ax, pnts3)
    >>> pg.viewer.mpl.drawSensors(ax, pnts4)
    >>> pg.viewer.mpl.drawSensors(ax, pnts5)
    """
    pnts = kwargs.pop('pnts', None)
    if pnts is not None:
        if len(pnts) == 1:
            return createRectangle(pos=pnts[0], size=[1, 1], **kwargs)
        if len(pnts) == 2:
            return createRectangle(start=pnts[0], end=pnts[1], **kwargs)

        minBB = kwargs.pop('minBB', False)
        minBBBoundary = kwargs.pop('minBBOffset', [1.0, 1.0])

        if not minBB:
            xMin = min(pg.x(pnts))
            xMax = max(pg.x(pnts))
            yMin = min(pg.y(pnts))
            yMax = max(pg.y(pnts))

            bb = np.asarray([[xMin, yMin], [xMax, yMax]])

            bbScale = (bb[1]-bb[0])*(minBBBoundary-np.asarray([1.0, 1.0]))
            return createRectangle(start=bb[0]-bbScale, end=bb[1]+bbScale,
                                   **kwargs)
        else:
            # create convex hull
            m = pg.meshtools.createMesh(pnts)
            if m.boundaryCount() == 0:
                # probably linear pnts so no convex hull
                for i in range(m.nodeCount()-1):
                    m.createBoundary([i, i+1])

                    bs = pg.meshtools.createLine(pnts[0], pnts[-1])
            else:
                bs = m.createSubMesh(m.boundaries([b.id()
                                     for b in m.boundaries() if b.outside()]))

            def getBB(m, off, rot):
                # Rotate hull and find bb
                m2 = pg.Mesh(m)
                m2.translate(-off)
                m2.transform(rot)
                m2.translate(off)

                # Increase bb if zero width or lenght
                bb = m2.bb()
                if (bb[1]-bb[0])[0] < 1e-12:
                    bb[1][0] += 0.5
                    bb[0][0] -= 0.5
                if (bb[1]-bb[0])[1] < 1e-12:
                    bb[1][1] += 0.5
                    bb[0][1] -= 0.5

                return bb

            off = 0
            minSize = [9e99, None, None, None]

            for b in bs.boundaries():
                # normalize to origin
                off = b.node(0).pos()
                # rotation to origin x axis
                rot = pg.core.getRotation((b.node(1).pos()-off), [1, 0])

                # get bounding box for normalized mesh
                bb = getBB(bs, off, rot)

                # compare size off bb and collect minimum size
                s = (bb[1]-bb[0]).abs()
                if s < minSize[0]:
                    minSize[0] = s
                    minSize[1] = bb
                    minSize[2] = b.node(1).pos()
                    minSize[3] = off

            # Rotate bb back and create rectangle
            bb = minSize[1]
            off = minSize[3]

            bbScale = (bb[1]-bb[0])*(minBBBoundary-np.asarray([1.0, 1.0]))
            r = createRectangle(start=bb[0]-bbScale, end=bb[1]+bbScale,
                                **kwargs)

            rot = pg.core.getRotation([1, 0], minSize[2]-off)
            r.translate(-off)
            r.transform(rot)
            r.translate(off)

            # pg.show(r)
            # pg.wait()
            return r

    if start is None:
        start = [-0.5, 0.5]
    if end is None:
        end = [0.5, -0.5]

    poly = pg.Mesh(dim=2, isGeometry=True)

    sPos = pg.Pos(start)
    ePos = pg.Pos(end)

    verts = [sPos, [sPos[0], ePos[1]], ePos, [ePos[0], sPos[1]]]

    # TODO refactor with polyCreatePolygon

    if kwargs.pop("leftDirection", False):
        for v in verts[::-1]:
            poly.createNode(v)
    else:
        for v in verts:
            poly.createNode(v)

    # Note that we do not support the usage of start/end AND size/pos. Only one
    # of the pairs. Otherwise strange things will happen with the region
    # markers!
    if size is not None:
        poly.scale(size)
    if pos is not None:
        poly.translate(pos)

    _polyCreateDefaultEdges(poly, **kwargs)

    sPos = poly.bb()[0]
    ePos = poly.bb()[1]

    kwargs['markerPosition'] = kwargs.pop('markerPosition',
                                          sPos + (ePos - sPos) * 0.2)

    setPolyRegionMarker(poly, **kwargs)

    return poly


def createWorld(start, end, marker=1, area=0., layers=None, worldMarker=True,
                **kwargs):
    """Create simple rectangular 2D or 3D world.

    Create simple rectangular [hexagonal] world with appropriate boundary
    conditions. Surface boundary is set  pg.core.MARKER_BOUND_HOMOGEN_NEUMANN,
    and inner subsurface is set to pg.core.MARKER_BOUND_MIXED, i.e., -2 OR
    Numbered: 1, 2, 3, 4, 5, 6 for left, right, bottom, top, front and back,
    if worldMarker is set to false and no layers are given. With layers, it is
    numbered in ascending order.

    TODO
    ----
        * 3D with layers

    Parameters
    ----------
    start: [x, y, [z]]
        Upper/Left/[Front] Corner
    end: [x, y, [z]]
        Lower/Right/[Back] Corner
    marker: int
        Marker for the resulting triangle cells after mesh generation.
    area: float | list
        Maximum cell size for resulting triangles after mesh generation.
        If area is a float set it global, if area is a list set it per layer.
    layers: [float] [None]
        List of depth coordinates for some layers.
    worldMarker: bool [True]
        Specify boundary markers:
        True: [-1, -2] for [surface, subsurface] boundaries
        False: ascending order [1, 2, 3, 4 ..]

    Other Parameters
    ----------------
    Forwarded to createCube

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import createWorld
    >>> from pygimli.viewer.mpl import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = createWorld(start=[-5, 0], end=[5, -5], layers=[-1,-2,-3])
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, world)
    >>> plt.show()
    """
    if len(start) == 3 and len(end) == 3:

        if layers is not None:
            pg.critical("3D with layers is not yet implemented.")

        world = createCube(size=pg.Pos(end)-pg.Pos(start),
                           pos=(pg.Pos(end)+pg.Pos(start))/2.0,
                           area=area, marker=marker, **kwargs)

        for i, b in enumerate(world.boundaries()):
            if worldMarker is True:
                if b.norm()[2] == 1.0:
                    b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
                else:
                    b.setMarker(pg.core.MARKER_BOUND_MIXED)
            else:
                if b.norm() == [-1, 0, 0]:
                    b.setMarker(1)
                elif b.norm() == [1, 0, 0]:
                    b.setMarker(2)
                elif b.norm() == [0, 0, -1]:
                    b.setMarker(3)
                elif b.norm() == [0, 0, 1]:
                    b.setMarker(4)
                elif b.norm() == [0, -1, 0]:
                    b.setMarker(5)
                elif b.norm() == [0, 1, 0]:
                    b.setMarker(6)

        return world

    z = np.array(start[1], dtype=float)
    if layers is not None:
        z = np.append(z, np.array(layers, dtype=float))
    z = np.append(z, end[1])
    # ensure - decreasing order if layers are out of bounding box
    z = np.sort(z)[::-1]

    poly = pg.Mesh(dim=2, isGeometry=True)

    if isinstance(area, float) or isinstance(area, int):
        area = np.ones(len(z)-1) * float(area)

    if len(area) < len(z) - 1:
        pg.warn('Missing {} area value, padding with zeros'.format(
            (len(z) - 1) - len(area)))
        _area = np.zeros(len(z)-1)
        _area[0:len(area)] = area
        area = _area

    for i, depth in enumerate(z):
        n = poly.createNode([start[0], depth])
        if i > 0:
            if len(z) == 2:
                poly.addRegionMarker(n.pos() + [0.2, 0.2],
                                     marker=marker, area=area[0])
            else:
                poly.addRegionMarker(n.pos() + [0.2, 0.2],
                                     marker=i, area=area[i - 1])

    for i, depth in enumerate(z[::-1]):
        poly.createNode([end[0], depth])

    _polyCreateDefaultEdges(poly,
                            boundaryMarker=range(1, poly.nodeCount() + 1))

    if worldMarker:
        for b in poly.boundaries():
            if b.norm()[1] == 1.0:
                b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
            else:
                b.setMarker(pg.core.MARKER_BOUND_MIXED)
    elif layers is None:
        for b in poly.boundaries():
            if b.norm() == [-1, 0]:
                b.setMarker(1)
            elif b.norm() == [1, 0]:
                b.setMarker(2)
            elif b.norm() == [0, -1]:
                b.setMarker(3)
            elif b.norm() == [0, 1]:
                b.setMarker(4)

    if layers is not None:
        for i in range(len(layers)):
            poly.createEdge(poly.node(i + 1),
                            poly.node(poly.nodeCount() - i - 2),
                            poly.boundaryCount() + 1)

    # pg.warnNonEmptyArgs(kwargs)
    return poly


def createCircle(pos=None, radius=1, nSegments=12, start=0, end=2.*math.pi,
                 **kwargs):
    """Create simple circle polygon.

    Create simple circle polygon with given attributes.

    Parameters
    ----------
    pos : [x, y] [[0.0, 0.0]]
        Center position
    radius : float | [a,b] [1]
        radius or halfaxes of the circle
    nSegments : int [12]
        Discrete amount of segments for the circle.
    start : double [0]
        Starting angle in radians
    end : double [2*pi]
        Ending angle in radians

    **kwargs:

        marker: int [1]
            Marker for the resulting triangle cells after mesh generation
        markerPosition : floats [x, y] [0.0, 0.0]
            Position of the marker (works for both regions and holes)
        area: float [0]
            Maximum cell size for resulting triangles after mesh generation
        isHole: bool [False]
            The polygon will become a hole instead of a triangulation
        boundaryMarker: int [1]
            Marker for the resulting boundary edges
        leftDirection: bool [True]
            Rotational direction
        isClosed: bool [True]
            Add closing edge between last and first node.

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>>  # no need to import matplotlib. pygimli's show does
    >>> import math
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawMesh
    >>> import pygimli.meshtools as mt
    >>> c0 = mt.createCircle(pos=(-5.0, 0.0), radius=2, nSegments=6)
    >>> c1 = mt.createCircle(pos=(-2.0, 2.0), radius=1, area=0.01, marker=2)
    >>> c2 = mt.createCircle(pos=(0.0, 0.0), nSegments=5, start=0, end=math.pi)
    >>> c3 = mt.createCircle(pos=(5.0, 0.0), nSegments=3, start=math.pi,
    ...                      end=1.5*math.pi, isClosed=False)
    >>> plc = mt.mergePLC([c0, c1, c2, c3])
    >>> fig, ax = pg.plt.subplots()
    >>> drawMesh(ax, plc, fillRegion=False)
    >>> pg.wait()
    """
    pg.renameKwarg('segments', 'nSegments', kwargs, '1.2')  # 20210312
    nSegments = kwargs.pop('nSegments', nSegments)

    # TODO refactor with polyCreatePolygon
    if pos is None:
        pos = [0.0, 0.0]

    poly = pg.Mesh(dim=2, isGeometry=True)

    dPhi = (end - start) / (nSegments)
    nPhi = nSegments + 1

    if abs((end % (2. * math.pi) - start)) < 1e-6:
        nPhi = nSegments

    for i in range(0, nPhi):
        if kwargs.pop('leftDirection', True):
            phi = start + i * dPhi
        else:
            phi = start - i * dPhi

        xp = np.cos(phi)
        yp = np.sin(phi)
        poly.createNode([xp, yp])

    if hasattr(radius, '__len__'):
        poly.scale(radius)
    else:
        poly.scale([radius, radius])
    poly.translate(pos)

    _polyCreateDefaultEdges(poly, **kwargs)

    if kwargs.pop('isClosed', True):
        setPolyRegionMarker(poly, **kwargs)

    # need a better way mess with these or wrong kwargs
    # pg.warnNonEmptyArgs(kwargs)

    return poly


def createLine(start, end, nSegments=1, **kwargs):
    """Create simple line polygon.

    Create simple line polygon from start to end.

    Parameters
    ----------
    start : [x, y]
        start position
    end : [x, y]
        end position
    nSegments : int
        Discrete amount of segments for the line

    Keyword Arguments
    -----------------
    boundaryMarker : int [1]
        Marker for the resulting boundary edges
    leftDirection : bool [True]
        Rotational direction

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>>  # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>>
    >>> w = mt.createWorld(start=[0, 0], end=[3, 3])
    >>> l1 = mt.createLine(start=[1, 1], end=[1, 2], nSegments=1,
    ...                    leftDirection=False)
    >>> l2 = mt.createLine(start=[1, 1], end=[2, 1], nSegments=20,
    ...                    leftDirection=True)
    >>>
    >>> ax, _ = pg.show(mt.createMesh([w, l1, l2,]))
    >>> ax, _ = pg.show([w, l1, l2,], ax=ax, fillRegion=False)
    >>> pg.wait()
    """
    pg.renameKwarg('segments', 'nSegments', kwargs, '1.2')  # 20210312
    nSegments = kwargs.pop('nSegments', nSegments)

    # TODO refactor with polyCreatePolygon
    poly = pg.Mesh(dim=2, isGeometry=True)
    startPos = pg.RVector3(start)
    endPos = pg.RVector3(end)
    a = endPos - startPos

    dt = 1. / nSegments
    left = kwargs.pop('leftDirection', True)

    for i in range(0, nSegments + 1):
        if left:
            p = startPos + a * (dt * i)
        else:
            p = endPos - a * (dt * i)

        poly.createNode(p)

    _polyCreateDefaultEdges(poly, isClosed=False, **kwargs)
    return poly


def createPolygon(verts, isClosed=False, addNodes=0, interpolate='linear',
                  **kwargs):
    """Create a polygon from a list of vertices.

    All vertices need to be unique and duplicate vertices will be ignored.
    If you want the polygon be a closed region you can set the 'isClosed' flag.
    Closed region can be attributed by assigning a region marker.
    The automatic region marker is placed in the center of all vertices.

    Parameters
    ----------
    verts : []
        * List of x y pairs [[x0, y0], ... ,[xN, yN]]

    isClosed : bool [True]
        Add closing edge between last and first node.

    addNodes : int [1], iterable
        Constant or (for each) Number of additional nodes to be added,
        equidistant between sensors.

    interpolate : str ['linear']
        Interpolation rule for addNodes. 'linear' or 'spline'. TODO 'harmfit'

    **kwargs:

        marker : int [None]
            Marker for the resulting triangle cells after mesh generation.
        markerPosition : floats [x, y] [0.0, 0.0]
            Position (absolute) of the marker (works for both regions and
            holes)
        area : float [0]
            Maximum cell size for resulting triangles after mesh generation
        isHole : bool [False]
            The polygon will become a hole instead of a triangulation
        boundaryMarker : int [1]
            Marker for the resulting boundary edges
        leftDirection : bool [True]
            Rotational direction

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> # no need to import matplotlib, pygimli show does.
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> p1 = mt.createPolygon([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
    ...                       isClosed=True, marker=3, area=0.1)
    >>> p2 = mt.createPolygon([[0.3, 0.15], [0.85, 0.15], [0.85, 0.7]],
    ...                       isClosed=True, isHole=True)
    >>> p3 = mt.createPolygon([[-0.1, 0.2], [-1.1, 0.2], [-1.1, 1.2], [-0.1, 1.2]],
    ...                       isClosed=True, addNodes=3, marker=2)
    >>> p4 = mt.createPolygon([[-0.1, 0.2], [-1.1, 0.2], [-1.1, 1.2], [-0.1, 1.2]],
    ...                       isClosed=True, addNodes=5, interpolate='spline',
    ...                       marker=4)
    >>> ax, _ = pg.show(mt.mergePLC([p1, p2, p3, p4]), showNodes=True)
    >>> pg.wait()
    """
    poly = pg.Mesh(dim=2, isGeometry=True)

    if hasattr(addNodes, '__iter__') or addNodes > 0:
        if isClosed:
            verts = np.array(verts)
            verts = np.vstack([verts, verts[0]])

        tV = pg.utils.cumDist(verts)

        if isinstance(addNodes, int) and addNodes > 0:
            addNodes = np.full(len(tV)-1, addNodes)

        if len(addNodes) != len(tV)-1:
            print(addNodes)
            pg.error('Amount of addNodes does not match needed length:',
                     len(tV)-1)

        tI = []

        for i, t in enumerate(tV[0:len(tV)-1]):
            tI.append(t)
            for j in range(addNodes[i]):
                dt = (tV[i+1]-tV[i]) / (addNodes[i]+1)
                tI.append(tV[i] + dt*(j+1))

        if not isClosed:
            tI.append(tV[-1])

        verts = pg.meshtools.interpolate(verts, tI,
                                         method=interpolate,
                                         periodic=isClosed)

    if kwargs.pop("leftDirection", False):
        for v in verts[::-1]:
            if isinstance(v, float) or isinstance(v, int):
                poly.createNodeWithCheck([v, 0], warn=True)
            else:
                poly.createNodeWithCheck(v, warn=True)
    else:
        for v in verts:
            if isinstance(v, float) or isinstance(v, int):
                poly.createNodeWithCheck([v, 0], warn=True)
            else:
                poly.createNodeWithCheck(v, warn=True)

    _polyCreateDefaultEdges(poly, isClosed=isClosed,
                            boundaryMarker=kwargs.pop('boundaryMarker', 1))

    if isClosed:
        setPolyRegionMarker(poly, **kwargs)

    return poly


def merge(*args, **kwargs):
    """Little syntactic sugar to merge.

    All args are forwarded to mergeMeshes if isGeometry is not set.
    Otherwise it considers the mesh as PLC to merge.

    Args
    ----
    List of meshes or comma separated list of meshes that will be forwarded to
    mergeMeshes or meshPLC.
    """
    if len(args) == 1 and isinstance(args[0], list):
        return merge(*args[0], **kwargs)

    for arg in args:
        if hasattr(arg, 'isGeometry'):
            if arg.isGeometry():
                return mergePLC([*args], **kwargs)

    return pg.meshtools.mergeMeshes([*args], **kwargs)


def mergePLC(plcs, tol=1e-3):
    """Merge multiply polygons.

    Merge multiply polygons into a single polygon.
    Common nodes and common edges will be checked and removed.
    When a node touches an edge, the edge will be splited.

    3D only OOC with polytools

    TODO:
        * Crossing or Node/Edge intersections will NOT be
        recognized yet.
        * Edge on Node touch

    Parameters
    ----------
    plcs: [:gimliapi:`GIMLI::Mesh`]
        List of PLC that want to be merged into one new PLC

    tol : double
        Tolerance to check for duplicated nodes. [1e-3]

    Returns
    -------
    plc : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> from pygimli.viewer.mpl import drawMesh
    >>> world = mt.createWorld(start=[-10, 0], end=[10, -10], marker=1)
    >>> c1 = mt.createCircle([-1, -4], radius=1.5, area=0.1,
    ...                       marker=2, nSegments=5)
    >>> c2 = mt.createCircle([-6, -5], radius=[1.5, 3.5], isHole=1)
    >>> r1 = mt.createRectangle(pos=[3, -5], size=[2, 2], marker=3)
    >>> r2 = mt.createRectangle(start=[4, -4], end=[6, -6],
    ...                          marker=4, area=0.1)
    >>> plc = mt.mergePLC([world, c1, c2, r1, r2])
    >>> fig, ax = pg.plt.subplots()
    >>> drawMesh(ax, plc)
    >>> drawMesh(ax, mt.createMesh(plc))
    >>> pg.wait()
    """
    if plcs[0].dim() == 3:
        return mergePLC3D(plcs, tol)

        # tmp = pg.optImport('tempfile')
        # names = []
        # for p in plcs:
        #     _, namePLC = tmp.mkstemp(suffix='.poly')
        #     pg.meshtools.exportPLC(p, namePLC)
        #     names.append(namePLC)

        # for n in names[1:]:
        #     syscal = 'polyMerge {0} {1} {0}'.format(names[0], n)
        #     pg.debug(syscal)
        #     os.system(syscal)

        # plc = readPLC(names[0])

        # for n in names:
        #     try:
        #         pg.debug('Remove:', n)
        #         os.remove(n)
        #     except Exception:
        #         print("can't remove:", n)

        # return plc

    # handle 2D geometries
    plc = pg.Mesh(dim=2, isGeometry=True)

    for p in plcs:
        nodes = []
        for n in p.nodes():
            nn = plc.createNodeWithCheck(n.pos(), tol,
                                         warn=False, edgeCheck=True)
            if n.marker() != 0:
                nn.setMarker(n.marker())
            nodes.append(nn)

        for e in p.boundaries():
            plc.createEdge(nodes[e.node(0).id()], nodes[e.node(1).id()],
                           e.marker())

        if len(p.regionMarkers()) > 0:
            for rm in p.regionMarkers():
                plc.addRegionMarker(rm)

        if len(p.holeMarker()) > 0:
            for hm in p.holeMarker():
                plc.addHoleMarker(hm)

    return plc


def mergePLC3D(plcs, tol=1e-3):
    """Merge a list of 3D PLC into one.

    Experimental replacement for polyMerge. Don't expect too much.

    Works if:
        * all plcs are free and does not have any contact to each other
        * contact of two facets if the second is completely within the first
    """
    if len(plcs) < 2:
        pg.critical("Give at least 2 PLCs.")

    if plcs[0].dim() != 3:
        pg.warn("2D poly found. redirect to mergePLC")
        return mergePLC(plcs, tol)

    p0 = pg.Mesh(plcs[0])

    for p in plcs[1:]:
        for b in p.boundaries():
            p0.copyBoundary(b)

        if len(p.regionMarkers()) > 0:
            for rm in p.regionMarkers():
                p0.addRegionMarker(rm)

        for hm in p.holeMarker():
            p0.addHoleMarker(hm)

    return p0


def createParaDomain2D(*args, **kwargs):
    """API change here .. use createParaMeshPLC instead."""
    pg.deprecated("use createParaMeshPLC")
    return createParaMeshPLC(*args, **kwargs)


def createParaMeshPLC(sensors, paraDX=1, paraDepth=-1, paraBoundary=2,
                      paraMaxCellSize=0.0, boundary=-1, boundaryMaxCellSize=0,
                      balanceDepth=True,
                      isClosed=False, addNodes=1, **kwargs):
    """Create a geometry (PLC) for an inversion parameter mesh.

    Create an inversion mesh geometry (PLC) for a given list of
    sensor positions. Sensor positions are assumed to be on the surface and
    must be unique and sorted along x coordinate.

    You can create a parameter mesh without sensors if you just set [xMin,
    xMax] as sensors.

    The PLC is a :gimliapi:`GIMLI::Mesh` and contain nodes, edges and two
    region markers, one for the parameters domain (marker=2) and a larger
    boundary around the outside (marker=1).

    TODO:
        * additional topographic points
        * spline interpolation between sensorpoints or addpoints for non closed
        * subsurface sensors (partly .. see example)

    Parameters
    ----------
    sensors : [RVector3] | DataContainer with sensorPositions() | [xMin, xMax]
        Sensor positions. Must be sorted and unique in positive x direction.
        Depth need to be y-coordinate.

    paraDX : float [1]
        Relative distance for refinement nodes between two sensors (1=none),
        e.g., 0.5 means 1 additional node between two neighboring sensors
        e.g., 0.33 means 2 additional equidistant nodes between two sensors

    paraDepth : float[-1], optional
        Maximum depth in m for parametric domain.
        Automatic (<=0) results in 0.4 * maximum sensor span range in m

    balanceDepth: bool [True]
        Equal depth for the parametric domain.

    paraBoundary : float, optional
        Margin for parameter domain in absolute sensor distances. 2 (default).

    paraMaxCellSize: double, optional
        Maximum cell size for parametric region in m²

    boundaryMaxCellSize: double, optional
        Maximum cells size in the boundary region in m²

    boundary : float, optional
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lover 0 force the boundary to be 4 times para domain width.

    isClosed : bool [False]
        Create a closed geometry from sensor positions.
        Region marker is 1. Boundary marker is -1 (homogeneous Neumann)

    addNodes : int [1]
        Number of additional nodes to be added equidistant between sensors.

    trapRatio : float [0]
        Form a trapezoidal shape instead of a rectangle.
        The value is a ratio of the total length to put inside at depth.

    Returns
    -------
    poly: :gimliapi:`GIMLI::Mesh`
        Piecewise linear complex (PLC) containing nodes and edges

    Examples
    --------
    >>> # no need to import matplotlib, pygimli show does.
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> # Create the simplest paramesh PLC with a para box of 10 m without
    >>> # sensors
    >>> p = mt.createParaMeshPLC([0,10])
    >>> # you can add subsurface sensors now with
    >>> for z in range(1,4):
    ...     n = p.createNode((5,-z), -99)
    >>> ax,_ = pg.show(p)
    """
    if isClosed:
        plc = createPolygon(sensors, isClosed=True, addNodes=addNodes,
                            boundaryMarker=-1, marker=1,
                            area=paraMaxCellSize, **kwargs)
        return plc

    noSensors = False
    if hasattr(sensors, 'sensorPositions'):  # obviously a DataContainer type
        sensors = sensors.sensorPositions()
    elif isinstance(sensors, np.ndarray):
        if sensors.ndim == 1:
            sensors = [pg.RVector3(s, 0) for s in sensors]
        else:  # assume 2d array with 2 or 3 values per item
            sensors = [pg.RVector3(s) for s in sensors]
    elif isinstance(sensors, list):
        if len(sensors) == 2:
            # guess we have just a desired Pbox with
            sensors = [pg.RVector3(sensors[0], 0.0),
                       pg.RVector3(sensors[1], 0.0)]
            noSensors = True
            paraBoundary = 0

    eSpacing = kwargs.pop('eSpacing', sensors[0].distance(sensors[1]))

    iz = 1
    xMin, yMin, zMin = sensors[0][0], sensors[0][1], sensors[0][2]
    xMax, yMax, zMax = xMin, yMin, zMin
    for e in sensors:
        xMin = min(xMin, e[0])
        xMax = max(xMax, e[0])
        yMin = min(yMin, e[1])
        yMax = max(yMax, e[1])
        zMin = min(zMin, e[2])
        zMax = max(zMax, e[2])

    if abs(yMin) < 1e-8 and abs(yMax) < 1e-8:
        iz = 2

    paraBound = eSpacing * paraBoundary

    if paraDepth <= 0.0:
        paraDepth = 0.4 * (xMax - xMin)

    poly = pg.Mesh(dim=2, isGeometry=True)
    # define para domain without surface
    n1 = poly.createNode([xMin - paraBound, sensors[0][iz]])
    xStart = xMin - paraBound
    xEnd = xMax + paraBound
    trapRatio = np.minimum(kwargs.pop("trapRatio", 0), 0.45)
    dxTrap = (xEnd - xStart) * trapRatio

    if balanceDepth:
        bD = min(sensors[0][iz] - paraDepth, sensors[-1][iz] - paraDepth)
        n2 = poly.createNode([xStart + dxTrap, bD])
        n3 = poly.createNode([xEnd - dxTrap, bD])
    else:
        n2 = poly.createNode([xStart + dxTrap, sensors[0][iz] - paraDepth])
        n3 = poly.createNode([xEnd - dxTrap, sensors[-1][iz] - paraDepth])

    n4 = poly.createNode([xMax + paraBound, sensors[-1][iz]])

    if boundary < 0:
        boundary = 4

    bound = abs(xMax - xMin) * boundary
    if bound > paraBound:
        # define world without surface
        n11 = poly.createNode(n1.pos() - [bound, 0.])
        n12 = poly.createNode(n11.pos() - [0., bound + paraDepth])
        n14 = poly.createNode(n4.pos() + [bound, 0.])
        n13 = poly.createNode(n14.pos() - [0., bound + paraDepth])

        poly.createEdge(n1, n11, pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
        poly.createEdge(n11, n12, pg.core.MARKER_BOUND_MIXED)
        poly.createEdge(n12, n13, pg.core.MARKER_BOUND_MIXED)
        poly.createEdge(n13, n14, pg.core.MARKER_BOUND_MIXED)
        poly.createEdge(n14, n4, pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
        poly.addRegionMarker(n12.pos() + [1e-3, 1e-3], 1, boundaryMaxCellSize)

    poly.createEdge(n1, n2, 1)
    poly.createEdge(n2, n3, 1)
    poly.createEdge(n3, n4, 1)
    poly.addRegionMarker(n2.pos() + [1e-3, 1e-3], 2, paraMaxCellSize)

    # define surface
    nSurface = []
    nSurface.append(n1)
    if paraDX == 0.0:
        paraDX = 1.0

    if not noSensors:
        for i, e in enumerate(sensors):
            if iz == 2:
                e.rotateX(-math.pi / 2)

            # nSurface.append(poly.createNode(e, pg.core.MARKER_NODE_SENSOR))
            if addNodes > 1:
                nSurface.append(poly.createNode(e, pg.core.MARKER_NODE_SENSOR))
                if i < len(sensors) - 1:
                    e1 = sensors[i + 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)

                    for j in range(addNodes):
                        nSurface.append(poly.createNode(
                            e + (e1 - e) * (j+1)/(addNodes+1)))
            elif paraDX >= 0.5:
                nSurface.append(poly.createNode(e, pg.core.MARKER_NODE_SENSOR))
                if i < len(sensors) - 1:
                    e1 = sensors[i + 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)
                    nSurface.append(poly.createNode((e + e1) * 0.5))

            elif paraDX < 0.5:
                if i > 0:
                    e1 = sensors[i - 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)

                    nSurface.append(poly.createNode(e - (e - e1) * paraDX))

                nSurface.append(poly.createNode(e, pg.core.MARKER_NODE_SENSOR))
                if i < len(sensors) - 1:
                    e1 = sensors[i + 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)

                    nSurface.append(poly.createNode(e + (e1 - e) * paraDX))

    nSurface.append(n4)

    ### remove duplicates
    nSurface = [ni for i, ni in enumerate(nSurface) if ni not in nSurface[:i]]
    ### sort again .. just in case
    nSurface.sort(key=lambda n: n.pos()[0])

    for i in range(len(nSurface) - 1, 0, -1):
        poly.createEdge(nSurface[i], nSurface[i - 1],
                        pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

    return poly


def createParaMeshSurface(sensors, paraBoundary=None, boundary=-1,
                          surfaceMeshQuality=30, surfaceMeshArea=0,
                          addTopo=None):
    r"""Create surface mesh for an 3D inversion parameter mesh.

        Topographic information (non-zero z-coodinate) can be from sensors
        together with addTopo, or in addTopo alone if provided.
        Outside boundary corners are set to median of all topography.

    Args
    ----
    sensors: DataContainer with sensorPositions()
        Sensor positions.

    paraBoundary: [float, float] [1.1, 1.1]
        Margin for parameter domain in relative extend.

    boundary: [float, float] [10., 10.]
        Boundary width to be appended for domain prolongation in relative
        para domain size.

    surfaceMeshQuality: float [30]
        Quality of the surface mesh.

    surfaceMeshArea: float [0]
        Max cell Size for parametric domain.

    addTopo: [[x,y,z],]
        Number of additional nodes for topography.

    Returns
    -------
    surface: :gimliapi:`GIMLI::Mesh`
        3D Surface mesh

    Examples
    --------
    >>> # no need to import matplotlib, pygimli show does.
    >>> import numpy as np
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> # very simple design: 10 sensors on 1D profile in 3D topography
    >>> x = np.linspace(-10, 10, 10)
    >>> topo = [[15, -15, 10], [-15, 15, -10]]
    >>> surface = mt.createParaMeshSurface(np.asarray([x, x, x*0]).T,
    ...                                    paraBoundary=[1.2, 1.2],
    ...                                    boundary=[2, 2],
    ...                                    surfaceMeshQuality=30,
    ...                                    addTopo=topo)
    >>> _ = pg.show(surface, showMesh=True, color='white')
    """
    if hasattr(sensors, 'sensors'):
        sensors = sensors.sensors()

    sensors = np.asarray(sensors)

    if paraBoundary is None:
        paraBoundary = [1.1, 1.1]
    elif isinstance(paraBoundary, (int, float)):
        paraBoundary = [paraBoundary, paraBoundary]

    if boundary is None:
        boundary = [10.0, 10.0]

    # find maximum extent
    boundaryRect = pg.meshtools.createRectangle(pnts=sensors[:, 0:2],
                                                minBB=False,
                                                minBBOffset=boundary)
    for i in range(4):
        boundaryRect.boundary(i).setMarker(i+1)
        boundaryRect.node(i).setMarker((i % 4 + 1) * 10)

    # collect all pnts with topography
    if addTopo is not None:
        if min(sensors[:, 2]) != max(sensors[:, 2]) and sensors[0][2] != 0.0:
            pnts = np.vstack((sensors, addTopo))
        else:
            pnts = np.asarray(addTopo)
    else:
        pnts = np.array(sensors)

    # add maximal extent corners to median topo
    boundaryRect.translate([0, 0, np.median(pnts[:, 2])])
    pnts = np.vstack((boundaryRect.positions(), pnts))

    # create mesh for topo interpolation
    pntsSurface = pg.meshtools.createMesh(pnts[:, 0:2])

    # find parameter extent
    paraRect = pg.meshtools.createRectangle(pnts=sensors[:, 0:2],
                                            minBB=True,
                                            minBBOffset=paraBoundary,
                                            area=surfaceMeshArea,
                                            )
    for i in range(4):
        paraRect.boundary(i).setMarker(i + 5)
        paraRect.node(i).setMarker((i % 4 + 5) * 10)

    # create surface mesh of sensors and with maximal and parameter extent
    if boundary == 0:
        surfacePLC = paraRect + sensors[:, 0:2]
    else:
        surfacePLC = boundaryRect + paraRect + sensors[:, 0:2]

    surface = pg.meshtools.createMesh(surfacePLC, quality=surfaceMeshQuality)

    # interpolate topography to surface
    # ax, _ = pg.show(pntsSurface)
    # p_ = np.asarray(list(zip(pg.x(surface.positions()),
    #                          pg.y(surface.positions()))))
    # ax.scatter(*p_.T)
    sZ = pg.interpolate(pntsSurface, pnts[:, 2], surface.positions())
    for n in surface.nodes():
        n.translate(0, 0, sZ[n.id()])

    # create 3D surfacemesh
    s = pg.meshtools.createSurface(
        surface, boundaryMarker=pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

    return s
    # pg.show(surface, showMesh=True)


def createParaMeshPLC3D(sensors, paraDX=0, paraDepth=-1, paraBoundary=None,
                        paraMaxCellSize=0.0, boundary=None,
                        boundaryMaxCellSize=0,
                        surfaceMeshQuality=30, surfaceMeshArea=0,
                        addTopo=None, isClosed=False, **kwargs):
    r"""Create a geometry (PLC) for an 3D inversion parameter mesh.

    Todo
    ----
        * 3D without TOPO
        * better paraBoundary scale for with and height

    Args
    ----
    sensors: Sensor list or pg.DataContainer with .sensors()
        Sensor positions.

    paraDX : float [1]
        Absolute distance for node refinement (0=none).
        Refinement node will be placed below the surface.

    paraDepth : float[-1], optional
        Maximum depth in m for parametric domain.
        Automatic (<=0) results in 0.4 * maximum sensor span range in m.
        Depth is set to median sensors depth + paraDepth.

    paraBoundary: [float, float] | float [1.1, 1.1]
        Margin for parameter domain in relative extend.

    paraMaxCellSize: double, optional
        Maximum cell size for parametric region in m³

    boundaryMaxCellSize: double, optional
        Maximum cells size in the boundary region in m³

    boundary: [float, float] | float [10., 10.]
        Boundary width to be appended for domain prolongation in relative
        para domain size.

    surfaceMeshQuality: float [30]
        Quality of the surface mesh.

    surfaceMeshArea: float [0]
        Max boundary size for surface area in parametric region.

    addTopo: [[x,y,z],]
        Number of additional nodes for topography.

    Returns
    -------
    poly: :gimliapi:`GIMLI::Mesh`
        Piecewise linear complex (PLC) containing nodes and edges
    """
    if hasattr(sensors, 'sensors'):
        sensors = sensors.sensors()

    sensors = np.asarray(sensors)

    if boundary is None:
        boundary = [10.0, 10.0]
    elif isinstance(boundary, (float, int)):
        boundary = [boundary, boundary]

    surface = pg.meshtools.createParaMeshSurface(
        sensors, paraBoundary=paraBoundary, boundary=boundary,
        surfaceMeshQuality=surfaceMeshQuality,
        surfaceMeshArea=surfaceMeshArea,
        addTopo=addTopo)

    # find depth and paradepth
    xSpan = (max(sensors[:, 0]) - min(sensors[:, 0]))
    ySpan = (max(sensors[:, 1]) - min(sensors[:, 1]))

    if paraDepth == -1:
        paraDepth = (0.4*(max(xSpan, ySpan)))

    paraDepth = np.median(sensors[:, 2]) - paraDepth
    depth = paraDepth - max(boundary[0]*xSpan, boundary[1]*ySpan)/2

    def sortP(p):
        base = pg.core.Line(p[0], p[1]).at(-1e7)

        def cmp_(p1, p2):
            if p1.distSquared(base) < p2.distSquared(base):
                return -1
            else:
                return 1

        p.sort(key=functools.cmp_to_key(cmp_))

    bounds = [surface]
    # close outer surfaces
    bttm = []
    for i in range(4):

        p = [n.pos() for n in surface.nodes() if n.marker() == i+1]
        p.append(surface.nodes(
            surface.nodeMarkers() == (i % 4 + 1) * 10)[0].pos())
        p.append(surface.nodes(
            surface.nodeMarkers() == ((i + 1) % 4 + 1) * 10)[0].pos())
        sortP(p)

        p0 = pg.Pos(p[-1])
        p0[2] = depth
        p.append(p0)

        p1 = pg.Pos(p[0])
        p1[2] = depth
        p.append(p1)

        m = pg.meshtools.createPolygon(p, isClosed=True)
        f = pg.meshtools.createFacet(m, boundaryMarker=-2)

        bounds.append(f)
        bttm.append(p0)
        bttm.append(p1)

    m = pg.meshtools.createRectangle(pnts=bttm, minBB=True)
    m.translate([0, 0, depth])
    bttmA = pg.meshtools.createFacet(m, boundaryMarker=-2)
    bounds.append(bttmA)

    # close para surfaces
    bttm = []
    for i in range(4):
        p = [n.pos() for n in surface.nodes() if n.marker() == i+5]
        p.append(surface.nodes(
            surface.nodeMarkers() == (i % 4 + 5) * 10)[0].pos())
        p.append(surface.nodes(
            surface.nodeMarkers() == ((i + 1) % 4 + 5) * 10)[0].pos())
        sortP(p)

        p0 = pg.Pos(p[-1])
        p0[2] = paraDepth
        p.append(p0)

        p1 = pg.Pos(p[0])
        p1[2] = paraDepth
        p.append(p1)

        m = pg.meshtools.createPolygon(p[::-1], isClosed=True)
        f = pg.meshtools.createFacet(m, boundaryMarker=1)

        bounds.append(f)
        bttm.append(p0)
        bttm.append(p1)

    m = pg.meshtools.createRectangle(pnts=bttm, minBB=True)
    m.translate([0, 0, paraDepth])
    bttmP = pg.meshtools.createFacet(m, boundaryMarker=1)
    bounds.append(bttmP)

    pdPLC = pg.meshtools.mergePLC(bounds)

    if paraDX > 0:
        for s in sensors:
            pdPLC.createNode(s - [0.0, 0.0, paraDX])

    pdPLC.addRegionMarker(pg.center(bttmA.positions()) + [0.0, 0.0, 0.1],
                          marker=1, area=boundaryMaxCellSize)
    pdPLC.addRegionMarker(pg.center(bttmP.positions()) + [0.0, 0.0, 0.1],
                          marker=2, area=paraMaxCellSize)

    return pdPLC


def readPLC(filename, comment='#'):
    r"""Read in a piece-wise linear complex object (PLC) from .poly file.

    A PLC is a pyGIMLi geometry, e.g., created using `mt.exportPLC`.

    Read 2D :term:`Triangle` or 3D :term:`Tetgen` PLC files.

    Parameters
    ----------
    filename: string
        Filename *.poly

    comment: string ('#')
        String containing all characters that define a comment line.
        Identified lines will be ignored during import.

    Returns
    -------
    poly :
        :gimliapi:`GIMLI::Mesh`

    See Also
    --------
    exportPLC
    """
    with open(filename, 'r') as fi:
        content = fi.readlines()

    # Filter comment lines
    comment_lines = []
    for i, line in enumerate(content):
        if line[0] in comment:
            comment_lines.append(i)
    for j in comment_lines[::-1]:
        del(content[j])

    # Read header
    headerLine = content[0].split('\r\n')[0].split()

    if len(headerLine) != 4:
        raise Exception("Format unknown! header size != 4", headerLine)

    fromOne = 0
    nVerts = int(headerLine[0])
    dimension = int(headerLine[1])
    nPointsAttributes = int(headerLine[2])
    haveNodeMarker = int(headerLine[3])

    poly = pg.Mesh(dim=dimension, isGeometry=False)
    # isGeometry forces expensive checks: we assume the plc is valid so we set
    # this flag in the end

    # Nodes section
    for i in range(nVerts):
        row = content[1 + i].split('\r\n')[0].split()

        if len(row) == (1 + dimension + nPointsAttributes + haveNodeMarker):
            if i == 0:
                fromOne = int(row[0])
            if dimension == 2:
                n = poly.createNode((float(row[1]), float(row[2])))
            elif dimension == 3:
                n = poly.createNode((float(row[1]), float(row[2]),
                                     float(row[3])))
            if haveNodeMarker:
                n.setMarker(int(row[-1]))

        else:
            print(i, len(row), row,
                  (1 + dimension + nPointsAttributes + haveNodeMarker))
            raise Exception("Poly file seams corrupt: node section line: " +
                            content[1 + i])

    # Segment section
    row = content[1 + nVerts].split()

    if len(row) != 2:
        raise Exception("Format unknown for segment section " + row)

    nSegments = int(row[0])
    haveBoundaryMarker = int(row[1])

    if dimension == 2:
        for i in range(nSegments):
            row = content[2 + nVerts + i].split()

            if len(row) == (3 + haveBoundaryMarker):
                marker = 0
                if haveBoundaryMarker:
                    marker = int(row[3])

                poly.createEdge(
                    poly.node(int(row[1]) - fromOne),
                    poly.node(int(row[2]) - fromOne), marker)
    else:
        segment_offset = 0
        for i in range(nSegments):
            row = content[2 + nVerts + i + segment_offset].split()
            numBounds = int(row[0])
            numHoles = int(row[1])
            # if numHoles != '0':
            #     pg.error("Can't handle 3D faces with holes yet")
            marker = 0
            if haveBoundaryMarker:
                marker = int(row[2])

            face = None
            for k in range(numBounds):
                boundRow = content[2 + nVerts + i + segment_offset + 1]\
                    .split()
                # nNodes = int(boundRow[0])
                nodeIdx = [int(_b) for _b in boundRow[1:]]

                if k == 0:
                    face = poly.createPolygonFace(poly.nodes(nodeIdx),
                                                  marker=marker, check=True)
                else:
                    if len(nodeIdx) == 2:
                        if nodeIdx[0] == nodeIdx[1]:
                            face.addSecondaryNode(poly.node(nodeIdx[0]))
                    else:
                        face.addSubface(nodeIdx)

                segment_offset += 1

            for k in range(numHoles):
                r = content[2 + nVerts + i + segment_offset + 1]\
                    .split()
                face.addHoleMarker([float(hm) for hm in r[1:]])

                segment_offset += 1

        nSegments += segment_offset

    # Hole section
    row = content[2 + nVerts + nSegments].split()

    if len(row) != 1:
        raise Exception("Format unknown for hole section " + row)

    nHoles = int(row[0])
    for i in range(nHoles):
        row = content[3 + nVerts + nSegments + i].split()
        if len(row) == 3:
            poly.addHoleMarker([float(row[1]), float(row[2])])
        elif len(row) == 4 and dimension == 3:
            poly.addHoleMarker([float(row[1]), float(row[2]), float(row[3])])
        else:
            raise Exception("Poly file seams corrupt: hole section line (3):" +
                            row + " : " + str(i) + " " + str(len(row)))

    if (3 + nVerts + nSegments + nHoles) < len(content):
        # Region section
        row = content[3 + nVerts + nSegments + nHoles].split()

        if len(row) != 1:
            raise Exception("Format unknown for region section " + row)

        nRegions = int(row[0])

        for i in range(nRegions):
            row = content[4 + nVerts + nSegments + nHoles + i].split()
            if len(row) == 5:
                poly.addRegionMarker([float(row[1]), float(row[2])],
                                     marker=int(float(row[3])),
                                     area=float(row[4]))
            elif len(row) == 6 and dimension == 3:
                poly.addRegionMarker([float(row[1]), float(row[2]),
                                      float(row[3])],
                                     marker=int(float(row[4])),
                                     area=float(row[5]))
            else:
                raise Exception("Poly file seams corrupt: region section " +
                                "line (5): " + str(i) + " " + str(len(row)))

    poly.setGeometry(True)
    return poly


def exportPLC(poly, fname, **kwargs):
    r"""Export a piece-wise linear complex (PLC) to a .poly file (2D or 3D).

    Chooses from poly.dimension() and forwards accordingly to
    :gimliapi:`GIMLI::Mesh::exportAsTetgenPolyFile`
    or :py:mod:`pygimli.meshtools.writeTrianglePoly`

    Parameters
    ----------
    poly : :gimliapi:`GIMLI::Mesh`
        The polygon to be written.

    fname : string
        Filename of the file to write (\\*.n, \\*.e).

    Examples
    --------
    >>> import pygimli as pg
    >>> import tempfile, os
    >>> fname = tempfile.mktemp() + '.poly' # Create temporary filename.
    >>> world2d = pg.meshtools.createWorld(start=[-10, 0], end=[20, -10])
    >>> pg.meshtools.exportPLC(world2d, fname)
    >>> read2d = pg.meshtools.readPLC(fname)
    >>> print(read2d)
    Mesh: Nodes: 4 Cells: 0 Boundaries: 4
    >>> world3d = pg.createGrid([0, 1], [0, 1], [-1, 0])
    >>> pg.meshtools.exportPLC(world3d, fname)
    >>> os.remove(fname)

    See Also
    --------
    readPLC
    """
    if poly.dimension() == 2:
        exportTrianglePoly(poly, fname, **kwargs)
    else:
        exportTetgenPoly(poly, fname, **kwargs)


def exportTrianglePoly(poly, fname, float_format='.15e'):
    r"""Write :term:`Triangle` poly.

    Write :term:`Triangle` :cite:`Shewchuk96b` ASCII file.
    See: ://www.cs.cmu.edu/~quake/triangle.html

    Parameters
    ----------
    poly : :gimliapi:`GIMLI::Mesh`
        mesh PLC holding nodes, edges, holes & regions

    fname : string
        Target filename *.poly

    float_format : string
        format string for floats according to str.format()

    verbose : boolean [False]
        Be verbose during import.
    """
    if fname.rfind('.poly') == -1:
        fname = fname + '.poly'

    if float_format[0] != '{':
        pfmt = '{:' + float_format + '}'
    else:
        pfmt = float_format
    with open(fname, 'w') as fid:
        fid.write('{:d}\t2\t0\t1\n'.format(poly.nodeCount()))
        nm = poly.nodeMarkers()
        bm = poly.boundaryMarkers()

        fmt = '{:d}' + ('\t' + pfmt) * 2 + '\t{:d}\n'
        for i, p in enumerate(poly.positions()):
            fid.write(fmt.format(i, p.x(), p.y(), nm[i]))
        fid.write('{:d}\t1\n'.format(poly.boundaryCount()))

        for i, b in enumerate(poly.boundaries()):
            fid.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(i, b.node(0).id(),
                                                        b.node(1).id(), bm[i]))
        fid.write('{:d}\n'.format(len(poly.holeMarker())))

        fmt = '{:d}' + ('\t' + pfmt) * 2 + '\n'
        for i, h in enumerate(poly.holeMarker()):
            fid.write(fmt.format(i, h.x(), h.y()))
        fid.write('{:d}\n'.format(len(poly.regionMarkers())))

        fmt = '{:d}' + ('\t' + pfmt) * 3 + '\t{:.15e}\n'
        for i, r in enumerate(poly.regionMarkers()):
            fid.write(fmt.format(i, r.x(), r.y(), r.marker(), r.area()))

    return


def writeTrianglePoly(*args, **kwargs):
    """Backward compatibility.

    Please use :py:mod:`pygimli.meshtools.exportTrianglePoly`.
    """
    return exportTrianglePoly(*args, **kwargs)


def exportTetgenPoly(poly, filename, float_format='.12e', **kwargs):
    r"""Export PLC as tetgen poly file.

    Write given piecewise linear complex (mesh/poly) into Ascii file in
    :term:`Tetgen` .poly format.

    Parameters
    ----------
    filename: string
        Name in which the result will be written. The recommended file
        ending is '.poly'.

    poly: :gimliapi:`GIMLI::Mesh`
        Piecewise linear complex as :gimliapi:`GIMLI::Mesh` to be exported.

    float_format: format string ('.12e')
        Format that will be used to write float values in the Ascii file.
        Default is the exponential float form with a precision of 12 digits.

    kwargs:
        * extraBoundaries:
            Add additional polygons (#c42 still needed?)

    """
    if filename[-5:] != '.poly':
        filename = filename + '.poly'
    polytxt = ''
    sep = '\t'  # standard tab seperated file
    linesep = '\n'  # os.linesep does not work in mingwshell, testit!!
    assert poly.dim() == 3, 'Exit, only for 3D meshes.'
    boundary_marker = 1
    attribute_count = 0

    # Part 1/4: node list
    # intro line
    # <nodecount> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    polytxt += '{0}{5}{1}{5}{2}{5}{3}{4}'.format(poly.nodeCount(), 3,
                                                 attribute_count,
                                                 boundary_marker,
                                                 linesep, sep)
    # loop over positions, attributes and marker(node)
    # <point idx> <x> <y> <z> [attributes] [boundary marker]
    point_str = '{:d}'  # index of the point
    for i in range(3):
        # coords as float with given precision
        point_str += sep + '{:%s}' % (float_format)
    point_str += sep + '{:d}' + linesep  # node marker
    for j, node in enumerate(poly.nodes()):
        fill = [node.id()]
        fill.extend([pos for pos in node.pos()])
        fill.append(node.marker())
        polytxt += point_str.format(*fill)

    # Part 2/4: boundary list
    # intro line
    # <# of facets> <boundary markers (0 or 1)>
    nBoundaries = poly.boundaryCount()
    # look for extra boundaries present in either the PLC or in kwargs
    extraBoundaries = []
    if 'extraBoundaries' in kwargs:
        extraBoundaries += kwargs.pop('extraBoundaries', [])

    if hasattr(poly, 'extraBoundaries'):
        extraBoundaries += poly.extraBoundaries

    if len(extraBoundaries) > 0:
        print("Detected ", len(extraBoundaries), " extra boundaries!")

    nBoundaries += len(extraBoundaries)
    polytxt += '{0:d}{2}1{1}'.format(nBoundaries, linesep, sep)
    # loop over facets, each facet can contain an arbitrary number of holes
    # and polygons, in our case, there is always one polygon per facet.

    hole_str = '{:d}'
    for m in range(3):
        hole_str += sep + '{:%s}' % float_format

    hole_str += linesep

    for bound in poly.boundaries():
        # one line per facet
        # <# of polygons> [# of holes] [boundary marker]
        try:
            nSubs = bound.subfaceCount()
        except BaseException:
            nSubs = 0
        try:
            nHoles = len(bound.holeMarkers())
        except BaseException:
            nHoles = 0

        npolys = 1 + nSubs + len(bound.secondaryNodes())
        polytxt += '{3}{2}{4}{2}{0:d}{1}'.format(bound.marker(), linesep,
                                                 sep, npolys, nHoles)
        # inner loop over polygons
        # <# of corners> <corner 1> <corner 2> ... <corner #>
        for k in range(1):
            poly_str = '{:d}'.format(bound.nodeCount())
            poly_str += sep + sep.join(['{:d}'.format(n) for n in bound.ids()])
            polytxt += '{0}{1}'.format(poly_str, linesep)

        # loop over subfaces
        for k in range(nSubs):
            sub = bound.subface(k)
            poly_str = '{:d}'.format(len(sub))
            poly_str += sep + sep.join(['{:d}'.format(n.id()) for n in sub])
            polytxt += '{0}{1}'.format(poly_str, linesep)

        # inner loop over holes
        if nHoles > 0:
            for n, hole in enumerate(bound.holeMarkers()):
                polytxt += hole_str.format(n, *hole)

        # not necessary yet ?! why is there an extra hole section?
        # because this is for 2D holes in facets only

        # loop over secondaryNodes add them as single points
        for k in range(len(bound.secondaryNodes())):
            ind = bound.secondaryNodes()[k].id()
            poly_str = '{:d}'.format(2)
            poly_str += sep + '{0:d} {0:d}'.format(ind)
            polytxt += '{0}{1}'.format(poly_str, linesep)

    # part 2b: extra boundaries that cannot be part of mesh class
    for nodes in extraBoundaries:
        # <# of polygons> [# of holes] [boundary marker]
        npolys = 1
        polytxt += '1{2}0{2}{0:d}{1}'.format(111, linesep, sep)
        # <# of corners> <corner 1> <corner 2> ... <corner #>
        poly_str = '{:d}'.format(len(nodes))
        for ind in nodes:
            poly_str += sep + '{:d}'.format(ind)

        polytxt += '{0}{1}'.format(poly_str, linesep)

    # part 3/4: hole list
    # intro line
    # <# of holes>
    holes = poly.holeMarker()
    polytxt += '{:d}{}'.format(len(holes), linesep)
    # loop over hole markers
    # <hole #> <x> <y> <z>

    for n, hole in enumerate(holes):
        polytxt += hole_str.format(n, *hole)

    # part 4/4: region attributes and volume constraints (optional)
    # intro line
    # <# of regions>
    regions = poly.regionMarkers()
    polytxt += '{:d}{}'.format(len(regions), linesep)
    # loop over region markers
    # <region #> <x> <y> <z> <region number> <region attribute>
    region_str = '{:d}'
    for o in range(3):
        region_str += sep + '{:%s}' % (float_format)

    region_str += sep + '{:d}%s{:%s}' % (sep, float_format) + linesep
    for p, region in enumerate(regions):
        polytxt += region_str.format(p, region.x(), region.y(), region.z(),
                                     region.marker(),
                                     region.area())

    # writing file
    with open(filename, 'w') as out:
        out.write(polytxt)


def syscallTetgen(filename, quality=1.2, area=0, preserveBoundary=False,
                  verbose=False, tetgen='tetgen'):
    """Create a mesh from a PLC by system-calling :term:`Tetgen`.

    Create a :term:`Tetgen` :cite:`Si2004` mesh from a PLC.

    Parameters
    ----------
    filename: str

    quality: float [1.2]
        Refines mesh (to improve mesh quality). [1.1 ... ]

    area: float [0.0]
        Maximum cell size (m³)

    preserveBoundary: bool [False]
        Preserve PLC boundary mesh

    verbose: bool [False]
        be verbose

    tetgen: str | path ['tetgen']
        Binary for tetgen. Given as complete path or simple the binary name
        if its known in the system path.

    Returns
    -------
    mesh : :gimliapi:`GIMLI::Mesh`
    """
    filebody = filename.replace('.poly', '')

    # tetgen -pazVAC -q1.2 $MESH
    # test -O2
    syscal = tetgen + ' -pzAC'

    if area > 0:
        syscal += 'a' + str(area)
    else:
        syscal += 'a'

    syscal += 'q' + str(quality)

    if not verbose:
        syscal += 'Q'
    else:
        pass
        # syscal += 'V'

    if preserveBoundary:
        syscal += 'Y'

    syscal += ' ' + filebody + '.poly'

    if verbose:
        print(syscal)

    pg.debug(syscal)

    system(syscal)

    mesh = None
    if os.path.isfile(filebody + '.1.node'):
        # system('meshconvert -it -BD -o ' + filebody + ' ' + filebody + '.1')
        mesh = pg.meshtools.readTetgen(filebody + '.1')
        try:
            os.remove(filebody + '.1.node')
            os.remove(filebody + '.1.ele')
            os.remove(filebody + '.1.face')
        except BaseException as e:
            print(e)
    else:
        # system('meshconvert -it -BD -o ' + filebody + ' ' + filebody + '-1')
        mesh = pg.meshtools.readTetgen(filebody + '-1')
        try:
            os.remove(filebody + '-1.node')
            os.remove(filebody + '-1.ele')
            os.remove(filebody + '-1.face')
        except BaseException as e:
            print(e)

    # mesh = pg.Mesh(filebody)
    return mesh


def polyCreateWorld(filename, x=None, depth=None, y=None, marker=0,
                    maxCellSize=0, verbose=True):
    """Create the PLC of a default world.

    Out-of-core wrapper for dcfemlib::polytools::polyCreateWorld

    Todo
    * needs to be converted to the Python-only tools.

    Parameters
    ----------
    filename : str
        file name
    x : float
        x dimension
    y : float
        y dimension
    depth : float
        z dimension
    marker : int
        region marker
    maxCellSize : float
        maximum cell size
    verbose : bool
        be verbose
    """
    if depth is None:
        print("Please specify worlds depth.")
        return

    if x is None:
        print("Please specify worlds x dimension.")
        return

    dimension = 3
    z = depth

    if y is None:
        dimension = 2

    syscal = 'polyCreateWorld -d ' + str(dimension) \
        + ' -x ' + str(x) \
        + ' -y ' + str(y) \
        + ' -z ' + str(z) \
        + ' -m ' + str(marker) \

    if maxCellSize > 0:
        syscal += " -a " + str(maxCellSize)

    syscal = syscal + ' ' + filename

    if verbose:
        print(syscal)

    os.system(syscal)


def createSurface(mesh, boundaryMarker=None, verbose=True):
    """Convert a 2D mesh into a 3D surface mesh.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        The 2D input mesh.
    boundaryMarker: int[0]
        Boundary marker for the resulting faces.
        If None the cell markers of the mesh are taken.

    Returns
    -------
    :gimliapi:`GIMLI::Mesh`
        The 3D surface mesh.
    """
    if mesh.dimension() != 2:
        pg.error("Need two dimensional mesh")
    if mesh.cellCount() == 0:
        pg.error("Need a two dimensional mesh with cells")

    # think to use this
    # s = surface.createHull()

    surface = pg.Mesh(dim=3, isGeometry=True)

    [surface.createNode(n.pos(), n.marker()).id() for n in mesh.nodes()]

    for c in mesh.cells():
        surface.createBoundary(c.ids(), marker=c.marker())

    if boundaryMarker is not None:
        surface.setBoundaryMarkers(np.full(surface.boundaryCount(),
                                           boundaryMarker))

    return surface


def createFacet(mesh, boundaryMarker=None):
    """Create coplanar PLC from a 2d mesh or PLC.

    Parameters
    ----------
    mesh : pg.Mesh
        2D mesh or PLC to be converted to a 3D facet
    boundaryMarker : int
        boundary marker for the facet, otherwise taken from 2d region markers

    Returns
    -------
    plc : pyGIMLi mesh
        plc of the created facet

    Todo
    ----
    * mesh with cell into plc with boundaries
    * poly account for inner edges
    """
    if mesh.dimension() != 2:
        pg.error("need two dimensional mesh or poly")

    if mesh.cellCount() > 0:
        pg.critical("Implementme")

    poly = pg.Mesh(dim=3, isGeometry=True)

    nodes = [poly.createNode(n.pos()).id() for n in mesh.nodes()]

    if boundaryMarker is None:
        for rm in mesh.regionMarkers():
            boundaryMarker = rm.marker()
            continue

    b = poly.createBoundary(nodes, marker=boundaryMarker or 0)

    for h in mesh.holeMarker():
        b.addHoleMarker(h)

    return poly


def createCube(size=None, pos=None, start=None, end=None,
               rot=None, boundaryMarker=0, **kwargs):
    """Create cube PLC as geometrie definition.

    Create cube PLC as geometrie definition.
    You can either give size and center position or start and end position.

    Parameters
    ----------
    size: [x, y, z]
        x, y, and z-size of the cube. Default = [1.0, 1.0, 1.0] in m
    pos: [x, y, z]
        The center position, default is at the origin.
    start: [x, y, z]
        Left Front Bottom corner.
    end: [x, y, z]
        Right Back Top corner.
    rot: pg.Pos [None]
        Rotate on the center.
    boundaryMarker: int[0]
        Boundary marker for the resulting faces.

    ** kwargs:
        Marker related arguments:
        See :py:mod:`pygimli.meshtools.polytools.setPolyRegionMarker`

    Examples
    --------
    >>> import pygimli.meshtools as mt
    >>> cube = mt.createCube()
    >>> print(cube)
    Mesh: Nodes: 8 Cells: 0 Boundaries: 6
    >>> cube = mt.createCube([10, 10, 1])
    >>> print(cube.bb())
    [RVector3: (-5.0, -5.0, -0.5), RVector3: (5.0, 5.0, 0.5)]
    >>> cube = mt.createCube([10, 10, 1], pos=[-4.0, 0.0, 0.0])
    >>> print(pg.center(cube.positions()))
    RVector3: (-4.0, 0.0, 0.0)

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.
    """
    if size is None:
        size = [1.0, 1.0, 1.0]

    if start is not None and end is not None:
        size = pg.Pos(end) - pg.Pos(start)
        pos = pg.Pos(start) + pg.Pos(size)/2

    poly = pg.Mesh(3, isGeometry=True)

    for y in [-0.5, 0.5]:
        poly.createNode(-0.5, y, -0.5)
        poly.createNode(+0.5, y, -0.5)
        poly.createNode(+0.5, y, +0.5)
        poly.createNode(-0.5, y, +0.5)

    faces = [[4, 5, 1, 0],
             [5, 6, 2, 1],
             [6, 7, 3, 2],
             [7, 4, 0, 3],
             [0, 1, 2, 3],
             [7, 6, 5, 4], ]

    if isinstance(boundaryMarker, list):
        for i, f in enumerate(faces):
            poly.createPolygonFace(poly.nodes(f), marker=boundaryMarker[i])
    else:
        for f in faces:
            poly.createPolygonFace(poly.nodes(f), marker=boundaryMarker)

    poly.scale(size)

    if rot is not None:
        poly.rotate(rot)

    if pos is not None:
        poly.translate(pos)

    setPolyRegionMarker(poly, **kwargs)

    return poly


def createSphere(size=None, pos=None, nSegments=20, nRings=10,
                 rot=None, boundaryMarker=0, var='uvsphere',
                 triFaces=True, **kwargs):
    """Create sphere PLC as geometrie definition.

    Create sphere PLC as geometrie definition.
    You can either give size and center position.

    Arguments
    ---------
    size: [x, y, z]
        x, y, and z-size of the cube. Default = [1.0, 1.0, 1.0] in m
    pos: [x, y, z]
        The center position, default is at the origin.
    nSegments: int[0]
        Discrete number of segments along the horizontal axis
    nRings: int[0]
        Discrete number of segments along the vertical axis
    rot: pg.Pos [None]
        Rotate on the center.
    boundaryMarker: int[0]
        Boundary marker for the resulting faces.
    var: str='uvsphere'
        Variance of sphere: 'uvsphere', 'qsphere' or 'icosphere'.
    triFaces: bool [True]
        Output triangular faces instead of quadrilateral faces. Only for qsphere
        at the moment.

    Keyword Args
    ------------
    ** kwargs:
        * Marker related arguments:
        See :py:mod:`pygimli.meshtools.polytools.setPolyRegionMarker`
        * refine: int [3]
            Number of refinements to create a smoother sphere.
            Default `refine=3` for variance='qsphere' and `refine=2`
            for 'icosphere'

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> import pygimli.meshtools as mt
    >>> sphere = mt.createSphere()
    >>> print(sphere)
    Mesh: Nodes: 182 Cells: 0 Boundaries: 360
    >>> sphere = mt.createSphere()
    >>> print(sphere.bb())
    [RVector3: (-0.5, -0.5, -0.5), RVector3: (0.5, 0.5, 0.5)]
    >>> uvs = mt.createSphere([6, 5, 2], nSegments=25, nRings=15)
    >>> pg.show(uvs, showMesh=True, markers=True) # doctest: +ELLIPSIS
    (<pyvista...
    >>> uv = pg.meshtools.createSphere(var='uvsphere', pos=[3,0,0])
    >>> qs4 = pg.meshtools.createSphere(var='qsphere', pos=[2,0,0],
    ...                                 refine=3, triFaces=False)
    >>> qs3 = pg.meshtools.createSphere(var='qsphere', pos=[1,0,0],
    ...                                 refine=3, triFaces=True)
    >>> ico = pg.meshtools.createSphere(var='icosphere', pos=[0,0,0],
    ...                                 refine=2)
    >>> pg.show([uv, qs4, qs3, ico], showMesh=True, markers=True) # doctest: +ELLIPSIS
    (<pyvista...
    """
    if size is None:
        size = [1.0, 1.0, 1.0]

    poly = pg.Mesh(3, isGeometry=True)

    if var == 'uvsphere':
        if nSegments < 3:
            nSegments = 3

        for i in range(nRings+1):
            theta = np.pi * (i / nRings)
            for j in range(nSegments):
                phi = 2 * np.pi * (j / (nSegments))
                x = 0.5 * np.sin(theta) * np.cos(phi)
                y = 0.5 * np.sin(theta) * np.sin(phi)
                z = 0.5 * np.cos(theta)
                poly.createNode(x, y, z)
                if i == 0 or i == nRings:
                    break

        num = len(poly.nodes())

        for f in range(nSegments-1):
            poly.createBoundary(poly.nodes([0, 1+f, 2+f]),
                                marker=boundaryMarker)
            poly.createBoundary(poly.nodes([num-1, num-2-f, num-3-f]),
                                marker=boundaryMarker)

        poly.createBoundary(poly.nodes([0, 1, nSegments]),
                            marker=boundaryMarker)
        poly.createBoundary(poly.nodes([num-1, num-2, num-1-nSegments]),
                            marker=boundaryMarker)

        for b in range(nRings-2):
            b_l = b*nSegments
            for f in range(1, nSegments):
                poly.createBoundary(poly.nodes([f+b_l,
                                                f+b_l+nSegments,
                                                f+b_l+nSegments+1]),
                                    marker=boundaryMarker)
                poly.createBoundary(poly.nodes([num-1-f-b_l,
                                                num-1-f-b_l-nSegments,
                                                num-2-f-b_l-nSegments]),
                                    marker=boundaryMarker)

            poly.createBoundary(poly.nodes([b_l+nSegments,
                                            b_l+nSegments*2,
                                            b_l+nSegments+1]),
                                marker=boundaryMarker)
            poly.createBoundary(poly.nodes([num-1-b_l-nSegments,
                                            num-1-b_l-nSegments*2,
                                            num-2-b_l-nSegments]),
                                marker=boundaryMarker)
    elif var == 'qsphere':
        ## trivial version of quad sphere
        poly = createCube()
        refine = kwargs.pop('refine', 3)
        for r in range(refine):
            poly = poly.createH2()
            for n in poly.nodes():
                r = pg.Line([0.0, 0.0, 0.0], n.pos())
                # normalizing with project on the unit sphere
                n.setPos(r.at(0.5/n.pos().abs()))
                ### alternatively could be smooth with neighbors, and
                ### post normalize the end result, TODO test and compare

        if triFaces is True:
            poly = pg.meshtools.refineQuad2Tri(poly)

    elif var == 'icosphere':
        ## trivial version of quad sphere

        poly = pg.Mesh(3, isGeometry=True)
        # Create an icosahedron (20 faces, 12 vertices)
        phi = (1 + math.sqrt(5)) / 2  # golden ratio

        # Vertices of a unit icosahedron centered at origin
        verts = [
            [-1,  phi, 0], [ 1,  phi, 0], [-1, -phi, 0], [ 1, -phi, 0],
            [0, -1,  phi], [0,  1,  phi], [0, -1, -phi], [0,  1, -phi],
            [ phi, 0, -1], [ phi, 0,  1], [-phi, 0, -1], [-phi, 0,  1],]

        # Normalize to radius 0.5
        verts = np.array(verts)
        verts = 0.5 * verts / np.linalg.norm(verts[0])

        nodes = [poly.createNode(*v) for v in verts]

        # Faces of the icosahedron (each as 3 indices into verts)
        faces = [
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],]

        for f in faces:
            poly.createBoundary([nodes[i].id() for i in f], marker=boundaryMarker)

        refine = kwargs.pop('refine', 2)
        for r in range(refine):
            poly = poly.createH2()
            for n in poly.nodes():
                r = pg.Line([0.0, 0.0, 0.0], n.pos())
                # normalizing with project on the unit sphere
                n.setPos(r.at(0.5/n.pos().abs()))
                # alternative could be smooth with neighbors, and
                # post normalize with the end result, test and compare



    poly.scale(size)

    if rot is not None:
        poly.rotate(rot)

    if pos is not None:
        poly.translate(pos)

    setPolyRegionMarker(poly, **kwargs)

    return poly


def extrude(p2, z=-1.0, boundaryMarker=0, **kwargs):
    """Create 3D body by extruding a closed 2D poly into z direction.

    Parameters
    ----------
    p2 : :gimliapi:`GIMLI::Mesh`
        2D geometry

    z : float [-1.0]
        2D geometry

    Keyword Arguments
    -----------------
    ** kwargs:
        Marker related arguments:
        See :py:mod:`pygimli.meshtools.polytools.setPolyRegionMarker`

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.
    """
    if p2.dimension() != 2:
        pg.error("need two dimensional mesh or poly")

    if p2.cellCount() > 0:
        pg.critical("Implementme")

    poly = pg.Mesh(3, isGeometry=True)
    top = []
    for n in p2.nodes():
        top.append(poly.createNode(n.pos()).id())

    bot = []
    for n in p2.nodes():
        bot.append(poly.createNode(n.pos() + [0.0, 0.0, z]).id())
    N = len(top)

    poly.createPolygonFace(poly.nodes(top), marker=boundaryMarker)
    poly.createPolygonFace(poly.nodes(bot[::-1]), marker=boundaryMarker)

    for i in range(len(top)):
        poly.createPolygonFace(poly.nodes([i, N + i,
                                           N + (i + 1) % N,
                                           (i+1) % N]),
                               marker=boundaryMarker)

    setPolyRegionMarker(poly, **kwargs)

    return poly


def createCylinder(radius=1, height=1, nSegments=8,
                   pos=None, rot=None, boundaryMarker=0, **kwargs):
    """Create PLC of a cylinder.

    Out of core wrapper for dcfemlib::polytools.

    Note, there is a bug in the old polytools which ignores the area settings
    for marker == 0.

    Parameters
    ----------
    radius : float
        Radius of the cylinder.

    height : float
        Height of the cylinder

    nSegments : int [8]
        Number of segments of the cylinder.

    pos : pg.Pos [None]
        The center position, default is at the origin.

    Keyword Arguments
    ----------------
    ** kwargs:
        Marker related arguments:
        See :py:mod:`pygimli.meshtools.polytools.setPolyRegionMarker`

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.
    """
    circ = createCircle(radius=radius, nSegments=nSegments)
    poly = extrude(circ, z=height, boundaryMarker=boundaryMarker, **kwargs)
    # move it to z=0
    poly.translate([0.0, 0.0, -height/2])

    if rot is not None:
        c = pg.center(poly.positions())
        poly.translate(-c)
        poly.rotate(rot)
        poly.translate(c)

    if pos is not None:
        poly.translate(pos)

    return poly


def boundaryPlaneIntersectionLines(boundaries, plane):
    """Create Lines from boundaries that intersect a plane."""
    lines = []

    for b in boundaries:
        ps = []
        for i, n in enumerate(b.shape().nodes()):
            line = pg.Line(n.pos(), b.shape().node(
                (i + 1) % b.shape().nodeCount()).pos())
            p = plane.intersect(line, 1e-8, True)
            if p.valid():
                ps.append(p)

        if len(ps) == 2:
            lines.append(list(zip([ps[0].x(), ps[1].x()],
                                  [ps[0].z(), ps[1].z()])))
    return lines


if __name__ == "__main__":
    pass
