#!/usr/bin/env python3
"""Mesh based data transformation and mapping, interpolation, extrapolation."""

import numpy as np
import pygimli as pg


def nodeDataToCellData(mesh, data):
    """Convert node data to cell data.

    Convert node data to cell data via interpolation to cell centers.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    data : iterable [float]
        Data of len mesh.nodeCount().
        TODO complex, R3Vector, ndarray
    """
    if len(data) != mesh.nodeCount():
        raise BaseException(
            "Dimension mismatch, expecting nodeCount(): " +
            str(mesh.nodeCount()) + "got: " + str(len(data)),
            str(len(data[0])))

    return pg.interpolate(mesh, data, mesh.cellCenters())


def cellDataToNodeData(mesh, data, style='mean'):
    """Convert cell data to node data.

    Convert cell data to node data via non-weighted averaging (mean) of common
    cell data.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    data : iterable [float]
        Data of len mesh.cellCount().
        TODO complex, R3Vector, ndarray

    style : str ['mean']
        Interpolation style.
        * 'mean' : non-weighted averaging
        TODO harmonic averaging
        TODO weighted averaging (mean, harmonic)
        TODO interpolation via cell centered mesh

    Examples
    --------
    >>> import pygimli as pg
    >>> grid = pg.createGrid(x=(1,2,3),y=(1,2,3))
    >>> celldata = np.array([1, 2, 3, 4])
    >>> nodedata = pg.meshtools.cellDataToNodeData(grid, celldata)
    >>> print(nodedata.array())
    [1.  1.5 2.  2.  2.5 3.  3.  3.5 4. ]
    """
    if len(data) != mesh.cellCount():
        raise BaseException(
            "Dimension mismatch, expecting cellCount(): " +
            str(mesh.cellCount()) + "got: " + str(len(data)),
            str(len(data[0])))

    if style == 'mean':

        if np.ndim(data) == 1:
            return pg.core.cellDataToPointData(mesh, data)

        if mesh.dim() == 1:
            return pg.core.cellDataToPointData(mesh, data)
        elif mesh.dim() == 2:
            return np.array([
                pg.core.cellDataToPointData(mesh, data[:, 0]),
                pg.core.cellDataToPointData(mesh, data[:, 1])
            ]).T
        elif mesh.dim() == 3:
            return np.array([
                pg.core.cellDataToPointData(mesh, data[:, 0]),
                pg.core.cellDataToPointData(mesh, data[:, 1]),
                pg.core.cellDataToPointData(mesh, data[:, 2])
            ])
    else:
        raise BaseException("Style '" + style + "'not yet implemented."
                            "Currently styles available are: 'mean, '")


def nodeDataToBoundaryData(mesh, data):
    """Convert node data to boundary data.

    Assuming [NodeCount, dim] data
    DOCUMENT_ME
    """
    if len(data) != mesh.nodeCount():
        raise BaseException(
            "Dimension mismatch, expecting nodeCount(): " +
            str(mesh.nodeCount()) + " got: " + str(len(data)),
            str(len(data[0])))

    if isinstance(data, pg.PosVector):
        ret = np.zeros((mesh.boundaryCount(), 3))
        for b in mesh.boundaries():
            ret[b.id()] = sum(data[b.ids()]) / b.nodeCount()
        return ret

    dim = len(data[0])
    ret = np.zeros((mesh.boundaryCount(), dim))
    if dim == 1:
        for b in mesh.boundaries():
            ret[b.id()] = b.pot(b.center(), data)  # / b.nodeCount()
    elif dim == 2:
        for b in mesh.boundaries():
            if b.nodeCount() > 2:
                print(b)
                raise BaseException("implement me")

            ret[b.id()] = sum(data[b.ids()]) / b.nodeCount()
            continue
            # v = b.data(b.center(), data)
            # interpolation is hell slow here .. check!!!!!!
            v2 = (data[b.node(0).id()] + data[b.node(1).id()]) * 0.5
            # print(v -v2)
            ret[b.id()] = [v2[0], v2[1]]
    else:
        for b in mesh.boundaries():
            ret[b.id()] = sum(data[b.ids()]) / b.nodeCount()
            continue
            raise BaseException("don't use this until further checking")
            ret[b.id()] = b.vec(b.center(), data)  # / b.nodeCount()

    return ret


def cellDataToBoundaryData(mesh, data):
    """Convert cell data to boundary data."""
    if len(data) != mesh.cellCount():
        raise BaseException(
            "Dimension mismatch, expecting cellCount(): " +
            str(mesh.cellCount()) + "got: " + str(len(data)),
            str(len(data[0])))

    CtB = mesh.cellToBoundaryInterpolation()

    if isinstance(data, pg.PosVector()):
        return np.array([CtB * pg.x(data), CtB * pg.y(data),
                         CtB * pg.z(data)]).T
    else:
        return CtB * data


def fillEmptyToCellArray(mesh, vals, slope=True):
    """Prolongate empty cell values to complete cell attributes.

    It is possible to have zero values that are filled with appropriate
    attributes. This function tries to fill empty values successively by
    prolongation of the non-zeros.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        For each cell of mesh a value will be returned.

    vals : array
        Array of size cellCount().

    Returns
    -------
    atts : array
        Array of length mesh.cellCount()

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>>
    >>> # Create a mesh with 3 layers and an outer region for extrapolation
    >>> layers = mt.createWorld([0,-50],[100,0], layers=[-15,-35])
    >>> inner = mt.createMesh(layers, area=3)
    >>> mesh = mt.appendTriangleBoundary(inner, xbound=120, ybound=50,
    ...                                  area=20, marker=0)
    >>>
    >>> # Create data for the inner region only
    >>> layer_vals = [20,30,50]
    >>> data = np.array(layer_vals)[inner.cellMarkers() - 1]
    >>>
    >>> # The following fails since len(data) != mesh.cellCount(), extrapolate
    >>> # pg.show(mesh, data)
    >>>
    >>> # Create data vector, where zeros fill the outer region
    >>> data_with_outer = np.array([0] + layer_vals)[mesh.cellMarkers()]
    >>>
    >>> # Actual extrapolation
    >>> extrapolated_data = mt.fillEmptyToCellArray(mesh,
    ...                                  data_with_outer, slope=False)
    >>> extrapolated_data_with_slope = mt.fillEmptyToCellArray(mesh,
    ...                                 data_with_outer, slope=True)
    >>>
    >>> # Visualization
    >>> fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(10,8), sharey=True)
    >>> _ = pg.show(mesh, data_with_outer, ax=ax1, cMin=0)
    >>> _ = pg.show(mesh, extrapolated_data, ax=ax2, cMin=0)
    >>> _ = pg.show(mesh, extrapolated_data_with_slope, ax=ax3, cMin=0)
    >>> _ = ax1.set_title("Original data")
    >>> _ = ax2.set_title("Extrapolated with slope=False")
    >>> _ = ax3.set_title("Extrapolated with slope=True")
    """
    # atts = pg.Vector(mesh.cellCount(), 0.0)  # not used
    # oldAtts = mesh.cellAttributes()  # not used
    mesh.setCellAttributes(vals)
    mesh.createNeighborInfos()
    # std::vector< Cell * >
    # empties = []

    if slope:
        # search all cells with empty neighbors
        ids = pg.find(mesh.cellAttributes() != 0.0)

        for c in mesh.cells(ids):
            for i in range(c.neighborCellCount()):
                nc = c.neighborCell(i)

                if nc:
                    if nc.attribute() == 0.0:
                        # c.setAttribute(99999)

                        b = pg.core.findCommonBoundary(c, nc)
                        # search along a slope
                        pos = b.center() - b.norm() * 1000.
                        sf = pg.Vector()
                        startCell = c

                        while startCell:

                            startCell.shape().isInside(pos, sf, False)
                            nextC = startCell.neighborCell(sf)
                            if nextC:
                                if nextC.attribute() == 0.0:
                                    nextC.setAttribute(c.attribute())
                                else:
                                    break

                            startCell = nextC

    vals = mesh.cellAttributes()
    mesh.prolongateEmptyCellsValues(vals, background=-9e99)
    mesh.setCellAttributes(vals)
    return vals


def interpolateAlongCurve(curve, t, **kwargs):
    """Interpolate along curve.

    Return curve coordinates for a piecewise linear curve
    :math:`C(t) = {x_i,y_i,z_i}` at positions :math:`t`.
    Curve and :math:`t` values are expected to be sorted along distance from
    the origin of the curve.

    Parameters
    ----------
    curve : [[x,z]] | [[x,y,z]] | [:gimliapi:`GIMLI::Pos`,]
        Discrete curve for 2D :math:`x,z` curve=[[x,z]], 3D :math:`x,y,z`

    t: 1D iterable
        Query positions along the curve in absolute distance

    kwargs :
        If kwargs are given, an additional curve smoothing is applied using
        :py:mod:`pygimli.meshtools.interpolate`. The kwargs will be delegated.

        periodic : bool [False]
            Curve is periodic.
            Usefull for closed parametric spline interpolation.

    TODO
    ----
        * 'linear'

    Returns
    -------
    p : np.array
        Curve positions at query points :math:`t`.
        Dimension of p match the size of curve the coordinates.

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> fig, axs = pg.plt.subplots(1,3, figsize=(9,4))
    >>> topo = np.array([[-2., 0.], [-1., 0.], [0.5, 0.], [3., 8.], [4., 8.], [10., 1.], [13., 1.], [14., 1.]])
    >>> t = np.arange(26.0)
    >>> p = mt.interpolateAlongCurve(topo, t)
    >>> _= axs[0].plot(p[:,0], p[:,1], 'o', color='red', label='linear') #doctest: +ELLIPSIS
    >>>
    >>> p = mt.interpolateAlongCurve(topo, t, method='spline')
    >>> _= axs[1].plot(p[:,0], p[:,1], 'o', color='black', label='spline') #doctest: +ELLIPSIS
    >>>
    >>> p = mt.interpolateAlongCurve(topo, t, method='harmonic', nc=2)
    >>> _= axs[2].plot(p[:,0], p[:,1], 'o', color='green', label='harmonic') #doctest: +ELLIPSIS
    >>>
    >>> for a in axs:
    ...     _= a.plot(topo[:,0], topo[:,1], '-', lw=1, label='curve')
    ...     _= a.legend()
    ...     _= a.set_aspect(1)
    >>> pg.plt.show()
    """
    xC = np.zeros(len(curve))
    yC = np.zeros(len(curve))
    zC = np.zeros(len(curve))

    tCurve = kwargs.pop('tCurve', None)
    if tCurve is None:
        tCurve = pg.utils.cumDist(curve)
    dim = 3

    # extrapolate starting overlaps
    if min(t) < min(tCurve):
        d = pg.Pos(curve[1]) - pg.Pos(curve[0])
        # d[2] = 0.0
        d.normalise()
        curve = np.insert(curve, [0], [
            curve[0] - np.array(d * (min(tCurve) - min(t)))[0:curve.shape[1]]
        ], axis=0)
        tCurve = np.insert(tCurve, 0, min(t), axis=0)

    # extrapolate ending overlaps
    if max(t) > max(tCurve):
        d = pg.Pos(curve[-2]) - pg.Pos(curve[-1])
        # d[2] = 0.0
        d.normalise()
        curve = np.append(curve, [
            curve[-1] - np.array(d * (max(t) - max(tCurve)))[0:curve.shape[1]]
        ], axis=0)
        tCurve = np.append(tCurve, max(t))

    if isinstance(curve, pg.PosVector):
        xC = pg.x(curve)
        yC = pg.y(curve)
        zC = pg.z(curve)
    else:
        curve = np.array(curve)

        if curve.shape[1] == 2:
            xC = curve[:, 0]
            yC = curve[:, 1]
            dim = 2
        else:
            xC = curve[:, 0]
            yC = curve[:, 1]
            zC = curve[:, 2]

    if len(kwargs.keys()) > 0 and (kwargs.get('method', 'linear') != 'linear'):

        # interpolate more curve points to get a smooth line, guarantee to keep
        # original positions
        ti = np.array([np.linspace(tCurve[i], tCurve[i+1], 20)[:-1]
                       for i in range(len(tCurve)-1)]).flatten()
        ti = np.append(ti, tCurve[-1])

        xC = pg.interpolate(ti, tCurve, xC, **kwargs)
        yC = pg.interpolate(ti, tCurve, yC, **kwargs)

        if dim == 3:
            zC = pg.interpolate(ti, tCurve, zC, **kwargs)
        tCurve = ti

    xt = interpolate(t, tCurve, xC)
    yt = interpolate(t, tCurve, yC)

    if dim == 2:
        return np.vstack([xt, yt]).T

    zt = interpolate(t, tCurve, zC)

    return np.vstack([xt, yt, zt]).T


def tapeMeasureToCoordinates(tape, pos):
    """Interpolate 2D tape measured topography to 2D Cartesian coordinates.

    Tape and pos value are expected to be sorted along distance to the origin.

    DEPRECATED will be removed, use
    :py:mod:`pygimli.meshtools.interpolateAlongCurve` instead

    TODO optional smooth curve with harmfit
    TODO parametric
    TODO parametric + Topo: 3d

    Parameters
    ----------
    tape : [[x,z]] | [Pos] | R3Vector
        List of tape measured topography points with measured distance (x)
        from origin and height (z)

    pos : iterable
        Array of query positions along the tape measured profile t[0 ..

    Returns
    -------
    res : ndarray(N, 2)
        Same as pos but with interpolated height values.
        The Distance between pos points and res (along curve) points remains.

    """
    pg.deprecated("tapeMeasureToCoordinates", "interpolateAlongCurve")
    return interpolateAlongCurve(tape, pos)


def interpolate(*args, **kwargs):
    r"""Interpolate different kind of data.

    Convenience function for different interpolation schemes:

    **1. Interpolate mesh based data from one mesh to another**
      syntactic sugar for the core based `pg.core.interpolate`

    Parameters
    ----------
        args: :gimliapi:`GIMLI::Mesh`, :gimliapi:`GIMLI::Mesh`, values
            `outData = pg.core.interpolate(outMesh, inMesh, values)`

            Interpolated values based on inMesh to outMesh.

            Values can be of length `inMesh.cellCount()` interpolated to
            `outMesh.cellCenters()`, or `inMesh.nodeCount()` which are
            interpolated to `outMesh.positions()`.

    Returns
    -------
        Interpolated values.

    **2. Mesh based values to arbitrary points, based on finite elements**
      syntactic sugar for the core based `pg.core.interpolate`

    Parameters
    ----------
        args: :gimliapi:`GIMLI::Mesh`, ...
            Arguments forwarded to :gimliapi:`GIMLI::interpolate`
        kwargs:
            Arguments forwarded to :gimliapi:`GIMLI::interpolate`

            `interpolate(srcMesh, destMesh)`
            All data from inMesh are interpolated to outMesh

    Returns
    -------
        Interpolated values

    **3. Interpolate along curve.**
      Forwarded to :py:mod:`pygimli.meshtools.interpolateAlongCurve`

    Parameters
    ----------
        args: curve, t

        kwargs:
            Arguments forwarded to
            :py:mod:`pygimli.meshtools.interpolateAlongCurve`

        periodic: bool [False]
            Curve is periodic.
            Useful for closed parametric spline interpolation.

    **4. 1D interpolation**
      1D point set :math:`u(x)` for ascending :math:`x`.
      Find interpolation function :math:`I = u(x)` and
      returns :math:`u_{\text{i}} = I(x_{\text{i}})`
      (interpolation methods are [**linear** via numpy,
      cubic **spline** via scipy, fit **harmonic** functions' via pygimli])
      Note, for 'linear' and 'spline' the interpolate contains all original
      coordinates while 'harmonic' returns an approximate best fit.
      The amount of harmonic coefficients can be specified by the 'nc' keyword.

    Parameters
    ----------
        args: xi, x, u
            * :math:`x_{\text{i}}` - target sample points
            * :math:`x` - function sample points
            * :math:`u` - function values

        kwargs:
            * method : string ['linear']
                Specify interpolation method 'linear, 'spline', 'harmonic'
            * nc : int
                Number of harmonic coefficients for the 'harmonic' method.
            * periodic : bool [False]
                Curve is periodic.
                Useful for closed parametric spline interpolation.

    Returns
    -------
        ui: array of length xi
            :math:`u_{\text{i}} = I(x_{\text{i}})`, with :math:`I = u(x)`

    To use the core functions :gimliapi:`GIMLI::interpolate` start with a
    mesh instance as first argument or use the appropriate keyword arguments.

    TODO
    ----
        * 2D parametric to points (method=['linear, 'spline', 'harmonic'])
        * 2D/3D point cloud to points/grids
            ('Delauney', 'linear, 'spline', 'harmonic')
        * Mesh to points based on nearest neighbor values (pgcore)

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> fig, ax = pg.plt.subplots(1, 1, figsize=(10, 5))
    >>> u = np.array([1.0, 12.0, 3.0, -4.0, 5.0, 6.0, -1.0])
    >>> xu = np.array(range(len(u)))
    >>> xi = np.linspace(xu[0], xu[-1], 1000)
    >>> _= ax.plot(xu, u, 'o')
    >>> _= ax.plot(xi, pg.interpolate(xi, xu, u, method='linear'),
    ...         color='blue', label='linear')
    >>> _= ax.plot(xi, pg.interpolate(xi, xu, u, method='spline'),
    ...            color='red', label='spline')
    >>> _= ax.plot(xi, pg.interpolate(xi, xu, u, method='harmonic'),
    ...         color='green', label='harmonic')
    >>> _= ax.legend()
    >>>
    >>> pg.plt.show()
    """
    fallback = kwargs.pop('fallback', 0.0)
    verbose = kwargs.pop('verbose', False)

    if 'srcMesh' in kwargs and len(args) == 0:
        return pg.core.interpolate(**kwargs,
                                   fillValue=fallback,
                                   verbose=verbose)

    elif len(args) > 0:
        # pg._b([type(a) for a in args])

        if isinstance(args[0], pg.Mesh):
            if len(args) == 2 and isinstance(args[1], pg.Mesh):
                return pg.core.interpolate(args[0], args[1],
                                           fillValue=fallback,
                                           verbose=verbose)

        if len(args) == 4 and isinstance(args[3], str):
            return pg.interpolate(args[0], args[1], args[2], method=args[3],
                                  **kwargs)

            # if args[1].ndim == 2:  # outData = (inMesh, mat(dim>1), vR3)
            #     outMat = pg.Matrix()
            #     pg.core.interpolate(args[0], inMat=np.array(args[1]),
            #                         destPos=args[2], outMat=outMat,
            #                         fillValue=fallback,
            #                         verbose=verbose)
            #     return np.array(outMat)

        if len(args) == 4:  # args: (inMesh, inData(dim==1), outPos, outData)

            if args[1].ndim == 1 and args[2].ndim == 1 and args[3].ndim == 1:
                return pg.core.interpolate(args[0], inVec=args[1],
                                           x=args[2], y=args[3],
                                           fillValue=fallback,
                                           verbose=verbose)

            if isinstance(args[1], pg.Matrix) and \
               isinstance(args[3], pg.Matrix):
                return pg.core.interpolate(args[0], inMat=args[1],
                                           destPos=args[2],
                                           outMat=args[3],
                                           fillValue=fallback,
                                           verbose=verbose)
            if isinstance(args[1], pg.Vector) and \
               isinstance(args[3], pg.Vector):
                return pg.core.interpolate(args[0], inVec=args[1],
                                           destPos=args[2],
                                           outVec=args[3],
                                           fillValue=fallback,
                                           verbose=verbose)

        if len(args) == 5: # (mesh, Vec, x=Vec, y=Vec, z=Vec)
            if args[1].ndim == 1 and args[2].ndim == 1 and \
               args[3].ndim == 1 and args[4].ndim == 1:
                return pg.core.interpolate(args[0], inVec=args[1],
                                           x=args[2], y=args[3],
                                           z=args[4],
                                           fillValue=fallback,
                                           verbose=verbose)


        if len(args) == 3 and pg.isPosList(args[2]):
            # args: (inMesh, inData(dim==1), posList)
            if isinstance(args[1], list):
                inMat = np.asarray(args[1])
            else:
                inMat = args[1]

            if inMat.ndim == 2:
                outMat = pg.Matrix()
                pg.core.interpolate(args[0],
                                    inMat=inMat,
                                    destPos=args[2],
                                    outMat=outMat,
                                    fillValue=fallback,
                                    verbose=verbose)
                return outMat
            else:
                return pg.core.interpolate(args[0], inMat, destPos=args[2],
                                           fillValue=fallback,
                                           verbose=verbose)

        try:
            # check for Boost.Python.ArgumentError
            return pg.core.interpolate(*args, **kwargs,
                                        fillValue=fallback,
                                        verbose=verbose)
        except Exception as e:

            pass


    # end if pg.core:

    if len(args) == 3:

        if isinstance(args[0], pg.Mesh):  # args: (outMesh, inMesh, data)
            outMesh = args[0]
            inMesh = args[1]
            data = args[2]

            if isinstance(data, pg.PosVector):
                x = pg.interpolate(outMesh, inMesh, pg.x(data))
                y = pg.interpolate(outMesh, inMesh, pg.y(data))
                z = pg.interpolate(outMesh, inMesh, pg.z(data))
                return np.vstack([x, y, z]).T

            if isinstance(data, np.ndarray):
                if data.ndim == 2 and data.shape[1] == 3:
                    x = pg.interpolate(outMesh, inMesh, data[:, 0])
                    y = pg.interpolate(outMesh, inMesh, data[:, 1])
                    z = pg.interpolate(outMesh, inMesh, data[:, 2])
                    return np.vstack([x, y, z]).T

            if len(data) == inMesh.cellCount():
                return pg.interpolate(srcMesh=inMesh, inVec=data,
                                      destPos=outMesh.cellCenters())
            elif len(data) == inMesh.nodeCount():
                return pg.interpolate(srcMesh=inMesh, inVec=data,
                                      destPos=outMesh.positions())
            else:
                print(inMesh)
                print(outMesh)
                raise Exception("Don't know how to interpolate data of size",
                                str(len(data)))

            print("data: ", data)
            raise Exception("Cannot interpret data: ", str(len(data)))

        else:  # args: xi, x, u

            xi = args[0]
            x = args[1]
            u = args[2]

            method = kwargs.pop('method', 'linear')

            if 'linear' in method:
                return np.interp(xi, x, u)

            if 'harmonic' in method:
                coeff = kwargs.pop('nc', int(np.ceil(np.sqrt(len(x)))))
                from pygimli.frameworks import harmfitNative
                return harmfitNative(u, x=x, nc=coeff, xc=xi, err=None)[0]

            if 'spline' in method:
                if pg.optImport("scipy", requiredFor="Spline interpolation."):
                    from scipy import interpolate
                    tck = interpolate.splrep(x, u, s=0, k=min(3, (len(x)-1)),
                                             per=kwargs.pop('periodic', False))
                    return interpolate.splev(xi, tck, der=0)
                else:
                    return xi * 0.

    if len(args) == 2:  # args curve, t
        curve = args[0]
        t = args[1]
        return interpolateAlongCurve(curve, t, **kwargs)


def extract2dSlice(mesh, origin=None, normal=None, angle=0, dip=0, **kwargs):
    """Extract slice from 3D mesh as triangle mesh.

    Parameters
    ----------
    mesh : pg.Mesh
        Input mesh
    origin : [float, float, float]
        origin to be shifted [x, y, z]
    normal : [float, float, float] | str
        normal vector for extracting plane, or
        "x", "y", "z" equal to "yz", "xz", "yz", OR
    angle : float [0]
        azimuth of plane in the xy plane (0=x, 90=y)
    dip : float [0]
        angle to be tilted into the x'z plane (0=vertical)

    Optional keywords
    -----------------
    x,y,z : float [None]
        axis-parallel flices, determines normal and origin

    Returns
    -------
    2d triangular pygimli mesh with all data fields
    """
    from pygimli.viewer.pv import pgMesh2pvMesh
    # from pygimli.meshtools import convertPVPolyData  # to be written yet
    meshtmp = pg.Mesh(mesh)
    if origin is None and normal is None:  # use x/y/z
        if "z" in kwargs:
            origin = [0, 0, kwargs["z"]]
            normal = "z"
        elif "y" in kwargs:
            origin = [0, kwargs["y"], 0]
            normal = "y"
        elif "x" in kwargs:
            origin = [kwargs["x"], 0, 0]
            normal = "x"
        else:
            normal = "y"

    if origin == "center":
        bb = mesh.boundingBox()
        origin = [(bb.xMin()+bb.xMax())/2, (bb.yMin()+bb.yMax())/2, 0]

    if origin is not None:
        meshtmp.translate(-pg.Pos(origin))

    if isinstance(normal, str):  # "x", "yz" etc.
        if normal == "z" or normal == "xy":
            normal = [0, 0, 1]
        elif normal == "y" or normal == "xz":
            normal = [0, 1, 0]
        elif normal == "x" or normal == "yz":
            normal = [1, 0, 0]

    if angle:
        meshtmp.rotate(pg.Pos(0, 0, np.deg2rad(-angle)))
    if dip:
        meshtmp.rotate(pg.Pos(0, np.deg2rad(-dip), 0))

    pvmesh = pgMesh2pvMesh(meshtmp)

    pvs = pvmesh.slice(normal=normal, origin=[0, 0, 0],
                       generate_triangles=True)
    # return convertPVPolyData(pvs)  # that's the better way
    tri = pvs.faces.reshape((-1, 4))[:, 1:]
    mesh2d = pg.Mesh(dim=2)
    for point in pvs.points:
        mesh2d.createNode(point)

    for face in tri:
        mesh2d.createTriangle(*[mesh2d.node(p) for p in face])

    if normal == [1, 0, 0]:
        mesh2d.swapCoordinates(0, 1)

    for key in mesh.dataKeys():
        mesh2d[key] = pvs[key]

    mesh2d.setCellMarkers(pvs["Cell Marker"])
    mesh2d.swapCoordinates(1, 2)
    return mesh2d


if __name__ == '__main__':
    # no need to import matplotlib. pygimli's show does
    # import pygimli as pg
    import pygimli.meshtools as mt
    elec = np.arange(11.)
    topo = np.array([[0., 0.], [3., 2.], [4., 2.], [6., 1.], [10., 1.]])
    pg.plt.plot(topo[:, 0], topo[:, 1])
    p = mt.tapeMeasureToCoordinates(topo, elec)
    pg.plt.plot(p[:, 0], p[:, 1], 'o')
    pg.plt.gca().set_aspect(1)
    pg.wait()
