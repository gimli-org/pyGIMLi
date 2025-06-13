#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generic mesh visualization tools."""

import os
import sys
import time
import traceback

import numpy as np
from pygimli.viewer.mpl.colorbar import setMappableData

from .. core.logger import renameKwarg

import pygimli as pg
from .showmatrix import showMatrix
from .mpl import drawMesh, drawModel, drawField, drawSensors, drawStreams
from .mpl import drawSelectedMeshBoundaries
from .mpl import addCoverageAlpha
from .mpl import updateAxes
from .mpl import createColorBar, updateColorBar
from .mpl import CellBrowser
from .mpl.colorbar import cmapFromName


def show(obj=None, data=None, **kwargs):
    """ Mesh and model visualization.

    Syntactic sugar to show a obj with data. Forwards to
    a known visualization for obj. Typical is
    :py:mod:`pygimli.viewer.showMesh` or
    :py:mod:`pygimli.viewer.pv.showMesh3D` to show most of the typical 2D
    and 3D content.
    See tutorials and examples for usage hints. An empty show
    call creates an empty ax window.

    Parameters
    ----------
    obj: obj
        obj can be so far.
        * None (for empty axes)
        * int, int (for apropriate subplots)
        * :gimliapi:`GIMLI::Mesh` or list of meshes
        * DataContainer
        * pg.core.Sparse[Map]Matrix
        * [Pos,] -> scatter

    data: iterable
        Optionally data to visualize. See appropriate show function.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forward to appropriate show functions.

        * ax : axe [None]
            Matplotlib axes object. Create a new if necessary.
        * fitView : bool [True]
            Scale x and y limits to match the view.

    Returns
    -------
    Return the results from the showMesh* functions. Usually the axe object
    and a colorbar.

    See Also
    --------
    showMesh
    """
    def _removeFigHeader(ax):
        ## remove annoying 'Figure Nr.' for ipympl widgets
        if isinstance(ax, np.ndarray):
            return _removeFigHeader(ax[0])
        ax.figure.canvas.header_visible = False

    if "axes" in kwargs:  # remove me in 1.2 #20200515
        pg.critical("Deprecation Warning: Please use "
                    "keyword `ax` instead of `axes`")
        kwargs['ax'] = kwargs.pop('axes', None)

    ax = None
    # Empty call just to create an mpl axes
    # if obj is None and 'mesh' not in kwargs:
    if obj is None and 'mesh' not in kwargs.keys():
        ax = kwargs.pop('ax', None)
        if ax is None:
            ax = pg.plt.subplots(figsize=kwargs.pop('figsize', None))[1]

        _removeFigHeader(ax)
        return ax, None

    # create table of axes aka mpl.subplots
    if isinstance(obj, int):
        nrows = obj
        ncols = 1
        if isinstance(data, int):
            ncols = data

        fig, ax = pg.plt.subplots(nrows=nrows, ncols=ncols,
                          figsize=kwargs.pop('figsize', None))

        _removeFigHeader(ax)
        return ax, None

    ### obj contains a mesh
    if hasattr(obj, 'mesh'):
        ### data has values
        if hasattr(obj, 'values'):
            if obj.mesh.dim() > 1:
                return pg.show(obj.mesh, obj.values,
                           label=kwargs.pop('label', obj.name),
                            **kwargs)
            else:
                return show1D(obj.mesh, obj, **kwargs)
        ### data need to be evaluated
        if hasattr(obj, 'eval'):
            if obj.mesh.dim() > 1:
                return pg.show(obj.mesh, obj.eval(),
                               **kwargs)
            else:
                return show1D(obj.mesh, obj, **kwargs)

        ### try to evaluate obj later
        return pg.show(obj.mesh, obj, **kwargs)

    ### try to interpret obj as ERT Data
    if isinstance(obj, pg.DataContainerERT):
        from pygimli.physics.ert import showERTData
        return showERTData(obj, vals=kwargs.pop('vals', data), **kwargs)

    if isinstance(obj, pg.DataContainer):  # some other
        if obj.isSensorIndex("s") and obj.isSensorIndex("g"):
            from pygimli.viewer.mpl import showDataContainerAsMatrix
            if data is None and obj.haveData("t"):
                data = "t"

            return showDataContainerAsMatrix(obj, "s", "g", data)

    # try to interpret obj as matrices
    if isinstance(obj, pg.core.MatrixBase) or (isinstance(obj, np.ndarray) and
                                               obj.ndim == 2):
        return showMatrix(obj, **kwargs)

    try:
        from scipy.sparse import spmatrix
        if isinstance(obj, spmatrix):
            return showMatrix(obj, **kwargs)
    except ImportError:
        pass

    ### obj is list of Pos -> draw dots
    if pg.isPosList(obj):
        ax = kwargs.pop('ax', None)
        if ax is None:
            ax, cbar = pg.show()
        ax.scatter(pg.x(obj), pg.y(obj),  **kwargs)
        _removeFigHeader(ax)
        return ax, None


    # try to interprete obj as mesh or list of meshes
    mesh = kwargs.pop('mesh', obj)

    fitView = kwargs.get('fitView', True)

    if isinstance(mesh, list):
        ax = kwargs.pop('ax', None)
        label = kwargs.pop('label', None)
        if not isinstance(label, list):
            label = [label]* len(mesh)

        ax, cBar = show(mesh[0], data, hold=True, fitView=fitView,
                        ax=ax, label=label[0],
                        **kwargs)

        for i, m in enumerate(mesh[1:]):
            ax, cBar = show(m, data, ax=ax, hold=True, fitView=False,
                            label=label[i], **kwargs)

        if fitView is not False:
            ax.autoscale(enable=True, axis='both', tight=True)
            ax.set_aspect('equal')

        _removeFigHeader(ax)
        return ax, cBar

    if isinstance(mesh, pg.Mesh):
        if mesh.dim() == 1:
            return show1D(mesh, data, **kwargs)

        elif mesh.dim() == 2:
            if pg.zero(pg.y(mesh)):
                pg.info("swap z<->y coordinates for visualization.")
                meshSwap = pg.Mesh(mesh)
                for n in meshSwap.nodes():
                    n.pos()[1] = n.pos()[2]
                return showMesh(meshSwap, data, **kwargs)

            return showMesh(mesh, data, **kwargs)

        elif mesh.dim() == 3:

            from .pv import showMesh3D
            return showMesh3D(mesh, data, **kwargs)
        else:
            pg.error("ERROR: Mesh not valid.", mesh)

    if isinstance(obj, pg.core.Boundary):
        ax = kwargs.pop('ax', None)
        drawSelectedMeshBoundaries(ax, [obj], **kwargs)

        _removeFigHeader(ax)
        return ax, None

    pg.error("Can't interprete obj: {0} to show.".format(obj))
    return None, None


def showMesh(mesh, data=None, block=False, colorBar=None,
             label=None, coverage=None, ax=None, savefig=None,
             showMesh=False, showBoundary=None, factor=1,
             markers=False, **kwargs):
    """2D Mesh visualization.

    Create an axis object and plot a 2D mesh with given node or cell data.
    Returns the axis and the color bar. The type of data determines the
    appropriate draw method.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh
    data: iterable [None]
        Optionally data to visualize.

        . None (draw mesh only)
            forward to :py:mod:`pygimli.viewer.mpl.drawMesh`
            or if no cells are given:
            forward to :py:mod:`pygimli.viewer.mpl.drawPLC`

        . [[marker, value], ...]
            List of Cellvalues per cell marker
            forward to :py:mod:`pygimli.viewer.mpl.drawModel`

        . float per cell -- model, patch
            forward to :py:mod:`pygimli.viewer.mpl.drawModel`

        . float per node -- scalar field
            forward to :py:mod:`pygimli.viewer.mpl.drawField`

        . iterable of type [float, float] -- vector field
            forward to :py:mod:`pygimli.viewer.mpl.drawStreams`

        . pg.PosVector -- vector field
            forward to :py:mod:`pygimli.viewer.mpl.drawStreams`

        . pg.core.stdVectorRVector3 -- sensor positions
            DEPRECATED from remove of pg.core.stdVectorRVector3
            forward to :py:mod:`pygimli.viewer.mpl.drawSensors`
    block: bool [False]
        Force to open the Figure of your content and blocks the script until
        you close the current figure. Same like pg.show(); pg.wait()
    colorBar: bool [None], Colorbar
        Create and show a colorbar. If colorBar is a valid colorbar then only
        its values will be updated.
    label: str
        Set colorbar label. If set colorbar is toggled to True. [None]
    coverage: iterable [None]
        Weight data by the given coverage array and fadeout the color.
    ax: matplotlib.Axes [None]
        Instead of creating a new and empty ax, just draw into the given one.
        Useful to combine multiple plots into one figure.
    savefig: string
        Filename for a direct save to disc.
    showMesh: bool [False]
        Shows the mesh itself additional.
    showBoundary: bool [None]
        Highlight all boundaries with marker != 0. None means automatic.
        True for cell data and False for node data.
    marker: bool [False]
        Show cell markers and boundary marker.
    boundaryMarkers: bool [False]
        Highlight boundaries with marker !=0 and add Marker annotation.
        Applies :py:mod:`pygimli.viewer.mpl.drawBoundaryMarkers`.
        Dictionary "boundaryProps" can be added and will be forwarded to
        :py:mod:`pygimli.viewer.mpl.drawBoundaryMarkers`.

    Keyword Arguments
    -----------------
    xl: str ["$x$ in m"]
        Add label to the x axis. Default is '$x$ in m'

    yl: str [None]
        Add label to the y axis. Default is '$y$ in m' or 'Depth in m' with
        world boundary markers.

    fitView: bool
        Fit the axes limits to the all content of the axes. Default True.

    boundaryProps: dict
        Arguments for plot boundary.

    hold: bool [pg.hold()]
        Holds back the opening of the Figure.
        If set to True [default] nothing happens until you either force another
        show with hold=False or block=True or call pg.wait() or pg.plt.show().
        If hold is set to False your script will open the figure and continue
        working. You can change global hold with pg.hold(bool).

    axisLabels: bool [True]
        Set x/yLabels for ax. X will be "$x$ in m" and "$y$ in m".
        Y ticks change to depth values for a mesh with world
        boundary markers and the label becomes "Depth in m".

    All remaining will be forwarded to the draw functions
    and matplotlib methods, respectively.

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> world = mt.createWorld(start=[-10, 0], end=[10, -10],
    ...                        layers=[-3, -7], worldMarker=True)
    >>> mesh = mt.createMesh(world, quality=32, area=0.2, smooth=[1, 10])
    >>> _ = pg.viewer.showMesh(mesh, markers=True, xl='$x$-coordinate')

    Returns
    -------
    ax : matplotlib.axes

    cBar : matplotlib.colorbar
    """
    renameKwarg('cmap', 'cMap', kwargs)
    cMap = kwargs.pop('cMap', 'viridis')
    cBarOrientation = kwargs.pop('orientation', 'horizontal')
    replaceData = kwargs.pop('replaceData', False)
    axisLabels = kwargs.pop('axisLabels', True)
    xl = kwargs.pop('xl', "$x$ in m")
    yl = kwargs.pop('yl', None)

    if ax is None:
        ax, _ = pg.show(figsize=kwargs.pop('figsize', None), **kwargs)

    # adjust limits only when axis is empty
    fitViewDefault = True
    # if (ax.lines or ax.collections or ax.patches):
    #     fitViewDefault = False
    # else:

    # plt.subplots() resets locale setting to system default .. this went
    # horrible wrong for german 'decimal_point': ','
    pg.checkAndFixLocaleDecimal_point(verbose=False)

    from pygimli.viewer.mpl.utils import __holdAxes__
    hold = kwargs.pop('hold', __holdAxes__)
    # hold = kwargs.pop('hold', pg.viewer.mpl.utils.__holdAxes__)

    if block is True:
        hold = True

    # lastHoldStatus = pg.viewer.mpl.utils.__holdAxes__
    lastHoldStatus = __holdAxes__
    pg.viewer.mpl.hold(val=hold)

    gci = None
    validData = False

    uniquemarkers = None
    ### if levels equals unqiue cell data, plot as markers
    if kwargs.get('levels', None) is not None:
        ud = pg.unique(data)
        levs = kwargs['levels']
        if len(ud) == len(levs) and ud == levs:
            kwargs["boundaryMarkers"] = kwargs.get("boundaryMarkers", False)

            uniquemarkers, uniqueidx = np.unique(np.array(data),
                                                 return_inverse=True)

            markers = True

    if markers:
        kwargs["boundaryMarkers"] = kwargs.get("boundaryMarkers", True)

        if mesh.cellCount() > 0:
            if uniquemarkers is None:
                uniquemarkers, uniqueidx = np.unique(
                                                np.array(mesh.cellMarkers()),
                                                return_inverse=True)
                label = "Cell markers"

            if cMap == 'viridis':
                cMap = "Set3"
            cMap = cmapFromName("Set3", ncols=len(uniquemarkers))

            kwargs["logScale"] = False
            kwargs["cMin"] = -0.5
            kwargs["cMax"] = len(uniquemarkers) - 0.5

            data = np.arange(len(uniquemarkers))[uniqueidx]

    if isinstance(data, str):
        #if data in mesh.dataKeys():
        #     data = mesh[data]
        #     # elif 0:  # maybe check x, y, z, cellMarker etc.
        # else:
        #     pg.error(f"Could not retrieve data from key {data}")
        #     return None, None
        if mesh.haveData(data):
            #print(factor)
            data = mesh[data] * factor
        else:
            raise IndexError("Mesh does not contain field ", data)
    elif callable(data):
        data = data(mesh.positions())

    if data is None:
        showMesh = True
        mesh.createNeighborInfos()

        if showBoundary is None:
            showBoundary = True

    elif pg.isPosList(data):
        drawStreams(ax, mesh, data, label=label, **kwargs)
    else:
        # check for map like data=[[marker, val], ....]
        if isinstance(data, list) and \
                isinstance(data[0], list) and isinstance(data[0][0], int):
            data = pg.solver.parseMapToCellArray(data, mesh)

        #pg._g(data)
        # if hasattr(data[0], '__len__'):
        #     print(hasattr(data[0], '__len__'),
        #       isinstance(data, np.ma.core.MaskedArray),
        #       isinstance(data[0], str),
        #       len(data), len(data[0]))

        if hasattr(data[0], '__len__') and not \
                isinstance(data, np.ma.core.MaskedArray) and not \
                isinstance(data[0], str) and not \
                (len(data) == 1 or len(data[0] == 0)):

            #pg._y(data)
            data = np.asarray(data)

            if len(data) == 2:  # [u,v] x N
                data = np.asarray(data).T

            if len(data) == 1: # 1xM matrix
                return showMesh(np.ravel(data), **kwargs)
            if data.shape[1] == 1:  # Mx1 matrix
                return showMesh(np.ravel(data), **kwargs)
            if data.shape[1] == 2:  # N x [u,v]
                drawStreams(ax, mesh, data, label=label, **kwargs)
            elif data.shape[1] == 3:  # N x [u,v,w]
                # if sum(data[:, 0]) != sum(data[:, 1]):
                # drawStreams(ax, mesh, data, **kwargs)
                drawStreams(ax, mesh, data[:, :2], label=label, **kwargs)
            else:

                ### Try animation frames x N
                if data.ndim == 2:
                    if data.shape[1] == mesh.cellCount() or \
                       data.shape[1] == mesh.nodeCount():

                        return showAnimation(mesh, data, cMap=cMap,
                                             ax=ax, **kwargs)

                pg.warn("No valid stream data:", data.shape, data.ndim)
                showMesh = True

        else:
            data = np.asarray(data)
            if bool(colorBar) is not False:
                colorBar = True

            if kwargs.pop("contour", False):
                data = pg.meshtools.cellDataToNodeData(mesh, data)
                kwargs.setdefault("nLevs", 11)

            if len(data) == mesh.cellCount():
                if showBoundary is None:
                    showBoundary = True

            def _drawField(ax, mesh, data, kwargs):  # like view.mpl.drawField?
                # kwargs as reference here to set defaults valid outside too

                validData = True
                if len(data) == mesh.cellCount():
                    kwargs['nCols'] = kwargs.pop('nCols', 256)
                    gci = drawModel(ax, mesh, data, **kwargs)

                elif len(data) == mesh.nodeCount():
                    kwargs['nLevs'] = kwargs.pop('nLevs', 5)
                    kwargs['nCols'] = kwargs.pop('nCols', kwargs['nLevs']-1)

                    gci = drawField(ax, mesh, data, symLevels='RdBu' in cMap,
                                    **kwargs)
                else:
                    pg.error("Data size invalid")
                    print("Data: ", data.shape)
                    print("Data: ", len(data), min(data), max(data),
                          pg.core.haveInfNaN(data))
                    print("Mesh: ", mesh)
                    validData = False
                    drawMesh(ax, mesh)

                return gci, validData


            try:
                if label is None:
                    label = ""

                if replaceData and hasattr(mesh, 'gci') and ax in mesh.gci:
                    gci = mesh.gci[ax]

                    if 'TriContourSet' in str(type(gci)):
                        ax.clear()
                        gci, validData = _drawField(ax, mesh, data, kwargs)
                        updateAxes(ax, force=True)
                    else:
                        setMappableData(gci, data,
                                        cMin=kwargs.get('cMin', None),
                                        cMax=kwargs.get('cMax', None),)
                        updateAxes(ax, force=True)
                        return ax, gci.colorbar
                else:
                    gci, validData = _drawField(ax, mesh, data, kwargs)

                # Cache mesh and scalarmappable to make replaceData work
                if not hasattr(mesh, 'gci'):
                    mesh.gci = {}

                mesh.gci[ax] = gci

                if cMap is not None and gci is not None:
                    gci.set_cmap(cmapFromName(cMap))
                    # gci.cmap.set_under('k')

            except BaseException as ex:  # super ugly!
                print(ex)
                traceback.print_exc(file=sys.stdout)

    if mesh.cellCount() == 0:
        showMesh = False
        if mesh.boundaryCount() == 0:
            pg.viewer.mpl.drawPLC(ax, mesh, showNodes=True,
                                  fillRegion=False,
                                  showBoundary=False,
                                  **kwargs)
            showBoundary = False
            # ax.plot(pg.x(mesh), pg.y(mesh), '.', color='black')
        else:
            kwargs['orientation'] = cBarOrientation
            pg.viewer.mpl.drawPLC(ax, mesh, **kwargs)

    if showMesh:
        if gci is not None and hasattr(gci, 'set_antialiased'):
            gci.set_antialiased(True)
            gci.set_linewidth(0.3)
            gci.set_edgecolor(kwargs.pop('color', "0.1"))
            #drawMesh(ax, mesh, lw=0.3, **kwargs)
        #else:

        if mesh.dim() == 1:
            return show1D(mesh, None, showMesh=True, ax=ax)

        drawMesh(ax, mesh, lw=0.3, **kwargs)
        # pg.viewer.mpl.drawSelectedMeshBoundaries(ax,
        #         mesh.boundaries(),
        #         color=kwargs.pop('color', "0.1"),
        #         linewidth=kwargs.pop('lw', 0.3))

    if bool(showBoundary) is True:
        b = mesh.boundaries(mesh.boundaryMarkers() != 0)
        pg.viewer.mpl.drawSelectedMeshBoundaries(ax, b,
                                                 color=(0.0, 0.0, 0.0, 1.0),
                                                 linewidth=1.4)

    if kwargs.pop("boundaryMarkers", False):
        pg.viewer.mpl.drawBoundaryMarkers(ax, mesh,
                    clipBoundaryMarkers=kwargs.pop('clipBoundaryMarkers', False),
                    bc=kwargs.pop('bc', None),
                    **kwargs.pop('boundaryProps', {}) )

    fitView = kwargs.pop('fitView', fitViewDefault)

    if fitView is not False:
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_aspect(kwargs.pop('aspect', 'equal'))

    cBar = None

    if label is not None and colorBar is None:
        colorBar = True

    # pg._r(validData, gci)
    if validData:
        labels = ['cMin', 'cMax', 'nCols', 'nLevs', 'logScale', 'levels']
        subkwargs = {key: kwargs[key] for key in labels if key in kwargs}

        subkwargs['cMap'] = cMap

        if isinstance(colorBar, bool):

            if colorBar is True:
                subkwargs['label'] = label
                subkwargs['orientation'] = cBarOrientation

            cBar = createColorBar(gci,
                                  size=kwargs.pop('size', 0.2),
                                  pad=kwargs.pop('pad', None),
                                  **subkwargs,
                                  onlyColorSet=not colorBar)

        elif colorBar is not False:
            cBar = updateColorBar(colorBar, **subkwargs)

        if markers and cBar is not None:
            ticks = np.arange(len(uniquemarkers))
            cBar.set_ticks(ticks)
            labels = []
            for marker in uniquemarkers:
                labels.append(pg.pf(marker, mathtex=True))
            cBar.set_ticklabels(labels)

    if coverage is not None:
        if isinstance(coverage, (float, int)):
            gci.set_alpha(coverage)
        elif len(data) == len(coverage) == mesh.cellCount():
            addCoverageAlpha(gci, coverage,
                             dropThreshold=kwargs.pop('dropThreshold', 0.4))
        else:
            pg.error('Coverage needs to be either of type float or an array',
                     'with the same length as data and mesh.cellCount().')
            # addCoverageAlpha(gci, pg.core.cellDataToPointData(mesh,
            #                                                   coverage))

    if not hold or block is not False and \
            pg.plt.get_backend().lower() != "agg":
        if data is not None:
            if len(data) == mesh.cellCount():
                CellBrowser(mesh, data, ax=ax)

        pg.plt.show(block=block)
        try:
            pg.plt.pause(0.01)
        except BaseException:
            pass

    if axisLabels == True and mesh.dim() == 2:

        try:
            useDepth = min(mesh.boundaryMarkers()) < 0 and max(pg.y(mesh)) <= 0
            pg.viewer.mpl.adjustWorldAxes(ax, useDepth=useDepth, xl=xl, yl=yl)
        except BaseException:
            pass
    else:
        pg.viewer.mpl.updateAxes(ax)

    pg.viewer.mpl.hold(val=lastHoldStatus)

    if savefig:
        print('saving: ' + savefig + ' ...', end="")
        if '.' not in savefig:
            savefig += '.pdf'

        ax.figure.savefig(savefig, bbox_inches='tight')
        # rc params savefig.format=pdf
        print('.. done')

    return ax, cBar


def showBoundaryNorm(mesh, normMap=None, **kwargs):
    """Show mesh boundaries normals.

    Show the mesh and draw a black line along the normal direction of all
    boundaries. If you provide a boundary marker vs. norm direction map,
    then only these norms are drawn.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    normMap : list
        list of [boundary marker, [norm]] pairs. e.g. [[1, [0.0,1.0]], ... ]

    **kwargs :
        Will be forwarded to the draw functions and matplotlib methods,
        respectively.

    Returns
    -------
    ax : matplotlib.ax
    """
    ax = kwargs.pop('ax', None)

    col = kwargs.pop('color', 'Black')

    if normMap:
        for pair in normMap:
            bounds = mesh.findBoundaryByMarker(pair[0])

            for b in bounds:
                c1 = b.center()

                if (pair[1][0] != 0) or (pair[1][1] != 0):
                    ax.arrow(c1[0], c1[1], pair[1][0], pair[1][1],
                             head_width=0.1, head_length=0.3, color=col,
                             **kwargs)
                else:
                    ax.plot(c1[0], c1[1], 'o', color=col)
        return

    ax = show(mesh, hold=True, ax=ax)[0]
    for b in mesh.boundaries():
        c1 = b.center()
        c2 = c1 + b.norm()
        ax.plot([c1[0], c2[0]], [c1[1], c2[1]], color=col, **kwargs)

    time.sleep(0.05)

    return ax


def show1D(mesh, obj, **kwargs):
    """Show simple plot for 1D modelling results
    """
    kwargs.pop('hold', None)
    kwargs.pop('fitView', None)

    ax = kwargs.pop('ax', None)
    newAxe = False

    if ax is None:
        newAxe = True
        ax = pg.show()[0]

    if hasattr(obj, 'eval'):
        x = pg.sort(pg.x(mesh))
        v = obj(x)

        if hasattr(v, 'ndim') and v.ndim == 2 and v.shape[0] == mesh.nodeCount():
        # Vector Field -- show x--component
            v = v[:,0]

    elif hasattr(obj, 'values'):
        pg._r(kwargs)
        pg._r(mesh)
        pg._r(obj)
        pg.critical('implementme')

    elif pg.isArray(obj, mesh.nodeCount()):
        x = pg.sort(pg.x(mesh))
        v = obj

    elif isinstance(obj, list):
        return show1D(mesh, np.array(obj), ax=ax, **kwargs)

    elif hasattr(obj, 'ndim') and obj.ndim == 2 and pg.isArray(obj[0], mesh.nodeCount()):
        # list of values for animation
        return showAnimation(mesh, obj, ax=ax, **kwargs)

    elif obj is None:
        showMesh = kwargs.pop('showMesh', False)

        if showMesh:
            pg.viewer.mpl.drawSelectedMeshBoundaries(ax, mesh.cells(),
                                                     color='k', linewidth=0.3,
                                                     linestyle="-", **kwargs)


        pg.viewer.mpl.drawPLC(ax, mesh, showNodes=True,
                              fillRegion=False,
                              showBoundary=False,
                                  **kwargs)
        return ax, None
    else:
        pg._r(kwargs)
        pg._r(mesh)
        pg._r(obj)
        #pg._r(obj.shape)
        pg.critical('implementme')

    swapAxes = kwargs.pop('swapAxes', False)
    label = kwargs.pop('label', None)

    if label is None and hasattr(obj, '_OP'):
        label = obj.__str__()

    grid = kwargs.pop('grid', True)
    xLabel = kwargs.pop('xl', 'Depth in m')

    if swapAxes is True:
        ax.set_ylabel(xLabel)
        pg.viewer.mpl.renameDepthTicks(ax)
    else:
        ax.set_xlabel(xLabel)

    if 'yl' in kwargs:
        if not isinstance(obj, (list, np.ndarray)):
            yLabel = kwargs.pop('yl', str(obj))
        else:
            yLabel = kwargs.pop('yl', None)

        if swapAxes is True:
            ax.set_xlabel(yLabel)
        else:
            ax.set_ylabel(yLabel)

    if swapAxes is True:
        curve = ax.plot(v, x, label=label, **kwargs)
    else:
        curve = ax.plot(x, v, label=label, **kwargs)

    if label is not None and label != '':
        ax.legend()
    ax.grid(grid)

    return ax, curve


__Animation_Keeper__ = None

def showAnimation(mesh, data, ax=None, **kwargs):
    """Show timelapse mesh data.

    Time will be annotated if the mesh contains a valid 'times' data array.
    Note, there can be only one animation per time.

    Best viewed in a notebook, because of caching and better animation control
    elements.

    TODO
    ----
        * 3D
        * allow for multiple animations per script

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        2D GIMLi mesh
    data: [NxM] iterable
        Data matrix for N frames with suitable data size.

    Keyword Args
    ------------
    plc: Mesh
        Overlay plc mesh.
    dpi: int[96]
        Movie resolution.
    interval: int [20]
        Forwarded to Matplotlib animation plotter.
    ffmpeg: str [None]
        Write ffmpeg movie with name *.mp4. Need ffmpeg package installed.
    **kwargs:
        Forwarded to :py:func:pygimli.viewer.show:

    """
    import matplotlib.animation
    plt = pg.plt

    plt.rcParams["animation.html"] = "jshtml"
    plt.rcParams['figure.dpi'] = kwargs.pop('dpi', 96)
    plt.rcParams['animation.embed_limit'] = 50

    figsize = kwargs.pop("figsize", None)
    flux = kwargs.pop('flux', None)
    plc = kwargs.pop('plc', None)
    interval = kwargs.pop('interval', 20)
    swapAxes = kwargs.get('swapAxes', False)
    ffmpeg = kwargs.pop('ffmpeg', None)

    plt.ioff()

    if mesh.dim() == 1:
        ax, curve = pg.show(mesh, data[0], ax=ax, **kwargs)
        # ensure p2 meshs are sorted ascending in x
        sortX = np.argsort(pg.x(mesh))

        if swapAxes is True:
            ax.set_xlim(min(data.flatten()), max(data.flatten()))
        else:
            ax.set_ylim(min(data.flatten()), max(data.flatten()))
    else:
        ax = pg.show(mesh, data[0], ax=ax, **kwargs)[0]
        if flux is not None:
            pg.show(mesh, flux[0], ax=ax)

        ax.figure.tight_layout()

    try:
        times = mesh['times']
    except Exception as e:
        times = None

    p = pg.utils.ProgressBar(len(data))

    def animate(t):
        p.update(t)
        if mesh.dim() == 1:
            if swapAxes is True:
                curve[0].set_xdata(data[t][sortX])
            else:
                curve[0].set_ydata(data[t][sortX])
        else:
            ax.clear()
            pg.show(mesh, data[t], ax=ax, **kwargs)

            if flux is not None:
                try:
                    pg.show(mesh, flux[t], ax=ax)
                except:
                    pass

        if plc is not None:
            pg.viewer.mpl.drawMesh(ax, plc, fillRegion=False, fitView=False)

        if times is not None and len(times) > t:
            # ax.text(0.02, 0.02, f't={pg.pf(times[t])}',
            ax.text(0.01, 1.01, f't={pg.utils.prettyTime(times[t])}',
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    transform=ax.transAxes, color='k', fontsize=8)

    if pg.isNotebook() is False:
        global __Animation_Keeper__

    __Animation_Keeper__ = matplotlib.animation.FuncAnimation(ax.figure,
                                                              animate,
                                                              interval=interval,
                                                              frames=len(data))

    if ffmpeg is not None:
        __Animation_Keeper__.save(ffmpeg +'.mp4',
                              writer=matplotlib.animation.FFMpegWriter(
                                    fps=14, bitrate=5000, codec='h264'),
                              dpi=plt.rcParams['figure.dpi']
                              )

    return __Animation_Keeper__

