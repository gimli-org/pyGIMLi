# -*- coding: utf-8 -*-
"""Draw mesh/model/fields with matplotlib."""

import textwrap
import numpy as np
import pygimli as pg

from pygimli.utils import streamline

from .colorbar import autolevel, createColorBar, updateColorBar, cmapFromName
from .utils import updateAxes as updateAxes_


class CellBrowserCacheSingleton(object):
    __instance = None
    cbCache_ = []

    def __new__(cls):
        if CellBrowserCacheSingleton.__instance is None:
            CellBrowserCacheSingleton.__instance = object.__new__(cls)
        return CellBrowserCacheSingleton.__instance

    def add(self, c):
        self.cbCache_.append(c)

    def remove(self, c):
        self.cbCache_.remove(c)


# We only want one instance of this global cache so its a singleton class
__CBCache__ = CellBrowserCacheSingleton()


class CellBrowser(object):
    """Interactive cell browser on current or specified ax for a given mesh.

    Cell information can be displayed by mouse picking. Arrow keys up and down
    can be used to scroll through the cells, while ESC closes the cell
    information window.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        The plotted 2D mesh to browse through.
    data : iterable
        Cell data.
    ax : mpl axis instance, optional
        Axis instance where the mesh is plotted (default is current axis).

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawModel
    >>> from pygimli.viewer.mpl import CellBrowser
    >>>
    >>> mesh = pg.createGrid(range(5), range(5))
    >>> fig, ax = plt.subplots()
    >>> plc = drawModel(ax, mesh, mesh.cellMarkers())
    >>> browser = CellBrowser(mesh)
    >>> browser.connect()
    """

    def __init__(self, mesh, data=None, ax=None):
        """ Construct CellBrowser on a specific `mesh`.
        """
        if ax:
            self.ax = ax
        else:
            self.ax = pg.plt.gca()

        self._connected = False

        self.fig = self.ax.figure
        self.mesh = None
        self.data = None
        self.highLight = None
        self.text = None

        self.cellID = None
        self.event = None
        self.artist = None
        self.pid = None
        self.kid = None
        self.text = None

        self.setMesh(mesh)
        self.setData(data)
        self.connect()

    def __del__(self):
        """ Deregister if the cellBrowser has been deleted.
        """
        self.disconnect()

    def connect(self):
        """ Connect to matplotlib figure canvas.
        """
        if not self._connected:
            self.pid = self.fig.canvas.mpl_connect('pick_event', self.onPick)
            self.kid = self.fig.canvas.mpl_connect('key_press_event',
                                                   self.onPress)
            __CBCache__.add(self)
            self._connected = True

    def disconnect(self):
        """ Disconnect from matplotlib figure canvas.
        """
        if self._connected:
            __CBCache__.remove(self)
            self.fig.canvas.mpl_disconnect(self.pid)
            self.fig.canvas.mpl_disconnect(self.kid)
            self._connected = False

    def initText(self):
        """ Initialize hint text properties.
        """
        import matplotlib as mpl
        bbox = dict(boxstyle='round, pad=0.5', fc='w', alpha=0.5)
        arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0.5')
        kwargs = dict(fontproperties='monospace', visible=False,
                      fontsize=mpl.rcParams['font.size'] - 2, weight='bold',
                      xytext=(50, 20), arrowprops=arrowprops,
                      textcoords='offset points', bbox=bbox, va='center')

        self.text = self.ax.annotate(None, xy=(0, 0), **kwargs)

    def setMesh(self, mesh):
        """ Set mesh.
        """
        self.mesh = mesh

    def setData(self, data=None):
        """ Set data, if not set look for the artist array data.
        """
        self.hide()
        if data is not None:
            if len(data) == self.mesh.cellCount():
                self.data = data
            elif len(data) == self.mesh.nodeCount():
                self.data = pg.meshtools.nodeDataToCellData(self.mesh, data)
            else:
                pg.warn('Data length mismatch mesh.cellCount(): ' +
                        str(len(data)) + "!=" + str(self.mesh.cellCount()) +
                        ". Mapping data to cellMarkers().")
                self.data = data[self.mesh.cellMarkers()]

    def hide(self):
        """ Hide info window.
        """
        self.cellID = -1

        if self.text is not None:
            self.text.set_visible(False)

        self.removeHighlightCell()

        self.fig.canvas.draw()

    def removeHighlightCell(self):
        """ Remove cell highlights.
        """
        if self.highLight is not None:
            if self.highLight in self.ax.collections:
                self.highLight.remove()
            self.highLight = None

    def highlightCell(self, cell):
        """ Highlight selected cell.
        """
        import matplotlib as mpl
        self.removeHighlightCell()
        self.highLight = mpl.collections.PolyCollection(
            [_createCellPolygon(cell)])
        self.highLight.set_edgecolor('0')
        self.highLight.set_linewidth(1.5)
        self.highLight.set_facecolor([0.9, 0.9, 0.9, 0.4])
        self.ax.add_collection(self.highLight)

    def onPick(self, event):
        """ Call `self.update()` on mouse pick event.
        """
        self.event = event
        self.artist = event.artist

        if self.data is None:
            self.data = self.artist.get_array()
            # self.edgeColors = self.artist.get_edgecolors()

        if 'mouseevent' in event.__dict__.keys():
            # print(event.__dict__.keys())
            # print(event.mouseevent)
            if (event.mouseevent.xdata is not None and
                    event.mouseevent.ydata is not None and
                    event.mouseevent.button == 1):
                c = self.mesh.findCell((event.mouseevent.xdata,
                                        event.mouseevent.ydata))
                if c and self.cellID != c.id():
                    self.cellID = c.id()
                else:
                    self.cellID = -1

                self.update()
        else:  # variant before (seemed inaccurate)
            self.cellID = event.ind[0]

    def onPress(self, event):
        """ Call `self.update()` if up, down, or escape keys are pressed.
        """
        # print(event, event.key)
        if self.data is None:
            return
        if event.key not in ('up', 'down', 'escape'):
            return
        if event.key == 'up':
            if self.cellID is not None:
                self.cellID += 1
        elif event.key == 'down':
            if self.cellID is not None:
                self.cellID -= 1
        else:
            self.hide()
            return

        if self.cellID is not None:
            self.cellID = int(np.clip(self.cellID, 0,
                                      self.mesh.cellCount() - 1))
            self.update()

    def update(self):
        """ Update the information window.

        Hide the information window for self.cellID == -1
        """
        try:
            if self.cellID > -1:
                cell = self.mesh.cell(self.cellID)
                center = cell.center()
                x, y = center.x(), center.y()
                marker = cell.marker()
                data = self.data[self.cellID]
                header = "Cell %d:\n" % self.cellID
                header += "-" * (len(header) - 1)
                istr = "\nx: {:.2f}\n y: {:.2f}\n data: {:.2e}\n marker: {:d}"
                info = istr.format(x, y, data, marker)
                text = header + textwrap.dedent(info)

                if self.text is None or self.text not in self.ax.texts:
                    self.initText()

                self.text.set_text(text)
                self.text.xy = x, y
                self.text.set_visible(True)
                self.highlightCell(cell)
                self.fig.canvas.draw()
            else:
                self.hide()

        except BaseException as e:
            print(e)


def drawMesh(ax, mesh, fitView=True, **kwargs):
    """ Draw a 2d mesh into a given ax.

    Set the limits of the ax tor the mesh extent.

    Parameters
    ----------
    ax : mpl axe instance
        Axis instance where the mesh is plotted.
    mesh : :gimliapi:`GIMLI::Mesh`
        The 2D mesh which will be drawn.
    fitView : bool [True]
        Adjust ax limits to mesh bounding box.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forward to drawPLC or drawMeshBoundaries.

    %(drawMeshBoundaries)

    %(drawPLC)

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawMesh
    >>> n = np.linspace(1, 2, 10)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, mesh)
    >>> plt.show()
    """
    if mesh.cellCount() == 0 or mesh.dim() == 1:
        pg.viewer.mpl.drawPLC(ax, mesh, **kwargs)
    else:
        pg.viewer.mpl.drawMeshBoundaries(ax, mesh, **kwargs)

    if fitView is True:
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_aspect('equal')

    updateAxes_(ax)


def drawModel(ax, mesh, data=None, tri=False, rasterized=False,
              cMin=None, cMax=None, logScale=False,
              xlabel=None, ylabel=None, fitView=True, verbose=False,
              **kwargs):
    """ Draw a 2d mesh and color the cell by the data.

    Parameters
    ----------
    ax : mpl axis instance, optional
        Axis instance where the mesh is plotted (default is current axis).
    mesh : :gimliapi:`GIMLI::Mesh`
        The plotted mesh to browse through.
    data : array, optional
        Data to draw. Should either equal numbers of cells or nodes of the
        corresponding `mesh`.
    tri : boolean, optional
        use MPL tripcolor (experimental)
    rasterized : boolean, optional
        Rasterize mesh patches to reduce file size and avoid zooming artifacts
        in some PDF viewers.
    fitView : bool [True]
        Adjust ax limits to mesh bounding box.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forwarded to the draw functions and mpl methods,
        respectively.

    Returns
    -------
    gci : matplotlib graphics object

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawModel
    >>> n = np.linspace(0, -2, 11)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> mx = pg.x(mesh.cellCenter())
    >>> my = pg.y(mesh.cellCenter())
    >>> data = np.cos(1.5 * mx) * np.sin(1.5 * my)
    >>> fig, ax = plt.subplots()
    >>> drawModel(ax, mesh, data)
    <matplotlib.collections.PolyCollection object at ...>
    """
    # deprecated .. remove me
    if 'cMap' in kwargs or 'cmap' in kwargs:
        pg.warn('cMap|cmap argument is deprecated for draw functions. ' +
                'Please use show or customize a colorbar.')
    # deprecated .. remove me

    if mesh.nodeCount() == 0:
        pg.error("drawModel: The mesh is empty.", mesh)

    if tri or 'shading' in kwargs:
        gci = drawField(ax, mesh, data,
                        cMin=cMin, cMax=cMax, logScale=logScale,
                        **kwargs)
    else:
        gci = pg.viewer.mpl.createMeshPatches(ax, mesh, rasterized=rasterized,
                                              verbose=verbose)
        ax.add_collection(gci)

        if data is None:
            data = pg.Vector(mesh.cellCount())

        if len(data) != mesh.cellCount():
            pg.info("drawModel have wrong data length .. " +
                    " indexing data from cellMarkers()")
            viewdata = data[mesh.cellMarkers()]
        else:
            viewdata = data

        pg.viewer.mpl.setMappableData(gci, viewdata,
                                      cMin=cMin, cMax=cMax, logScale=logScale,
                                      **kwargs)

    gci.set_antialiased(True)
    gci.set_linewidths(0.1)
    gci.set_edgecolors("face")

    if xlabel is not None:
        ax.set_xlabel(xlabel)

    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if fitView is True:
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_aspect('equal')

    updateAxes_(ax)
    return gci


def drawSelectedMeshBoundaries(ax, boundaries, color=None, linewidth=1.0,
                               linestyle="-", **kwargs):
    """Draw mesh boundaries into a given axes.

    Parameters
    ----------
    ax : matplotlib axes
        axes to plot into
    boundaries : :gimliapi:`GIMLI::Mesh` boundary vector
        collection of boundaries to plot
    color : matplotlib color |str [None]
        matching color or string, else colors are according to markers
    linewidth : float [1.0]
        line width
    linestyles : linestyle for line collection, i.e. solid or dashed

    Returns
    -------
    lco : matplotlib line collection object
    """
    import matplotlib as mpl
    drawAA = True
    lines = []

    if hasattr(boundaries, '__len__'):
        if len(boundaries) == 0:
            return

    for bound in boundaries:
        lines.append(list(zip([bound.node(0).x(), bound.node(1).x()],
                              [bound.node(0).y(), bound.node(1).y()])))

    lineCollection = mpl.collections.LineCollection(lines,
                                                    antialiaseds=drawAA,
                                                    **kwargs)

    if color is None:
        viewdata = [b.marker() for b in boundaries]
        pg.viewer.mpl.setMappableValues(lineCollection, viewdata,
                                        logScale=False)
    else:
        lineCollection.set_color(color)

    lineCollection.set_linewidth(linewidth)
    lineCollection.set_linestyle(linestyle)
    ax.add_collection(lineCollection)

    updateAxes_(ax)

    return lineCollection


def drawSelectedMeshBoundariesShadow(ax, boundaries, first='x', second='y',
                                     color=(0.3, 0.3, 0.3, 1.0)):
    """Draw mesh boundaries as shadows into a given axes.

    Parameters
    ----------
    ax : matplotlib axes
        axes to plot into
    boundaries : :gimliapi:`GIMLI::Mesh` boundary vector
        collection of boundaries to plot
    first / second : str ['x' / 'y']
        attribute names to retrieve from nodes
    color : matplotlib color |str [None]
        matching color or string, else colors are according to markers
    linewidth : float [1.0]
        line width

    Returns
    -------
    lco : matplotlib line collection object
    """
    import matplotlib as mpl
    polys = []

    for cell in boundaries:
        polys.append(list(zip([getattr(cell.node(0), first)(),
                               getattr(cell.node(1), first)(),
                               getattr(cell.node(2), first)()],
                              [getattr(cell.node(0), second)(),
                               getattr(cell.node(1), second)(),
                               getattr(cell.node(2), second)()])))

    collection = mpl.collections.PolyCollection(polys, antialiaseds=True)

    collection.set_color(color)
    collection.set_edgecolor(color)
    collection.set_linewidth(0.2)
    ax.add_collection(collection)

    updateAxes_(ax)
    return collection


def drawBoundaryMarkers(ax, mesh, clipBoundaryMarkers=False, **kwargs):
    """Draw boundary markers for mesh.boundaries with marker != 0

    Args
    ----
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh that have the boundary markers.

    clipBoundaryMarkers: bool [False]
        Clip boundary marker to the axes limits if needed.

    Keyword Arguments
    ----------------
    **kwargs
        Forwarded to plot

    Examples
    --------
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> c0 = mt.createCircle(pos=(0.0, 0.0), radius=1, nSegments=4)
    >>> l0 = mt.createPolygon([[-0.5, 0.0], [.5, 0.0]], boundaryMarker=2)
    >>> l1 = mt.createPolygon([[-0.25, -0.25], [0.0, -0.5], [0.25, -0.25]],
    ...                       interpolate='spline', addNodes=4,
    ...                       boundaryMarker=3)
    >>> l2 = mt.createPolygon([[-0.25, 0.25], [0.0, 0.5], [0.25, 0.25]],
    ...                       interpolate='spline', addNodes=4,
    ...                       isClosed=True, boundaryMarker=3)
    >>> mesh = mt.createMesh([c0, l0, l1, l2], area=0.01)
    >>> ax, _ = pg.show(mesh)
    >>> pg.viewer.mpl.drawBoundaryMarkers(ax, mesh)
    """
    ms = pg.unique(pg.sort(mesh.boundaryMarkers()[mesh.boundaryMarkers()!=0]))
    bc = kwargs.pop('bc', None)

    #cMap = plt.cm.get_cmap("Set3", len(ms))
    kwargs['lw'] = kwargs.pop('lw', 4)

    for i, m in enumerate(ms):
        bs = mesh.findBoundaryByMarker(m)
        paths = mesh.findPaths(bs)
        col = 'C' + str(i)

        for p in paths:

            xs = pg.x(mesh.nodes(p))
            ys = pg.y(mesh.nodes(p))
            path = np.array([xs, ys]).T

            ax.plot(xs, ys, color=col, **kwargs)

            center = pg.meshtools.interpolateAlongCurve(
                path, [pg.utils.cumDist(path)[-1]/2])[0]

            x = center[0]
            y = center[1]
            bbox_props = dict(boxstyle="circle,pad=0.2", fc="w", ec=col)
            txt = ax.text(x, y, str(m), color=col, va="center", ha="center",
                          zorder=20, bbox=bbox_props, fontsize=9,
                          fontdict={'weight': 'bold'})

            # clipping avoid visuablity outside axes.
            # Needed if the axes limits do not match mesh size.
            txt.set_clip_on(clipBoundaryMarkers)

            ax.plot(xs[0], ys[0], 'o', color='k')
            ax.plot(xs[-1], ys[-1], 'o', color='k')

    if bc is not None:
        pg._r(bc)

        # pos = ax.get_position()
        # print(pos)
        # print(pos.height)
        # bca = ax.figure.add_axes([0.8, 0.75, 0.1, pos.height])
        bcstr = 'BC\n'
        bcstr += '1: $p=1$\n'
        bcstr += '2: $p=0$\n'
        bcstr += '3: $v=$no-flow\n'
        bcstr += '4: $v=$no-flow\n'

        ax.text(1.02, 0.0, bcstr,
                 #transform=ax.figure.transFigure,
                #  bbox=dict(boxstyle="square",
                #    ec=(1., 0.5, 0.5),
                #    fc=(1., 0.8, 0.8),
                #    )
                 )


    # for b in mesh.boundaries():
    #     if b.marker() != 0:
    #         x = b.center()[0]
    #         y = b.center()[1]
    #         bbox_props = dict(boxstyle="circle,pad=0.1", fc="w", ec="k")
    #         ax.text(x, y, str(b.marker()), color="k", va="center",
    #                 ha="center", zorder=20, bbox=bbox_props, fontsize=9)


def drawMeshBoundaries(ax, mesh, hideMesh=False, useColorMap=False,
                       fitView=True, lw=None, color=None, **kwargs):
    """Draw mesh on ax with boundary conditions colorized.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`

    hideMesh : bool [False]
        Show only the boundary of the mesh and omit inner edges that
        separate the cells.
    useColorMap : bool[False]
        Apply the default colormap to boundaries with marker values > 0
    fitView : bool [True]
        Adjust ax limits to mesh bounding box.
    lw : float [None]
        Linewidth. When set to None then lw depends on boundary marker.
        Linewidth [0.3] for edges with marker == 0 if hideMesh is False.
    color : None
        Color for special lines. If set to None automatic "black".

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawMeshBoundaries
    >>> n = np.linspace(0, -2, 11)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> for bound in mesh.boundaries():
    ...     if not bound.rightCell():
    ...         bound.setMarker(pg.core.MARKER_BOUND_MIXED)
    ...     if bound.center().y() == 0:
    ...         bound.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
    >>> fig, ax = plt.subplots()
    >>> drawMeshBoundaries(ax, mesh)
    """
    if not mesh:
        raise Exception("drawMeshBoundaries(ax, mesh): invalid mesh")

    if mesh.dimension() != 2:
        raise Exception("No 2d mesh: dim = ", mesh.dimension())

    if mesh.nodeCount() < 2:
        raise Exception("drawMeshBoundaries(ax, mesh): to few nodes",
                        mesh.nodeCount())

    if fitView is True:
        ax.autoscale(enable=True, axis='both', tight=True)

    mesh.createNeighborInfos()

    if not hideMesh:
        drawSelectedMeshBoundaries(ax,
                                   mesh.findBoundaryByMarker(0),
                                   color=color or (0.0, 0.0, 0.0, 1.0),
                                   linewidth=lw or 0.3)

    drawSelectedMeshBoundaries(
        ax, mesh.findBoundaryByMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN),
                                      color=(0.0, 1.0, 0.0, 1.0),
                                      linewidth=lw or 1.0)
    drawSelectedMeshBoundaries(
        ax, mesh.findBoundaryByMarker(pg.core.MARKER_BOUND_MIXED),
                                      color=(1.0, 0.0, 0.0, 1.0),
                                      linewidth=lw or 1.0)

    col = color

    b0 = [b for b in mesh.boundaries() if b.marker() > 0]
    if useColorMap:
        drawSelectedMeshBoundaries(ax, b0, color=None,
                                   linewidth=lw or 1.5)
    else:
        drawSelectedMeshBoundaries(ax, b0,
                                   color=col or (0.0, 0.0, 0.0, 1.0),
                                   linewidth=lw or 1.5)

    b4 = [b for b in mesh.boundaries() if b.marker() < -4]
    drawSelectedMeshBoundaries(ax, b4,
                               color=col or (0.0, 0.0, 0.0, 1.0),
                               linewidth=lw or 1.5)

    updateAxes_(ax)


def drawPLC(ax, mesh, fillRegion=True, regionMarker=True,
            boundaryMarkers=False, showNodes=False, fitView=True, **kwargs):
    """Draw 2D PLC into given axes.

    Parameters
    ----------
    ax : mpl axe

    mesh : :gimliapi:`GIMLI::Mesh`

    fillRegion: bool [True]
        Fill the regions with default colormap.
    regionMarker: bool [True]
        Show region marker.
    boundaryMarkers: bool [False]
        Show boundary marker.
    showNodes: bool [False]
        Draw all nodes as little dots.
    fitView : bool [True]
        Adjust ax limits to mesh bounding box.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forwarded to the draw functions and mpl methods,
        respectively.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> # Create geometry definition for the modelling domain
    >>> world = mt.createWorld(start=[-20, 0], end=[20, -16],
    ...                        layers=[-2, -8], worldMarker=False)
    >>> # Create a heterogeneous block
    >>> block = mt.createRectangle(start=[-6, -3.5], end=[6, -6.0],
    ...                            marker=10,  boundaryMarker=10, area=0.1)
    >>> fig, ax = plt.subplots()
    >>> geom = world + block
    >>> _ = pg.viewer.mpl.drawPLC(ax, geom)
    """
    #    eCircles = []
    cbar = None
    if fillRegion and mesh.boundaryCount() > 2:
        tmpMesh = pg.meshtools.createMesh(mesh, quality=20, area=0)
        if tmpMesh.cellCount() == 0:
            gci = None
        else:
            markers = np.array(tmpMesh.cellMarkers())
            uniquemarkers, uniqueidx = np.unique(markers, return_inverse=True)
            gci = drawModel(ax=ax,
                            data=np.arange(len(uniquemarkers))[uniqueidx],
                            mesh=tmpMesh,
                            alpha=1,
                            linewidth=0.0,
                            tri=True,
                            snap=True,
                            )

            if regionMarker is True:
                orient = kwargs.pop('orientation', 'horizontal')
                cbar = createColorBar(gci, orientation=orient,
                                      label="Region markers")

                updateColorBar(cbar,
                               cMap=cmapFromName("Set3",
                                                 ncols=len(uniquemarkers)),
                               cMin=-0.5,
                               cMax=len(uniquemarkers) - 0.5)
                ticks = np.arange(len(uniquemarkers))

                cbar.set_ticks(ticks)
                areas = {}
                for reg in mesh.regionMarkers():
                    areas[reg.marker()] = reg.area()
                    # if kwargs.get("regionMarkers", True):
                    if regionMarker:
                        ax.plot(reg.x(), reg.y(), "mx", alpha=0.5)
                        ax.text(reg.x(), reg.y(), str(reg.marker()),
                                color="m", ha="center", va="center")

                labels = []
                for marker in uniquemarkers:
                    label = "{:d}".format(marker)
                    if marker in areas and areas[marker] > 0:
                        label += "\n$A$={:g}".format(areas[marker])
                        # label += "\n(area: %s)" % areas[marker]
                    labels.append(label)
                cbar.set_ticklabels(labels)

    else:
        gci = None
        if kwargs.pop('showBoundary', True):
            drawMeshBoundaries(ax, mesh, **kwargs)

    if showNodes:
        for n in mesh.nodes():
            col = (0.0, 0.0, 0.0, 0.5)

            if n.marker() == pg.core.MARKER_NODE_SENSOR:
                col = (0.0, 0.0, 0.0, 1.0)

            # ms = kwargs.pop('markersize', 5)
            ax.plot(n.pos()[0], n.pos()[1], 'o',
                    color=col, zorder=10, **kwargs)

    #        eCircles.append(mpl.patches.Circle((n.pos()[0], n.pos()[1])))
    #        eCircles.append(mpl.patches.Circle((n.pos()[0], n.pos()[1]), 0.1))
    #        cols.append(col)
    #    p = mpl.collections.PatchCollection(eCircles, color=cols)
    #    ax.add_collection(p)

    if regionMarker:
        for hole in mesh.holeMarker():
            ax.text(hole[0], hole[1], 'H', color='black',
                    va="center", ha="center")

    if fitView:
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_aspect('equal')

    updateAxes_(ax)
    # if cbar is None:
    #     return gci
    return gci, cbar


def _createCellPolygon(cell):
    """Utility function to polygon for cell shape to be used by MPL."""
    if cell.shape().nodeCount() == 3:
        return list(zip([cell.node(0).x(), cell.node(1).x(),
                         cell.node(2).x()],
                        [cell.node(0).y(), cell.node(1).y(),
                         cell.node(2).y()]))
    elif cell.shape().nodeCount() == 4:
        return list(zip([cell.node(0).x(), cell.node(1).x(),
                         cell.node(2).x(), cell.node(3).x()],
                        [cell.node(0).y(), cell.node(1).y(),
                         cell.node(2).y(), cell.node(3).y()]))

    pg.warn("Unknown shape to patch: ", cell)


def createMeshPatches(ax, mesh, rasterized=False, verbose=True):
    """Utility function to create 2d mesh patches within a given ax."""
    import matplotlib as mpl
    if not mesh:
        pg.error("drawMeshBoundaries(ax, mesh): invalid mesh:", mesh)
        return

    if mesh.nodeCount() < 2:
        pg.error("drawMeshBoundaries(ax, mesh): to few nodes:", mesh)
        return

    pg.tic()
    polys = [_createCellPolygon(c) for c in mesh.cells()]

    patches = mpl.collections.PolyCollection(polys, picker=True,
                                             rasterized=rasterized)

    if verbose:
        pg.info("Creation of mesh patches took = ", pg.toc())

    return patches


def createTriangles(mesh):
    """Generate triangle objects for later drawing.

    Creates triangle for each 2D triangle cell or 3D boundary.
    Quads will be split into two triangles.
    Result will be cached into mesh._triData.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        2D mesh or 3D mesh

    Returns
    -------
    x : numpy array
        x position of nodes
    y : numpy array
        x position of nodes
    triangles : numpy array Cx3
        cell indices for each triangle, quad or boundary face
    z : numpy array
        z position for given indices
    dataIdx : list of int
        List of indices for a data array
    """
    if hasattr(mesh, '_triData'):
        if hash(mesh) == mesh._triData[0]:
            return mesh._triData[1:]

    x = pg.x(mesh)
    y = pg.y(mesh)
    z = pg.z(mesh)
    #    x.round(1e-1)
    #    y.round(1e-1)

    if mesh.dim() == 2:
        ents = mesh.cells()
    else:
        ents = mesh.boundaries(mesh.boundaryMarkers() != 0)
        if len(ents) == 0:
            for b in mesh.boundaries():
                if b.leftCell() is None or b.rightCell() is None:
                    ents.append(b)

    triangles = []
    dataIdx = []

    for c in ents:
        triangles.append([c.node(0).id(), c.node(1).id(), c.node(2).id()])
        dataIdx.append(c.id())

        if c.shape().nodeCount() == 4:
            triangles.append([c.node(0).id(), c.node(2).id(), c.node(3).id()])
            dataIdx.append(c.id())

    mesh._triData = [hash(mesh), x, y, triangles, z, dataIdx]
    return x, y, triangles, z, dataIdx


def drawField(ax, mesh, data=None, levels=None, nLevs=5,
              cMin=None, cMax=None, nCols=None, logScale=False, fitView=True,
              **kwargs):
    """Draw mesh with scalar field data.

    Draw scalar field into MPL axes using matplotlib triplot.
    Only for triangle/quadrangle meshes currently

    Parameters
    ----------
    ax : Matplotlib axis object

    mesh : :gimliapi:`GIMLI::Mesh`
        2D mesh
    data: iterable
        Scalar field values. Can be of length mesh.cellCount()
        or mesh.nodeCount().
    levels : iterable of type float
        Values for contour lines. If empty auto generated from nLevs.
    nLevs : int
        Number of contour levels based on cMin, cMax and logScale.
    nCols : int
        Number of color, can be different from nLevs.
    cMin : float [None]
        Minimal contour value. If None min(data).
    cMax : float [None]
        Maximal contour value. If None max(data).
    logScale : bool [False]
        Levels and colors distributes with logarithmic scale.
    fitView : bool [True]
        Adjust ax limits to mesh bounding box.

    Keyword Arguments
    -----------------
    shading: 'flat' | 'gouraud'

    fillContour: [True]

    contourLines: [True]
    **kwargs
        Additional kwargs forwarded to ax.tripcolor,
        ax.tricontour, ax.tricontourf

    symLevels: bool [False]
        Symmetric levels if no cMap is given.

    Returns
    -------
    gci : image object
        The current image object useful for post color scaling

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawField
    >>> n = np.linspace(0, -2, 11)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> nx = pg.x(mesh.positions())
    >>> ny = pg.y(mesh.positions())
    >>> data = np.cos(1.5 * nx) * np.sin(1.5 * ny)
    >>> fig, ax = plt.subplots()
    >>> tri = drawField(ax, mesh, data)
    """
    x, y, triangles, _, dataIndex = createTriangles(mesh)
    if len(data) == mesh.cellCount():
        z = data[dataIndex]
    else:
        z = data

    gci = None

    if levels is None:
        levels = autolevel(data, nLevs,
                           zMin=cMin, zMax=cMax, logScale=logScale,
                           symLevels=kwargs.pop('symLevels', False))
    else:
        nCols = len(levels)-1

    if len(z) == len(triangles) and len(data) != mesh.nodeCount():
        shading = kwargs.pop('shading', 'flat')

        # bounds = np.linspace(levels[0], levels[-1], nLevs)
        # norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        if shading == 'gouraud':
            z = pg.meshtools.cellDataToNodeData(mesh, data)
            gci = ax.tripcolor(x, y, triangles, z,
                               shading=shading, **kwargs)
        else:
            gci = ax.tripcolor(x, y, triangles, facecolors=z,
                               shading=shading, **kwargs)

    elif len(z) == mesh.nodeCount():
        shading = kwargs.pop('shading', None)

        if shading is not None:
            gci = ax.tripcolor(x, y, triangles, z, shading=shading, **kwargs)
        else:

            fillContour = kwargs.pop('fillContour', True)
            contourLines = kwargs.pop('contourLines', True)

            if fillContour is True:
                # add outer climits to fill lower and upper too
                levs = np.array(levels)

                # if min(z) < min(levels):
                #     levs = np.hstack([min(z), levs])

                # if max(z) > max(levels):
                #     levs = np.hstack([levs, max(z)])

                if nCols is not None:
                    if logScale:
                        levs = np.geomspace(min(levels), max(levels), nCols+1)
                    else:
                        levs = np.linspace(min(levels), max(levels), nCols+1)

                gci = ax.tricontourf(x, y, triangles, z,
                                     # antialiased=True, # not allways nice
                                     levels=levs, **kwargs)

            if contourLines is True:
                ax.tricontour(x, y, triangles, z, levels=levels,
                              colors=kwargs.pop('colors', ['0.5']),
                              alpha=kwargs.pop('alpha', 0.5),
                              linewidths=kwargs.pop('lw', 1.0),
                              **kwargs)
    else:
        gci = None
        raise Exception("Data size does not fit mesh size: ", len(z),
                        mesh.cellCount(), mesh.nodeCount())

    # we should ne adapt cols here at all -- test remove
    # if gci and cMin and cMax:
    #     gci.set_clim(cMin, cMax)

    if fitView is True:
        ax.autoscale(enable=True, axis='both', tight=True)
        ax.set_aspect('equal')

    updateAxes_(ax)
    return gci


def drawStreamLines(ax, mesh, u, nx=25, ny=25, **kwargs):
    """Draw streamlines for the gradients of field values u on a mesh.

    The matplotlib routine streamplot needs equidistant spacings so
    we interpolate first on a grid defined by nx and ny nodes.
    Additionally arguments are piped to streamplot.

    This works only for rectangular regions.
    You should use pg.viewer.mpl.drawStreams, which is more comfortable and
    more flexible.

    Parameters
    ----------
    ax : mpl axe

    mesh : :gimliapi:`GIMLI::Mesh`
        2D mesh
    u : iterable float
        Scalar data field.

    """
    X, Y = np.meshgrid(
        np.linspace(mesh.xmin(), mesh.xmax(), nx),
        np.linspace(mesh.ymin(), mesh.ymax(), ny))

    U = X.copy()
    V = X.copy()

    for i, row in enumerate(X):
        for j in range(len(row)):
            p = [X[i, j], Y[i, j]]
            gr = [0.0, 0.0]
            c = mesh.findCell(p)
            if c:
                gr = c.grad(p, u)

            U[i, j] = -gr[0]
            V[i, j] = -gr[1]

    gci = ax.streamplot(X, Y, U, V, **kwargs)

    updateAxes_(ax)
    return gci


def drawStreamLine(ax, mesh, c, data, dataMesh=None, linewidth=1.0,
                    dropTol=0.0, **kwargs):
    """Draw a single streamline.

    Draw a single streamline into a given mesh for given data stating at
    the center of cell c.
    The Streamline will be enlarged until she reached a cell that
    already contains a streamline.

    TODO
        linewidth and color depends on absolute velocity
        or background color saturation

    Parameters
    ----------
    ax : matplotlib.ax
        ax to draw into
    mesh : :gimliapi:`GIMLI::Mesh`
        2d mesh
    c : :gimliapi:`GIMLI::Cell`
        Start point is c.center()
    data : iterable float | [float, float]
        If data is an array (per cell or node) gradients are calculated
        otherwise the data will be interpreted as vector field per nodes or
        cell centers.
    dataMesh : :gimliapi:`GIMLI::Mesh` [None]
        Optional mesh for the data. If you want high resolution
        data to plot on coarse draw mesh.
    linewidth : float [1.0]
        Streamline linewidth
    dropTol : float [0.0]
        Don't draw stream lines with velocity lower than drop tolerance.

    Keyword Arguments
    -----------------
    **kwargs
        arrowSize: int
            Size of the arrow's head.
        arrowColor: str
            Color of the arrow's head.
        Additional kwargs are forwarded to mpl.LineCollection, mpl.Polygon
    """
    import matplotlib as mpl
    x, y, v = streamline(mesh, data, startCoord=c.center(), dLengthSteps=5,
                         dataMesh=dataMesh, maxSteps=10000, verbose=False,
                         coords=[0, 1])

    if 'color' not in kwargs:
        kwargs['color'] = 'black'

    arrowSize = kwargs.pop('arrowSize', 12)
    arrowColor = kwargs.pop('arrowColor', kwargs.get('color'))

    lines = None

    if len(x) > 2:
        points = np.array([x, y]).T.reshape(-1, 1, 2)

        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lwidths = pg.Vector(len(v), kwargs.pop('lw', linewidth))
        lwidths[pg.find(pg.Vector(v) < dropTol)] = 0.0

        lines = mpl.collections.LineCollection(segments,
                                               linewidths=lwidths, **kwargs)
        ax.add_collection(lines)

        # probably the limits are wrong without plot call
        # lines = ax.plot(x, y, **kwargs)
        # updateAxes_(ax, lines)
        # ax.plot(x, y, '.-', color='black', **kwargs)
    if len(x) > 3:
        xmid = int(len(x) / 2)
        ymid = int(len(y) / 2)
        dx = x[xmid + 1] - x[xmid]
        dy = y[ymid + 1] - y[ymid]
        c = mesh.findCell([x[xmid], y[ymid]])

        if v[xmid] > dropTol:

            absArrowSize = True
            if absArrowSize:
                ax.annotate('', xytext=(x[xmid]-dx, y[ymid]-dy),
                            xy=(x[xmid], y[ymid]),
                            arrowprops=dict(arrowstyle="-|>",
                                            color=arrowColor, lw=0),
                            size=arrowSize, **kwargs)
            else:
                ax.arrow(x[xmid], y[ymid], dx, dy,
                         shape='full', lw=0,
                         length_includes_head=True,
                         fc=arrowColor,
                         head_width=.35, **kwargs)

            # dx90 = -dy
            # dy90 = dx
            # aLen = 3
            # aWid = 1
            # xy = list(zip([x[xmid] + dx90*aWid, x[xmid] + dx*aLen,
            #                x[xmid] - dx90*aWid],
            #               [y[ymid] + dy90*aWid, y[ymid] + dy*aLen,
            #                y[ymid] - dy90*aWid]))

            # arrow = mpl.patches.Polygon(xy, ls=None, lw=0, closed=True,
            #                             **kwargs)
            # ax.add_patch(arrow)

    return lines


def drawStreams(ax, mesh, data, startStream=3, coarseMesh=None, quiver=False,
                fitView=True, **kwargs):
    """Draw streamlines based on an unstructured mesh.

    Every cell contains only one streamline and every new stream line
    starts in the center of a cell. You can alternatively provide a second mesh
    with coarser mesh to draw streams for.

    Parameters
    ----------
    ax : matplotlib.ax
        ax to draw into
    mesh : :gimliapi:`GIMLI::Mesh`
        2d mesh
    data : iterable float | [float, float] | pg.PosVector
        If data is an array (per cell or node) gradients are calculated
        otherwise the data will be interpreted as vector field per nodes or
        cell centers.
    startStream : int
        variate the start stream drawing, try values from 1 to 3 what every
        you like more.
    coarseMesh : :gimliapi:`GIMLI::Mesh`
        Instead of draw a stream for every cell in mesh, draw a streamline
        segment for each cell in coarseMesh.
    quiver : bool [False]
        Draw arrows instead of streamlines.

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forwarded to axe.quiver, drawStreamLine

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawStreams
    >>> n = np.linspace(0, 1, 10)
    >>> mesh = pg.createGrid(x=n, y=n)
    >>> nx = pg.x(mesh.positions())
    >>> ny = pg.y(mesh.positions())
    >>> data = np.cos(1.5 * nx) * np.sin(1.5 * ny)
    >>> fig, ax = plt.subplots()
    >>> drawStreams(ax, mesh, data, color='red')
    >>> drawStreams(ax, mesh, data, dropTol=0.9)
    >>> drawStreams(ax, mesh, pg.solver.grad(mesh, data),
    ...             color='green', quiver=True)
    >>> ax.set_aspect('equal')
    >>> pg.wait()
    """
    viewMesh = None
    dataMesh = None

    if quiver:

        x = None
        y = None
        u = None
        v = None

        if len(data) == mesh.nodeCount():
            x = pg.x(mesh.positions())
            y = pg.y(mesh.positions())
        elif len(data) == mesh.cellCount():
            x = pg.x(mesh.cellCenters())
            y = pg.y(mesh.cellCenters())
        elif len(data) == mesh.boundaryCount():
            x = pg.x(mesh.boundaryCenters())
            y = pg.y(mesh.boundaryCenters())

        if isinstance(data, pg.PosVector):
            u = pg.x(data)
            v = pg.y(data)
        else:
            u = data[:, 0]
            v = data[:, 1]

        q = ax.quiver(x, y, u, v, **kwargs)

        label = kwargs.pop('label', None)
        if label is not None:
            ax.quiverkey(q, 0.05, -0.1, 0.1, label, labelpos='E',
                         coordinates='axes')

        updateAxes_(ax)
        return

    if coarseMesh is not None:
        viewMesh = coarseMesh
        dataMesh = mesh
        dataMesh.createNeighborInfos()
    else:
        viewMesh = mesh

    viewMesh.createNeighborInfos()

    for c in viewMesh.cells():
        c.setValid(True)

    if startStream == 1:
        # start a stream from each boundary cell
        for y in np.linspace(viewMesh.ymin(), viewMesh.ymax(), 100):
            c = viewMesh.findCell(
                [(viewMesh.xmax() - viewMesh.xmax()) / 2.0, y])
            if c is not None:
                if c.valid():
                    drawStreamLine(ax, viewMesh, c, data, dataMesh, **kwargs)

    elif startStream == 2:
        # start a stream from each boundary cell
        for x in np.linspace(viewMesh.xmin(), viewMesh.xmax(), 100):
            c = viewMesh.findCell(
                [x, (viewMesh.ymax() - viewMesh.ymax()) / 2.0])
            if c is not None:
                if c.valid():
                    drawStreamLine(ax, viewMesh, c, data, dataMesh, **kwargs)

    elif startStream == 3:
        # start a stream from each boundary cell
        for b in viewMesh.findBoundaryByMarker(1, 99):
            c = b.leftCell()
            if c is None:
                c = b.rightCell()

            if c.valid():
                drawStreamLine(ax, viewMesh, c, data, dataMesh, **kwargs)

    # start a stream from each unused cell
    for c in viewMesh.cells():
        if c.valid():
            drawStreamLine(ax, viewMesh, c, data, dataMesh, **kwargs)

    for c in viewMesh.cells():
        c.setValid(True)

    label = kwargs.pop('label', None)

    if label is not None:
        from matplotlib.lines import Line2D
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        arrowSize = kwargs.pop('arrowSize', 12)
        arrowColor = kwargs.pop('arrowColor', kwargs.get('color'))

        yOff = -26
        ## streamline
        ax.annotate('', xy=(3, yOff), xycoords='axes points',
                    xytext=(25, yOff),
                    arrowprops=dict(arrowstyle="->",
                                    color=arrowColor,
                                    lw=1),
                    size=0 )

        ## arrow head
        ax.annotate('', xy=(4, yOff), xycoords='axes points',
                    xytext=(20, yOff),
                    arrowprops=dict(arrowstyle="<|-",
                                    color=arrowColor,
                                    lw=0),
                    size=arrowSize )
        ## text
        ax.annotate(label, xy=(30, yOff),
                    xycoords='axes points',
                    horizontalalignment='left',
                    verticalalignment='center')

    updateAxes_(ax)


def drawSensors(ax, sensors, diam=None, coords=None, **kwargs):
    """Draw sensor positions as black dots with a given diameter.

    Parameters
    ----------
    ax : mpl axe instance

    sensors : vector or list of Pos
        List of positions to plot.
    diam : float [None]
        Diameter (absolute in m) of circles (None leads to point distance by 4).
    coords: (int, int) [0, 1]
        Coordinates to take (usually x and y).

    Keyword Arguments
    -----------------
    **kwargs
        Additional kwargs forwarded to mpl.PatchCollection, mpl.Circle

    sensorMarker:
        Also 'sm'. Set marker style: 'o' Circle, 'v' Triangle pointing down.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawSensors
    >>> sensors = np.random.rand(5, 2)
    >>> fig, ax = pg.plt.subplots()
    >>> drawSensors(ax, sensors, diam=0.02, coords=[0, 1])
    >>> ax.set_aspect('equal')
    >>> pg.wait()
    """
    import matplotlib as mpl
    if coords is None:
        coords = [0, 2]

        if len(sensors[0]) == 2 or \
                not pg.core.zVari(sensors) and sensors[0][2] == 0.0:
            coords = [0, 1]


    if diam is None:
        eSpacing = pg.Pos(sensors[0]).distance(sensors[1])
        diam = eSpacing / 2.5

    eSensors = []
    sm = kwargs.pop('sm', kwargs.pop('sensorMarker', 'o'))
    for e in sensors:
        x = e[coords[0]]
        y = e[coords[1]]
        if sm == 'o':
            eSensors.append(mpl.patches.Circle((x,y), diam/2,
                                               **kwargs))
        else:
            eSensors.append(mpl.patches.Polygon(([x, y],
                                                 [x+diam/1.5/1.4, y+diam/1.5],
                                                 [x-diam/1.5/1.4, y+diam/1.5]
                                                 ),
                                                closed=True,
                                                **kwargs))

    p = mpl.collections.PatchCollection(eSensors, clip_on=False, **kwargs)
    p.set_zorder(100)
    ax.add_collection(p)

    updateAxes_(ax)
    #ax.set(ylim=[None, y+diam])


def _createParameterContraintsLines(mesh, cMat, cWeights=None):
    """Create line segments representing constrains.
    """
    C = None

    if isinstance(cMat, pg.matrix.SparseMapMatrix):
        # TODO super hackish .. clean me up!!

        tmp = pg.optImport('tempfile')
        _, tmpFile = tmp.mkstemp(suffix='.matrix')
        C = pg.Matrix()
        cMat.save(tmpFile)
        pg.core.loadMatrixCol(C, tmpFile)

        try:
            import os
            os.remove(tmpFile)
        except Exception as e:
            pg.error(e)
            print("can't remove:", tmpFile)

    else:
        pg.critical('implementme')
        C = cMat

    cellList = dict()
    for c in mesh.cells():
        pID = c.marker()
        if pID not in cellList:
            cellList[pID] = []
        cellList[pID].append(c)

    paraCenter = dict()
    for pID, vals in list(cellList.items()):
        p = pg.Pos(0.0, 0.0, 0.0)
        for c in vals:
            p += c.center()
        p /= float(len(vals))
        paraCenter[pID] = p

    nConstraints = cMat.rows()

    if cWeights is None:
        cWeights = np.ones(nConstraints)

    start = []
    end = []
    #    swatch = pg.core.Stopwatch(True)  # not used
    i = -1
    while i < C[0].size():
        cID = int(C[0][i])

        a = C[1][i]
        b = None

        if i < C[0].size()-1:
            if C[0][i+1] == cID:
                b = C[1][i+1]
                i += 1
        if b is not None:
            if cWeights[cID] > 0:
                p1 = paraCenter[a]
                p2 = paraCenter[b]

                if cWeights is not None:
                    pa = pg.Pos(p1 + (p2 - p1)/2. * (1. - cWeights[cID]))
                    pb = pg.Pos(p2 + (p1 - p2)/2. * (1. - cWeights[cID]))
                else:
                    pa = p1
                    pb = p2

                start.append(pa)
                end.append(pb)
        else:
            start.append(paraCenter[a])
            end.append(paraCenter[a])

        i += 1

    return start, end


def drawParameterConstraints(ax, mesh, cMat, cWeights=None):
    """Draw inter parameter constraints between cells.

    Parameters
    ----------
    ax : MPL axes

    mesh : :gimliapi:`GIMLI::Mesh`
        2d mesh
    cMat : :gimliapi:`GIMLI::SparseMapMatrix`
        ConstraintsMatrix
    cWeights : iterable float
        Weights for all constraints. Need to have a lengths == cMat.rows()
    """
    import matplotlib as mpl
    start, end = _createParameterContraintsLines(mesh, cMat, cWeights)

    lines = []
    colors = []
    linewidths = []
    col = (0.0, 0.0, 1.0, 1.0)
    for i, _ in enumerate(start):
        if start[i] == end[i]:
            ax.plot(start[i].x(), end[i].y(), '.', color=col, markersize=2)
        else:
            lines.append(list(zip([start[i].x(), end[i].x()],
                                  [start[i].y(), end[i].y()])))

            linewidths.append(0.5)
        colors.append(col)

    lc = mpl.collections.LineCollection(lines, antialiaseds=True)
    lc.set_color(colors)
    lc.set_linewidth(linewidths)
    ax.add_collection(lc)

    updateAxes_(ax)
