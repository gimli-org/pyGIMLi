#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib
import pygimli as pg
from .utils import pgMesh2pvMesh

from pygimli.viewer.mpl.colorbar import cmapFromName

pv = pg.optImport('pyvista', requiredFor="properly visualize 3D data")


def drawMesh(ax, mesh, notebook=False, **kwargs):
    """Draw a mesh into a given plotter.

    Parameters
    ----------
    ax: pyvista.Plotter [optional]
        The plotter to draw everything. If none is given, one will be created.
    mesh: pg.Mesh
        The mesh to show.
    notebook: bool [False]
        Sets the plotter up for jupyter notebook/lab.
    cMap: str ['viridis']
        The colormap string.
    bc: pyvista color ['#EEEEEE']
        Background color.
    style: str['surface']
        Possible options: "surface", "wireframe", "points"
    label: str
        Data to be plotted. If None the first is taken.

    Returns
    -------
    ax: pyvista.Plotter [optional]
        The plotter
    """
    # sort out a few kwargs to not confuse the plotter initialization
    opacity = kwargs.pop('alpha', kwargs.pop('opacity', 1))
    cMap = kwargs.pop('cMap', None)
    color = kwargs.pop('color', None)
    style = kwargs.pop('style', 'surface')
    returnActor = kwargs.pop('returnActor', False)
    showMesh = kwargs.pop('showMesh', False)
    grid = kwargs.pop('grid', False)
    colorBar = kwargs.pop('colorBar', style != 'wireframe')

    if pv.BUILDING_GALLERY:
        bc = "#ffffff"
    else:
        bc = kwargs.pop('bc', '#EEEEEE')  # background color

    lw = kwargs.pop('line_width', 0.1)
    filt = kwargs.pop('filter', {})
    log_scale = kwargs.pop("logScale", False)
    clim = None

    if "cMin" in kwargs and "cMax" in kwargs:
        clim = [kwargs.pop("cMin"), kwargs.pop("cMax")]

    # show_edges = kwargs.pop('show_edges', True)  # not used
    # name = kwargs.pop('name', 'Mesh')  # not used

    if isinstance(mesh, pg.Mesh):
        mesh = pgMesh2pvMesh(mesh)

    dataName = kwargs.pop('label', list(mesh.cell_data.keys())[0])

    try:
        theme = pv.themes.Theme()
    except:  # older pv versions
        theme = pv.themes.DefaultTheme()

    theme.background = bc
    # seems to be broken .. results on pure black screens on some machines
    # theme.antialiasing = True

    theme.font.color = 'k'

    if ax is None:
        ax = pv.Plotter(notebook=notebook, theme=theme,
                        **kwargs)

    if grid is True:
        pass  # implementme

    ax.show_bounds(all_edges=True, minor_ticks=True)
    ax.add_axes()

    for k, fi in filt.items():
        if k.lower() == 'clip':
            if isinstance(mesh, pv.core.pointset.UnstructuredGrid):
                fi.setdefault('crinkle', True)

            mesh = mesh.clip(**fi)
        elif k.lower() == 'threshold':
            mesh = mesh.threshold(**fi)
        elif k.lower() == 'slice':
            mesh = mesh.slice(**fi)
        else:
            pg.error('filter:', k, 'not yet implemented')

    _actor = ax.add_mesh(mesh,  # type: pv.UnstructuredGrid
                         scalars=dataName,
                         cmap=cMap,
                         color=color,
                         style=style,
                         show_edges=showMesh,
                         line_width=lw,
                         # edge_color='white',
                         show_scalar_bar=colorBar,
                         opacity=opacity,
                         clim=clim,
                         log_scale=log_scale
                         )

    if returnActor:
        return ax, _actor
    else:
        return ax


def drawModel(ax=None, mesh=None, data=None, **kwargs):
    """Draw a mesh with given data.

    Parameters
    ----------
    ax: pyvista.Plotter [None]
        Pyvista's basic Plotter to add the mesh to.
    mesh: pg.Mesh
        The Mesh to plot.
    data: iterable
        Data that should be displayed with the mesh.

    Returns
    -------
    ax: pyvista.Plotter [optional]
        The plotter
    """
    defaultCMap = kwargs.pop('cMap', 'viridis')
    dataName = kwargs.pop('label', None)

    if all(v is None for v in [ax, mesh, data]):
        pg.critical("At least mesh or data should not be None")
        return None

    if kwargs.pop('markers', False) is True:
        # show boundary mesh with markers
        data = mesh.boundaryMarkers()
        defaultCMap = cmapFromName("Set3", ncols=max(1, len(pg.unique(data))))
        dataName = 'Boundary Marker'
        mesh = pgMesh2pvMesh(mesh, data, dataName, boundaries=True)
    else:

        if data is not None or len(mesh.dataMap()) != 0:
            kwargs.setdefault('style', 'surface')
            kwargs['color'] = None
        if dataName is None and data is not None:
            if len(data) == mesh.cellCount():
                dataName = 'Cell data'
            elif len(data) == mesh.nodeCount():
                dataName = 'Node data'

        mesh = pgMesh2pvMesh(mesh, data, dataName)

    kwargs['cMap'] = defaultCMap
    kwargs['label'] = dataName

    return drawMesh(ax, mesh, **kwargs)


def drawSensors(ax, sensors, diam=0.01, color='grey', **kwargs):
    """
    Draw the sensor positions to given mesh or the the one in given plotter.

    Parameters
    ----------
    ax: pyvista.Plotter
        The plotter to draw everything. If none is given, one will be created.
    sensors: iterable
        Array-like object containing tuple-like (x, y, z) positions.
    diam: float [0.01]
        Radius of sphere markers.
    color: str ['grey']
        Color of sphere markers.

    Returns
    -------
    ax: pyvista.Plotter
        The plotter containing the mesh and drawn electrode positions.
    """
    for pos in sensors:
        s = pv.Sphere(radius=diam / 2, center=pos)
        ax.add_mesh(s, color=color, **kwargs)

    return ax


def drawSlice(ax, mesh, normal=[1, 0, 0], **kwargs):
    """Draw a slice in a 3D mesh for given pygimli mesh.

    Parameters
    ----------
    ax: pyvista.Plotter
        The Plotter to draw on.
    mesh: pg.Mesh
        The mesh to take the slice out of.
    normal: list [[1, 0, 0]]
        Coordinates to orientate the slice.

    Returns
    -------
    ax: pyvista.Plotter
        The plotter containing the mesh and drawn electrode positions.

    Keyword arguments passed to pyvista.slice
    -----------------------------------------
    normal: [float, float, float] | str
        normal vector constructing the slice
    origin: [float, float, float]
        origin for the slice (by default mesh center)
    generate_triangles: bool [False]
        generate triangle mesh
    contour: bool [False]
        draw contours instead

    Keyword arguments passed to pyvista.add_mesh
    --------------------------------------------
    cmap|cMap : str [None]
        colormap
    log_scale|logScale : bool [False]
        use logarithmic colormap scaling
    clim : [float, float]
        color limits as tuple/list
    cMin, cMax : float
        color limits in pg style

    More information at
    https://docs.pyvista.org/api/core/_autosummary/pyvista.CompositeFilters.slice.html
    """
    label = kwargs.pop('label', None)
    data = kwargs.pop('data', None)
    origin = kwargs.pop('origin', None)
    generate_triangles = kwargs.pop('generate_triangles', False)
    contour = kwargs.pop('contour', False)
    kwargs.setdefault('cmap', kwargs.pop('cMap', None))
    kwargs.setdefault('log_scale', kwargs.pop('logScale', False))
    if 'cMin' in kwargs and 'cMax' in kwargs:
        kwargs.setdefault('clim', [kwargs.pop('cMin'), kwargs.pop('cMax')])
    pvmesh = pgMesh2pvMesh(mesh, data, label)

    try:
        single_slice = pvmesh.slice(normal, origin=origin, contour=contour,
                                    generate_triangles=generate_triangles)
    except AssertionError as e:
        # 'contour' kwarg only works with point data and breaks execution
        pg.error(e)
    else:
        # REVIEW: bounds and axes might be confused with the outline..?!
        outline = pvmesh.outline()
        ax.add_mesh(outline, color="k")
        ax.add_mesh(single_slice, **kwargs)

    return ax


def drawStreamLines(ax, mesh, data, label=None, radius=0.01, **kwargs):
    """
    Draw streamlines of given data.

    PyVista streamline needs a vector field of gradient data per cell.

    Parameters
    ----------
    ax: pyvista.Plotter [None]
        The plotter that should be used for visualization.
    mesh: pyvista.UnstructuredGrid|pg.Mesh [None]
        Structure to plot the streamlines in to.
        If pv grid a check is performed if the data set is already contained.
    data: iterable [None]
        Values used for streamlining.
    label: str
        Label for the data set. Will be searched for within the data.
    radius: float [0.01]
        Radius for the streamline tubes.

    Note
    ----
    All kwargs will be forwarded to pyvistas streamline filter:
    https://docs.pyvista.org/api/core/_autosummary/pyvista.DataSetFilters.streamlines.html
    """
    if label is None:
        label = 'grad'

    if isinstance(mesh, pg.Mesh):

        # create gradient of cell data if not provided
        if np.ndim(data) == 1:
            grad = pg.solver.grad(mesh, data)  # should be -grad
        else:
            grad = data

        # ensure that it's point/node data in the mesh
        if len(grad) == mesh.cellCount():
            grad = pg.meshtools.cellDataToNodeData(mesh, grad)

        # add data to the mesh and convert to pyvista grid
        mesh = pgMesh2pvMesh(mesh, grad.T, label)

    elif isinstance(mesh, pv.UnstructuredGrid):
        if label not in mesh.point_arrays:  # conversion needed
            mesh.cell_data_to_point_data()

    if label is None:
        label = list(mesh.point_arrays.keys())[0]

    # kwargs['vectors'] = label

    streams = mesh.streamlines(vectors=label, **kwargs)
    ax.add_mesh(streams.tube(radius=radius), show_scalar_bar=False)
