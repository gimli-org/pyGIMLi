#!/usr/bin/env python
"""Plot 3D mesh using pyvista."""
import pygimli as pg

pyvista = pg.optImport("pyvista", requiredFor="properly visualize 3D data")
trame = pg.optImport(
    "trame",
    requiredFor="use interactive 3D visualizations within Jupyter notebooks",
)

if pyvista is None:
    view3Dcallback = "showMesh3DFallback"
else:
    view3Dcallback = "showMesh3DVista"
    vers_users = pyvista.__version__
    vers_userf = float(pyvista.__version__[::-1].replace(".", "", 1)[::-1])
    vers_needs = "0.34.0"
    vers_needf = 0.340

    if vers_userf < vers_needf:
        pg.warn(f"Please consider updating PyVista to at least {vers_needs}")

    from .draw import drawModel


def showMesh3D(mesh, data, **kwargs):
    """Call the defined function to show the 3D object."""
    if pg.rc["view3D"] == "fallback":
        return showMesh3DFallback(mesh, data, **kwargs)

    return globals()[view3Dcallback](mesh, data, **kwargs)


def showMesh3DFallback(mesh, data, **kwargs):
    """Plot the 3D object sketchy."""
    ax = kwargs.pop("ax", None)
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if ax is None or not isinstance(ax, Axes3D):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d', proj_type="persp")
        #ax = fig.add_subplot(projection='3d', proj_type="ortho")

    if mesh.boundaryCount() > 0:
        x, y, tri, z, dataIndex = pg.viewer.mpl.createTriangles(mesh)
        ax.plot_trisurf(x, y, tri, z)
    else:
        if mesh.nodeCount() < 1e4:
            x = pg.x(mesh.positions())
            y = pg.y(mesh.positions())
            z = pg.z(mesh.positions())
            ax.scatter(x, y, z, "ko")
    ax.set_title("Fallback, install pyvista for proper 3D visualization")

    return ax, None


def showMesh3DVista(mesh, data=None, **kwargs):
    """Make use of the actual 3D visualization tool kit.

    Parameters
    ----------
    data: pg.Vector or np.ndarray
        Dictionary of cell values, sorted by key. The values need to be
        numpy arrays.

    Returns
    -------
    plotter: pyvista.Plotter
        The plotter from pyvista.
    gui: Show3D [None]
        The small gui based on pyvista. Note that this is returned as 'None'
        if gui is passed as 'False'.

    Note
    ----
    Not having PyQt5 installed results in displaying the first key
    (and values) from the dictionary.
    """
    # for compatibility remove show kwargs that are not needed
    kwargs.pop("figsize", False)

    hold = kwargs.pop("hold", False)
    cMap = kwargs.pop("cMap", "viridis")
    notebook = kwargs.pop("notebook", pg.isNotebook())

    # For sphinx builds and non-interactive agg backend
    if not notebook and not pg.viewer.mpl.isInteractive():
        notebook = False
        kwargs["backend"] = None

    gui = kwargs.pop("gui", False)

    # GUI tmp deactivated
    if gui:
        pg.error("pyqt show gui currently not maintained")
        return None

    plotter = drawModel(kwargs.pop("ax", None),
                        mesh, data, notebook=notebook, cMap=cMap,
                        #offScreen=True,
                        **kwargs)

    # seems top be broken on some machines
    if kwargs.get("aa", False):
        plotter.enable_anti_aliasing()

    if notebook and trame and not pyvista.BUILDING_GALLERY:
        # monkeypatch show of this plotter instance so we can use multiple
        # backends and only plotter.show() .. whoever this needs.
        plotter.__show = plotter.show
        plotter.show = lambda *args, **kwargs: plotter.__show(
            *args, jupyter_backend="trame", **kwargs
        )
    elif pyvista.BUILDING_GALLERY:
        plotter.__show = plotter.show
        plotter.show = lambda *args, **kwargs: plotter.__show(
            *args, jupyter_backend="html", **kwargs
        )

    if hold is False:
        plotter.show()

    # , None to keep compatibility
    return plotter, None
