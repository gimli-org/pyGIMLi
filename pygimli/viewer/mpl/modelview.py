"""Model viewer functions."""
import numpy as np

import pygimli as pg
from pygimli.utils import rndig

from .colorbar import setMappableData
from .utils import updateAxes as updateAxes_


def showModel1D(*args, **kwargs):
    """Show 1D model."""
    if "ax" in kwargs:
        ax = kwargs.pop("ax")
    else:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots()

    drawModel1D(ax, *args, **kwargs)
    ax.invert_yaxis()
    return ax, _

def drawModel1D(ax, thickness=None, values=None, model=None, depths=None,
                plot='plot',
                xlabel=r'Resistivity $(\Omega$m$)$', zlabel='Depth (m)',
                z0=0, zmax=None,
                **kwargs):
    """Draw 1d block model into axis ax.

    Draw 1d block model into axis ax defined by values and thickness vectors
    using plot function.
    For log y cases, z0 should be set > 0 so that the default becomes 1.

    Parameters
    ----------
    ax : mpl axes
        Matplotlib Axes object to plot into.

    values : iterable [float]
        [N] Values for each layer plus lower background.

    thickness : iterable [float]
        [N-1] thickness for each layer. Either thickness or depths must be set.

    depths : iterable [float]
        [N-1] Values for layer depths (positive z-coordinates).
        Either thickness or depths must be set.

    model : iterable [float]
        Shortcut to use default model definition.
        thks = model[0:nLay]
        values = model[nLay:]

    plot : string
        Matplotlib plotting function.
        'plot', 'semilogx', 'semilogy', 'loglog'

    xlabel : str
        Label for x axis.

    ylabel : str
        Label for y axis.

    z0 : float
        Starting depth in m

    **kwargs : dict()
        Forwarded to the plot routine

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import pygimli as pg
    >>> # plt.style.use('ggplot')
    >>> thk = [1, 4, 4]
    >>> res = np.array([10., 5, 15, 50])
    >>> fig, ax = plt.subplots()
    >>> pg.viewer.mpl.drawModel1D(ax, values=res*5, depths=np.cumsum(thk),
    ...                          plot='semilogx', color='blue')
    >>> pg.viewer.mpl.drawModel1D(ax, values=res, thickness=thk, z0=1,
    ...                          plot='semilogx', color='red')
    >>> pg.wait()
    """
    if model is not None:
        nLayers = (len(model)-1)//2
        thickness = model[:nLayers]
        values = model[nLayers:]

    if thickness is None and depths is None:
        raise AttributeError("Either thickness or depths must be given.")

    nLayers = len(values)
    px = np.zeros(nLayers * 2)
    pz = np.zeros(nLayers * 2)

    z1 = np.cumsum(thickness) + z0 if thickness is not None else depths

    for i in range(nLayers):
        px[2 * i] = values[i]
        px[2 * i + 1] = values[i]

        if i == nLayers - 1:
            pz[2 * i + 1] = zmax or z1[i - 1] * 1.2
        else:
            pz[2 * i + 1] = z1[i]
            pz[2 * i + 2] = z1[i]

    if plot == 'loglog' or plot == 'semilogy':
        if z0 == 0:
            pz[0] = z1[0] / 2.
        else:
            pz[0] = z0

    try:
        plot = getattr(ax, plot)
        plot(px, pz+z0, **kwargs)
    except BaseException as e:
        print(e)

    ax.set_ylabel(zlabel)
    ax.set_xlabel(xlabel)
    # assume positive depths pointing upward
    if plot in ['plot', 'semilogx']:
        ax.set_ylim(pz[-1], pz[0])

    ax.grid(True)


showModel1D.__doc__ += drawModel1D.__doc__


def draw1DColumn(ax, x, val, thk, width=30, ztopo=0, cMin=1, cMax=1000,
                 cMap=None, name=None, textoffset=0.0, **kwargs):
    """Draw a 1D column (e.g., from a 1D inversion) on a given ax.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from pygimli.viewer.mpl import draw1DColumn
    >>> thk = [1, 2, 3, 4]
    >>> val = thk
    >>> fig, ax = plt.subplots()
    >>> draw1DColumn(ax, 0.5, val, thk, width=0.1, cmin=1, cmax=4, name="VES")
    <matplotlib.collections.PatchCollection object at ...>
    >>> _ = ax.set_ylim(-np.sum(thk), 0)
    """
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import LogNorm

    z = -np.hstack([0., np.cumsum(thk), np.sum(thk) * 1.5]) + ztopo
    recs = []
    for i in range(len(val)):
        recs.append(Rectangle((x - width / 2., z[i]), width, z[i + 1] - z[i]))

    pp = PatchCollection(recs)
    col = ax.add_collection(pp)

    pp.set_edgecolor(None)
    pp.set_linewidth(0.0)

    if cMap is not None:
        if isinstance(cMap, str):
            pp.set_cmap(pg.viewer.mpl.cmapFromName(cMap))
        else:
            pp.set_cmap(cMap)

    if kwargs.pop("logScale", True):
        pp.set_norm(LogNorm(cMin, cMax))

    pp.set_array(np.array(val))
    pp.set_clim(cMin, cMax)
    if name:
        ax.text(x+textoffset, ztopo, name, ha='center', va='bottom')

    updateAxes_(ax)
    return col


def showmymatrix(mat, x, y, dx=2, dy=1, xlab=None, ylab=None, cbar=None):
    """What is this good for?."""
    pg.critical('who use this?')
    plt = pg.plt
    plt.imshow(mat, interpolation='nearest')
    plt.xticks(np.arange(0, len(x), dx), ["%g" % rndig(xi, 2) for xi in x])
    plt.yticks(np.arange(0, len(y), dy), ["%g" % rndig(yi, 2) for yi in y])
    plt.ylim((len(y) - 0.5, -0.5))

    if xlab is not None:
        plt.xlabel(xlab)

    if ylab is not None:
        plt.ylabel(ylab)

    plt.axis('auto')
    if cbar is not None:
        plt.colorbar(orientation=cbar)

    return


def draw1dmodelErr(x, xL, xU=None, thk=None, xcol='g', ycol='r', **kwargs):
    """TODO."""
    pg.critical("in use?")
    if thk is None:
        nlay = (len(x) + 1) / 2
        thk = np.array(x)[:nlay - 1]
        x = np.asarray(x)[nlay - 1:nlay * 2 - 1]
        thkL = np.array(xL)[:nlay - 1]
        thkU = np.array(xU)[:nlay - 1]
        xL = np.asarray(xL)[nlay - 1:nlay * 2 - 1]
        xU = np.asarray(xU)[nlay - 1:nlay * 2 - 1]

#    thk0 = np.hstack((thk, 0.))
#    thkL0 = np.hstack((thkL, 0.))
#    thkU0 = np.hstack((thkU, 0.))
    zm = np.hstack((np.cumsum(thk) - thk / 2, np.sum(thk) * 1.2))  # midpoint
    zc = np.cumsum(thk)  # cumulative
    draw1dmodel(x, thk, **kwargs)
    plt = pg.plt
    plt.xlim(min(xL) * 0.95, max(xU) * 1.05)
    plt.ylim(zm[-1] * 1.1, 0.)
    plt.errorbar(
        x, zm, fmt='.', xerr=np.vstack(
            (x - xL, xU - x)), ecolor=xcol, **kwargs)
    plt.errorbar((x[:-1] + x[1:]) / 2, zc, fmt='.',
                 yerr=np.vstack((thk - thkL, thkU - thk)), ecolor=ycol)


def showStitchedModels(models, ax=None, x=None, cMin=None, cMax=None, thk=None,
                       logScale=True, title=None, zMin=0, zMax=0, zLog=False,
                       **kwargs):
    """Show several 1d block models as (stitched) section.

    Parameters
    ----------
    model : iterable of iterable (np.ndarray or list of np.array)
        1D models (consisting of thicknesses and values) to plot
    ax : matplotlib axes [None - create new]
        axes object to plot in
    x : iterable
        positions of individual models
    cMin/cMax : float [None - autodetection from range]
        minimum and maximum colorscale range
    logScale : bool [True]
        use logarithmic color scaling
    zMin/zMax : float [0 - automatic]
        range for z (y axis) limits
    zLog : bool
        use logarithmic z (y axis) instead of linear
    topo : iterable
        vector of elevation for shifting
    thk : iterable
        vector of layer thicknesses for all models

    Returns
    -------
    ax : matplotlib axes [None - create new]
        axes object to plot in
    """
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import LogNorm

    if x is None:
        x = np.arange(len(models))

    topo = kwargs.pop('topo', np.zeros_like(x))

    fig = None
    if ax is None:
        fig, ax = pg.plt.subplots()

    dxmed2 = kwargs.pop("width", np.median(np.diff(x))) / 2.
    patches = []
    zMinLimit = 9e99
    zMaxLimit = 0

    if thk is not None:
        nlay = len(models[0])
    else:
        nlay = int(np.floor((len(models[0]) + 1) / 2.))

    vals = np.zeros((len(models), nlay))
    for i, imod in enumerate(models):
        if thk is not None:  # take only resistivity from model
            vals[i, :] = imod
            thki = thk
        else:  # extract thickness from model vector
            if isinstance(imod, pg.Vector):
                vals[i, :] = imod[nlay - 1:2 * nlay - 1]
                thki = np.asarray(imod[:nlay - 1])
            else:
                vals[i, :] = imod[nlay - 1:2 * nlay - 1]
                thki = imod[:nlay - 1]

        if zMax > 0:
            z = np.hstack((0., np.cumsum(thki), zMax))
        else:
            thki = np.hstack((thki, thki[-1]*3))
            z = np.hstack((0., np.cumsum(thki)))

        z = topo[i] - z
        zMinLimit = min(zMinLimit, z[-1])
        zMaxLimit = max(zMaxLimit, z[0])

        for j in range(nlay):
            rect = Rectangle((x[i] - dxmed2, z[j]),
                             dxmed2 * 2, z[j+1]-z[j])
            patches.append(rect)

    p = PatchCollection(patches)  # , cmap=cmap, linewidths=0)
    if cMin is not None:
        p.set_clim(cMin, cMax)

    setMappableData(p, vals.ravel(), logScale=logScale)
    ax.add_collection(p)

    if logScale:
        norm = LogNorm(cMin, cMax)
        p.set_norm(norm)

    if 'cMap' in kwargs:
        p.set_cmap(kwargs['cMap'])

#    ax.set_ylim((zMaxLimit, zMin))
    ax.set_ylim((zMinLimit, zMaxLimit))

    if zLog:
        ax.set_yscale("log", nonposy='clip')

    ax.set_xlim((min(x) - dxmed2, max(x) + dxmed2))

    if title is not None:
        ax.set_title(title)

    if kwargs.pop('colorBar', True):
        cb = pg.viewer.mpl.createColorBar(p, cMin=cMin, cMax=cMax, nLevs=5)
#    cb = plt.colorbar(p, orientation='horizontal',aspect=50,pad=0.1)
        if 'cticks' in kwargs:
            xt = np.unique(np.clip(kwargs['cticks'], cMin, cMax))
            cb.set_ticks(xt)
            cb.set_ticklabels([str(xti) for xti in xt])
        if 'label' in kwargs:
            cb.set_label(kwargs['label'])

    pg.plt.draw()
    return ax  # maybe return cb as well?


def draw1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0):
    """DEPRECATED."""
    pg.critical("STYLE_WARNING!!!!!!! don't use this call. "
          "WHO use this anymore?? Use show1dmodel or drawModel1D instead.")
    show1dmodel(x, thk, xlab, zlab, islog, z0)


def show1dmodel(x, thk=None, xlab=None, zlab="z in m", islog=True, z0=0,
                **kwargs):
    """Show 1d block model defined by value and thickness vectors."""
    pg.error("rename after naming convention, don't use plt")
    pg.critical("STYLE_WARNING!!!!!!! don't use this call. "
          "WHO use this anymore??.")

    if xlab is None:
        xlab = "$\\rho$ in $\\Omega$m"

    if thk is None:  # gimli blockmodel (thk+x together) given
        nl = int(np.floor((len(x) - 1) / 2.)) + 1
        thk = np.asarray(x)[:nl - 1]
        x = np.asarray(x)[nl - 1:nl * 2 - 1]

    z1 = np.concatenate(([0], np.cumsum(thk))) + z0
    z = np.concatenate((z1, [z1[-1] * 1.2]))
    nl = len(x)  # x.size()
    px = np.zeros((nl * 2, 1))
    pz = np.zeros((nl * 2, 1))
    for i in range(nl):
        px[2 * i] = x[i]
        px[2 * i + 1] = x[i]
        pz[2 * i + 1] = z[i + 1]
        if i < nl - 1:
            pz[2 * i + 2] = z[i + 1]

    plt = pg.plt
    #    plt.cla()
    if islog:
        plt.semilogx(px, pz, **kwargs)
    else:
        plt.plot(px, pz, **kwargs)

    plt.ion()
    plt.grid(which='both')
    plt.xlim((np.min(x) * 0.9, np.max(x) * 1.1))
    plt.ylim((max(z1) * 1.15, 0.))
    plt.xlabel(xlab)
    plt.ylabel(zlab)
    plt.show()
    return


def showfdemsounding(freq, inphase, quadrat, response=None, npl=2):
    """Show FDEM sounding as real(inphase) and imaginary (quadrature) fields.

    Show FDEM sounding as real(inphase) and imaginary (quadrature) fields
        normalized by the (purely real) free air solution.
    """
    pg.error("rename after naming convention, don't use plt")
    pg.critical("STYLE_WARNING!!!!!!! don't use this call. "
          "WHO use this anymore??.")
    nf = len(freq)
    plt = pg.plt
    fig = plt.figure(1)
    fig.clf()
    ax1 = fig.add_subplot(1, npl, npl - 1)
    plt.semilogy(inphase, freq, 'x-')
    if response is not None:
        plt.semilogy(np.asarray(response)[:nf], freq, 'x-')

    plt.grid(which='both')
    ax2 = fig.add_subplot(1, npl, npl)
    plt.semilogy(quadrat, freq, 'x-')
    if response is not None:
        plt.semilogy(np.asarray(response)[nf:], freq, 'x-')
    plt.grid(which='both')
    fig.show()
    ax = [ax1, ax2]
    if npl > 2:
        ax3 = fig.add_subplot(1, npl, 1)
        ax.append(ax3)

    return ax
