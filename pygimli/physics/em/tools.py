"""Tools for electromagnetics."""
import numpy as np
import pygimli as pg

def cmapDAERO():
    """Standardized colormap from A-AERO projects (purple=0.3 to red=500)."""
    from matplotlib.colors import LinearSegmentedColormap
    CMY = np.array([
        [127, 255, 31], [111, 255, 47], [95, 255, 63], [79, 255, 79],
        [63, 255, 95], [47, 255, 111], [31, 255, 127], [16, 255, 159],
        [0, 223, 159], [0, 191, 159], [0, 159, 207], [0, 127, 175],
        [0, 95, 175], [0, 63, 175], [0, 47, 175], [0, 31, 191], [0, 0, 255],
        [0, 0, 159], [15, 0, 127], [47, 0, 143], [79, 0, 143], [111, 0, 143],
        [143, 0, 127], [159, 31, 63], [175, 47, 31], [207, 63, 0],
        [223, 111, 0], [231, 135, 0], [239, 159, 0], [255, 191, 47],
        [239, 199, 63], [223, 207, 79], [207, 239, 111]], dtype=float)
    RGB = 1.0 - CMY/255
    return LinearSegmentedColormap.from_list('D-AERO', RGB)

def registerDAEROcmap():
    """Standardized colormap from A-AERO projects (purple=0.3 to red=500).

    Example
    -------
    >>> import pygimli as pg
    >>> cmap = pg.physics.em.hemmodelling.registerDAEROcmap()
    >>> mesh = pg.createGrid(20,2)
    >>> data = pg.x(mesh.cellCenters())
    >>> _ = pg.show(mesh, data, cMap=cmap)
    """
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib as mpl

    CMY = np.array([
        [127, 255, 31], [111, 255, 47], [95, 255, 63], [79, 255, 79],
        [63, 255, 95], [47, 255, 111], [31, 255, 127], [16, 255, 159],
        [0, 223, 159], [0, 191, 159], [0, 159, 207], [0, 127, 175],
        [0, 95, 175], [0, 63, 175], [0, 47, 175], [0, 31, 191], [0, 0, 255],
        [0, 0, 159], [15, 0, 127], [47, 0, 143], [79, 0, 143], [111, 0, 143],
        [143, 0, 127], [159, 31, 63], [175, 47, 31], [207, 63, 0],
        [223, 111, 0], [231, 135, 0], [239, 159, 0], [255, 191, 47],
        [239, 199, 63], [223, 207, 79], [207, 239, 111]], dtype=float)
    RGB = 1.0 - CMY/255
    daero = LinearSegmentedColormap.from_list('D-AERO', RGB)
    mpl.colormaps.register(name='daero', cmap=daero)
    return daero

def xfplot(ax, DATA, x, freq, everyx=5, orientation='horizontal', aspect=30,
           label=None, cMap="Spectral_r"):
    """Plot a matrix according to x and frequencies."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    nt = list(range(0, len(x), everyx))
    im = ax.matshow(DATA.T, interpolation='nearest', cmap=cMap)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_xticks(nt)
    ax.set_xticklabels([f"{xi:g}" for xi in x[nt]])
    ax.set_yticks(list(range(0, len(freq) + 1, 2)))
    ax.set_yticklabels([f"{freq[i]:g}" for i in range(0, len(freq), 2)])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('f [Hz]')
    ax.xaxis.set_label_position('top')
    ax.set_aspect("auto")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='5%', pad=0.3)
    pg.plt.colorbar(im, ax=ax, cax=cax, orientation=orientation, aspect=aspect)
    if label is not None:
        cax.set_title(label)
    # plt.colorbar(im, ax=ax, orientation=orientation, aspect=aspect)
    return im

class FDEMsystems:
    """FDEM system frequency/distance definitions."""

    def __init__(self):
        """Initialize FDEM system definitions."""
        self.frequencies = {
            'MaxMin10': 2**np.arange(10) * 110.,
            'MaxMin8': 2**np.arange(8) * 110.,
            'ResolveHCP': np.array([387., 1820., 8330., 41500., 133400.]),
            'ResolveVCX': 5410.,
            'ResolveHCPOld': np.array([380., 1770., 8300., 41000., 129500.]),
            'BKS36a': np.array([386, 1817, 8360, 41420, 133200, 5390]),
            'BKS60': np.array([380, 1773, 8300, 41000, 129500, 5410])
        }
        self.distances = {
            'MaxMin10': np.array([3.66, 4.57, 5.49, 6.40, 7.32, 8.23, 9.14,
                                  10.06, 10.97, 11.89]),
            'MaxMin8': np.array([3.66, 5.49, 7.32, 9.14, 10.97,
                                 12.80, 14.63, 16.46]),
            'ResolveHCP': np.array([7.938, 7.931, 7.925, 7.912, 7.918]),
            'ResolveVCX': 9.055,
            'ResolveHCPOld': np.array([7.918, 7.918, 7.957, 8.033, 7.906]),
            'BKS36a': np.array([7.938, 7.931, 7.925, 7.912, 7.918, 9.055]),
            'BKS60': np.array([7.918, 7.918, 7.957, 8.033, 7.906, 9.042])
        }

    def getFrequenciesDistances(self, system):
        """Get frequencies for a given system."""
        return self.frequencies.get(system, None), \
            self.distances.get(system, None)
