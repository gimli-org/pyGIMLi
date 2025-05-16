# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 09:33:14 2015

@author: Marcus Wennermark
"""
import numpy as np

from .modelview import draw1DColumn
from .colorbar import cmapFromName


def create_legend(ax, cmap, ids, classes):
    """
    Create a list of patch objects that can be used for borehole legends.
    """
    import matplotlib.patches as mpatches

    patches = [mpatches.Patch(color=cmap(i), label=classes[i]) for i in ids]
    return patches


class BoreHole(object):
    r"""Class for handling (structural) borehole data for inclusion in plots.

    Each row in the data file must contain a start depth [m], end depth and a
    classification. The values should be separated by whitespace.
    The first line should contain the inline position (x and z), text ID and an
    offset for the text (for plotting).

    The format is very simple, and a sample file could look like this:

        16.0 0.0 BOREHOLEID_TEXTOFFSET \n
        0.0 1.6 clti \n
        1.6 10.0 shale
    """
    def __init__(self, fname):
        self._fname = fname
        self._load()

    def __repr__1(self):
        return self.__class__.__name__ + '("{}")'.format(self._fname)

    def __repr__(self):
        out = 'Borehole id: {}\n Inline position (x, z): {}\n Layers:'.format(
            self.borehole_id, self.inline_pos)

        for layer in self.data[1:]:
            out = ''.join([out, '\n ', str(layer)])

        return out

    def _load(self):
        """Loads the data file."""
        self.data = np.genfromtxt(self._fname, dtype=None)
        if len(self.data) > 1:
            header = self.data[0][2].decode('UTF-8').split('_')
            self.borehole_id = header[0]
            if len(header) > 1:
                self._textoffset = float(header[1])

            self.inline_pos = (self.data[0][0], self.data[0][1])
            self.classes = [d[-1].decode('UTF-8') for d in self.data[1:]]
            self.unique_classes, self.class_id = \
                np.unique(self.classes, return_inverse=True)
        else:
            raise Warning('File "{}" contains no layers!'.format(self._fname))

    def plot(self, ax, plot_thickness=1.0, cmin=None, cmax=None, cm=None,
             do_legend=True, **legend_kwargs):
        """Plots the data on the specified axis."""
        start_depths = np.asarray([d[0] for d in self.data[1:]])
        end_depths = np.asarray([d[1] for d in self.data[1:]])
        thickness = end_depths - start_depths

        if cmin is None or cmax is None:
            cmin = min(self.class_id)
            cmax = max(self.class_id)

        if cm is None:
            cm = cmapFromName("Set3", len(self.unique_classes))
            # cm = plt.get_cmap('jet', len(self.unique_classes))

        draw1DColumn(ax, self.inline_pos[0], self.class_id, thickness,
                     ztopo=self.inline_pos[1], width=plot_thickness, cMin=cmin,
                     cMax=cmax, name=self.borehole_id, cMap=cm,
                     textoffset=self._textoffset, logScale=False)

        if do_legend:
            self.add_legend(ax, cm, **legend_kwargs)

    def add_legend(self, ax, cmap, **legend_kwargs):
        """Add a legend to the plot."""

        leg = create_legend(ax, cmap, self.class_id, self.unique_classes)
        ax.legend(handles=leg, **legend_kwargs)


class BoreHoles(object):
    """Class to load and handle several boreholes belonging to one profile.

    It makes the color coding consistent across boreholes.
    """

    def __init__(self, fnames):
        """Load a list of bore hole from filenames."""
        if isinstance(fnames, str):
            if fnames.find("*") >= 0:
                from glob import glob
                fnames = glob(fnames)
        self._fnames = fnames
        if len(fnames) > 0:
            self.boreholes = [BoreHole(f) for f in fnames]
        else:
            raise Warning('No filenames specified!')

    def __repr__1(self):
        return self.__class__.__name__ + '({})'.format(self._fnames)

    def __repr__(self):
        out = '{} borehole files loaded:'.format(len(self._fnames))

        for b in self.boreholes:
            out = ''.join([out, str(b)])
        return out

    def _build_common_colormap(self):
        """Create a common colormap for all boreholes.

        Such that a certain classification has the same color on all boreholes.
        """
        import matplotlib.pyplot as plt
        self.common_unique, rev_idx = np.unique(
            np.hstack([b.classes for b in self.boreholes]),
            return_inverse=True)

        self.class_id = rev_idx

        start_idx = 0
        for b in self.boreholes:
            diff_idx = len(b.class_id)
            end_idx = start_idx + diff_idx
            b.class_id = self.class_id[start_idx:end_idx]
            b.unique_classes = self.common_unique
            start_idx += diff_idx

        self.cm = plt.get_cmap('jet', len(self.common_unique))
        self.cmin = min(self.class_id)
        self.cmax = max(self.class_id)

    def plot(self, ax, plot_thickness=1.0, do_legend=True, **kwargs):
        """Plot the boreholes on the specified axis.
        
        Parameters
        ----------
        ax : matplotlib.Axes
        plot_thickness : float [1.0]
            width (in m) for borehole columns
        do_legend : bool
            draw a legend for the geological units
        **kwargs : keyword arguments passed to legend
            fontsize, markerscale, (s. matplotlib.legend)
        """
        self._build_common_colormap()

        for b in self.boreholes:
            b.plot(ax, plot_thickness=plot_thickness, do_legend=False,
                   cmin=self.cmin, cmax=self.cmax, cm=self.cm)

        if do_legend:
            self.add_legend(ax, self.cm, **kwargs)

    def add_legend(self, ax, cmap, **legend_kwargs):
        """Add a legend to the plot.
        
        Parameters
        ----------
        ax : matplotlib.axes
        cmap : matplotlib.colormap
        
        Keyword arguments passed to ax.legend:
        fontsize, markerscale
        """
        leg = create_legend(ax, cmap, np.arange(cmap.N), self.common_unique)

        extra = dict(bbox_to_anchor=(0.9, 0.05, 0.1, 0.1), ncol=cmap.N / 2,
                     mode='scale', fancybox=False, shadow=False, fontsize=3.0,
                     columnspacing=1.0, loc='center right', markerscale=0.7,
                     framealpha=1.0, borderpad=0.5, handleheight=0.5,
                     frameon=True)
        
        extra.update(legend_kwargs)
        ax.legend(handles=leg, **extra)

