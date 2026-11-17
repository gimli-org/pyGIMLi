"""Matplotlib style parameters used for the pyGIMLi documentation."""

import math

phi = (math.sqrt(5) + 1) / 2
fontscale = 1.0
dark_gray = ".15"
light_gray = ".8"

plot_rcparams = {
    "font.size": 12 * fontscale,
    "axes.labelsize": 11 * fontscale,
    "axes.titlesize": 12 * fontscale,
    "xtick.labelsize": 10 * fontscale,
    "ytick.labelsize": 10 * fontscale,
    "legend.fontsize": 10 * fontscale,
    "grid.linewidth": 1,
    "lines.linewidth": 1.75,
    "patch.linewidth": 0.3,
    "lines.markersize": 7,
    "lines.markeredgewidth": 0,
    "xtick.major.width": 1,
    "ytick.major.width": 1,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.pad": 7,
    "ytick.major.pad": 7,
    "text.usetex": False,
    "figure.subplot.bottom": 0.2,
    "figure.subplot.left": 0.2,
    "figure.subplot.right": 0.9,
    "figure.subplot.top": 0.85,
    "figure.subplot.wspace": 0.4,
    "figure.figsize": (5 * phi, 5),
    "figure.dpi": 96,
    "figure.facecolor": "white",
    "text.color": dark_gray,
    "axes.labelcolor": dark_gray,
    "legend.frameon": False,
    "legend.numpoints": 1,
    "legend.scatterpoints": 1,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.color": dark_gray,
    "ytick.color": dark_gray,
    "axes.axisbelow": True,
    "image.cmap": "viridis",
    "font.family": ["sans-serif"],
    "font.sans-serif": [
        "Arial",
        "DejaVu Sans",
        "Liberation Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ],
    "grid.linestyle": "-",
    "lines.solid_capstyle": "round",
    "axes.facecolor": "white",
    "axes.edgecolor": dark_gray,
    "axes.linewidth": 1.25,
    "grid.color": light_gray,
}
