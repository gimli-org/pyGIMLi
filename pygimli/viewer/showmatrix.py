# -*- coding: utf-8 -*-
"""Generic matrix visualization tools."""

# import matplotlib as mpl
import numpy as np

import pygimli as pg

from .mpl import createColorBar  # , updateColorBar
from .mpl.matrixview import drawBlockMatrix, drawSparseMatrix
from .mpl.colorbar import cmapFromName


def showMatrix(mat, ax=None, **kwargs):
    """Show various pyGIMLi matrices using matplotlib.

    Args
    ----
    mat: matrix

    ax: mpl.axes

    Keyword Args
    ------------
        **kwargs : forwarded to mpl plotting commands

    Returns
    -------
        mpl.axes, Colorbar
    """
    cBar = None
    if ax is None:
        ax = pg.show()[0]

    try:
        from scipy.sparse import spmatrix

        if isinstance(mat, spmatrix):
            kwargs.setdefault('markersize', 1)
            gci = drawSparseMatrix(ax, mat, **kwargs)
            return ax, None
    except ImportError:
        pass

    if isinstance(mat, (pg.core.RSparseMapMatrix, pg.core.RSparseMatrix)):
        gci = drawSparseMatrix(ax, mat, **kwargs)
    elif isinstance(mat, pg.matrix.BlockMatrix):
        gci, cBar = drawBlockMatrix(ax, mat, **kwargs)

        if cBar is None:
            uniqueIDs = pg.unique([e.matrixID for e in mat.entries()])
            cMap = cmapFromName("Set3", ncols=len(uniqueIDs))
            sm = pg.plt.cm.ScalarMappable(cmap=cMap)

            cBar = createColorBar(sm, ax=ax, label="Matrix ID",
                                  cMin=-0.5, cMax=len(uniqueIDs)-0.5)
            ticks = np.arange(len(uniqueIDs))
            cBar.set_ticks(ticks)

            labels = []
            for ID in uniqueIDs:
                label = "{:d}".format(ID)
                labels.append(label)
            cBar.set_ticklabels(labels)
    elif isinstance(mat, pg.matrix.Matrix):
        from pygimli.utils import gmat2numpy
        gci = ax.matshow(gmat2numpy(mat))
        cBar = ax.figure.colorbar(gci)
    elif isinstance(mat, pg.matrix.RealNumpyMatrix):
        gci = ax.matshow(mat.M)
        cBar = ax.figure.colorbar(gci)
    elif isinstance(mat, (pg.matrix.IdentityMatrix,
                          pg.matrix.DiagonalMatrix)):
        x = [0, mat.cols()]
        ax.plot(x, x)
        ax.set_ylim(x[::-1])
    else:
        nC = mat.cols()
        nR = mat.rows()
        ax.fill([0, nC, nC, 0, 0], [0, 0, nR, nR, 0], hatch='/')
        ax.set_ylim(nR, 0)

    return ax, cBar
