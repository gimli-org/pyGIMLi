#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Functions to draw various pygimli matrices with matplotlib."""

import numpy as np
# try to avoid plt impoert #
# import matplotlib.pyplot as plt

import pygimli as pg
from .colorbar import cmapFromName


def drawSparseMatrix(ax, mat, **kwargs):
    """Draw a view of a matrix into the axes.

    Parameters
    ----------
    ax : mpl axis instance, optional
        Axis instance where the matrix will be plotted.

    mat: pg.matrix.SparseMatrix or pg.matrix.SparseMapMatrix

    Returns
    -------
    mpl.lines.line2d

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawSparseMatrix
    >>> A = pg.randn((10,10), seed=0)
    >>> SM = pg.matrix.SparseMapMatrix()
    >>> for i in range(10):
    ...     SM.setVal(i, i, 5.0)
    >>> fig, (ax1, ax2) = pg.plt.subplots(1, 2, sharey=True, sharex=True)
    >>> _ = drawSparseMatrix(ax1, A, colOffset=5, rowOffset=5, color='blue')
    >>> _ = drawSparseMatrix(ax2, SM, color='green')
    """
    row = kwargs.pop('rowOffset', 0)
    col = kwargs.pop('colOffset', 0)
    color = kwargs.pop('color', None)

    mat = pg.utils.sparseMatrix2coo(mat)
    mat.row += row
    mat.col += col
    gci = ax.spy(mat, color=color, **kwargs)

    ax.autoscale(enable=True, axis='both', tight=True)
    return gci

def drawBlockMatrix(ax, mat, **kwargs):
    """Draw a view of a matrix into the axes.

    Arguments
    ---------

    ax : mpl axis instance, optional
        Axis instance where the matrix will be plotted.

    mat: pg.Matrix.BlockMatrix

    Keyword Arguments
    -----------------
    spy: bool [False]
        Draw all matrix entries instead of colored blocks

    Returns
    -------
    ax:

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> I = pg.matrix.IdentityMatrix(10)
    >>> SM = pg.matrix.SparseMapMatrix()
    >>> for i in range(10):
    ...     SM.setVal(i, 10 - i, 5.0)
    ...     SM.setVal(i,  i, 5.0)
    >>> B = pg.matrix.BlockMatrix()
    >>> B.add(I, 0, 0)
    0
    >>> B.add(SM, 10, 10)
    1
    >>> print(B)
    pg.matrix.BlockMatrix of size 20 x 21 consisting of 2 submatrices.
    >>> fig, (ax1, ax2) = pg.plt.subplots(1, 2, sharey=True)
    >>> _ = pg.show(B, ax=ax1)
    >>> _ = pg.show(B, spy=True, ax=ax2)
    """
    if kwargs.pop('spy', False):
        gci = []
        ids = pg.unique([e.matrixID for e in mat.entries()])
        cMap = cmapFromName("Set3", ncols=len(ids))

        for e in mat.entries():
            mid = e.matrixID

            mati = mat.mat(mid)
            if isinstance(mati, pg.matrix.IdentityMatrix):
                mati = np.eye(mati.size())

            gci.append(drawSparseMatrix(ax, mati,
                                        rowOffset=e.rowStart,
                                        colOffset=e.colStart,
                                        color=cMap(mid)))

        return gci, None
    else:
        plcs = []
        for e in mat.entries():
            mid = e.matrixID
            widthy = mat.mat(mid).rows() - 0.1 # to make sure non-matrix regions are not connected in the plot
            widthx = mat.mat(mid).cols() - 0.1
            plc = pg.meshtools.createRectangle(
                [e.colStart, e.rowStart],
                [e.colStart + widthx, e.rowStart + widthy],
                marker=mid)
            plcs.append(plc)

        bm = pg.meshtools.mergePLC(plcs)
        gci, cBar = pg.viewer.mpl.drawPLC(ax, bm, fitView=False, **kwargs)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        if cBar is not None:
            cBar.set_label("Matrix ID")

        if len(mat.entries()) > 10:
            gci.set_cmap("viridis")

        return gci, cBar
