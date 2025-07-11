#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions to convert pg SparseMatrices from and to numpy objects"""

import numpy as np
import pygimli as pg


def toSparseMatrix(A):
    """Convert any matrix type to pg.SparseMatrix and return copy of it.

    No conversion if A is a SparseMatrix already
    Arguments
    ---------
    A: pg or scipy matrix

    Returns
    -------
        pg.SparseMatrix
    """
    if isinstance(A, pg.matrix.BlockMatrix):
        return pg.matrix.SparseMatrix(A.sparseMapMatrix())

    if isinstance(A, pg.matrix.CSparseMapMatrix):
        return pg.matrix.CSparseMatrix(A)

    if isinstance(A, pg.matrix.CSparseMatrix):
        return A

    if isinstance(A, pg.matrix.SparseMatrix):
        return A
        #return pg.matrix.SparseMatrix(A)

    if isinstance(A, pg.matrix.SparseMapMatrix):
        return pg.matrix.SparseMatrix(A)

    from scipy.sparse import csr_matrix

    if isinstance(A, csr_matrix):
        return pg.SparseMatrix(A.indptr, A.indices, A.data)

    from scipy.sparse import coo_matrix
    if isinstance(A, coo_matrix):
        pg.critical('implement me')
        return pg.matrix.SparseMatrix(A)

    return toSparseMatrix(csr_matrix(A))


def toSparseMapMatrix(A):
    """Convert any matrix type to pg.SparseMatrix and return copy of it.

    Arguments
    ---------
    A: pg or scipy matrix

    Returns
    -------
        pg.SparseMatrix
    """
    if isinstance(A, pg.matrix.SparseMapMatrix):
        return pg.matrix.SparseMapMatrix(A)

    if isinstance(A, pg.matrix.BlockMatrix):
        return A.sparseMapMatrix()

    if isinstance(A, pg.matrix.SparseMatrix):
        return pg.matrix.SparseMapMatrix(A)

    from scipy.sparse import csr_matrix
    if isinstance(A, csr_matrix):
        pg.critical('implement me')

        return pg.matrix.SparseMapMatrix(A)

    from scipy.sparse import coo_matrix
    if isinstance(A, coo_matrix):
        pg.critical('implement me')
        return pg.matrix.SparseMapMatrix(A)

    return toSparseMapMatrix(csr_matrix(A))

def toCSR(A):
    return sparseMatrix2csr(A)

def toCOO(A):
    return sparseMatrix2coo(A)

def sparseMatrix2csr(A):
    """Convert SparseMatrix to scipy.csr_matrix.

    Compressed Sparse Row matrix, i.e., Compressed Row Storage (CRS)

    Parameters
    ----------
    A: pg.matrix.SparseMapMatrix | pg.matrix.SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.csr_matrix
        Matrix to convert into.
    """
    #optImport(scipy.sparse, requiredFor="toCRS_matrix")
    shape=(A.rows(), A.cols())

    from scipy.sparse import csr_matrix
    if isinstance(A, pg.matrix.CSparseMapMatrix):
        C = pg.matrix.CSparseMatrix(A)
        return csr_matrix((C.vecVals().array(),
                           C.vecRowIdx(),
                           C.vecColPtr()), shape=shape, dtype=complex)
    if isinstance(A, pg.matrix.SparseMapMatrix):
        C = pg.matrix.SparseMatrix(A)
        return csr_matrix((C.vecVals().array(),
                           C.vecRowIdx(),
                           C.vecColPtr()), shape=shape)
    if isinstance(A, pg.matrix.SparseMatrix):
        return csr_matrix((A.vecVals().array(),
                           A.vecRowIdx(),
                           A.vecColPtr()), shape=shape)
    if isinstance(A, pg.matrix.CSparseMatrix):
        csr = csr_matrix((A.vecVals().array(),
                           A.vecRowIdx(),
                           A.vecColPtr()), shape=shape, dtype=complex)
        return csr
    if isinstance(A, pg.matrix.BlockMatrix):
        M = A.sparseMapMatrix()
        return sparseMatrix2csr(M)

    return csr_matrix(A)


def sparseMatrix2coo(A, rowOffset=0, colOffset=0):
    """Convert SparseMatrix to scipy.coo_matrix.

    Parameters
    ----------
    A: pg.matrix.SparseMapMatrix | pg.matrix.SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.coo_matrix
        Matrix to convert into.
    """
    from scipy.sparse import coo_matrix

    vals = pg.Vector()
    rows = pg.core.IndexArray([0])
    cols = pg.core.IndexArray([0])

    if isinstance(A, pg.matrix.SparseMatrix):
        C = pg.matrix.SparseMapMatrix(A)
        C.fillArrays(vals=vals, rows=rows, cols=cols)
        rows += rowOffset
        cols += colOffset
        return coo_matrix((vals, (rows, cols)), shape=(A.rows(), A.cols()))

    elif isinstance(A, pg.matrix.SparseMapMatrix):
        A.fillArrays(vals, rows, cols)
        rows += rowOffset
        cols += colOffset

        return coo_matrix((vals, (rows, cols)), shape=(A.rows(), A.cols()))

    return coo_matrix(A)


def convertCRSIndex2Map(rowIdx, colPtr):
    """Converts CRS indices to uncompressed indices (row, col)."""
    ii = []
    jj = []
    for i in range(len(colPtr)-1):
        for j in range(colPtr[i], colPtr[i + 1]):
            ii.append(i)
            jj.append(rowIdx[j])
    return ii, jj


def sparseMatrix2Array(matrix, indices=True, getInCRS=True):
    """Extract indices and value from sparse matrix (SparseMap or CRS)

    Get python Arrays from SparseMatrix or SparseMapMatrix in either CRS
    convention (row index, column Start_End, values) or row index, column
    index, values.

    Parameters
    ----------

    matrix: pg.matrix.SparseMapMatrix or pg.matrix.SparseMatrix
        Input matrix to be transformed to numpy arrays.

    indices: boolean (True)
        Decides weather the indices of the matrix will be returned or not.

    getInCSR: boolean (True)
        If returned, the indices can have the format of a compressed row
        storage (CSR), the default or uncompressed lists with column and row
        indices.

    Returns
    -------

    vals: numpy.ndarray
        Entries of the matrix.

    indices: list, list
        Optional. Returns additional array with the indices for reconstructing
        the matrix in the defined format.

    """
    io_warn = 'Only working for pygimli.SparseMapMatrix or CSR shaped ' +\
              'pygimli.Sparse matrices. Import type is {}'
    assert isinstance(matrix, pg.matrix.SparseMapMatrix) or\
        isinstance(matrix, pg.matrix.SparseMatrix), io_warn.format(type(matrix))

    if not isinstance(matrix, pg.matrix.SparseMatrix):
        matrix = pg.matrix.SparseMatrix(matrix)

    vals = np.array(matrix.vecVals())
    if indices is True:
        rows = list(matrix.vecRowIdx())
        cols = list(matrix.vecColPtr())
        if getInCRS:
            return list(rows), list(cols), vals
        else:
            rr, cc = convertCRSIndex2Map(rows, cols)
            return rr, cc, vals
    else:
        return vals

def sparseMatrix2Dense(matrix):
    """Convert sparse matrix to dense ndarray"""

    rr, cc, vals = sparseMatrix2Array(matrix, indices=True, getInCRS=False)
    mat = np.zeros((matrix.rows(), matrix.cols()))

    for i, r in enumerate(rr):
        mat[r, cc[i]] = vals[i]

    return mat

if __name__ == '__main__':
    pass
# The End
