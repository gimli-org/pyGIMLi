# -*- coding: utf-8 -*-
"""Some matrix specialization."""

import time

import numpy as np
import pygimli as pg
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix

from .core import pgcore

from .core import (CMatrix, CSparseMapMatrix, CSparseMatrix,
                   RSparseMapMatrix, RSparseMatrix, ElementMatrix,
                   IVector, MatrixBase, R3Vector, RVector)

from .logger import critical, warn, error
from .base import isArray


def warn(*a):
    import pygimli as pg
    pg.warn(*a)

# make core matrices (now in pgcore, later pg.core) available here for brevity
# useful aliases
IdentityMatrix = pgcore.IdentityMatrix

IndexArray = pgcore.IndexArray
Vector = pgcore.RVector
Matrix = pgcore.RMatrix
SparseMatrix = pgcore.RSparseMatrix
SparseMapMatrix = pgcore.RSparseMapMatrix
BlockMatrix = pgcore.RBlockMatrix

# General Monkeypatch core classes
__Matrices = [pgcore.MatrixBase,
              #pgcore.RSparseMatrix,
              #pgcore.RSparseMapMatrix,
              #pgcore.CSparseMatrix,
              #pgcore.CSparseMapMatrix,
              #pgcore.RBlockMatrix,
              #pgcore.IdentityMatrix
              ]

def __Matrix_array_ufunc__(self, ufunc, method, *inputs, **kwargs):
    """Handle ufunc to communicate with numpy."""
    if ufunc == np.multiply:
        ## for np.float * self
        if id(self) == id(inputs[0]) and len(inputs) == 2:
            return self.__mul__(inputs[1])

    if ufunc == np.matmul:
        ## for np.mat * self
        if id(self) == id(inputs[0]) and len(inputs) == 2:
            return self.__mul__(inputs[1])
        # if isinstance(inputs[0], np.ndarray):
        #     return pg.core.mult(self, np.squeeze(inputs[0]))

    pg._r('self:', self)
    pg._r('self:', id(self))
    pg._r('ufunc:', ufunc)
    pg._r('method:', method)
    pg._r('inputs:', len(inputs))
    for i, ii in enumerate(inputs):
        pg._r(f'\t input {i}:', type(ii))
        pg._r(f'\t input {i}:', id(ii))
        pg._r(f'\t input {i}:', ii)

    pg._r('kwargs', kwargs)
    pg.critical('implementme.')


for m in __Matrices:
    m.ndim = 2
    m.__len__ = lambda self: self.rows()
    m.shape = property(lambda self: (self.rows(), self.cols()))
    ## To allow for ignoring np.*.__mul__ in case A has the __rmul__ function
    ## see test TestMatrixMethods.testSparseMatrixBasics, TestRVectorMethods.testUFunc
    m.__array_ufunc__ = __Matrix_array_ufunc__
    m.__rmul__ = lambda self: critical('not yet implemented')


pgcore.RMatrix.dtype = float
pgcore.CMatrix.dtype = complex#
pgcore.RDenseMatrix.dtype = float
pgcore.RSparseMapMatrix.dtype = float
pgcore.CSparseMapMatrix.dtype = complex
pgcore.RSparseMatrix.dtype = float
pgcore.CSparseMatrix.dtype = complex
pgcore.RBlockMatrix.dtype = float
pgcore.CBlockMatrix.dtype = complex


def __RMatrix_str(self):
    import pygimli as pg
    s = str(type(self)).split('.')[-1].split("'")[0] + ": " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 10:
        s += '\n'

        for row in self:
            for v in row:
                s += pg.pf(v).rjust(8)
            s+='\n'
    return s

def __CMatrix_str(self):
    s = "CMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 10:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__() + '\n'
    return s


def __ElementMatrix_str(self):
    """Show entries of an ElementMatrix."""
    import pygimli as pg

    self.integrate()

    if self.mat().cols() == 0 and self.mat().rows() == 0:
        return 'Empty ElementMatrix\n'

    maxRowID = int(np.log10(max(self.rowIDs())))+2

    sAdd = ''
    if hasattr(self, 'isIdentity'):
        sAdd +=f'I={self.isIdentity}'

    if hasattr(self, 'mulR') and self.mulR is not None:
        if (pg.isScalar(self.mulR) and self.mulR != 1.0) or \
            pg.isPos(self.mulR) or \
            pg.isArray(self.mulR) or\
            pg.isMatrix(self.mulR):
            s = '\n ' + f'multR = {self.mulR} (applied)' + '\n ' + ' ' * maxRowID
        elif pg.isScalar(self.mulR, 1.0):
            s = '\n ' + ' ' * maxRowID
        else:
            s = '\n ' + f'multR = {self.mulR} (unknown how to apply)' + '\n ' + ' ' * maxRowID
    else:
        self.mulR = 1
        s = '\n ' + ' ' * maxRowID

    # print(self.mat())
    # print(self.colIDs())
    # print(self.rowIDs())
    # if sAdd != '':
    #     s += f' {sAdd}\n'

    for i in range(self.mat().cols()):
        s += str(self.colIDs()[i]).rjust(9)
    s += '    ' + sAdd
    s += '\n'

    s += '  ' + '-' * self.mat().cols()*(9 + maxRowID) + '-\n'

    for i in range(self.mat().rows()):
        s += str(self.rowIDs()[i]).rjust(maxRowID) + " :"
        if isinstance(self.mulR, (int, float)):
            for v in self.row_RM(i):
                s += pg.pf(v*self.mulR).rjust(9)
        elif pg.isPos(self.mulR):
            if self.mat().cols() <= len(self.mulR):
                for v in self.row_RM(i) * self.mulR[0:self.mat().cols()]:
                    s += pg.pf(v).rjust(9)
            else:
                print(self.mat())
                print(self.row_RM(i))
                print(self.mulR)
                pg.critical('invalid element multR.')
        elif pg.isArray(self.mulR, self.mat().cols()):
            for v in self.row_RM(i)*self.mulR:
                s += pg.pf(v).rjust(9)

        elif isinstance(self.mulR, pg.solver.utils.ConstitutiveMatrix) \
            or self.mulR.shape[0] == self.mat().shape[1]:
            ## or ElasticityMatrix, or any other matrix with same cols as self.mat()
            for v in self.row_RM(i)@self.mulR:
                s += pg.pf(v).rjust(9)

        elif pg.isMatrix(self.mulR):
            if pg.isArray(self.mulR.flatten(), self.mat().cols()):
                for v in self.row(i)*self.mulR.flatten():
                    s += pg.pf(v).rjust(9)
            else:
                try:
                    for v in self.row(i):
                        s += pg.pf(v).rjust(9)
                        s += '\n'
                    s += f'mulR = {self.mulR} (unknown how to apply)\n'
                except Exception as e:
                    print(f"Error processing row {i}: {e}")
                # print('mat:', self.mat())

                # pg._r(self.mat())
                # pg._r(self.rows())
                # print(self.mulR)
                # pg.critical('invalid matrix element mulR.')

        elif pg.isArray(self.mulR):
            # per node vals
            # print(i, self.row_RM(i))
            # print(self.rowIDs())
            for v in (self.row_RM(i) * self.mulR[self.rowIDs()[i]]):
                s += pg.pf(v).rjust(9)

        elif self.mulR is not None:
            print('mat:', self.mat())
            print('multR:', self.mulR)
            try:
                print(mulE(self, f=self.mulR))
            except BaseException:
                pass
            warn('invalid element mulR. should be evaluated with mulE?')
            return ''


        else:
            for v in self.row(i):
                s += pg.pf(v).rjust(9)
        s += '\n'
    return s


pgcore.RMatrix.__repr__ = __RMatrix_str
pgcore.RDenseMatrix.__repr__ = __RMatrix_str
pgcore.CMatrix.__repr__ = __CMatrix_str
pgcore.ElementMatrix.__repr__ = __ElementMatrix_str


def __SparseMatrix_str(self):
    """Show entries of an ElementMatrix."""
    import pygimli as pg

    S = pg.matrix.asSparseMapMatrix(self)
    if S.cols() == 0 or S.rows() == 0:
        return 'Empty ElementMatrix\n'

    name = str(type(self)).split('.')[-1].split("'")[0]
    if name.startswith('R'):
        name = name[1:]
    s = f"{name}: {S.rows()} x {S.cols()}"
    pc = int(self.nVals()/self.cols()/self.rows()*1000) / 10
    s += " (nnz=" + str(self.nVals()) + " / " + str(pc)+ "%)\n"

    if S.cols() < 25:
        s += '\n'

        M = pg.matrix.asDense(self)
        for mi in M:
            for v in mi:
                if (abs(v) < 1e-12 and abs(v) > 0):
                    s += ('+'+pg.pf(v)).rjust(9)
                elif (abs(v)< 1e-12 and abs(v) < 0):
                    s += ('-'+pg.pf(v)).rjust(9)
                else:
                    s += pg.pf(v).rjust(9)
            s += '\n'
    return s
pgcore.RSparseMatrix.__repr__ =__SparseMatrix_str
pgcore.RSparseMapMatrix.__repr__ =__SparseMatrix_str


def __RVector_format(self, f):
    print(f)
    return str(self)
pgcore.RVector.__format__ = __RVector_format


# for copy by reference objects we need to keep the owner objects alive, i.e. increase the reference counter
__ElementMatrix_mat = pgcore.ElementMatrix.mat
def __ElementMatrix_mat_GC__(self):
    R = __ElementMatrix_mat(self)
    R.owner = self
    return R
pgcore.ElementMatrix.mat = __ElementMatrix_mat_GC__


# temporary object for debugging
# __TElementMatrix_mat = pgcore.TestEM.mat
# def __TElementMatrix_mat_GC__(self):
#     R = __TElementMatrix_mat(self)
#     R.owner = self
#     return R
# pgcore.TestEM.mat = __TElementMatrix_mat_GC__


__DenseMatrix_row = pgcore.RDenseMatrix.row
def __DenseMatrix_row_GC__(self, i):
    R = __DenseMatrix_row(self, i)
    R.owner = self
    return R
pgcore.RDenseMatrix.row = __DenseMatrix_row_GC__


def __CopyRMatrixTranspose__(self):
    return pgcore.RMatrix(np.array(self).T)
pgcore.RMatrix.T = property(__CopyRMatrixTranspose__)


def __CopyRDenseMatrixTranspose__(self):
    return pgcore.RDenseMatrix(np.array(self).T)
pgcore.RDenseMatrix.T = property(__CopyRDenseMatrixTranspose__)


## Overload some operators

__RSparseMapMatrix_Mul__orig = pgcore.RSparseMapMatrix.__mul__
def __RSparseMapMatrix_Mul__(self, b):
    """Multiply SparseMapMatrix with b."""
    if isinstance(b, (int, np.float64)):
        ## or this will be wrongly parsed into mul(self, RVector(b))
        return __RSparseMapMatrix_Mul__orig(self, float(b))
    return __RSparseMapMatrix_Mul__orig(self, b)
pgcore.RSparseMapMatrix.__mul__ = __RSparseMapMatrix_Mul__


__RSparseMatrix_Mul__orig = pgcore.RSparseMatrix.__mul__
def __RSparseMatrix_Mul__(self, b):
    """Multiply SparseMatrix with b."""
    # pg._b(b)
    if isinstance(b, (int, np.float64)):
        ## or this will be wrongly parsed into mul(self, RVector(b))
        return __RSparseMatrix_Mul__orig(self, float(b))

    if isinstance(b, np.ndarray):
        return __RSparseMatrix_Mul__orig(self, np.squeeze(b))

    return __RSparseMatrix_Mul__orig(self, b)
pgcore.RSparseMatrix.__mul__ = __RSparseMatrix_Mul__


def __RSparseMapMatrix_Neg__(self):
    """Multiply SparseMapMatrix with b."""
    return self * -1.0
pgcore.RSparseMapMatrix.__neg__ = __RSparseMapMatrix_Neg__
pgcore.RSparseMatrix.__neg__ = __RSparseMapMatrix_Neg__


def __stdVectorRSparseMapMatrix_Mult__(self, b):
    ret = None
    if isinstance(b, list):
        # ret = np.array([pgcore.mult(self, bi) for bi in b]).T
        # print(ret)

        br = b[0]
        for i in range(1, len(b)):
            br = pgcore.cat(br, b[i])

        ret = pgcore.stdVectorR3Vector()
        pgcore.mult(self, br, ret)

    elif isinstance(b, np.ndarray) and b.ndim == 2:
        ret = pgcore.stdVectorR3Vector()
        pgcore.mult(self, b.flatten('F'), ret)

    elif isArray(b):
        ret = pgcore.stdVectorRVector()
        pgcore.mult(self, b, ret)

    else:
        print(b)
        pg.critical("Don't know how to convert b")

    return ret
pgcore.stdVectorRSparseMapMatrix.__mul__ = __stdVectorRSparseMapMatrix_Mult__


# Special Monkeypatch core classes
# keep original method
__BlockMatrix_addMatrix__ = pgcore.RBlockMatrix.addMatrix

def __BlockMatrix_addMatrix_happy_GC__(self, M, row=None, col=None,
                                       scale=1.0, transpose=False):
    """Add an existing matrix to this block matrix and return a unique index.

    As long row and col are None, the Matrix will not be used until a matrix
    entry is added.

    Monkeypatched version to increase the reference counter of M to keep the
    garbage collector happy.

    TODO
    ----
        * Add numpy matrices or convertable
        * Transpose is only for 1d arrays. Needed for matrices?

    Parameters
    ----------
    M: Matrix | pg.Vector | 1d iterable | scipy.sparse
        Matrix to add to the block.
    row: long
        Starting row index.
    col: long
        Starting column index.
    scale: float[1.0]
        Scale all matrix entries.
    transpose: bool [False]
        Transpose the matrix.
    """
    if M.ndim == 1:
        if transpose is False:
            _M = SparseMapMatrix(list(range(len(M))), [0]*len(M), M)
        else:
            warn('BlockMatrix add (transpose==True) ... Move me to core')
            _M = SparseMapMatrix([0]*len(M), list(range(len(M))), M)
        M = _M
    else:
        if transpose is True:
            if isinstance(M, pgcore.RSparseMapMatrix):
                warn('BlockMatrix add (transpose==True) ... Move me to core')
                v = pgcore.RVector()
                i = pgcore.IndexArray([0])
                j = pgcore.IndexArray([0])
                M.fillArrays(v, i, j)
                M = SparseMapMatrix(j, i, v)
            else:
                critical("don't know yet how to add transpose matrix of type",
                         type(M))

        from scipy.sparse import spmatrix

        if isinstance(M, spmatrix):
            M = toSparseMatrix(M)

    if not hasattr(self, '__mats__'):
        self.__mats__ = []
    self.__mats__.append(M)

    matrixID = __BlockMatrix_addMatrix__(self, M)
    # print('add:', matrixID, row, col)

    if row is not None and col is not None:
        self.addMatrixEntry(matrixID, row, col, scale)

    # for e in self.entries():
    #     print(e)
    return matrixID


def __BlockMatrix_str__(self):
    string = ("BlockMatrix of size %d x %d consisting of %d "
              "submatrices.")
    return string % (self.rows(), self.cols(), len(self.entries()))


pgcore.RBlockMatrix.addMatrix = __BlockMatrix_addMatrix_happy_GC__
pgcore.RBlockMatrix.add = __BlockMatrix_addMatrix_happy_GC__
pgcore.RBlockMatrix.__repr__ = __BlockMatrix_str__
pgcore.RBlockMatrix.ndim = 2
# pgcore.CBlockMatrix.addMatrix = __BlockMatrix_addMatrix_happy_GC__
# pgcore.CBlockMatrix.add = __BlockMatrix_addMatrix_happy_GC__


def __SparseMatrixEqual__(self, T):
    """Compare two SparseMatrices"""
    from .matrix import sparseMatrix2Array

    if self.shape[0] != T.shape[0] or self.shape[1] != T.shape[1]:
        warn(f'Compare sizes invalid {self.shape} vs. {T.shape}: ')
        return False

    rowsA, colsA, valsA = sparseMatrix2Array(self, indices=True, getInCRS=True)
    rowsB, colsB, valsB = sparseMatrix2Array(T, indices=True, getInCRS=True)

    # print(rowsA, colsA, valsA)
    # print(rowsB, colsB, valsB)

    if len(valsA) != len(valsB):
        warn(f'Compare value sizes invalid: {len(valsA)}, {len(valsB)}')
        return False

    meanA = np.mean(abs(valsA))
    nAB = np.linalg.norm(valsA-valsB)

    if np.isnan(nAB):
        print(valsA)
        print(valsB)
        error('norm(A-B) is nan')

    if pgcore.deepDebug() == -1:
        print(self, T)
        print('diff rows:', np.linalg.norm(np.array(rowsA)-np.array(rowsB)))
        print('diff cols:', np.linalg.norm(np.array(colsA)-np.array(colsB)))

        print(f'|A-B|={nAB}', f'mean(A)={meanA} mean(B)={np.mean(abs(valsB))}')

        print(f'nAB/meanA = {nAB/meanA}')

    if meanA > 1e-10:
        return np.all(rowsA == rowsB) and np.all(colsA == colsB) and nAB/meanA < 1e-12
    else:
        return nAB < 1e-12
pgcore.RSparseMatrix.__eq__ = __SparseMatrixEqual__
pgcore.RSparseMapMatrix.__eq__ = __SparseMatrixEqual__


def __SparseMatrixCopy__(self):
    """Create a copy."""
    return pgcore.RSparseMatrix(self)
pgcore.RSparseMatrix.copy = __SparseMatrixCopy__


def __SparseMapMatrixCopy__(self):
    """Create a copy."""
    return pgcore.RSparseMapMatrix(self)
pgcore.RSparseMapMatrix.copy = __SparseMapMatrixCopy__


def toSparseMatrix(A):
    """Convert any matrix type to pg.SparseMatrix and return copy of it.

    No conversion if A is a SparseMatrix already.

    Args
    ----
    A: pg or scipy matrix

    Returns
    -------
        pg.SparseMatrix
    """
    if isinstance(A, BlockMatrix):
        S = SparseMatrix()
        S.copy_(A.sparseMapMatrix())
        return S

    if isinstance(A, CSparseMapMatrix):
        S = CSparseMatrix()
        S.copy_(A)
        return S

    if isinstance(A, CSparseMatrix):
        return A

    if isinstance(A, SparseMatrix):
        return A

    if isinstance(A, SparseMapMatrix):
        S = SparseMatrix()
        S.copy_(A)
        return S

    if isinstance(A, csr_matrix):
        #pg.core.setDeepDebug(1)
        if len(A.data) == 0:
            print(A.indptr, A.indices, A.data)
            pg.error('Empty matrix conversion')
            return SparseMatrix(A.shape[0], A.shape[1])

        r, c = A.shape[0], A.shape[1]
        with pg.tictoc('ToSparseMatrix'):

            A = SparseMatrix(np.asarray(A.indptr, dtype=np.uint64),
                             np.asarray(A.indices, dtype=np.uint64),
                             A.data)

        ## sometimes zero diagnonals are skipped
        A.resize(r, c)
        #pg.core.setDeepDebug(0)
        return A

    return asSparseMatrix(asCSR(A))


def toSparseMapMatrix(A):
    """Convert any matrix type to pg.SparseMatrix and return copy of it.

    TODO
    ----
        * Always copy?  why?

    Arguments
    ---------
    A: pg or scipy matrix

    Returns
    -------
        pg.SparseMatrix
    """
    if isinstance(A, SparseMapMatrix):
        return SparseMapMatrix(A)

    if isinstance(A, BlockMatrix):
        return A.sparseMapMatrix()

    if isinstance(A, SparseMatrix):
        S = SparseMapMatrix()
        S.copy_(A)
        return S

    if isinstance(A, CSparseMatrix):
        S = CSparseMapMatrix()
        S.copy_(A)
        return A

    if isinstance(A, csr_matrix):
        return toSparseMapMatrix(toSparseMatrix(A))

    if isinstance(A, csc_matrix):
        return toSparseMapMatrix(A.tocsr())

    if isinstance(A, coo_matrix):
        return SparseMapMatrix(np.asarray(A.row, dtype=np.uint64),
                               np.asarray(A.col, dtype=np.uint64),
                               A.data)

    return toSparseMapMatrix(toCOO(A))

def toCSR(A):
    return sparseMatrix2csr(A)

def toCSC(A):
    return sparseMatrix2csc(A)

def toCOO(A):
    return sparseMatrix2coo(A)


def asCSR(A):
    return sparseMatrix2csr(A)

def asCSC(A):
    return sparseMatrix2csc(A)

def asCOO(A):
    return sparseMatrix2coo(A)

def asSparseMapMatrix(A):
    return toSparseMapMatrix(A)

def asSparseMatrix(A):
    return toSparseMatrix(A)


def sparseMatrix2csr(A):
    """Convert SparseMatrix to scipy.csr_matrix.

    Compressed Sparse Row matrix, i.e., Compressed Row Storage (CRS)

    Parameters
    ----------
    A: SparseMapMatrix | SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.csr_matrix
        Matrix to convert into.
    """
    #optImport(scipy.sparse, requiredFor="toCRS_matrix")


    if isinstance(A, CSparseMapMatrix):
        return toCOO(A).tocsr()

    if isinstance(A, SparseMapMatrix):
        return toCOO(A).tocsr()
        # C = pg.utils.toSparseMatrix(A)
        # return csr_matrix((C.vecVals().array(),
        #                    C.vecRowIdx(),
        #                    C.vecColPtr()))
    elif isinstance(A, SparseMatrix):
        return csr_matrix((A.vecVals().array(),
                           A.vecRowIdx().array(),
                           A.vecColPtr().array()), shape=A.shape)

    elif isinstance(A, CSparseMatrix):
        return csr_matrix((A.vecVals().array(),
                           A.vecRowIdx().array(),
                           A.vecColPtr().array()), shape=A.shape, dtype=complex)
        return csr
    elif isinstance(A, BlockMatrix):
        M = A.sparseMapMatrix()
        warn('bad efficency BlockMatrix->csr')
        return sparseMatrix2csr(M)

    return csr_matrix(A)

def sparseMatrix2csc(A):
    """Convert SparseMatrix to scipy.csr_matrix.

    Compressed Sparse Column matrix, i.e., Compressed Column Storage (CCS)

    Parameters
    ----------
    A: pg.matrix.SparseMapMatrix | pg.matrix.SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.csc_matrix
        Matrix to convert into.
    """
    #optImport(scipy.sparse, requiredFor="toCRC_matrix")

    if isinstance(A, CSparseMapMatrix):
        return toCOO(A).tocsc()
    if isinstance(A, SparseMapMatrix):
        return toCOO(A).tocsc()
    elif isinstance(A, SparseMatrix):
        return toCSR(A).tocsc()
    elif isinstance(A, CSparseMatrix):
        return toCSR(A).tocsc()

    elif isinstance(A, BlockMatrix):
        M = A.sparseMapMatrix()
        warn('bad efficency BlockMatrix->csc')
        return sparseMatrix2csc(M)

    return csc_matrix(A)


def sparseMatrix2coo(A, rowOffset=0, colOffset=0):
    """Convert SparseMatrix to scipy.coo_matrix.

    Parameters
    ----------
    A: SparseMapMatrix | SparseMatrix
        Matrix to convert from.

    Returns
    -------
    mat: scipy.coo_matrix
        Matrix to convert into.
    """
    vals = Vector()
    rows = IndexArray([0])
    cols = IndexArray([0])

    if isinstance(A, SparseMatrix):
        C = toSparseMapMatrix(A)
        C.fillArrays(vals=vals, rows=rows, cols=cols)
        rows += rowOffset
        cols += colOffset
        return coo_matrix((vals, (rows, cols)), shape=(A.rows(), A.cols()))

    elif isinstance(A, SparseMapMatrix):
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

    matrix: SparseMapMatrix or SparseMatrix
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
    if getInCRS is True:
        mat = asCSR(matrix)
        rows = mat.indices
        cols = mat.indptr
        vals = mat.data
    else:
        mat = asCOO(matrix)
        rows = mat.row
        cols = mat.col
        vals = mat.data

    if indices is True:
        return rows, cols, vals
    else:
        return vals


def sparseMatrix2Dense(matrix):
    """Convert sparse matrix to dense ndarray"""

    rr, cc, vals = sparseMatrix2Array(matrix, indices=True, getInCRS=False)
    mat = np.zeros((matrix.rows(), matrix.cols()))

    for i, r in enumerate(rr):
        mat[r, cc[i]] = vals[i]

    return mat
toDense = sparseMatrix2Dense
asDense = sparseMatrix2Dense


def removeEntries(A, rows=None, cols=None):
    """Remove rows and columns from a sparse matrix.

    Convert matrix into SCR

    """
    A = toCSR(A)
    rMask = np.ones(A.shape[0], dtype=bool)
    cMask = np.ones(A.shape[1], dtype=bool)

    if rows is not None:
        rows = list(rows)
    else:
        rows = []

    if cols is not None:
        cols = list(cols)
    else:
        cols = []

    if len(rows) > 0 and len(cols) > 0:
        rMask[rows] = False
        cMask[cols] = False
        return A[rMask][:,cMask]
    elif len(rows) > 0:
        rMask[rows] = False
        return A[rMask]
    elif len(cols) > 0:
        cMask[cols] = False
        return A[:,cMask]

    return A


def reduceEntries(A, idx):
    """Remove all values of A in row[idx] and col[idx] and add 1 to diag """
    import pygimli as pg
    debug = False
    if debug:
        pg.tic()

    if 1:
        try:
            with pg.tictoc('reduce'):
                #raise Exception('not yet implemented')
                A.reduce(idx, True)
        except:
            with pg.tictoc('manual reduce'):
                for i, ix in enumerate(idx):
                    with pg.tictoc('rows'):
                        A.cleanRow(ix)
                    with pg.tictoc('cols'):
                        A.cleanCol(ix)

        with pg.tictoc('diag'):
            for i, ix in enumerate(idx):
                A.setVal(ix, ix, 1.0)
    else:
        with pg.tictoc('->CSC'):
            csc = toCSC(A)

        with pg.tictoc('->col*=0'):
            csc[:,idx] *= 0

        with pg.tictoc('->CSR'):
            csr = csc.tocsr()

        with pg.tictoc('->row*=0'):
            csr[idx] *= 0

        # with pg.tictoc('->zero=0'):
        #     csr.eliminate_zeros()

        with pg.tictoc('<-copy'):
            A.copy_(asSparseMatrix(csr))

        with pg.tictoc('diag'):
            for i, ix in enumerate(idx):
                A.setVal(ix, ix, 1.0)
