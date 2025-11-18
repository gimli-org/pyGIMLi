# -*- coding: utf-8 -*-
"""Some matrix specialization."""

import time
import numpy as np

from .core import pgcore

from .core import (CMatrix, CSparseMapMatrix, CSparseMatrix,
                   RSparseMapMatrix, RSparseMatrix, ElementMatrix,
                   IVector, MatrixBase, R3Vector, RVector)

from .logger import critical, warn, info


# make pygimli core (pgcore) matrices available here for clarity
# useful aliases
IdentityMatrix = pgcore.IdentityMatrix

Matrix = pgcore.RMatrix
SparseMatrix = pgcore.RSparseMatrix
SparseMapMatrix = pgcore.RSparseMapMatrix
BlockMatrix = pgcore.RBlockMatrix

# General Monkeypatch core classes
__Matrices = [pgcore.MatrixBase,
              pgcore.RSparseMatrix,
              pgcore.RSparseMapMatrix,
              pgcore.CSparseMatrix,
              pgcore.CSparseMapMatrix,
              pgcore.RBlockMatrix,
              pgcore.IdentityMatrix]

for m in __Matrices:
    m.ndim = 2
    m.__len__ = lambda self: self.rows()
    m.shape = property(lambda self: (self.rows(), self.cols()))
    ## To allow for ignoring np.*.__mul__ in case A has the __rmul__ function
    ## see test TestMatrixMethods.testSparseMatrixBasics, TestRVectorMethods.testUFunc
    m.__array_ufunc__ = None
    m.__rmul__ = lambda self: critical('not yet implemented')


pgcore.RMatrix.dtype = float
pgcore.CMatrix.dtype = complex
pgcore.RSparseMapMatrix.dtype = float
pgcore.CSparseMapMatrix.dtype = complex
pgcore.RSparseMatrix.dtype = float
pgcore.CSparseMatrix.dtype = complex


def __RMatrix_str(self):
    s = "RMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__() + '\n'
    return s


def __CMatrix_str(self):
    s = "CMatrix: " + str(self.rows()) + " x " + str(self.cols())

    if self.rows() < 6:
        s += '\n'
        for v in range(self.rows()):
            s += self[v].__str__() + '\n'
    return s


def __ElementMatrix_str(self):
    """Show entries of an ElementMatrix."""
    import pygimli as pg
    if self.mat().cols() == 0 and self.mat().rows() == 0:
        return 'Empty ElementMatrix\n'

    maxRowID = int(np.log10(max(self.rowIDs())))+2

    s = '\n ' + ' ' * maxRowID
    for i in range(self.mat().cols()):
        s += str(self.colIDs()[i]).rjust(9)
    s += '\n'

    s += '  ' + '-'*self.mat().cols()*(9 + maxRowID) + '-\n'

    for i in range(self.mat().rows()):
        s += str(self.rowIDs()[i]).rjust(maxRowID) + " :"
        for v in self.row(i):
            s += pg.pf(v).rjust(9)
        s += '\n'
    return s


pgcore.RMatrix.__repr__ = __RMatrix_str
pgcore.CMatrix.__repr__ = __CMatrix_str
pgcore.ElementMatrix.__repr__ = __ElementMatrix_str


def __RSparseMapMatrix_str(self):
    s = "SparseMapMatrix: " + str(self.rows()) + " x " + str(self.cols())
    if self.rows() * self.cols() > 0:
        pc = int(self.nVals()/self.cols()/self.rows()*1000) / 10
        s += " (nnz=" + str(self.nVals()) + " / " + str(pc)+ "%)"

    return s


pgcore.RSparseMapMatrix.__repr__ = __RSparseMapMatrix_str

def __RSparseMatrix_str(self):
    s = "SparseMatrix: " + str(self.rows()) + " x " + str(self.cols())
    if self.rows() * self.cols() > 0:
        pc = int(self.nVals()/self.cols()/self.rows()*1000) / 10
        s += " (nnz=" + str(self.nVals()) + " / " + str(pc)+ "%)"

    return s


pgcore.RSparseMatrix.__repr__ = __RSparseMatrix_str


def __RVector_format(self, f):
    print(f)
    return str(self)

pgcore.RVector.__format__ = __RVector_format

# Special Monkeypatch core classes
__BlockMatrix_addMatrix__ = pgcore.RBlockMatrix.addMatrix


def __BlockMatrix_addMatrix_happy_GC__(self, M, row=None, col=None,
                                       scale=1.0, transpose=False):
    """Add an existing matrix to this block matrix and return a unique index.

    As long row and col are None, the Matrix will not be used until a matrix
    entry is has been added.

    Monkeypatched version to increase the reference counter of M to keep the
    garbage collector happy.

    TODO
    ----
        * Add numpy matrices or convertable
        * Transpose is only for 1d arrays. Needed for matrices?

    Parameters
    ----------
    M: pg.matrix.Matrix | pg.Vector | 1d iterable
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

    if not hasattr(self, '__mats__'):
        self.__mats__ = []
    self.__mats__.append(M)

    matrixID = __BlockMatrix_addMatrix__(self, M)

    if row is not None and col is not None:
        self.addMatrixEntry(matrixID, row, col, scale)

    return matrixID


def __BlockMatrix_str__(self):
    string = ("pg.matrix.BlockMatrix of size %d x %d consisting of %d "
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
    from pygimli.utils import sparseMatrix2Array
    if self.rows() != T.rows() or self.cols() != T.cols():
        warn("Compare sizes invalid {0},{1} vs. {2},{3}: ".format(
            self.rows(), self.cols(), T.rows(), T.cols()))
        return False

    rowsA, colsA, valsA = sparseMatrix2Array(self, indices=True)
    rowsB, colsB, valsB = sparseMatrix2Array(T, indices=True)

    if len(valsA) != len(valsB):
        warn("Compare value sizes invalid: ", len(valsA), len(valsB))
        return False

    # print(self, T)
    # print('rows:', np.linalg.norm(np.array(rowsA)-np.array(rowsB)))
    # print('cols:', np.linalg.norm(np.array(colsA)-np.array(colsB)))

    # print(np.linalg.norm(valsA-valsB), np.mean(abs(valsA)),
    #     np.mean(abs(valsB)))
    # print(np.linalg.norm(valsA-valsB)/np.mean(abs(valsA)))

    return rowsA == rowsB and \
        colsA == colsB and \
        np.linalg.norm(valsA-valsB)/np.mean(abs(valsA)) < 1e-14


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


# class MultMatrix(pgcore.MatrixBase):
#     """Base Matrix class for all matrix types holding a matrix."""

#     def __init__(self, A, verbose=False):
#         self._A = A
#         self.ndim = self._A.ndim
#         super(MultMatrix, self).__init__(verbose)

#     @property
#     def A(self):
#         return self._A

#     @A.setter
#     def A(self, A):
#         self._A = A

#     def rows(self):
#         """Return number of rows (using underlying matrix)."""
#         return self.A.rows()  # this should be _A

#     def cols(self):
#         """Return number of columns (using underlying matrix)."""
#         return self.A.cols()  # this should be _A

#     def save(self, filename):
#         """So it can be used in inversion with dosave flag"""
#         pass


# class MultLeftMatrix(MultMatrix):
#     """Matrix consisting of actual RMatrix and lef-side vector."""

#     def __init__(self, A, left, verbose=False):
#         """Constructor saving matrix and vector."""
#         if A.rows() != len(left):
#             raise Exception("Matrix columns do not fit vector length!")
#         super(MultLeftMatrix, self).__init__(A, verbose)

#         self._l = left

#     @property
#     def l(self):  # better use left and right instead (pylint E743)?
#         return self._l

#     @l.setter
#     def r(self, l):
#         self._l = l

#     def mult(self, x):
#         """Multiplication from right-hand-side (dot product A*x)."""
#         return self.A.mult(x) * self.l

#     def transMult(self, x):
#         """Multiplication from right-hand-side (dot product A.T * x)"""
#         return self.A.transMult(x * self.l)


# LMultRMatrix = MultLeftMatrix  # alias for backward compatibility


# class MultRightMatrix(MultMatrix):
#     """Some Matrix, multiplied with a right hand side vector r."""

#     def __init__(self, A, r=None, verbose=False):
#         super(MultRightMatrix, self).__init__(A, verbose)

#         if r is None:
#             self._r = pgcore.RVector(self.cols(), 1.0)
#         else:
#             self._r = r

#     @property
#     def r(self):
#         return self._r

#     @r.setter
#     def r(self, r):
#         self._r = r

#     def mult(self, x):
#         """Return M*x = A*(r*x)"""
#         if hasattr(x, '__len__') and hasattr(self.r, '__len__'):
#             if len(x) != len(self.r):
#                 # assuming A was complex
#                 # warn('need to double x')
#                 # print('mult:', self.A.rows(), " x " , self.A.cols(),
#                 #        'x:', len(x), 'r:', len(self.r))
#                 # print(self.perm)
#                 return self.A.mult(x[self.perm] * self.r)
#                 # return self.A.mult(pgcore.cat(x, x) * self.r)
#         return self.A.mult(x * self.r)

#     def transMult(self, x):
#         """Return M.T*x=(A.T*x)*r"""
#         # print('transmult', self.A.rows(), " x " , self.A.cols(), x, self.r, )
#         return self.A.transMult(x) * self.r


# RMultRMatrix = MultRightMatrix  # alias for backward compatibility


# class MultLeftRightMatrix(MultMatrix):
#     """Matrix consisting of actual RMatrix and left-hand-side vector."""

#     def __init__(self, A, left, right, verbose=False):
#         """Constructor saving matrix and vector."""
#         if A.cols() != len(right):
#             raise Exception("Matrix columns do not fit right vector length!")
#         if A.rows() != len(left):
#             raise Exception("Matrix rows do not fit left vector length!")

#         super(MultLeftRightMatrix, self).__init__(A, verbose)
#         self._r = right
#         self._l = left

#     @property
#     def l(self):
#         return self._l

#     @l.setter
#     def l(self, l):
#         self._l = l

#     @property
#     def r(self):
#         return self._r

#     @r.setter
#     def r(self, r):
#         self._r = r

#     def mult(self, x):
#         """Multiplication from right-hand-side (dot product A*x)."""
#         return self.A.mult(x * self._r) * self._l

#     def transMult(self, x):
#         """Multiplication from right-hand-side (dot product A.T*x)."""
#         return self.A.transMult(x * self._l) * self._r


# LRMultRMatrix = MultLeftRightMatrix  # alias for backward compatibility


# class Add2Matrix(pgcore.MatrixBase):
#     """Matrix by adding two matrices."""

#     def __init__(self, A, B):
#         super().__init__()
#         self.A = A
#         self.B = B
#         assert A.rows() == B.rows()
#         assert A.cols() == B.cols()

#     def mult(self, x):
#         """Return M*x = A*(r*x)"""
#         return self.A.mult(x) + self.B.mult(x)

#     def transMult(self, x):
#         """Return M.T*x=(A.T*x)*r"""
#         return self.A.transMult(x) + self.B.transMult(x)

#     def cols(self):
#         """Number of columns."""
#         return self.A.cols()

#     def rows(self):
#         """Number of rows."""
#         return self.A.rows()


# class Mult2Matrix(pgcore.MatrixBase):
#     """Matrix by multiplying two matrices.
#         M*x = A * (B*x)
#         M.T*x = (A*x) * B
#     """

#     def __init__(self, A, B):
#         super().__init__()
#         self.A = A
#         self.B = B
#         assert A.cols() == B.rows()

#     def mult(self, x):
#         """Return M*x = A*(B*x)"""
#         return self.A.mult(self.B.mult(x))

#     def transMult(self, x):
#         """Return M.T*x=(A.T*x)*B"""
#         return self.B.transMult(self.A.transMult(x))

#     def cols(self):
#         """Number of columns."""
#         return self.B.cols()

#     def rows(self):
#         """Number of rows."""
#         return self.A.rows()


# class DiagonalMatrix(pgcore.MatrixBase):
#     """Square matrix with a vector on the main diagonal."""

#     def __init__(self, d):
#         super().__init__()
#         self.d = d

#     def mult(self, x):
#         """Return M*x = d*x (element-wise)"""
#         return x * self.d

#     def transMult(self, x):
#         """Return M.T*x = M*x"""
#         return x * self.d

#     def cols(self):
#         """Number of columns (length of diagonal)."""
#         return len(self.d)

#     def rows(self):
#         """Number of rows (length of diagonal)."""
#         return len(self.d)

# @pg.cache
# def createCm05(A):
#     """Globally cached helper function to create Cm05Matrix."""
#     return pg.matrix.Cm05Matrix(A, verbose=True)

# class Cm05Matrix(pgcore.MatrixBase):
#     """Matrix implicitly representing the inverse square-root."""

#     def __init__(self, A, verbose=False):
#         """Constructor saving matrix and vector.

#         Parameters
#         ----------
#         A : ndarray
#             numpy type (full) matrix
#         """
#         super().__init__(verbose)  # only in Python 3
#         self._mul = None

#         if isinstance(A, str):
#             self.load(A)
#         else:
#             from scipy.linalg import eigh  # , get_blas_funcs

#             if A.shape[0] != A.shape[1]:  # rows/cols for pgcore matrix
#                 raise Exception("Matrix must by square (and symmetric)!")

#             if verbose:
#                 t = time.perf_counter()
#             self.ew, self.EV = eigh(A)

#             if verbose:
#                 info('(C) Time for eigenvalue decomposition: {:.1f}s'.format(
#                     time.perf_counter()-t))

#             #self.A = A

#     @property
#     def mul(self):
#         if self._mul is None:
#             self._mul = np.sqrt(1./self.ew)
#         return self._mul

#     def save(self, fileName):
#         """Save the content of this matrix. Used for caching until pickling is possible for this class"""
#         np.save(fileName, dict(ew=self.ew, EV=self.EV), allow_pickle=True)

#     def load(self, fileName):
#         """Load the content of this matrix. Used for caching until pickling is possible for this class"""
#         d = np.load(fileName + '.npy', allow_pickle=True).tolist()
#         self.ew = d['ew']
#         self.EV = d['EV']

#     def rows(self):
#         """Return number of rows (using underlying matrix)."""
#         return len(self.ew)

#     def cols(self):
#         """Return number of columns (using underlying matrix)."""
#         return self.row()

#     def mult(self, x):
#         """Multiplication from right-hand side (dot product)."""
#         part1 = (np.dot(np.transpose(x), self.EV).T*self.mul).reshape(-1, 1)
#         return self.EV.dot(part1).reshape(-1,)
# #        return self.EV.dot((x.T.dot(self.EV)*self.mul).T)

#     def transMult(self, x):
#         """Multiplication from right-hand side (dot product)."""
#         return self.mult(x)  # matrix is symmetric by definition


# class NDMatrix(BlockMatrix):
#     """Diagonal block (block-Jacobi) matrix derived from pg.matrix.BlockMatrix.

#     (to be moved to a better place at a later stage)
#     """

#     def __init__(self, num, nrows, ncols):
#         super(NDMatrix, self).__init__()  # call inherited init function
#         self.Ji = []  # list of individual block matrices
#         for i in range(num):
#             self.Ji.append(pgcore.Matrix())
#             self.Ji[-1].resize(nrows, ncols)
#             n = self.addMatrix(self.Ji[-1])
#             self.addMatrixEntry(n, nrows * i, ncols * i)

#         self.recalcMatrixSize()


# class GeostatisticConstraintsMatrix(pgcore.MatrixBase):
#     """Geostatistic constraints matrix

#     Uses geostatistical operators described by Jordi et al. (2018),
#     however corrects for the remaining non-smooth (damping) part by
#     correcting for the spur of the inverse root matrix.

#     Jordi, C., Doetsch, J., Günther, T., Schmelzbach, C. & Robertsson, J.O.A.
#     (2018): Geostatistical regularisation operators for geophysical inverse
#     problems on irregular meshes. Geoph. J. Int. 213, 1374-1386,
#     doi:10.1093/gji/ggy055.
#     """

#     def __init__(self, CM=None, mesh=None, **kwargs):
#         """Initialize by computing the covariance matrix & its inverse root.

#         Parameters
#         ----------
#         CM : pg.Matrix or pg.SparseMapMatrix
#             covariance matrix, if not given, use mesh and I
#         mesh : pg.Mesh
#             mesh of which the cell midpoints are used for covariance
#         I : float | iterable of floats
#             axis correlation length (isotropic) or lengths (anisotropic)
#         dip : float [0]
#             angle of main axis corresponding to I[0] (2D) or I[0]&I[1] (3D)
#         strike : float [0]
#             angle of main axis corresponding to I[0] versus I[1] (3D)
#         withRef : bool [False]
#             neglect spur (reference model effect) that is otherwise corrected
#         """
#         super().__init__(kwargs.pop('verbose', False))
#         self.withRef = kwargs.pop('withRef', False)
#         self._spur = None

#         if isinstance(CM, str):
#             self.load(CM)
#         else:
#             from pygimli.utils.geostatistics import covarianceMatrix

#             if isinstance(CM, pgcore.Mesh):
#                 CM = covarianceMatrix(CM, **kwargs)

#             if CM is None and mesh is not None:
#                 CM = covarianceMatrix(mesh, **kwargs)
#             else:
#                 pg.critical('Give either CM or mesh')

#             self.Cm05 = createCm05(CM)

#     @property
#     def spur(self):
#         if self._spur is None:
#             if self.withRef is True:
#                 self._spur = np.zeros(self.rows())
#             else:
#                 self._spur = self.Cm05 * pgcore.RVector(self.rows(), 1.0)
#         return self._spur

#     @property
#     def nModel(self):
#         try:
#             return self.Cm05.size()
#         except Exception as e:
#             return 0

#     def save(self, fileName):
#         """Save the content of this matrix. Used for caching until pickling is possible for this class
#         """
#         self.Cm05.save(fileName + '-Cm05')
#         np.save(fileName, dict(verbose=self.verbose(),
#                                withRef=self.withRef,
#                                Cm05=fileName +'-Cm05'),
#                         allow_pickle=True)


#     def load(self, fileName):
#         """Load the content of this matrix. Used for caching until pickling is possible for this class
#         """
#         d = np.load(fileName + '.npy', allow_pickle=True).tolist()
#         self.setVerbose(d['verbose'], )
#         self.withRef = d['withRef']
#         self.Cm05 = Cm05Matrix(d['Cm05'])


#     def mult(self, x):
#         return self.Cm05.mult(x) - self.spur * x

#     def transMult(self, x):
#         return self.Cm05.transMult(x) - self.spur * x

#     def cols(self):
#         return self.nModel

#     def rows(self):
#         return self.nModel

#     def clear(self):
#         self.Cm05 = None
#         self._spur = None
