"""Some special and usefull matrices."""

import numpy as np

import pygimli as pg
import pygimli.core as pgcore

from pygimli.core import (CMatrix, CSparseMapMatrix, CSparseMatrix,
                          RSparseMapMatrix, RSparseMatrix, ElementMatrix,
                          MatrixBase)

IdentityMatrix = pgcore.IdentityMatrix
Matrix = pgcore.RMatrix
SparseMatrix = pgcore.RSparseMatrix
SparseMapMatrix = pgcore.RSparseMapMatrix
BlockMatrix = pgcore.RBlockMatrix


class ScaledMatrix(MatrixBase):
    """Super-simple matrix with scaling factor B=A*s."""

    def __init__(self, A, scale=1.0, verbose=False):
        super().__init__(verbose)
        self._A = A
        self._scale = scale

    def rows(self):
        """Return number of rows (forwarded)."""
        return self._A.rows()

    def cols(self):
        """Return number of cols (forwarded)."""
        return self._A.cols()

    def mult(self, x):
        """Multiplication from right-hand-side (A.T*x)."""
        return self._A.mult(x) * self._scale

    def transMult(self, x):
        """Multiplication from right-hand-side (A*x)."""
        return self._A.transMult(x) * self._scale


class TransposedMatrix(MatrixBase):
    """Wrapper for transposed matrix (of any kind)."""

    def __init__(self, A, verbose=False):
        super().__init__(verbose)
        self._A = A

    def rows(self):
        """Return number of rows (cols of underlying matrix)."""
        return self._A.cols()

    def cols(self):
        """Return number of cols (rows of underlying matrix)."""
        return self._A.rows()

    def mult(self, x):
        """Multiplication from right-hand-side (A.T*x)."""
        return self._A.transMult(x)

    def transMult(self, x):
        """Multiplication from right-hand-side (A*x)."""
        return self._A.mult(x)


class SquaredMatrix(MatrixBase):
    """Wrapper for squared (least-squares) matrix A^T*A (of any kind)."""

    def __init__(self, A, verbose=False):
        super().__init__(verbose)
        self._A = A

    def rows(self):
        """Return number of rows (cols of underlying matrix)."""
        return self._A.cols()

    def cols(self):
        """Return number of cols (cols of underlying matrix)."""
        return self._A.cols()

    def mult(self, x):
        """Multiplication from right-hand-side (A.T*A*x)."""
        return self._A.transMult(self._A.mult(x))

    def transMult(self, x):
        """Multiplication from right-hand-side ((A^T*A)^T*x)."""
        return self._A.transMult(self._A.mult(x))


class SquaredTransposeMatrix(MatrixBase):
    """Wrapper for squared (least-squares) matrix A*A^T (of any kind)."""

    def __init__(self, A, verbose=False):
        super().__init__(verbose)
        self._A = A

    def rows(self):
        """Return number of rows (rows of underlying matrix)."""
        return self._A.rows()

    def cols(self):
        """Return number of cols (rows of underlying matrix)."""
        return self._A.rows()

    def mult(self, x):
        """Multiplication from right-hand-side (A.T*A*x)."""
        return self._A.mult(self._A.transMult(x))

    def transMult(self, x):
        """Multiplication from right-hand-side ((A^T*A)^T*x)."""
        return self._A.mult(self._A.transMult(x))


class RealNumpyMatrix(MatrixBase):
    """Matrix Wrapper for for ndarrays, providing syntax for pygimli c++ core.

    Holds reference to a real matrix, providing the correct
    multiplication algorithms for the pygimli inversion process.
    """

    def __init__(self, mat, copy=False):
        super().__init__()

        self.ndim = 2

        # print('real Matrix')
        if isinstance(mat, str):
            self.M = np.load(mat)
        else:
            if copy is True:
                self.M = np.copy(mat)
            else:
                self.M = mat

    def __repr__(self):
        """Return string representation of the matrix."""
        return self.M.__repr__()

    def rows(self):
        """Return number of rows."""
        return np.shape(self.M)[0]

    def cols(self):
        """Return number of columns."""
        return np.shape(self.M)[1]

    def save(self, name):
        """Save the content of this matrix."""
        np.save(name, self.M)

    def mult(self, vector):
        """Multiplication from right-hand side."""
        return self.M.dot(vector)

    def transMult(self, vector):
        """Multiplication from left-hand side."""
        return np.dot(vector, self.M)


def ComplexNumpyMatrix(mat, copy=False):
    """Return numpy matrix with real and imaginary parts stacked vertically."""
    return RealNumpyMatrix(np.row_stack([mat.real, mat.imag]), copy=copy)


def NumpyMatrix(mat, copy=False):
    """Return numpy-based pyGIMLi matrix (real or complex)."""
    if isinstance(mat, str):
        mat = np.load(mat)

    if isinstance(mat[0][0], complex):
        return ComplexNumpyMatrix(mat, copy=copy)
    else:
        return RealNumpyMatrix(mat, copy=copy)


class MultMatrix(MatrixBase):
    """Base Matrix class for all matrix types holding a matrix."""

    def __init__(self, A, verbose=False):
        self._A = A
        self.ndim = self._A.ndim
        super().__init__(verbose)

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self, A):
        self._A = A

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return self.A.rows()  # this should be _A

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self.A.cols()  # this should be _A

    def save(self, filename):
        """Save it to be used in inversion with dosave flag."""
        pass


class MultLeftMatrix(MultMatrix):
    """Matrix consisting of actual RMatrix and lef-side vector."""

    def __init__(self, A, left, verbose=False):
        """Initialize by saving matrix and vector."""
        if A.rows() != len(left):
            raise ValueError("Matrix columns do not fit vector length!")
        super().__init__(A, verbose)

        self._l = left

    @property
    def l(self):  # better use left and right instead (pylint E743)?
        return self._l

    @l.setter
    def r(self, l):
        self._l = l

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self.A.mult(x) * self.l

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T * x)."""
        return self.A.transMult(x * self.l)


LMultRMatrix = MultLeftMatrix  # alias for backward compatibility


class MultRightMatrix(MultMatrix):
    """Some Matrix, multiplied with a right hand side vector r."""

    def __init__(self, A, r=None, verbose=False):
        super().__init__(A, verbose)

        if r is None:
            self._r = pgcore.RVector(self.cols(), 1.0)
        else:
            self._r = r

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        self._r = r

    def mult(self, x):
        """Return M*x = A*(r*x)."""
        if hasattr(x, '__len__') and hasattr(self.r, '__len__') and \
            len(x) != len(self.r):
                # assuming A was complex
                # warn('need to double x')
                # print('mult:', self.A.rows(), " x " , self.A.cols(),
                #        'x:', len(x), 'r:', len(self.r))
                # print(self.perm)
                return self.A.mult(x[self.perm] * self.r)
                # return self.A.mult(pgcore.cat(x, x) * self.r)
        return self.A.mult(x * self.r)

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r."""
        # print('transmult', self.A.rows(), " x " , self.A.cols(), x, self.r, )
        return self.A.transMult(x) * self.r


RMultRMatrix = MultRightMatrix  # alias for backward compatibility


class MultLeftRightMatrix(MultMatrix):
    """Matrix consisting of actual RMatrix and left-hand-side vector."""

    def __init__(self, A, left, right, verbose=False):
        """Initialize matrix by matrix and vectors.

        Parameters
        ----------
        A : pg.matrix.MatrixBase derived matrix
            main matrix
        left, right : pg.Vector
            left and right side vectors to weight individual rows & columns
        """
        if A.cols() != len(right):
                raise IndexError("Matrix columns don't fit right vector length!",
                                A.cols(), len(right))
        if A.rows() != len(left):
            raise IndexError("Matrix rows don't fit left vector length!",
                            A.rows(), len(left))

        super().__init__(A, verbose)
        self._r = right
        self._l = left

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, l):
        self._l = l

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, r):
        self._r = r

    def mult(self, x):
        """Multiplication from right-hand-side (dot product A*x)."""
        return self.A.mult(x * self._r) * self._l

    def transMult(self, x):
        """Multiplication from right-hand-side (dot product A.T*x)."""
        return self.A.transMult(x * self._l) * self._r


LRMultRMatrix = MultLeftRightMatrix  # alias for backward compatibility


class Add2Matrix(MatrixBase):
    """Matrix by addition of two matrices implicitly.

    The matrix holds two matrices and distributes multiplication in 2 parts

    M = A + B
    M * x = A * x + B * x
    M.T * y = B^T * y + A.T^*y
    """

    def __init__(self, A, B):
        """Initialize matrix.

        Parameters
        ----------
        A, B : any Matrix derived from pg.core.MatrixBase
            submatrices to be multiplied (A/B.cols() and rows() need to match)
        """
        super().__init__()
        self.A = A
        self.B = B
        assert A.rows() == B.rows()
        assert A.cols() == B.cols()

    def mult(self, x):
        """Return M*x = A*(r*x)."""
        return self.A.mult(x) + self.B.mult(x)

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*r."""
        return self.A.transMult(x) + self.B.transMult(x)

    def cols(self):
        """Return number of columns."""
        return self.A.cols()

    def rows(self):
        """Return number of rows."""
        return self.A.rows()


class Mult2Matrix(MatrixBase):
    """Matrix by multipication of two matrices implicitly.

    The matrix holds two matrices and distributes multiplication in 2 parts

        M = A * B
        M * x = A * (B*x)
        M.T * y = B^T * (A.T^*y) = (y^T * A * B)^T
    """

    def __init__(self, A, B):
        """Initialize matrix.

        Parameters
        ----------
        A, B : any Matrix derived from pg.core.MatrixBase
            submatrices to be multiplied (A.cols() must equal B.rows())
        """
        super().__init__()
        self.A = A
        self.B = B
        assert A.cols() == B.rows()

    def mult(self, x):
        """Return M*x = A*(B*x)."""
        return self.A.mult(self.B.mult(x))

    def transMult(self, x):
        """Return M.T*x=(A.T*x)*B."""
        return self.B.transMult(self.A.transMult(x))

    def cols(self):
        """Return number of columns."""
        return self.B.cols()

    def rows(self):
        """Return number of rows."""
        return self.A.rows()


class DiagonalMatrix(MatrixBase):
    """Square matrix with a vector on the main diagonal."""

    def __init__(self, d):
        """Initialize matrix.

        Parameters
        ----------
        d : Vector
            vector holding diagonal elements
        """
        super().__init__()
        self.d = d

    def mult(self, x):
        """Return M*x = d*x (element-wise)."""
        return x * self.d

    def transMult(self, x):
        """Return M.T*x = M*x."""
        return x * self.d

    def cols(self):
        """Return number of columns (length of diagonal)."""
        return len(self.d)

    def rows(self):
        """Return number of rows (length of diagonal)."""
        return len(self.d)


@pg.cache
def createCm05(A):
    """Globally cached helper function to create Cm05Matrix."""
    return pg.matrix.Cm05Matrix(A, verbose=True)


class Cm05Matrix(MatrixBase):
    """Matrix implicitly representing the inverse square-root."""

    def __init__(self, A, trsh=1e-4, verbose=False):
        """Initialize by storing matrix and vector.

        Parameters
        ----------
        A : ndarray
            numpy type (full) matrix
        """
        super().__init__(verbose)
        self._mul = None
        self._trsh = trsh

        if isinstance(A, str):
            self.load(A)
        else:
            import scipy
            from packaging import version

            if A.shape[0] != A.shape[1]:  # rows/cols for pgcore matrix
                raise ValueError("Matrix must by square (and symmetric)!")

            if verbose:
                pg.tic(key='init cm05')

            eigkw = {}
            trsh = self._trsh
            if version.parse(scipy.__version__) >= version.parse("1.5"):
                eigkw["subset_by_value"] = [trsh, np.inf]
                eigkw["driver"] = "evr"
                self.ew, self.EV = scipy.linalg.eigh(A, **eigkw)
            else:
                self.ew, self.EV = scipy.linalg.eigh(A)
                self.EV = self.EV[:, self.ew > trsh]
                self.ew = self.ew[self.ew > trsh]

            if verbose:
                pg.info('(C) Time for eigenvalue decomposition {:.1f}s'.format(
                        pg.dur(key='init cm05')))

    @property
    def mul(self):
        """Scaling vector for A**(-0.5), 1/sqrt(eigenvalues)."""
        if self._mul is None:
            self._mul = np.sqrt(1./self.ew)
        return self._mul

    def save(self, fileName):
        """Save the content of this matrix.

        Used for caching until pickling is possible for this class
        """
        np.save(fileName, dict(ew=self.ew, EV=self.EV), allow_pickle=True)

    def load(self, fileName):
        """Load the content of this matrix.

        Used for caching until pickling is possible for this class
        """
        d = np.load(fileName + '.npy', allow_pickle=True).tolist()
        self.ew = d['ew']
        self.EV = d['EV']

    def rows(self):
        """Return number of rows (using underlying matrix)."""
        return len(self.ew)

    def cols(self):
        """Return number of columns (using underlying matrix)."""
        return self.rows()

    def mult(self, x):
        """Multiplication from right-hand side (dot product)."""
        return self.EV.dot(np.dot(np.transpose(x), self.EV).T*self.mul)

    def transMult(self, x):
        """Multiplication from right-hand side (dot product)."""
        return self.mult(x)  # matrix is symmetric by definition


class RepeatVMatrix(BlockMatrix):
    """Matrix repeating a base matrix N times vertically. Only A is stored.

    M = | A |
        | A |
        | A |
    """

    def __init__(self, A, num):
        """Init matrix by specify matrix and number of repetitions.

        Parameters
        ----------
        A : pg.MatrixBase (or derived)
            matrix for constraining a single frame (e.g. 1st order smoothness)
        num : int
            number of repetitions
        """
        super().__init__()
        self.A_ = A
        self.Aid = self.addMatrix(self.A_)
        nr = 0
        for _ in range(num):
            self.addMatrixEntry(self.Aid, nr, 0)
            nr += A.rows()

        self.recalcMatrixSize()


class RepeatHMatrix(BlockMatrix):
    """Matrix repeating a base matrix N times horizontally. Only A is stored.

    M = [ A A A ]
    """

    def __init__(self, A, num):
        """Init matrix by specify matrix and number of repetitions.

        Parameters
        ----------
        A : pg.MatrixBase (or derived)
            matrix for constraining a single frame (e.g. 1st order smoothness)
        num : int
            number of repetitions
        """
        super().__init__()
        self.A_ = A
        self.Aid = self.addMatrix(self.A_)
        nc = 0
        for _ in range(num):
            self.addMatrixEntry(self.Aid, 0, nc)
            nc += A.cols()

        self.recalcMatrixSize()


class RepeatDMatrix(BlockMatrix):
    """Matrix repeating a base matrix N times diagonally. Only A is stored.

    M = | A     |
        |   A   |
        |     A |
    """

    def __init__(self, A, num):
        """Init matrix by specify matrix and number of repetitions.

        Parameters
        ----------
        A : pg.MatrixBase (or derived)
            matrix for constraining a single frame (e.g. 1st order smoothness)
        num : int
            number of repetitions
        """
        super().__init__()
        self.A_ = A
        self.Aid = self.addMatrix(self.A_)
        nc = 0
        nr = 0
        for _ in range(num):
            self.addMatrixEntry(self.Aid, nr, nc)
            nc += A.cols()
            nr += A.rows()

        self.recalcMatrixSize()


class FrameConstraintMatrix(RepeatDMatrix):
    """Matrix constraining meshes inside and between frames of same mesh.

    M = |  A        |
        |     A     |
        |        A  |
        | -I +I     |
        |    -I +I  |
    """

    def __init__(self, A, num, scale=1.0):
        """Init matrix by specifying number of frames.

        Parameters
        ----------
        A : pg.MatrixBase (or derived)
            matrix for constraining a single frame (e.g. 1st order smoothness)
        num : int
            number of frames
        scale : float | iterable of size num - 1
        """
        super().__init__(A, num)
        self.nm = A.cols()  # model parameter per frame
        self.nb = A.rows()  # boundaries per frame
        if isinstance(scale, float):
            self.Iminus_ = pgcore.IdentityMatrix(self.nm, val=-scale)
            self.Iplus_ = pgcore.IdentityMatrix(self.nm, val=scale)
            self.Im = self.addMatrix(self.Iminus_)
            self.Ip = self.addMatrix(self.Iplus_)
            nc = self.rows()
            for i in range(num-1):
                self.addMatrixEntry(self.Im, nc, i*self.nm)
                self.addMatrixEntry(self.Ip, nc, (i+1)*self.nm)
                nc += self.nm
        elif hasattr(scale, "__iter__"):
            assert len(scale) == self.nm-1, "scaling vector length != num-1"
            self.diag_ = np.repeat(scale, self.nm)
            self.Iminus_ = pg.matrix.DiagonalMatrix(-self.diag_)
            self.Iplus_ = pg.matrix.DiagonalMatrix(self.diag_)
            nc = self.rows()
            self.addMatrix(self.Iminus_, nc, 0)
            self.addMatrix(self.Iplus_, nc, self.nm)
        else:
            raise TypeError("Could not use scale")

        self.recalcMatrixSize()


class NDMatrix(pgcore.BlockMatrix):
    """Diagonal block (block-Jacobi) matrix ."""

    def __init__(self, num, nrows, ncols):
        super().__init__()
        self.Ji = []
        for i in range(num):
            self.Ji.append(pgcore.Matrix())
            self.Ji[-1].resize(nrows, ncols)
            n = self.addMatrix(self.Ji[-1])
            self.addMatrixEntry(n, nrows * i, ncols * i)

        self.recalcMatrixSize()


class KroneckerMatrix(MatrixBase):
    """Memory-saving implementation of Kronecker matrix.

    The Kronecker matrix consists of repetitions of an inner
    matrix I, multiplied with elements of an outer matrix O.
        | O_11 I  O_12 I ... ]
    A = | O_21 I  O_22 I
        | ...
    """

    def __init__(self, outer, inner, verbose=False):
        """Initialize with outer and inner matrix."""
        super().__init__(verbose)
        self._I = inner
        self._O = outer
        self.ni = inner.rows()
        self.no = outer.rows()
        self.mi = inner.cols()
        self.mo = outer.cols()

    def rows(self):
        """Return number of rows (rows(I)*rows(O))."""
        return self._I.rows() * self._O.rows()

    def cols(self):
        """Return number of cols (cols(I)*cols(O))."""
        return self._I.cols() * self._O.cols()

    def mult(self, x):
        """Multiplication from right-hand-side (A.T*x)."""
        xx = np.reshape(x, [self.mo, self.mi])
        yy = np.zeros([self.no, self.ni])
        for i, xi in enumerate(xx):
            Ixi = self._I.mult(xi)
            for j in range(self.no):
                yy[j] += Ixi * self._O[j, i]

        return yy.ravel()

    def transMult(self, x):
        """Multiplication from right-hand-side (A*x)."""
        xx = np.reshape(x, [self.no, self.ni])
        yy = np.zeros([self.mo, self.mi])
        for i, xi in enumerate(xx):
            Ixi = self._I.transMult(xi)
            for j in range(self.mo):
                yy[j] += Ixi * self._O[i, j]

        return yy.ravel()

class GeostatisticConstraintsMatrix(MatrixBase):
    """Geostatistic constraints matrix.

    Uses geostatistical operators described by Jordi et al. (2018),
    however corrects for the remaining non-smooth (damping) part by
    correcting for the spur of the inverse root matrix.

    Jordi, C., Doetsch, J., Günther, T., Schmelzbach, C. & Robertsson, J.O.A.
    (2018): Geostatistical regularisation operators for geophysical inverse
    problems on irregular meshes. Geoph. J. Int. 213, 1374-1386,
    doi:10.1093/gji/ggy055.
    """

    def __init__(self, CM=None, mesh=None, **kwargs):
        """Initialize by computing the covariance matrix & its inverse root.

        Parameters
        ----------
        CM : pg.Matrix or pg.SparseMapMatrix
            covariance matrix, if not given, use mesh and I
        mesh : pg.Mesh
            mesh of which the cell midpoints are used for covariance
        I : float | iterable of floats
            axis correlation length (isotropic) or lengths (anisotropic)
        dip : float [0]
            angle of main axis corresponding to I[0] (2D) or I[0]&I[1] (3D)
        strike : float [0]
            angle of main axis corresponding to I[0] versus I[1] (3D)
        withRef : bool [False]
            neglect spur (reference model effect) that is otherwise corrected
        """
        super().__init__(kwargs.pop('verbose', False))
        self.withRef = kwargs.pop('withRef', False)
        self._spur = None
        self.Cm05 = None

        if isinstance(CM, str):
            self.load(CM)
        else:
            from pygimli.utils.geostatistics import covarianceMatrix

            if isinstance(CM, pgcore.Mesh):
                CM = covarianceMatrix(CM, **kwargs)
            elif CM is None and mesh is not None:
                CM = covarianceMatrix(mesh, **kwargs)
            else:
                pg.critical('Give either CM or mesh')

            self.Cm05 = createCm05(CM)

    def __repr__(self):
        """Return string representation of the matrix."""
        return "Geostatistical constraints matrix of size " + \
            f"{self.rows():d}x{self.cols():d} using {len(self.Cm05.ew):d} eigenvalues"

    @property
    def spur(self):
        """Return spur vector correcting for reference model effect."""
        if self._spur is None:
            if self.withRef is True:
                self._spur = np.zeros(self.rows())
            else:
                self._spur = self.Cm05 * pgcore.RVector(self.rows(), 1.0)
        return self._spur

    @property
    def nModel(self):
        """Return number of model parameters (size of Cm05)."""
        return self.Cm05.EV.shape[0] if hasattr(self.Cm05, "EV") else 0

    def save(self, fileName):
        """Save content of this matrix.

        Used for caching until pickling is possible for this class
        """
        self.Cm05.save(fileName + '-Cm05')
        np.save(fileName, dict(verbose=self.verbose(),
                               withRef=self.withRef,
                               Cm05=fileName + '-Cm05'),
                allow_pickle=True)

    def load(self, fileName):
        """Load the content of this matrix.

        Used for caching until pickling is possible for this class
        """
        d = np.load(fileName + '.npy', allow_pickle=True).tolist()
        self.setVerbose(d['verbose'], )
        self.withRef = d['withRef']
        self.Cm05 = Cm05Matrix(d['Cm05'])

    def mult(self, x):
        """Multiplication from right-hand side (A.T*x)."""
        return self.Cm05.mult(x) - self.spur * x

    def transMult(self, x):
        """Multiplication from left-hand side (A*x)."""
        return self.Cm05.transMult(x) - self.spur * x

    def cols(self):
        """Return number of columns (number of model parameters)."""
        return self.nModel

    def rows(self):
        """Return number of rows (number of constraints)."""
        return self.nModel

    def clear(self):
        """Clear cached matrices."""
        self.Cm05 = None
        self._spur = None


def hstack(mats):
    """Syntactic sugar function to horizontally stacked matrix.

    Parameters
    ----------
    mats : iterable
        list of matrices (any pyGIMLi matrix)

    Returns
    -------
    mat : :gimliapi:`GIMLI::Blockmatrix`
        pyGIMLi block matrix

    Examples
    --------
    >>> import pygimli as pg
    >>> A = pg.Matrix(4, 3)
    >>> B = pg.Matrix(4, 2)
    >>> H = pg.matrix.hstack([A, B])
    >>> ax, _ = pg.show(H)
    >>> print(H)
    pg.matrix.BlockMatrix of size 4 x 5 consisting of 2 submatrices.
    """
    assert not np.any(np.diff([m.rows() for m in mats])), \
        "Matrix row numbers do not match!"
    A = pgcore.BlockMatrix()
    icol = 0
    for mat in mats:
        A.addMatrix(mat, 0, icol)
        icol += mat.cols()

    A.recalcMatrixSize()
    return A


def vstack(mats):
    """Syntactic sugar function to vertically stacked matrix.

    Parameters
    ----------
    mats : iterable
        list of matrices (any pyGIMLi matrix)

    Returns
    -------
    mat : :gimliapi:`GIMLI::Blockmatrix`
        pyGIMLi block matrix

    Examples
    --------
    >>> import pygimli as pg
    >>> A = pg.Matrix(3, 4)
    >>> B = pg.Matrix(2, 4)
    >>> V = pg.matrix.vstack([A, B])
    >>> ax, _ = pg.show(V)
    >>> print(V)
    pg.matrix.BlockMatrix of size 5 x 4 consisting of 2 submatrices.
    """
    assert not np.any(np.diff([m.cols() for m in mats])), \
        "Matrix column numbers do not match!"
    A = pgcore.BlockMatrix()
    irow = 0
    for mat in mats:
        A.addMatrix(mat, irow, 0)
        irow += mat.rows()

    A.recalcMatrixSize()
    return A


def dstack(mats):
    """Syntactic sugar function to diagonally stacked matrix.

    Parameters
    ----------
    mats : iterable
        list of matrices (any pyGIMLi matrix)

    Returns
    -------
    mat : :gimliapi:`GIMLI::Blockmatrix`
        pyGIMLi block matrix

    Examples
    --------
    >>> import pygimli as pg
    >>> A = pg.Matrix(2, 3)
    >>> B = pg.Matrix(4, 5)
    >>> D = pg.matrix.dstack([A, B])
    >>> ax, _ = pg.show(D)
    >>> print(D)
    pg.matrix.BlockMatrix of size 6 x 8 consisting of 2 submatrices.
    """
    A = pgcore.BlockMatrix()
    irow = 0
    icol = 0
    for mat in mats:
        A.addMatrix(mat, irow, icol)
        irow += mat.rows()
        icol += mat.cols()

    A.recalcMatrixSize()
    return A


def matrixRow(A, n):
    """Return matrix row of arbitrary matrix."""
    assert n < A.rows(), "number exceeds number of rows in matrix"
    one = pg.Vector(A.rows())
    one[n] = 1.0
    return A.transMult(one)

def matrixColumn(A, n):
    """Return matrix column of arbitrary matrix."""
    assert n < A.cols(), "number exceeds number of columns in matrix"
    one = pg.Vector(A.cols())
    one[n] = 1.0
    return A.mult(one)

def complexMatrix(R, I, scaleR=1.0, scaleI=1.0):
    """Create a complex-valued matrix from two real-valued ones."""
    assert R.cols() == I.cols(), "Number of columns need to match"
    assert R.rows() == I.rows(), "Number of columns need to match"
    if isinstance(R, SparseMapMatrix) and isinstance(I, SparseMapMatrix):
        assert R.vecColPtr() == I.vecColPtr(), "sparsity structure differs"
        assert R.vecRowIdx() == I.vecRowIdx(), "sparsity structure differs"
        return pg.core.CSparseMatrix(R.vecColPtr(), R.vecRowIdx(),
                                     pg.core.toComplex(R.vecVals()*scaleR,
                                     I.vecVals()*scaleI))



if __name__ == "__main__":
    Amat = pg.Matrix(3, 4)
    Bmat = TransposedMatrix(Amat)
    xvec = pg.Vector(3, 1.0)
    print(Bmat*xvec)
    yvec = pg.Vector(4, 1.0)
    Cmat = SquaredMatrix(Amat)
    print(Cmat*yvec)
