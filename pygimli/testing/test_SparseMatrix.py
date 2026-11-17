#!/usr/bin/env python
"""Test the sparse matrix classes and their conversions."""

# write a correct test!
import unittest

import pygimli as pg
import numpy as np


class TestSparseMatrix(unittest.TestCase):
    """Test the sparse matrix classes and their conversions."""

    def test_Misc(self):
        """Test miscellaneous previously failing operations."""
        D = pg.matrix.asSparseMapMatrix(np.ones((3, 4)))

        np.testing.assert_allclose(D.col(2), pg.Vector(D.rows(), 1.0))
        np.testing.assert_allclose(D.row(2), pg.Vector(D.cols(), 1.0))

        D.cleanRow(1)
        np.testing.assert_allclose(D.col(2), [1.0, 0.0, 1.0])

        D.cleanCol(1)
        np.testing.assert_allclose(D.row(2), [1.0, 0.0, 1.0, 1.0])

        D = pg.matrix.asSparseMatrix(D)
        np.testing.assert_allclose(D.values().array(), np.ones(6))
        np.testing.assert_allclose(D.vecRowIdx().array(), [0, 2, 3, 0, 2, 3])

        D = pg.SparseMapMatrix(3, 4)
        for i in range(D.rows()):
            for j in range(D.cols()):
                D.setVal(i, j, 1.0)

    def test_Convert(self):
        """Test Sparse matrix conversion from * to *."""

        def _cmp_(A, B, timings=False):

            At = pg.matrix.asCSR(A)
            if timings:
                pg.toc("cmp A:", reset=True)
            Bt = pg.matrix.asCSR(B)
            if timings:
                pg.toc("cmp B:", reset=True)
            # np.testing.assert_equal(pg.utils.toSparseMatrix(A)
            #                         == pg.utils.toSparseMatrix(B), True)

            # A.indptr, A.indices, A.data
            # B.indptr, B.indices, B.data

            np.testing.assert_equal(At.indptr, Bt.indptr)
            np.testing.assert_equal(At.indices, Bt.indices)
            np.testing.assert_allclose(At.data, Bt.data)
            np.testing.assert_equal(np.sum(np.abs(At.data)) > 0, True)

        def _test_(A, timings=False):
            """Test all conversions for matrix A with timings."""
            if timings:
                print(f"A.shape: {A.shape}")

            pg.tic()
            T = pg.matrix.asSparseMatrix(A)
            if timings:
                pg.toc(f"{type(A)} -> RSparseMatrix:", reset=True)
            _cmp_(A, T, timings)

            pg.tic()
            T = pg.matrix.asSparseMapMatrix(A)
            if timings:
                pg.toc(f"{type(A)} -> RSparseMapMatrix:", reset=True)
            _cmp_(A, T, timings)

            pg.tic()
            T = pg.matrix.asCOO(A)
            if timings:
                pg.toc(f"{type(A)} -> COO", reset=True)
            _cmp_(A, T, timings)

            pg.tic()
            T = pg.matrix.asCSR(A)
            if timings:
                pg.toc(f"{type(A)} -> CSR", reset=True)
            _cmp_(A, T, timings)

            pg.tic()
            T = pg.matrix.asCSC(A)
            if timings:
                pg.toc(f"{type(A)} -> CSC", reset=True)
            _cmp_(A, T, timings)

        grid = pg.createGrid(20, 20, 20)

        # alpha = pg.math.toComplex(np.ones(grid.cellCount()),
        #                           np.ones(grid.cellCount())*1.0)

        A = pg.solver.createStiffnessMatrix(
            grid, a=np.ones(grid.cellCount()) * 3.14
        )

        _test_(A, timings=False)
        _test_(pg.matrix.asSparseMapMatrix(A), timings=False)
        _test_(pg.matrix.asCOO(A), timings=False)
        _test_(pg.matrix.asCSR(A), timings=False)

        colIds = range(10)
        rowIds = range(10)
        vals = np.ones(10)

        # Construct SparseMap Matrix from python arrays
        A = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)

        # SparseMap -> CRS (compressed row storage)
        S = pg.matrix.asSparseMatrix(A)

        # CRS -> SparseMap (COO)
        A2 = pg.matrix.asSparseMapMatrix(S)

        # all should be identity matrix
        np.testing.assert_equal(A2.getVal(1, 1), 1.0)
        np.testing.assert_equal(sum(S * np.ones(S.cols())), S.rows())
        np.testing.assert_equal(sum(A2 * np.ones(A2.cols())), A2.rows())

        MAP1 = pg.matrix.SparseMapMatrix(r=3, c=15)
        CSR = pg.matrix.asSparseMatrix(MAP1)
        MAP2 = pg.matrix.asSparseMapMatrix(CSR)

        v3 = pg.Vector(3)
        v15 = pg.Vector(15)

        np.testing.assert_equal((MAP1 * v15).size(), 3)
        np.testing.assert_equal((MAP1.transMult(v3)).size(), 15)

        np.testing.assert_equal((CSR * v15).size(), 3)
        np.testing.assert_equal((CSR.transMult(v3)).size(), 15)

        np.testing.assert_equal(MAP1.cols(), MAP2.cols())
        np.testing.assert_equal(CSR.cols(), MAP1.cols())
        np.testing.assert_equal(CSR.rows(), MAP1.rows())
        np.testing.assert_equal(MAP1.rows(), MAP2.rows())

        # testing SparseMatrix to Numpy
        mm = pg.matrix.SparseMapMatrix(r=4, c=5)
        check_rows = [0, 0, 1, 2, 3]
        check_cols = [0, 1, 2, 3, 4]
        check_vals = np.array([1.0, 3, np.pi, 1e-12, -1.12345e13])

        for i in range(len(check_rows)):
            mm.addVal(check_rows[i], check_cols[i], check_vals[i])

        # pg.solver.showSparseMatrix(mm, full=True)

        check_csr_rows = [0, 1, 2, 3, 4]
        check_csr_colPtr = [0, 2, 3, 4, 5]

        r1, c1, v1 = pg.matrix.sparseMatrix2Array(mm)
        np.testing.assert_allclose(r1, check_csr_rows)
        np.testing.assert_allclose(c1, check_csr_colPtr)
        np.testing.assert_allclose(v1, check_vals)

        sciA1 = pg.matrix.asCSR(pg.matrix.asSparseMatrix(mm))
        np.testing.assert_equal(sciA1.indices, check_csr_rows)
        np.testing.assert_equal(sciA1.indptr, check_csr_colPtr)

        sciA1 = pg.matrix.asCSR(mm)
        np.testing.assert_equal(sciA1.indices, check_csr_rows)
        np.testing.assert_equal(sciA1.indptr, check_csr_colPtr)

        r2, c2, v2 = pg.matrix.sparseMatrix2Array(
            pg.matrix.asSparseMatrix(mm), getInCRS=False
        )
        np.testing.assert_allclose(r2, check_rows)
        np.testing.assert_allclose(c2, check_cols)
        np.testing.assert_allclose(v2, check_vals)

        A1 = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)
        A2 = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)
        A1 += A2

    def test_Access(self):
        """Test matrix access (currently a placeholder)."""
        # addVal(0, 1, 1.2)
        pass

    def test_Operators(self):
        """Test the arithmetic operators of the sparse matrices."""

        def _test(A):
            R = pg.matrix.asCOO(A)

            ## needed? useful?
            # np.testing.assert_equal((A+2).values(), A.values()+2)
            # np.testing.assert_equal((A-2).values(), A.values()-2)
            # np.testing.assert_equal((A+2.).values(), A.values()+2.)
            # np.testing.assert_equal((A-2.).values(), A.values()-2.)
            # np.testing.assert_equal((2.+A).values(), 2+A.values())
            # np.testing.assert_equal((2.-A).values(), 2-A.values())
            # np.testing.assert_equal((2.+A).values(), 2.+A.values())
            # np.testing.assert_equal((2.-A).values(), 2.-A.values())

            np.testing.assert_equal(A * 2.0 == R * 2.0, True)
            np.testing.assert_equal(A / 2.0 == R / 2.0, True)
            np.testing.assert_equal(2.0 * A == 2.0 * R, True)
            np.testing.assert_equal(A * 2 == R * 2, True)
            np.testing.assert_equal(A / 2 == R / 2, True)
            np.testing.assert_equal(2 * A == 2 * R, True)
            np.testing.assert_equal(A * np.float64(0.1) == R * 0.1, True)

            np.testing.assert_equal(-A == -R, True)
            np.testing.assert_equal(A + A == R + R, True)
            np.testing.assert_equal(A + A * 2 == R + R * 2, True)
            np.testing.assert_equal(A - A * 2 == R - R * 2, True)
            np.testing.assert_equal(-A + A * 2 == -R + R * 2, True)

        A = np.linspace(0, 15, num=16).reshape(4, 4)
        _test(pg.matrix.asSparseMapMatrix(A))
        _test(pg.matrix.asSparseMatrix(A))

    def test_ComplexMatrix(self):
        """Test complex-valued sparse matrices (currently deactivated)."""
        verbose = False
        grid = pg.createGrid(3, 3)
        # print(grid)

        alpha = pg.math.toComplex(
            np.ones(grid.cellCount()), np.ones(grid.cellCount()) * 1.0
        )

        A = pg.solver.createStiffnessMatrix(grid, a=alpha)

        return  # test deactivated until we need complex matrices again
        pg.solver.solver.applyDirichlet(A, None, [0], [0.0])
        # pg.solver.showSparseMatrix(A)
        # pg.solver.assembleDirichletBC(A, [[grid.boundary(0), 0.0]])

        b = pg.math.toComplex(np.ones(A.rows()), np.ones(A.rows()) * 0.0)
        x = pg.solver.linSolve(A, b, verbose=verbose, solver="pg")
        np.testing.assert_allclose(A.mult(x), b, rtol=1e-10)

        x2 = pg.solver.linSolve(A, b, verbose=verbose, solver="scipy")
        np.testing.assert_allclose(x2, x, rtol=1e-10)

        # result not used here, kept until complex matrices are needed again
        pg.core.SparseMatrix(
            A.vecColPtr(), A.vecRowIdx(), pg.core.real(A.vecVals()), A.stype()
        )
        x3 = pg.solver.linSolve(
            pg.utils.squeezeComplex(A),
            pg.utils.squeezeComplex(b),
            verbose=verbose,
            solver="pg",
        )

        np.testing.assert_allclose(pg.utils.toComplex(x3), x, rtol=1e-10)

    def test_BlockMatrix(self):
        """Test the block matrix assembly."""
        A = pg.SparseMapMatrix(2, 2)
        A.setVal(0, 0, 1.0)

        B = pg.BlockMatrix()
        B.add(A, 0, 0)

        np.testing.assert_allclose(B.row(0), [1.0, 0.0], rtol=1e-10)
        B.add(A, 0, 0)
        np.testing.assert_allclose(B.row(0), [2.0, 0.0], rtol=1e-10)

        C = B.sparseMapMatrix()
        np.testing.assert_allclose(C.row(0), [2.0, 0.0], rtol=1e-10)

        B.add(A, 10, 10)

    def test_ReduceEntries(self):
        """Test the reduction of entries in a sparse matrix."""
        A = np.ones((5, 5))

        A1 = pg.matrix.asSparseMapMatrix(A)
        A2 = pg.matrix.asSparseMapMatrix(A)
        idx = [1, 2]
        pg.matrix.reduceEntries(A1, idx)

        for _i, ix in enumerate(idx):
            A2.cleanRow(ix)
            A2.cleanCol(ix)
            A2.setVal(ix, ix, 1.0)

        np.testing.assert_equal(A1 == A2, True)

    def test_Masks(self):
        """Test mask operations on sparse matrices."""
        A = np.linspace(1, 25, num=25).reshape(5, 5)
        A1 = pg.matrix.asSparseMatrix(A)

        diag = A1.createDiagonalMask()
        np.testing.assert_equal(A1.values(diag), [1, 7, 13, 19, 25])
        A1.setMaskValues(diag, 2.0)

        ## all diagonal values set to 2.0
        np.testing.assert_equal(A1.values(diag), [2.0] * A1.size())

        ## all values in row 1 and 3 set to 3.0
        mask = A1.createReduceMask([1, 3])
        ## check if mask is sorted and unique
        np.testing.assert_equal(mask, np.unique(mask))

        A1.setMaskValues(mask, 0.0)
        np.testing.assert_equal(sum(A1.col(1)), 0)
        np.testing.assert_equal(sum(A1.row(3)), 0)

        maskDiag = mask[np.nonzero(np.isin(mask, diag))[0]]
        A1.setMaskValues(maskDiag, 1.0)

        # manual reduce
        A2 = pg.matrix.asSparseMatrix(A)
        A2.setMaskValues(diag, 2.0)
        pg.matrix.reduceEntries(A2, [1, 3])
        A2.setVal(1, 1, 1.0)
        A2.setVal(3, 3, 1.0)

        ## check if mask reduce and manual reduce are equal
        np.testing.assert_equal(A1 == A2, True)


if __name__ == "__main__":
    unittest.main()
