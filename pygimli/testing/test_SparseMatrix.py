#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestSparseMatrix(unittest.TestCase):

    def test_Convert(self):
        """
        """
        colIds = range(10)
        rowIds = range(10)
        vals = np.ones(10)

        # Construct SparseMap Matrix from python arrays
        A = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)

        # Construct SparseMap -> CRS (compressed row storage)
        S = pg.matrix.SparseMatrix(A)

        # Construct CRS -> SparseMap
        A2 = pg.matrix.SparseMapMatrix(S)

        # all should by identity matrix
        np.testing.assert_equal(A2.getVal(1, 1), 1.0)
        np.testing.assert_equal(sum(S * np.ones(S.cols())), S.rows())
        np.testing.assert_equal(sum(A2 * np.ones(A2.cols())), A2.rows())

        MAP1 = pg.matrix.SparseMapMatrix(r=3, c=15)
        CSR = pg.matrix.SparseMatrix(MAP1)
        MAP2 = pg.matrix.SparseMapMatrix(CSR)

        v3 = pg.Vector(3)
        v15 = pg.Vector(15)

        np.testing.assert_equal((MAP1*v15).size(), 3)
        np.testing.assert_equal((MAP1.transMult(v3)).size(), 15)

        np.testing.assert_equal((CSR*v15).size(), 3)
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

        #pg.solver.showSparseMatrix(mm, full=True)

        check_csr_rows = [0, 1, 2, 3, 4]
        check_csr_colPtr = [0, 2, 3, 4, 5]

        check_csc_cols = [0, 0, 1, 2, 3]
        check_csc_rowptr = [0, 1, 2, 3, 4, 5]

        r1, c1, v1 = pg.utils.sparseMatrix2Array(mm)
        np.testing.assert_allclose(r1, check_csr_rows)
        np.testing.assert_allclose(c1, check_csr_colPtr)
        np.testing.assert_allclose(v1, check_vals)

        sciA1 = pg.utils.sparseMatrix2csr(pg.matrix.SparseMatrix(mm))
        np.testing.assert_equal(sciA1.indices, check_csr_rows)
        np.testing.assert_equal(sciA1.indptr, check_csr_colPtr)

        sciA1 = pg.utils.sparseMatrix2csr(mm)
        np.testing.assert_equal(sciA1.indices, check_csr_rows)
        np.testing.assert_equal(sciA1.indptr, check_csr_colPtr)

        r2, c2, v2 = pg.utils.sparseMatrix2Array(pg.matrix.SparseMatrix(mm),
                                                 getInCRS=False)
        np.testing.assert_allclose(r2, check_rows)
        np.testing.assert_allclose(c2, check_cols)
        np.testing.assert_allclose(v2, check_vals)

        A1 = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)
        A2 = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)
        A1 += A2

        sciA1 = pg.utils.sparseMatrix2csr(pg.matrix.SparseMatrix(mm))
        sciA2 = pg.utils.sparseMatrix2csr(mm)
        np.testing.assert_equal(len(sciA1.data), mm.size())
        np.testing.assert_equal(sciA1.data, sciA2.data)
        np.testing.assert_equal(sciA1.indices, sciA2.indices)
        np.testing.assert_equal(sciA1.indptr, sciA2.indptr)

        sciA1 = pg.utils.sparseMatrix2coo(pg.matrix.SparseMatrix(mm))
        sciA2 = pg.utils.sparseMatrix2coo(mm)
        np.testing.assert_equal(len(sciA1.data), mm.size())
        np.testing.assert_equal(sciA1.data, sciA2.data)
        np.testing.assert_equal(sciA1.row, sciA2.row)
        np.testing.assert_equal(sciA1.col, sciA2.col)

        ### toSparseMatrix

        sciCSR = pg.utils.sparseMatrix2csr(pg.matrix.SparseMatrix(mm))
        np.testing.assert_equal(pg.utils.toSparseMatrix(sciCSR) == mm, True)



    def test_Access(self):
        #addVal(0, 1, 1.2) kommt nach der konvertierung auch wieder [0], [1], [1.2]
        pass

    def test_Operators(self):
        colIds = range(10)
        rowIds = range(10)
        vals = np.ones(10)
        A = pg.matrix.SparseMapMatrix(colIds, rowIds, vals)
        S = pg.matrix.SparseMatrix(A)

        S2 = S + S * 0.1 * 0.3

    def test_ComplexMatrix(self):
        verbose = False
        grid = pg.createGrid(3, 3)
        # print(grid)

        alpha = pg.math.toComplex(np.ones(grid.cellCount()),
                                  np.ones(grid.cellCount())*1.0)

        A = pg.solver.createStiffnessMatrix(grid, a=alpha)
        pg.solver.solver.applyDirichlet(A, None, [0], [0.0])
        #pg.solver.showSparseMatrix(A)
        #pg.solver.assembleDirichletBC(A, [[grid.boundary(0), 0.0]])

        b = pg.math.toComplex(np.ones(A.rows()), np.ones(A.rows())*0.0)
        x = pg.solver.linSolve(A, b, verbose=verbose, solver='pg')
        np.testing.assert_allclose(A.mult(x), b, rtol=1e-10)

        x2 = pg.solver.linSolve(A, b, verbose=verbose, solver='scipy')
        np.testing.assert_allclose(x2, x, rtol=1e-10)

        x3 = pg.solver.linSolve(pg.utils.squeezeComplex(A),
                                pg.utils.squeezeComplex(b),
                                verbose=verbose, solver='pg')

        np.testing.assert_allclose(pg.utils.toComplex(x3), x, rtol=1e-10)


    def test_BlockMatrix(self):
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
        print(B)


    def test_Misc(self):
        D = pg.SparseMapMatrix(3, 4)
        for i in range(D.rows()):
            for j in range(D.cols()):
                D.setVal(i, j, 1.0)

        np.testing.assert_allclose(D.col(2), pg.Vector(D.rows(), 1.0))
        np.testing.assert_allclose(D.row(2), pg.Vector(D.cols(), 1.0))

        D.cleanRow(1)
        np.testing.assert_allclose(D.col(2), [1.0, 0.0, 1.0])

        D.cleanCol(1)
        np.testing.assert_allclose(D.row(2), [1.0, 0.0, 1.0, 1.0])




if __name__ == '__main__':
    unittest.main()
