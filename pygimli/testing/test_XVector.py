#!/usr/bin/env python3
"""Test for core.Vector<*> conversion and operators."""

import unittest
import numpy as np
import pygimli as pg


class TestRVectorMethods(unittest.TestCase):
    """Vector methods."""

    def test_XVectorBasics(self):
        """Test basic operators for Vector, CVector, IVector."""

        def testVector(v):

            for t in v:
                self.assertEqual(t, 1.0)

            self.assertEqual(sum(v), 5.0)
            self.assertFalse(pg.core.haveInfNaN(v))

            v[1] = 0
            self.assertEqual(v[1], 0)

            v[1] = 1
            self.assertEqual(v[1], 1)
            #print(v/v)

        testVector(pg.Vector(5, 1.0))
        testVector(pg.CVector(5, 1.0))
        testVector(pg.BVector(5, True))
        testVector(pg.IVector(5, 1))


    def test_XVectorSetVal(self):
        """Test setVal for Vector, CVector, IVector."""
        #
        vec = pg.Vector(5, 1.0)
        self.assertEqual(vec[0], 1.0)

        vec[0] += 1.0
        self.assertEqual(vec[0], 2.0)

        vec[0] = vec[0] + 1.0
        self.assertEqual(vec[0], 3.0)

        vec.setVal(4.0, 0)
        self.assertEqual(vec[0], 4.0)

        vec.setVal(5.0, 0, 1)
        self.assertEqual(vec[0], 5.0)

        vec.setVal(6.0)
        self.assertEqual(vec[0], 6.0)

        vec.setVal(7.0, vec > 0.0)
        self.assertEqual(vec[0], 7.0)


    def test_IVectorOP(self):
        """Test operators for IVector."""
        v = pg.IVector(5, 1)

        # print(v + 2)
        # print(v - 2)
        # print(2 * v)
        #print(v * 2)
        #self.assertEqual(sum(v * 2), 10)

        self.assertEqual(sum(v + 1), 10)
        self.assertEqual(sum(v - 2), -5)
        self.assertEqual(sum(v / 1), 5)
        self.assertEqual(sum(1 + v), 10)
        self.assertEqual(sum(-1 - v), -10)
        self.assertEqual(sum(1 / v), 5)
        # no clue why this doesn't work .. we might could hack them if
        # someone need it
        #self.assertEqual(sum(v * 2), 10)

        self.assertEqual(sum(v + v), 10)
        self.assertEqual(sum(v * v), 5)
        self.assertEqual(sum(v - v), 0)
        self.assertEqual(sum(v / v), 5)
        self.assertEqual(sum(2 * v), 10)


    def test_IndexArray(self):
        """Test IndexArray."""
        v = pg.core.IndexArray([0,1,2,3])
        self.assertEqual(sum(v), 6)
        np.testing.assert_array_equal(v + 1, [1, 2, 3, 4])


    def test_RVectorOP(self):
        """Test operators for Vector."""
        v = pg.Vector(5, 1.0)

        self.assertEqual(sum(v + 1), 10)
        self.assertEqual(sum(v - 2), -5)
        self.assertEqual(sum(v * 2), 10)
        self.assertEqual(sum(v / 1), 5)
        self.assertEqual(sum(1 + v), 10)
        self.assertEqual(sum(-1 - v), -10)
        self.assertEqual(sum(2 * v), 10)
        self.assertEqual(sum(1 / v), 5)
        self.assertEqual(sum(v + 1.0), 10)
        self.assertEqual(sum(v - 2.0), -5)
        self.assertEqual(sum(v * 2.0), 10)
        self.assertEqual(sum(v / 1.0), 5)
        self.assertEqual(sum(1.0 + v), 10)
        self.assertEqual(sum(-1.0 - v), -10)
        self.assertEqual(sum(2.0 * v), 10)
        self.assertEqual(sum(1.0 / v), 5)

        v2 = np.ones(len(v))* 0.01
        # check pg * np
        self.assertEqual(sum(v * v2), 5*0.01)
        # check np * pg
        self.assertEqual(sum(v2 * v), 5*0.01)

        v = pg.Vector(3, -1.0)
        np.testing.assert_array_equal(v, [-1]*len(v))
        np.testing.assert_array_equal(abs(v), [1]*len(v))

        #v = pg.CVector(5, 1.0)

        #self.assertEqual(sum(v + 1), 10)
        #self.assertEqual(sum(v - 1), 0)
        #self.assertEqual(sum(v * 1), 5)
        #self.assertEqual(sum(v / 1), 5)
        #self.assertEqual(sum(1 + v), 10)
        #self.assertEqual(sum(1 - v), 0)
        #self.assertEqual(sum(1 * v), 5)
        #self.assertEqual(sum(1 / v), 5)
        #self.assertEqual(sum(v + 1.0), 10)
        #self.assertEqual(sum(v - 1.0), 0)
        #self.assertEqual(sum(v * 1.0), 5)
        #self.assertEqual(sum(v / 1.0), 5)
        #self.assertEqual(sum(1.0 + v), 10)
        #self.assertEqual(sum(1.0 - v), 0)
        #self.assertEqual(sum(1.0 * v), 5)
        #self.assertEqual(sum(1.0 / v), 5)


    def test_RVectorIndexRW(self):
        """Test index access for Vector."""
        v = pg.Vector(5, 2.0)
        np.testing.assert_array_equal(v, [2, 2, 2, 2, 2])

        v += 1.0
        np.testing.assert_array_equal(v, [3, 3, 3, 3, 3])

        v += 1
        np.testing.assert_array_equal(v, [4, 4, 4, 4, 4])

        v[1] = 1.0
        np.testing.assert_array_equal(v, [4, 1, 4, 4, 4])

        v[1] += 1.0
        np.testing.assert_array_equal(v, [4, 2, 4, 4, 4])

        v[[1,2]] = 2.0
        np.testing.assert_array_equal(v, [4, 2, 2, 4, 4])

        v[pg.IVector(1,3)] = 3.0
        np.testing.assert_array_equal(v, [4, 2, 2, 3, 4])

        v[pg.IVector(5,2)] = 1.0
        np.testing.assert_array_equal(v, [4, 2, 1, 3, 4])

        v[pg.find(v==4.0)] = 5.0
        np.testing.assert_array_equal(v, [5, 2, 1, 3, 5])

        v[v==5.0] = 4.0
        np.testing.assert_array_equal(v, [4, 2, 1, 3, 4])

        v[v==4.0] = 5.0
        np.testing.assert_array_equal(v, [5, 2, 1, 3, 5])

        #this will work only if we overwrite __iadd__
        #v[v==4.0] += 1.0
        #np.testing.assert_array_equal(v, [6, 2, 1, 3, 6])

        v.setVal(1.0, 1)
        np.testing.assert_array_equal(v, [5, 1, 1, 3, 5])


    def test_RVectorFuncts(self):
        """Test universal functions for Vector."""
        v = pg.Vector(5, 2.0)
        self.assertEqual(sum(pg.math.pow(v, 2)), 20)
        self.assertEqual(sum(pg.math.pow(v, 2.0)), 20)
        self.assertEqual(sum(pg.math.pow(v, v)), 20)


    def test_Pos(self):
        """Test operators for Pos."""
        p = pg.Pos(1.0, 0.0)

        np.testing.assert_array_equal(p, [1.0, 0.0, 0.0])
        np.testing.assert_array_equal(p + p, [2.0, 0.0, 0.0])
        np.testing.assert_array_equal(p + [0., 1, .0], [1.0, 1.0, 0.0])
        np.testing.assert_array_equal(p + [0, 1, 0], [1.0, 1.0, 0.0])
        np.testing.assert_array_equal([0., 1, 0] + p, [1.0, 1.0, 0.0])
        np.testing.assert_array_equal([0, 1, 0] + p, [1.0, 1.0, 0.0])
        np.testing.assert_array_equal([0., 1, 0] - p, [-1.0, 1.0, 0.0])
        np.testing.assert_array_equal([0, 1, 0] - p, [-1.0, 1.0, 0.0])
        np.testing.assert_array_equal(p + p/2, p + p/2.)
        np.testing.assert_array_equal(p + p+2, p + p+2.)
        np.testing.assert_array_equal(p + p*2, p + p*2.)
        np.testing.assert_array_equal(p + p-2, p + p-2.)

        p = pg.Pos(1.0, 1.0, 1.0)
        #np.testing.assert_equal(abs(p), p.abs()) # not implemented


    def test_PosVectorOP(self):
        """Test operators for [Pos, ...]."""
        r3 = pg.core.PosVector(10) + 2.0

        np.testing.assert_array_equal(r3, [pg.Pos(2, 2, 2)]*10)
        np.testing.assert_array_equal(r3 * 2.0, [pg.Pos(4, 4, 4)]*10)
        np.testing.assert_array_equal(r3 * np.long(2), [pg.Pos(4, 4, 4)]*10)

        #TODO native RValue suuport for int->double conversion to avoid this hack
        #pg.core.setDeepDebug(1)
        #r3*2
        #halt

        np.testing.assert_array_equal(r3 * 2, [pg.Pos(4, 4, 4)]*10)
        np.testing.assert_array_equal(r3 / 2.0, [pg.Pos(1.0, 1.0, 1.0)]*10)
        np.testing.assert_array_equal(1 / r3, [pg.Pos(0.5, 0.5, 0.5)]*10)
        np.testing.assert_array_equal(-1 - r3, [pg.Pos(-3.0, -3.0, -3.0)]*10)
        np.testing.assert_array_equal(abs(r3 / abs(r3)), [1.0]*10)


    def test_stdVectorRVectorOP(self):
        """Test operators for sr=[RVector, ...].
            BINOP(sr, float) -> sr
            BINOP(float, sr) -> sr
            BINOP(np.array, sr) -> sr
            BINOP(sr, np.array) -> sr
        """
        sr = pg.core.stdVectorRVector()
        sr.append(pg.core.RVector(5))
        sr.append(pg.core.RVector(10))

        for i, v in enumerate(sr):
            np.testing.assert_array_equal(v, [0]*len(v))

        sr += 1.0
        for i, v in enumerate(sr):
            np.testing.assert_array_equal(v, [1]*len(v))

        sr -= 2
        for i, v in enumerate(sr):
            np.testing.assert_array_equal(v, [-1]*len(v))

        sr *= -2
        for i, v in enumerate(sr):
            np.testing.assert_array_equal(v, [2]*len(v))

        sr /= 4.0
        for i, v in enumerate(sr):
            np.testing.assert_array_equal(v, [0.5]*len(v))

        for i, v in enumerate(sr*3):
            np.testing.assert_array_equal(v, sr[i]*3)

        for i, v in enumerate(1/sr):
            np.testing.assert_array_equal(v, 1/sr[i])

        for i, v in enumerate(sr*abs(-1*sr)):
            np.testing.assert_array_equal(v, sr[i]**2)

        for i, v in enumerate(sr + 1.0):
            np.testing.assert_array_equal(v, sr[i]+1)

        #
        ###  BINOP(ndarray, sr) -> sr
        #
        a = np.full(len(sr), 2.0)

        for i, v in enumerate(sr + a):
            np.testing.assert_array_equal(v, sr[i] + a[i])

        for i, v in enumerate(a + sr):
            np.testing.assert_array_equal(v, a[i] + sr[i])

        for i, v in enumerate(sr - a):
            np.testing.assert_array_equal(v, sr[i] - a[i])

        for i, v in enumerate(a - sr):
            np.testing.assert_array_equal(v, a[i] - sr[i])

        for i, v in enumerate(sr * a):
            np.testing.assert_array_equal(v, sr[i] * a[i])

        for i, v in enumerate(a * sr):
            np.testing.assert_array_equal(v, a[i] * sr[i])

        for i, v in enumerate(sr / a):
            np.testing.assert_array_equal(v, sr[i] / a[i])

        for i, v in enumerate(a / sr):
            np.testing.assert_array_equal(v, a[i] / sr[i])

        #
        ### unary operators
        #
        for i, v in enumerate(-sr):
            np.testing.assert_array_equal(v, -sr[i])

        for i, v in enumerate(sr**-1):
            np.testing.assert_array_equal(v, sr[i]**-1)

        for i, v in enumerate(sr**a):
            np.testing.assert_array_equal(v, sr[i]**a[i])



    def test_stdVectorPosVectorOP(self):
        """Test operators for [PosVector, ...]."""
        r3 = pg.core.stdVectorPosVector()
        r3.append(pg.core.PosVector(7) + 2.0)
        r3.append(2 + pg.core.PosVector(4))
        r3.append(3.15 + pg.core.PosVector(10))
        r3.append(pg.core.PosVector(9) + 4)

        t = abs(r3) ## perfcheck candidate
        np.testing.assert_equal(len(t), len(r3))

        for i in range(len(t)):
            np.testing.assert_array_equal(t[i], abs(r3[i]))

        t = r3/abs(r3)  ## perfcheck candidate
        np.testing.assert_equal(len(t), len(r3))

        for i in range(len(t)):
            np.testing.assert_array_equal(abs(t[i]), abs(r3[i]/abs(r3[i])))

        t = t * 2.0
        for i in range(len(t)):
            np.testing.assert_array_equal(abs(t[i]), [2]*len(r3[i]))


    def test_PosVectorIndex(self):
        """Test index access for [Pos, ...]."""
        r3 = pg.core.PosVector(10)

        self.assertEqual(r3[0], pg.Pos(0, 0, 0))
        np.testing.assert_array_equal(r3[0], pg.Pos(0, 0, 0))

        r3[1] = pg.Pos(0.0, 1.0, 0.0)
        np.testing.assert_array_equal(r3[1], pg.Pos(0.0, 1.0, 0.0))

        r3[2] = (0.0, 2.0, 0.0)
        np.testing.assert_array_equal(r3[2], pg.Pos(0.0, 2.0, 0.0))

        r3[3] = (0.0, 3.0, 0.0)
        np.testing.assert_array_equal(r3[3], pg.Pos(0.0, 3.0, 0.0))

        d = pg.utils.dist(r3)
        self.assertEqual(sum(d), 1+2+3)


    def test_Slices(self):
        """Test slice access for Vector."""
        a = pg.Vector(np.arange(10.))

        np.testing.assert_array_equal(a[:], np.arange(10.)[:])
        np.testing.assert_array_equal(a[:, np.newaxis],
                                      np.array(a)[:, np.newaxis])
        np.testing.assert_array_equal(a[::], np.arange(10.)[::])
        np.testing.assert_array_equal(a[:-1], np.arange(10.)[:-1])
        np.testing.assert_array_equal(a[:-9], np.arange(10.)[:-9])

        np.testing.assert_array_equal(a[:-10], np.arange(10.)[:-10])

        np.testing.assert_array_equal(a[::1], np.arange(10.)[::1])
        np.testing.assert_array_equal(a[::-1], np.arange(10.)[::-1])

        np.testing.assert_array_equal(a[0:3:1], np.arange(10.)[0:3:1])
        np.testing.assert_array_equal(a[0:3:2], np.arange(10.)[0:3:2])

        np.testing.assert_array_equal(a[3:0:-1], np.arange(10.)[3:0:-1])
        np.testing.assert_array_equal(a[3:0:-2], np.arange(10.)[3:0:-2])

        np.testing.assert_array_equal(a[0:3:-1], np.arange(10.)[0:3:-1])
        np.testing.assert_array_equal(a[0:3:-2], np.arange(10.)[0:3:-2])

        np.testing.assert_array_equal(sum(a[:-1]), sum(a[:len(a)-1]))


    def test_IndexAccess(self):
        """Test index access for Vector."""
        # (double) array/vector
        an = np.arange(10.)
        ag = pg.Vector(an)

        # bn = nd.array(bool)
        bn = (an > 4.)
        self.assertEqual(type(bn), np.ndarray)
        self.assertEqual(bn.dtype, 'bool')
        self.assertEqual(sum(bn), 5)

        # bg = BVector
        bg = (ag > 4.)
        self.assertEqual(type(bg), pg.BVector)
        self.assertEqual(sum(bg), 5)

        # BVector(nd.array(bool))
        self.assertEqual(len(bg), len(pg.BVector(bn)))
        self.assertEqual(sum(bg), sum(pg.BVector(bn)))
        self.assertEqual(bg[0], pg.BVector(bn)[0])
        np.testing.assert_array_equal(bg, pg.BVector(bn))

        # In = nd.array(int)
        In = np.nonzero(bn)[0]
        self.assertEqual(type(In), np.ndarray)
        self.assertEqual(In.dtype, 'int64')
        self.assertEqual(len(In), 5)
        self.assertEqual(In[0], 5)

        # np.nonzero(bg)
        np.testing.assert_array_equal(In, np.nonzero(bg)[0])

        # Ig = IndexArray
        Ig = pg.find(bg)
        self.assertEqual(type(Ig), pg.core.IndexArray)
        self.assertEqual(len(Ig), 5)
        self.assertEqual(Ig[0], 5)

        # pg.find(nd.array(bool))
        np.testing.assert_array_equal(Ig, pg.find(bn))

        ## Indexoperators ##
        # ndarray [nd.array(bool)] == ndarray [nd.array(int)]
        np.testing.assert_equal(an[bn], an[In])
        self.assertEqual(len(an[bn]), 5)
        self.assertEqual(an[bn][0], 5)

        # ndarray[IndexArray] == ndarray [nd.array(int)]
        np.testing.assert_equal(an[Ig], an[In])

        # ndarray[BVector] == ndarray [nd.array(bool)]
        np.testing.assert_array_equal(an[np.array(bg, dtype='bool')], an[bn])
        np.testing.assert_array_equal(an[np.array(bg)], an[bn])
        np.testing.assert_array_equal(an[bg.array()], an[bn])
        np.testing.assert_array_equal(an[an>5], [6, 7, 8, 9])

        np.testing.assert_array_equal(ag[bg], ag[Ig])
        self.assertEqual(len(ag[bg]), 5)
        self.assertEqual(ag[bg][0], 5)

        # RVector [BVector] ==  RVector [nd.array(bool)]
        np.testing.assert_array_equal(ag[bg], ag[bn])
        np.testing.assert_equal(sum(ag[bg]), sum(ag[bn]))

        # RVector [IndexArray] ==  RVector [nd.array(int)]
        np.testing.assert_array_equal(ag[Ig], ag[In])
        # RVector(BVector) ==  RVector(nd.array(bool))
        # RVector(IndexArray) ==  RVector(nd.array(int))

        I = pg.core.IndexArray([0,1,1,0])
        np.testing.assert_array_equal(pg.sum(I), 2)
        np.testing.assert_array_equal(sum(I), 2)
        np.testing.assert_array_equal(np.sum(I), 2)


    def testIndexGetter(self):
        """Test index getter for Vector."""
        s = pg.core.RVector(range(5))
        np.testing.assert_array_equal(s, s[range(5)])

        s = pg.IVector(range(5))
        np.testing.assert_array_equal(s, s[range(5)])
        np.testing.assert_array_equal(s, s[s])

        s = pg.core.IndexArray(range(5))
        np.testing.assert_array_equal(s, s[range(5)])
        np.testing.assert_array_equal(s, s[s])


    def testComparison(self):
        """Test comparison operators for Vector."""
        a = pg.Vector(10, 1)
        b = pg.Vector(10, 2)

        np.testing.assert_array_equal(a < 1, [False]*10)
        np.testing.assert_array_equal(a == 1, [True]*10)
        np.testing.assert_array_equal(b > 1, [True]*10)
        np.testing.assert_array_equal(a < b, [True]*10)
        np.testing.assert_array_equal(a > b, [False]*10)
        np.testing.assert_equal(len(a > b), 10)

        a = pg.IVector(10, 1)
        np.testing.assert_array_equal(a == 1, [True]*10)
        np.testing.assert_array_equal(a != 1, [False]*10)
        np.testing.assert_array_equal(a < 1, [False]*10)

        # print(a < 1, np.array([False]*10))
        # np.testing.assert_equal(len(a < 1), 10)


    def testUFunc(self):
        """Test universal functions for Vector."""
        t = pg.Vector(5, 1)

        np.testing.assert_equal(t*2, np.array(t)*2.0)
        np.testing.assert_equal(2*t, np.array(t)*2.0)

        np.testing.assert_equal(t*np.float64(2), np.array(t)*2.0)
        np.testing.assert_equal(np.float64(2)*t, np.array(t)*2.0)
        np.testing.assert_equal(np.int64(2)*t, np.array(t)*2.0)
        np.testing.assert_equal(t*np.int64(2), np.array(t)*2.0)

        M = pg.Matrix(3, 3) + 1.0
        np.testing.assert_equal(np.float64(2)*M, 2*np.array(M))
        np.testing.assert_equal(M*np.float64(2), 2*np.array(M))
        np.testing.assert_equal(np.int64(2)*M, 2*np.array(M))
        np.testing.assert_equal(M*np.int64(2), 2*np.array(M))
        np.testing.assert_equal(np.int64(2) + M, 2 + np.array(M))
        np.testing.assert_equal(M + np.int64(2), 2 + np.array(M))
        np.testing.assert_equal(np.int64(2) - M, 2 - np.array(M))
        np.testing.assert_equal(M - np.int64(2), np.array(M)-2)


class TestMatrixMethods(unittest.TestCase):
    """Test matrix methods."""

    def testSparseMatrixBasics(self):
        """Test basic operators for SparseMapMatrix."""
        A = pg.SparseMapMatrix(3,3)
        for i in range(A.size()):
            A.setVal(i,i, 1.0)

        np.testing.assert_equal(A + 2.0 * A ==
                                A + A * 2.0, True)

        ## check to ignore np.*.__mul__ in case A has __rmul__ function
        np.testing.assert_equal(A + np.float64(2.0) * A ==
                                A + A * np.float64(2.0), True)


    def testRMatrixIndex(self):
        """Test index access for Matrix."""
        A = pg.Matrix(3,4)
        A[0] = pg.Vector(4,1)

        np.testing.assert_equal(sum(A[0]), 4)

        A[1,2] = 2.0
        # np.testing.assert_equal(sum(A[1]), 2)
        np.testing.assert_equal(A[1,2], 2)

        np.testing.assert_equal(A[:,2], A.col(2))
        np.testing.assert_equal(A[2], A.row(2))

        ## will not work because A[2] refer to A[2]__getItem__ which only can
        # return a const reference. use the tuple idx above
        # A[2][2] = 2.0
        # np.testing.assert_equal(sum(A[2]), 2)


    def testRDensMatrixAccess(self):
        """Test access for RDenseMatrix."""
        #A = pg.core.RMatrix(2,2)
        #A = pg.core.RDenseMatrix(2,2)
        # print(type(A[0] + 1.0))
        # print(A[0] + 1)
        # print(A[0] + 1)
        # for v in A[0] + 1:
        #     print(v)
        # for v in A[0] + 1:
        #     print(v)


    def testRDensMatrixRefCounting(self):
        """Test reference counting for RDenseMatrix."""
        def createE():
            E = pg.core.ElementMatrix()
            E.setIds([0, 1], [0,0])
            E.multR = 1.0
            E.resize(2, 2)
            A = pg.core.RMatrix(2,2)
            A = pg.core.RDenseMatrix(A+1.0)
            E.setMat(A)
            return E

            ## for debugging only
            # E = pg.core.TestEM()
            # A = pg.core.RDenseMatrix(2, 2)
            # A[0].fill(1.0)
            # E.setMat(A)
            # return E

        E = createE()
        np.testing.assert_equal(E.mat().row(0), [1.0, 1.0])

        ### Without GC patch () this fails since mat() is copy
        # by reference, and E is destructed after createE().mat()
        # has been called so mat is not more valid
        np.testing.assert_equal(createE().mat().row(0), [1.0, 1.0])


if __name__ == '__main__':

    # pg.setDeepDebug(1)
    # t = TestRVectorMethods()

    # # t.test_IVectorOP()
    # t.test_Slices()
    # t.testRMatrixIndex()

    unittest.main()
