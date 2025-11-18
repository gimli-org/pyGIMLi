#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

import numpy as np

import pygimli as pg


class TestRVectorMethods(unittest.TestCase):

    def test_XVectorBasics(self):

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
        """
        """
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
        # no clue why this doesnt work .. we might could hack them if someone need it
        #self.assertEqual(sum(v * 2), 10)

        self.assertEqual(sum(v + v), 10)
        self.assertEqual(sum(v * v), 5)
        self.assertEqual(sum(v - v), 0)
        self.assertEqual(sum(v / v), 5)
        self.assertEqual(sum(2 * v), 10)

    def test_IndexArray(self):
        v = pg.core.IndexArray([0,1,2,3])
        self.assertEqual(sum(v), 6)
        np.testing.assert_array_equal(v + 1, [1, 2, 3, 4])

    def test_RVectorOP(self):
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
        v = pg.Vector(5, 2.0)
        self.assertEqual(sum(pg.math.pow(v, 2)), 20)
        self.assertEqual(sum(pg.math.pow(v, 2.0)), 20)
        self.assertEqual(sum(pg.math.pow(v, v)), 20)

    def test_R3VectorIndex(self):
        r3 = pg.core.R3Vector(10)

        self.assertEqual(r3[0], pg.RVector3(0, 0, 0))
        np.testing.assert_array_equal(r3[0], pg.RVector3(0, 0, 0))

        r3[1] = pg.RVector3(0.0, 1.0, 0.0)
        np.testing.assert_array_equal(r3[1], pg.RVector3(0.0, 1.0, 0.0))

        r3[2] = (0.0, 2.0, 0.0)
        np.testing.assert_array_equal(r3[2], pg.RVector3(0.0, 2.0, 0.0))

        r3[3] = (0.0, 3.0, 0.0)
        np.testing.assert_array_equal(r3[3], pg.RVector3(0.0, 3.0, 0.0))

        d = pg.utils.dist(r3)
        self.assertEqual(sum(d), 1+2+3)

    def test_Slices(self):

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

    def test_IndexAccess(self):
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
        s = pg.core.RVector(range(5))
        np.testing.assert_array_equal(s, s[range(5)])
                
        s = pg.IVector(range(5))
        np.testing.assert_array_equal(s, s[range(5)])
        np.testing.assert_array_equal(s, s[s])

        s = pg.core.IndexArray(range(5))
        np.testing.assert_array_equal(s, s[range(5)])
        np.testing.assert_array_equal(s, s[s])

    def testComparison(self):
        a = pg.Vector(10, 1)
        b = pg.Vector(10, 2)

        np.testing.assert_equal(len(a < 1), 10)
        np.testing.assert_equal(len(a > 2), 10)

        np.testing.assert_equal(len(a < b), 10)
        np.testing.assert_equal(len(a > b), 10)

    def testUFunc(self):
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
        np.testing.assert_equal(np.int64(2) + M, 2+np.array(M))
        np.testing.assert_equal(M + np.int64(2), 2+np.array(M))


    def testRMatrixIndex(self):
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

if __name__ == '__main__':

    # pg.setDeepDebug(1)
    # t = TestRVectorMethods()

    # # t.test_IVectorOP()
    # t.test_Slices()
    # t.testRMatrixIndex()

    unittest.main()
