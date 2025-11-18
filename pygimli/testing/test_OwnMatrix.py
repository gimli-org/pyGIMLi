#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg

from pygimli.core.matrix import RVector, MatrixBase

class Matrix(pg.core.MatrixBase):

    def __init__(self):
        pg.MatrixBase.__init__(self)

    def rows(self): return 1

    def cols(self): return 1

    def mult(self, b):
        ret = pg.Vector(self.rows())
        print("TestMatrix::mult")
        return ret

    def transMult(self, b):
        ret = pg.Vector(self.cols())
        print("TestMatrix::transMult")
        return ret

    def save(self, name):
        print("TestMatrix::save", name)


class Modelling(pg.core.ModellingBase):

    def __init__(self):
        pg.core.ModellingBase.__init__(self, True)
        self.regionManager().setParameterCount(1)
        self.mat = Matrix()

    def response(self, model):
        print("TestModelling::response")
        res = pg.Vector(1, 1.0)
        return res

    def jacobian(self):
        print("TestModelling::jacobian()")
        return self.mat

    def createJacobian(self, model):
        print("TestModelling::createJacobian")


class TestOwnMatrix(unittest.TestCase):

    def test_runFOP(self):
        #F = TestModelling()

        #dat = pg.Vector(1, 1)
        #err = pg.Vector(1, 0.00001)

        #inv = pg.Inversion(dat, F, True, True)
        #inv.setError(err)

        #inv.run()
        pass

if __name__ == '__main__':

    unittest.main()
