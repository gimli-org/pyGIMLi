#!/usr/bin/env python
"""Test user-defined matrix and modelling classes."""

# write a correct test!
import unittest

import pygimli as pg


class Matrix(pg.core.MatrixBase):
    """Provide a minimal matrix implementation for testing."""

    def __init__(self):
        pg.MatrixBase.__init__(self)

    def rows(self):
        """Return the number of rows."""
        return 1

    def cols(self):
        """Return the number of columns."""
        return 1

    def mult(self, b):
        """Multiply the matrix with the vector b."""
        ret = pg.Vector(self.rows())
        print("TestMatrix::mult")
        return ret

    def transMult(self, b):
        """Multiply the transposed matrix with the vector b."""
        ret = pg.Vector(self.cols())
        print("TestMatrix::transMult")
        return ret

    def save(self, name):
        """Write the matrix to the file given by name."""
        print("TestMatrix::save", name)


class Modelling(pg.core.ModellingBase):
    """Provide a minimal modelling operator for testing."""

    def __init__(self):
        pg.core.ModellingBase.__init__(self, True)
        self.regionManager().setParameterCount(1)
        self.mat = Matrix()

    def response(self, model):
        """Return the model response for the given model."""
        print("TestModelling::response")
        res = pg.Vector(1, 1.0)
        return res

    def jacobian(self):
        """Return the Jacobian matrix."""
        print("TestModelling::jacobian()")
        return self.mat

    def createJacobian(self, model):
        """Create the Jacobian matrix for the given model."""
        print("TestModelling::createJacobian")


class TestOwnMatrix(unittest.TestCase):
    """Test user-defined matrices and modelling classes."""

    def test_runFOP(self):
        """Run a forward operator using the custom matrix."""
        # F = TestModelling()

        # dat = pg.Vector(1, 1)
        # err = pg.Vector(1, 0.00001)

        # inv = pg.Inversion(dat, F, True, True)
        # inv.setError(err)

        # inv.run()
        pass


if __name__ == "__main__":
    unittest.main()
