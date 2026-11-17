#!/usr/bin/env python
"""Test the pygimli iterator reference handling."""

# write a correct test!
import unittest

import pygimli as pg


class TestIterBug(unittest.TestCase):
    """Test that iteration keeps the referenced objects alive."""

    def test_MissingRefCounter(self):
        """Test slicing a vector and iterating over the result."""
        a = pg.Vector(10, 1)

        # this fails because of a missing reference counter: the iterator uses
        # the temporary object a(0, 9) whose reference count is not increased
        # I think our pygimli objects need a reference counter
        # when they are created from Python
        # print((a(0, 9).beginPyIter()[0], "!=", a[0]))
        # this might goes wrong but we cannot test it
        # self.assertNotEqual(a(0, 9).beginPyIter()[0], a[0])

        # this works because we define our own iterator that keeps the
        # reference
        # see: class VectorIter: in pygimli.__init__.py
        i = 0
        for ai in a[0:9]:
            self.assertEqual(ai, a[i])
            i = i + 1


if __name__ == "__main__":
    unittest.main()
