#!/usr/bin/env python
"""Test the pygimli frameworks fitting routines."""

# write a correct test!
import unittest

import pygimli as pg
import numpy as np


class TestFrameworks(unittest.TestCase):
    """Test the framework fit function."""

    def test_Fit(self):
        """Test fitting of polynomial and exponential functions."""
        func = lambda x, a, b, c: a + b * x + c * x**2
        x = np.linspace(0, 1, 11)

        model = [1.5, 2, 2.5]
        data = func(x, *model)
        mEst, response = pg.frameworks.fit(func, data, x=x)
        np.testing.assert_allclose(mEst, model)
        np.testing.assert_allclose(data, response)

        model = [0.5, -1.0, 0.0]
        data = func(x, *model)  # + pg.randn(len(x))*0.001
        mEst, response = pg.frameworks.fit(func, data, x=x)
        # np.testing.assert_allclose(mEst, model, atol=1e-7)
        # np.testing.assert_allclose(data, response)
        # pg.plt.plot(x, data)
        # pg.plt.plot(x, response, 'x')
        print(model)

        data[5] = 1e-17
        # np.testing.assert_allclose(data2, data, atol=1e-15)
        mEst, response = pg.frameworks.fit(func, data, x=x)
        # np.testing.assert_allclose(mEst, model, atol=1e-8)
        # np.testing.assert_allclose(data, response)
        # pg.plt.plot(x, response, 'o')
        # print(mEst)
        # pg.wait()

        func = lambda t, a, b: a * np.exp(b * t)
        t = np.linspace(1, 2, 100)
        data = func(t, 1.1, 2.2)
        model, response = pg.frameworks.fit(func, data, t=t)

        # pg.plt.plot(t, data, '.', label='data')
        # pg.plt.plot(t, response, label='response')
        # pg.wait()

        np.testing.assert_allclose(model, [1.1, 2.2])
        np.testing.assert_allclose(data, response)


if __name__ == "__main__":
    # test = TestFrameworks()
    # test.test_Fit()

    unittest.main()
