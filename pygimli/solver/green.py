#!/usr/bin/env python

"""Greens function for some analytical stuff."""

import numpy as np


def greenDiffusion1D(x, t=0, a=1, dim=1):
    r"""Greens function for diffusion operator.

    Provides the elementary solution for:

    .. math::
        g(x,t) = \partial t + a \Delta

    To find a solution for:

    .. math::
        u(x,t) \quad\text{for}\quad\frac{\partial u(x,t)}{\partial t} +
        a \Delta u(x,t) = f(x)

    .. math:: x = [-x, 0, x]

    .. math:: u(x,t) = g(x,t) * f(x)

    Parameters
    ----------
    x : array_like
        Discrete spatial coordinates. Should better be symmetric [-x, 0, x].
    t : float
        Time
    a : float, optional
        Scale for Laplacian operator
    dim : int, optional
        Spatial dimension [1]

    Returns
    -------
    g : array_like
        Discrete Greens'function

    Examples
    --------
    >>> import numpy as np
    >>> from pygimli.solver import greenDiffusion1D
    >>> dx = 0.001
    >>> x = np.arange(0, 0.1+dx, dx)
    >>> g = greenDiffusion1D(np.hstack((-x[:0:-1], x)), t=1.0, a=1e-4)
    >>> g *= dx
    >>> f = np.zeros(len(x))
    >>> f[int(len(x)/2)] = 1.
    >>> u = np.convolve(g, f)[len(x)-1:2*len(x)-1]

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> _ = ax.plot(x, u, label="u(x)")
    >>> _ = ax.plot(x, g[::2], label="g(x)")
    >>> _ = ax.legend()
    >>> fig.show()
    """
    return 1. / (4. * np.pi * a * t)**(dim / 2.0) * \
            np.exp(-(x**2) / (4. * a * t))
