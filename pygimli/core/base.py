# -*- coding: utf-8 -*-
"""
Misc stuff also needed for core imports and monkey patching
"""
import numpy as np

from .core import (RVector3, R3Vector, RMatrix)

def isInt(v, val=None):
    """Check if v is int , i.e. int, np.int32, np.int64.

    Examples
    --------
    >>> import pygimli as pg
    >>> print(pg.isInt(0))
    True
    >>> print(pg.isInt(np.int32(1)))
    True
    >>> print(pg.isInt(1.0))
    False
    """
    if val is None:
        return isinstance(v, (int, np.int32, np.int64))
    return isInt(v) and v == val

def isScalar(v, val=None):
    """Check if v is scalar, i.e. int, float or complex.

    Optional compare with val.

    Examples
    --------
    >>> import pygimli as pg
    >>> print(pg.isScalar(0))
    True
    >>> print(pg.isScalar(1.0))
    True
    >>> print(pg.isScalar(1.0, 0.0))
    False
    >>> print(pg.isScalar(1.0, 1.0))
    True
    >>> print(pg.isScalar(1+1j))
    True
    >>> print(pg.isScalar([0.0, 1.0]))
    False
    """
    if val is None:
        return isinstance(v, (int, float, complex, np.integer, np.complex128))
    # maybe add some tolerance check
    return isinstance(v, (int, float, complex, np.integer, np.complex128)) and v == val


def isIterable(v, N=None):
    """Check if `v` is iterable with optional size `N`.

    Examples
    --------
    >>> import pygimli as pg
    >>> import numpy as np
    >>> print(pg.isIterable([0, 1]))
    True
    >>> print(pg.isIterable([pg.Vector(5), pg.Vector(2)], N=2))
    True
    """
    if N is None:
        return hasattr(v, '__iter__')

    return isIterable(v) and len(v) == N


def isArray(v, N=None):
    """Check if `v` is a 1D array or a vector, with optional size `N`.

    Examples
    --------
    >>> import pygimli as pg
    >>> import numpy as np
    >>> print(pg.isArray([0, 1]))
    True
    >>> print(pg.isArray(np.ones(5)))
    True
    >>> print(pg.isArray(pg.Vector(5)))
    True
    >>> print(pg.isArray(pg.Vector(5), N=5))
    True
    >>> print(pg.isArray(pg.Vector(5), N=2))
    False
    >>> print(pg.isArray('foo'))
    False
    """
    if N is None:
        
        if isinstance(v, (tuple, list)):
            return isScalar(v[0])
            
        return (hasattr(v, '__iter__') and \
            not isinstance(v, (str))) and v.ndim == 1

    return isArray(v) and len(v) == N

def isComplex(vals):
    """Check numpy or pg.Vector if have complex data type"""
    if isScalar(vals):
        if isinstance(vals, (np.complex128, complex)):
            return True
    elif isArray(vals):
        return isComplex(vals[0])
    return False

def isPos(v):
    """Check if v is an array of size(3), [x,y,z], or pg.Pos.

    Examples
    --------
    >>> import pygimli as pg
    >>> print(pg.isPos([0.0, 0.0, 1.]))
    True
    >>> print(pg.isPos(pg.Pos(0.0, 0.0, 0.0)))
    True
    >>> print(pg.isPos(np.ones(3)))
    True
    >>> print(pg.isPos(np.ones(4)))
    False
    """
    return isArray(v, 2) or isArray(v, 3) or isinstance(v, RVector3)

def isR3Array(v, N=None):
    """Check if v is an array of size(N,3), a R3Vector or a list of pg.Pos.

    Examples
    --------
    >>> import pygimli as pg
    >>> print(pg.isR3Array([[0.0, 0.0, 1.], [1.0, 0.0, 1.]]))
    True
    >>> print(pg.isR3Array(np.ones((33, 3)), N=33))
    True
    >>> print(pg.isR3Array(pg.meshtools.createGrid(2,2).positions()))
    True
    """
    if N is None:
        return isinstance(v, R3Vector) or \
               (    isinstance(v, list) and isPos(v[0])) or \
               (not isinstance(v, list) and hasattr(v, '__iter__') and \
                not isinstance(v, (str)) and v.ndim == 2 and isPos(v[0]))
    return isR3Array(v) and len(v) == N

isPosList = isR3Array

def isMatrix(v, shape=None):
    """Check is v has ndim=2 or is comparable list"""
    if shape is None:
        return isinstance(v, RMatrix) or \
               hasattr(v, 'ndim') and v.ndim == 2 or \
                isinstance(v, list) and isArray(v[0])
    return isMatrix(v) and (hasattr(v, 'shape') and v.shape == shape)
