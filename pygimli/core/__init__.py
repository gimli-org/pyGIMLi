# -*- coding: utf-8 -*-
"""
Imports and extensions of the C++ bindings.
"""
import sys
import numpy as np

from .core import pgcore
from .core import *

from .logger import error, critical
from .base import (isInt, isScalar, isArray, isPos,
                   isPosList)

# #######################################
# ###  Global convenience functions #####
# #######################################

pgcore.load = None

###########################
# print function for GIMLi stuff
############################


def __RVector_str(self, valsOnly=False):
    """String representation of RVector."""
    s = str()

    if not valsOnly:
        s = str(self.size())

    if len(self) == 0:
        return s
    else:
        s += " ["

    if len(self) < 101:
        for i in range(0, len(self) - 1):
            s = s + str(self[i]) + ", "

        s = s + str(self[len(self) - 1]) + "]"
        return s
    return (str(self.size()) + " [" + str(self[0]) + ",...," + str(
        self[self.size() - 1]) + "]")


def __Pos_str(self):
    """String representation of Pos."""
    return ("Pos: (" + str(self.x()) + ", " + str(self.y()) + ", " + str(
        self.z()) + ")")


def __R3Vector_str(self):
    """String representation of R3Vector."""
    if self.size() < 20:
        return str(self.array())

    return "R3Vector: n=" + str(self.size())


def __Line_str(self):
    """String representation of Line."""
    return "Line: " + str(self.p0()) + "  " + str(self.p1())


def __BoundingBox_str(self):
    """String representation of BoundingBox."""
    s = ''
    s += f"BoundingBox [{self.min()}, {self.max()}]"
    return s


pgcore.RVector.__repr__ = __RVector_str
pgcore.CVector.__repr__ = __RVector_str
pgcore.BVector.__repr__ = __RVector_str
pgcore.IVector.__repr__ = __RVector_str
pgcore.IndexArray.__repr__ = __RVector_str
pgcore.Pos.__repr__ = __Pos_str
pgcore.R3Vector.__repr__ = __R3Vector_str

pgcore.Line.__repr__ = __Line_str
pgcore.BoundingBox.__repr__ = __BoundingBox_str

############################
# compatibility stuff
############################


def nonzero_test(self):
    """Function to throw a warning if any vector is used as bool vector."""
    raise BaseException("Warning! there is no 'and' and 'or' for "
                        "BVector and RVector. " +
                        "Use binary operators '&' or '|' instead. " +
                        "If you looking for the nonzero test, use len(v) > 0")


def np_round__(self, r):
    """Make numpy.round also work for pg.Vector."""
    return np.round(self.array(), r)


pgcore.RVector.__bool__ = nonzero_test
pgcore.R3Vector.__bool__ = nonzero_test
pgcore.BVector.__bool__ = nonzero_test
pgcore.CVector.__bool__ = nonzero_test
pgcore.IVector.__bool__ = nonzero_test
pgcore.IndexArray.__bool__ = nonzero_test

pgcore.RVector.__nonzero__ = nonzero_test
pgcore.R3Vector.__nonzero__ = nonzero_test
pgcore.BVector.__nonzero__ = nonzero_test
pgcore.CVector.__nonzero__ = nonzero_test
pgcore.IVector.__nonzero__ = nonzero_test
pgcore.IndexArray.__nonzero__ = nonzero_test

pgcore.RVector.__round__ = np_round__


def _invertBVector_(self):
    """Invert (toggle) the BVector."""
    return pgcore.inv(self)
pgcore.BVector.__invert__ = _invertBVector_
pgcore.BVector.__inv__ = _invertBVector_


def _lowerThan_(self, v2):
    """Test whether any vector is lower than another.

    Overwrite BVector = v1 < v2 since there is a wrong operator due to the
    boost binding generation
    """
    if isinstance(v2, int):
        return self < float(v2)
    return pgcore.inv(self >= v2)
pgcore.RVector.__lt__ = _lowerThan_
pgcore.R3Vector.__lt__ = _lowerThan_
pgcore.BVector.__lt__ = _lowerThan_
pgcore.CVector.__lt__ = _lowerThan_
# pgcore.IVector.__lt__ = _lowerThan_
# pgcore.IndexArray.__lt__ = _lowerThan_


__pgcore_RVector___gt__ = pgcore.RVector.__gt__
def _greaterThan_(self, v2):
    """ Overwrite BVector = v1 > v2 since there is a wrong default conversion
    from v2 of type int -> v2 (RVector(v2))
    boost binding generation
    """
    if isinstance(v2, int):
        return self > float(v2)
    return __pgcore_RVector___gt__(self, v2)
pgcore.RVector.__gt__ = _greaterThan_

######################
# special constructors
######################

# Overwrite constructor for IndexArray
# This seams ugly but necessary until we can recognize numpy array in
# custom_rvalue
__origIndexArrayInit__ = pgcore.IndexArray.__init__
def __newIndexArrayInit__(self, arr=None, val=None):
    """Initialize Index array."""
    if arr is None:
        __origIndexArrayInit__(self)
    elif hasattr(arr, 'dtype') and hasattr(arr, '__iter__'):
        __origIndexArrayInit__(self, np.asarray(arr, dtype=np.uint64))
        #__origIndexArrayInit__(self, [int(a) for a in arr])
    else:
        if val:
            __origIndexArrayInit__(self, arr, val)
        else:
            __origIndexArrayInit__(self, arr)
pgcore.IndexArray.__init__ = __newIndexArrayInit__


# Overwrite constructor for BVector
# This seams ugly but necessary until we can recognize numpy array in
# custom_rvalue
__origBVectorInit__ = pgcore.BVector.__init__
def __newBVectorInit__(self, arr=None, val=None):
    """Initialize Bool array."""
    if arr is None:
        __origIndexArrayInit__(self)
    elif hasattr(arr, 'dtype') and hasattr(arr, '__iter__'):
        # this is hell slow .. better in custom_rvalue.cpp or in
        # vector.h directly from python object
        __origBVectorInit__(self, len(arr))
        for i, a in enumerate(arr):
            self.setVal(bool(a), i)
    else:
        if val:
            __origBVectorInit__(self, arr, val)
        else:
            __origBVectorInit__(self, arr)
pgcore.BVector.__init__ = __newBVectorInit__

######################
# special overwrites
######################

# RVector + int fails .. so we need to tweak this command
__oldRVectorAdd__ = pgcore.RVector.__add__
def __newRVectorAdd__(a, b):
    if isinstance(b, np.ndarray) and b.dtype == complex:
        return __oldRVectorAdd__(a, pgcore.CVector(b))
    if isInt(b):
        return __oldRVectorAdd__(a, float(b))
    if isInt(a):
        return __oldRVectorAdd__(float(a), b)
    return __oldRVectorAdd__(a, b)
pgcore.RVector.__add__ = __newRVectorAdd__


__oldRVectorSub__ = pgcore.RVector.__sub__
def __newRVectorSub__(a, b):
    if isInt(b):
        return __oldRVectorSub__(a, float(b))
    if isInt(a):
        return __oldRVectorSub__(float(a), b)
    return __oldRVectorSub__(a, b)
pgcore.RVector.__sub__ = __newRVectorSub__


__oldRVectorMul__ = pgcore.RVector.__mul__
def __newRVectorMul__(a, b):
    if isInt(b):
        return __oldRVectorMul__(a, float(b))
    if isInt(a):
        return __oldRVectorMul__(float(a), b)
    return __oldRVectorMul__(a, b)
pgcore.RVector.__mul__ = __newRVectorMul__

try:
    __oldRVectorTrueDiv__ = pgcore.RVector.__truediv__

    def __newRVectorTrueDiv__(a, b):
        if isInt(b):
            return __oldRVectorTrueDiv__(a, float(b))
        if isInt(a):
            return __oldRVectorTrueDiv__(float(a), b)
        return __oldRVectorTrueDiv__(a, b)

    pgcore.RVector.__truediv__ = __newRVectorTrueDiv__
except BaseException:
    __oldRVectorTrueDiv__ = pgcore.RVector.__div__

    def __newRVectorTrueDiv__(a, b):
        if isInt(b):
            return __oldRVectorTrueDiv__(a, float(b))
        if isInt(a):
            return __oldRVectorTrueDiv__(float(a), b)
        return __oldRVectorTrueDiv__(a, b)
    pgcore.RVector.__div__ = __newRVectorTrueDiv__


__oldRMatMul__ = pgcore.RMatrix.__mul__
def __newRMatMul__(a, b):
    if isInt(b):
        return __oldRMatMul__(a, float(b))
    return __oldRMatMul__(a, b)
pgcore.RMatrix.__mul__ = __newRMatMul__


__oldRMatAdd__ = pgcore.RMatrix.__add__
def __newRMatAdd__(a, b):
    if isInt(b):
        return __oldRMatAdd__(a, float(b))
    return __oldRMatAdd__(a, b)
pgcore.RMatrix.__add__ = __newRMatAdd__


###############################################################################
# override wrong default conversion from int to IndexArray(int) for setVal    #
###############################################################################
__origRVectorSetVal__ = pgcore.RVector.setVal
def __newRVectorSetVal__(self, *args, **kwargs):
    # print('__newRVectorSetVal__', *args, **kwargs)
    if len(args) == 2:
        if isinstance(args[1], int):
            if args[1] < 0:
                return __origRVectorSetVal__(self, args[0],
                                             i=len(self) + args[1])
            else:
                return __origRVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], pgcore.BVector):
            return __origRVectorSetVal__(self, args[0], bv=args[1])
    return __origRVectorSetVal__(self, *args, **kwargs)
pgcore.RVector.setVal = __newRVectorSetVal__


__origR3VectorSetVal__ = pgcore.R3Vector.setVal
def __newR3VectorSetVal__(self, *args, **kwargs):
    # print('__newRVectorSetVal__', *args, **kwargs)
    if len(args) == 2:
        if isinstance(args[1], int):
            return __origR3VectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], pgcore.BVector):
            return __origR3VectorSetVal__(self, args[0], bv=args[1])
    return __origR3VectorSetVal__(self, *args, **kwargs)
pgcore.R3Vector.setVal = __newR3VectorSetVal__


__origBVectorSetVal__ = pgcore.BVector.setVal
def __newBVectorSetVal__(self, *args, **kwargs):
    if len(args) == 2:
        if isinstance(args[1], int):
            return __origBVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], pgcore.BVector):
            return __origBVectorSetVal__(self, args[0], bv=args[1])
    return __origBVectorSetVal__(self, *args, **kwargs)
pgcore.BVector.setVal = __newBVectorSetVal__


__origCVectorSetVal__ = pgcore.CVector.setVal
def __newCVectorSetVal__(self, *args, **kwargs):
    if len(args) == 2:
        if isinstance(args[1], int):
            return __origCVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], pgcore.BVector):
            return __origCVectorSetVal__(self, args[0], bv=args[1])
    return __origCVectorSetVal__(self, *args, **kwargs)
pgcore.CVector.setVal = __newCVectorSetVal__


__origIVectorSetVal__ = pgcore.IVector.setVal
def __newIVectorSetVal__(self, *args, **kwargs):
    if len(args) == 2:
        if isinstance(args[1], int):
            return __origIVectorSetVal__(self, args[0], i=args[1])
        if isinstance(args[1], pgcore.BVector):
            return __origIVectorSetVal__(self, args[0], bv=args[1])
    return __origIVectorSetVal__(self, *args, **kwargs)
pgcore.IVector.setVal = __newIVectorSetVal__


__origIndexArraySetVal__ = pgcore.IndexArray.setVal
def __newIndexArraySetVal__(self, *args, **kwargs):
    if len(args) == 2:
        if isinstance(args[1], int):
            return __origIndexArraySetVal__(self, args[0], i=args[1])
        if isinstance(args[1], pgcore.BVector):
            return __origIndexArraySetVal__(self, args[0], bv=args[1])
    return __origIndexArraySetVal__(self, *args, **kwargs)
pgcore.IndexArray.setVal = __newIndexArraySetVal__


############################
# Indexing [] operator for RVector, CVector, IndexArray,
#                          Pos, R3Vector, RMatrix, CMatrix
############################
def __getVal(self, idx):
    """Get vector value at index. Hell slow."""
    if isinstance(idx, pgcore.BVector):
        return self.get_(idx)
    elif isinstance(idx, pgcore.IVector):
        return self.getVSI_(idx)
    elif isinstance(idx, pgcore.IndexArray):
        return self.getVUI_(idx)
    elif isinstance(idx, slice):

        s = idx.start
        e = idx.stop
        if s is None:
            s = 0
        if e is None:
            e = len(self)

        if idx.step is None:
            return self.getVal(int(s), int(e))
        else:
            step = idx.step
            if step < 0 and idx.start is None and idx.stop is None:
                ids = range(e - 1, s - 1, idx.step)
            else:
                ids = range(s, e, idx.step)

            if len(ids):
                return self.getVSI_(ids)

            return self.get_(0)

    elif isinstance(idx, list) or hasattr(idx, '__iter__'):
        if isinstance(idx[0], int):
            return self.getVSI_(idx)

        if hasattr(idx[0], 'dtype'):
            # print("numpy: ", idx[0].dtype.str, idx[0].dtype ,type(idx[0]))
            if idx[0].dtype == 'bool':
                return self.getVUI_([i for i, x in enumerate(idx) if x])
                # return self[np.nonzero(idx)[0]]
        elif isinstance(idx[0], slice):  # try fixing newaxis
            # probably the call x = x[:, np.newaxis]
            # so we return np.array here
            return np.array(self)[idx]
        # elif isinstance(idx[0], None) and isinstance(idx[1], slice):
            # return self[idx[1]]

        return self.getVSI_([int(a) for a in idx])

    elif idx < 0:
        idx = len(self) + idx

    return self.getVal(int(idx))


def __setVal(self, idx, val):
    """Index write access (x[i]=y)."""
    # print("__setVal", self, 'idx', idx, 'val:', val)
    if isinstance(idx, slice):
        if idx.step is None:
            if int(idx.stop) >= self.size():
                self.resize(idx.stop)

            if idx.start is None:
                self.setVal(val, 0, int(idx.stop))
            else:
                self.setVal(val, int(idx.start), int(idx.stop))
            return
        else:
            critical("not yet implemented for slice:", slice)
    elif isinstance(idx, tuple):
        # print("tuple", idx, type(idx))
        if isinstance(self, pgcore.RMatrix):
            self.rowRef(int(idx[0])).setVal(val, int(idx[1]))
            return
        else:
            error("Can't set index with tuple", idx, "for", self)
            return
    # if isinstance(idx, pgcore.BVector):
    # print("__setVal", self, idx, 'val:', val)
    # self.setVal(val, bv=idx)
    # return
    if isinstance(val, complex):
        if isinstance(idx, int):
            return self.setVal(val=val, id=idx)
        else:
            return self.setVal(val=val, ids=idx)

    if isinstance(self, pgcore.RMatrix):
        self.setVal(idx, val)
    else:
        self.setVal(val, idx)


def __getValR3(self, idx):
    """Get R3 value at index."""
    if isinstance(idx, slice):
        s = idx.start
        e = idx.stop
        if s is None or s < 0:
            s = 0
        if e is None or e > 3:
            e = 3

        if idx.step is not None:
            print('slice step not supported here', idx)
        return np.array(self)[s:e]

    return self.getVal(int(idx))


def __getValMatrix(self, idx):
    """Get matrix value at index."""
    #pg.error('in use?')
    #    print(idx, type(idx))
    if isinstance(idx, slice):
        step = idx.step
        if step is None:
            step = 1
        start = idx.start
        if start is None:
            start = 0
        stop = idx.stop
        if stop is None:
            stop = len(self)

        return [self.rowRef(i) for i in range(start, stop, step)]

    elif isinstance(idx, tuple):
        # print(idx, type(idx))
        if isinstance(idx[0], slice):
            if isinstance(idx[1], int):
                tmp = self[idx[0]]
                ret = pgcore.RVector(len(tmp))
                for i, t in enumerate(tmp):
                    ret[i] = t[idx[1]]
                return ret
        else:
            return self.row(int(idx[0]))[idx[1]]

    if idx == -1:
        idx = len(self) - 1

    return self.row(idx)

pgcore.RMatrix.__setitem__ = __setVal
pgcore.RDenseMatrix.__setitem__ = __setVal
pgcore.RMatrix.__getitem__ = __getValMatrix # in use? (slow)-- inline is better

pgcore.RVector.__setitem__ = __setVal
pgcore.RVector.__getitem__ = __getVal  # very slow -- inline is better

pgcore.CVector.__setitem__ = __setVal
pgcore.CVector.__getitem__ = __getVal  # very slow -- inline is better

pgcore.BVector.__setitem__ = __setVal
pgcore.BVector.__getitem__ = __getVal  # very slow -- inline is better

pgcore.IVector.__setitem__ = __setVal
pgcore.IVector.__getitem__ = __getVal  # very slow -- inline is better

pgcore.R3Vector.__setitem__ = __setVal
pgcore.R3Vector.__getitem__ = __getVal  # very slow -- inline is better

pgcore.IndexArray.__setitem__ = __setVal
pgcore.IndexArray.__getitem__ = __getVal  # very slow -- inline is better

pgcore.Pos.__setitem__ = __setVal
pgcore.Pos.__getitem__ = __getValR3 # support slice


# len(RVector), RMatrix
_vec = [pgcore.RVector,
         pgcore.BVector,
         pgcore.CVector,
         pgcore.IVector,
         pgcore.IndexArray]

for v in _vec:
    v.ndim = 1
    v.__len__ = lambda self: self.size()
    v.shape = property(lambda self: (self.size(),))
    # if hasattr(v, '__call__') and callable(getattr(v, '__call__')):
    try:
        del v.__call__
    except AttributeError:
        pass

pgcore.RVector.dtype = float
pgcore.BVector.dtype = bool
pgcore.CVector.dtype = complex
pgcore.IVector.dtype = int
pgcore.IndexArray.dtype = np.uint

pgcore.Pos.dtype = float
pgcore.Pos.__len__ = lambda self: 3
pgcore.Pos.ndim = 1
pgcore.Pos.shape = (3,)

pgcore.R3Vector.dtype = float
pgcore.R3Vector.__len__ = lambda self: self.size()
pgcore.R3Vector.ndim = 2
pgcore.R3Vector.shape = property(lambda self: (self.size(), 3))

# remove me 23_08_28
# pgcore.stdVectorPos.ndim = 2

############################
# abs(RVector), RMatrix
############################
pgcore.RVector.__abs__ = pgcore.fabs
pgcore.CVector.__abs__ = pgcore.mag
pgcore.R3Vector.__abs__ = pgcore.absR3

############################
# __hash__ settings
############################
pgcore.RVector.__hash__ = pgcore.RVector.hash
pgcore.CVector.__hash__ = pgcore.CVector.hash
pgcore.IVector.__hash__ = pgcore.IVector.hash
pgcore.IndexArray.__hash__ = pgcore.IndexArray.hash
pgcore.R3Vector.__hash__ = pgcore.R3Vector.hash
pgcore.Pos.__hash__ = pgcore.Pos.hash
pgcore.DataContainer.__hash__ = pgcore.DataContainer.hash
pgcore.DataContainerERT.__hash__ = pgcore.DataContainerERT.hash
pgcore.Mesh.__hash__ = pgcore.Mesh.hash


def _stdVectorR3Vector_hash_(self):
    """Hash function for stdVectorR3Vector."""
    h = 0
    for rv in self:
        h = h ^ rv.__hash__()
    return h
pgcore.stdVectorR3Vector.__hash__ = _stdVectorR3Vector_hash_

############################
# Iterator support for RVector allow to apply python build-ins
############################


class VectorIter:
    """Iterator for RVector, CVector, BVector, IVector, IndexArray."""
    def __init__(self, vec):
        """Initialize iterator."""
        self.it = vec.beginPyIter()
        self.vec = vec

    def __iter__(self):
        """Return iterator."""
        return self

    # this is for python < 3
    def next(self):
        """Return next element."""
        return self.it.nextForPy()

    # this is the same but for python > 3
    def __next__(self):
        """Return next element."""
        return self.it.nextForPy()


def __VectorIterCall__(self):
    """Return iterator."""
    return VectorIter(self)
    # don't use pygimli iterators here until the reference for temporary
    # vectors are collected
    # return pgcore.RVectorIter(self.beginPyIter())

pgcore.RVector.__iter__ = __VectorIterCall__
pgcore.R3Vector.__iter__ = __VectorIterCall__
pgcore.BVector.__iter__ = __VectorIterCall__
pgcore.IVector.__iter__ = __VectorIterCall__
pgcore.IndexArray.__iter__ = __VectorIterCall__
pgcore.CVector.__iter__ = __VectorIterCall__


class DefaultContainerIter:
    """Iterator for RMatrix, RDenseMatrix, CMatrix."""
    def __init__(self, vec):
        """Initialize iterator."""
        self.vec = vec
        self.length = len(vec)
        self.pos = -1

    def __iter__(self):
        """Return iterator."""
        return self

    def next(self):
        """Return next element."""
        return self.__next__()

    # this is the same but for python > 3
    def __next__(self):
        """Return next element."""
        self.pos += 1
        if self.pos == self.length:
            raise StopIteration()
        else:
            return self.vec[self.pos]

def __MatIterCall__(self):
    """Return iterator."""
    return DefaultContainerIter(self)

pgcore.RMatrix.__iter__ = __MatIterCall__
pgcore.RDenseMatrix.__iter__ = __MatIterCall__
pgcore.CMatrix.__iter__ = __MatIterCall__


class PosIter():
    """Simple iterator for Pos/PosVector.

    Because it lacks the core function .beginPyIter()
    """
    def __init__(self, vec):
        """Initialize iterator."""
        self.vec = vec
        self.length = 3
        self.pos = -1

    def __iter__(self):
        """Return iterator."""
        return self

    def next(self):
        """Return next element."""
        return self.__next__()

    def __next__(self):
        """Return next element."""
        self.pos += 1
        if self.pos == self.length:
            raise StopIteration()
        else:
            return self.vec[self.pos]

def __PosIterCall__(self):
    return PosIter(self)
pgcore.Pos.__iter__ = __PosIterCall__

## there are weird exposed value expressions for POS OP with int argument
## so we need to overwrite them here until its known how to avoid the exposure
## 231006


__POS_orig__add__ = getattr(pgcore.Pos, '__add__')
def __POS_new__add__(self, b):
    if isinstance(b, int):
        return __POS_orig__add__(self, float(b))
    return __POS_orig__add__(self, b)
setattr(pgcore.Pos, '__add__', __POS_new__add__)


__POS_orig__radd__ = getattr(pgcore.Pos, '__radd__')
def __POS_new__radd__(self, b):
    if isinstance(b, int):
        return __POS_orig__radd__(self, float(b))
    elif isinstance(b, list):
        return pgcore.Pos(b) + self
    else:
        return __POS_orig__radd__(self, b)
setattr(pgcore.Pos, '__radd__', __POS_new__radd__)


__POS_orig__sub__ = getattr(pgcore.Pos, '__sub__')
def __POS_new__sub__(self, b):
    if isinstance(b, int):
        return __POS_orig__sub__(self, float(b))
    return __POS_orig__sub__(self, b)
setattr(pgcore.Pos, '__sub__', __POS_new__sub__)


__POS_orig__rsub__ = getattr(pgcore.Pos, '__rsub__')
def __POS_new__rsub__(self, b):
    if isinstance(b, int):
        return __POS_orig__rsub__(self, float(b))
    elif isinstance(b, list):
        return pgcore.Pos(b) - self
    else:
        return __POS_orig__rsub__(self, b)
setattr(pgcore.Pos, '__rsub__', __POS_new__rsub__)


__POS_orig__mul__ = getattr(pgcore.Pos, '__mul__')
def __POS_new__mul__(self, b):
    if isinstance(b, int):
        return __POS_orig__mul__(self, float(b))
    return __POS_orig__mul__(self, b)
setattr(pgcore.Pos, '__mul__', __POS_new__mul__)


__POS_orig__truediv__ = getattr(pgcore.Pos, '__truediv__')
def __POS_new__truediv__(self, b):
    if isinstance(b, int):
        return __POS_orig__truediv__(self, float(b))
    return __POS_orig__truediv__(self, b)
setattr(pgcore.Pos, '__truediv__', __POS_new__truediv__)


# ######### c to python converter ######
# default converter from Pos to numpy array
def __PosArrayCall__(self, dtype=None, **kwargs):
    # if idx:
    # print(self)
    # print(idx)
    # raise Exception("we need to fix this")
    import numpy as np
    return np.array([self.getVal(0), self.getVal(1), self.getVal(2)])

# default converter from RVector to numpy array


def __RVectorArrayCall__(self, dtype=None, **kwargs):
    #
    # This wrapper is needed to catch explicit type conversion from np.asarray
    #
    #if idx and not isinstance(idx, numpy.dtype):
    #print("self:", self)
    #print("idx:", idx, type(idx) )
    #raise Exception("we need to fix this")
    # probably fixed!!!
    # import numpy as np
    # we need to copy the array until we can handle increasing the reference
    # counter in self.array() else it leads to strange behavior
    # test in testRValueConverter.py:testNumpyFromRVec()
    # return np.array(self.array())
    return self.array()


def __CVectorArrayCall__(self, dtype=None, **kwargs):
    #
    # This wrapper is needed to catch explicit type conversion from np.asarray
    #
    # #if idx and not isinstance(idx, numpy.dtype):
    #print("self:", self)
    #print("idx:", idx, type(idx) )
    #raise Exception("we need to fix this")
    # probably fixed!!! or not!!
    # import numpy as np
    # we need to copy the array until we can handle increasing the reference
    # counter in self.array() else it leads to strange behavior
    # test in testRValueConverter.py:testNumpyFromRVec()
    # return np.array(self.array())
    return self.array()


# default converter from RVector to numpy array
pgcore.RVector.__array__ = __RVectorArrayCall__
# not yet ready handmade_wrappers.py
pgcore.BVector.__array__ = __RVectorArrayCall__
pgcore.IndexArray.__array__ = __RVectorArrayCall__
pgcore.stdVectorI.__array__ = __RVectorArrayCall__
pgcore.R3Vector.__array__ = __RVectorArrayCall__

# should work with handmade wrapper -- test me!
pgcore.Pos.__array__ = __PosArrayCall__

pgcore.R3Vector.append = pgcore.R3Vector.push_back

# see bug description
pgcore.CVector.__array__ = __CVectorArrayCall__


##################################
# stdVectorRVector operators
# Check why they are not exposed by generator!!
##################################

__IOP__ = ['__iadd__',
           '__isub__',
           '__imul__',
           #'__itruediv__',
           ]

__BINOP__ = ['__add__',
             '__sub__',
             '__mul__',
             '__truediv__',
             '__radd__',
             '__rsub__',
             '__rmul__',
             '__rtruediv__',
             ]

def __stdVectorRVector_NEG__(self):
    """Unary minus operator."""
    ret = pgcore.stdVectorRVector()
    for i, ai in enumerate(self):
        ret.append(-ai)
    return ret
pgcore.stdVectorRVector.__neg__ = __stdVectorRVector_NEG__


def __stdVectorRVector_POW__(self, exp):
    """Power operator."""
    ret = pgcore.stdVectorRVector()
    for i, ai in enumerate(self):
        ret.append(ai**exp)
    return ret
pgcore.stdVectorRVector.__pow__ = __stdVectorRVector_POW__


def __stdVectorRVector_ABS__(self):
    """Absolute value operator."""
    ret = pgcore.stdVectorRVector()
    for i, ai in enumerate(self):
        ret.append(pgcore.fabs(ai))
    return ret
pgcore.stdVectorRVector.__abs__ = __stdVectorRVector_ABS__


def __stdVectorRVector_IOP__(self, a, OP):
    if isScalar(a):
        for i in range(len(self)):
            self[i] = getattr(self[i], OP)(a)
    elif len(self) == len(a):
        for i in range(len(self)):
            getattr(self[i], OP)(a)
    else:
        critical(f'Cannot {OP} stdVectorRVector with different lengths. {len(self)} != {len(a)}')

    return self

def __stdVectorRVector_BIOP__(self, b, OP):

    ret = pgcore.stdVectorRVector()
    if isScalar(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b))
    elif len(self) == len(b):
        for i, ai in enumerate(self):
            #print(i, ai, b[i], OP)
            ret.append(getattr(ai, OP)(b[i]))
    else:
        critical(f'Cannot {OP} stdVectorRVector with different lengths. {len(self)} != {len(b)}')

    return ret

# for _OP in __IOP__:
#     def _closure(OP):
#         return lambda a, b: __stdVectorRVector_IOP__(a, b, OP)
#     setattr(pgcore.stdVectorRVector, _OP, _closure(_OP))

for _OP in __BINOP__:
    def _closure(OP):
        return lambda a, b: __stdVectorRVector_BIOP__(a, b, OP)
    setattr(pgcore.stdVectorRVector, _OP, _closure(_OP))


def __stdVectorRVector__array_ufunc__(self, ufunc, method, *inputs, **kwargs):
    """
    For numpy operators
    """
    ###  more advanced behaviour possible ... instead of __array_ufunc__ = None
    ### https://docs.scipy.org/doc/numpy-1.13.0/neps/ufunc-overrides.html
    if ufunc == np.power:
        vRv = inputs[0]
        exp = inputs[1]
        # print('*'*80, vRv)
        # print('-'*80, exp)
        return vRv ** exp

    ### rmul
    if ufunc == np.multiply:
        if isinstance(inputs[0], pgcore.stdVectorRVector):
            return inputs[0] * inputs[1]
        else:
            return inputs[1] * inputs[0]

    ## default apply numpy function
    ret = pgcore.stdVectorRVector()
    for i, ai in enumerate(self):
        ret.append(ufunc(ai))
    return ret


    pg._r(self)
    pg._y(f'ufunc: {ufunc}')
    pg._y(f'method: {method}')
    pg._y(f'inputs: {inputs}')
    pg._y(f'kwargs: {kwargs}')
    pg.critical('implementme')

pgcore.stdVectorRVector.__array_ufunc__ = __stdVectorRVector__array_ufunc__
#pgcore.stdVectorRVector.__array_ufunc__ = None

##################################
# stdVectorR3Vector operators
##################################

def __stdVectorR3Vector_BINARY_OP__(self, b, OP):

    ret = pgcore.stdVectorR3Vector()

    if isScalar(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b))
        return ret

    elif len(self) == len(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b[i]))
    else:
        pg.critical(f'Cannot {OP} stdVectorR3Vector with different lengths:'
                    f'{len(self)}{len(b)}')

    return ret

def __stdVectorR3Vector_ABS__(self):
    ret = pgcore.stdVectorRVector()
    for i, ai in enumerate(self):
        ret.append(abs(ai))
    return ret

pgcore.stdVectorR3Vector.__abs__ = __stdVectorR3Vector_ABS__


for _OP in __BINOP__:
    def _closure(OP):
        return lambda a, b: __stdVectorR3Vector_BINARY_OP__(a, b, OP)
    setattr(pgcore.stdVectorR3Vector, _OP, _closure(_OP))


def __stdVectorRMatrix_BINARY_OP__(self, b, OP):

    ret = pgcore.stdVectorRMatrix()

    if isScalar(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b))
        return ret

    elif len(self) == len(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b[i]))
    else:
        pg.critical(f'Cannot {OP} stdVectorRMatrix with different lengths:'
                    f'{len(self)}{len(b)}')

    return ret


for _OP in __BINOP__:
    def _closure(OP):
        return lambda a, b: __stdVectorRMatrix_BINARY_OP__(a, b, OP)
    setattr(pgcore.stdVectorRMatrix, _OP, _closure(_OP))


def __stdVectorMatrixVector_BINARY_OP__(self, b, OP):

    ret = pgcore.stdVectorMatrixVector()

    if isScalar(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b))
        return ret

    elif len(self) == len(b):
        for i, ai in enumerate(self):
            ret.append(getattr(ai, OP)(b[i]))
    else:
        pg.critical(f'Cannot {OP} stdVectorMatrixVector with different lengths:'
                    f'{len(self)}{len(b)}')

    return ret


for _OP in __BINOP__:
    def _closure(OP):
        return lambda a, b: __stdVectorMatrixVector_BINARY_OP__(a, b, OP)
    setattr(pgcore.stdVectorMatrixVector, _OP, _closure(_OP))


def __RMatrix_NEG__(self):
    return pgcore.RMatrix(self*-1.0)
pgcore.RMatrix.__neg__ = __RMatrix_NEG__

def __stdVectorRMatrix_NEG__(self):
    ret = pgcore.stdVectorRMatrix()
    for i, ai in enumerate(self):
        ret.append(-ai)
    return ret
pgcore.stdVectorRMatrix.__neg__ = __stdVectorRMatrix_NEG__


def __stdVectorMatrixVector_NEG__(self):
    ret = pgcore.stdVectorMatrixVector()
    for i, ai in enumerate(self):
        ret.append(-ai)
    return ret
pgcore.stdVectorMatrixVector.__neg__ = __stdVectorMatrixVector_NEG__


##################################
# custom rvalues for special cases
##################################

def find(v):
    """Find a specific entry in vector."""
    if hasattr(v, 'dtype') and hasattr(v, '__iter__'):
        # print('new find', v, pgcore.BVector(v))
        return pgcore.find(pgcore.BVector(v))
    else:
        # print('orig find')
        return pgcore.find(v)


def pow(v, p):
    """Power function.

    pow(v, int) is misinterpreted as pow(v, rvec(int)), so we need to fix this
    """
    if isinstance(p, int):
        return pgcore.pow(v, float(p))
    return pgcore.pow(v, p)


def __RVectorPower(self, m):
    return pow(self, m)


pgcore.RVector.__pow__ = __RVectorPower

##################################
# usefull aliases
##################################

Vector = pgcore.RVector
Pos = pgcore.Pos
PosVector = pgcore.R3Vector
PosList = PosVector
Inversion = pgcore.RInversion


############################
# non automatic exposed functions
############################
__PY_ABS__ = abs

def abs(v):
    """Create abs in the sense of distance instead of just vanishing the sign.

    Create abs in the sense of distance instead of vanishing the sign. Used
    to calculate the length of coordinates, or anything that can be interpreted
    as coordinate.

    Args
    ----
    v: iterable of float, complex, or :gimliapi:`GIMLI::Pos`

    Returns
    -------
    length: iterable or scalar
        Array of lenghts.

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> pg.abs([1.0, 1.0, 1.0])
    1.7320508075688772
    >>> pg.abs(np.array([1.0, 1.0, 1.0]))
    1.7320508075688772
    >>> pg.abs(np.array([1.0, 1.0]))
    1.4142135623730951
    >>> pg.abs([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    2 [1.7320508075688772, 1.7320508075688772]
    >>> pg.abs(np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]))
    2 [1.7320508075688772, 1.7320508075688772]
    >>> # Note, this will be interpreted as 3 2Dim Pos
    >>> pg.abs(np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]).T)
    3 [1.4142135623730951, 1.4142135623730951, 1.4142135623730951]
    >>> pg.abs(pg.PosList([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]))
    2 [1.7320508075688772, 1.7320508075688772]
    """
    if isinstance(v, pgcore.CVector):
        return pgcore.mag(v)
    elif isPos(v):
        return pgcore.Pos(v).abs()
    elif isPosList(v):
        return pgcore.absR3(v)
    elif isinstance(v, list):
        # possible [x,y,[z]] or [pos, ...]
        try:
            return pgcore.Pos(v).abs()
        except:
            return pgcore.absR3(np.array(v).T)
    elif isinstance(v, pgcore.R3Vector):
        return pgcore.absR3(v)
    elif isinstance(v, np.ndarray):
        if v.ndim == 1:
            return np.abs(v)
        if v.shape[0] == 2 or v.shape[0] == 3:
            return pgcore.absR3(v.T)
        else:
            return pgcore.absR3(v)
    elif isinstance(v, pgcore.RMatrix):
        raise BaseException("IMPLEMENTME")
        for i in range(len(v)):
            v[i] = pgcore.abs(v[i])
        return v
    elif hasattr(v, 'values'):
        #import pygimli as pg
        return abs(v.values)
    elif hasattr(v, 'vals'):
        return abs(v.vals)
    elif isinstance(v, (int, float)):
        return __PY_ABS__(v)

    try:
        return __PY_ABS__(v)
    except:
        return pgcore.fabs(v)


# default BVector operator == (RVector, int) will be casted to
# BVector operator == (RVector, RVector(int)) and fails
# this needs a monkey patch for BVector operator == (RVector, int)
pgcore.__EQ_RVector__ = pgcore.RVector.__eq__


def __EQ_RVector__(self, val):
    if isinstance(val, int):
        val = float(val)
    return pgcore.__EQ_RVector__(self, val)


pgcore.RVector.__eq__ = __EQ_RVector__


############################
# usefull stuff
############################
def toIVector(v):
    print("do not use toIVector(v) use ndarray directly .. "
          "this method will be removed soon")
    ret = pgcore.IVector(len(v), 0)
    for i, r in enumerate(v):
        ret[i] = int(r)
    return ret


############################
# for pickleing
############################

def __ByteBuffer_setstate(self, state):
    """Recover bytestream data from state['_bytestream_']"""
    self.fill(state['_bytestream_'])

def __ByteBuffer_reduce(self):
    """Create bytestream data"""
    return (pgcore.ByteBuffer, (),
        dict(_bytestream_=self.array()))

pgcore.ByteBuffer.__setstate__ = __ByteBuffer_setstate
pgcore.ByteBuffer.__reduce__ = __ByteBuffer_reduce


def __FEAFunction_reduce(self):
    return (pgcore.FEAFunction, (self.valueSize(), self.getEvalOrder()),
                    )
pgcore.FEAFunction.__reduce__ = __FEAFunction_reduce


def __Pos_reduce(self):
    return (pgcore.Pos, (self.x(), self.y(), self.z()),
                    )
pgcore.Pos.__reduce__ = __Pos_reduce


def __enablePickle(cls):
    """Enable pickling for classes that support serialization.
    """
    def _setstate(self, state):
        self.deserialize(state['_bytestream_'])
    def _reduce(self):
        return (cls, (),
                dict(_bytestream_=self.serialize()))

    cls.__setstate__ = _setstate
    cls.__reduce__ = _reduce

__enablePickle(pgcore.Mesh)
__enablePickle(pgcore.RVector)
__enablePickle(pgcore.IVector)
__enablePickle(pgcore.BVector)
__enablePickle(pgcore.IndexArray)



# DEPRECATED for backward compatibility should be removed
def asvector(array):
    """Convert numpy array into vector (not to be used anymore!)."""
    print("do not use asvector(ndarray) use ndarray directly .. "
          "this method will be removed soon")
    return pgcore.RVector(array)


# ##########################
# We want ModellingBase with multi threading jacobian brute force
# ##########################


def __GLOBAL__response_mt_shm_(fop, model, shm, i):
    resp = fop.response_mt(model, i)

    for j in range(len(resp)):
        shm[j] = resp[j]


def __ModellingBase__createJacobian_mt__(self, model, resp):
    from math import ceil
    from multiprocessing import Process, Array
    import numpy as np

    nModel = len(model)
    # nData = len(resp)  # not used

    fak = 1.05

    dModel = pgcore.RVector(len(model))
    nProcs = self.multiThreadJacobian()

    if sys.platform == 'win32' or sys.platform == 'darwin':
        # strange pickle problem: see  python test_PhysicsManagers.py ves
        from .logger import warn
        warn('Multiprocess Jacobian currently unavailable for Win32 and Mac.')
        nProcs = 1

    if nProcs == 1:
        self.createJacobian(model, resp)
        return

    shm = []

    oldBertThread = self.threadCount()
    self.setThreadCount(1)

    # print("Model/Data/nProcs", nModel, nData, nProcs,
    #       int(ceil(float(nModel)/nProcs)))
    for pCount in range(int(ceil(float(nModel) / nProcs))):
        procs = []
        if self.verbose():
            print("Jacobian MT:(", pCount * nProcs, "--",
                  (pCount + 1) * nProcs, ") /", nModel, '... ')

        for i in range(int(pCount * nProcs), int((pCount + 1) * nProcs)):
            if i < nModel:
                modelChange = pgcore.RVector(model)
                modelChange[i] *= fak
                dModel[i] = modelChange[i] - model[i]

                shm.append(Array('d', len(resp)))
                procs.append(
                    Process(target=__GLOBAL__response_mt_shm_,
                            args=(self, modelChange, shm[i], i)))

        for i, p in enumerate(procs):
            p.start()

        for i, p in enumerate(procs):
            p.join()

        # if self.verbose():
        #     print(dur(), 's')
    self.setThreadCount(oldBertThread)

    for i in range(nModel):
        dData = np.array(shm[i]) - resp
        self._J.setCol(i, dData / dModel[i])


def __ModellingBase__responses_mt__(self, models, respos):

    nModel = len(models)
    nProcs = self.multiThreadJacobian()

    if nProcs == 1:
        for i, m in enumerate(models):
            respos[i] = self.response_mt(m, i)
        return

    from math import ceil
    from multiprocessing import Process, Array
    import numpy as np

    if models.ndim != 2:
        raise BaseException("models need to be a matrix(N, nModel):" +
                            str(models.shape))
    if respos.ndim != 2:
        raise BaseException("respos need to be a matrix(N, nData):" +
                            str(respos.shape))

    nData = len(respos[0])
    shm = []

    oldBertThread = self.threadCount()
    self.setThreadCount(1)

    # print("*"*100)
    # print(nModel, nProcs)
    # print("*"*100)
    for pCount in range(int(ceil(nModel / nProcs))):
        procs = []
        if self.verbose():
            print(pCount * nProcs, "/", nModel)
        for i in range(int(pCount * nProcs), int((pCount + 1) * nProcs)):

            if i < nModel:
                shm.append(Array('d', nData))
                procs.append(
                    Process(target=__GLOBAL__response_mt_shm_,
                            args=(self, models[i], shm[i], i)))

        for i, p in enumerate(procs):
            p.start()

        for i, p in enumerate(procs):
            p.join()

    self.setThreadCount(oldBertThread)

    for i in range(nModel):
        resp = np.array(shm[i])
        respos[i] = resp


class ModellingBaseMT__(pgcore.ModellingBase):

    def __init__(self, mesh=None, dataContainer=None, verbose=False):
        if mesh and dataContainer:
            pgcore.ModellingBase.__init__(
                self, mesh=mesh, dataContainer=dataContainer, verbose=verbose)
        elif isinstance(mesh, pgcore.Mesh):
            pgcore.ModellingBase.__init__(self, mesh=mesh, verbose=verbose)
        elif dataContainer:
            pgcore.ModellingBase.__init__(self, dataContainer=dataContainer,
                                          verbose=verbose)
        else:
            pgcore.ModellingBase.__init__(self, verbose=verbose)

        self._J = pgcore.RMatrix()
        self.setJacobian(self._J)


ModellingBaseMT__.createJacobian_mt = __ModellingBase__createJacobian_mt__
ModellingBaseMT__.responses = __ModellingBase__responses_mt__

ModellingBase = ModellingBaseMT__

###########################
# unsorted stuff
###########################

# DEPRECATED
# pgcore.interpolate = pgcore.interpolate_GILsave__

############################
# some backward compatibility
############################


def __getCoords(coord, dim, ent):
    """Syntactic sugar to find all x-coordinates of a given entity."""

    if isScalar(ent):
        return ent
    if isPos(ent):
        return ent[dim]
    if isinstance(ent, np.ndarray) and (ent.shape[1] >= dim and ent.shape[1] <= 3):
        return ent[:, dim]
    if isinstance(ent, (pgcore.R3Vector)):
        return getattr(pgcore, coord)(ent)
    if isinstance(ent, list) and isinstance(ent[0], pgcore.Pos):
        return getattr(pgcore, coord)(ent)
    if isinstance(ent, list) and isPos(ent[0]):
        return getattr(pgcore, coord)(ent)
    if isinstance(ent, DataContainer):
        return getattr(pgcore, coord)(ent.sensorPositions())
    if isinstance(ent, Mesh):
        return getattr(pgcore, coord)(ent.positions())
    if isinstance(ent, pgcore.stdVectorNodes):
        return np.array([n.pos()[dim] for n in ent])
    if isinstance(ent, pgcore.Node):
        return ent.pos()[dim]
    if isinstance(ent, list) and isinstance(ent[0], pgcore.Node):
        return [n.pos()[dim] for n in ent]

    if hasattr(ent, 'ndim') and ent.ndim == 2:
        if ent.shape[0] > 3 and ent.shape[1] > dim:
            ## for sure (N x [x, y, z])
            return ent[:, dim]
        if ent.shape[1] > 3 and ent.shape[0] > dim:
            ## for sure ([x, y, z] x N)
            return ent[dim]

        if hasattr(ent, 'flags') and ent.flags['F_CONTIGUOUS'] == True or \
            (ent.shape[0] >= dim and ent.shape[0] <= 3) or \
            ent.shape[1] > dim:
            ## (N x [x, y, z]).T || ([x, y, z] x N)
            return ent[dim]

        if ent.shape[0] > dim or \
            (ent.shape[1] >= dim and ent.shape[1] <= 3):
            ## assuming (Nx[x, y, z])
            return ent[:, dim]

    # use logger here
    raise Exception(
        "Don't know how to find the " + coord + "-coordinates of entity:", ent)


def x(instance):
    """Syntactic sugar to find all x-coordinates of a given class instance.

    Convenience function to return all associated x-coordinates
    of a given class instance. If instance is 1D array its returned itself.

    Parameters
    ----------
    instance : DataContainer, Mesh, R3Vector, np.array, list(Pos)
        Return the associated coordinate positions for given class instance.

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> pg.x([[1.0, 1.0, 1.0]])
    1 [1.0]
    >>> pg.x([[0, 0], [1, 0]])
    2 [0.0, 1.0]
    """
    ## test atoms first (perfomance)
    if isPos(instance):
        return instance[0]
    if isArray(instance) and not isPos(instance):
        return instance
    return __getCoords('x', 0, instance)


def y(instance):
    """Syntactic sugar to find all y-coordinates of a given class instance.

    Convenience function to return all associated x-coordinates
    of a given class instance.

    Parameters
    ----------
    instance : DataContainer, Mesh, R3Vector, np.array, list(Pos)
        Return the associated coordinate positions for given class instance.
    """
    ## test atoms first (perfomance)
    if isPos(instance):
        return instance[1]
    if isinstance(instance, (int,float)):
        return 0
    return __getCoords('y', 1, instance)


def z(instance):
    """Syntactic sugar to find all z-coordinates of a given class instance.

    Convenience function to return all associated x-coordinates
    of a given class instance.

    Parameters
    ----------
    instance : DataContainer, Mesh, R3Vector, np.array, list(Pos)
        Return the associated coordinate positions for given class instance.
    """
    ## test atoms first (perfomance)
    if isPos(instance):
        return instance[2]
    if isinstance(instance, (int,float)):
        return 0
    return __getCoords('z', 2, instance)


def search(what):
    """Utility function to search docstrings for string `what`."""
    np.lookfor(what, module="pygimli", import_modules=False)


from .base import (isInt, isScalar, isIterable, isArray, isPos, isR3Array,
                   isVecField, isPosList, isComplex, isMatrix)

# Import from submodules at the end
from .mesh import Mesh, MeshEntity, Node
from .datacontainer import DataContainer, DataContainerERT
from .trans import *  # why do we need that?

# from .matrix import (Cm05Matrix, LMultRMatrix, LRMultRMatrix, MultLeftMatrix,
#                      MultLeftRightMatrix, MultRightMatrix, RMultRMatrix)
from .matrix import (BlockMatrix, SparseMatrix, SparseMapMatrix, IdentityMatrix,
                     Matrix)

from .matrix import (asCSC, asCSR, asCOO, asSparseMapMatrix, asSparseMatrix,
                     asDense, reduceEntries, removeEntries,
                     )
