#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Collection of several utility functions."""

from __future__ import print_function
import sys
from math import floor, sqrt

import numpy as np
import pygimli as pg

# scooby is a soft dependency.
try:
    from scooby import Report as ScoobyReport
except ImportError:
    class ScoobyReport:
        """Local scooby reporting class."""

        def __init__(self, *args, **kwargs):
            """Do nothing."""
            pass

        def __repr__(self):
            """Representation."""
            message = (
                "`Report` requires `scooby`. Install via `pip install scooby` "
                "or `conda install -c conda-forge scooby`."
            )
            return message

        def to_dict(self):
            """Dictionary representation (empty for now)."""
            return {}


class ProgressBar(object):
    """Animated text-based progress bar.

    Animated text-based progressbar for intensive loops. Should work in the
    console. In IPython Notebooks a 'tqdm' progressbar instance is created and
    can be configured with appropriate keyword arguments.
    """

    def __init__(self, its, width=80, sign=":", **kwargs):
        """Create animated text-based progress bar.

        Todo
        ----
        * optional: 'estimated time' instead of 'x of y complete'

        Parameters
        ----------
        its : int
            Number of iterations of the process.
        width : int
            Width of the ProgressBar, default is 80.
        sign : str
            Sign used to fill the bar.

        Additional Args
        ---------------
        Forwarded to create the tqdm progressbar instance. See
        https://tqdm.github.io/docs/tqdm/

        Examples
        --------
        >>> from pygimli.utils import ProgressBar
        >>> pBar = ProgressBar(its=20, width=40, sign='+')
        >>> pBar.update(5)
        \r[+++++++++++       30%                 ] 6 of 20 complete
        """
        self.its = int(its)
        self.width = width
        self.sign = sign[0]  # take first character only if sign is longer
        self.pBar = "[]"
        self._amount(0)
        self.nbProgress = None
        self._iter = -1

        if pg.isNotebook():
            tqdm = pg.optImport(
                'tqdm', requiredFor="use nice progressbar in jupyter notebook")
            if tqdm is not None:
                from tqdm.notebook import tqdm
                fmt = kwargs.pop(
                    'bar_format',
                    '{desc} {percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt}' +
                    ' [{elapsed} < {remaining}]')
                self.nbProgress = tqdm(total=its, bar_format=fmt, **kwargs)

    def __call__(self, it, msg=""):
        """Update progress."""
        self.update(it, msg)

    def update(self, iteration, msg=""):
        """Update by iteration number starting at 0 with optional message."""
        if self.nbProgress is not None:
            # TODO maybe catch if someone doesn't call with iteration steps=1
            self.nbProgress.update(n=iteration-self._iter)
        else:
            self._setbar(iteration + 1)
            if len(msg) >= 1:
                self.pBar += " (" + msg + ")"
            print("\r" + self.pBar, end="")
            sys.stdout.flush()

        # last iteration here
        if iteration == self.its-1:
            if self.nbProgress is not None:
                self.nbProgress.close()
            else:
                print()

        self._iter = iteration

    def _setbar(self, elapsed_it):
        """Reset pBar based on current iteration number."""
        self._amount((elapsed_it / float(self.its)) * 100.0)
        self.pBar += " %d of %s complete" % (elapsed_it, self.its)

    def _amount(self, new_amount):
        """Calculate amount by which to update the pBar."""
        pct_done = int(round((new_amount / 100.0) * 100.0))
        full_width = self.width - 2
        num_signs = int(round((pct_done / 100.0) * full_width))
        self.pBar = "[" + self.sign * num_signs + \
            " " * (full_width - num_signs) + "]"
        pct_place = (len(self.pBar) // 2) - len(str(pct_done))
        pct_string = " %d%% " % pct_done
        self.pBar = self.pBar[0:pct_place] + \
            (pct_string + self.pBar[pct_place + len(pct_string):])


def boxprint(s, width=80, sym="#"):
    """Print string centered in a box.

    Examples
    --------
    >>> from pygimli.utils import boxprint
    >>> boxprint("This is centered in a box.", width=40, sym='+')
    ++++++++++++++++++++++++++++++++++++++++
    +      This is centered in a box.      +
    ++++++++++++++++++++++++++++++++++++++++
    """
    row = sym * width
    centered = str(s).center(width - 2)
    print("\n".join((row, centered.join((sym, sym)), row)))


def trimDocString(docstring):
    """Return properly formatted docstring.

    From: https://www.python.org/dev/peps/pep-0257/

    Examples
    --------
    >>> from pygimli.utils import trimDocString
    >>> docstring = '    This is a string with indention and whitespace.   '
    >>> trimDocString(docstring).replace('with', 'without')
    'This is a string without indention and whitespace.'
    """
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = 2**16 - 1
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < 2**16 - 1:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a single string:
    return '\n'.join(trimmed)


def unicodeToAscii(text):
    """TODO DOCUMENTME."""
    if isinstance(text, str):
        return text.encode("iso-8859-1", "ignore")
    else:
        return text


def logDropTol(p, dropTol=1e-3):
    """Create logarithmic scaled copy of p.

    Examples
    --------
    >>> from pygimli.utils import logDropTol
    >>> x = logDropTol((-10, -1, 0, 1, 100))
    >>> print(x.array())
    [-4. -3.  0.  3.  5.]
    """
    tmp = pg.Vector(p)

    tmp = pg.abs(tmp / dropTol)
    tmp.setVal(1.0, pg.find(tmp < 1.0))

    tmp = pg.log10(tmp)
    tmp *= pg.math.sign(p)
    return tmp


def prettify(value, roundValue=False):
    """Return prettified string for value .. if possible."""
    if isinstance(value, dict):
        import json
        # class CustomEncoder(json.JSONEncoder):
        #     def __init__(self, *args, **kwargs):
        #         super().__init__(*args, **kwargs)

        #     def _iterencode(self, o):
        #         try:
        #             return super()._iterencode(o)
        #         except:
        #             return "{0} is not JSON serializable".format(type(o))

        try:
            return json.dumps(value, indent=4)
        except Exception as e:
            pg.warning('prettify fails:', e)
            return str(value)
    elif pg.isScalar(value):
        return prettyFloat(value, roundValue)
    pg.warn("Don't know how to prettify the string representation for: ",
            type(value), value)
    return value


def prettyFloat(value, roundValue=None):
    """Return prettified string for a float value.

    Todo
    ----
        add number for round to
        add test
    """
    # test-cases:
    # if change things her, look that they are still good (mod-dc-2d)

    if (isinstance(roundValue, int) and abs(round(value)-value) < 1e-4 and
            abs(value) < 1e3 and 0):
        string = str(int(round(value, roundValue)))
    elif abs(value) < 1e-14:
        string = "0"
    elif abs(value) > 1e4 or abs(value) <= 1e-3:
        string = str("%.1e" % value)
    elif abs(value) < 1e-2:
        string = str("%.4f" % round(value, 4))
    # max two values after comma
    # elif abs(value) < 1e-1:
    #     string = str("%.3f" % round(value, 3))
    elif abs(value) < 1e0:
        string = str("%.2f" % round(value, 2))
    elif abs(value) < 1e1:
        string = str("%.2f" % round(value, 2))
    elif abs(value) < 1e2:
        string = str("%.2f" % round(value, 2))
    else:
        string = str("%.0f" % round(value, 2))

    # pg._y(string)
    # print(string.endswith("0") and string[-2] == '.')
    if string.endswith(".0"):
        # pg._r(string.replace(".0", ""))
        return string.replace(".0", "")
    elif string.endswith(".00"):
        return string.replace(".00", "")
    elif '.' in string and 'e' not in string and string.endswith("00"):
        return string[0:len(string)-2]
    elif '.' in string and 'e' not in string and string.endswith("0"):
        # pg._r(string[0:len(string)-1])
        return string[0:len(string)-1]
    else:
        return string


def prettyTime(t):
    """Return prettified time in seconds as string. No months, no leap year.

    TODO
    ----
        * weeks (needed)
        * > 1000 years

    Args
    ----
    t: float
        Time in seconds, should be > 0

    Examples
    --------
    >>> from pygimli.utils import prettyTime
    >>> print(prettyTime(1))
    1 s
    >>> print(prettyTime(3600*24))
    1 day
    >>> print(prettyTime(2*3600*24))
    2 days
    >>> print(prettyTime(365*3600*24))
    1 year
    >>> print(prettyTime(3600))
    1 hour
    >>> print(prettyTime(2*3600))
    2 hours
    >>> print(prettyTime(3660))
    1h1m
    >>> print(prettyTime(1e-3))
    1 ms
    >>> print(prettyTime(1e-6))
    1 µs
    >>> print(prettyTime(1e-9))
    1 ns
    """
    if abs(t) > 1:
        seconds = int(t)
        years, seconds = divmod(seconds, 365*86400)
        days, seconds = divmod(seconds, 86400)
        hours, seconds = divmod(seconds, 3600)
        minutes, seconds = divmod(seconds, 60)
        if years > 0:
            if days >= 1:
                return '%dy%dd' % (years, days)
            else:
                if years > 1:
                    return '%d years' % (years,)
                else:
                    return '%d year' % (years,)
        elif days > 0:
            if hours >= 1:
                return '%dd%dh' % (days, hours)
            else:
                if days > 1:
                    return '%d days' % (days,)
                else:
                    return '%d day' % (days,)
        elif hours > 0:
            if minutes >= 1:
                return '%dh%dm' % (hours, minutes)
            else:
                if hours > 1:
                    return '%d hours' % (hours)
                else:
                    return '%d hour' % (hours)
        elif minutes > 0:
            if seconds >= 1:
                return '%dm%ds' % (minutes, seconds)
            else:
                if minutes > 1:
                    return '%d minutes' % (minutes)
                else:
                    return '%d minute' % (minutes)
        else:
            return '%d s' % (seconds,)
    else:
        if abs(t) >= 1e-3 and abs(t) <= 0.1:
            return prettyFloat(t*1e3) + " ms"
        elif abs(t) >= 1e-6 and abs(t) <= 1e-3:
            return prettyFloat(t*1e6) + " µs"
        elif abs(t) >= 1e-9 and abs(t) <= 1e-6:
            return prettyFloat(t*1e9) + " ns"
        return prettyFloat(t) + " s"


def niceLogspace(vMin, vMax, nDec=10):
    """Nice logarithmic space from decade < vMin to decade > vMax.

    Parameters
    ----------
    vMin : float
        lower limit need to be > 0
    vMax : float
        upper limit need to be >= vMin
    nDec : int
        Amount of logarithmic equidistant steps for one decade

    Examples
    --------
    >>> from pygimli.utils import niceLogspace
    >>> v1 = niceLogspace(vMin=0.1, vMax=0.1, nDec=1)
    >>> print(v1)
    [0.1 1. ]
    >>> v1 = niceLogspace(vMin=0.09, vMax=0.11, nDec=1)
    >>> print(v1)
    [0.01 0.1  1.  ]
    >>> v1 = niceLogspace(vMin=0.09, vMax=0.11, nDec=10)
    >>> print(len(v1))
    21
    >>> print(v1)
    [0.01       0.01258925 0.01584893 0.01995262 0.02511886 0.03162278
     0.03981072 0.05011872 0.06309573 0.07943282 0.1        0.12589254
     0.15848932 0.19952623 0.25118864 0.31622777 0.39810717 0.50118723
     0.63095734 0.79432823 1.        ]
    """
    if vMin > vMax or vMin < 1e-12:
        print("vMin:", vMin, "vMax", vMax)
        raise Exception('vMin > vMax or vMin <= 0.')

    vMin = 10**np.floor(np.log10(vMin))
    vMax = 10**np.ceil(np.log10(vMax))

    if vMax == vMin:
        vMax *= 10

    n = np.log10(vMax / vMin) * nDec + 1

    q = 10.**(1. / nDec)

    return vMin * q**np.arange(n)


def grange(start, end, dx=0, n=0, log=False):
    """Create array with possible increasing spacing.

    Create either array from start step-wise filled with dx until end reached
    [start, end] (like np.array with defined end).
    Fill the array from start to end with n steps.
    [start, end] (like np.linespace)
    Fill the array from start to end with n steps but logarithmic increasing,
    dx will be ignored.

    Parameters
    ----------
    start: float
        First value of the resulting array
    end: float
        Last value of the resulting array
    dx: float
        Linear step length, n will be ignored
    n: int
        Amount of steps
    log: bool
        Logarithmic increasing range of length = n from start to end.
        dx will be ignored.

    Examples
    --------
    >>> from pygimli.utils import grange
    >>> v1 = grange(start=0, end=10, dx=3)
    >>> v2 = grange(start=0, end=10, n=3)
    >>> print(v1)
    4 [0.0, 3.0, 6.0, 9.0]
    >>> print(v2)
    3 [0.0, 5.0, 10.0]

    Returns
    -------
    ret: :gimliapi:`GIMLI::RVector`
        Return resulting array
    """
    s = float(start)
    e = float(end)
    d = float(dx)

    if dx != 0 and not log:
        if end < start and dx > 0:
            # print("grange: decreasing range but increasing dx, swap dx sign")
            d = -d
        if end > start and dx < 0:
            # print("grange: increasing range but decreasing dx, swap dx sign")
            d = -d
        ret = pg.Vector(range(int(floor(abs((e - s) / d)) + 1)))
        ret *= d
        ret += s
        return ret

    elif n > 0:
        if not log:
            return grange(start, end, dx=(e - s) / (n - 1))
        else:
            return pg.core.increasingRange(start, end, n)[1:]
    else:
        raise Exception('Either dx or n have to be given.')


def diff(v):
    """Calculate approximate derivative.

    Calculate approximate derivative from v as d = [v_1-v_0, v2-v_1, ...]

    Parameters
    ----------
    v : array(N) | pg.PosVector(N)
        Array of double values or positions

    Returns
    -------
    d : [type(v)](N-1) |
        derivative array

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.utils import diff
    >>> p = pg.PosVector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> print(diff(p)[0])
    RVector3: (0.0, 1.0, 0.0)
    >>> print(diff(p)[1])
    RVector3: (0.0, -1.0, 0.0)
    >>> print(diff(p)[2])
    RVector3: (0.0, 0.0, 0.0)
    >>> p = pg.Vector(3)
    >>> p[0] = 0.0
    >>> p[1] = 1.0
    >>> p[2] = 2.0
    >>> print(diff(p))
    2 [1.0, 1.0]
    """
    d = None

    if isinstance(v, np.ndarray):
        if v.ndim == 2:
            if v.shape[1] < 4:
                # v = pg.PosVector(v.T)
                vt = v.copy()
                v = pg.PosVector(len(vt))
                for i, vi in enumerate(vt):
                    v.setVal(pg.RVector3(vi), i)
            else:
                v = pg.PosVector(v)
        else:
            v = pg.Vector(v)
    elif isinstance(v, list):
        v = pg.PosVector(v)

    if isinstance(v, pg.PosVector) or isinstance(v, pg.core.stdVectorRVector3):
        d = pg.PosVector(len(v) - 1)
    else:
        d = pg.Vector(len(v) - 1)

    for i, _ in enumerate(d):
        d[i] = v[i + 1] - v[i]
    return d


def dist(p, c=None):
    """Calculate the distance for each position in p relative to pos c(x,y,z).

    Parameters
    ----------
    p : ndarray(N,2) | ndarray(N,3) | pg.PosVector

        Position array
    c : [x,y,z] [None]
        relative origin. default = [0, 0, 0]

    Returns
    -------
    d : ndarray(N)
        Distance array

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.utils import dist
    >>> import numpy as np
    >>> p = pg.PosVector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> print(dist(p))
    [0. 1. 0. 0.]
    >>> x = pg.Vector(4, 0)
    >>> y = pg.Vector(4, 1)
    >>> print(dist(np.array([x, y]).T))
    [1. 1. 1. 1.]
    """
    if c is None:
        c = pg.RVector3(0.0, 0.0, 0.0)
    d = np.zeros(len(p))
    pI = None
    for i, _ in enumerate(p):
        if isinstance(p[i], pg.RVector3):
            pI = p[i]
        elif isinstance(p[i], float):
            pI = pg.RVector3(p[i], 0)
        else:
            pI = pg.RVector3(p[i])
        d[i] = (pI - c).abs()

    return d


def cumDist(p):
    """The progressive, i.e., cumulative length for a path p.

    d = [0.0, d[0]+ | p[1]-p[0] |, d[1] + | p[2]-p[1] | + ...]

    Parameters
    ----------
    p : ndarray(N,2) | ndarray(N,3) | pg.PosVector
        Position array

    Returns
    -------
    d : ndarray(N)
        Distance array

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.utils import cumDist
    >>> import numpy as np
    >>> p = pg.PosVector(4)
    >>> p[0] = [0.0, 0.0]
    >>> p[1] = [0.0, 1.0]
    >>> p[2] = [0.0, 1.0]
    >>> p[3] = [0.0, 0.0]
    >>> print(cumDist(p))
    [0. 1. 1. 2.]
    """
    d = np.zeros(len(p))
    d[1:] = np.cumsum(dist(diff(p)))
    return d


def cut(v, n=2):
    """Cut array v into n parts."""
    N = len(v)
    Nc = N//n
    cv = [v[i*Nc:(i+1)*Nc] for i in range(n)]
    return cv


def randn(n, seed=None):
    """Create n normally distributed random numbers with optional seed.

    Parameters
    ----------
    n: long
        length of random numbers array.
    seed: int[None]
        Optional seed for random number generator

    Returns
    -------
    r: np.array
        Random numbers.

    Examples
    --------
    >>> import numpy as np
    >>> from pygimli.utils import randn
    >>> a = randn(5, seed=1337)
    >>> b = randn(5)
    >>> c = randn(5, seed=1337)
    >>> print(np.array_equal(a, b))
    False
    >>> print(np.array_equal(a, c))
    True
    """
    if seed is not None:
        np.random.seed(seed)

    if isinstance(n, tuple):
        return np.random.randn(n[0], n[1])
    return np.random.randn(n)


def rand(n, minVal=0.0, maxVal=1.0, seed=None):
    """Create RVector of length n with normally distributed random numbers."""
    if seed is not None:
        np.random.seed(seed)
    return np.random.rand(n) * (maxVal - minVal) + minVal


def getIndex(seq, f):
    """TODO DOCUMENTME."""
    pg.error('getIndex in use?')
    # DEPRECATED_SLOW
    idx = []
    if isinstance(seq, pg.Vector):
        for i, _ in enumerate(seq):
            v = seq[i]
            if f(v):
                idx.append(i)
    else:
        for i, d in enumerate(seq):
            if f(d):
                idx.append(i)
    return idx


def filterIndex(seq, idx):
    """TODO DOCUMENTME."""
    pg.error('filterIndex in use?')
    if isinstance(seq, pg.Vector):
        # return seq(idx)
        ret = pg.Vector(len(idx))
    else:
        ret = list(range(len(idx)))

    for i, ix in enumerate(idx):
        ret[i] = seq[ix]

    return ret


def findNearest(x, y, xp, yp, radius=-1):
    """TODO DOCUMENTME."""
    idx = 0
    minDist = 1e9
    startPointDist = pg.Vector(len(x))
    for i, _ in enumerate(x):
        startPointDist[i] = sqrt((x[i] - xp) * (x[i] - xp) + (y[i] - yp) * (y[
            i] - yp))

        if startPointDist[i] < minDist and startPointDist[i] > radius:
            minDist = startPointDist[i]
            idx = i
    return idx, startPointDist[idx]


def unique_everseen(iterable, key=None):
    """Return iterator of unique elements ever seen with preserving order.

    Return iterator of unique elements ever seen with preserving order.

    From: https://docs.python.org/3/library/itertools.html#itertools-recipes

    Examples
    --------
    >>> from pygimli.utils import unique_everseen
    >>> s1 = 'AAAABBBCCDAABBB'
    >>> s2 = 'ABBCcAD'
    >>> list(unique_everseen(s1))
    ['A', 'B', 'C', 'D']
    >>> list(unique_everseen(s2, key=str.lower))
    ['A', 'B', 'C', 'D']

    See Also
    --------
    unique, unique_rows
    """
    try:
        from itertools import ifilterfalse
    except BaseException:
        from itertools import filterfalse

    seen = set()
    seen_add = seen.add
    if key is None:
        try:
            for element in ifilterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield element
        except BaseException:
            for element in filterfalse(seen.__contains__, iterable):
                seen_add(element)
                yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element


def unique(a):
    """Return list of unique elements ever seen with preserving order.

    Examples
    --------
    >>> from pygimli.utils import unique
    >>> unique((1,1,2,2,3,1))
    [1, 2, 3]

    See Also
    --------
    unique_everseen, unique_rows
    """
    return list(unique_everseen(a))


def uniqueRows(data, precision=2):
    """Equivalent of Matlabs unique(data, 'rows') with tolerance check.

    Additionally returns forward and reverse indices

    Examples
    --------
    >>> from pygimli.utils.utils import uniqueRows
    >>> import numpy as np
    >>> A = np.array(([1,2,3],[3,2,1],[1,2,3]))
    >>> unA, ia, ib = uniqueRows(A)
    >>> bool(np.all(A[ia] == unA))
    True
    """
    pg.deprecated(hint="Please use np.unique(axis=0) instead.")
    fak = 100**precision
    dFix = np.fix(data * fak) / fak + 0.0
    dtype = np.dtype((np.void, dFix.dtype.itemsize * dFix.shape[1]))
    b = np.ascontiguousarray(dFix).view(dtype)
    _, ia, ib = np.unique(b, return_index=True, return_inverse=True)
    return np.unique(b).view(dFix.dtype).reshape(-1, dFix.shape[1]), ia, ib


def uniqueAndSum(indices, to_sum, return_index=False, verbose=False):
    """Sum double values found by indices in a various number of arrays.

    Returns the sorted unique elements of a column_stacked array of indices.
    Another column_stacked array is returned with values at the unique
    indices, while values at double indices are properly summed.

    Parameters
    ----------
    ar : array_like
        Input array. This will be flattened if it is not already 1-D.
    to_sum : array_like
        Input array to be summed over axis 0. Other existing axes will be
        broadcasted remain untouched.
    return_index : bool, optional
        If True, also return the indices of `ar` (along the specified axis,
        if provided, or in the flattened array) that result in the unique
        array.

    Returns
    -------
    unique : ndarray
        The sorted unique values.
    summed_array : ndarray
        The summed array, whereas all values for a specific index is the sum
        over all corresponding nonunique values.
    unique_indices : ndarray, optional
        The indices of the first occurrences of the unique values in the
        original array. Only provided if `return_index` is True.

    Examples
    --------
    >>> import numpy as np
    >>> from pygimli.utils import uniqueAndSum
    >>> idx1 = np.array([0, 0, 1, 1, 2, 2])
    >>> idx2 = np.array([0, 0, 1, 2, 3, 3])
    >>> # indices at positions 0 and 1 and at positions 5 and 6 are not unique
    >>> to_sort = np.column_stack((idx1, idx2))
    >>> # its possible to stack more than two array
    >>> # you need for example 3 array to find unique node positions in a mesh
    >>> values = np.arange(0.1, 0.7, 0.1)
    >>> print(values)
    [0.1 0.2 0.3 0.4 0.5 0.6]
    >>> # some values to be summed together (for example attributes of nodes)
    >>> unique_idx, summed_vals = uniqueAndSum(to_sort, values)
    >>> print(unique_idx)
    [[0 0]
     [1 1]
     [1 2]
     [2 3]]
    >>> print(summed_vals)
    [0.3 0.3 0.4 1.1]
    """
    flag_mult = len(indices) != indices.size
    if verbose:
        print('Get {} indices for sorting'.format(np.shape(indices)))
    if flag_mult:
        ar = indices.ravel().view(
            np.dtype((np.void,
                      indices.dtype.itemsize * indices.shape[1]))).flatten()
    else:
        ar = np.asanyarray(indices).flatten()

    to_sum = np.asanyarray(to_sum)

    if ar.size == 0:
        ret = (ar, )
        ret += (to_sum)
        if return_index:
            ret += (np.empty(0, np.bool), )
        return ret
    if verbose:
        print('Performing argsort...')
    perm = ar.argsort(kind='mergesort')
    aux = ar[perm]
    flag = np.concatenate(([True], aux[1:] != aux[:-1]))
    if flag_mult:
        ret = (indices[perm[flag]], )

    else:
        ret = (aux[flag], )  # unique indices
    if verbose:
        print('Identified {} unique indices'.format(np.shape(ret)))
    if verbose:
        print('Performing reduceat...')
    summed = np.add.reduceat(to_sum[perm], np.nonzero(flag)[0])

    ret += (summed, )  # summed values

    if return_index:
        ret += (perm[flag], )  # optional: indices

    return ret


def filterLinesByCommentStr(lines, comment_str='#'):
    """Filter lines from file.readlines() beginning with symbols in comment."""
    comment_line_idx = []
    for i, line in enumerate(lines):
        if line[0] in comment_str:
            comment_line_idx.append(i)
    for j in comment_line_idx[::-1]:
        del lines[j]
    return lines


class Report(ScoobyReport):
    r"""Report date, time, system, and package version information.

    Use ``scooby`` to report date, time, system, and package version
    information in any environment, either as html-table or as plain text.

    Parameters
    ----------
    additional : {package, str}, default: None
        Package or list of packages to add to output information (must be
        imported beforehand or provided as string).

    """

    def __init__(self, additional=None, **kwargs):
        """Initialize a scooby.Report instance."""
        # Mandatory packages.
        core = ['pygimli', 'pgcore', 'numpy', 'matplotlib']
        # Optional packages.
        optional = ['scipy', 'tqdm', 'IPython', 'meshio', 'tetgen', 'pyvista']
        inp = {
            'additional': additional,
            'core': core,
            'optional': optional,
            **kwargs  # User input overwrites defaults.
        }

        super().__init__(**inp)
