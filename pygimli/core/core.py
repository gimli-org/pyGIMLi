"""Core module of pyGIMLi.

Import the compiled extension module (in-place build or installed package) and
provide small helpers on top of it.
"""

import sys
import traceback

# if sys.platform == 'win32':
#    os.environ['PATH'] = os.path.abspath(__file__) + ';' + os.environ['PATH']

pgcore = None

__PY_MAX__ = max

try:
    # from . import _pygimli_  # if it works: as pgcore, replace all _pygimli_
    # from ._pygimli_ import *  # check if . can be omitted
    # pgcore = _pygimli_

    # print("importing inplace core")
    from .libs import _pygimli_ as pgcore  # noqa: F401  # re-exported
    from .libs._pygimli_ import *  # noqa: F403  # C++ bindings star import
    # print("inplace core")

except ImportError:
    # print("did not find in-place pg core, try import pgcore")
    # import pgcore as _pygimli_
    # from pgcore import _pygimli_  # check version compatibility
    try:
        from pgcore import *  # noqa: F403  # C++ bindings star import
        # pass
    except ImportError:
        traceback.print_exc(file=sys.stdout)
        sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")
        raise


# keep original python max
def max(a, b):
    """Return the maximum of a and b (C++ core version if available)."""
    try:
        return pg.core.max(a, b)  # noqa: F405  # from C++ star import
    except Exception:
        return __PY_MAX__(a, b)
