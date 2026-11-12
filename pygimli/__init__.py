"""pyGIMLi - open-source library for modelling and inversion in geophysics."""

import sys
import locale
from pathlib import Path

from . import core
from .core.decorators import (
    renamed,
    singleton,
    moduleProperty,
    skipOnDefaultTest,
    deprecate,
)

# Import everything that should be accessible through main namespace.
from .core import BVector, CVector, IVector, RVector, Vector

from .core import (
    DataContainer,
    DataContainerERT,
    Line,
    Mesh,
    Plane,
    Pos,
    PosList,
    PosVector,
    Stopwatch,
    Swatches,
    abs,
    cat,
    center,
    exp,
    find,
    log,
    log10,
    logDropTol,
    max,
    mean,
    median,
    min,
    search,
    sort,
    sum,
    trans,
    versionStr,
    x,
    y,
    z,
    zero,
)

from .core import (
    isInt,
    isScalar,
    isIterable,
    isArray,
    isPos,
    isR3Array,
    isPosList,
    isVecField,
    isComplex,
    isMatrix,
)

from .core.base import isSquareMatrix

from .core.threads import (
    setThreadCount,
    threadsInfo,
    getOMPThreadCount,
    getBLASThreadCount,
    setOMPThreadCount,
    setBLASThreadCount,
)

from .core import math  # alias all from .core.math.* to pg.math.*

# from .core import matrix # alias all from .core.matrix.* to pg.matrix.*
from .core.matrix import BlockMatrix, Matrix, SparseMapMatrix, SparseMatrix

from .core.logger import (
    _,
    _d,
    _g,
    _r,
    _y,
    _b,
    critical,
    d,
    deprecated,
    renameKwarg,
    renameArg,
    error,
    info,
    debug,
    setDebug,
    setLogLevel,
    setVerbose,
    v,
    verbose,
    warn,
)

warning = warn  # conveniance

from .core.config import getConfigPath, rc, getCPUCount

from .meshtools import createGrid, interpolate
from .solver import solve
from .utils import boxprint, cache, cut, unique, unit, cmap, randn
from .utils import prettify as pf
from .utils.utils import Report, Table


def ppf(value):
    """Print *value* with the pretty formatter."""
    print(pf(value))


from .viewer import show, wait, noShow, hold

from .frameworks import fit
from .frameworks import Modelling
from .frameworks import Inversion
from .testing import test  # , setTestingMode, testingMode

from .math import matrix  # alias all from .core.matrix.* to pg.matrix.*
from .core.load import (
    load,
    optImport,
    getCachePath,
    getExampleFile,
    getExampleData,
)


def checkAndFixLocaleDecimal_point(verbose=False):  # noqa: F811
    """Switch the numeric locale to a decimal point if required."""
    if locale.localeconv()["decimal_point"] == "," and verbose:
        print(
            "Found locale decimal_point ',' and change it to: decimal point '.'"
        )
    try:
        locale.localeconv()["decimal_point"]
        locale.setlocale(locale.LC_NUMERIC, "C")
    except Exception as e:
        print(e)
        print("cannot set locale to decimal point")

    # LC_CTYPE should be something with UTF-8
    # export LC_CTYPE="de_DE.UTF-8"
    # python -c 'import sys; print(sys.stdout.encoding)'


checkAndFixLocaleDecimal_point(verbose=True)
# print(locale.localeconv()['decimal_point'])
# if locale.localeconv()['decimal_point'] == ',':
#   print("Found locale decimal_point ',' and change it to: decimal point '.'")
# try:
#    locale.localeconv()['decimal_point']
#    locale.setlocale(locale.LC_NUMERIC, 'C')
# except:
#    print('cannot set locale to decimal point')


if "--debug" in sys.argv or "-d" in sys.argv:
    setDebug(True)
else:
    setDebug(False)

if "--verbose" in sys.argv or "-v" in sys.argv:
    setVerbose(True)
else:
    setVerbose(False)

# if '--test' in sys.argv or '-t' in sys.argv:
#     setTestingMode(True)
# else:
#     setTestingMode(False)


###############################################################################
# Please leave this block here until the following issue is fixed:
# https://github.com/ContinuumIO/anaconda-issues/issues/1068
# if "conda" in __path__[0]:
#     try:
#         import PyQt5
#         import matplotlib
#         matplotlib.use("qt5agg", warn=False)
#     except ImportError:
#         pass
###############################################################################


def isNotebook():
    """Determine if run inside Jupyter notebook or Spyder."""
    import sys

    return "ipykernel_launcher.py" in sys.argv[0] or "ipykernel" in sys.modules


def isIPyTerminal():
    """Determine if run inside ipython terminal, e.g., sphinx-gallery."""
    import sys

    return "IPython" in sys.modules


class SWatches:
    """Singleton class to access Stopwatch instances."""

    @staticmethod
    def __getitem__(key: str):
        """Return Stopwatch with specific key name."""
        if "*" in key:
            allKeys = SWatches().keys()
            import fnmatch

            matches = fnmatch.filter(allKeys, key)
            if len(matches) == 0:
                raise KeyError(f"No stopwatch found for: {key}")
            elif len(matches) > 1:
                return [SWatches()[k] for k in matches]
            else:
                key = matches[0]

        # return SWatches()['/'+key]
        return Swatches.instance()[key]

    @staticmethod
    def keys():
        """Return list of all stopwatch keys."""
        return list(Swatches.instance().keys())

    @staticmethod
    def remove(key: str, recursive: bool = False) -> None:
        """Remove and delete the stopwatch with the given key.

        Arguments
        ---------
        key: str
            Root swatch key name to remove.
        recursive: bool [False]
            If True, remove all timings for sub-keys as well.
        """
        return Swatches.instance().remove(key, recursive)


def tic(msg=None, key=""):
    """Start global timer. Print elapsed time with `toc()`.

    You can start multiple stopwatches with optional identifier.

    Parameters
    ----------
    msg : string, optional
        Print message string just before starting the timer.
    key: identifier
        Identifier for your Stopwatch.
    """
    if msg:
        print(msg)

    SWatches()["/" + key].start()


def toc(msg=None, box=False, stop=False, reset=False, key=""):
    """Print elapsed time since global timer was started with `tic()`.

    Arguments
    ---------
    msg: string [None]
        Print message string just after printing the elapsed time. If box is
        True, then embed msg into its own box
    box: bool [False]
        Embed the time in an ascii box
    stop: bool [False]
        Stops the stopwatch.
    reset: bool [False]
        Reset timer to 0.0 but don't stop it. Empties stored values.
    key: identifier
        Identifier for your Stopwatch.
    """
    if msg:
        if box is True:
            boxprint(msg)
        else:
            print(msg, end=" ")

    seconds = dur(key, stop=stop, reset=reset)

    ## refactor with prettyTime
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if h <= 0 and m <= 0:
        time = pf(s)
    elif h <= 0:
        if m == 1.0:
            time = f"{int(m)} minute and {s:.2f}"
        else:
            time = f"{int(m)} minutes and {s:.2f}"
    elif h == 1.0:
        time = f"{int(h)} hour, {int(m)} minutes and {s:.2f}"
    else:
        time = f"{int(h)} hours, {int(m)} minutes and {s:.2f}"
    p = print if not box else boxprint

    if len(SWatches().keys()) and key != "":
        p(f"Elapsed time ({key}) is {time} seconds.")
    else:
        p(f"Elapsed time is {time} seconds.")


def dur(key="", stop=False, reset=False):
    """Return time in seconds since global timer was started with `tic()`.

    Arguments
    ---------
    key: identifier
        Identifier for your Stopwatch.
    stop: bool [False]
        Stops the stopwatch.
    reset: bool [False]
        Reset timer to 0.0 but don't stop it. Empties stored values.
    """
    if isinstance(stop, str):
        key = stop

    if stop is True:
        SWatches()["/" + key].stop()
    return SWatches()["/" + key].duration(restart=reset)


def store(key="", stop=True):
    """Store the current time since the global timer was started.

    Arguments
    ---------
    stop: bool [True]
        Reset the stopwatch.
    key: identifier
        Identifier for your Stopwatch.
    """
    if stop is True:
        SWatches()["/" + key].stop()

    return SWatches()["/" + key].store()


class tictoc:
    """Context manager for timing code blocks."""

    def __init__(self, key: str):
        """Initialize the tictoc context manager."""
        self._tt = core.TicToc(key)

    def __enter__(self):
        """Start the timer when entering the context."""
        return self

    def __exit__(self, type, value, traceback):
        """Stop the timer when exiting the context."""
        ## no explicit stop needed, destructor will handle it
        del self._tt


def swatch(key):
    """Return Stopwatch with specific key name.

    Arguments
    ---------
    key: identifier
        Identifier for your Stopwatch.
        Single wildcards '*' are supported.

    Returns
    -------
    Stopwatch instance or list of Stopwatch instances if multiple matches found.
    """
    if key.startswith("/"):
        key = key[1:]
    return SWatches()["/" + key]


def timings(name="/"):
    """Return table of timings for a given root swatch key name."""
    import numpy as np

    class TTree:
        def __init__(self, parent=None):
            self.parent = parent
            self.name = None
            self.childs = {}
            self.data = None

        def __getitem__(self, k):
            names = k.split("/")

            if self.name is None:
                self.name = names[0]
            else:
                if names[0] != self.name:
                    error(f"Wrong tree({self.name}) for {k}")

            if len(names) > 1:
                if names[1] not in self.childs:
                    self.childs[names[1]] = TTree(parent=self)

                return self.childs[names[1]]["/".join(names[1:])]

            return self

        @property
        def fullname(self):
            ps = "" if self.name is None else self.name

            p = self.parent
            while 1:
                try:
                    ps = p.name + "/" + ps
                    p = p.parent
                except BaseException:
                    break

            return ps + ":" + str(self.data)

        def __str__(self):
            s = self.fullname + "\n"
            for _n, t in self.childs.items():
                s += str(t)

            return s

    tree = TTree()
    header = ["", "single", "count", "sum", "rel.(%)", "uncov.(%)"]
    table = []

    # if len(SWatches().items()) == 0:
    #     pg.error('')
    # _g(SWatches().keys())
    maxTime = 0

    if not name.startswith("/"):
        name = "/" + name

    for k in list(SWatches().keys()):
        s = SWatches()[k]

        if isinstance(k, str) and k.startswith(name):
            ts = s.stored()
            if ts is None:
                print("no swatch for: ", k)
                table.append([k, 0.0, 0, 0.0, "", None])
            elif len(ts) == 0:
                print("no stored times for: ", s.duration(), k)
                table.append([k, s.duration(), 1, s.duration(), "", None])
            else:
                #     sts = sum(ts)
                maxTime = max(float(maxTime), sum(ts))

                perc = 0
                try:
                    perc = int(sum(ts) / maxTime * 100)
                except ZeroDivisionError:
                    pass

                table.append(
                    [
                        k,
                        np.mean(ts),
                        len(ts),
                        sum(ts),
                        str(perc).rjust(3, "-").rjust(3 * (k.count("/")), "-"),
                        None,
                    ]
                )

                tree[k].data = sum(ts)

        # print(f'\t{k}: {pg.pf(sts)}s ({len(ts)}x{pg.pf(np.mean(ts)*1000)}ms)')

    for row in table:
        if len(list(tree[row[0]].childs.keys())) > 0:
            row[-1] = row[-3]
            for n, tc in tree[row[0]].childs.items():
                try:
                    if tc.data is not None:
                        row[-1] -= tc.data
                except BaseException as e:
                    print("error for: ", n, tc.fullname, tc.data)
                    print(e)
                    pass

            row[-1] = f"{pf(row[-1] / maxTime * 100)}"

    if len(table) == 0:
        error(f"No timings for: {name}")
        return

    return Table(table, header, align="lrcrlr", transpose=False)


def removeTimings(name: str, recursive: bool = False) -> None:
    """Remove all timings for a given root swatch key name.

    Arguments
    ---------
    name: str
        Root swatch key name to remove.
    recursive: bool [False]
        If True, remove all timings for subkeys as well.
    """
    if not name.startswith("/"):
        name = "/" + name
    SWatches().remove(name, recursive=recursive)


# special shortcut pg.plt with lazy evaluation
__MPL_PLT__ = None


@moduleProperty
def _plt():
    # import time
    # t0 = time.time()
    global __MPL_PLT__, rc

    if __MPL_PLT__ is None:
        if rc["matplotlib"] is not None:
            try:
                get_ipython().run_line_magic("matplotlib", rc["matplotlib"])
                debug("matplotlib notebook backend set to: ", rc["matplotlib"])
            except NameError:
                pass
            except BaseException as e:
                info(
                    "matplotlib notebook backend set to "
                    f"{rc['matplotlib']} failed: ",
                    e,
                )
        # tic()
        import matplotlib.pyplot as plt

        # if isNotebook():
        #     pass
        # else:
        #     import matplotlib
        #     matplotlib.use('qtagg')
        #     print('############### importing plt took ', dur())

        #     print('############### backend:', plt.get_backend())

        from .viewer.mpl import registerShowPendingFigsAtExit, hold

        registerShowPendingFigsAtExit()

    return plt


def findVersion():
    """Find current version generated by versioneer."""
    # setDebug(False)
    root = Path(__file__).resolve().parent
    gitPath = root / ".git"

    debug("Fetching version info.")
    from ._version import get_versions

    _versions = get_versions()

    if "error" in _versions and _versions["error"] is not None:
        from importlib.metadata import version

        return version("pygimli")

    _version = _versions["version"]

    def _get_branch():
        """Get current git branch."""
        from shutil import which
        from subprocess import check_output

        if not gitPath.exists():
            return None
        git = which("git")
        if git is None:
            return None
        out = check_output(  # noqa: S603 - path resolved via shutil.which
            [
                git,
                "--git-dir",
                str(gitPath),
                "rev-parse",
                "--abbrev-ref",
                "HEAD",
            ]
        ).decode("utf8")

        branch = out.split("\n")[0]
        if "HEAD" not in branch:
            return branch
        return None

    # def _get_latest_tag():
    #     from os.path import exists

    #     if exists(gitPath):
    #         from subprocess import check_output
    #         out = check_output(["git", "--git-dir", gitPath,
    #             "describe", "--tag"]).decode("utf8")

    #         tag = out.split("\n")[0].split('-')[0]
    #         return tag
    #     return None

    _branch = _get_branch()

    if _versions["dirty"]:
        _version = _version.replace(".dirty", " (with local changes")

        if _branch:
            _version += f" on {_branch} branch)"
        else:
            _version += ")"
    elif _branch and "+" in _version:
        _version += f" ({_branch})"

    return _version


def version():
    """Shortcut to show and return current version."""
    pg = sys.modules[__name__]
    vers = pg.__version__  # triggers lazy evaluation
    info("Version: " + vers + " core:" + versionStr())
    return vers


# call once to get version from cache, setup or _version.py
# patch __version__ into the main module class for lazy evaluation
class pygimli(sys.modules[__name__].__class__):
    """Module class with a lazily evaluated ``__version__``."""

    __version__ = property(lambda self: self.findVersion())


sys.modules[__name__].__class__ = pygimli
