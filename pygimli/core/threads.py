#!/usr/bin/env python3
"""Process wide control of the BLAS and OpenMP thread pools.

The C++ core can only address the OpenBLAS instance that ``libgimli`` itself
is linked against. numpy and scipy usually bring their own private copies
(loaded with ``RTLD_LOCAL``), so ``openblas_set_num_threads()`` from the core
never reaches them. A program that imports numpy *before* pygimli therefore
keeps using all cores although ``pg.setThreadCount(1)`` was called.

This module closes that gap from Python: it enumerates every shared object
that is currently loaded into the process, probes the known BLAS/OpenMP
thread entry points via :mod:`ctypes` and calls them on each instance that
is found. No third party dependency and no rebuild of the C++ core needed.

Notes
-----
* ``OPENBLAS_NUM_THREADS``, ``MKL_NUM_THREADS``, ``BLIS_NUM_THREADS`` and
  ``OMP_NUM_THREADS`` stay authoritative: :func:`applyEnvironment` never
  overrides them, ``GIMLI_NUM_THREADS`` is only used as a fallback and the
  default determined by ``__setTC__()`` as the last resort.
* That default (``min(8, ncpu - 2)``) is applied to every instance and
  exported to the environment, for unset variables only. numpy and scipy
  read their variables while their private copy is dlopen'ed on first use,
  so they start with the same value. More threads have to be requested
  explicitly.
* Libraries loaded *after* a call are picked up by the next call and get the
  last explicitly requested value applied automatically.
* Platforms without a way to enumerate loaded objects (e.g. Windows) fall
  back to the C++ implementation of the core plus a one time warning.

"""

import ctypes
import os
import re
import sys
import threading
import warnings
from pathlib import Path

__all__ = [
    "blasRuntimes",
    "ompRuntimes",
    "setBLASThreadCount",
    "getBLASThreadCount",
    "setOMPThreadCount",
    "getOMPThreadCount",
    "setThreadCount",
    "threadCount",
    "threadsInfo",
    "applyEnvironment",
]

#: Objects whose basename matches are probed first (informational only).
BLAS_NAME_PATTERN = re.compile(r"openblas|blis|mkl|blas|lapack", re.I)
OMP_NAME_PATTERN = re.compile(r"lib(?:g?omp|iomp5?)[.\-]", re.I)

BLAS_OBJECT_PATTERN = re.compile(r"(openblas|flexiblas|libblis|mkl_)", re.I)
_SHARED_OBJECT_PATTERN = re.compile(r"\.(?:so|dylib)(?:\.[\w.]+)?$")

#: (kind, setter names, getter names) -- order defines the precedence.
BLAS_FAMILIES = (
    (
        "scipy_openblas",
        ("scipy_openblas_set_num_threads64_", "scipy_openblas_set_num_threads"),
        ("scipy_openblas_get_num_threads64_", "scipy_openblas_get_num_threads"),
    ),
    (
        "openblas",
        (
            "openblas_set_num_threads64_",
            "openblas_set_num_threads",
            "openblas_set_num_threads_",
        ),
        (
            "openblas_get_num_threads64_",
            "openblas_get_num_threads",
            "openblas_get_num_threads_",
        ),
    ),
    ("blis", ("bli_thread_set_num_threads",), ("bli_thread_get_num_threads",)),
    (
        "mkl",
        ("mkl_set_num_threads", "mkl_set_num_threads_local"),
        ("mkl_get_max_threads",),
    ),
)

#: Setters that take the thread count *by value* although their name looks
#: like a Fortran alias (trailing ``64_``). Belongs to the scipy-openblas
#: wheel, which describes its ABI as ``scipy_openblas_set_num_threads64_``
#: taking an ``int64`` (measured: a pointer is accepted but ignored).
VALUE_SETTERS = frozenset(
    ("scipy_openblas_set_num_threads64_", "scipy_openblas_set_num_threads")
)

OMP_SYMBOLS = (
    "omp_set_num_threads",
    "omp_get_max_threads",
    "omp_get_thread_limit",
    "omp_get_num_procs",
    "omp_in_parallel",
    "omp_get_num_threads",
)

_lock = threading.RLock()
_handles = {}
_blasRuntimes = []
_ompRuntimes = []
_probedBlas = set()
_probedOmp = set()
_lastBlasTarget = None
_lastOmpTarget = None
_threadCount = None
_warnedBlas = False
_warnedOmp = False
_warnedPlatform = False


#####################
# optional core access
#####################
def _core():
    """Return the C++ extension module or None if not available."""
    try:
        from .core import pgcore
    except Exception:  # pragma: no cover - only without C++ core
        return None
    return pgcore


def _coreCall(name, *args):
    """Call ``pgcore.<name>(*args)``; return None if that is not possible."""
    core = _core()
    if core is None:
        return None
    func = getattr(core, name, None)
    if func is None:
        return None
    try:
        return func(*args)
    except Exception:  # pragma: no cover - defensive
        return None


def _logger():
    """Return the pygimli logger module or None when running standalone."""
    try:
        from . import logger
    except Exception:
        return None
    return logger if hasattr(logger, "info") else None


def _warn(message):
    logger = _logger()
    if logger is not None:
        try:
            logger.warn(message)
        except Exception:  # pragma: no cover - defensive
            pass
        else:
            return
    warnings.warn(message, RuntimeWarning, stacklevel=3)


def _logInfo(message):
    logger = _logger()
    if logger is not None:
        try:
            logger.info(message)
        except Exception:  # pragma: no cover - defensive
            pass
        else:
            return
    print(message)


###########################
# enumeration of shared objects
###########################
def _pathsLinux():
    paths = []
    seen = set()
    try:
        with Path("/proc/self/maps").open(errors="replace") as fp:
            for line in fp:
                fields = line.split()
                if len(fields) < 6:
                    continue
                path = fields[-1]
                if not path.startswith("/") or path in seen:
                    continue
                if path.endswith(" (deleted)"):
                    continue
                if not _SHARED_OBJECT_PATTERN.search(path):
                    continue
                seen.add(path)
                paths.append(path)
    except OSError:  # pragma: no cover - no procfs
        return []
    return paths


def _pathsDarwin():
    try:
        libc = ctypes.CDLL(None)
        count = libc._dyld_image_count
        count.restype = ctypes.c_uint32
        name = libc._dyld_get_image_name
        name.restype = ctypes.c_char_p
        name.argtypes = [ctypes.c_uint32]
    except (AttributeError, OSError):  # pragma: no cover - unexpected
        return []
    paths = []
    for i in range(count()):
        raw = name(i)
        if not raw:
            continue
        path = raw.decode("utf-8", "replace")
        if path not in paths:
            paths.append(path)
    return paths


def sharedObjectPaths():
    """Return the paths of all shared objects loaded into this process.

    Uses ``/proc/self/maps`` on Linux and ``_dyld_image_count`` /
    ``_dyld_get_image_name`` on macOS, the Python pendant of
    ``dl_iterate_phdr``. An empty list means the platform offers no
    enumeration, the callers then fall back to the C++ core.
    """
    global _warnedPlatform
    if sys.platform.startswith("linux"):
        return _pathsLinux()
    if sys.platform == "darwin":
        return _pathsDarwin()
    if not _warnedPlatform:
        _warnedPlatform = True
        _warn(
            "cannot enumerate loaded libraries on this platform, "
            "thread count for BLAS is limited to the C++ core instance"
        )
    return []


###################
# symbol handling
###################
def _setterSpec(name):
    """Return ``(ctype, passByRef)`` for a BLAS thread setter.

    Names with a trailing underscore are Fortran aliases and take a pointer
    to the thread count instead of the value itself. :data:`VALUE_SETTERS`
    lists the exceptions, they take the count by value although their name
    has the Fortran suffix.
    """
    if name in VALUE_SETTERS:
        return ctypes.c_longlong, False
    if name.endswith("64_"):
        return ctypes.c_longlong, True
    if name.endswith("_"):
        return ctypes.c_int, True
    if name == "bli_thread_set_num_threads":
        return ctypes.c_longlong, False
    return ctypes.c_int, False


def _getterSpec(name):
    """Return the C return type for a BLAS thread getter."""
    if name.endswith("64_"):
        return ctypes.c_longlong
    return ctypes.c_int


def _handle(path):
    """Return a cached ctypes handle for *path* (None on failure)."""
    if path in _handles:
        return _handles[path]
    handle = None
    try:
        mode = getattr(os, "RTLD_LOCAL", 0)
        handle = ctypes.CDLL(path, mode=mode)
    except (OSError, TypeError):
        handle = None
    _handles[path] = handle
    return handle


def _bind(handle, name, argTypes, restype=ctypes.c_int):
    """Bind ``handle.<name>`` with explicit signature, None if missing."""
    try:
        func = getattr(handle, name)
    except AttributeError:
        return None
    if func is None:
        return None
    func.argtypes = list(argTypes)
    func.restype = restype
    return func


##############
# runtime types
##############
class BlasRuntime:
    """A BLAS thread runtime (OpenBLAS, BLIS or MKL) in this process."""

    __slots__ = (
        "path",
        "kind",
        "setterNames",
        "getterName",
        "_setters",
        "_getter",
    )

    def __init__(self, path, kind, getterName, setters, getter):
        self.path = path
        self.kind = kind
        self.getterName = getterName
        self.setterNames = tuple(name for name, _, _, _ in setters)
        self._setters = setters
        self._getter = getter

    def setNumThreads(self, nThreads):
        """Apply *nThreads* to this runtime."""
        for _, func, argType, byRef in self._setters:
            if byRef:
                # keep the value alive until the call returned
                value = argType(nThreads)
                func(ctypes.byref(value))
            else:
                func(nThreads)

    def numThreads(self):
        """Return the current thread count or None if unknown."""
        if self._getter is None:
            return None
        try:
            value = int(self._getter())
        except Exception:  # pragma: no cover - defensive
            return None
        return value

    def __repr__(self):
        return f"<BlasRuntime {self.kind} {self.path}={self.numThreads()}>"


class OmpRuntime:
    """An OpenMP runtime (libgomp, libomp, libiomp5) in this process."""

    __slots__ = ("path", "_functions")

    def __init__(self, path, functions):
        self.path = path
        self._functions = functions

    def _call(self, name, default=None):
        func = self._functions.get(name)
        if func is None:
            return default
        try:
            value = int(func())
        except Exception:  # pragma: no cover - defensive
            return default
        return value

    def setNumThreads(self, nThreads):
        func = self._functions.get("omp_set_num_threads")
        if func is not None:
            func(nThreads)

    def maxThreads(self):
        return self._call("omp_get_max_threads")

    def threadLimit(self):
        return self._call("omp_get_thread_limit")

    def numProcs(self):
        return self._call("omp_get_num_procs")

    def inParallel(self):
        return self._call("omp_in_parallel")

    def numThreads(self):
        return self._call("omp_get_num_threads")

    def __repr__(self):
        return f"<OmpRuntime {self.path} max_threads={self.maxThreads()}>"


##################
# probing
##################
def _probeBlas(path, handle=None):
    """Return the :class:`BlasRuntime` of *path*, None if it is no BLAS.

    *handle* replaces the real ``ctypes`` handle, the tests use that to
    inject a stub object instead of a loaded library.
    """
    if handle is None:
        handle = _handle(path)
    if handle is None:
        return None
    for kind, setterNames, getterNames in BLAS_FAMILIES:
        setters = []
        for name in setterNames:
            argType, byRef = _setterSpec(name)
            argTypes = [ctypes.POINTER(argType)] if byRef else [argType]
            func = _bind(handle, name, argTypes)
            if func is not None:
                setters.append((name, func, argType, byRef))
        if not setters:
            continue
        getterName = None
        getter = None
        for name in getterNames:
            getter = _bind(handle, name, [], _getterSpec(name))
            if getter is not None:
                getterName = name
                break
        return BlasRuntime(path, kind, getterName, setters, getter)
    return None


def _probeOmp(path):
    handle = _handle(path)
    if handle is None:
        return None
    functions = {}
    for name in OMP_SYMBOLS:
        argTypes = [ctypes.c_int] if name == "omp_set_num_threads" else []
        func = _bind(handle, name, argTypes)
        if func is not None:
            functions[name] = func
    if "omp_set_num_threads" not in functions:
        return None
    return OmpRuntime(path, functions)


def _sortKey(path):
    """Probe name matched objects first, all others afterwards."""
    base = Path(path).name
    return (0 if BLAS_NAME_PATTERN.search(base) else 1, path)


def _refreshBlas():
    """Probe newly loaded objects, return all known BLAS runtimes."""
    found = []
    for path in sorted(sharedObjectPaths(), key=_sortKey):
        if not BLAS_OBJECT_PATTERN.search(Path(path).name):
            continue
        if path in _probedBlas:
            continue
        _probedBlas.add(path)
        runtime = _probeBlas(path)
        if runtime is not None:
            found.append(runtime)
    if found:
        if _lastBlasTarget is not None:
            for runtime in found:
                runtime.setNumThreads(_lastBlasTarget)
        _blasRuntimes.extend(found)
    return _blasRuntimes


def _refreshOmp():
    """Probe newly loaded objects, return all known OMP runtimes."""
    found = []
    for path in sharedObjectPaths():
        if path in _probedOmp:
            continue
        _probedOmp.add(path)
        if not OMP_NAME_PATTERN.search(Path(path).name):
            continue
        runtime = _probeOmp(path)
        if runtime is not None:
            found.append(runtime)
    if found:
        if _lastOmpTarget is not None:
            for runtime in found:
                runtime.setNumThreads(_lastOmpTarget)
        _ompRuntimes.extend(found)
    return _ompRuntimes


##################
# public API
##################
def blasRuntimes():
    """Return all BLAS thread runtimes found in this process."""
    with _lock:
        return list(_refreshBlas())


def ompRuntimes():
    """Return all OpenMP runtimes found in this process."""
    with _lock:
        return list(_refreshOmp())


def setBLASThreadCount(nThreads):
    """Limit the threads of all BLAS libraries loaded in this process.

    Priority
    --------
    1. The last explicit call of this function (or of :func:`setThreadCount`),
       it is also applied to libraries that get loaded later
    2. ``OPENBLAS_NUM_THREADS``, ``MKL_NUM_THREADS``, ``BLIS_NUM_THREADS``,
       read once by :func:`applyEnvironment` while importing
    3. pygimli's own default: ``__setTC__()`` in ``core/src/gimli.cpp`` runs
       whenever ``libgimli`` is loaded and calls
       ``setBLASThreadCount(min(8, numberOfCPU() - 2))``
    4. the default of the BLAS library itself

    ``__setTC__()`` is a static initialiser and ends with an unconditional
    ``setThreadCount()``. ``OPENBLAS_NUM_THREADS`` therefore only survives
    because :func:`applyEnvironment` re-applies it after the import.
    ``OMP_NUM_THREADS`` survives anyway, ``setThreadCount()`` touches BLAS
    only.

    The default of step 3 is also what :func:`applyEnvironment` applies to
    every other instance and exports to the environment, so a plain
    ``import pygimli`` never uses more than ``min(8, ncpu - 2)`` threads.

    Parameters
    ----------
    nThreads : int
        Maximum amount of threads, values below 1 are clipped to 1.
    """
    global _lastBlasTarget, _warnedBlas
    nThreads = max(1, int(nThreads))
    with _lock:
        _lastBlasTarget = nThreads
        runtimes = _refreshBlas()
        for runtime in runtimes:
            runtime.setNumThreads(nThreads)
        if not runtimes and not _warnedBlas:
            _warnedBlas = True
            _warn(
                "no BLAS runtime found in this process, "
                f"cannot set the BLAS thread count to {nThreads}"
            )


def getBLASThreadCount():
    """Return the largest thread count of all BLAS runtimes (0 if none)."""
    with _lock:
        runtimes = _refreshBlas()
        counts = [
            n for n in (rt.numThreads() for rt in runtimes) if n is not None
        ]
    if counts:
        return max(counts)
    value = _coreCall("getBLASThreadCount")
    return 0 if value is None else int(value)


def setOMPThreadCount(nThreads):
    """Limit the threads of all OpenMP runtimes in this process."""
    global _lastOmpTarget, _warnedOmp
    nThreads = max(1, int(nThreads))
    _coreCall("setOMPThreadCount", nThreads)
    with _lock:
        _lastOmpTarget = nThreads
        runtimes = _refreshOmp()
        for runtime in runtimes:
            runtime.setNumThreads(nThreads)
        if not runtimes and not _warnedOmp:
            _warnedOmp = True
            _warn(
                "no OpenMP runtime found in this process, "
                f"cannot set the OMP thread count to {nThreads}"
            )


def getOMPThreadCount():
    """Return the OpenMP thread limit of the process.

    Several OpenMP runtimes can be loaded at once (libgomp for the C++ core,
    libomp for a wheel). The largest ``omp_get_max_threads()`` is returned,
    that is the one that actually limits the slowest pool.
    """
    with _lock:
        runtimes = _refreshOmp()
        counts = [n for n in (rt.maxThreads() for rt in runtimes) if n]
    if counts:
        return max(counts)
    value = _coreCall("getOMPThreadCount")
    return 1 if value is None else max(1, int(value))


def setThreadCount(nThreads):
    """Set the OpenMP *and* the BLAS thread count of the whole process."""
    global _threadCount
    nThreads = max(1, int(nThreads))
    _coreCall("setThreadCount", nThreads)
    setOMPThreadCount(nThreads)
    setBLASThreadCount(nThreads)
    _threadCount = nThreads


def threadCount():
    """Return the thread count requested from gimli."""
    if _threadCount is not None:
        return _threadCount
    value = _coreCall("threadCount")
    return 0 if value is None else int(value)


def _openblasConfig(runtime):
    """Return the ``openblas_get_config()`` string of *runtime* or None."""
    handle = _handle(runtime.path)
    if handle is None:
        return None
    func = _bind(handle, "openblas_get_config", [], ctypes.c_char_p)
    if func is None:
        return None
    try:
        config = func()
    except Exception:  # pragma: no cover - defensive
        return None
    if not config:
        return None
    return config.decode("utf-8", "replace")


def _implementationKey(runtime):
    """Return a key identifying the thread pool behind *runtime*.

    FlexiBLAS is loaded as a front end plus one plugin per kernel, all of
    them resolve ``openblas_get_num_threads()`` to the same code address.
    Two runtimes with the same address are the same pool.
    """
    func = runtime._getter
    if func is not None:
        try:
            return ("code", ctypes.cast(func, ctypes.c_void_p).value)
        except (TypeError, ValueError):  # pragma: no cover - stub in tests
            pass
    return ("path", runtime.path)


def _uniqueRuntimes(runtimes):
    """Collapse runtimes that resolve to the same implementation."""
    unique = {}
    for runtime in runtimes:
        key = _implementationKey(runtime)
        known = unique.get(key)
        if known is None or len(runtime.path) < len(known.path):
            unique[key] = runtime
    return list(unique.values())


def threadsInfo():
    """Log the current thread configuration of the whole process."""
    with _lock:
        omp = list(_refreshOmp())
        blas = list(_refreshBlas())

    if omp:
        runtime = omp[0]
        _logInfo(f"omp  max_threads  : {runtime.maxThreads()}")
        _logInfo(f"omp  thread_limit : {runtime.threadLimit()}")
        _logInfo(f"omp  num_procs    : {runtime.numProcs()}")
        _logInfo(f"omp  in_parallel  : {runtime.inParallel()}")
        _logInfo(f"omp  runtime      : {runtime.path}")
    else:
        _logInfo(f"omp  max_threads  : {getOMPThreadCount()} (C++ core)")

    if not blas:
        _logInfo("blas no runtime found")
        return
    count = getBLASThreadCount()
    for runtime in _uniqueRuntimes(blas):
        threads = runtime.numThreads()
        if threads is None:
            threads = count
        getter = runtime.getterName or "unknown"
        _logInfo(f"blas {runtime.kind} ({getter}): {threads}   {runtime.path}")
        config = _openblasConfig(runtime)
        if config is not None and "USE_OPENMP" in config:
            _warn(
                f"{runtime.path} is built with OpenMP, "
                "use omp_set_num_threads()/OMP_NUM_THREADS there"
            )


#########################
# environment handling
#########################
#: Variables read by the BLAS libraries while they are loaded.
BLAS_ENV_VARS = ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")
#: Variables read by the OpenMP runtimes while they are loaded.
OMP_ENV_VARS = ("OMP_NUM_THREADS",)


def _envInt(*names):
    """Return the first integer found in ``os.environ[names]``."""
    for name in names:
        value = os.environ.get(name)
        if value is None or value == "":
            continue
        try:
            return int(value)
        except ValueError:
            _warn(f"ignoring non integer value {name}={value!r}")
    return None


def _exportDefaults(names, value):
    """Export *value* for libraries that get loaded later.

    numpy and scipy read ``OPENBLAS_NUM_THREADS`` while their private
    OpenBLAS copy is dlopen'ed, which usually happens at the first matrix
    product, long after this module was imported. Only unset variables are
    written, a value chosen by the user always wins. Child processes inherit
    the defaults.
    """
    for name in names:
        if not os.environ.get(name):
            os.environ[name] = str(value)


def applyEnvironment():
    """Apply the thread counts of the process.

    ``OPENBLAS_NUM_THREADS``, ``MKL_NUM_THREADS``, ``BLIS_NUM_THREADS`` and
    ``OMP_NUM_THREADS`` win over ``GIMLI_NUM_THREADS``, which wins over the
    default that ``__setTC__()`` determined while ``libgimli`` was loaded
    (``min(8, ncpu - 2)``). That default is applied to every instance and
    exported to the environment, so libraries that are dlopen'ed later pick
    it up as well -- more threads have to be requested explicitly. A value
    set by the user is never overwritten and an explicit
    :func:`setBLASThreadCount`/:func:`setOMPThreadCount` call afterwards
    replaces the result.
    """
    gimli = _envInt("GIMLI_NUM_THREADS")
    blas = _envInt(*BLAS_ENV_VARS)
    omp = _envInt(*OMP_ENV_VARS)
    if blas is None:
        blas = gimli
    if omp is None:
        omp = gimli
    if blas is None:
        # value that __setTC__() gave the instance libgimli is linked with
        blas = _coreCall("getBLASThreadCount") or None
    if omp is None:
        omp = _coreCall("threadCount") or None
    if blas is not None:
        _exportDefaults(BLAS_ENV_VARS, blas)
        setBLASThreadCount(blas)
    if omp is not None:
        _exportDefaults(OMP_ENV_VARS, omp)
        setOMPThreadCount(omp)


applyEnvironment()
