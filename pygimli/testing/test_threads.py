#!/usr/bin/env python
"""Regression tests for the process wide thread control.

Bug: ``pg.setThreadCount(1)`` / ``pg.core.setBLASThreadCount(1)`` only
affected the OpenBLAS copy that ``libgimli`` is linked against. A script that
imported numpy *before* pygimli kept using all cores, because numpy brings its
own (``RTLD_LOCAL``) OpenBLAS instance that was never touched.
``OPENBLAS_NUM_THREADS=1`` in the environment worked around it.

The tests import numpy first on purpose: that is the worst case reported.
"""

import os
import subprocess
import sys
import time
import unittest
from unittest import mock

import numpy  # noqa: F401  # must stay in front of pygimli (see bug above)

import pygimli as pg
from pygimli.core import threads

#: Variables that have to be cleared for a truly bare child process.
THREAD_VARS = (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "OMP_NUM_THREADS",
)


def _hasBlas():
    return len(pg.core.blasRuntimes()) > 0


def _runPython(code, **env):
    """Run *code* in a fresh interpreter, return the ``RESULT`` fields.

    A separate process is needed because the thread counts are fixed while
    the libraries are loaded, and because this module exports its defaults
    into ``os.environ``.
    """
    environment = dict(os.environ)
    for name in THREAD_VARS:
        environment.pop(name, None)
    for name, value in env.items():
        environment[name] = str(value)
    result = subprocess.run(  # noqa: S603 - fixed code, no user input
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        check=True,
        timeout=300,
        env=environment,
    )
    lines = [
        line for line in result.stdout.splitlines() if line.startswith("RESULT")
    ]
    return lines[-1].split()[1:]


def _cpuWallRatio(size=1500):
    """Return cpu time divided by wall time of one matrix product."""
    values = []
    for _ in range(2):
        a = numpy.random.random((size, size))
        cpu = time.process_time()
        wall = time.perf_counter()
        a @ a
        values.append(
            (time.process_time() - cpu) / (time.perf_counter() - wall)
        )
    return min(values)


class TestThreadControl(unittest.TestCase):
    """Thread counts have to reach all runtimes, not only the core one."""

    def setUp(self):
        """Remember the thread counts so tearDown can restore them."""
        self.oldBlas = pg.core.getBLASThreadCount()
        self.oldOmp = pg.core.getOMPThreadCount()

    def tearDown(self):
        """Restore the thread counts of the test process."""
        pg.core.setBLASThreadCount(max(1, self.oldBlas))
        pg.core.setOMPThreadCount(max(1, self.oldOmp))

    def test_setBLASThreadCount(self):
        """Bug: BLAS thread count was a silent no-op / hit the wrong lib."""
        if not _hasBlas():
            self.skipTest("no BLAS runtime in this process")
        pg.core.setBLASThreadCount(1)
        self.assertEqual(pg.core.getBLASThreadCount(), 1)
        for runtime in pg.core.blasRuntimes():
            self.assertEqual(runtime.numThreads(), 1, runtime.path)

    def test_setBLASThreadCountIsIdempotent(self):
        """Repeated calls keep the requested value."""
        if not _hasBlas():
            self.skipTest("no BLAS runtime in this process")
        for _ in range(3):
            pg.core.setBLASThreadCount(1)
            self.assertEqual(pg.core.getBLASThreadCount(), 1)

    def test_setOMPThreadCount(self):
        """Bug: getOMPThreadCount returned omp_get_num_threads() -> 1."""
        pg.core.setOMPThreadCount(2)
        self.assertEqual(pg.core.getOMPThreadCount(), 2)
        pg.core.setOMPThreadCount(1)
        self.assertEqual(pg.core.getOMPThreadCount(), 1)

    def test_setThreadCountSetsOmpAndBlas(self):
        """gimli.h documents setThreadCount as OMP *and* BLAS."""
        pg.core.setThreadCount(2)
        self.assertEqual(pg.core.threadCount(), 2)
        self.assertEqual(pg.core.getOMPThreadCount(), 2)
        if _hasBlas():
            self.assertEqual(pg.core.getBLASThreadCount(), 2)

    def test_envGimliNumThreadsIsRespected(self):
        """Bug: GIMLI_NUM_THREADS was overwritten by the hard coded default."""
        with mock.patch.dict(
            os.environ, {"GIMLI_NUM_THREADS": "1"}, clear=False
        ):
            for name in THREAD_VARS:
                os.environ.pop(name, None)
            pg.core.applyEnvironment()
        self.assertEqual(pg.core.getOMPThreadCount(), 1)
        if _hasBlas():
            self.assertEqual(pg.core.getBLASThreadCount(), 1)

    def test_envBlasNumThreadsIsAuthoritative(self):
        """OPENBLAS_NUM_THREADS wins over GIMLI_NUM_THREADS."""
        if not _hasBlas():
            self.skipTest("no BLAS runtime in this process")
        with mock.patch.dict(
            os.environ,
            {"GIMLI_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "3"},
            clear=False,
        ):
            pg.core.applyEnvironment()
        self.assertEqual(pg.core.getBLASThreadCount(), 3)

    def test_runtimeLoadedAfterSettingIsCovered(self):
        """A BLAS lib imported later (e.g. scipy) gets the value as well."""
        if not _hasBlas():
            self.skipTest("no BLAS runtime in this process")
        pg.core.setBLASThreadCount(1)
        try:
            import scipy.linalg  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed")
        self.assertEqual(pg.core.getBLASThreadCount(), 1)
        for runtime in pg.core.blasRuntimes():
            self.assertEqual(runtime.numThreads(), 1, runtime.path)

    def test_scipyOpenblasSetterIsCalledByValue(self):
        """Regression: scipy_openblas_set_num_threads64_ takes int64 by value.

        The ``64_`` suffix normally marks a Fortran alias that wants a
        pointer, so the setter was called with ``byref(c_longlong(1))``.
        That is accepted but ignored, the pool kept using every core.
        """
        calls = []
        fake = mock.Mock()
        # a plain function: _bind() sets argtypes/restype on it
        fake.scipy_openblas_set_num_threads64_ = lambda value: calls.append(
            value
        )
        for name in (
            "openblas_set_num_threads64_",
            "openblas_set_num_threads",
            "openblas_set_num_threads_",
            "bli_thread_set_num_threads",
            "mkl_set_num_threads",
            "mkl_set_num_threads_local",
            "scipy_openblas_set_num_threads",
            "scipy_openblas_get_num_threads64_",
        ):
            setattr(fake, name, None)
        fake.scipy_openblas_get_num_threads64_ = lambda: 1
        runtime = threads._probeBlas("<fake>", handle=fake)
        self.assertIsNotNone(runtime)
        self.assertEqual(runtime.kind, "scipy_openblas")
        runtime.setNumThreads(1)
        self.assertEqual(calls, [1])

    def test_blasThreadCountHasEffect(self):
        """Regression: 30 threads despite setThreadCount(1).

        Only the real effect counts: limiting the BLAS pool to one thread
        has to drop the cpu/wall ratio of a matrix product.
        """
        if os.cpu_count() < 4:
            self.skipTest("needs a multi core machine")
        if not _hasBlas():
            self.skipTest("no BLAS runtime in this process")
        old = _cpuWallRatio()
        pg.core.setBLASThreadCount(1)
        new = _cpuWallRatio()
        pg.core.setBLASThreadCount(max(1, int(round(old))))
        self.assertLess(
            new,
            max(2.0, old / 2),
            f"BLAS ignores the thread count (cpu/wall {old:.1f} -> {new:.1f})",
        )

    def test_blasRuntimesSkipForeignObjects(self):
        """Regression: libgimli/cholmod/umfpack were reported as BLAS.

        All three re-export the OpenBLAS symbols through their dependency
        chain, they must not show up as a BLAS instance of their own.
        """
        names = [runtime.path for runtime in pg.core.blasRuntimes()]
        for needle in ("libgimli", "_pygimli_", "cholmod", "umfpack"):
            self.assertFalse(
                [n for n in names if needle in n],
                f"{needle} reported as BLAS: {names}",
            )

    def test_plainRunUsesTheCoreDefault(self):
        """Regression: numpy's own OpenBLAS ignored the pygimli default.

        ``__setTC__()`` caps the instance libgimli is linked with, but the
        wheel brings a private copy that is dlopen'ed at the first matrix
        product. It has to start with the same default, otherwise a plain
        ``import pygimli`` still uses every core.
        """
        core, blas = (
            int(value)
            for value in _runPython(
                "import numpy, pygimli as pg\n"
                "a = numpy.random.random((1200, 1200))\n"
                "a @ a\n"
                "print('RESULT', pg.core.pgcore.threadCount(),\n"
                "      pg.core.getBLASThreadCount())\n"
            )
        )
        self.assertGreater(core, 0)
        self.assertLessEqual(blas, core)

    def test_gimliNumThreadsReachesTheCore(self):
        """Regression: __setTC__() discarded GIMLI_NUM_THREADS."""
        value = _runPython(
            "import pygimli as pg\n"
            "print('RESULT', pg.core.pgcore.threadCount())\n",
            GIMLI_NUM_THREADS=1,
        )
        self.assertEqual(int(value[0]), 1)

    def test_exportedDefaultKeepsTheUserValue(self):
        """The exported default must not replace an explicit setting."""
        value = _runPython(
            "import os, pygimli as pg\n"
            "print('RESULT', os.environ['OPENBLAS_NUM_THREADS'])\n",
            OPENBLAS_NUM_THREADS=3,
        )
        self.assertEqual(int(value[0]), 3)

    def test_threadsInfo(self):
        """Check that threadsInfo reports the filtered runtimes only."""
        a = numpy.random.random((200, 200))
        a @ a
        messages = []
        with mock.patch.object(threads, "_logInfo", messages.append):
            pg.core.threadsInfo()
        self.assertTrue(messages)
        for message in messages:
            for needle in ("libgimli", "_pygimli_", "cholmod", "umfpack"):
                self.assertNotIn(needle, message)
        blasLines = [m for m in messages if m.startswith("blas ")]
        if _hasBlas():
            self.assertTrue(blasLines, messages)


if __name__ == "__main__":
    unittest.main()
