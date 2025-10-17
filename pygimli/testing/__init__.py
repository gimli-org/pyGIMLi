#!/usr/bin/env python3
"""Testing utilities.

In Python you can call `pygimli.test()` to run all docstring
examples.

Writing tests for pyGIMLi
-------------------------

Please check: https://docs.pytest.org/en/latest/
"""
import sys
from os.path import join, realpath

import numpy as np
import pygimli as pg
import warnings

__devTests__ = True


def setDevTests(mode):
    """Set pygimli testing mode.

    Testing mode ensures a constant seed for the random generator if you use
    pg.randn().
    """
    global __devTests__
    __devTests__ = mode


def devTests():
    """Determine pygimli testing mode.

    Returns True if pygimli is in testing mode.
    """
    import os
    if os.getenv('DEVTESTS') == '1':
        return True
    if os.getenv('DEVTESTS') == '0':
        return False

    global __devTests__
    return __devTests__


def test(target=None, show=False, onlydoctests=False, coverage=False,
         htmlreport=False, abort=False, verbose=True, devTests=False,
         exitOnError=False):
    """Run docstring examples and additional tests.

    Examples
    --------
    >>> import pygimli as pg
    >>> # Run the whole test suite.
    >>> pg.test() # doctest: +SKIP
    >>> # Test a single function by a string.
    >>> pg.test("utils.boxprint", verbose=False) # doctest: +SKIP
    >>> # The target argument can also be the function directly
    >>> from pygimli.utils import boxprint
    >>> pg.test(boxprint, verbose=False) # doctest: +SKIP
    >>> # Use some logical expressions
    >>> pg.test("draw and not drawMesh") # doctest: +SKIP

    Parameters
    ----------
    target : function or string or pattern (-k flag in pytest), optional
        Function or method to test. By default everything is tested.
    show : boolean, optional
        Show viewer windows during test run. They will be closed
        automatically.
    onlydoctests : boolean, optional
        Run test files in testing as well.
    coverage : boolean, optional
        Create a coverage report. Requires the pytest-cov plugin.
    htmlreport : str, optional
        Filename for HTML report such as www.pygimli.org/build_tests.html.
        Requires pytest-html plugin.
    exitOnError : boolean, optional
        Run sys.exit(1) if pytest detects a test error.
    abort : boolean, optional
        Deprecated alias for exitOnError.
    devTests: boolean[False]
        Don't skip special tests marked for development, only with the
        @pg.skipOnDefaultTest decorator. Can be overwritten by env DEVTESTS.

    Returns
    -------
    exitcode: integer
        Return code of pytest, usually 0 for a successful test run.

    """
    setDevTests(devTests)

    try:
        import pytest
    except ImportError:
        raise ImportError("pytest is required to run test suite. "
                          "Try 'pip install pytest'.")

    exitOnError = exitOnError | abort

    # Remove figure warnings
    np.random.seed(1337)
    plt = pg.plt
    plt.rcParams["figure.max_open_warning"] = 1000
    warnings.filterwarnings("ignore", category=UserWarning,
                            message='Matplotlib is currently using agg, a '
                                    'non-GUI backend, so cannot show figure.')

    printopt = np.get_printoptions()


    # Numpy compatibility (array string representation has changed)
    if np.__version__[:4] == "1.14":
        pg.warn("Some doctests will fail due to old numpy version.",
                "Consider upgrading to numpy >= 1.15")

    old_backend = plt.get_backend()
    old_backend_pv = pg.rc["pyvista.backend"]


    if not show:
        import pyvista
        pyvista.OFF_SCREEN=True

        plt.switch_backend("Agg")
        pg.rc["pyvista.backend"] = None
    else:
        plt.ion()

    cwd = join(realpath(__path__[0]), '..')

    excluded = [
        "gui", "physics/traveltime/example.py", "physics/em/fdemexample.py"
    ]

    if onlydoctests:
        excluded.append("testing")

    cmd = (["--color", "yes", "--doctest-modules", "-p", "no:warnings"])

    string = f"pygimli {pg.__version__}"

    target_source = False
    if target:
        if not isinstance(target, str):
            import inspect
            target_source = inspect.getsourcefile(target)
            target = target.__name__
        else:
            target = target.replace("pg.", "")
            target = target.replace("pygimli.", "")

        cmd.extend(["-k", target, "--no-header", "--doctest-report", "udiff"])
        if not verbose:
            cmd.extend(["-qq", "-rN"])

        if show:  # Keep figure opened if single function is tested
            plt.ioff()

        string = f"'{target}' from {string}"

    if verbose:
        cmd.extend(["-v", "--durations", "5"])
        pg.boxprint(f"Testing {string}", sym="+", width=90)


    if coverage:
        pc = pg.optImport("pytest_cov", "create a code coverage report")
        if pc:
            cmd.extend(["--cov", "pygimli"])
            cmd.extend(["--cov-report", "term"])

    if htmlreport:
        ph = pg.optImport("pytest_html", "create a html report")
        if ph:
            cmd.extend(["--html", htmlreport])

    for directory in excluded:
        cmd.extend(["--ignore", join(cwd, directory)])

    plt.close("all")
    if target_source:
        cmd.extend([target_source])
    else:
        cmd.extend([cwd])

    exitcode = pytest.main(cmd)

    plt.switch_backend(old_backend)
    pg.rc["pyvista.backend"] = old_backend_pv

    np.set_printoptions(**printopt)

    # if abort is True:
    #     if verbose:
    #         print("Exiting with exitcode", exitcode)
    #     if exitcode != 0:
    #         pg.critical("some tests failed, exiting with exitcode", exitcode)
    if exitOnError == pytest.ExitCode.TESTS_FAILED:
        if verbose:
            print("Exiting with exitcode", exitcode)
        sys.exit(exitcode)

    return exitcode
