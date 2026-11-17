"""Demo of the pyGIMLi logging facilities.

Run with ``python pygimli/testing/test_Logger.py``. Importing this module must
not log or change the global log settings, because the module is imported by
the API documentation (autosummary) and by ``pg.test``.
"""

import pygimli as pg


def main():
    """Exercise the pyGIMLi log levels, decorators and colored output."""
    pg.version()

    # test pygimli log
    pg.info("Start numeric log test." + str(pg.log(pg.Vector(1, 1.0))))
    pg.setVerbose(True)
    pg.verbose("some verbose notes")
    pg.warn("Start warning test.")

    def testTraceback1():
        """Log an error from a nested call to show the call location."""

        def testTraceback2():
            """Do the actual logging."""
            pg.error("Start error test.: int", 1, " vec", pg.Vector(2))

        testTraceback2()

    testTraceback1()

    @pg.v
    def testVerboseDecorator1():
        """Be verbose even if verbose is switched off."""
        pg.verbose(
            "testVerboseDecorator1 should be seen even if verbose is false"
        )

    def testVerboseDecorator2():
        """Be verbose only if verbose is switched on."""
        pg.verbose(
            "testVerboseDecorator2 should be seen even if verbose is true"
        )

    pg.setVerbose(False)
    testVerboseDecorator1()

    pg.setVerbose(False)
    testVerboseDecorator2()
    pg.setVerbose(True)
    testVerboseDecorator2()

    @pg.d
    def testDebugDecorator():
        """Be debug even if debug is switched off."""
        pg.debug(
            "testDebugDecorator should always be seen even if debug is false"
        )

    testDebugDecorator()

    pg.d(testDebugDecorator())

    pg.debug("debug 0")

    pg.setDebug(1)
    pg.debug("debug ON")
    pg.setThreadCount(2)

    pg.debug("debug with a trace", withTrace=True)

    # should not be printed
    pg.setDebug(0)
    pg.debug("debug OFF")
    pg.setThreadCount(2)

    # test core log (should not be used outside the core)
    pg.log(pg.core.Info, "core log ")
    pg.log(pg.core.Warning, "core log ")
    pg.log(pg.core.Error, "core log ")
    # a core Critical only logs and never raises (see GIMLI::log), it is
    # therefore not demonstrated here: it would abort the documentation build
    # (see doc/_sphinx-ext/criticalguard.py)

    def testMethod(**kwargs):
        """Warn about an empty keyword-argument call."""
        pg.core.logger.warnNonEmptyArgs(kwargs)

    testMethod(a=1, b="foo")

    # teste colored output
    print(pg._("Green", c="g"), pg._("Red", c="r"), pg._("Yellow", c="y"))
    print(pg._("more", "then", "one", c="6;30;42"))


if __name__ == "__main__":
    main()
