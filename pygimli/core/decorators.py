"""Useful utility decorators."""

import functools
import sys
import unittest


class renamed:  # noqa: N801  public API name (pg.renamed)
    """Mark functions that are renamed and will be removed later.

    ```
    @pg.renamed(newname, "1.2")
    def oldname(args, kwargs):
        pass
    ```
    """

    def __init__(self, newFunc, removed=""):
        self.newFunc = newFunc
        self.removed = removed

    def __call__(self, func):
        """Return the wrapper that forwards to the new function."""

        def wrapper(*args, **kwargs):
            import pygimli as pg

            pg.warning(
                func.__name__
                + " had been renamed to "
                + self.newFunc.__name__
                + " and will be removed in: "
                + self.removed
            )
            return self.newFunc(*args, **kwargs)

        return wrapper


class deprecate:  # noqa: N801  public API name (pg.deprecate)
    """Mark a function as deprecated.

    ```
    @pg.deprecate(newname, 1.5)
    def test(args, kwargs):
        pass
    ```
    """

    def __init__(self, newFunc, removedIn=""):
        self.newFunc = newFunc
        self.removed = str(removedIn)

    def __call__(self, func):
        """Return the wrapper that warns and forwards to the new function."""

        def wrapper(*args, **kwargs):
            import pygimli as pg

            pg.deprecated(
                f"{func.__name__}({args[1:]}, {kwargs})"
                + " had been marked as deprecated"
                + " and will be removed in: "
                + self.removed,
                hint=f"\nUse: {self.newFunc}",
            )

            return self.newFunc(*args, **kwargs)

        return wrapper


def skipOnDefaultTest(func):
    """Skip a test in the default test run.

    The test is reported as *skipped* (not silently passed) and runs only with
    ``pg.test(devTests=True)`` or ``DEVTESTS=1``.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        from pygimli.testing import devTests

        if devTests() is True:
            return func(*args, **kwargs)

        raise unittest.SkipTest(f"development test {func.__name__}")

    return wrapper


def singleton(cls):
    """Make a class a Singleton class (only one instance)."""

    @functools.wraps(cls)
    def wrapper(*args, **kwargs):
        if wrapper.instance is None:
            wrapper.instance = cls(*args, **kwargs)
        return wrapper.instance

    wrapper.instance = None
    return wrapper


# Lazy evaluation of expensive modules (eg. mpl.pylab)
# found: https://stackoverflow.com/questions/880530/can-modules-have-properties-the-same-way-that-objects-can
class AttrGetter:
    """Forward attribute access to a lazy module getter."""

    def __new__(cls, gt):
        """Return an instance around the getter *gt*."""
        if isinstance(gt, cls):
            return gt
        else:
            o = super().__new__(cls)
            o.oldgetattr = gt
            o.funcmap = {}
            return o

    def __call__(self, name):
        """Return the lazily evaluated function or the original attribute."""
        name2 = "_" + name
        if name2 in self.funcmap:
            return self.funcmap[name2]()
        else:
            return self.oldgetattr(name)

    def add(self, func):
        """Register *func* for lazy evaluation."""
        self.funcmap[func.__name__] = func


def moduleProperty(func):
    """Turn module functions into a property.

    Function names must be prefixed with an underscore.
    """
    module = sys.modules[func.__module__]

    def base_getattr(name):
        raise AttributeError(
            f"module '{module.__name__}' has no attribute '{name}'"
        )

    ag = AttrGetter(getattr(module, "__getattr__", base_getattr))
    module.__getattr__ = ag
    ag.add(func)
    return func
