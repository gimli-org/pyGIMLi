
# -*- coding: utf-8 -*-
"""
Useful utility decorators.
"""
import sys


class renamed:
    """Decorator to mark functions are renamed and will be removed some later
    ```
    @pg.renamed(newname, '1.2')
    def oldname(args, kwargs):
        pass
    ```
    """
    def __init__(self, newFunc, removed=''):
        self.newFunc = newFunc
        self.removed = removed

    def __call__(self, func):
        def wrapper(*args, **kwargs):
            import pygimli as pg
            pg.warning(func.__name__ + ' had been renamed to ' + \
                       self.newFunc.__name__ + \
                       ' and will be removed in: ' + self.removed)
            return self.newFunc(*args, **kwargs)
        return wrapper


def skipOnDefaultTest(func):
    """Decorator to mark a test to be skipped with default testing
    ```
    @pg.skipOnDefaultTest()
    def test(args, kwargs):
        pass
    ```
    """
    def wrapper(*args, **kwargs):
        import pygimli as pg
        from pygimli.testing import devTests

        if devTests() == True:
            return func(*args, **kwargs)
        else:
            pg.info('Skipped test:', func)

    return wrapper


import functools
def singleton(cls):
    """Make a class a Singleton class (only one instance)"""
    @functools.wraps(cls)
    def wrapper(*args, **kwargs):
        if wrapper.instance is None:
            wrapper.instance = cls(*args, **kwargs)
        return wrapper.instance
    wrapper.instance = None
    return wrapper


# Lazy evalation of expensive modules (eg. mpl.pylab)
# found: https://stackoverflow.com/questions/880530/can-modules-have-properties-the-same-way-that-objects-can
class AttrGetter:
    def __new__(cls, gt):
        if isinstance(gt, cls):
            return gt
        else:
            o = super().__new__(cls)
            o.oldgetattr = gt
            o.funcmap = {}
            return o

    def __call__(self, name):
        name2 = "_" + name
        if name2 in self.funcmap:
            return self.funcmap[name2]()
        else:
            return self.oldgetattr(name)

    def add(self, func):
        self.funcmap[func.__name__] = func


def moduleProperty(func):
    """Decorator to turn module functions into properties.
    Function names must be prefixed with an underscore."""
    module = sys.modules[func.__module__]
    def base_getattr(name):
        raise AttributeError(
            f"module '{module.__name__}' has no attribute '{name}'")
    ag = AttrGetter(getattr(module, '__getattr__', base_getattr))
    module.__getattr__ = ag
    ag.add(func)
    return func
