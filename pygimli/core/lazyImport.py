#!/usr/bin/env python3
"""Lazy import utility class."""

import importlib


class LazyImport:
    """Ensure lazy imports on first use only.

    Idea: https://stackoverflow.com/questions/4177735/\
        best-practice-for-lazy-loading-python-modules
    """

    def __init__(self, module) :
        self._module  = module
        self._loaded  = None

    def __getattr__ (self, attr) :
        """Load on first access."""
        try:
            return getattr(self._loaded, attr)
        except Exception as e:
            print(e)
            if self._loaded is None:
                # module not loaded -> load it
                self._loaded = importlib.import_module(self._module)
                return getattr(self._loaded, attr)
            else:
                # module is loaded but got any problems
                raise e
