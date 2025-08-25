#!/usr/bin/env python3
"""Lazy import utility class."""

import importlib


class LazyImport:
    """Ensure lazy imports on first use only.

    Idea: https://stackoverflow.com/questions/4177735/\
        best-practice-for-lazy-loading-python-modules
    """

    def __init__(self, moduleName) :
        self._moduleName  = moduleName
        self._loaded  = None

    def __getattr__ (self, attr) :
        """Load on first access."""
        try:
            return getattr(self._loaded, attr)
        except AttributeError:
            self._loaded = importlib.import_module(self._moduleName)
            return getattr(self._loaded, attr)

