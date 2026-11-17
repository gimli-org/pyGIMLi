#!/usr/bin/env python3
"""Distribution setup of the pybtex style plugin for the pyGIMLi docs."""

from setuptools import setup

setup(
    name="pybtex-pgstyle",
    version="0.1.2",
    entry_points={
        "pybtex.style.formatting": ["pgstyle = pybtexPGstylePlugin:PGStyle"],
        "pybtex.style.labels": ["alpha = pybtexPGstylePlugin:Alpha"],
    },
    py_modules=["pybtexPGstylePlugin"],
)
