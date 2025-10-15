#!/usr/bin/env python3
from setuptools import setup

setup(name='pybtex-pgstyle',
      version='0.1.0',
      entry_points={'pybtex.style.formatting': ['pgstyle = pybtexPGstylePlugin:PGStyle'],
                    'pybtex.style.labels': ['alpha = pybtexPGstylePlugin:Alpha']},
      py_modules=['pybtexPGstylePlugin'])
