# -*- coding: utf-8 -*-
"""
Basic import of gimli core extension.
"""
import sys, os

if sys.platform == 'win32':
    os.environ['PATH'] = os.path.join(os.path.dirname(__file__), 'libs') + ':' + os.environ['PATH']
    #print(os.environ['LD_LIBRARY_PATH'])
elif sys.platform == 'linux':
    if os.getenv('LD_LIBRARY_PATH'):
        os.environ['LD_LIBRARY_PATH'] = os.path.join(os.path.dirname(__file__), 'libs') + ':' + os.environ['LD_LIBRARY_PATH']
    #print(os.environ['LD_LIBRARY_PATH'])

_pygimli_ = None

try:
    # from . import _pygimli_
    # from ._pygimli_ import *

    from .libs import _pygimli_ as pgcore  
    from .libs._pygimli_ import *  

except ImportError as e:
    import traceback
    print(e)
    traceback.print_exc(file=sys.stdout)
    sys.stderr.write("ERROR: cannot import the library '_pygimli_'.\n")

