# -*- coding: utf-8 -*-
"""Utility functions for downloading, caching, and importing."""

import sys
import os.path

from importlib import import_module
import contextlib
from urllib.request import urlopen

import numpy as np
import pygimli as pg
from pygimli.meshtools import (readFenicsHDF5Mesh, readGmsh, readPLC, readSTL,
                               readMeshIO)
from pygimli.utils import readGPX
# from pygimli.utils import cache  # not used yet
from pygimli.physics.traveltime import load as loadTT


__gimliExampleDataRepo__ = 'gimli-org/example-data/'
__gimliExampleDataBase__ = 'example-data'


def optImport(module, requiredFor="use the full functionality"):
    """Import and return module only if it exists.

    If `module` cannot be imported, a warning is printed followed by the
    `requiredFor` string. Otherwise, the imported `module` will be returned.
    This function should be used to import optional dependencies in order to
    avoid repeated try/except statements.

    Parameters
    ----------
    module : str
        Name of the module to be imported.
    requiredFor : str, optional
        Info string for the purpose of the dependency.

    Examples
    --------
    >>> import pygimli as pg
    >>> pg = pg.optImport("pygimli")
    >>> pg.__name__
    'pygimli'
    >>> pg.optImport("doesNotExist", requiredFor="do something special")
    No module named 'doesNotExist'.
    You need to install this optional dependency to do something special.
    """
    # set default message for common imports
    if not requiredFor and "matplotlib" in module:
        requiredFor = "visualize 2D content"

    if module.count(".") > 2:
        raise ImportError("Can only import modules and sub-packages.")

    mod = None
    try:
        mod = import_module(module)
    except ImportError:
        msg = ("No module named \'%s\'.\nYou need to install this optional "
               "dependency to %s.")
        print(msg % (module, requiredFor))
        # exception would be better her but then the test fails
        # raise Exception(msg % (module, requiredFor))

    return mod


def opt_import(*args, **kwargs):
    """Optionally import module (deprecated)."""
    pg.deprecated()  # last vis: 20190903
    return optImport(*args, **kwargs)


def getCachePath():
    r"""Return the user space path for cache files.

    - Windows: 'C:\Documents and Settings\username\Appdata\pygimli\Cache'
    - Linux: '~/.cache/pygimli' (if not overwritten by $XDG_CONFIG_HOME)
    - Mac: '~/Library/Caches/pygimli'

    See Also
    --------
    pygimli.getConfigPath
    """
    system = sys.platform
    configpath = pg.core.config.getConfigPath()
    if system == "win32":
        path = os.path.join(configpath, "Cache")
    elif system == "darwin":
        path = configpath.replace("Preferences", "Caches")
    else:
        path = configpath.replace(".config", ".cache")
    return path


def load(fname, verbose=False, testAll=True, realName=None):
    """General import function to load data and meshes from file.

    Parameters
    ----------
    fname : string
        Filename or folder of files to load.
    testAll : bool [True]
        Test all filter when file suffix is unknown or loading fails.
    realName : str [None]
        Real file name.
        When fname is generic (i.e. without suffix) we can test the
        realName for type info.
    verbose : bool
        Be verbose.

    Examples
    --------
    >>> import os, tempfile
    >>> import pygimli as pg
    >>> fname = tempfile.mktemp(suffix=".bms")
    >>> pg.createGrid(range(3), range(3)).save(fname)
    1
    >>> mesh = pg.load(fname)
    >>> os.remove(fname)
    >>> mesh.cellCount()
    4
    """
    ImportFilter = {
        # maybe inflate the importer list from the submodules itself.
        # Data
        ".dat": pg.DataContainerERT,
        ".data": pg.DataContainerERT,
        ".ohm": pg.DataContainerERT,  # BERT compatibility
        ".shm": pg.DataContainerERT,  # BERT compatibility
        ".sgt": loadTT,
        ".gtt": loadTT,
        ".tom": loadTT,
        ".collect": pg.core.DataMap,
        # Vectors
        # ".dat": pg.Vector,
        ".vector": pg.Vector,
        ".vec": pg.Vector,
        ".idx": pg.IVector,
        # Matrices
        ".bmat": pg.Matrix,
        ".mat": pg.Matrix,
        ".matrix": pg.matrix.SparseMapMatrix,
        # Meshes
        ".poly": readPLC,
        ".bms": pg.Mesh,
        ".msh": readGmsh,
        ".mod": pg.Mesh,
        ".vtk": [pg.Mesh, readMeshIO],
        ".vtu": [pg.Mesh, readMeshIO],
        ".stl": readSTL,
        ".h5": readFenicsHDF5Mesh,  # fenics specs as default
        # Misc
        ".gpx": readGPX,  # read gpx waypoints
        ".xy": np.loadtxt,  #
        ".txt": np.loadtxt,  #
        ".npy": np.load,  # numpy binary format (important for pg.getExample*)
        ".npz": np.load,  #
    }

    if fname.startswith('https://') or fname.startswith('http://'):
        return getExampleData(fname)

    if not os.path.exists(fname):
        raise Exception("File or directory named %s does not exist." % (fname))

    # recursive function call if fname is a folder of files
    if os.path.isdir(fname):
        files = os.listdir(fname)
        if verbose:
            print("Reading %s with %d files..." % (fname, len(files)))

        return [load(os.path.join(fname, f)) for f in files]

    suffix = None
    if realName is not None:
        suffix = os.path.splitext(realName)[1]
    else:
        suffix = os.path.splitext(fname)[1]

    if suffix in ImportFilter:
        importTrys = ImportFilter[suffix]
        if not isinstance(importTrys, list):
            importTrys = [importTrys]

        for importer in importTrys:
            try:
                if verbose:
                    pg.info("Reading {0} ({1})".format(fname, importer))
                return importer(fname)
            except Exception as e:
                if verbose or pg.core.debug():
                    import sys
                    import traceback
                    traceback.print_exc(file=sys.stdout)
                    pg.warn(e)
        if verbose:
            pg.warn("File extension %s seems to be not correct. "
                    "Trying auto-detect." % suffix)
    else:
        if verbose:
            print("File extension {0} is unknown. Trying auto-detect.".format(
                suffix))

    if testAll:
        for routine in ImportFilter.values():
            try:
                return routine(fname)
            except Exception:
                # print(e)
                pass

    raise Exception("File type of {0} is unknown or file does not exist "
                    "and could not be imported.".format(suffix))


def getMD5(fileName):
    """Return md5 checksum for given fileName."""
    import hashlib
    md5 = hashlib.md5()

    with open(fileName, "rb") as fi:
        content = fi.read()
        md5.update(content)

    return md5.hexdigest()


def getUrlFile(url, fileName, timeout=10, verbose=False):
    """Write file from url. Path will be created."""
    import hashlib
    md5_hash = hashlib.md5()

    with contextlib.closing(urlopen(url, timeout=timeout)) as fp:
        header = dict(fp.getheaders())
        # print(pg.pf(header))
        length = int(header['Content-Length'])

        if verbose:
            p = pg.utils.ProgressBar(length)

        blockSize = 1024 * 8
        block = fp.read(blockSize)
        blockCounter = 1
        if block:
            os.makedirs(os.path.dirname(fileName), exist_ok=True)
            with open(fileName, 'wb') as outFile:
                if verbose:
                    p(min(length-1, blockCounter*blockSize))

                md5_hash.update(block)
                outFile.write(block)
                while True:
                    block = fp.read(blockSize)
                    blockCounter += 1
                    if not block:
                        break
                    if verbose:
                        p(min(length-1, blockCounter*blockSize))
                    md5_hash.update(block)
                    outFile.write(block)
        else:
            pg.critical('File or connection error:', url, fileName)
    if verbose:
        digest = md5_hash.hexdigest()
        print('md5:', digest)
        # print('md5:', getMD5(fileName))


def getExampleFile(path, load=False, force=False, verbose=False, **kwargs):
    """Download and return a filename to the example repository.

    Content will not be downloaded if already exists.
    TODO:
        * checksum or hash test for the content.
        * provide release or github version for specific repo

    Args
    ----
    path: str
        Path to the remote repo
    load: bool [False]
        Try to load the file and return the relating object.
    force: bool [False]
        Don't use cached file and get it in every case.

    Keyword Args
    ------------
    githubRepo: str
        Third party github data repository.
    branch: str ['master']
        branchname

    Returns
    -------
    filename: str
        Local file name to the data content.
    data: obj
        Content of the file if 'load' is set True.
    """
    repo = kwargs.pop('githubRepo', __gimliExampleDataRepo__)
    branch = kwargs.pop('branch', 'master')

    fileName = ''
    if not path.startswith('http://'):
        url = '/'.join(('https://raw.githubusercontent.com/',  # RAW files
                        repo, branch, path))

        pg.info(f'Looking for {path} in {repo}')

        fileName = os.path.join(getCachePath(),
                                __gimliExampleDataBase__,
                                repo, branch, path)
    else:
        url = path
        fileName = os.path.join(getCachePath(),
                                __gimliExampleDataBase__,
                                *url.split('http://')[1].split('/'),
                                )

        if verbose is True:
            pg.info(f'Looking for {url}')

    if not os.path.exists(fileName) or force is True:
        if verbose:
            pg.info(f'Getting: {fileName} from {url}')

        getUrlFile(url, fileName, verbose=verbose)
    else:
        if verbose:
            pg.info("File already exists:", fileName)

    if load is True:
        return pg.load(fileName, verbose=verbose)

    return fileName


def getExampleData(path, verbose=False):
    """Shortcut to load example data."""
    return getExampleFile(path, load=True, verbose=verbose)
