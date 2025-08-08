#!/usr/bin/env python
"""Caching manager with function decorator.

Input supports python base types and all `pg.core` objects that are hashable
with the `.hash()` method.
Output supports everything that's pickle-able.
Drop a ticket if some `pg` types need hash or pickle support.

To use just add the decorator.

```
@pg.cache
def myLongRunningStuff(*args, **kwargs):
    #...
    return results
```

To use the cache without the decorator, you can call it also like this:
`pg.cache(myLongRunningStuff)(*args, **kwargs)`
"""
import sys
import os
from pathlib import Path
import traceback
import inspect
import hashlib
import json
import time
import pickle

import numpy as np
import pygimli as pg


class CacheUseSettings:
    """Class to manage cache settings."""

    __NO_CACHE__ = False

    @classmethod
    def setNoCache(cls, c:bool=True):
        """Set the caching to noCache mode.

        This will disable the caching mechanism and all decorated functions.
        """
        cls.__NO_CACHE__ = c


    @classmethod
    def isNoCache(cls):
        """Check if caching is disabled."""
        return cls.__NO_CACHE__


    @classmethod
    def reset(cls):
        """Reset the caching settings to default."""
        cls.__NO_CACHE__ = False


    @classmethod
    def inUse(cls, sysArgs:list, manualSkip:bool=False):
        """Check if caching is in use.

        Or disabled by command line arguments, prior global setting using
        `pygimli.utils.cache.CacheUseSettings.setNoCache()`
        or environment variable `SKIP_CACHE`.

        Arguments
        ---------
        sysArgs: list
            The command line arguments to check for cache disabling flags are:
            '--noCache', '--skipCache', '-N'

        manualSkip: bool
            If True, the cache will be skipped regardless of command line
            arguments.
        """
        return not any(('--noCache' in sysArgs,
                        '--skipCache' in sysArgs,
                        '-N' in sysArgs,
                        os.getenv('SKIP_CACHE'),
                        manualSkip is True,
                        cls.__NO_CACHE__))


def strHash(s: str) -> int:
    """Create a hash value for the given string.

    Uses sha224 to create a 16 byte hash value.

    Arguments
    ---------
    s: str
        The string to hash.

    Returns
    -------
    hash: int
        The hash value of the string.
    """
    return int(hashlib.sha224(s.encode()).hexdigest()[:16], 16)


def valHash(a:any, verbose:bool=False)-> int:
    """Create a hash value for the given value.

    Arguments
    ---------
    a: any
        The value to hash. Can be a string, int, list, numpy array or any
        other object. Logs an error if the type is not supported.

    verbose: bool
        If True, additional logging will be performed.

    Returns
    -------
    hash: int
        The hash value of the value.
    """
    if verbose is True:
        pg._b(type(a))

    if isinstance(a, np.ndarray):
        if a.ndim == 1:
            return hash(pg.Vector(a))
        elif a.ndim == 2:
            # convert to RVector to use mem copy
            return hash(pg.Vector(a.reshape((1,a.shape[0]*a.shape[1]))[0]))
        else:
            print(a)
            pg.error('no hash for numpy array')

    if isinstance(a, list | tuple):
        hsh = 1
        for i, item in enumerate(a):
            if hasattr(item, '__hash__') and not isinstance(item, list |np.ndarray):
                if verbose is True:
                    pg._g(i, type(item))
                h = valHash(str(i)) ^ hash(item)
            else:
                h = valHash(str(i)+str(item), verbose=verbose)\
                  ^ valHash(item, verbose=verbose)
                if verbose is True:
                    pg._y(i, type(item), h)
            hsh = hsh ^ h
        if verbose is True:
            pg._y(hsh)

        return hsh
    elif isinstance(a, str):
        return strHash(a)
    elif isinstance(a, int):
        return a
    elif pg.isScalar(a):
        return valHash(str(a))
    elif isinstance(a, dict):
        hsh = 0
        for k, v in a.items():
            if verbose is True:
                pg._y(k, valHash(k))
                pg._y(v, valHash(v))
            hsh = hsh ^ valHash(k) ^ valHash(v)

        return hsh
    elif isinstance(a, pg.core.stdVectorNodes):
        ### cheap hash .. we assume the nodes are part of a mesh which
        ### is hashed anyways
        # pg._r('pg.core.stdVectorNodes: ', a, hash(len(a)))
        return hash(len(a))
    elif isinstance(a, np.ndarray):
        if a.ndim == 1:
            return hash(pg.Vector(a))
        elif a.ndim == 2:
            # convert to RVector to use mem copy
            return hash(pg.Vector(a.reshape((1,a.shape[0]*a.shape[1]))[0]))
        else:
            print(a)
            pg.error('no hash for numpy array')
    elif hasattr(a, '__hash__'):
        return hash(a)
    elif isinstance(a, pg.DataContainer):
        # not used
        return hash(a)
    elif callable(a):
        try:
            if hasattr(a, '_func'):
                ## FEAFunctions or any other wrapper containing lambda as _func
                # pg._g(inspect.getsource(a._func))
                return strHash(inspect.getsource(a._func))
            # for lambdas
            # pg._r('callable: ', inspect.getsource(a))
            else:
                return strHash(inspect.getsource(a))
        except BaseException:
            return valHash(str(a))

    pg.critical('cannot find hash for:', a)
    return hash(a)


class Cache:
    """Cache class to store and restore data.

    This class is used to store and restore data in a cache.
    """

    def __init__(self, hashValue:int):
        """Initialize the cache with a hash value.

        Arguments
        ---------
        hashValue: int
            The hash value of the function and its arguments.
        """
        self._value = None
        self._hash = hashValue
        self._name = str(CacheManager().cachingPath(str(self._hash)))
        self._info = None
        self.restore()


    @property
    def info(self):
        """Return the cache info dictionary.

        This dictionary contains information about the cache like type, file,
        date, duration, restored count, code info, version, args and kwargs.
        """
        if self._info is None:
            self._info = {'type': '',
                          'file': '',
                          'date': 0,
                          'dur': 0.0,
                          'restored': 0,
                          'codeinfo': '',
                          'version': '',
                          'args': '',
                          'kwargs': {},
                          }
        return self._info


    @info.setter
    def info(self, i):
        """Set the cache info dictionary.

        Arguments
        ---------
        i: dict
            The cache info dictionary to set.
        """
        self._info = i


    @property
    def value(self):
        """Return the cached value."""
        return self._value


    @value.setter
    def value(self, v):
        """Set the cached value and store it in the cache.

        Arguments
        ---------
        v: any
            The value to cache. Can be a DataContainerERT, Mesh, RVector,
            ndarray or any other object with either a save method or can be
            pickled.
        """
        self.info['type'] = str(type(v).__name__)

        # if len(self.info['type']) != 1:
        #     pg.error('only single return caches supported for now.')
        #     return

        self.info['file'] = self._name

        # pg._r(self.info)

        if self.info['type'] == 'Mesh':
            pg.info('Save Mesh binary')
            v.save(self._name)
        elif self.info['type'] == 'RVector':
            pg.info('Save RVector binary')
            v.save(self._name, format=pg.core.Binary)
        elif self.info['type'] == 'ndarray':
            pg.info('Save ndarray')
            np.save(self._name, v, allow_pickle=True)
        elif hasattr(v, 'save') and hasattr(v, 'load'):
            v.save(self._name)
        else:
            self.info['type'] = 'pickle'
            with Path(self._name + '.pkl').open('wb') as f:
                pickle.dump(v, f)

        self.updateCacheInfo()
        self._value = v
        pg.info('Cache stored:', self._name)


    def updateCacheInfo(self):
        """Update the cache info dictionary and save it to a json file."""
        with Path(self._name + '.json').open('w', encoding='utf-8') as of:
            json.dump(self.info, of, sort_keys=False,
                      indent=4, separators=(',', ': '))

    def restore(self):
        """Restore cache from json infos."""
        if Path(self._name + '.json').exists():

            # mpl kills locale setting to system default .. this went
            # horrible wrong for german 'decimal_point': ','
            pg.checkAndFixLocaleDecimal_point(verbose=False)

            try:
                with Path(self._name + '.json').open(encoding='utf-8') as file:
                    self.info = json.load(file)

                if self.info['type'] == 'DataContainerERT':
                    self._value = pg.DataContainerERT(self.info['file'],
                                                      removeInvalid=False)
                    # print(self._value)
                elif self.info['type'] == 'RVector':
                    self._value = pg.Vector()
                    self._value.load(self.info['file'], format=pg.core.Binary)
                elif self.info['type'] == 'Mesh':
                    self._value = pg.Mesh()
                    self._value.load(self.info['file'] + '.bms')
                    pg.debug("Restoring cache took:", pg.dur(), "s")
                elif self.info['type'] == 'ndarray':
                    self._value = np.load(self.info['file'] + '.npy',
                                          allow_pickle=True)
                elif self.info['type'] == 'Cm05Matrix':
                    self._value = pg.matrix.Cm05Matrix(self.info['file'])
                elif self.info['type'] == 'GeostatisticConstraintsMatrix':
                    self._value = pg.matrix.GeostatisticConstraintsMatrix(
                                                            self.info['file'])
                else:
                    with Path(self.info['file'] + '.pkl').open('rb') as f:
                        self._value = pickle.load(f)

                if self.value is not None:
                    self.info['restored'] = self.info['restored'] + 1
                    self.updateCacheInfo()
                    pg.info(f'Cache {self.info["codeinfo"]} restored '
                            f'({round(self.info["dur"], 1)}s x '
                            f'{self.info["restored"]}): {self._name}')
                else:
                    # default try numpy
                    pg.warn('Could not restore cache of type '
                            f'{self.info["type"]}.')

                pg.debug("Restoring cache took:", pg.dur(), "s")
            except BaseException as e:
                traceback.print_exc(file=sys.stdout)
                print(self.info)
                pg.error('Cache restoring failed:', e)


#TODO unify singleton handling
class CacheManager:
    """Cache manager to handle caching of functions and data.

    This class is a singleton and should be accessed via the instance method.
    It provides methods to create unique cache paths, hash functions and
    cache function calls.

    TODO
    ----
        * Unify singleton handling
    """

    __instance = None
    __has_init = False

    def __new__(cls):
        """Create a new instance of the CacheManager."""
        if cls.__instance is None:
            cls.__instance = object.__new__(cls)
        return cls.__instance


    def __init__(self):
        """Initialize the CacheManager just once."""
        if not self.__has_init:
            self._caches = {}
            self.__has_init = True


    @classmethod
    def instance(cls):
        """Get the singleton instance of the CacheManager."""
        return cls.__instance__


    def cachingPath(self, fName:str):
        """Create a full path name for the cache.

        Arguments
        ---------
        fName: str
            The name of the file to cache.

        Returns
        -------
        path: str
            The full path to the cache file.
        """
        path = pg.getCachePath() if pg.rc["globalCache"] else ".cache"

        if not Path(path).exists():
            Path(path).mkdir()

        return Path(path) / fName


    def funcInfo(self, func):
        """Return unique info string about the called function.

        Arguments
        ---------
        func: function
            The function to get the info from.

        Returns
        -------
        info: str
            A string containing the file name and the qualified name of the
            function.
        """
        return func.__code__.co_filename + ":" + func.__qualname__


    def hash(self, func, *args, **kwargs):
        """Create a hash value.

        Arguments
        ---------
        func: function
            The function to hash.
        *args: any
            The positional arguments of the function.
        **kwargs: any
            The keyword arguments of the function.

        Returns
        -------
        hash: int
            A unique hash value for the function and its arguments.
        """
        with pg.tictoc('hash'):
            funcInfo = self.funcInfo(func)
            funcHash = strHash(funcInfo)
            versionHash = strHash(pg.versionStr())
            codeHash = strHash(inspect.getsource(func))

            #pg._b('fun:', funcHash, 'ver:', versionHash, 'code:', codeHash)
            argHash = valHash(args, verbose=False)
            kwargHash = valHash(kwargs, verbose=False)
            #pg._b('argHash:', argHash, kwargHash)
        return funcHash ^ versionHash ^ codeHash ^ argHash ^ kwargHash


    def cache(self, func, *args, **kwargs):
        """Create a unique cache.

        Arguments
        ---------
        func: function
            The function to cache.
        *args: any
            The positional arguments of the function.
        **kwargs: any
            The keyword arguments of the function.

        Returns
        -------
        c: Cache
            A Cache object containing the cached value, info and hash value.
        """
        hashVal = self.hash(func, *args, **kwargs)

        c = Cache(hashVal)
        c.info['codeinfo'] = self.funcInfo(func)
        c.info['version'] = pg.versionStr()
        c.info['args'] = str(args)

        kw = dict(kwargs)
        for k, v in kw.items():
            if isinstance(v, np.ndarray):
                kw[k] = 'ndarray'
        c.info['kwargs'] = str(kw)

        return c


def cache(func):
    """Cache decorator.

    This decorator caches the return value of the function and stores it in a
    Cache object. If the function is called again with the same arguments,
    the cached value is returned instead of calling the function again.
    If the cache is not found, the function is called and the result is stored
    in the cache.

    This can be used without using the decorator by calling:
    `pg.cache(func)(*args, **kwargs)`

    Arguments
    ---------
    func: function
        The function to cache.

    Returns
    -------
    wrapper: function
        A wrapper function that caches the return value of the function.
    """
    def wrapper(*args, **kwargs):

        if not CacheUseSettings.inUse(sys.argv,
                                    manualSkip=kwargs.pop('skipCache', False)):

            return func(*args, **kwargs)

        c = CacheManager().cache(func, *args, **kwargs)
        if c.value is not None:
            return c.value

        # pg.tic will not work because there is only one global __swatch__
        sw = pg.Stopwatch(True)
        rv = func(*args, **kwargs)
        c.info['date'] = time.time()
        c.info['dur'] = sw.duration()
        try:
            c.value = rv
        except Exception as e:
            print(e)
            pg.warn("Can't cache:", rv)
        return rv

    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper
