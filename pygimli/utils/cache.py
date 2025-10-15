#!/usr/bin/env python
"""Caching manager with function decorator.

Input supports python base types and all pg.core objects with .hash() method.
Output supports DataContainerERT, ...

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
import traceback
import inspect
import hashlib
import json
import time

import numpy as np
import pygimli as pg


__NO_CACHE__ = False

def noCache(c:bool=True):
    """ Set the caching to noCache mode.

    This will disable the caching mechanism and all decorated functions
    """
    global __NO_CACHE__
    __NO_CACHE__ = c


def strHash(s: str) -> int:
    """ Create a hash value for the given string.

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


def valHash(a:any)-> int:
    """ Create a hash value for the given value.

    Arguments
    ---------
    a: any
        The value to hash. Can be a string, int, list, numpy array or any
        other object. Logs an error if the type is not supported.
    Returns
    -------
    hash: int
        The hash value of the value.
    """
    if isinstance(a, str):
        return strHash(a)
    elif isinstance(a, int):
        return a
    elif isinstance(a, list):
        hsh = 0
        for i, item in enumerate(a):
            hsh = hsh ^ valHash(str(i)+str(item))
        return hsh
    elif isinstance(a, np.ndarray):
        if a.ndim == 1:
            return hash(pg.Vector(a))
        elif a.ndim == 2:
            # convert to RVector to use mem copy
            return hash(pg.Vector(a.reshape((1,a.shape[0]*a.shape[1]))[0]))
        else:
            print(a)
            pg.error('no hash for numpy array')

    return hash(a)


class Cache(object):
    """ Cache class to store and restore data.

    This class is used to store and restore data in a cache.
    """
    def __init__(self, hashValue:int):
        """ Initialize the cache with a hash value.

        Arguments
        ---------
        hashValue: int
            The hash value of the function and its arguments.
        """
        self._value = None
        self._hash = hashValue
        self._name = CacheManager().cachingPath(str(self._hash))
        self._info = None
        self.restore()


    @property
    def info(self):
        """ Return the cache info dictionary.

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
        """ Set the cache info dictionary.

        Arguments
        ---------
        i: dict
            The cache info dictionary to set.
        """
        self._info = i


    @property
    def value(self):
        """ Return the cached value.
        """
        return self._value


    @value.setter
    def value(self, v):
        """ Set the cached value and store it in the cache.

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

        self.updateCacheInfo()

        if self.info['type'] == 'Mesh':
            pg.info('Save Mesh binary v2')
            v.saveBinaryV2(self._name)
        elif self.info['type'] == 'RVector':
            pg.info('Save RVector binary')
            v.save(self._name, format=pg.core.Binary)
        elif self.info['type'] == 'ndarray':
            pg.info('Save ndarray')
            np.save(self._name, v, allow_pickle=True)
        elif hasattr(v, 'save') and hasattr(v, 'load'):
            v.save(self._name)
        else:
            np.save(self._name, v, allow_pickle=True)
            # pg.warn('ascii save of type', self.info['type'],
            # 'might by dangerous')
            # v.save(self._name)

        self._value = v
        pg.info('Cache stored:', self._name)


    def updateCacheInfo(self):
        """ Update the cache info dictionary and save it to a json file.
        """
        with open(self._name + '.json', 'w') as of:
            json.dump(self.info, of, sort_keys=False,
                      indent=4, separators=(',', ': '))

    def restore(self):
        """ Restore cache from json infos.
        """
        if os.path.exists(self._name + '.json'):

            # Fricking mpl kills locale setting to system default .. this went
            # horrible wrong for german 'decimal_point': ','
            pg.checkAndFixLocaleDecimal_point(verbose=False)

            try:
                with open(self._name + '.json') as file:
                    self.info = json.load(file)

                # if len(self.info['type']) != 1:
                #     pg.error('only single return caches supported for now.')

                #pg._y(pg.pf(self.info))

                if self.info['type'] == 'DataContainerERT':
                    self._value = pg.DataContainerERT(self.info['file'],
                                                      removeInvalid=False)
                    # print(self._value)
                elif self.info['type'] == 'RVector':
                    self._value = pg.Vector()
                    self._value.load(self.info['file'], format=pg.core.Binary)
                elif self.info['type'] == 'Mesh':
                    pg.tic()
                    self._value = pg.Mesh()
                    self._value.loadBinaryV2(self.info['file'] + '.bms')
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
                    self._value = np.load(self.info['file'] + '.npy',
                                          allow_pickle=True)

                if self.value is not None:
                    self.info['restored'] = self.info['restored'] + 1
                    self.updateCacheInfo()
                    pg.info('Cache {3} restored ({1}s x {0}): {2}'.\
                        format(self.info['restored'],
                               round(self.info['dur'], 1),
                               self._name, self.info['codeinfo']))
                else:
                    # default try numpy
                    pg.warn('Could not restore cache of type '
                            f'{self.info["type"]}.')

                pg.debug("Restoring cache took:", pg.dur(), "s")
            except BaseException as e:
                traceback.print_exc(file=sys.stdout)
                print(self.info)
                pg.error('Cache restoring failed:', e)


class CacheManager(object):
    """ Cache manager to handle caching of functions and data.

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
        """ Create a new instance of the CacheManager.
        """
        if cls.__instance is None:
            cls.__instance = object.__new__(cls)
        return cls.__instance


    def __init__(self):
        """ Initialize the CacheManager just once.
        """
        if not self.__has_init:
            self._caches = {}
            self.__has_init = True


    @staticmethod
    def instance(cls):
        """ Get the singleton instance of the CacheManager.
        """
        return cls.__instance__


    def cachingPath(self, fName:str):
        """ Create a full path name for the cache.

        Arguments
        ---------
        fName: str
            The name of the file to cache.

        Returns
        -------
        path: str
            The full path to the cache file.
        """
        if pg.rc["globalCache"]:
            path = pg.getCachePath()
        else:
            path = ".cache"
        if not os.path.exists(path):
            os.mkdir(path)
        return os.path.join(path, fName)


    def funcInfo(self, func):
        """ Return unique info string about the called function.

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
        """ Create a hash value.

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
        pg.tic()
        funcInfo = self.funcInfo(func)
        funcHash = strHash(funcInfo)
        versionHash = strHash(pg.versionStr())
        codeHash = strHash(inspect.getsource(func))

        argHash = 0
        for i, a in enumerate(args):
            if pg.isScalar(a):
                argHash = argHash ^ valHash(str(i) + str(a))
            else:
                argHash = argHash ^ (valHash(i) ^ valHash(a))

        for k, v in kwargs.items():
            if pg.isScalar(v):
                argHash = argHash ^ (valHash(k + str(v)))
            else:
                argHash = argHash ^ valHash(k) ^ valHash(v)

        pg.debug("Hashing took:", pg.dur(), "s")
        return funcHash ^ versionHash ^ codeHash ^ argHash


    def cache(self, func, *args, **kwargs):
        """ Create a unique cache.

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
        c.info['kwargs'] = str(kwargs)

        return c


def cache(func):
    """ Cache decorator.

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

        nc = kwargs.pop('skipCache', False)

        if any(('--noCache' in sys.argv,
                '-N' in sys.argv, nc is True, __NO_CACHE__)):

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
