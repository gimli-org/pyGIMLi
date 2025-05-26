#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Caching manager with function decorator.

Input supports python base types and all pg.core objects with .hash() method.
Output supports DataContainerERT, ...

TODO:

  *   Output types:
        numpy.ndarray, pg.Mesh. pg.Vector, pg.Matrix

To use just add the decorator.

@pg.cache
def myLongRunningStuff(*args, **kwargs):
    #...
    return results
"""
import sys
import os
import inspect
import hashlib
import json
import time

import numpy as np
import pygimli as pg
import pickle


__NO_CACHE__ = False

def noCache(c):
    """Set the cache to noCache mode."""
    global __NO_CACHE__
    __NO_CACHE__ = c


def strHash(string):
    """Create a hash value for the given string."""
    return int(hashlib.sha224(string.encode()).hexdigest()[:16], 16)


def valHash(a):
    """Create a hash value for the given value."""
    if isinstance(a, str):
        return strHash(a)
    elif isinstance(a, int):
        return a
    elif isinstance(a, dict):
        hsh = 0
        for k, v in a.items():
            # pg._y(k, valHash(k))
            # pg._y(v, valHash(v))
            hsh = hsh ^ valHash(k) ^ valHash(v)

        return hsh
    elif isinstance(a, list):
        hsh = 0
        for i, item in enumerate(a):
            hsh = hsh ^ valHash(str(i)+str(item))
        return hsh
    elif isinstance(a, pg.core.stdVectorNodes):
        ### cheap hash .. we assume the nodes are part of a mesh which is hashed anyways
        # pg._r('pg.core.stdVectorNodes: ', a, hash(len(a)))
        return hash(len(a))
    elif isinstance(a, np.ndarray):
        if a.ndim == 1:
            return hash(pg.Vector(a))
        elif a.ndim == 2:
            # convert to RVector to use memcopy
            return hash(pg.Vector(a.reshape((1,a.shape[0]*a.shape[1]))[0]))
        else:
            print(a)
            pg.error('no hash for numpy array')
    elif hasattr(a, '__hash__') and not callable(a):
        # pg._r('has hash: ', a, hash(a))
        return hash(a)
    elif isinstance(a, pg.DataContainer):
        return hash(a)
    elif callable(a):
        try:
            if hasattr(a, '_funct'):
                ## FEAFunctions or any other wrapper containing lambda as _funct
                # pg._g(inspect.getsource(a._funct))
                return strHash(inspect.getsource(a._funct))
            # for lambdas
            # pg._r('callable: ', inspect.getsource(a))
            else:
                return strHash(inspect.getsource(a))
        except:
            return valHash(str(a))

    pg.critical('cannot find hash for:', a)
    return hash(a)


class Cache(object):
    """Cache class to store and restore data.

    This class is used to store and restore data in a cache.
    It is used to store and restore data in a cache.
    """
    def __init__(self, hashValue):
        self._value = None
        self._hash = hashValue
        self._name = CacheManager().cachingPath(str(self._hash))
        self._info = None
        self.restore()

    @property
    def info(self):
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
        self._info = i

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, v):
        self.info['type'] = str(type(v).__name__)

        # if len(self.info['type']) != 1:
        #     pg.error('only single return caches supported for now.')
        #     return

        self.info['file'] = self._name

        # pg._r(self.info)

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
        elif isinstance(v, tuple) and 0:
            pass
            # pg._r(v)
            # with open(self._name + '.npy', 'wb') as f:
            #     for i in range(len(v)):
            #         pg._b(i, len(v), '#'*60)
            #         pg._g(v[i])
            #         np.save(f, v[i])

        else:
            self.info['type'] = 'pickle'
            with open(self._name + '.pkl', 'wb') as f:
                pickle.dump(v, f)
            #np.save(self._name, np.array(v, dtype=object), allow_pickle=True)
            # pg.warn('ascii save of type', self.info['type'], 'might by dangerous')
            # v.save(self._name)

        self.updateCacheInfo()
        self._value = v
        pg.info('Cache stored:', self._name)

    def updateCacheInfo(self):
        with open(self._name + '.json', 'w') as of:
            json.dump(self.info, of, sort_keys=False,
                      indent=4, separators=(',', ': '))

    def restore(self):
        """Read data from json infos"""

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
                elif self.info['type'] == 'tuple' and 0:
                    r = []
                    with open(self.info['file'] + '.npy', 'rb') as f:
                        v = np.load(f)
                        r.append(v)
                    self._value = (*r,)

                else:
                    with open(self.info['file'] + '.pkl', 'rb') as f:
                        self._value = pickle.load(f)

                if self.value is not None:
                    self.info['restored'] = self.info['restored'] + 1
                    self.updateCacheInfo()
                    pg.info('Cache {3} restored ({1}s x {0}): {2}'.\
                        format(self.info['restored'],
                               round(self.info['dur'], 1),
                               self._name, self.info['codeinfo']))
                else:
                    # default try numpy
                    pg.warn('Could not restore cache of type {0}.'.format(self.info['type']))

                pg.debug("Restoring cache took:", pg.dur(), "s")
            except Exception as e:
                import traceback
                traceback.print_exc(file=sys.stdout)
                print(self.info)
                pg.error('Cache restoring failed.')


#@pg.singleton
class CacheManager(object):
    __instance = None
    __has_init = False

    def __new__(cls):
        if cls.__instance is None:
            cls.__instance = object.__new__(cls)
        return cls.__instance

    def __init__(self):
        if not self.__has_init:
            self._caches = {}
            self.__has_init = True

    @staticmethod
    def instance(cls):
        return cls.__instance__

    def cachingPath(self, fName):
        """Create a path name for the cache"""
        if pg.rc["globalCache"]:
            path = pg.getCachePath()
        else:
            path = ".cache"
        if not os.path.exists(path):
            os.mkdir(path)
        return os.path.join(path, fName)

    def functInfo(self, funct):
        """Return unique info string about the called function."""
        return funct.__code__.co_filename + ":" + funct.__qualname__

    def hash(self, funct, *args, **kwargs):
        """"Create a hash value"""
        pg.tic()
        funcInfo = self.functInfo(funct)
        funcHash = strHash(funcInfo)
        versionHash = strHash(pg.versionStr())
        codeHash = strHash(inspect.getsource(funct))

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

    def cache(self, funct, *args, **kwargs):
        """ Create a unique cache """
        hashVal = self.hash(funct, *args, **kwargs)

        cached = Cache(hashVal)
        cached.info['codeinfo'] = self.functInfo(funct)
        cached.info['version'] = pg.versionStr()
        cached.info['args'] = str(args)
        #cached.info['kwargs'] = pg.pf(kwargs)

        kw = dict(kwargs)
        for k, v in kw.items():
            if isinstance(v, np.ndarray):
                kw[k] = 'ndarray'
        cached.info['kwargs'] = str(kw)

        return cached


def cache(funct):
    """Cache decorator."""
    def wrapper(*args, **kwargs):

        nc = kwargs.pop('skipCache', False)

        if any(('--noCache' in sys.argv, '--skipCache' in sys.argv,
                os.getenv('SKIP_CACHE'),
                '-N' in sys.argv, nc is True, __NO_CACHE__)):

            return funct(*args, **kwargs)

        cache = CacheManager().cache(funct, *args, **kwargs)

        if cache.value is not None:
            return cache.value
        else:
            # pg.tic will not work because there is only one global __swatch__
            sw = pg.Stopwatch(True)
            rv = funct(*args, **kwargs)
            cache.info['date'] = time.time()
            cache.info['dur'] = sw.duration()
            try:
                cache.value = rv
            except Exception as e:
                print(e)
                pg.warn("Can't cache:", rv)
            return rv
    return wrapper
