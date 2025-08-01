#!/usr/bin/env python
"""Useful utility functions."""

from .base import (rms, rmsWithErr, nanrms, createDateTimeString,
                   createResultPath, createPath, createPathTimePath,
                   getSavePath, gmat2numpy, interperc, interpExtrap, inthist,
                   nanrms, num2str, numpy2gmat, rrms, rms, rndig, saveResult,
                   chi2)

mkdir = createPath

# compatibility for dev #can be removed? (20200515)
from .base import rmsWithErr as rmswitherr

from .complex import (isComplex, toComplex, toPolar, squeezeComplex,
                      toRealMatrix, KramersKronig)

from .cache import (cache, strHash, CacheUseSettings, valHash)
from .geostatistics import (computeInverseRootMatrix, covarianceMatrix,
                            generateGeostatisticalModel)
from .gps import GKtoUTM, findUTMZone, getProjection, getUTMProjection, readGPX
from .hankel import hankelFC

from .postinversion import iterateBounds, modelCovariance, modelResolutionMatrix

from .units import (unit, cmap)
from .units import cmap as cMap # for compatibility (will be removed)

from .utils import (ProgressBar, Table, boxprint, cumDist, cut, diff, dist, rate,
                    filterIndex, filterLinesByCommentStr, findNearest,
                    getIndex, grange, logDropTol, niceLogspace,
                    prettify, prettyFloat, prettyTime,
                    randn, rand, trimDocString, unicodeToAscii, unique,
                    unique_everseen, uniqueAndSum)

from .streams import streamline, streamlineDir
from .dem import DEM


__all__ = [name for name in dir() if '_' not in name]
