"""Pygimli base functions.

Some needs to be sorted. Need to fit nameing conventions!
"""
import os.path
import time

import numpy as np
import pygimli as pg


def rms(v, axis=None):
    """Compute the root mean square."""
    # abs for complex values
    return np.sqrt(np.mean(np.abs(v)**2, axis))


def rrms(a, b, axis=None):
    """Compute the relative (regarding a) root mean square."""
    # abs for complex values
    return rms(np.abs(a-b)/np.abs(a), axis)


def nanrms(v, axis=None):
    """Compute the root mean square excluding nan values."""
    # abs for complex values
    return np.sqrt(np.nanmean(np.abs(v)**2, axis))


def rmsWithErr(a, b, err, errtol=1):
    """Compute (abs-)root-mean-square of values with error above threshold."""
    fi = pg.find(err < errtol)
    return rms(a[fi] - b[fi])


def chi2(a, b, err, trans=None):
    """Return chi square value."""
    if trans is None:
        trans = pg.trans.Trans()

    d = (trans(a) - trans(b)) / trans.error(a, err)
    return pg.math.dot(d, d) / len(d)


# fc_cleaning compatibilty to bert
rmswitherr = rmsWithErr


def rrmsWithErr(a, b, err, errtol=1):
    """Compute root mean square of values with error above a threshold."""
    fi = pg.find(err < errtol)
    return rms((a[fi]-b[fi])/a[fi])


def gmat2numpy(mat):
    """Convert pygimli matrix into numpy.array.

    TODO implement correct rval
    """
    nmat = np.zeros((len(mat), len(mat[0])))
    for i, row in enumerate(mat):
        nmat[i] = row
    return nmat


def numpy2gmat(nmat):
    """Convert numpy.array into pygimli RMatrix.

    TODO implement correct rval
    """
    gmat = pg.Matrix()
    for arr in nmat:
        gmat.push_back(arr)  # pg.asvector(arr))
    return gmat


def rndig(a, ndig=3):
    """Round float using a number of counting digits."""
    if np.abs(a) < 1e-4:
        return a
    else:
        return np.around(a, ndig - int(np.ceil(np.log10(np.abs(a) + 1e-4))))


def num2str(a, fmtstr='%g'):
    """List of strings (deprecated, for backward-compatibility)."""
    return [fmtstr % rndig(ai) for ai in a]


def inthist(a, vals, bins=None, islog=False):
    """Return point of integral (cumulative) histogram.

    E.g. inthist(a, [25, 50, 75]) provides quartiles and median of an array"""
    if bins is None:
        bins = int(np.min((np.round(len(a) / 20), 10)))

    if islog:
        hists, edges = np.histogram(np.log(a), bins=bins)
    else:
        hists, edges = np.histogram(a, bins=bins)

    cums = np.cumsum(np.hstack((0., hists))) / np.sum(hists) * 100.
    out = np.interp(vals, cums, edges)
    if islog:
        return np.exp(out)
    else:
        return out


def interperc(a, trimval=3.0, islog=False, bins=None):
    """Return symmetric interpercentiles for alpha-trim outliers.

    E.g. interperc(a, 3) returns range of inner 94% (3 to 97%)
    which is particularly useful for colorscales)."""
    return inthist(a, np.array([trimval, 100. - trimval]),
                   bins=bins, islog=islog)


def interpExtrap(x, xp, yp):
    """numpy.interp interpolation function extended by linear extrapolation."""
    y = np.interp(x, xp, yp)
    y = np.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    return np.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2]) /
                    (xp[-1]-xp[-2]), y)


def saveResult(fname, data, rrms=None, chi2=None, mode='w'):
    """Save rms/chi2 results into filename."""
    pg.warn("utils.saveResult .. in use? (debug)")
    with open(fname, mode) as f:
        np.savetxt(f, data)
        if rrms is not None:
            f.write('\nrrms:{}\n'.format(rrms))
        if chi2 is not None:
            f.write('\nchi2:{}\n'.format(chi2))


def createDateTimeString(now=None):
    """Return datetime as string (e.g. for saving results)."""
    if now is None:
        now = time.localtime()
    return str(now.tm_year) + str(now.tm_mon).zfill(2) + \
        str(now.tm_mday).zfill(2) + '-' + \
        str(now.tm_hour).zfill(2) + '.' + \
        str(now.tm_min).zfill(2)


def getSavePath(folder=None, subfolder='', now=None):
    """TODO."""
    if folder is None:
        path = createResultPath(subfolder, now=now)
    else:
        path = createPath([folder, subfolder])
    return path


def createResultPath(subfolder, now=None):
    """Create a result Folder."""
    result = createDateTimeString(now)
    return createPath(['.', result, subfolder])


def createPath(pathList):
    """Create the path structure specified by list.

    Parameters
    ----------
    pathList: str | list(str)
        Create Path with option subpaths
    """
    if hasattr(pathList, '__iter__') and not isinstance(pathList, str):
        path = os.path.join('', *pathList)
    else:
        path = os.path.join('', pathList)

    try:
        os.makedirs(path)
    except FileExistsError:
        print(f'Path {path} already exists. Skipping')
    except OSError as e:
        pg.error(f'Unable to create path "{path}".')
        raise e
    return path


@pg.renamed(createResultPath, '1.2')  # 20200515
def createResultFolder(subfolder, now=None):
    pass


@pg.renamed(createPath, '1.2')  # 20200515
def createfolders(foldername_list):
    pass


@pg.renamed(createPath, '1.2')  # 20200515
def createFolders(pathList):
    pass
