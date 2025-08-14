"""Utility functions for ERT data processing."""
import numpy as np
# from numpy import ma

import pygimli as pg
from .ert import createGeometricFactors


def uniqueERTIndex(data, nI=0, reverse=False, unify=True):
    """Generate unique index from sensor indices A/B/M/N for matching

    Parameters
    ----------
    data : DataContainerERT
        data container holding a b m n field registered as indices (int)
    nI : int [0]
        index to generate (multiply), by default (0) sensorCount
        if two data files with different sensorCount are compared make sure
        to use the same nI for both
    reverse : bool [False]
        exchange current (A, B) with potential (M, N) for reciprocal analysis
    unify : bool [True]
        sort A/B and M/N so that bipole orientation does not matter
    """
    if nI == 0:
        nI = data.sensorCount() + 1

    if unify:
        normABMN = {'a': np.minimum(data['a'], data['b']) + 1,
                    'b': np.maximum(data['a'], data['b']) + 1,
                    'm': np.minimum(data['m'], data['n']) + 1,
                    'n': np.maximum(data['m'], data['n']) + 1}
    else:
        normABMN = {tok: data[tok] + 1 for tok in "abmn"}

    abmn = "abmn"
    if reverse:
        abmn = "mnab"

    ind = np.zeros(data.size(), dtype=np.int64)
    for el in abmn:
        ind = ind * nI + normABMN[el]

    return ind  # np.array(ind, dtype=int)

def generateDataFromUniqueIndex(ind, data=None, nI=None):
    """Generate data container from unique index."""
    scheme = pg.DataContainerERT()
    if isinstance(data, pg.DataContainer):
        scheme = pg.DataContainerERT(data)
    elif isinstance(data, pg.PosVector):
        scheme.setSensorPositions(data)
    elif isinstance(data, int):  # check for positions
        for i in range(data):
            scheme.createSensor([i, 0, 0])

    nI = nI or scheme.sensorCount() + 1
    scheme.resize(0)  # make sure all data are deleted
    scheme.resize(len(ind))
    nmba = np.zeros([len(ind), 4], dtype=int)
    ind1 = ind.copy()
    for i in range(4):
        col = ind1 % nI
        ind1 -= col
        ind1 = ind1 // nI
        nmba[:, i] = col

    for i, tok in enumerate("nmba"):
        scheme[tok] = pg.IVector(nmba[:, i]) - 1

    scheme["valid"] = 1
    return scheme

def reciprocalIndices(data, onlyOnce=False, unify=True):
    """Return indices for reciprocal data.

    Parameters:
    data : DataContainerERT
        data containing reciprocal data
    onlyOnce : bool [False]
        return every pair only once

    Returns
    -------
    iN, iR : np.array(dtype=int)
        indices into the data container for normal and reciprocals
    """
    unF = uniqueERTIndex(data, unify=unify)
    unB = uniqueERTIndex(data, unify=unify, reverse=True)
    iF, iB = [], []
    for i, u in enumerate(unF):
        ii = np.nonzero(unB == u)[0]
        if len(ii) > 0:
            iF.append(i)
            iB.append(ii[0])

    iF = np.array(iF, dtype=int)
    iB = np.array(iB, dtype=int)
    if onlyOnce:
        return iF[iF < iB], iB[iF < iB]
    else:
        return iF, iB

def getResistance(data):
    """Return data resistance."""
    if data.allNonZero('r'):
        R = data['r']
    elif data.haveData('u') and data.haveData('i'):
        R = data['u'] / data['i']
    else:  # neither R or U/I
        if not data.allNonZero('k'):
            data['k'] = createGeometricFactors(data)

        if data.haveData('rhoa'):
            R = data['rhoa'] / data['k']
        else:
            pg.error('Either R, U/I or rhoa must be defined!')

    return R

def fitReciprocalErrorModel(data, nBins=None, show=False, rel=False):
    """Fit an error by statistical normal-reciprocal analysis.

    Identify normal reciprocal pairs and fit error model to it

    Parameters
    ----------
    data : pg.DataContainerERT
        input data
    nBins : int
        number of bins to subdivide data (4 < data.size()//30 < 30)
    rel : bool [False]
        fit relative instead of absolute errors
    show : bool [False]
        generate plot of reciprocals, mean/std bins and error model

    Returns
    -------
    ab : [float, float]
        indices of a+b*R (rel=False) or a+b/R (rec=True)
    """
    R = getResistance(data)
    iF, iB = reciprocalIndices(data, True)
    n30 = len(iF) // 30
    nBins = nBins or np.maximum(np.minimum(n30, 30), 4)
    RR = np.abs(R[iF] + R[iB]) / 2
    sInd = np.argsort(RR)
    RR = RR[sInd]
    dR = (R[iF] - R[iB])[sInd]
    inds = np.linspace(0, len(RR), nBins+1, dtype=int)
    stdR = np.zeros(nBins)
    meanR = np.zeros(nBins)
    for b in range(nBins):
        ii = range(inds[b], inds[b+1])
        stdR[b] = np.std(dR[ii])
        meanR[b] = np.mean(RR[ii])

    G = np.ones([len(meanR), 2]) # a*x+b
    w = np.reshape(np.isfinite(meanR), [-1, 1])
    meanR[np.isnan(meanR)] = 1e-16
    stdR[np.isnan(stdR)] = 0
    if rel:
        G[:, 1] = 1 / meanR
        ab, *_ = np.linalg.lstsq(w*G, stdR/meanR, rcond=None)
    else:
        G[:, 1] = meanR
        ab, *_ = np.linalg.lstsq(w*G, stdR, rcond=None)

    if show:
        x = np.logspace(np.log10(min(meanR)), np.log10(max(meanR)), 30)
        _, ax = pg.plt.subplots()
        if rel:
            ax.semilogx(RR, dR/RR, '.', label='data')  # /RR
            ii = meanR > 1e-12
            x = x[x > 1e-12]
            eModel = ab[0] + ab[1] / x
            ax.plot(meanR[ii], stdR[ii]/meanR[ii], '*', label='std dev') # /meanR
            ax.set_title(r'$\delta$R/R={:.6f}+{:.6f}/|R|'.format(*ab))
            ax.set_ylabel(r'$\delta$R/R')
        else:
            eModel = ab[0] + ab[1] * x
            ax.semilogx(RR, dR, '.', label='data')  # /RR
            ax.plot(meanR, stdR, '*', label='std dev') # /meanR
            ax.set_title(r'$\delta$R={:.6f}+{:.6f}*|R|'.format(*ab))
            ax.set_ylabel(r'$\delta$R ($\Omega$)')

        ax.plot(x, eModel, '-', label='error model') # /x
        ax.grid(which='both')
        ax.set_xlabel(r'R ($\Omega$)')
        ax.legend()
        return ab, ax
    else:
        return ab

def reciprocalProcessing(data, rel=True, maxrec=0.2, maxerr=0.2):
    """Reciprocal data analysis and error estimation.

    Carry out workflow for reciprocal
    * identify normal reciprocal pairs
    * fit error model to it
    * estimate error for all measurements
    * remove one of the pairs (and duplicates) by averaging
    * filter the data for maximum reciprocity and error

    Parameters
    ----------
    out : pg.DataContainerERT
        input data container with possible
    rel : bool
        base analysis on relative instead of absolute errors
    maxrec : float [0.2]
        maximum reciprocity allowed in data
    maxerr : float [0.2]
        maximum error estimate

    Returns
    -------
    out : pg.DataContainerERT
        filtered data container with pairs removed/averaged
    """
    out = data.copy()
    out.averageDuplicateData()
    iF, iB = reciprocalIndices(out, True)
    if len(iF) == 0:
        return data

    R = getResistance(out)
    meanR = np.abs(R[iF]+R[iB]) / 2
    dRbyR = np.abs(R[iF]-R[iB]) / meanR
    out['rec'] = 0
    out['rec'][iF] = dRbyR
    out.markInvalid(out['rec'] > maxrec)
    ab = fitReciprocalErrorModel(out, rel=rel)
    if rel:
        out["err"] = ab[0] + ab[1] / R
    else:
        out["err"] = ab[0] + ab[1] * R

    out['r'][iF] = meanR
    out.markInvalid(iB)
    out.markInvalid(out['err'] > maxerr)
    out.removeInvalid()
    return out

def getReciprocals(data, change=False, remove=False):
    """Compute data reciprocity from forward and backward data.

    The reciprocity (difference between forward and backward array divided by
    their mean) is computed and saved under the dataContainer field 'rec'

    Parameters
    ==========
    data : pg.DataContainerERT
        input data container to be changed inplace
    change : bool [True]
        compute current-weighted mean of forward and backward values
    remove : bool [False]
        remove backward data that are present as forward data
    """
    if not data.allNonZero('r'):
        data['r'] = data['u'] / data['i']

    unF = uniqueERTIndex(data)
    unB = uniqueERTIndex(data, reverse=True)
    rF, rB = [], []
    rec = np.zeros(data.size())
    data['rec'] = 0
    for iB in range(data.size()):
        if unB[iB] in unF:
            iF = int(np.nonzero(unF == unB[iB])[0][0])
            rF.append(data['r'][iF])
            rB.append(data['r'][iB])
            rec[iB] = (rF[-1]-rB[-1]) / (rF[-1]+rB[-1]) * 2
            data['rec'][iF] = rec[iB]
            IF, IB = data['i'][iF], data['i'][iB]  # use currents for weighting
            if change and data['valid'][iF]:
                data['r'][iF] = (rF[-1] * IF + rB[-1] * IB) / (IF + IB)
                data['i'][iF] = (IF**2 + IB**2) / (IF + IB)  # according weight
                data['u'][iF] = data['r'][iF] * data['i'][iF]
                if remove:
                    data['valid'][iB] = 0  # for adding all others later on

    print(len(rF), "reciprocals")
    if remove:
        data.removeInvalid()


def removeDuplicates(data, mode:str="average"):
    """Remove duplicate rows from an ERT datacontainer.
    
    Parameters
    ----------
    data : ERTDataContainer
        The data container from which duplicates should be removed.
    mode : str, optional
        The method to handle duplicates, choose between
        - "average"
        - "first"
        - "last"
        - "minerr"
        - "weighted"
    """
    ind = uniqueERTIndex(data)
    tokens = ["u", "i", "r", "rhoa", "err", "ip", "iperr"]
    for i in np.unique(ind):
        fi = np.nonzero(ind == i)[0]
        if len(fi) > 1:
            f0 = fi[0]
            data["valid"][f0] = True
            for ff in fi[1:]:
                data["valid"][ff] = False

            if mode == "average":
                for tok in tokens:
                    data[tok][f0] = np.mean(data[tok][fi])
            else:
                ff = f0
                if mode == "last":
                    ff = fi[-1]
                elif mode == "minerr":
                    ff = fi[np.argmin(data["err"][fi])]
                for tok in tokens:
                    data[tok][f0] = data[tok][ff]
    
    data.removeInvalid()  


def extractReciprocals(fwd, bwd):
    """Extract reciprocal data from forward/backward DataContainers."""
    nMax = max(fwd.sensorCount(), bwd.sensorCount())
    unF = uniqueERTIndex(fwd, nI=nMax)
    unB = uniqueERTIndex(bwd, nI=nMax, reverse=True)
    rF, rB = [], []
    rec = np.zeros(bwd.size())
    both = pg.DataContainerERT(fwd)
    both.set('rec', pg.Vector(both.size()))
    back = pg.DataContainerERT(bwd)
    back.set('rec', pg.Vector(back.size()))
    for iB in range(bwd.size()):
        if unB[iB] in unF:
            iF = int(np.nonzero(unF == unB[iB])[0][0])
            rF.append(fwd('r')[iF])
            rB.append(bwd('r')[iB])
            rec[iB] = (rF[-1]-rB[-1]) / (rF[-1]+rB[-1]) * 2
            both('rec')[iF] = rec[iB]
            IF, IB = fwd('i')[iF], bwd('i')[iB]  # use currents for weighting
            both('r')[iF] = (rF[-1] * IF + rB[-1] * IB) / (IF + IB)
            both('i')[iF] = (IF**2 + IB**2) / (IF + IB)  # according to weight
            both('u')[iF] = fwd('r')[iF] * fwd('i')[iF]
            back('valid')[iB] = 0  # for adding all others later on
    print(len(rF), "reciprocals")
    back.removeInvalid()
    both.add(back)
    return rec, both

def combineMultipleData(DATA):
    """Combine multiple data containers into data/err matrices.

    Generates a unified data container and a data matrices
    Non-existing data are filled with NaN values.
    
    Parameters
    ----------
    DATA : list of DataContainerERT
        data contea

    Returns
    -------
    data : pg.DataContainerERT
        ERT data container with quadrupols from all input data
    RHOA : np.array of size data.size() x len(DATA)
        apparent resistivity data matrix
    ERR : np.array (same size)
        relative error matrix
    IP : np.array(same size)
        induced polarization data matrix
    IPERR : np.array(same size)
        induced polarization error matrix
    """
    assert hasattr(DATA, '__iter__'), "DATA should be DataContainers or str!"
    if isinstance(DATA[0], str):  # read in if strings given
        DATA = [pg.DataContainerERT(data) for data in DATA]

    nEls = [data.sensorCount() for data in DATA]
    assert max(np.abs(np.diff(nEls))) == 0, "Electrodes not equal"
    uIs = [uniqueERTIndex(data) for data in DATA]
    uI = np.unique(np.hstack(uIs))
    scheme = generateDataFromUniqueIndex(uI, DATA[0])
    uI = uniqueERTIndex(scheme)   #, unify=False)
    R = np.ones([scheme.size(), len(DATA)]) * np.nan
    IP = R.copy()
    ERR = np.zeros_like(R)
    IPERR = R.copy()
    if not scheme.haveData('k'):  # just do that only once
        scheme['k'] = createGeometricFactors(scheme)  # check numerical

    for i, di in enumerate(DATA):
        ii = np.searchsorted(uI, uIs[i])
        if not di.haveData('r'):
            if di.allNonZero('u') and di.allNonZero('i'):
                di['r'] = di['u']/di['i']
            elif di.allNonZero('rhoa'):
                di['r'] = di['rhoa'] / scheme['k'][ii]

        R[ii, i] = di['r']
        ERR[ii, i] = di['err']
        IP[ii, i] = di['ip']
        IPERR[ii, i] = di['iperr']

    RHOA = np.abs(np.reshape(scheme['k'], [-1, 1]) * R)
    return scheme, RHOA, ERR, IP, IPERR
