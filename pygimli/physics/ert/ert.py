# -*- coding: utf-8 -*-
"""Electrical resistivity tomography"""

import numpy as np

import pygimli as pg

from .ertModelling import ERTModelling
from .ertScheme import createData
createERTData = createData  # backward compatibility


def simulate(mesh, scheme, res, **kwargs):
    """Simulate an ERT measurement.

    Perform the forward task for a given mesh, resistivity distribution &
    measuring scheme and return data (apparent resistivity) or potentials.

    For complex resistivity, the data contains an apparent phase or the
    returned potentials are complex.

    The forward operator itself only calculates potential values for the
    electrodes in the given data scheme.
    To calculate apparent resistivities, geometric factors (k) are needed.
    If there are no values k in the DataContainerERT scheme, the function
    tries to calculate them, either analytically or numerically by using a
    p2-refined version of the given mesh.

    TODO
    ----
    * 2D + Complex + SR

    Args
    ----
    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D Mesh to calculate for.

    res : float, array(mesh.cellCount()), array(N, mesh.cellCount()), list
        Resistivity distribution for the given mesh cells can be:
        . float for homogeneous resistivity (e.g. 1.0)
        . single array of length mesh.cellCount()
        . matrix of N resistivity distributions of length mesh.cellCount()
        . resistivity map as [[regionMarker0, res0],
                              [regionMarker0, res1], ...]

    scheme : :gimliapi:`GIMLI::DataContainerERT`
        Data measurement scheme.

    Keyword Args
    ------------
    verbose: bool[False]
        Be verbose. Will override class settings.
    calcOnly: bool [False]
        Use fop.calculate instead of fop.response. Useful if you want
        to force the calculation of impedances for homogeneous models.
        No noise handling. Solution is put as token 'u' in the returned
        DataContainerERT.
    noiseLevel: float [0.0]
        add normally distributed noise based on
        scheme['err'] or on noiseLevel if error>0 is not contained
    noiseAbs: float [0.0]
        Absolute voltage error in V
    returnArray: bool [False]
        Returns an array of apparent resistivities instead of
        a DataContainerERT
    returnFields: bool [False]
        Returns a matrix of all potential values (per mesh nodes)
        for each injection electrodes.
    sr : bool
        use secondary field (singularity removal)
    seed : int
        numpy.random seed for repeatable noise in synthetic experiments
    phiErr : float|iterable
        absolute phase error, if not given, data['iperr'] or noiseLevel is used
    contactImpedances float|iterables
        contact impedances for being used with CEM model
    current : float
        current to be assumed

    Returns
    -------
    DataContainerERT | array(data.size()) | array(N, data.size()) |
    array(N, mesh.nodeCount()):
        Data container with resulting apparent resistivity data['rhoa'] and
        errors (if noiseLevel or noiseAbs is set).
        Optionally return a Matrix of rhoa values
        (for returnArray==True forces noiseLevel=0).
        In case of complex-valued resistivity, phase values are contained in
        data['phia'] or returned as additionally returned array.

    Examples
    --------
    >>> from pygimli.physics import ert
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> world = mt.createWorld(start=[-50, 0], end=[50, -50],
    ...                        layers=[-1, -5], worldMarker=True)
    >>> scheme = ert.createData(
    ...                     elecs=pg.utils.grange(start=-10, end=10, n=21),
    ...                     schemeName='dd')
    >>> for pos in scheme.sensorPositions():
    ...     _= world.createNode(pos)
    ...     _= world.createNode(pos + [0.0, -0.1])
    >>> mesh = mt.createMesh(world, quality=34)
    >>> rhomap = [
    ...    [1, 100. + 0j],
    ...    [2, 50. + 0j],
    ...    [3, 10.+ 1j],
    ... ]
    >>> data = ert.simulate(mesh, res=rhomap, scheme=scheme, verbose=True)
    """
    verbose = kwargs.pop('verbose', True)
    calcOnly = kwargs.pop('calcOnly', False)
    returnFields = kwargs.pop("returnFields", False)
    returnArray = kwargs.pop('returnArray', False)
    noiseLevel = kwargs.pop('noiseLevel', 0.0)
    noiseAbs = kwargs.pop('noiseAbs', 1e-4)
    seed = kwargs.pop('seed', None)
    sr = kwargs.pop('sr', True)
    returnFOP = kwargs.pop("returnFOP", False)

    fop = ERTModelling(sr=sr, verbose=verbose)
    # fop = self.createForwardOperator(useBert=True, sr=sr, verbose=verbose)
    fop.data = scheme
    fop.setMesh(mesh, ignoreRegionManager=True)
    cI = kwargs.pop("contactImpedances", None)
    if cI is not None:
        if isinstance(cI, float):
            cI = pg.Vector(scheme.sensorCount(), cI)

        fop._core.setContactImpedances([1e-3, 1e-4, 1e-5, 1e-6])

    rhoa = None
    phia = None

    isArrayData = False
    # parse the given res into mesh-cell-sized array
    if isinstance(res, (int, float)):
        res = np.ones(mesh.cellCount()) * float(res)
    elif isinstance(res, complex):
        res = np.ones(mesh.cellCount()) * res
    elif hasattr(res[0], '__iter__'):  # ndim == 2
        if len(res[0]) == 2:  # res seems to be a res map
            # check if there are markers in the mesh that are not defined
            # the rhomap. better signal here before it results in errors
            meshMarkers = list(set(mesh.cellMarkers()))
            mapMarkers = [m[0] for m in res]
            if any([mark not in mapMarkers for mark in meshMarkers]):
                left = [m for m in meshMarkers if m not in mapMarkers]
                pg.critical("Mesh contains markers without assigned "
                            "resistivities {}. Please fix given "
                            "rhomap.".format(left))
            res = pg.solver.parseArgToArray(res, mesh.cellCount(), mesh)
        else:  # probably nData x nCells array
            # better check for array data here
            isArrayData = True

    if isinstance(res[0], complex) or isinstance(res, pg.CVector):
        pg.info("Complex resistivity values found.")
        fop.setComplex(True)
    else:
        fop.setComplex(False)

    if not isinstance(scheme, pg.DataContainerERT):
        raise TypeError("Scheme must be DataContainerERT!")

    if not scheme.allNonZero('k') and not calcOnly:
        if verbose:
            pg.info('Calculate geometric factors.')
        scheme.set('k', fop.calcGeometricFactor(scheme))

    ret = pg.DataContainerERT(scheme)
    # just to be sure that we don't work with artifacts
    ret['u'] = 0.0
    ret['i'] = 0.0
    ret['r'] = 0.0

    if isArrayData:
        rhoa = np.zeros((len(res), scheme.size()))
        for i, r in enumerate(res):
            rhoa[i] = fop.response(r)
            if verbose:
                print(i, "/", len(res), " : ", pg.dur(), "s",
                      "min r:", min(r), "max r:", max(r),
                      "min r_a:", min(rhoa[i]), "max r_a:", max(rhoa[i]))
    else:  # res is single resistivity array
        if len(res) == mesh.cellCount():

            if calcOnly:
                fop.mapERTModel(res, 0)

                dMap = pg.core.DataMap()
                fop.calculate(dMap)
                if fop.complex():
                    pg.critical('Implement me')
                else:
                    ret['r'] = dMap.data(scheme)
                    ret['i'] = kwargs.pop("current", 1.0)
                    ret['u'] = ret['r'] * ret['i']

                if returnFields:
                    return pg.Matrix(fop.solution())

                return ret
            else:
                if fop.complex():
                    res = pg.utils.squeezeComplex(res)

                resp = fop.response(res)

                if fop.complex():
                    rhoa, phia = pg.utils.toPolar(resp)
                else:
                    rhoa = resp
        else:
            print(mesh)
            print("res: ", res)
            raise BaseException(
                "Simulate called with wrong resistivity array.")

    if not isArrayData:
        ret['rhoa'] = rhoa

        if phia is not None:
            ret['phia'] = phia
    else:
        ret['rhoa'] = rhoa[0]
        if phia is not None:
            ret['phia'] = phia[0]

    if returnFields:
        return pg.Matrix(fop.solution())

    if noiseLevel > 0:  # if errors in data noiseLevel=1 just triggers
        if not ret.allNonZero('err'):
            # 1A  and #100µV
            ret['err'] = estimateError(ret,
                                         relativeError=noiseLevel,
                                         absoluteUError=noiseAbs,
                                         absoluteCurrent=1)
            if verbose:
                pg.info("Data error estimate (min:max) ",
                      min(ret['err']), ":", max(ret['err']))

        rhoa *= 1. + pg.randn(ret.size(), seed=seed) * ret('err')
        ret['rhoa'] = rhoa

        ipError = None
        if phia is not None:
            if scheme.allNonZero('iperr'):
                ipError = scheme['iperr']
            else:
                # np.abs(self.data("phia") +TOLERANCE) * 1e-4absoluteError
                if noiseLevel > 0.5:
                    noiseLevel /= 100.

                if 'phiErr' in kwargs:
                    ipError = np.ones(ret.size()) * kwargs.pop('phiErr') / 1000
                else:
                    ipError = abs(ret['phia']) * noiseLevel

                if verbose:
                    print("Data IP abs error estimate (min:max) ",
                          min(ipError), ":", max(ipError))

            phia += pg.randn(ret.size(), seed=seed) * ipError
            ret['iperr'] = ipError
            ret['phia'] = phia

    # check what needs to be setup and returned

    if ret.allNonZero('k'):  # also provide r if user changes k later
        ret['r'] = ret['rhoa'] / ret['k']
        ret['i'] = kwargs.pop("current", 1.0)
        ret['u'] = ret['r'] * ret['i']

    if returnArray:
        if phia is not None:
            return rhoa, phia
        else:
            return rhoa

    ret['valid'] = 1
    if returnFOP:
        return ret, fop
    else:
        return ret


def simulateOld(mesh, scheme, res, sr=True, useBert=True,
                verbose=False, **kwargs):
    """ERT forward calculation.

    Convenience function to use the ERT modelling operator
    if you like static functions.

    See :py:mod:`pygimli.ert.ERTManager.simulate` for description
    of the arguments.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh` | str
        Modelling domain. Mesh can be a file name here.
    scheme: :gimliapi:`GIMLI::DataContainerERT` | str
        Data configuration. Scheme can be a file name here.
    res: see :py:mod:`pygimli.ert.ERTManager.simulate`
        Resistivity distribution.
    sr: bool [True]
        Use singularity removal technique.
    useBert: bool [True]
        Use Bert forward operator instead of the reference implementation.
    **kwargs:
        Forwarded to :py:mod:`pygimli.ert.ERTManager.simulate`
    """
    from .ertManager import ERTManager
    ert = ERTManager(useBert=useBert, sr=sr, verbose=verbose)

    if isinstance(mesh, str):
        mesh = pg.load(mesh)

    if isinstance(scheme, str):
        scheme = pg.physics.ert.load(scheme)

    return ert.simulate(mesh=mesh, res=res, scheme=scheme,
                        verbose=verbose, **kwargs)


@pg.cache
def createGeometricFactors(scheme, numerical=None, mesh=None, dim=3,
                           h2=True, p2=True, forceFlatEarth=False, verbose=False):
    """Create geometric factors for a given data scheme.

    Create geometric factors for a data scheme with and without topography.
    Calculation will be done analytical (only for half space geometry)
    or numerical.

    This function caches the result depending on scheme, mesh and pg.version()

    Parameters
    ----------
    scheme: :gimliapi:`GIMLI::DataContainerERT`
        Datacontainer of the scheme.
    numerical: bool | None [False]
        If numerical is None, False is assumed, we try to guess topography
        and warn if we think we found them.
        If set to True or False, numerical calculation will used or not.
    mesh: :gimliapi:`GIMLI::Mesh` | str
        Mesh for numerical calculation. If not given, analytical geometric
        factors for halfspace earth are guessed or a default mesh will be
        created (and h/p refined according to h2/p2). If given topo is set to
        True. If the numerical effort is to high or the accuracy to low
        you should consider calculating the factors manually.
    h2: bool [True]
        Default spatial refinement to achieve high accuracy
    p2: bool [True]
        Default polynomial refinement to achieve high accuracy
    forceFlatEarth : bool [False]
        Set z values of electrodes to zero for geometric factor (usful only for topography).
    verbose: bool
        Give some output.
    """
    if numerical is None:
        numerical = False
        # if min(pg.z(scheme)) != max(pg.z(scheme)):
        #     verbose = True
        #     pg.warn('Sensor z-coordinates not equal. Is there topography?')

    if numerical is False and mesh is None:
        if verbose:
            pg.info('Calculate analytical flat-earth geometric factors.')

        # better to be replaced by a Python function analyticalGeometricFactors
        return pg.core.geometricFactors(scheme, forceFlatEarth=forceFlatEarth, dim=dim)

    if mesh is None:
        pg.info('Create default mesh for geometric factor calculation.')
        m = createInversionMesh(scheme)
    else:
        m = mesh

    if verbose:
        pg.info('mesh', m)

    if h2 is True:
        m = m.createH2()
        if verbose:
            pg.info('h2 refine', m)

    if p2 is True:
        m = m.createP2()
        if verbose:
            pg.info('p2 refine', m)

    if verbose:
        pg.info('Calculate numerical geometric factors.')

    d = simulate(m, res=1.0, scheme=scheme, sr=False, useBert=True,
                 calcOnly=True, verbose=True)
    return 1./d['u']


def createInversionMesh(data, **kwargs):
    """Create default mesh for ERT inversion.

    Parameters
    ----------
    data: :gimliapi:`GIMLI::DataContainerERT`
        Data Container needs at least sensors to define the geometry of the
        mesh.

    Other Parameters
    ----------------
        Forwarded to :py:mod:`pygimli.meshtools.createParaMesh`

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
        Inversion mesh with default marker (1 for background,
        2 parametric domain)
    """
    mesh = pg.meshtools.createParaMesh(data.sensors(), **kwargs)
    return mesh


def createERTDataNotUsedAnymore(elecs, schemeName='none', **kwargs):
    """Create data scheme for compatibility (advanced version in BERT).

    Parameters
    ----------
    sounding : bool [False]
        Create a 1D VES Schlumberger configuration.
        elecs need to be an array with elecs[0] = mn/2 and elecs[1:] = ab/2.
    """
    if kwargs.pop('sounding', False):
        data = pg.DataContainerERT()
        data.setSensors(pg.cat(-elecs[::-1], elecs))

        nElecs = len(elecs)
        for i in range(nElecs-1):
            data.createFourPointData(i, i, 2*nElecs-i-1, nElecs-1, nElecs)

        return data

    if schemeName != "dd":
        return createData(elecs, schemeName, **kwargs)

    isClosed = kwargs.pop('closed', False)

    data = pg.DataContainerERT()
    data.setSensors(elecs)

    nElecs = len(elecs)
    a = []
    b = []
    m = []
    n = []
    eb = 0
    for i in range(nElecs):
        for j in range(eb + 2, nElecs):
            ea = i
            eb = ea + 1
            em = j
            en = em + 1

            if isClosed:
                en = en % nElecs

            if en < nElecs and en != ea:
                a.append(ea)
                b.append(eb)
                m.append(em)
                n.append(en)

    data.resize(len(a))
    data.add('a', a)
    data.add('b', b)
    data.add('m', m)
    data.add('n', n)
    data.set('valid', np.ones(len(a)))

    return data


def estimateError(data, relativeError=0.03, absoluteUError=None,
                  absoluteError=0, absoluteCurrent=0.1):
    """Estimate relative error based on relative and absolute parts.

    Parameters
    ----------
    relativeError : float [0.03]
        relative error level in %/100 for u, R or rhoa

    absoluteUError : float [None]
        Absolute potential error in V. Needs 'u' values in data. Otherwise
        calculate them from 'rhoa', 'k' and absoluteCurrent if no 'i' given.

    absoluteError : float [0.0]
        Absolute data error in Ohm m. Needs 'R' or 'rhoa' values.

    absoluteCurrent : float [0.1]
        Current level in A for reconstruction for absolute potential V

    Returns
    -------
    error : Array
    """
    if relativeError >= 0.5:
        print("relativeError set to a value > 0.5 .. assuming this "
              "is a percentage Error level dividing them by 100")
        relativeError /= 100.0

    if absoluteUError is None:
        if not data.allNonZero('rhoa'):
            pg.critical("We need apparent resistivity values "
                        "(rhoa) in the data to estimate a "
                        "data error.")
        error = relativeError + pg.abs(absoluteError / data['rhoa'])
    else:
        u = None
        i = absoluteCurrent
        if data.haveData('i'):
            i = data['i']

        if data.haveData('u'):
            u = data['u']
        else:
            if data.haveData('r'):
                u = data['r'] * i
            elif data.haveData('rhoa'):
                if data.haveData('k'):
                    u = data['rhoa'] / data['k'] * i
                else:
                    pg.critical("We need (rhoa) and (k) in the"
                                "data to estimate data error.")

            else:
                pg.critical("We need apparent resistivity values "
                            "(rhoa) or impedances (r) "
                            "in the data to estimate data error.")

        error = pg.abs(absoluteUError / u) + relativeError

    return error


def __DataContainerERT_createGeometricFactors(self, *args,**kwargs):
    self['k'] = createGeometricFactors(self, *args, **kwargs)

pg.DataContainerERT.createGeometricFactors = __DataContainerERT_createGeometricFactors
pg.DataContainerERT.createGeometricFactors.__doc__ = createGeometricFactors.__doc__


def __DataContainerERT_estimateError(self, *args,**kwargs):
    if not self.haveData('k'):
        self.createGeometricFactors()

    self['err'] = estimateError(self, *args, **kwargs)

pg.DataContainerERT.estimateError = __DataContainerERT_estimateError
pg.DataContainerERT.estimateError.__doc__ = estimateError.__doc__.replace(
    "Estimate ", "Set ")[:-4]


if __name__ == "__main__":
    pass
