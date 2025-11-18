# -*- coding: utf-8 -*-

import numpy as np
# from numpy import ma

import pygimli as pg
from pygimli.physics import ert
from pygimli.viewer.mpl.colorbar import createColorBarOnly


def createData(elecs, schemeName='none', **kwargs):
    """ Utility one-liner to create a BERT datafile

    Parameters
    ----------
    elecs : int | list[pos] | array(x)
        Number of electrodes or electrode positions or x-positions

    schemeName : str ['none']
        Name of the configuration. If you provide an unknown scheme name, all
        known schemes ['wa', 'wb', 'pp', 'pd', 'dd', 'slm', 'hw', 'gr'] listed.

    **kwargs :

        Arguments that will be forwarded to the scheme generator.

        * inverse : bool
            interchange AB MN with MN AB
        * reciprocity : bool
            interchange AB MN with BA NM
        * addInverse : bool
            add additional inverse measurements
        * spacing : float [1]
            electrode spacing in meters
        * closed : bool
            Close the chain. Measure from the end of the array to the first
            electrode.

    Returns
    -------
    data : DataContainerERT

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from pygimli.physics import ert
    >>>
    >>> schemes = ['wa', 'wb', 'pp', 'pd', 'dd', 'slm', 'hw', 'gr']
    >>> fig, ax = plt.subplots(3,3)
    >>>
    >>> for i, schemeName in enumerate(schemes):
    ...     s = ert.createData(elecs=41, schemeName=schemeName)
    ...     k = ert.createGeometricFactors(s)
    ...     _ = ert.show(s, vals=k, ax=ax.flat[i], label='k - ' + schemeName)
    >>>
    >>> plt.show()
    """
    if kwargs.pop('sounding', False):
        data = pg.DataContainerERT()
        data.setSensors(pg.cat(-elecs[::-1], elecs))

        nElecs = len(elecs)
        for i in range(nElecs-1):
            data.createFourPointData(i, i, 2*nElecs-i-1, nElecs-1, nElecs)

        return data

    mg = DataSchemeManager()

    if schemeName == "none":
        pg.error('argument "schemeName" not set. Valid schemeNames are:')
        for i in mg.schemes():
            print(i, "scheme: " + mg.scheme(i).prefix)

    scheme = mg.scheme(schemeName)

    scheme.setInverse(kwargs.pop('inverse', False))
    scheme.addInverse(kwargs.pop('addInverse', False))
    scheme._closed = kwargs.pop('closed', False)

    if isinstance(elecs, int):
        data = scheme.create(nElectrodes=elecs,
                             electrodeSpacing=kwargs.pop('spacing', 1),
                             **kwargs)
    elif hasattr(elecs, '__iter__'):
        if isinstance(elecs[0], (float, int)):
            data = scheme.create(nElectrodes=len(elecs), **kwargs)
            data.setSensors(elecs)
        else:
            data = scheme.create(sensorList=elecs, **kwargs)
    else:
        print(elecs)
        pg.critical("Can't interpret elecs")

    data['k'] = ert.createGeometricFactors(data)
    return data


def createDataVES(ab2, mn2):
    """ Utility one-liner to create a BERT datafile for Schlumberger 1D VES

    Parameters
    ----------
    ab2: array
        Half distance between current electrodes

    mn2: float
        Half distance between measurement electrodes

    Returns
    -------
    data : DataContainerERT
    """

    data = pg.DataContainerERT()

    if isinstance(mn2, (float, int)):
        mn2 = [mn2]

    count = 0
    for mn in mn2:
        emID = data.createSensor([-mn, 0.0, 0.0])
        enID = data.createSensor([mn, 0.0, 0.0])

        for x in ab2:
            eaID = data.createSensor([-x, 0.0, 0.0])
            ebID = data.createSensor([x, 0.0, 0.0])
            data.createFourPointData(count, eaID, ebID, emID, enID)
            count += 1

    data.fitFillSize()

    return data


class Pseudotype:
    unknown = 0
    A_M = 1
    AB_MN = 2
    AB_M = 3
    AB_N = 4
    DipoleDipole = 5
    Schlumberger = 6
    WennerAlpha = 7
    WennerBeta = 8
    Gradient = 9
    PoleDipole = 10
    HalfWenner = 11
    PolePole = 12
    Test = 99


class DataSchemeManager(object):
    """Data scheme manager."""
    def __init__(self):
        """Initialize all the data schemes."""
        self.schemes_ = {}
        self.addScheme(DataSchemeBase())
        self.addScheme(DataSchemeWennerAlpha())
        self.addScheme(DataSchemeWennerBeta())
        self.addScheme(DataSchemeDipoleDipole())
        self.addScheme(DataSchemeSchlumberger())
        self.addScheme(DataSchemePolePole())
        self.addScheme(DataSchemePoleDipole())
        self.addScheme(DataSchemeHalfWenner())
        self.addScheme(DataSchemeMultipleGradient())

        self.addScheme(DataSchemeBase(typ=Pseudotype.A_M, name='A_M'))
        self.addScheme(DataSchemeBase(typ=Pseudotype.AB_MN, name='AB_MN'))
        self.addScheme(DataSchemeBase(typ=Pseudotype.AB_M, name='AB_M'))
        self.addScheme(DataSchemeBase(typ=Pseudotype.AB_N, name='AB_N'))

    def addScheme(self, scheme):
        """A a scheme from given name."""
        self.schemes_[scheme.name] = scheme

    def scheme(self, name):
        """
            Return DataScheme for a given name if registered.

        Parameters
        ----------
        name : str | int

            Name or prefix name of a known data scheme. If the name is unknown
            all known data schemes are listed.
            Name can be a integer number that
            represents the internal Pseudotype.

        Return
        ------

        scheme : DataScheme

        """
        if isinstance(name, int):
            s = self.schemeFromTyp(name)
            if s:
                return s

        elif isinstance(name, str):  # or type(name) == unicode: (always in Py3)
            s = self.schemeFromPrefix(name)
            if s:
                return s

            if name in self.schemes_:
                return self.schemes_[name]

            print('Unknown scheme name:', name)
            print('-----------------------')
            print('Valid names or prefixes')
            print('-----------------------')
            for s in self.schemes_.values():

                print(s.name, ': ', s.prefix)
            raise Exception("No scheme known for name: ", name)

        return DataSchemeBase()

    def schemeFromPrefix(self, prefix):
        """ Return DataScheme for a given prefix name.
        """
        for s in list(self.schemes_.values()):
            if s.prefix == prefix:
                return s
        return None

    def schemeFromTyp(self, typ):
        for s in list(self.schemes_.values()):
            if s.type == typ:
                return s
        return None

    def schemes(self):
        '''.'''
        return list(self.schemes_.keys())


class DataSchemeBase(object):
    """Base class for ERT data schemes

    Attributes
    ----------
    closed : bool
        Close the chain. Measure from the end of the array to the first
        electrode.

    """
    def __init__(self, typ=Pseudotype.unknown, name="unknown", prefix='uk'):
        self.name = name
        self.prefix = prefix
        self.type = typ
        self.data_ = None
        self.inverse_ = False
        self.addInverse_ = False
        self.reciprocity = False
        self.nElectrodes_ = 0
        self.maxSeparation = 1e99
        self._closed = False

    @property
    def closed(self):
        return self._closed

    def create(self, nElectrodes=24, electrodeSpacing=1, sensorList=None,
               **kwargs):
        """."""
        self.createElectrodes(nElectrodes, electrodeSpacing, sensorList)
        self.setMaxSeparation(kwargs.pop("maxSeparation", 999))
        self.createData(**kwargs)

        if self.addInverse_:
            out = pg.DataContainerERT(self.data_)
            self.setInverse(not self.inverse_)
            self.createData(**kwargs)
            self.data_.add(out)
            self.data_.removeInvalid()
            self.data_.sortSensorsIndex()

        if kwargs.values():
            print("Warning! DataSchemeBase::create has unhandled arguments")
            print(kwargs)

        return self.data_

    def createElectrodes(self, nElectrodes=24, electrodeSpacing=1,
                         sensorList=None):
        self.data_ = pg.DataContainerERT()

        if sensorList is not None:
            for p in sensorList:
                if isinstance(p, float):
                    self.data_.createSensor((p, 0.))
                else:
                    self.data_.createSensor(p)
        else:
            for i in range(nElectrodes):
                self.data_.createSensor(pg.Pos(float(i) *
                                        electrodeSpacing, 0.0))

        self.nElectrodes_ = self.data_.sensorCount()

    def createData(self, **kwargs):
        print('*'*100)

    def setInverse(self, inverse=False):
        self.inverse_ = inverse

    def addInverse(self, addInverse=False):
        """Add inverse value to create a full dataset."""
        self.addInverse_ = addInverse

    def setMaxSeparation(self, maxSep):
        if maxSep > 0.0:
            self.maxSeparation = maxSep
        else:
            self.maxSeparation = 1e99

    def createDatum_(self, a, b, m, n, count):
        if a < self.nElectrodes_ and b < self.nElectrodes_ and \
                m < self.nElectrodes_ and n < self.nElectrodes_:
            if self.inverse_:
                self.data_.createFourPointData(count, m, n, a, b)
            else:
                self.data_.createFourPointData(count, a, b, m, n)

            count += 1
        return count


class DataSchemePolePole(DataSchemeBase):
    """Pole-Pole data scheme."""
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Pole Pole (C-P)"
        self.prefix = "pp"
        self.type = Pseudotype.PolePole

    def createData(self, **kwargs):
        """
        Create a Pole-Pole dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='pp', **kwargs) instead.
        """
        nElectrodes = self.nElectrodes_
        # reserve a couple more than nesseccary ###
        self.data_.resize((nElectrodes) * (nElectrodes))

        count = 0
        # enlargeEverySep = 0  # not used

        b = -1
        n = -1

        for a in range(0, nElectrodes):
            for m in range(a + 1, nElectrodes):
                if m - a > self.maxSeparation:
                    break

                count = self.createDatum_(a, b, m, n, count)

        self.data_.removeInvalid()
        return self.data_


class DataSchemeDipoleDipole(DataSchemeBase):
    """Dipole-dipole data scheme. """
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Dipole Dipole (CC-PP)"
        self.prefix = "dd"
        self.type = Pseudotype.DipoleDipole
        self.enlargeEverySep = 0
        self.spacings = [1]

    def createData(self, **kwargs):
        """
        Create a Dipole-Dipole dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='dd', **kwargs) instead.

        Parameters
        ----------

        **kwargs:
            * complete : bool
                Add reciprocity measurements.
            * enlarge : int
                Enlarge dipole length every n dipole separations.
            * spacings : array[int]
                vector of spacings (dipole lengths) to use
        """
        nElectrodes = self.nElectrodes_
        complete = kwargs.pop('complete', False)

        if complete:
            self._closed = True

        self.enlargeEverySep = kwargs.pop('enlarge', 0)
        self.spacings = kwargs.pop('spacings', self.spacings)
        # self.createElectrodes(nElectrodes, electrodeSpacing)

        # reserve a couple more than necessary ###
        nElectrodes = self.nElectrodes_
        self.data_.resize(nElectrodes * nElectrodes)

        count = 0

        if self.closed:
            space = 1
            for i in range(nElectrodes):
                a = i
                b = (a + space) % nElectrodes

                for j in range(nElectrodes):
                    m = (j) % nElectrodes
                    n = (m + space) % nElectrodes

                    if not complete:
                        if j <= i:
                            continue

                    if a != m and a != n and b != m and b != n:
                        count = self.createDatum_(a, b, m, n, count)

        else:
            for space in self.spacings:
                maxSep = nElectrodes - space
                maxInj = nElectrodes - space

                if self.maxSeparation < maxSep:
                    maxSep = self.maxSeparation
                for sep in range(1, maxSep + 1):

                    if self.enlargeEverySep > 0:
                        if (sep-1) % self.enlargeEverySep == 0:
                            space += 1

                    for i in range(maxInj - sep):
                        a = i
                        b = (a + space) % nElectrodes
                        m = (b + sep) % nElectrodes
                        n = (m + space) % nElectrodes
                        if m + space < nElectrodes:
                            count = self.createDatum_(a, b, m, n, count)

        self.data_.removeInvalid()
        return self.data_
# class DataSchemeDipoleDipole

class DataSchemePoleDipole(DataSchemeBase):
    """Pole-dipole data scheme"""
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Pole Dipole (C-PP)"
        self.prefix = "pd"
        self.type = Pseudotype.PoleDipole
        self.spacings = [1]

    def createData(self, **kwargs):
        """
        Create a Pole-Dipole dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='pd', **kwargs) instead.

        Parameters
        ----------

        **kwargs:
            * enlarge : int
                Enlarge dipole length every n dipole separations.
            * spacings : array[int]
                vector of spacings (dipole lengths) to use
        """
        nElectrodes = self.nElectrodes_

        # self.createElectrodes(nElectrodes, electrodeSpacing)
        # reserve a couple more than nesseccary !!!
        self.data_.resize((nElectrodes) * (nElectrodes))

        count = 0
        self.enlargeEverySep = kwargs.pop('enlarge', 0)
        self.spacings = kwargs.pop('spacings', self.spacings)

        b = -1
        for a in range(0, nElectrodes):

            for m in range(a + 1, nElectrodes - 1):
                n = m + 1

                if m - a > self.maxSeparation:
                    break

                count = self.createDatum_(a, b, m, n, count)

        self.data_.removeInvalid()
        return self.data_
# class DataSchemePoleDipole

class DataSchemeHalfWenner(DataSchemeBase):
    """Pole-Dipole like Wenner Beta with increasing dipole distance"""
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Half Wenner (C-P-P)"
        self.prefix = "hw"
        self.type = Pseudotype.HalfWenner

    def createData(self, **kwargs):
        """
        Create a Half-Wenner dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='hw', **kwargs) instead.
        """
        nElectrodes = self.nElectrodes_

        # reserve a couple more than nesseccary !!!
        self.data_.resize((nElectrodes) * (nElectrodes))

        # print("create", self.maxSeparation)

        count = 0
        # space = 0  # not yet used
        # enlargeEverySep = 0  # not yet used

        b = -1
        for a in range(0, nElectrodes):
            inc = 1
            while True:
                m = a - inc
                n = m - inc

                if m < 0 or n < 0 or inc > self.maxSeparation:
                    break

                count = self.createDatum_(a, b, m, n, count)
                inc = inc + 1

            inc = 1
            while True:
                m = a + inc
                n = m + inc

                if m > nElectrodes or n > nElectrodes or \
                        inc > self.maxSeparation:
                    break

                count = self.createDatum_(a, b, m, n, count)
                inc = inc + 1

        self.data_.removeInvalid()
        self.data_.sortSensorsIndex()
        return self.data_
#class DataSchemeHalfWenner

class DataSchemeWennerAlpha(DataSchemeBase):
    """Wenner alpha (C--P--P--C) data scheme with equal distances. """
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Wenner Alpha (C-P-P-C)"
        self.prefix = "wa"
        self.type = Pseudotype.WennerAlpha

    def createData(self, **kwargs):
        """Create a Wenner-alpha dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='wa', **kwargs) instead.
        """
        nElectrodes = self.nElectrodes_

        maxSep = nElectrodes - 2
        if self.maxSeparation < maxSep:
            maxSep = self.maxSeparation

        # reserve a couple more than nesseccary !!!
        self.data_.resize(nElectrodes * nElectrodes)

        count = 0

        for sep in range(1, maxSep + 1):
            for i in range((nElectrodes - 2) - sep):
                a = i
                m = a + sep
                n = m + sep
                b = n + sep
                count = self.createDatum_(a, b, m, n, count)

        self.data_.removeInvalid()
        return self.data_
#class DataSchemeWennerAlpha

class DataSchemeWennerBeta(DataSchemeBase):
    """Wenner-beta (C--C--P--P) data scheme with equal distance."""
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Wenner Beta(C-C-P-P)"
        self.prefix = "wb"
        self.type = Pseudotype.WennerBeta

    def createData(self, **kwargs):
        """Create a Wenner-beta dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='wb', **kwargs) instead.
        """
        nElectrodes = self.nElectrodes_

        maxSep = nElectrodes - 2
        if self.maxSeparation < maxSep:
            maxSep = self.maxSeparation

        # reserve a couple more than nesseccary ###
        self.data_.resize((nElectrodes * nElectrodes))

        count = 0

        for sep in range(1, maxSep + 1):
            for i in range((nElectrodes - 2) - sep):
                a = i
                b = a + sep
                m = b + sep
                n = m + sep

                count = self.createDatum_(a, b, m, n, count)

        self.data_.removeInvalid()
        return self.data_
# class DataSchemeWennerBeta(...)


class DataSchemeSchlumberger(DataSchemeBase):
    """Wenner-Schlumberger (C--P-P--C) data scheme. """
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "Schlumberger(C-PP-C)"
        self.prefix = "slm"
        self.type = Pseudotype.Schlumberger

    def createData(self, **kwargs):
        """Create a full (Wenner-)Schlumberger dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='sl', **kwargs) instead.
        """
        nElectrodes = self.nElectrodes_

        maxSep = nElectrodes - 2
        if self.maxSeparation < maxSep:
            maxSep = self.maxSeparation

        self.data_.resize(nElectrodes * nElectrodes)

        count = 0

        for sep in range(1, maxSep + 1):
            for i in range((nElectrodes - 2) - sep):
                a = i
                m = a + sep
                n = m + 1
                b = n + sep
                count = self.createDatum_(a, b, m, n, count)

        self.data_.removeInvalid()
        return self.data_
# class DataSchemeSchlumberger(...)


class DataSchemeMultipleGradient(DataSchemeBase):
    """MultipleGradient (C---P-P--C) data scheme. """
    def __init__(self):
        DataSchemeBase.__init__(self)
        self.name = "MultipleGradient(C--P-P--C)"
        self.prefix = "gr"
        self.type = Pseudotype.Gradient

    def createData(self, **kwargs):
        """Create a multi-gradient dataset.

        Don't use directly .. call create from DataSchemeManager or
        ert.createData(elecs, schemeName='gr', **kwargs) instead.
        """
        nElectrodes = self.nElectrodes_

        ab_sep_base = 9  # number of channels + 2
        ev = 2
        takeevery = 2
        max_fak = int(np.ceil(nElectrodes / ab_sep_base))
        ab_space = [ii*ab_sep_base for ii in range(max_fak) if ii % ev == 1]
        mn_space = [ii for ii in range(max_fak) if ii % takeevery == 1]
        count = 0
        a, b, m, n = [], [], [], []
        for ab in range(len(ab_space)):  # ab spacings
            for aa in np.arange(1, nElectrodes-ab_space[ab]+1, 1):  # a index
                mn = mn_space[ab]
                for mm in np.arange(aa+mn, aa+ab_space[ab]-mn, mn):
                    count += 1
                    a.append(int(aa))
                    b.append(int(aa+ab_space[ab]))
                    m.append(int(mm))
                    n.append(int(mm+mn))

        self.data_.resize(count)
        self.data_.set('a', pg.Vector(a) - 1)
        self.data_.set('b', pg.Vector(b) - 1)
        self.data_.set('m', pg.Vector(m) - 1)
        self.data_.set('n', pg.Vector(n) - 1)
        self.data_.set('valid', pg.Vector(count, 1))

        return self.data_
# class DataSchemeMultipleGradient(...)


if __name__ == '__main__':
    schemes = ['wa', 'wb', 'pp', 'pd', 'dd', 'slm', 'gr', 'hw']
    fig, ax = pg.plt.subplots(3, 3)
    kw = dict(cMin=10, cMax=1000, logScale=True, colorBar=False, cMap="viridis")
    for it, scheme in enumerate(schemes):
        shm = ert.createData(elecs=41, schemeName=scheme)
        print(scheme, shm)
        k = ert.createGeometricFactors(shm)
        mgr = DataSchemeManager()
        longname = mgr.scheme(scheme).name
        ert.show(shm, vals=np.abs(k), ax=ax.flat[it], colorBar=1, logScale=0,
                label='k ' + longname + ')-' + scheme)

    createColorBarOnly(**kw, ax=ax.flat[-1], aspect=0.1)
    pg.plt.show()
