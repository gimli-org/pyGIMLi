#!/usr/bin/env python
"""Vertical electrical sounding (VES) manager class."""
import numpy as np
import pygimli as pg
from pygimli.frameworks import MethodManager1d
from .vesModelling import VESModelling, VESCModelling


class VESManager(MethodManager1d):
    r"""Vertical electrical sounding (VES) manager class.

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.physics import VESManager
    >>>
    >>> ab2 = np.logspace(np.log10(1.5), np.log10(100), 32)
    >>> mn2 = 1.0
    >>> # 3 layer with 100, 500 and 20 Ohmm
    >>> # and layer thickness of 4, 6, 10 m
    >>> # over a half-space of 800 Ohmm
    >>> synthModel = pg.cat([4., 6., 10.], [100., 5., 20., 800.])
    >>> ves = VESManager()
    >>> ra, err = ves.simulate(synthModel, ab2=ab2, mn2=mn2, noiseLevel=0.01)
    >>> ax = ves.showData(ra, error=err)
    >>> model = ves.invert(ra, err, nLayers=4, showProgress=0, verbose=0)
    >>> ax, _ = ves.showModel(synthModel)
    >>> _ = ves.showResult(ax=ax)
    """

    def __init__(self, **kwargs):
        """Initialize instance.

        Keyword Arguments
        -----------------
        **kwargs:
            complex: bool
                Accept complex resistivities.

        """
        self._complex = kwargs.pop('complex', False)

        super().__init__(**kwargs)

        self.inv.setDeltaPhiStop(1)

        self.dataTrans = None
        self.rhoaTrans = pg.trans.TransLog()
        self.phiaTrans = pg.trans.TransLin()
        if "ab2" in kwargs or "data" in kwargs or "dataContainer" in kwargs:
            self.fop.setDataSpace(**kwargs)

    @property
    def complex(self):
        """Return whether the computations are complex."""
        return self._complex

    @complex.setter
    def complex(self, c):
        self._complex = c
        self.reinitForwardOperator()

    def createForwardOperator(self, **kwargs):
        """Create Forward Operator.

        Create Forward Operator based on complex attribute.
        """
        if self.complex:
            return VESCModelling(**kwargs)
        else:
            return VESModelling(**kwargs)

    def simulate(self, model, ab2=None, mn2=None, **kwargs):
        """Simulate measurement data."""
        if ab2 is not None and mn2 is not None:
            self._fw.fop.setDataSpace(ab2=ab2, mn2=mn2)

        return super().simulate(model, **kwargs)

    def preErrorCheck(self, err, dataVals=None):
        """Fct to be called before the validity check of the error values."""
        err = np.atleast_1d(err)
        if self.complex:
            if len(err) == 2:
                nData = len(dataVals) // 2
                err = pg.cat(np.ones(nData)*err[0],
                             np.abs(err[1] / dataVals[nData:]))
        else:
            if len(err) == 1:
                err = np.ones(nData)*err[0]

        return err

    def invert(self, data=None, err=None, ab2=None, mn2=None, **kwargs):
        """Invert measured data.

        Parameters
        ----------
        data : iterable
            data vector
        err : iterable
            error vector
        ab2: iterable
            AB/2 vector (otherwise taken from data)
        mn2: iterable
            MN/2 vector (otherwise taken from data)

        Keyword Arguments
        ----------------
        **kwargs
            Additional kwargs inherited from %(MethodManager1d.invert) and
            %(Inversion.run)

        Returns
        -------
        model : pg.Vector
            inversion result
        """
        if ab2 is not None and mn2 is not None:
            self.fop.setDataSpace(ab2=ab2, mn2=mn2)

        if data is not None:
            if self.complex:
                nData = len(data)//2
                self.dataTrans = pg.trans.TransCumulative()
                self.dataTrans.add(self.rhoaTrans, nData)
                self.dataTrans.add(self.phiaTrans, nData)
            else:
                self.dataTrans = pg.trans.TransLog()

            self.inv.dataTrans = self.dataTrans

        if 'layerLimits' not in kwargs:
            kwargs['layerLimits'] = [min(self.fop.mn2)/5,
                                     max(self.fop.ab2)/2]

        if 'paraLimits' in kwargs and self.complex:
            pL = kwargs['paraLimits'][1]
            kwargs['paraLimits'][1] = [pL[0]/1000, pL[1]/1000]

        return super().invert(data=data, err=err, **kwargs)

    def loadData(self, fileName, **kwargs):
        """Load simple data matrix."""
        mat = np.loadtxt(fileName)
        if len(mat[0]) == 4:
            self.fop.setDataSpace(ab2=mat[:, 0], mn2=mat[:, 1])
            return mat.T
        if len(mat[0]) == 6:
            self.complex = True
            self.fop.setDataSpace(ab2=mat[:, 0], mn2=mat[:, 1])
            return (mat[:, 0], mat[:, 1],
                    np.array(pg.cat(mat[:, 2], mat[:, 4])),
                    np.array(pg.cat(mat[:, 3], mat[:, 5])))

    def exportData(self, fileName, data=None, error=None):
        """Export data into simple ascii matrix.

        Usefull?
        """
        mn2 = np.abs((self.fop.am - self.fop.an) / 2.)
        ab2 = (self.fop.am + self.fop.bm) / 2.
        mat = None
        if data is None:
            data = self.inv.dataVals

        if error is None:
            error = self.inv.errorVals

        if self.complex:
            nData = len(data)//2
            mat = np.array([ab2, mn2,
                            data[:nData], error[:nData],
                            data[nData:], error[nData:]
                            ]).T
            np.savetxt(fileName, mat,
                       header=r'ab/2\tmn/2\trhoa\terr\tphia\terrphi')
        else:
            mat = np.array([ab2, mn2, data, error]).T
            np.savetxt(fileName, mat, header=r'ab/2\tmn/2\trhoa\terr')


def VESManagerApp():
    """Call VESManager as console app."""
    parser = VESManager.createArgParser(dataSuffix='ves')
    options = parser.parse_args()

    verbose = not options.quiet
    if verbose:
        print("VES Manager console application.")
        print(options._get_kwargs())

    mgr = VESManager(verbose=verbose, debug=pg.debug())

    ab2, mn2, ra, err = mgr.loadData(options.dataFileName)

    mgr.showData(ra, err)
    mgr.invert(ra, err, ab2, mn2,
               maxIter=options.maxIter,
               lam=options.lam,
               )
    mgr.showResultAndFit()


if __name__ == '__main__':
    VESManagerApp()
