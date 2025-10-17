#!/usr/bin/env python
"""Vertical electrical sounding (VES) manager class."""
import numpy as np

import pygimli as pg
from pygimli.frameworks import Block1DModelling
from pygimli.frameworks.modelling import DEFAULT_STYLES


class VESModelling(Block1DModelling):
    """Vertical Electrical Sounding (VES) forward operator.

    Parameters
    ----------
    ab2 : array-like
        Half distance between the current electrodes A and B.
    mn2 : array-like
        Half distance between the potential electrodes M and N.
        Only used for input (feeding am etc.) or plotting.
    am : array-like
        Part of data basis. Distances between A and M electrodes.
        A is first current, M is first potential electrode.
    bm : array-like
        Part of data basis. Distances between B and M electrodes.
        B is second current, M is first potential electrode.
    an : array-like
        Part of data basis. Distances between A and N electrodes.
        A is first current, N is second potential electrode.
    bn : array-like
        Part of data basis. Distances between B and N electrodes.
        B is second current, N is second potential electrode.
    nLayers : int [4]
        number of layers, i.e. the parameter number is nLayers*2-1
    nPara : int []
        number of parameters per layer
    """

    def __init__(self, ab2=None, mn2=None, **kwargs):
        """Initialize with distances.

        Either with
        * all distances AM, AN, BM, BN using am/an/bm/bn
        * a dataContainerERT using data or dataContainerERT
        * AB/2 and (optionally) MN/2 distances using ab2/mn2

        nLayers : int [4]
            Number of layers.
        """
        self.am = kwargs.pop("am", None)
        self.bm = kwargs.pop("bm", None)
        self.an = kwargs.pop("an", None)
        self.bn = kwargs.pop("bn", None)
        super().__init__(**kwargs)
        self.ab2 = ab2
        self.mn2 = mn2


        if 'dataContainerERT' in kwargs or 'data' in kwargs:
            if 'data' in kwargs:
                data = kwargs['data']
            else:
                data = kwargs['dataContainerERT']
            if isinstance(data, pg.DataContainerERT):
                kwargs['am'] = [data.sensorPosition(data('a')[i]).distance(
                    data('m')[i]) for i in range(data.size())]
                kwargs['an'] = [data.sensorPosition(data('a')[i]).distance(
                    data('n')[i]) for i in range(data.size())]
                kwargs['bm'] = [data.sensorPosition(data('b')[i]).distance(
                    data('m')[i]) for i in range(data.size())]
                kwargs['bn'] = [data.sensorPosition(data('b')[i]).distance(
                    data('n')[i]) for i in range(data.size())]

        self.setDataSpace(ab2=ab2, mn2=mn2,
                          am=self.am, an=self.an, bm=self.bm, bn=self.bn)

    def createStartModel(self, rhoa):
        """Create starting model."""
        if self.nLayers == 0:
            pg.critical("Model space is not been initialized.")

        startThicks = np.logspace(np.log10(min(self.mn2)/2),
                                  np.log10(max(self.ab2)/5),
                                  self.nLayers - 1)
        startThicks = pg.utils.diff(pg.cat([0.0], startThicks))

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks, trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa), trans='log')

        return super().createStartModel()

    def setDataSpace(self, ab2=None, mn2=None,
                     am=None, bm=None, an=None, bn=None,
                     **kwargs):
        """Set data basis, i.e., arrays for all am, an, bm, bn distances.

        You can set either
        * AB/2 and (optionally) MN/2 spacings for a classical sounding, or
        * all distances AM, AN, BM, BN for arbitrary arrays

        Parameters
        ----------
        ab2 : iterable
            AB/2 distances
        mn2 : iterable | float
            MN/2 distance(s)
        am, an, bm, bn : iterable
            distances between current and potential electrodes
        """
        # Sometimes you don't have AB2/MN2 but provide am etc.
        self.am = am
        self.an = an
        self.bm = bm
        self.bn = bn
        if ab2 is not None:
            if mn2 is None:
                mn2 = np.array(ab2) / 3

            if isinstance(mn2, float):
                mn2 = np.ones(len(ab2))*mn2

            if len(ab2) != len(mn2):
                print("ab2", ab2)
                print("mn2", mn2)
                raise ValueError("length of ab2 is unequal length of mn2")

            self.am = ab2 - mn2
            self.an = ab2 + mn2
            self.bm = ab2 + mn2
            self.bn = ab2 - mn2

        elif (am is not None and bm is not None and an is not None and
              bn is not None):
            self.am = am
            self.bm = bm
            self.an = an
            self.bn = bn

        if self.am is not None and self.bm is not None:
            self.ab2 = (self.am + self.bm) / 2
            self.mn2 = abs(self.am - self.an) / 2

            self.k = (2.0 * np.pi) / (1.0 / self.am - 1.0 / self.an -
                                      1.0 / self.bm + 1.0 / self.bn)

    def response(self, par):
        """Model response."""
        return self.response_mt(par, 0)

    def response_mt(self, par, i=0):
        """Multi-threading model response."""
        if self.am is not None and self.bm is not None:
            nLayers = (len(par)+1) // 2
            fop = pg.core.DC1dModelling(nLayers,
                                        self.am, self.bm, self.an, self.bn)
        else:
            pg.critical("No data space defined don't know what to calculate.")

        return fop.response(par)

    def drawModel(self, ax, model, **kwargs):
        """Draw model as 1D block model."""
        pg.viewer.mpl.drawModel1D(ax=ax,
                                  model=model,
                                  plot=kwargs.pop('plot', 'loglog'),
                                  xlabel=r'Resistivity ($\Omega$m)',
                                  **kwargs)
        ax.set_ylabel('Depth in (m)')

        return ax, None  # should return gci and not ax

    def drawData(self, ax, data, error=None, label=None, **kwargs):
        r"""Draw modeled apparent resistivity data.

        Parameters
        ----------
        ax: axes
            Matplotlib axes object to draw into.

        data: iterable
            Apparent resistivity values to draw.

        error: iterable [None]
            Adds an error bar if you have error values.

        label: str ['$\rho_a$']
            Set legend label for the amplitude.

        Other Parameters
        ----------------
        ab2: iterable
            Override ab2 that fits data size.
        mn2: iterable
            Override mn2 that fits data size.
        plot: function name
            Matplotlib plot function, e.g., plot, loglog, semilogx or semilogy
        """
        ab2 = kwargs.pop('ab2', self.ab2)
        # mn2 = kwargs.pop('mn2', self.mn2)
        plot = getattr(ax, kwargs.pop('plot', 'loglog'))

        ra = data
        raE = error

        style = dict(pg.frameworks.modelling.DEFAULT_STYLES.get(
            label, pg.frameworks.modelling.DEFAULT_STYLES['Default']))
        style.update(kwargs)

        if label is None:
            label = r'$\rho_a$'

        plot(ra, ab2, label=label, **style)

        if raE is not None:
            raErr = np.array(ra * raE)

            if pg.isArray(raErr, len(ra)):
                ax.errorbar(ra, ab2,
                            xerr=raErr, barsabove=True,
                            **DEFAULT_STYLES.get('Error',
                                                 DEFAULT_STYLES['Default']),
                            label='_nolegend_')

        ax.set_ylim(max(ab2), min(ab2))
        ax.set_xlabel(r'Apparent resistivity ($\Omega$m)')
        ax.set_ylabel(r'AB/2 (m)')
        ax.grid(True)
        ax.legend()
        return ax, None  # should return gci and not ax&cb


class VESCModelling(VESModelling):
    """Vertical Electrical Sounding (VES) complex forward operator.

    Vertical Electrical Sounding (VES) forward operator for complex
    resistivity values. see: :py:mod:`pygimli.physics.ert.VESModelling`
    """

    def __init__(self, **kwargs):
        """Initialize with number of parameters and distances."""
        super().__init__(nPara=2, **kwargs)
        self.phiAxe = None

    def phaseModel(self, model):
        """Return the current phase model values."""
        nLay = (len(model) + 1) // 3
        return pg.cat(model[0:nLay-1], 1000. * model[nLay*2-1::])

    def resModel(self, model):
        """Return the resistivity model values."""
        nLay = (len(model) + 1) // 3
        return model[0:nLay*2-1]

    def createStartModel(self, rhoa):
        """Create starting model of nlay-1 thicknesses & nlay resistivities."""
        startDepths = np.logspace(np.log10(min(self.mn2)/2),
                                  np.log10(max(self.ab2)/5),
                                  self._nLayers-1)
        startThicks = pg.utils.diff(pg.cat([0.0], startDepths))

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks,
                                 trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa),
                                 trans='log')

        self.setRegionProperties(2, startModel=np.median(rhoa[len(rhoa)//2::]),
                                 trans='log')

        sm = self.regionManager().createStartModel()
        return sm

    def response_mt(self, par, i=0):
        """Multi-threaded model response for parametrization.

        Returns
        -------
        response : iterable
            [rhoa, +phi(rad)] for [thicks, res, phi(rad)]
        """
        if self.am is not None and self.bm is not None:
            nLayers = (len(par) + 1) // 3
            fop = pg.core.DC1dModellingC(
                nLayers, self.am, self.bm, self.an, self.bn)
        else:
            pg.critical("No data basis known.")

        return fop.response(par)

    def drawModel(self, ax, model, **kwargs):
        """Draw 1D VESC Modell."""
        a1 = ax
        a2 = pg.viewer.mpl.createTwinY(ax)

        super().drawModel(
            a1, model=self.resModel(model), **kwargs)

        plot = kwargs.pop('plot', 'semilogy')
        if plot == 'loglog':
            plot = 'semilogy'
        elif plot == 'semilogx':
            plot = 'plot'

        pg.viewer.mpl.drawModel1D(ax=a2,
                                  model=self.phaseModel(model),
                                  plot=plot,
                                  color='C2',
                                  xlabel='Phase (mrad)',
                                  **kwargs)

        a2.set_xlabel('neg. phase (mRad)', color='C2')

    def drawData(self, ax, data, error=None, labels=None, ab2=None, mn2=None,
                 **kwargs):
        r"""Draw modeled apparent resistivity and apparent phase data.

        Parameters
        ----------
        ax: axes
            Matplotlib axes object to draw into.

        data: iterable
            Apparent resistivity values to draw. [rhoa phia].

        error: iterable [None]
            Rhoa in Ohm m and phia in radiand.
            Adds an error bar if you have error values. [err_rhoas err_phia]
            The error of amplitudes are assumed to be relative and the error
            of the phases is assumed to be absolute in mrad.

        labels: str [r'$\varrho_a$', r'$\varphi_a$']
            Set legend labels for amplitude and phase.

        Other Parameters
        ----------------
        ab2: iterable
            Override ab2 that fits data size.
        mn2: iterable
            Override mn2 that fits data size.
        plot: function name
            Matplotlib plot function, e.g., plot, loglog, semilogx or semilogy
        """
        a1 = None
        a2 = None

        if hasattr(ax, '__iter__'):
            if len(ax) == 2:
                a1 = ax[0]
                a2 = ax[1]
        else:
            a1 = ax
            a2 = pg.viewer.mpl.createTwinY(ax)

        if ab2 is not None and mn2 is not None:
            self.setDataSpace(ab2=ab2, mn2=mn2)

        ra = data[0:len(data)//2]
        phi = data[len(data)//2::] * 1000.  # mRad

        phiE = None  # abs err
        raE = None  # rel err

        if error is not None:
            if type(error) is float:
                raE = np.ones(len(data)//2) * error
                phiE = np.ones(len(data)//2) * error
            else:
                raE = error[0:len(data)//2]
                phiE = error[len(data)//2::]

        if labels is None:
            labels = [r'$\varrho_a$', r'$\varphi_a$']

        label = kwargs.pop('label', 'Data')

        style = dict(pg.frameworks.modelling.DEFAULT_STYLES.get(
            label, pg.frameworks.modelling.DEFAULT_STYLES['Default']))
        style.update(kwargs)

        super().drawData(a1, ra, error=raE,
                                            label=labels[0], **style)

        style['color'] = 'C2'

        a2.semilogy(phi, self.ab2, label=labels[1], **style)

        if phiE is not None:
            a2.errorbar(phi, self.ab2,
                        xerr=phiE,
                        **DEFAULT_STYLES.get('Error',
                                             DEFAULT_STYLES['Default']),
                        barsabove=True,
                        label='_nolegend_'
                        )

        a2.set_ylim(max(self.ab2), min(self.ab2))
        a2.set_xlabel('Apparent neg. phase (mRad)', color='C2')
        a2.set_ylabel('AB/2 in (m)')
        a2.legend()
        a2.grid(True)


# class VESRhoModelling(VESModelling):  # not working due to Block1dModelling
class VESRhoModelling(pg.frameworks.MeshModelling):
    """Vertical electrical sounding (VES) modelling with fixed layers."""

    def __init__(self, thk, verbose=False, **kwargs):
        """Initialize modelling operator by passing model and data space.

        Parameters
        ----------
        thk : iterable, optional
            Thickness vector of the individual layers.
        verbose : bool, optional
            some output. The default is False.
        **kwargs : geometric definition of the sounding, either
            ab2 : iterable
                AB/2 distances
            mn2 : iterable
                MN/2 distances (if not specified, ab2/3 by default) OR
            am : iterable
                A-M distance AND
            an : iterable
                A-N distance AND
            bm : iterable
                N-M distance AND
            bn : iterable
                B-N distance OR
            dataContainer : pg.DataContainerERT
                ERT data container to determine the AM/AN/BM/BN distances
        """
        super().__init__(verbose=verbose)
        # better do the following in a function like setDataSpace/setModelSpace
        self.bfop = VESModelling(**kwargs)
        self.thk = thk
        self.fwd = pg.core.DC1dRhoModelling(thk, self.bfop.am, self.bfop.bm,
                                            self.bfop.an, self.bfop.bn,
                                            verbose=verbose)
        self.mesh_ = pg.meshtools.createMesh1D(len(thk)+1)
        self.setMesh(self.mesh_)
        # self.mesh = self.mesh_  # could work with MeshModelling parent

    def response(self, par):
        """Forward response (app. resistivity for given resistivity vector)."""
        return self.fwd.response(par)

    def response_mt(self, par):
        """Forward response."""
        fwd = pg.core.DC1dRhoModelling(self.thk, self.bfop.am, self.bfop.bm,
                                       self.bfop.an, self.bfop.bn, False)
        return fwd.response(par)

    def createStartModel(self, rhoa):
        """Create starting model."""
        return pg.Vector(len(self.thk)+1, np.median(rhoa))
