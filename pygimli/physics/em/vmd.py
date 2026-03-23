#!/usr/bin/env python
"""Electromagnetics for Vertical Magnetic Dipole (VMD) functions and class."""

from math import sqrt, pi

import numpy as np
import pygimli as pg

from pygimli.frameworks import Block1DModelling


class VMDModelling(Block1DModelling):
    r"""Modelling operator for a Vertical Magnetic Dipole (VMD).

    Modelling operator for a Vertical Magnetic Dipole (VMD) to calculate the
    electromagnetic response in cylindric coordinates
    ::math::`H_z, H_r, E_{\phi} \in I\!C`
    for a layered halfspace (1D) using Hankel transformation.

    The VMD is at the origin ::math::`r_s = (0.0)` and ::math::`z_s=0`.
    """

    def __init__(self, **kwargs):
        """Initialize forward operator.

        Initialize forward operator with optional halfspace geometry
        and measurement spezifications.

        Parameters
        ----------
        **kwargs : dict()
            Forward to ModellingBase
        """
        super().__init__(**kwargs)
        self.nLayers = kwargs.pop("nLayers", 0)

    def createStartModel(self, rhoa):
        r"""Create suitable starting model.

        Create suitable starting model based on median apparent
        resistivity values and skin depth approximation.
        """
        if self.nLayers == 0:
            pg.critical("Model space is not been initialized.")

        skinDepth = np.sqrt(max(self.t) * pg.math.median(rhoa)) * 500
        thk = np.arange(self.nLayers) / sum(np.arange(self.nLayers)) * \
            skinDepth / 2.
        startThicks = thk[1:]

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks, trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa), trans='log')
        return super().regionManager().createStartModel()

    def response_mt(self, par, i=0):
        """Compute response vector for a set of model parameter.

        Parameters
        ----------
        par : iterabale
            DOCUMENTME
        """
        raise pg.critical("IMPLEMENTME")

    def response(self, par):
        """Return model response."""
        return self.response_mt(par)

    def calcEPhiF(self, f, rho, d, rmin=1, nr=41, ze=0, zs=0, tm=1):
        """Compute radial E field from VMD source.

        Parameters
        ----------
        f : float
            Frequency
        rho : iterable
            resistivity vector
        d : iterable
            thickness vector
        rmin : float
            minimum radious to compute
        nr : int
            number of radii
        ze : float [0]
            z coordinate of receiver in Meter
        zs : float [ze]
            z coordinate of source in Meter
        zs    z-Koordinate des Senders in Meter
        """
        # filter coefficients
        q = 10**0.1
        rr = rmin * q**(np.array(nr) - 1)

        he = ze
        hs = zs
        # zm = hs - he  # not used
        zp = hs + he

        # rm = np.sqrt(rr**2 + zm**2)  # not used
        rp = np.sqrt(rr**2 + zp**2)

        fcJ1, nc0 = pg.utils.hankelFC(4)
        nc = len(fcJ1)

        ncnr = nc + nr - 1

        # Vakuumwerte, nicht normiert
        if ze <= 0:
            ePhi = rr / rp**3

        # Create Wavenumbers
        nu = np.arange(1, ncnr + 1)

        n = nc0 - nc + nu
        q = np.log(10) * 0.1
        k = np.exp(-(n-1) * q) / rmin

        if ze <= 0:
            ePhi = rr / rp**3.0

        # Admittanzen for halfspace borders for each Wavenumbers k
        bt = np.zeros(ncnr) * 1j
        aa = np.zeros(ncnr) * 1j
        for i in range(ncnr):
            if ze <= 0:
                bt[i] = self.btp(k[i], f, rho, d, type=1)
            else:
                raise RuntimeError('NeedTests')
                # not used (uncommented)
                # aa[i], aap[i], bt[i] = downward(k[i], f, rho, d, ze)

        # Kernel functions
        if ze <= 0:
            e = np.exp(k * ze) * np.exp(k * zs)
            delta = e * (bt - k) / (bt + k)
        else:
            raise RuntimeError('NeedTests')
            e = np.exp(k * zs)
            delta = 2. * k**2. / (bt + k)*e

        # convolution
        aux3 = np.zeros(nr) * 1j

        for n in range(nr):
            for nn in range(nc):
                nu = nn + n
                mn = nc - nn

                nnn = nc0 - nc + nu
                k = np.exp(-nnn * q) / rmin

                del0 = delta[nu]

                if ze <= 0:
                    del1 = del0 * k
                    aux3[n] = aux3[n] + del1 * fcJ1[mn-1]
                else:
                    raise RuntimeError('NeedTests')
                    del3 = del0 * aa(nu)
                    aux3[n] = aux3[n] + del3 * fcJ1[mn]

            if nr > 1:
                aux3[n] = aux3[n] / rr[n]
                raise ValueError("more than one r .. check code here")
            else:
                aux3[n] = aux3[n] / rr

        if ze <= 0:
            # Air
            ePhi = ePhi - aux3
        else:
            raise ValueError("more than one r .. check code here")
            # Halfspace
            ePhi = aux3

        # Normalization
        ePhi *= -tm * 1j * (2*pi*f) * pg.physics.constants.mu0 / (4*pi)
        return ePhi

    def btp(self, k, f, rho, d, type=1):
        """Admittance of a layered halfspace.

        for TE (TYP=1) and TM (TYP=2) mode.
        k: wave number (1/m)
        f: frequency (1/s)
        rho: layer resitivities
        d: layer thicknesses
        """
        nl = len(rho)
        c = 1j * pg.physics.constants.mu0 * 2 * pi * f
        b = np.sqrt(k**2 + c/rho[nl-1])

        if nl > 1:
            beta = 1
            for nn in range(nl-2, -1, -1):
                alpha = np.sqrt(k**2 + c / rho[nn])
                if type == 2:
                    beta = rho[nn] / rho[nn + 1]

                cth = np.exp(-2. * d[nn] * alpha)
                cth = (1-cth) / (1+cth)
                b = (b + alpha * beta * cth) / (beta + cth * b / alpha)

        return b


class VMDTimeDomainModelling(VMDModelling):
    """Vertical magnetic dipole (VMD) modelling."""

    def __init__(self, times, txArea, rxArea=None, **kwargs):
        """Set up class with necessary arguments."""
        super().__init__(**kwargs)

        self.t = times
        self.txArea = txArea

        if rxArea is None:
            self.rxArea = txArea
        else:
            self.rxArea = rxArea

    def createStartModel(self, rhoa, nLayers=None, thickness=None):
        """Create suitable starting model.

        Create suitable starting model based on median apparent resistivity
        values and skin depth approximation.
        """
        if nLayers is None:
            nLayers = self.nLayers

        res = np.ones(nLayers) * pg.math.median(rhoa)
        if thickness is None:
            skinDepth = np.sqrt(max(self.t) * pg.math.median(rhoa)) * 500
            thk = np.arange(nLayers) / sum(np.arange(nLayers)) * skinDepth / 2.
            thk = thk[1:]
        else:
            thk = np.ones(nLayers-1) * thickness

        self.setStartModel(pg.cat(thk, res))
        return self.startModel()

    def response_mt(self, par, i=0):
        """Multi-threading result.

        par = [thicknesses, res]
        """
        nLay = (len(par)-1)//2
        thk = par[0:nLay]
        res = par[nLay:]
        return self.calcRhoa(thk, res)

    def response(self, par):
        """Return forward response.

        par = [thicknesses(nLay), res(nlay + 1)]
        """
        return self.response_mt(par, 0)

    def calcRhoa(self, thk, res):
        """Compute apparent resistivity response."""
        a = sqrt(self.txArea / pi)  # TX coil radius

        ePhiTD, tD = self.calcEphiT(tMin=min(self.t), tMax=max(self.t),
                                    rho=res, d=thk, rMin=a, rMax=a, z=0,
                                    dipm=self.rxArea)

        ePhi = np.exp(np.interp(np.log(self.t),
                                np.log(tD), np.log(ePhiTD[:, 0])))

        tmp = a**(4./3) * self.rxArea**(2./3) * \
            pg.physics.constants.mu0**(5./3) / (20**(2./3) * pi**(1./3))
        rhoa = tmp / (self.t**(5./3) * (ePhi * 2 * pi * a)**(2./3))

        return rhoa

    def calcEphiT(self, tMin, tMax, rho, d, rMin, rMax, z, dipm):
        """Compute radial electric field."""
        # nl = len(rho)
        r = pg.utils.niceLogspace(rMin, rMax, nDec=10)
        r = [rMin]

        t = pg.utils.niceLogspace(tMin, tMax, nDec=10)

        fcS, nc0 = pg.utils.hankelFC(1)

        nr = len(r)
        nc = len(fcS)
        nt = len(t)
        ncnt = nc + nt

        fef = np.zeros((ncnt, nr), complex)
        ePhi = np.zeros((nt, nr))

        omega = 10. ** (0.1 * (1 - (-nc + nc0 + np.arange(1, ncnt)))) / t[0]

        for nn, o in enumerate(omega):
            f = o / (2. * pi)
            ePhiF = self.calcEPhiF(f, rho, d, ze=z, zs=0.,
                                   rmin=r[0], nr=len(r), tm=1.0)

            fef[nn, :] = ePhiF / sqrt(o)

            ita = max(0, nn-nc)
            ite = min(nt, nn+1)
            for it in range(ita, ite):
                itn = nc-nn+it-1
                ePhi[it] = ePhi[it] - np.real(fef[nn]) * fcS[itn]

        for it in range(nt):
            ePhi[it] = ePhi[it] * dipm * sqrt(2/pi / t[it])

        return ePhi, t
