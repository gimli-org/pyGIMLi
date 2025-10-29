#!/usr/bin/env python
"""Modelling helicopter electromagnetics (HEM) using VMD solvers."""
from math import sqrt, pi
import numpy as np

import pygimli as pg
from pygimli.physics.constants import Constants
from pygimli.frameworks import Block1DModelling, MeshModelling
from .tools import FDEMsystems

class HEMmodelling(Block1DModelling):
    """HEM Airborne modelling class based on the BGR RESOLVE system."""

    ep0 = Constants.e0
    mu0 = Constants.mu0
    c0 = sqrt(1. / ep0 / mu0)
    fdefault = np.array([387.0, 1821.0, 8388.0, 41460.0, 133300.0], float)
    rdefault = np.array([7.94, 7.93, 7.93, 7.91, 7.92], float)
    scaling = 1e6

    def __init__(self, nlay, height, f=None, r=None, **kwargs):
        """Initialize class with geometry.

        Parameters
        ----------
        nlay : int
            number of layers
        height : float
            helicopter
        f : array [BGR RESOLVE system 387Hz-133kHz]
            frequency vector
        r : array [BGR RESOLVE system 7.91-7.94]
            distance vector
        scaling : float
            scaling factor or string (ppm=1e6, percent=1e2)
        """
        self.nlay = nlay
        self.height = height
        syst = kwargs.pop('system', None)
        if isinstance(syst, str):
            f, r = FDEMsystems.getFrequenciesDistances(syst)

        self.f = np.asarray(f)
        if r is None:
            raise ValueError("Specify separation value or vector!")
        if 'scaling' in kwargs:
            if kwargs['scaling'] == 'ppm':
                self.scaling = 1e6
            elif kwargs['scaling'] in ['percent', '%']:
                self.scaling = 1e2
            else:
                self.scaling = kwargs['scaling']

        if self.f is None:
            self.f = self.fdefault
        if isinstance(r, (float, int)):
            self.r = np.ones_like(f, dtype=float) * r
        else:
            if len(r) == len(self.f):
                self.r = r
            else:
                raise ValueError('Length vector have to be matching!')
        if self.r is None:
            self.r = self.rdefault

        self.wem = (2.0 * pi * self.f) ** 2 * self.ep0 * self.mu0
        self.iwm = 1.0j * 2.0 * pi * self.f * self.mu0
        mesh = pg.meshtools.createMesh1DBlock(nlay)
        super().__init__()
        self.setMesh(mesh)

    def response(self, model):
        """Compute response vector by pasting in-phase and out-phase data."""
        ip, op = self.vmd_hem(self.height,
                              np.asarray(model[self.nlay-1:self.nlay*2-1]),
                              np.asarray(model[:self.nlay-1]))
        return pg.cat(ip, op)

    def response_mt(self, par, i=0):
        """Multi-threaded forward response."""
        return self.response(par)

    def calc_forward(self, x, h, rho, d, epr, mur, quasistatic=False):
        """Calculate forward response."""
        field = np.zeros((self.f.size, x.size), complex)
        # Forward calculation for background model
        if d.size:
            for m in range(x.size):
                field[:, m] = self.vmd_hem(np.array([h[m]], float),
                                           rho[:, m], d[:, m], epr[:, m],
                                           mur[:, m], quasistatic).T[:, 0]
        else:
            for n in range(self.f.size):
                for m in range(x.size):
                    field[n, m] = self.vmd_hem(np.array([h[n, m]], float),
                                               np.array([rho[n, m]], float),
                                               d,
                                               np.array([epr[n, m]], float),
                                               np.array([mur[n, m]], float),
                                               quasistatic).T[n, 0]
        return field

    def downward(self, rho, d, z, epr, mur, lam):
        """Downward continuation of fields."""
        nl = rho.size
        alpha = np.zeros((nl, lam.shape[1], self.f.size), complex)
        b = np.zeros((nl, lam.shape[1], self.f.size), complex)
        aa = np.zeros((nl, lam.shape[1], self.f.size), complex)
        aap = np.zeros((nl, lam.shape[1], self.f.size), complex)
        rho = rho[:, np.newaxis, np.newaxis] * np.ones(
            (rho.size, lam.shape[1], self.f.size), float)
        d = d[:, np.newaxis, np.newaxis] * np.ones(
            (d.size, lam.shape[1], self.f.size), float)
        h = np.insert(np.cumsum(d[:, 0, 0]), 0, 0)
        epr = epr[:, np.newaxis, np.newaxis] * np.ones(
            (epr.size, lam.shape[1], self.f.size), float)
        mur = mur[:, np.newaxis, np.newaxis] * np.ones(
            (mur.size, lam.shape[1], self.f.size), float)
        lam = np.tile(lam, (nl, 1, 1))
        # progression constant
        alpha = np.sqrt(lam ** 2 - np.tile(self.wem, (nl, lam.shape[1], 1)) *
                        epr * mur + np.tile(self.iwm, (nl, lam.shape[1], 1)) *
                        mur / rho)
        if nl == 1:  # homogenous halfspace
            b1 = alpha.copy()  # (1, 100, nfreq)
            a = np.exp(-alpha * z)
            ap = a.copy()
            return b1, a, ap
        elif nl > 1:  # multi-layer case
            # tanh num unstable tanh(x)=(exp(x)-exp(-x))/(exp(x)+exp(-x))
            ealphad = np.exp(-2.0 * alpha[0:-1, :, :] * d)
            talphad = (1.0 - ealphad) / (1.0 + ealphad)
            b[-1, :, :] = np.copy(alpha[-1, :, :])
            # recursive admittance computation at upper layer boundary
            # from bottom to top, for nl-1 layers
            for n in range(nl-2, -1, -1):
                b[n, :, :] = alpha[n, :, :] * \
                    (b[n+1, :, :] + alpha[n, :, :] * talphad[n, :, :]) / \
                    (alpha[n, :, :] + b[n+1, :, :] * talphad[n, :, :])
            # Impedance
            c = 1.0 / b
            # b1 == 1. row in b (nl, 100, nfreq)
            b1 = np.copy(b[0, :, :][np.newaxis, :, :])  # (1, 100, nfreq)
            # Variation from one layer boundary to the other
            for n in range(0, nl-1):
                aa[n, :, :] = (b[n, :, :] + alpha[n, :, :]) / (
                    b[n+1, :, :] + alpha[n, :, :]) * \
                    np.exp(-alpha[n, :, :] * d[n, :, :])
                aap[n, :, :] = (1.0 + alpha[n, :, :] * c[n, :, :]) / (
                    1.0 + alpha[n, :, :] * c[n+1, :, :]) * \
                    np.exp(-alpha[n, :, :] * d[n, :, :])
            # Determine layer Index where z is
            ind = nl - 1
            for n in range(0, nl-1):
                if np.logical_and(z >= h[n], z < h[n+1]):
                    ind = n

            if (ind + 1) < nl:
                a = np.prod(aa[:ind, :, :], 0) * 0.5 * \
                    (1.0 + b[ind, :, :] / alpha[ind, :, :]) * \
                    (np.exp(-alpha[ind, :, :] * (z - h[ind])) -
                     (b[ind+1, :, :] - alpha[ind, :, :]) /
                     (b[ind+1, :, :] + alpha[ind, :, :]) *
                     np.exp(-alpha[ind, :, :] *
                            (d[ind, :, :] + h[ind+1] - z)))
                ap = np.prod(aap[:ind, :, :], 0) * 0.5 * \
                    (1.0 + alpha[ind, :, :] * c[ind, :, :]) * \
                    (np.exp(-alpha[ind, :, :] * (z - h[ind])) +
                        (1.0 - alpha[ind, :, :] * c[ind+1, :, :]) /
                        (1.0 + alpha[ind, :, :] * c[ind+1, :, :]) *
                        np.exp(-alpha[ind, :, :] *
                               (d[ind, :, :] + h[ind+1] - z)))
            else:
                a = np.prod(aa, 0) * np.exp(-alpha[ind, :, :] * (z - h[ind]))
                ap = np.prod(aap, 0) * np.exp(-alpha[ind, :, :] * (z - h[ind]))
            a = a[np.newaxis, :, :]  # (1, 100, nfreq)
            ap = ap[np.newaxis, :, :]  # (1, 100, nfreq)
            return b1, a, ap

    def vmd_hem(self, h, rho, d, epr=1., mur=1., quasistatic=False):
        """Vertical magnetic dipole (VMD) response.

        Parameters
        ----------
        h : float
            flight height
        rho : array
            resistivity vector
        d : array
            thickness vector
        """
        # filter coefficients
        if isinstance(epr, float):
            epr = np.ones((len(rho),), float)*epr
        if isinstance(mur, float):
            mur = np.ones((len(rho),), float)*mur
        fc0, nc, nc0 = hankelfc(3)
        fc1, nc, nc0 = hankelfc(4)
        # allocate arrays
        nf = len(self.f)
        lam = np.zeros((1, nc, nf), float)
        alpha0 = np.zeros((1, nc, nf), complex)
        delta0 = np.zeros((1, nc, nf), complex)
        delta1 = np.zeros((1, nc, nf), complex)
        delta2 = np.zeros((1, nc, nf), complex)
        delta3 = np.zeros((1, nc, nf), complex)
        aux0 = np.zeros((1, nf), complex)
        aux1 = np.zeros((1, nf), complex)
        aux2 = np.zeros((1, nf), complex)
        aux3 = np.zeros((1, nf), complex)
        Z = np.zeros(nf, complex)
        # r0
        r0 = np.copy(self.r)
        # determine optimum r0 (shift nodes) for f > 1e4 and h > 100
        if quasistatic:
            index = np.zeros(self.f.shape, np.bool)
        else:
            index = np.logical_and(self.f >= 1e4, h >= 100.0)

        if np.any(index):
            opt = np.floor(10.0 * np.log10(
                self.r[index] * 2.0 * np.pi * self.f[index] / self.c0) + nc0)
            r0[index] = self.c0 / (2.0 * np.pi * self.f[index]) * 10.0 ** (
                (opt + 0.5 - nc0) / 10.0)
        # Wave numbers
        n = np.arange(nc0 - nc, nc0, 1, float)
        q = 0.1 * np.log(10)
        lam = np.reshape(np.exp(-n[np.newaxis, :, np.newaxis] * q) /
                         r0[np.newaxis, np.newaxis, :], (-1, nc, r0.size))
        # (1, 100, nfreq)
        # wave number in air, quasistationary approximation
        alpha0 = np.copy(lam) * complex(1, 0)  # (1, 100, nfreq)
        # wave number in air, full solution for f > 1e4
        index = np.zeros(self.f.shape, np.bool) if quasistatic \
            else self.f >= 10000.0
        if np.any(index):
            alpha0[:, :, index] = np.sqrt(
                lam[:, :, index]**2 - np.tile(self.wem[index], (1, nc, 1)) +
                np.tile(self.iwm[index], (1, nc, 1)) / 1e9)  # (1, 100 , nfreq)
        # Admittanzen an der Oberfläche eines geschichteten Halbraums
        b1, _, _ = self.downward(rho, d, 0.0, epr, mur, lam)
        # Kernel functions
        e = np.exp(-2.0 * h * alpha0)  # (1, 100, nfreq)
        delta0 = (b1 - alpha0 * mur[0]) / (b1 + alpha0 * mur[0]) * e
        delta1 = (2 * mur[0]) / (b1 + alpha0 * mur[0]) * e  # (1, 100, nfreq)
        delta2 = 1 / h * e  # (1, 100, nfreq)
        delta3 = 1 / (2 * h) * e  # (1, 100, nfreq)
        # convolution
        # quasistationary approximation
        aux0 = np.sum(delta0 * lam ** 3 / alpha0 *
                      np.tile(fc0[::-1].T[:, :, np.newaxis],
                              (1, 1, self.f.size)), 1, complex) / r0
        # full solution, partial integration
        if np.any(index):
            aux1 = np.sum(delta1 * lam ** 3 *
                          np.tile(fc0[::-1].T[:, :, np.newaxis],
                                  (1, 1, self.f.size)), 1, complex) / r0
            aux2 = np.sum(
                delta2 * lam * np.tile(fc0[::-1].T[:, :, np.newaxis],
                                       (1, 1, self.f.size)), 1, complex)/r0
            aux3 = np.sum(delta3 * lam ** 2 *
                          np.tile(fc1[::-1].T[:, :, np.newaxis],
                                  (1, 1, self.f.size)), 1, complex) / r0
        # normed secondary field
        # quasistationary approximation
        Z = self.r ** 3 * aux0 * self.scaling
        # full solution
        if np.any(index):
            Z[:, index] = (-self.r[index]**3 * aux1[:, index] +
                           self.r[index]**3 * aux2[:, index] -
                           self.r[index]**4 * aux3[:, index]) * self.scaling
        return np.real(Z[0]), np.imag(Z[0])

    def vmd_total_Ef(self, h, z, rho, d, epr, mur, tm):
        """VMD E-phi field (not used actively)."""
        # only halfspace
        # Filter coefficients
        fc1, nc, nc0 = hankelfc(4)
        lam = np.zeros((1, nc, self.f.size), float)
        alpha0 = np.zeros((1, nc, self.f.size), complex)
        delta = np.zeros((1, nc, self.f.size), complex)
        aux = np.zeros((1, self.f.size), complex)
        Ef = np.zeros(self.f.size, complex)
        r0 = np.copy(self.r)
        # wave numbers
        n = np.arange(nc0 - nc, nc0, 1, float)
        q = 0.1 * np.log(10)
        lam = np.reshape(np.exp(-n[np.newaxis, :, np.newaxis] * q) /
                         r0[np.newaxis, np.newaxis, :], (-1, nc, r0.size))
        # wave numbers in air, full solution
        alpha0 = np.sqrt(lam ** 2 - np.tile(self.wem, (1, nc, 1)) +
                         np.tile(self.iwm, (1, nc, 1)) / 1e9)
        # admittances on surface of layered halfspace
        b1, aa, _ = self.downward(rho, d, z, epr, mur, lam)
        # Kernel functions
        e = np.exp(-h * alpha0)  # (1, 100, nfreq)
        delta = 2.0 / (alpha0 + b1) * e  # (1, 100, nfreq)
        # convolution
        # quasistationary approximation
        aux = np.sum(delta*lam**2*aa*np.tile(fc1[::-1].T[:, :, np.newaxis],
                                             (1, 1, self.f.size)), 1,
                     complex) / r0  # (1, nfreq)
        # absolute fields, full solution
        Ef = -tm * self.iwm / (4.0 * np.pi) * aux
        return Ef


def hankelfc(order):
    """Filter coefficients for Hankel transformation."""
    if order == 1:  # sin
        fc = np.array([
            2.59526236e-7, 3.66544843e-7, 5.17830795e-7, 7.31340622e-7,
            1.03322805e-6, 1.45918500e-6, 2.06161065e-6, 2.91137793e-6,
            4.11357863e-6, 5.80876420e-6, 8.20798075e-6, 1.15895083e-5,
            1.63778560e-5, 2.31228459e-5, 3.26800649e-5, 4.61329334e-5,
            6.52101085e-5, 9.20390575e-5, 1.30122935e-4, 1.83620431e-4,
            2.59656626e-4, 3.66311982e-4, 5.18141184e-4, 7.30717340e-4,
            1.03392184e-3, 1.45742714e-3, 2.06292302e-3, 2.90599911e-3,
            4.11471902e-3, 5.79042763e-3, 8.20004722e-3, 1.15192930e-2,
            1.63039133e-2, 2.28257757e-2, 3.22249222e-2, 4.47864328e-2,
            6.27329625e-2, 8.57059100e-2, 1.17418314e-1, 1.53632655e-1,
            1.97717964e-1, 2.28849849e-1, 2.40311038e-1, 1.65409220e-1,
            2.84701476e-3, -2.88016057e-1, -3.69097406e-1, -2.50107514e-2,
            5.71811256e-1, -3.92261572e-1, 7.63280044e-2, 5.16233994e-2,
            -6.48012082e-2, 4.89047141e-2, -3.26936331e-2, 2.10539842e-2,
            -1.33862549e-2, 8.47124695e-3, -5.35123972e-3, 3.37796651e-3,
            -2.13174466e-3, 1.34513833e-3, -8.48749612e-4, 5.35531006e-4,
            -3.37898780e-4, 2.13200109e-4, -1.34520273e-4, 8.48765787e-5,
            -5.35535069e-5, 3.37899801e-5, -2.13200365e-5, 1.34520337e-5,
            -8.48765949e-6, 5.35535110e-6, -3.37899811e-6, 2.13200368e-6,
            -1.34520338e-6, 8.48765951e-7, -5.35535110e-7, 3.37899811e-7],
            float)
        nc = 80
        nc0 = 40
    elif order == 2:  # cos
        fc = np.array([
            1.63740363e-7, 1.83719709e-7, 2.06136904e-7, 2.31289411e-7,
            2.59510987e-7, 2.91176117e-7, 3.26704977e-7, 3.66569013e-7,
            4.11297197e-7, 4.61483045e-7, 5.17792493e-7, 5.80972733e-7,
            6.51862128e-7, 7.31401337e-7, 8.20645798e-7, 9.20779729e-7,
            1.03313185e-6, 1.15919300e-6, 1.30063594e-6, 1.45933752e-6,
            1.63740363e-6, 1.83719709e-6, 2.06136904e-6, 2.31289411e-6,
            2.59510987e-6, 2.91176117e-6, 3.26704977e-6, 3.66569013e-6,
            4.11297197e-6, 4.61483045e-6, 5.17792493e-6, 5.80972733e-6,
            6.51862128e-6, 7.31401337e-6, 8.20645798e-6, 9.20779729e-6,
            1.03313185e-5, 1.15919300e-5, 1.30063594e-5, 1.45933752e-5,
            1.63740363e-5, 1.83719709e-5, 2.06136904e-5, 2.31289411e-5,
            2.59510987e-5, 2.91176117e-5, 3.26704977e-5, 3.66569013e-5,
            4.11297197e-5, 4.61483045e-5, 5.17792493e-5, 5.80972733e-5,
            6.51862128e-5, 7.31401337e-5, 8.20645798e-5, 9.20779729e-5,
            1.03313185e-4, 1.15919300e-4, 1.30063594e-4, 1.45933752e-4,
            1.63740363e-4, 1.83719709e-4, 2.06136904e-4, 2.31289411e-4,
            2.59510987e-4, 2.91176117e-4, 3.26704976e-4, 3.66569013e-4,
            4.11297197e-4, 4.61483045e-4, 5.17792493e-4, 5.80972733e-4,
            6.51862127e-4, 7.31401337e-4, 8.20645797e-4, 9.20779730e-4,
            1.03313185e-3, 1.15919300e-3, 1.30063593e-3, 1.45933753e-3,
            1.63740362e-3, 1.83719710e-3, 2.06136901e-3, 2.31289411e-3,
            2.59510977e-3, 2.91176115e-3, 3.26704948e-3, 3.66569003e-3,
            4.11297114e-3, 4.61483003e-3, 5.17792252e-3, 5.80972566e-3,
            6.51861416e-3, 7.31400728e-3, 8.20643673e-3, 9.20777603e-3,
            1.03312545e-2, 1.15918577e-2, 1.30061650e-2, 1.45931339e-2,
            1.63734419e-2, 1.83711757e-2, 2.06118614e-2, 2.31263461e-2,
            2.59454421e-2, 2.91092045e-2, 3.26529302e-2, 3.66298115e-2,
            4.10749753e-2, 4.60613861e-2, 5.16081994e-2, 5.78193646e-2,
            6.46507780e-2, 7.22544422e-2, 8.03873578e-2, 8.92661837e-2,
            9.80670729e-2, 1.07049506e-1, 1.13757572e-1, 1.18327217e-1,
            1.13965041e-1, 1.00497783e-1, 6.12958082e-2, -1.61234222e-4,
            -1.11788551e-1, -2.27536948e-1, -3.39004453e-1, -2.25128800e-1,
            8.98279919e-2, 5.12510388e-1, -1.31991937e-1, -3.35136479e-1,
            3.64868100e-1, -2.34039961e-1, 1.32085237e-1, -7.56739672e-2,
            4.52296662e-2, -2.78297002e-2, 1.73727753e-2, -1.09136894e-2,
            6.87397283e-3, -4.33413470e-3, 2.73388730e-3, -1.72477355e-3,
            1.08821012e-3, -6.86602007e-4, 4.33213523e-4, -2.73338487e-4,
            1.72464733e-4, -1.08817842e-4, 6.86594042e-5, -4.33211523e-5,
            2.73337984e-5, -1.72464607e-5, 1.08817810e-5, -6.86593962e-6,
            4.33211503e-6, -2.73337979e-6, 1.72464606e-6, -1.08817810e-6,
            6.86593961e-7, -4.33211503e-7, 2.73337979e-7, -1.72464606e-7],
            float)
        nc = 164
        nc0 = 122
    elif order == 3:  # J0
        fc = np.array([
            2.89878288e-7, 3.64935144e-7, 4.59426126e-7, 5.78383226e-7,
            7.28141338e-7, 9.16675639e-7, 1.15402625e-6, 1.45283298e-6,
            1.82900834e-6, 2.30258511e-6, 2.89878286e-6, 3.64935148e-6,
            4.59426119e-6, 5.78383236e-6, 7.28141322e-6, 9.16675664e-6,
            1.15402621e-5, 1.45283305e-5, 1.82900824e-5, 2.30258527e-5,
            2.89878259e-5, 3.64935186e-5, 4.59426051e-5, 5.78383329e-5,
            7.28141144e-5, 9.16675882e-5, 1.15402573e-4, 1.45283354e-4,
            1.82900694e-4, 2.30258630e-4, 2.89877891e-4, 3.64935362e-4,
            4.59424960e-4, 5.78383437e-4, 7.28137738e-4, 9.16674828e-4,
            1.15401453e-3, 1.45282561e-3, 1.82896826e-3, 2.30254535e-3,
            2.89863979e-3, 3.64916703e-3, 4.59373308e-3, 5.78303238e-3,
            7.27941497e-3, 9.16340705e-3, 1.15325691e-2, 1.45145832e-2,
            1.82601199e-2, 2.29701042e-2, 2.88702619e-2, 3.62691810e-2,
            4.54794031e-2, 5.69408192e-2, 7.09873072e-2, 8.80995426e-2,
            1.08223889e-1, 1.31250483e-1, 1.55055715e-1, 1.76371506e-1,
            1.85627738e-1, 1.69778044e-1, 1.03405245e-1, -3.02583233e-2,
            -2.27574393e-1, -3.62173217e-1, -2.05500446e-1, 3.37394873e-1,
            3.17689897e-1, -5.13762160e-1, 3.09130264e-1, -1.26757592e-1,
            4.61967890e-2, -1.80968674e-2, 8.35426050e-3, -4.47368304e-3,
            2.61974783e-3, -1.60171357e-3, 9.97717882e-4, -6.26275815e-4,
            3.94338818e-4, -2.48606354e-4, 1.56808604e-4, -9.89266288e-5,
            6.24152398e-5, -3.93805393e-5, 2.48472358e-5, -1.56774945e-5,
            9.89181741e-6, -6.24131160e-6, 3.93800058e-6, -2.48471018e-6,
            1.56774609e-6, -9.89180896e-7, 6.24130948e-7, -3.93800005e-7,
            2.48471005e-7, -1.56774605e-7, 9.89180888e-8, -6.24130946e-8],
            float)
        nc = 100
        nc0 = 60
    elif order == 4:  # J1
        fc = np.array([
            1.84909557e-13, 2.85321327e-13, 4.64471808e-13, 7.16694771e-13,
            1.16670043e-12, 1.80025587e-12, 2.93061898e-12, 4.52203829e-12,
            7.36138206e-12, 1.13588466e-11, 1.84909557e-11, 2.85321327e-11,
            4.64471808e-11, 7.16694771e-11, 1.16670043e-10, 1.80025587e-10,
            2.93061898e-10, 4.52203829e-10, 7.36138206e-10, 1.13588466e-9,
            1.84909557e-9, 2.85321326e-9, 4.64471806e-9, 7.16694765e-9,
            1.16670042e-8, 1.80025583e-8, 2.93061889e-8, 4.52203807e-8,
            7.36138149e-8, 1.13588452e-7, 1.84909521e-7, 2.85321237e-7,
            4.64471580e-7, 7.16694198e-7, 1.16669899e-6, 1.80025226e-6,
            2.93060990e-6, 4.52201549e-6, 7.36132477e-6, 1.13587027e-5,
            1.84905942e-5, 2.85312247e-5, 4.64449000e-5, 7.16637480e-5,
            1.16655653e-4, 1.79989440e-4, 2.92971106e-4, 4.51975783e-4,
            7.35565435e-4, 1.13444615e-3, 1.84548306e-3, 2.84414257e-3,
            4.62194743e-3, 7.10980590e-3, 1.15236911e-2, 1.76434485e-2,
            2.84076233e-2, 4.29770596e-2, 6.80332569e-2, 9.97845929e-2,
            1.51070544e-1, 2.03540581e-1, 2.71235377e-1, 2.76073871e-1,
            2.16691977e-1, -7.83723737e-2, -3.40675627e-1, -3.60693673e-1,
            5.13024526e-1, -5.94724729e-2, -1.95117123e-1, 1.99235600e-1,
            -1.38521553e-1, 8.79320859e-2, -5.50697146e-2, 3.45637848e-2,
            -2.17527180e-2, 1.37100291e-2, -8.64656417e-3, 5.45462758e-3,
            -3.44138864e-3, 2.17130686e-3, -1.36998628e-3, 8.64398952e-4,
            -5.45397874e-4, 3.44122545e-4, -2.17126585e-4, 1.36997597e-4,
            -8.64396364e-5, 5.45397224e-5, -3.44122382e-5, 2.17126544e-5,
            -1.36997587e-5, 8.64396338e-6, -5.45397218e-6, 3.44122380e-6,
            -2.17126543e-6, 1.36997587e-6, -8.64396337e-7, 5.45397218e-7],
            float)
        nc = 100
        nc0 = 60
    return (np.reshape(fc, (-1, 1)), nc, nc0)  # (100,) -> (100, 1)


class HEMRhoModelling(HEMmodelling):
    """Airborne EM (HEM) Forward modelling class for Occam inversion."""

    def __init__(self, dvec, height=1., **kwargs):
        """Init class by specifying frequencies and distances (s. HEMMod)."""
        nlay = len(dvec) + 1
        self.dvec = np.asarray(dvec)
        self.zvec = np.hstack((0, np.cumsum(dvec)))
        HEMmodelling.__init__(self, nlay, height, **kwargs)
        self.mymesh = pg.meshtools.createMesh1D(nlay)
        self.setMesh(self.mymesh)  # only for inversion

    def response(self, model):
        """Forward response as combined in-phase and out-of-phase."""
        ip, op = self.vmd_hem(self.height, rho=model, d=self.dvec)
        return pg.cat(ip, op)


class FDEMResSusModelling(HEMmodelling):
    """FDEM block modelling class using both conductivity & permittivity."""

    def __init__(self, nlay, height=1, **kwargs):
        """Init class (see HEMmodelling)."""
        HEMmodelling.__init__(self, nlay, height, **kwargs)
        self.scaling = 1e2
        self.mymesh = pg.meshtools.createMesh1DBlock(nlay, nProperties=2)
        self.setMesh(self.mymesh)  # only for inversion
        # pg.core.ModellingBase.__init__(self, self.mymesh)

    def response(self, model):
        """Response vector as combined in-phase and out-phase data."""
        thk = np.asarray(model[:self.nlay-1], dtype=float)
        res = np.asarray(model[self.nlay-1:2*self.nlay-1], dtype=float)
        mur = np.asarray(model[2*self.nlay-1:3*self.nlay-1], dtype=float) + 1
        ip, op = self.vmd_hem(self.height, rho=res, d=thk, mur=mur)
        return pg.cat(ip, op)


class HEMRhoSusModelling(HEMmodelling):
    """Airborne EM (HEM) smooth forward modelling including susceptibility."""

    def __init__(self, dvec, *args, **kwargs):
        """Initialize (not yet working)."""
        self.nlay = len(dvec) + 1
        self.dvec = np.asarray(dvec)
        self.zvec = np.hstack((0, np.cumsum(dvec)))
        HEMmodelling.__init__(self, self.nlay, *args, **kwargs)
        self.mymesh = pg.meshtools.createMesh1D(self.nlay, nProperties=2)
        self.setMesh(self.mymesh)  # only for inversion

    def response(self, model):
        """Response vector as combined in-phase and out-phase data."""
        res = np.asarray(model[:self.nlay])
        mur = np.asarray(model[self.nlay:]) + 1
        ip, op = self.vmd_hem(self.height, rho=res, d=self.dvec, mur=mur)
        return pg.cat(ip, op)


class FDEMLCIFOP(pg.core.ModellingBase):
    """FDEM 2d-LCI modelling class based on Block matrices."""

    def __init__(self, data, nlay=2, verbose=False, f=None, r=None):
        """Parameters: FDEM data class and number of layers."""
        super().__init__(verbose)
        self.nlay = nlay
        self.FOP = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        self.np = 2 * nlay - 1
        self.mesh2d = pg.meshtools.createMesh2D(self.np, self.nx)
        self.mesh2d.rotate(pg.RVector3(0, 0, -np.pi/2))
        self.setMesh(self.mesh2d)

        self.J = pg.matrix.BlockMatrix()
        self.FOP1d = []
        for i in range(self.nx):
            self.FOP1d.append(HEMmodelling(nlay, data.z[i], f, r))
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, self.nf*2*i, self.np*i)

        self.setJacobian(self.J)

    def response(self, model):
        """Cut together forward responses of all soundings."""
        modA = np.reshape(model, (self.nx, self.np))
        resp = pg.Vector(0)
        for i, modi in enumerate(modA):
            resp = pg.cat(resp, self.FOP1d[i].response(modi))

        return resp

    def createJacobian(self, model):
        """Fill the individual blocks of the Block-Jacobi matrix."""
        modA = np.reshape(model, (self.nx, self.np))
        for i, modi in enumerate(modA):
            self.FOP1d[i].createJacobian(modi)


class FDEM2dFOP(pg.core.ModellingBase):
    """FDEM 2d-LCI modelling class based on BlockMatrices."""

    def __init__(self, data, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers."""
        super().__init__(verbose)
        self.nlay = nlay
        self.FOP = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        npar = 2 * nlay - 1
        self.mesh2d = pg.Mesh()
        self.mesh2d.create2DGrid(range(npar+1), range(self.nx+1))
        self.setMesh(self.mesh2d)

        self.J = pg.matrix.BlockMatrix()
        self.FOP1d = []
        for i in range(self.nx):
            self.FOP1d.append(HEMmodelling(nlay, data.z[i]))
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, self.nf*2*i, npar*i)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())
        self.setJacobian(self.J)

    def response(self, model):
        """Response as pasted forward responses from all soundings."""
        modA = np.reshape(model, (self.nx, self.nlay*2-1))
        resp = pg.Vector(0)
        for i, modi in enumerate(modA):
            resp = pg.cat(resp, self.FOP1d[i].response(modi))

        return resp

    def createJacobian(self, model):
        """Fill Jacobian (block) matrix by computing each block."""
        modA = np.reshape(model, (self.nx, self.nlay*2-1))
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class FDEMSmoothModelling(MeshModelling):
    """Occam-style (smooth) inversion."""

    def __init__(self, thk, **kwargs):
        super().__init__()
        self.thk_ = thk
        self.nlay_ = len(thk)+1
        self.core = HEMmodelling(**kwargs, nLayers=self.nlay_)
        self.mesh_ = pg.meshtools.createMesh1D(self.nlay_)
        self.setMesh(self.mesh_)

    def response(self, par):
        """Model response (forward modelling)."""
        return self.core.response(pg.cat(self.thk_, par))


if __name__ == '__main__':
    numlay = 3
    elevation = 30.0
    resistivity = np.array([1000.0, 100.0, 1000.0], float)
    thickness = np.array([22.0, 29.0], float)
    fop = HEMmodelling(numlay, elevation, r=10)  # frequency, separation)
    IP, OP = fop.vmd_hem(elevation, resistivity, thickness)
    pg.info(IP, OP)
