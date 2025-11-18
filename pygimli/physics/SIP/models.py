#!/usr/bin/env python
"""Spectral induced polarization (SIP) relaxations models."""

from math import pi
import numpy as np
import pygimli as pg


def ColeColeRho(f, rho, m, tau, c, a=1):
    """Backward compatibility."""
    pg.deprecated("Please use modelColeColeRho instead of ColeColeRho.")
    return modelColeColeRho(f, rho, m, tau, c, a)


def modelColeColeRho(f, rho, m, tau, c, a=1):
    r"""Frequency-domain Cole-Cole impedance model after Pelton et al. (1978).

    Frequency-domain Cole-Cole impedance model after Pelton et al. (1978)
    :cite:`PeltonWarHal1978`

    .. math::

        Z(\omega) & = \rho_0\left[1 - m \left(1 -
        \frac{1}{1+(\text{i}\omega\tau)^c}\right)\right] \\
        \quad \text{with}\quad m & = \frac{1}{1+\frac{\rho_0}{\rho_1}} \quad
        \text{and}\quad \omega =2\pi f

    * :math:`Z(\omega)` - Complex impedance per 1A current injection
    * :math:`f` - Frequency
    * :math:`\rho_0` -- Background resistivity states the unblocked pore path
    * :math:`\rho_1` -- Resistance of the solution in the blocked pore passages
    * :math:`m` -- Chargeability after Seigel (1959) :cite:`Seigel1959` as
      being the ratio of voltage immediately after, to the voltage immediately
      before cessation of an infinitely long charging current.
    * :math:`\tau` --
      'Time constant' relaxation time [s] for 1/e decay
    * :math:`c` - Rate of charge accumulation.
      Cole-Cole exponent typically [0.1 .. 0.6]

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.physics.SIP import modelColeColeRho
    >>> f = np.logspace(-2, 5, 100)
    >>> m = np.linspace(0.1, 0.9, 5)
    >>> tau = 0.01
    >>> fImMin = 1/(tau*2*np.pi)
    >>> fig, axs = pg.plt.subplots(1, 2)
    >>> ax1 = axs[0]
    >>> ax2 = axs[0].twinx()
    >>> ax3 = axs[1]
    >>> ax4 = axs[1].twinx()
    >>> for i in range(len(m)):
    ...     Z = modelColeColeRho(f, rho=1, m=m[i], tau=tau, c=0.5)
    ...     _= ax1.loglog(f, np.abs(Z), color='black')
    ...     _= ax2.loglog(f, -np.angle(Z)*1000, color='b')
    ...     _= ax3.loglog(f, Z.real, color='g')
    ...     _= ax4.semilogx(f, Z.imag, color='r')
    ...     _= ax4.plot([fImMin, fImMin], [-0.2, 0.], color='r')
    >>> _ = ax4.text(fImMin, -0.1, r"$f($min($Z''$))=$\frac{1}{2*\pi\tau}$", color='r')
    >>> _ = ax4.text(0.1, -0.17, r"$f($min[$Z''$])=$\frac{1}{2\pi\tau}$", color='r')
    >>> _ = ax1.set_ylabel('Amplitude $|Z(f)|$', color='black')
    >>> _ = ax1.set_xlabel('Frequency $f$ [Hz]')
    >>> _ = ax1.set_ylim(1e-2, 1)
    >>> _ = ax2.set_ylabel(r'- Phase $\varphi$ [mrad]', color='b')
    >>> _ = ax2.set_ylim(1, 1e3)
    >>> _ = ax3.set_ylabel('re $Z(f)$', color='g')
    >>> _ = ax4.set_ylabel('im $Z(f)$', color='r')
    >>> _ = ax3.set_xlabel('Frequency $f$ [Hz]')
    >>> _ = ax3.set_ylim(1e-2, 1)
    >>> _ = ax4.set_ylim(-0.2, 0)
    >>> pg.plt.show()
    """
    z = (1. - m * (1. - relaxationTerm(f, tau, c, a))) * rho
    if np.isnan(z).any():
        print(f, 'rho', rho, 'm', m, 'tau', tau, 'c', c)
        pg.critical(z)
    return z


def ColeColeRhoDouble(f, rho, m1, t1, c1, m2, t2, c2):
    """Backward compatibility."""
    pg.deprecated("Use modelColeColeRhoDouble instead of ColeColeRhoDouble.")
    return modelColeColeRhoDouble(f, rho, m1, t1, c1, m2, t2, c2)


def modelColeColeRhoDouble(f, rho, m1, t1, c1, m2, t2, c2, a=1, mult=False):
    """Frequency-domain double Cole-Cole resistivity (impedance) model.

    Frequency-domain Double Cole-Cole impedance model returns the sum of
    two Cole-Cole Models with a common amplitude.
    Z = rho * (Z1(Cole-Cole) + Z2(Cole-Cole))
    """
    Z1 = modelColeColeRho(f, rho=1, m=m1, tau=t1, c=c1, a=a)
    Z2 = modelColeColeRho(f, rho=1, m=m2, tau=t2, c=c2, a=a)
    if mult:
        return rho * Z1 * Z2
    else:
        return rho * (Z1 + Z2 - 1)


def ColeColeSigma(f, sigma, m, tau, c, a=1):
    """Backward compatibility."""
    pg.deprecated("Please use modelColeColeSigma instead of ColeColeSigma.")
    return modelColeColeSigma(f, sigma, m, tau, c, a)


def modelColeColeSigma(f, sigma, m, tau, c, a=1):
    """Complex-valued conductivity (admittance) Cole-Cole model."""
    return (1. + m / (1-m) * (1. - relaxationTerm(f, tau, c, a))) * sigma


def modelColeColeSigmaTauRho(f, sigma, m, tau, c, a=1):
    """Complex-valued conductivity (admittance) Cole-Cole model."""
    return (1. + m / (1-m) * (1. - relaxationTerm(f, tau, c, a)*(1-m))) * sigma


def modelColeColeSigmaDouble(f, sigma, m1, t1, c1, m2, t2, c2, a=1,
                             tauRho=True):
    """Complex-valued double added conductivity (admittance) model."""
    if tauRho:
        R1 = 1. - relaxationTerm(f, tau=t1, c=c1, a=a, p=1-m1)
        R2 = 1. - relaxationTerm(f, tau=t2, c=c2, a=a, p=1-m2)
        return sigma * (1. + m1 / (1 - m1) * R1 + m2 / (1 - m2) * R2)
    else:
        A1 = modelColeColeSigma(f, sigma=1, m=m1, tau=t1, c=c1, a=a)
        A2 = modelColeColeSigma(f, sigma=1, m=m2, tau=t2, c=c2, a=a)
        return (A1 + A2 - 1.) * sigma


def tauRhoToTauSigma(tRho, m, c):
    r"""Convert :math:`\tau_{\rho}` to :math:`\tau_{\sigma}` Cole-Cole model.

    .. math::

        \tau_{\sigma} = \tau_{\rho}/(1-m)^{\frac{1}{c}}

    Examples
    --------
    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.physics.SIP import modelColeColeRho, modelColeColeSigma
    >>> from pygimli.physics.SIP import tauRhoToTauSigma
    >>> tr = 1.
    >>> Z = modelColeColeRho(1e5, rho=10.0, m=0.5, tau=tr, c=0.5)
    >>> ts = tauRhoToTauSigma(tr, m=0.5, c=0.5)
    >>> S = modelColeColeSigma(1e5, sigma=0.1, m=0.5, tau=ts, c=0.5)
    >>> abs(1.0/S - Z) < 1e-12
    True
    >>> float(np.angle(1.0/S / Z)) < 1e-12
    True
    """
    return tRho * (1-m) ** (1/c)


def relaxationTerm(f, tau, c=1., a=1., p=1.):
    r"""Auxiliary function for Debye type relaxation term of the basic type.

    .. math::

        A(f,\tau,c)  = \frac{1}{1+(\jmath 2\pi f \tau)^c}

    or the more generalized term

    .. math::

        A(f,\tau,c,a,p)  = \frac{1}{(1+(\jmath 2\pi f \tau)^c p)^a}

    With c=1 and a one yields the Cole-Davidson model.
    With p=(1-m) one can account for the difference of sigma and rho tau.
    """
    return 1. / ((f * 2. * pi * tau * 1j)**c * p + 1.)**a


def DebyeRelaxation(f, tau, m):
    """Backward compatibility."""
    pg.deprecated("Please use relaxationDebye instead of DebyeRelaxation.")
    return relaxationDebye(f, tau, m)


def relaxationDebye(f, tau, m):
    """Complex-valued single Debye relaxation term with chargeability."""
    return 1. - (1. - relaxationTerm(f, tau)) * m


def WarbugRelaxation(f, tau, m):
    """Backward compatibility."""
    pg.deprecated("Please use relaxationWarbug instead of WarbugRelaxation.")
    return relaxationWarbug(f, tau, m)


def relaxationWarbug(f, tau, m):
    """Complex-valued single Debye relaxation term with chargeability."""
    return 1. - (1. - relaxationTerm(f, tau, c=0.5)) * m


def ColeColeEpsilon(f, e0, eInf, tau, alpha):
    """Backward compatibility."""
    pg.deprecated("Use modelColeColeEpsilon instead of ColeColeEpsilon.")
    return modelColeColeEpsilon(f, e0, eInf, tau, alpha)


def modelColeColeEpsilon(f, e0, eInf, tau, alpha):
    """Original complex-valued permittivity formulation (Cole&Cole, 1941)."""
    return (e0 - eInf) * relaxationTerm(f, tau, c=1./alpha) + eInf


def ColeCole(f, R, m, tau, c, a=1):
    """Backward compatibility."""
    pg.deprecated("Please use modelColeColeRho instead of ColeCole.")
    return modelColeColeRho(f, R, m, tau, c, a)


def ColeDavidson(f, R, m, tau, a=1):
    """Backward compatibility."""
    pg.deprecated("Please use modelColeDavidson instead of ColeDavidson.")
    return modelColeDavidson(f, R, m, tau, a)


def modelColeDavidson(f, R, m, tau, a=1):
    """For backward compatibility."""
    return ColeCole(f, R, m, tau, c=1, a=a)


class ColeColePhi(pg.core.ModellingBase):
    r"""Cole-Cole model with EM term after Pelton et al. (1978).

    Modelling operator for the Frequency Domain
    :py:mod:`Cole-Cole <pygimli.physics.SIP.modelColeColeRho>` impedance model
    using :py:mod:`pygimli.physics.SIP.modelColeColeRho` after
    Pelton et al. (1978) :cite:`PeltonWarHal1978`

    * :math:`\textbf{m} =\{ m, \tau, c\}`

        Modelling parameter for the Cole-Cole model with :math:`\rho_0 = 1`

    * :math:`\textbf{d} =\{\varphi_i(f_i)\}`

        Modelling response for all given frequencies as negative phase angles
        :math:`\varphi(f) = -tan^{-1}\frac{\text{Im}\,Z(f)}{\text{Re}\,Z(f)}`
        and :math:`Z(f, \rho_0=1, m, \tau, c) =` Cole-Cole impedance.
    """

    def __init__(self, f, verbose=False):
        """Set up class by specifying the frequency."""
        super().__init__(self, verbose)
        self.f_ = f
        self.setMesh(pg.meshtools.createMesh1D(1, 3))

    def response(self, par):
        """Phase angle of the model."""
        spec = modelColeColeRho(self.f_, 1.0, par[0], par[1], par[2])
        return -np.angle(spec)


class DoubleColeCole(pg.Modelling):
    r"""Complex double Cole-Cole model with EM term after Pelton et al. (1978).

    Modelling operator for the Frequency Domain
    :py:mod:`Cole-Cole <pygimli.physics.SIP.modelColeColeRho>` impedance model
    using :py:mod:`pygimli.physics.SIP.modelColeColeRho` after
    Pelton et al. (1978) :cite:`PeltonWarHal1978`

    * :math:`\textbf{m} =\{ m_1, \tau_1, c_1, m_2, \tau_2, c_2\}`

        Modelling parameter for the Cole-Cole model with :math:`\rho_0 = 1`

    * :math:`\textbf{d} =\{\varphi_i(f_i)\}`

        Modelling Response for all given frequencies as negative phase angles
        :math:`\varphi(f) = \varphi_1(Z_1(f))+\varphi_2(Z_2(f)) =
        -tan^{-1}\frac{\text{Im}\,(Z_1Z_2)}{\text{Re}\,(Z_1Z_2)}`
        and :math:`Z_1(f, \rho_0=1, m_1, \tau_1, c_1)` and
        :math:`Z_2(f, \rho_0=1, m_2, \tau_2, c_2)` ColeCole impedances.

    """

    def __init__(self, f, rho=True, tauRho=True, mult=False, aphi=True,
                 verbose=False):
        """Set up class by specifying the frequency."""
        super().__init__(verbose=verbose)
        self.f_ = f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 7))  # 4 single parameters
        self.rho = rho
        self.tauRho = tauRho
        self.mult = mult
        self.aphi = aphi

    def response(self, par):
        """Amplitude/phase or real/imag angle of the model."""
        if self.rho:
            out = modelColeColeRhoDouble(self.f_, *par, a=1, mult=self.mult)
        else:
            out = modelColeColeSigmaDouble(self.f_, *par, a=1,
                                           tauRho=self.tauRho)

        if self.aphi:
            return np.hstack((np.abs(out), np.abs(np.angle(out))))
        else:
            return np.hstack((out.real, np.abs(out.imag)))


class DoubleColeColePhi(pg.core.ModellingBase):
    r"""Double Cole-Cole model after Pelton et al. (1978).

    Modelling operator for the Frequency Domain - phase only
    :py:mod:`Cole-Cole <pygimli.physics.SIP.modelColeColeRho>` impedance model
    using :py:mod:`pygimli.physics.SIP.modelColeColeRho` after
    Pelton et al. (1978) :cite:`PeltonWarHal1978`

    * :math:`\textbf{m} =\{ m_1, \tau_1, c_1, m_2, \tau_2, c_2\}`

        Modelling parameter for the Cole-Cole model with :math:`\rho_0 = 1`

    * :math:`\textbf{d} =\{\varphi_i(f_i)\}`

        Modelling Response for all given frequencies as negative phase angles
        :math:`\varphi(f) = \varphi_1(Z_1(f))+\varphi_2(Z_2(f)) =
        -tan^{-1}\frac{\text{Im}\,(Z_1Z_2)}{\text{Re}\,(Z_1Z_2)}`
        and :math:`Z_1(f, \rho_0=1, m_1, \tau_1, c_1)` and
        :math:`Z_2(f, \rho_0=1, m_2, \tau_2, c_2)` ColeCole impedances.

    """

    def __init__(self, f, verbose=False):
        """Set up class by specifying the frequency."""
        super().__init__(verbose)
        self.f_ = f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 6))  # 4 single parameters

    def response(self, par):
        """Phase angle of the model."""
        spec1 = modelColeColeRho(self.f_, 1.0, par[0], par[1], par[2])
        spec2 = modelColeColeRho(self.f_, 1.0, par[3], par[4], par[5])

        return -np.angle(spec1) - np.angle(spec2)


class ColeColeAbs(pg.core.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)."""

    def __init__(self, f, verbose=False):
        super().__init__(verbose)
        self.f_ = f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 4))  # 3 single parameters

    def response(self, par):
        """Amplitude of the model."""
        spec = modelColeColeRho(self.f_, par[0], par[1], par[2], par[3])
        return np.abs(spec)


class ColeColeComplex(pg.core.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)."""

    def __init__(self, f, verbose=False):
        """Set up instance with frequency."""
        super().__init__(verbose)
        self.f_ = f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 4))  # 4 single parameters

    def response(self, par):
        """Phase angle of the model."""
        spec = modelColeColeRho(self.f_, *par)
        return pg.cat(np.abs(spec), -np.angle(spec))


class ColeColeComplexSigma(pg.core.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)."""

    def __init__(self, f, verbose=False):
        """Set up instance with frequency."""
        super().__init__(verbose)
        self.f_ = f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 4))  # 4 single parameters

    def response(self, par):
        """Phase angle of the model."""
        spec = modelColeColeSigma(self.f_, *par)
        return pg.cat(np.real(spec), np.imag(spec))


class PeltonPhiEM(pg.core.ModellingBase):
    """Cole-Cole model with EM term after Pelton et al. (1978)."""

    def __init__(self, f, verbose=False):
        """Set up instance with frequency."""
        super().__init__(verbose)
        self.f_ = f  # save frequencies
        self.setMesh(pg.meshtools.createMesh1D(1, 4))  # 4 single parameters

    def response(self, par):
        """Phase angle of the model."""
        spec = modelColeColeRho(self.f_, 1.0, par[0], par[1], par[2]) * \
            relaxationTerm(self.f_, par[3])  # pure EM has c=1
        return -np.angle(spec)


class DebyePhi(pg.Modelling):
    """Debye decomposition (smooth Debye relaxations) phase only."""

    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        """Init with frequency and tau vector."""
        super().__init__(verbose=verbose)
        self.f_ = fvec
        self.nf_ = len(fvec)
        self.t_ = tvec
        self.mesh = pg.meshtools.createMesh1D(len(tvec))  # standard 1d mesh
        self.setMesh(self.mesh)

    def response(self, par):
        """amplitude/phase spectra as function of spectral chargeabilities."""
        y = np.ones(self.nf_, dtype=np.complex)  # 1 -
        for (tau, mk) in zip(self.t_, par, strict=False):
            y -= (1. - relaxationTerm(self.f_, tau)) * mk

        return -np.angle(y)


class DebyeComplex(pg.Modelling):
    """Debye decomposition (smooth Debye relaxations) of complex data."""

    def __init__(self, fvec, tvec, verbose=False):  # save reference in class
        """Init with frequency and tau vector."""
        self.f = fvec
        self.nf = len(fvec)
        self.t = tvec
        self.nt = len(tvec)
        mesh = pg.meshtools.createMesh1D(len(tvec))  # standard 1d mesh
        super().__init__(mesh=mesh, verbose=verbose)
        T, W = np.meshgrid(tvec, fvec * 2. * pi)
        WT = W*T
        self.A = WT**2 / (WT**2 + 1)
        self.B = WT / (WT**2+1)
        self.J = pg.Matrix()
        self.J.resize(len(fvec)*2, len(tvec))
        for i in range(self.nf):
            wt = fvec[i] * 2.0 * pi * tvec
            wt2 = wt**2
            self.J[i] = wt2 / (wt2 + 1.0)
            self.J[i+self.nf] = wt / (wt2 + 1.0)

        self.setJacobian(self.J)

    def response(self, par):
        """amplitude/phase spectra as function of spectral chargeabilities."""
        return self.J * par

    def createJacobian(self, par):
        """Linear jacobian after Nordsiek&Weller (2008)."""
        pass


if __name__ == "__main__":
    pass
