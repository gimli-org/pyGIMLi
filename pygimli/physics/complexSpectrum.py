"""Frequency sounding."""
import numpy as np
import pygimli as pg


class ComplexSpectrum(pg.DataContainer):
    """Class for storing frequency sounding data."""

    def __init__(self, **kwargs):
        """Init with f and either real/imag or amp/phi."""
        super().__init__()
        self["f"] = kwargs.pop("f", 0)
        self.unit = kwargs.pop("unit", "")
        self["real"] = kwargs.pop("real", 0)
        self["imag"] = kwargs.pop("imag", 0)
        if "amp" in kwargs and "phi" in kwargs:
            self.setAmpPhi(kwargs["amp"], kwargs["phi"])

    @property
    def f(self):
        """Frequency."""
        return self["f"]

    @f.setter
    def f(self, f):
        """Set frequency."""
        self["f"] = f

    @property
    def real(self):
        """Real part."""
        return self["real"].array()

    @real.setter
    def real(self, re):
        """Set real part."""
        self["real"] = re

    @property
    def imag(self):
        """Imaginary part."""
        return self["imag"].array()

    @imag.setter
    def imag(self, im):
        """Set imaginary part."""
        self["imag"] = im

    @property
    def amp(self):
        """Amplitude."""
        return np.sqrt(self["real"]**2+self["imag"]**2)

    @amp.setter
    def amp(self, a):
        """Set amplitude (and keep phase)."""
        phi = np.angle(self.real + self.imag*1j)
        self["real"] = np.cos(phi) * a
        self["imag"] = np.sin(phi) * a

    @property
    def phi(self):
        """Amplitude."""
        return np.angle(self.real+self.imag*1j)

    @phi.setter
    def phi(self, p):
        """Set amplitude (and keep phase)."""
        amp = np.sqrt(self["real"]**2+self["imag"]**2)
        self["real"] = np.cos(p) * amp
        self["imag"] = np.sin(p) * amp

    def show(self, ax=None, aphi=False, p1kw={}, p2kw={}, **kwargs):
        """Show sounding."""
        if aphi:
            p1, p2 = self.amp, self.phi
        else:
            p1, p2 = self.real, self.imag

        vertical = kwargs.pop("vertical", False)
        if ax is None:
            if vertical:
                fig, ax = pg.plt.subplots(ncols=2, sharey=True)
            else:
                fig, ax = pg.plt.subplots(nrows=2, sharex=True)

        if vertical:
            ax[0].plot(p1, self.f, **p1kw)
            ax[1].plot(p2, self.f, **p2kw)
            ax[0].set_yscale("log")
            ax[0].set_ylim(min(self.f), max(self.f))
        else:
            ax[0].plot(self.f, p1, **p1kw)
            ax[1].plot(self.f, p2, **p2kw)
            ax[0].set_xscale("log")
            ax[0].set_xlim(min(self.f), max(self.f))

        for a in ax:
            a.grid(True)

        return ax


# %%
if __name__ == "__main__":
    f = 2.**np.arange(14)
    self = ComplexSpectrum(f=f, real=1, imag=0.1+np.log10(f)/10)
    self.show(vertical=0)
    self.show(aphi=True, vertical=True)
