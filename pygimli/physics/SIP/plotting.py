#!/usr/bin/env python
"""Spectral induced polarization (SIP) plotting tools."""

import pygimli as pg


def showAmplitudeSpectrum(*args, **kwargs):
    """Backward compatibility."""
    pg.deprecated('drawAmplitudeSpectrum')
    return drawAmplitudeSpectrum(*args, **kwargs)


def showPhaseSpectrum(*args, **kwargs):
    """Backward compatibility."""
    pg.deprecated('drawPhaseSpectrum')
    return drawPhaseSpectrum(*args, **kwargs)


def drawAmplitudeSpectrum(ax, freq, amp, ylabel=r'$\rho$ ($\Omega$m)',
                          grid=True, marker='+', ylog=True, **kwargs):
    """Show amplitude spectrum (resistivity as a function of f)."""
    if 'label' not in kwargs:
        kwargs['label'] = 'obs'

    gci = ax.semilogx(freq, amp, marker=marker, **kwargs)
    if ylog is None:
        ylog = (min(amp) > 0)
    if ylog:
        ax.set_yscale('log')
    # ax.set_ylim(min(amp) * .99, max(amp * 1.01))
    ax.set_xlabel('f (Hz)')
    ax.set_ylabel(ylabel)
    ax.grid(grid)
    ax.legend()
    return gci


def drawPhaseSpectrum(ax, freq, phi, ylabel=r'$-\phi$ (mrad)',
                      grid=True, marker='+', ylog=False, **kwargs):
    """Show phase spectrum (-phi as a function of f)."""
    if 'label' not in kwargs:
        kwargs['label'] = 'obs'

    gci = ax.semilogx(freq, phi, marker=marker, **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f (Hz)')
    ax.set_ylabel(ylabel)
    ax.grid(grid)
    ax.legend()
    return gci


def showSpectrum(freq, amp, phi, nrows=2, ylog=None, axs=None, **kwargs):
    """Show amplitude and phase spectra in two subplots."""
    if axs is None:
        fig, axs = pg.plt.subplots(nrows=nrows, sharex=(nrows == 2))
    else:
        fig = axs[0].figure
    drawAmplitudeSpectrum(axs[0], freq, amp, ylog=ylog, **kwargs)
    drawPhaseSpectrum(axs[1], freq, phi, ylog=ylog, **kwargs)
    return fig, axs


def plotSpectrum(ax, freq, vals, ylabel=r'$-\phi$ (mrad)',
                 grid=True, marker='+', ylog=True, **kwargs):
    """Plot some spectrum (redundant)."""
    pg.deprecated('drawSpectrum')
    if 'label' not in kwargs:
        kwargs['label'] = 'obs'
    ax.loglog(freq, vals, marker=marker, **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f (Hz)')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


if __name__ == "__main__":
    pass
