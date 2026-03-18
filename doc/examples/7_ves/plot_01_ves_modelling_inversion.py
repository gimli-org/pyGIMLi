#!/usr/bin/env python
# encoding: utf-8
r"""
VES modelling and inversion
===========================

Vertical electrical sounding (VES) is a classical 1D geophysical technique to
determine layer resistivities and thicknesses from surface measurements of
apparent resistivity as a function of electrode spacing.

In a Schlumberger array, the current electrodes A and B are moved symmetrically
around a fixed centre while the potential electrodes M and N remain close to
the centre. The measured quantity is the apparent resistivity :math:`\rho_a`
as a function of the AB/2 half-spacing.

We use the :py:class:`pygimli.physics.ves.VESManager` to create synthetic
data, add noise, and invert for a layered subsurface model.
"""

# sphinx_gallery_thumbnail_number = 3
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.physics import VESManager

# %%%
# Model setup
# -----------
# We define a four-layer earth model with resistivities in Ohm·m and layer
# thicknesses in metres. The deepest layer (half-space) has no thickness.
#

nlay = 4                             # number of layers
synthRes = [100., 500., 20., 800.]   # resistivities in Ohm·m
synthThk = [0.5, 3.5, 6.]           # layer thicknesses in metres

# %%%
# Electrode geometry
# ------------------
# AB/2 spacings are chosen on a logarithmic scale covering two decades.
# The MN/2 spacing is one third of AB/2, which is common for Schlumberger
# arrays.
#

ab2 = np.logspace(np.log10(1.5), np.log10(100.), 25)   # AB/2 in metres
mn2 = ab2 / 3.                                          # MN/2 in metres

# %%%
# Forward modelling and data generation
# --------------------------------------
# We instantiate the :py:class:`VESManager` and generate synthetic apparent
# resistivity data from the known model, then add 3% Gaussian noise to
# simulate field measurements.
#

ves = VESManager()
noiseLevel = 0.03   # 3 % relative noise
synthModel = pg.cat(synthThk, synthRes)   # thicknesses first, then resistivities
ra, err = ves.simulate(synthModel, ab2=ab2, mn2=mn2, noiseLevel=noiseLevel,
                       seed=42)

# %%%
# Display the synthetic sounding curve (apparent resistivity vs. AB/2).
# The data are plotted on a log–log scale with the depth axis pointing
# downward, which is the standard convention for VES data.
#

fig, ax = plt.subplots()
ax.errorbar(ab2, ra, yerr=err * ra, fmt='x-', label='Synthetic data')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('AB/2 (m)')
ax.set_ylabel(r'Apparent resistivity ($\Omega\cdot$m)')
ax.set_title('Synthetic VES sounding curve')
ax.legend()
ax.grid(True, which='both', ls='--', alpha=0.5)

# %%%
# Inversion
# ---------
# We invert the noisy synthetic data using the
# :py:class:`VESManager`.  The manager automatically sets up the forward
# operator, the log-transformed model space, and the Gauss–Newton
# inversion loop.  ``nLayers`` specifies the number of layers to invert for.
#

model = ves.invert(ra, err, ab2=ab2, mn2=mn2, nLayers=nlay, verbose=True)
print("Inversion chi² =", round(ves.inv.chi2(), 3))

# %%%
# Data fit
# --------
# We compare the measured (noisy) data with the model response to assess
# the quality of the fit.
#

fig, ax = plt.subplots()
ax.errorbar(ab2, ra, yerr=err * ra, fmt='x', label='Data', color='C0')
ax.loglog(ab2, ves.inv.response, '-', label='Model response', color='C1')
ax.set_xlabel('AB/2 (m)')
ax.set_ylabel(r'Apparent resistivity ($\Omega\cdot$m)')
ax.set_title('VES data fit')
ax.legend()
ax.grid(True, which='both', ls='--', alpha=0.5)

# %%%
# Model comparison
# ----------------
# Finally, we show the synthetic (true) model and the inverted model
# side by side.  The :py:meth:`VESManager.showResult` method draws the
# inverted model, and we overlay the true model for comparison.
#

ax, _ = ves.showModel(synthModel, label='True model')
ves.showResult(ax=ax, label='Inverted model')
ax.set_title('VES 1D model comparison')

pg.wait()
