Electromagnetics
================

This module contains different solutions for electromagnetic measurements over laterally layered (VTI) media in either the frequency domain (FD) or time domain (TD).

* FDEM sounding (ground or airborne) using VMD (or possible HMD) sources and receivers, e.g. RESOLVE airborne receiver or MaxMin ground system
* CSEM (FD) sounding with grounded bipole transmitter and magnetic (or optionally electric) field receivers
* TDEM sounding with loop responses (either impulse or step)

For 3D (or 2D) EM inversion, we refer to the packages custEM (Rochlitz et al., 2023) and SAEM.

References

Siemon, B. (2012): Accurate 1D forward and inverse modeling of high-frequency helicopter-borne electromagnetic data, Geophysics 77(4), WB71-87.

Werthmüller, D. (2017): An open-source full 3D electromagnetic modeler for 1D VTI media in Python: empymod, Geophysics 82(6), WB9-19.

Rochlitz, R., Becken, M. & Günther, T. (2023): Three-dimensional inversion of semi-airborne electromagnetic data with a second-order finite-element forward solver. Geophys. J. Int. 234(1), 528-545, doi:10.1093/gji/ggad056.
