# Software built upon pyGIMLi

This section highlights research software and tools that are built upon pyGIMLi. These projects demonstrate how the library or specific components of it are used in more domain-specific or user-friendly GUI-based software packages.

## custEM

An open-source Python toolkit for customizable 3D finite-element modeling of controlled-source electromagnetic (CSEM) data. Built on FEniCS, it integrates with pyGIMLi for 3D inversion of CSEM data and supports marine, land-based, airborne, and mixed EM scenarios.

> Rochlitz, R., Skibbe, N. & Günther, T. (2019). custEM: Customizable finite-element simulation of complex controlled-source electromagnetic data. *Geophysics*, 84(2), F17–F33. [DOI: 10.1190/geo2018-0208.1](https://doi.org/10.1190/geo2018-0208.1)

## formikoj

A flexible open-source library for managing and processing geophysical data in environmental and engineering investigations. Provides seismic waveform modeling and refraction processing capabilities, combining modeling and processing tasks in a single workflow with pyGIMLi integration for inversion.

> Steiner, M. & Flores Orozco, A. (2023). formikoj: A flexible library for data management and processing in geophysics—Application for seismic refraction data. *Computers & Geosciences*, 176, 105339. [DOI: 10.1016/j.cageo.2023.105339](https://doi.org/10.1016/j.cageo.2023.105339)

## four-phase-inversion

An open-source research code for petrophysical joint inversion of seismic refraction and electrical resistivity data to quantitatively image water, ice, air, and rock contents in permafrost systems. The inversion framework is entirely built on pyGIMLi.

> Wagner, F.M., Mollaret, C., Günther, T., Kemna, A. & Hauck, C. (2019). Quantitative imaging of water, ice, and air in permafrost systems through petrophysical joint inversion of seismic refraction and electrical resistivity data. *Geophysical Journal International*, 219(3), 1866–1875. [DOI: 10.1093/gji/ggz402](https://doi.org/10.1093/gji/ggz402)

## MC_RB_EM_1D

Software for estimation of electrical conductivity models using multi-coil rigid-boom electromagnetic induction measurements. Evaluates the well-posedness of inverse problems for layered conductivity models.

> Carrizo Mascarell, M., Werthmüller, D. & Slob, E. (2024). MC_RB_EM_1D: Estimation of electrical conductivity models using multi-coil rigid-boom electromagnetic induction measurements. *Computers & Geosciences*, 193, 105732. [DOI: 10.1016/j.cageo.2024.105732](https://doi.org/10.1016/j.cageo.2024.105732)

## PyMERRY

A post-processing tool for electrical resistivity tomography (ERT) that computes a data-coverage mask and evaluates resistivity uncertainties to support the reliable interpretation of inverted ERT images. Designed to work with 2D models produced by pyGIMLi and other ERT codes.

> Gautier, M., Gautier, S. & Cattin, R. (2023). PyMERRY: A Python solution for an improved interpretation of electrical resistivity tomography images. *Geophysics*, 89(1), F23–F39. [DOI: 10.1190/geo2023-0105.1](https://doi.org/10.1190/geo2023-0105.1)

## PyRefra

Open-source Python software for display, processing, and picking of near-surface refraction seismic data with tomographic modeling. Integrates with pyGIMLi for travel-time tomographic inversion.

> Zeyen, H. & Léger, E. (2024). PyRefra: Refraction seismic data treatment and inversion. *Computers & Geosciences*, 185, 105556. [DOI: 10.1016/j.cageo.2024.105556](https://doi.org/10.1016/j.cageo.2024.105556)

## Refrapy

An open-source Python package for seismic refraction data analysis. Provides basic waveform processing, first break picking, and inversion through time-terms analysis or travel-time tomography via GUI interaction.

> Guedes, V.J.C.B., Maciel, S.T.R. & Rocha, M.P. (2022). Refrapy: A Python program for seismic refraction data analysis. *Computers & Geosciences*, 159, 105020. [DOI: 10.1016/j.cageo.2021.105020](https://doi.org/10.1016/j.cageo.2021.105020)

## SardineReborn

A Python/PyQt5 GUI application for seismic refraction data display, picking, forward modeling, and travel-time tomography inversion via pyGIMLi. Oriented toward teaching and exploratory workflows; reads standard seismic formats (SEG2, SEGY) through ObsPy.

> Michel, H. (2025). Sardine Reborn. Zenodo. [DOI: 10.5281/zenodo.15089167](https://doi.org/10.5281/zenodo.15089167)

---

```{note}
Miss your software in this list? Send the reference to [mail@pygimli.org](mailto:mail@pygimli.org) or [add it directly here](https://github.com/gimli-org/gimli/edit/dev/doc/community/other-software.md).
```
