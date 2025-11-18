<!---
Readme for Github repository only. (Gets selected before *.rst file)
-->

<a href="https://www.pygimli.org">
  <img src="https://www.pygimli.org/_images/pg_logo.png" width="50%">
</a>

[![pyGIMLi CI](https://github.com/gimli-org/gimli/actions/workflows/main.yml/badge.svg)](https://github.com/gimli-org/gimli/actions/workflows/main.yml)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/license.svg)](https://pygimli.org/license.html)
[![release](https://img.shields.io/github/release/gimli-org/gimli.svg?style=flat-square)](https://github.com/gimli-org/gimli/releases/latest)
[![Github commits (since latest release)](https://img.shields.io/github/commits-since/gimli-org/gimli/latest.svg?style=flat-square)](https://github.com/gimli-org/gimli/tree/dev)
[![Slack](https://img.shields.io/badge/pyGIMLi%20chat%20-%20mattermost?style=flat&logo=mattermost&label=mattermost&link=https%3A%2F%2Fmattermost.softwareunderground.org%2Fswung%2Fchannels%2Fpygimli
)](https://mattermost.softwareunderground.org/swung/channels/pygimli)

[pyGIMLi](https://www.pygimli.org) is an open-source library for modelling and inversion and in geophysics. The object-oriented library provides management for structured and unstructured meshes in 2D and 3D, finite-element and finite-volume solvers, various geophysical forward operators, as well as Gauss-Newton based frameworks for constrained, joint and fully-coupled inversions with flexible regularization.

What is pyGIMLi suited for?

- analyze, visualize and invert geophysical data in a reproducible manner
- forward modelling of (geo)physical problems on complex 2D and 3D geometries
- inversion with flexible controls on a-priori information and regularization
- combination of different methods in constrained, joint and fully-coupled inversions
- teaching applied geophysics (e.g. in combination with [Jupyter notebooks])

What is pyGIMLi **NOT** suited for?

-   for people that expect a ready-made GUI for interpreting their data

[jupyter notebooks]: https://jupyter.org

##### Installation

[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/platforms.svg)](https://anaconda.org/gimli/pygimli)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/downloads.svg)](https://anaconda.org/gimli/pygimli)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/version.svg)](https://anaconda.org/gimli/pygimli)
[![Anaconda-Server Badge](https://anaconda.org/gimli/pygimli/badges/latest_release_date.svg)](https://anaconda.org/gimli/pygimli)

On all platforms, we recommend to install pyGIMLi via the conda package manager
contained in the Anaconda distribution. For details on how to install Anaconda,
we refer to: https://docs.anaconda.com/anaconda/install/

Note that Anaconda comes with many (great) packages, many of which you likely
will not use. If you want to save space, you can install the [light-weight
version Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/).

To avoid conflicts with other packages, we recommend to install pyGIMLi in a
separate environment. Here we call this environment pg, but you can give it any
name. Note that this environment has to be created only once.

``` bash
conda create -n pg -c gimli -c conda-forge "pygimli>=1.5.0"
```

If you are using Windows or Mac, a new environment named “pg” should be visible in the Anaconda Navigator. If you want to use pygimli from the command line, you have to activate the environment. You can put this line in your ~/.bashrc file so that it is activated automatically if you open a terminal.

``` bash
conda activate pg
```

See https://www.pygimli.org/installation.html for more information.

##### Import convention

```python
import pygimli as pg
print(pg.__version__)
```

Check www.pygimli.org for additional information, detailed installation
instructions and many examples.

#### Citing pyGIMLi

More information can be found in [this paper]. If you use pyGIMLi for your work, please cite as:

> Rücker, C., Günther, T., Wagner, F.M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences, 109, 106-123, doi: 10.1016/j.cageo.2017.07.011.

[this paper]: http://www.sciencedirect.com/science/article/pii/S0098300417300584/pdfft?md5=44253eaacd5490e3fb32210671672496&pid=1-s2.0-S0098300417300584-main.pdf

BibTeX code:

```sourceCode
@article{Ruecker2017,
  title = {{pyGIMLi}: An open-source library for modelling and inversion in geophysics},
  journal = {Computers and Geosciences},
  volume = {109},
  pages = {106--123},
  year = {2017},
  doi = {10.1016/j.cageo.2017.07.011},
  url = {https://www.sciencedirect.com/science/article/pii/S0098300417300584},
  author = {R\"ucker, C. and G\"unther, T. and Wagner, F. M.}
}
```

##### License

pyGIMLi is distributed under the terms of the **Apache 2.0** license. Details on
the license agreement can be found [here].

[here]: https://www.pygimli.org/license.html

#### Credits

We use or link some third-party software (beside the usual tool stack: cmake, gcc, boost, python, numpy, scipy, matplotlib) and are grateful for all the work made by the authors of these awesome open-source tools:

* libkdtree++: Maybe abandoned, mirror: https://github.com/nvmd/libkdtree

* meshio: https://github.com/nschloe/meshio

* pyplusplus: https://pypi.org/project/pyplusplus/

* pyvista: https://docs.pyvista.org/

* suitesparse, umfpack: https://people.engr.tamu.edu/davis/suitesparse.html

* Tetgen: http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1

* Triangle: https://www.cs.cmu.edu/~quake/triangle.html
