<a href="https://www.pygimli.org">
  <img src="https://www.pygimli.org/_images/pg_logo.png" width="50%">
</a>

[![Build Status](http://jenkins.pygimli.org/job/pyGIMLi_dev/badge/icon?style=flat-square)](http://jenkins.pygimli.org/job/pyGIMLi_dev/)
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

Before you start, considering its not a bad idea to use virtual environments, so give this a try:

``` bash
python -m venv pygimli
source pygimli/bin/activate
```

To install pygimli :

``` bash
python -m pip install pygimli
```

You might add the 'all' option to install also optional dependencies.

``` bash
python -m pip install pygimli['all']
```

You can see if the installation was successful:

``` bash
python -c 'import pygimli as pg; pg.version()'
```

For more information visit [pyGIMLi](https://www.pygimli.org).
