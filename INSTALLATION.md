(installation)=
# Installation

```{raw} html
<p style="height:22px">
  <a href="https://anaconda.org/gimli/pygimli" >
    <img src="https://anaconda.org/gimli/pygimli/badges/license.svg"/>
  </a>
  <a href="https://anaconda.org/gimli/pygimli" >
    <img src="https://anaconda.org/gimli/pygimli/badges/downloads.svg"/>
 <a href="https://anaconda.org/gimli/pygimli" >
    <img src="https://anaconda.org/gimli/pygimli/badges/version.svg?style=flat-square"/>
  </a>
 <a href="https://anaconda.org/gimli/pygimli" >
    <img src="https://anaconda.org/gimli/pygimli/badges/latest_release_date.svg?style=flat-square"/>
 </a>
 <a href="https://anaconda.org/gimli/pygimli" >
    <img src="https://anaconda.org/gimli/pygimli/badges/platforms.svg?style=flat-square"/>
 </a>
</p><br>
```

pyGIMLi is available for Windows, macOS, and Linux and can be installed via various package managers:

::::{tab-set}

:::{tab-item} pip

```bash
pip install pygimli
```

:::

:::{tab-item} uv

```bash
uv add pygimli
```

:::

:::{tab-item} conda

```bash
conda install -c gimli -c conda-forge pygimli
```

:::

:::{tab-item} pixi

```bash
pixi add --channel gimli --channel conda-forge pygimli
```

:::

::::

While we do not have a preference on the package manager you use, we recommend to avoid mixing them and generally encourage the use of separated environments, which will be discussed for `pip` and `conda` below.

## Installing via [`pip`](https://pip.pypa.io/en/stable/getting-started/)

`pip` is the default command to install Python packages.
On most systems `pip` ships together with Python, so installing a Python 3 is usually all you need:

::::{tab-set}

:::{tab-item} Windows
:sync: os-windows

Install Python with the official
[python.org installer](https://www.python.org/downloads/windows/) and tick
**"Add python.exe to PATH"** during setup, or install it from the Microsoft
Store.

:::

:::{tab-item} macOS
:sync: os-macos

Recent macOS versions no longer ship Python. Install it with the official
[python.org installer](https://www.python.org/downloads/macos/) or via
[Homebrew](https://brew.sh):

```bash
brew install python
```

:::

:::{tab-item} Linux
:sync: os-linux

Python is preinstalled on most distributions. Install Python, `pip` and `venv`
via your package manager, e.g.:

```bash
sudo apt install python3 python3-pip python3-venv    # Debian/Ubuntu
sudo dnf install python3 python3-pip                 # Fedora
sudo pacman -S python python-pip                     # Arch
```

:::

::::

Verify that everything is available (use `py` instead of `python3` on Windows):

```bash
python3 --version
python3 -m pip --version
```

To avoid conflicts with other packages we install pyGIMLi into a separate
environment, here called `pg` (you can pick any name; the environment only has to be created once). 

::::{tab-set}

:::{tab-item} Windows
:sync: os-windows

```bash
python -m venv .venv --prompt=pg
.venv/Scripts/activate
```

:::

:::{tab-item} macOS
:sync: os-macos

```bash
python -m venv .venv --prompt=pg
source .venv/bin/activate
```

:::


:::{tab-item} Linux
:sync: os-linux

```bash
python -m venv .venv --prompt=pg
source .venv/bin/activate
```

:::
::::

To install pygimli:

```bash
pip install pygimli
```

To update pygimli:

```bash
pip install -U pygimli
```

Find available versions:

```bash
pip index versions pygimli
```

To install a specific version:

```bash
pip install pygimli==x.x.x
```

## Installing via `conda/mamba`

We recommend to install the conda/mamba package managers via the lightweight and free Miniforge distribution (<https://conda-forge.org/download/>).

```bash
# These two lines are optional but recommended
conda create -n pg
conda activate pg

conda install -c gimli pygimli
```

Once the environment is activated you can use pyGIMLi from the command line with
your editor of choice. To activate it automatically in every new terminal, add
the activation command (`conda activate pg`) to
your shell startup file, e.g. `~/.bashrc`.

## Using pyGIMLi with Spyder or JupyterLab

Depending on your preferences, you can also install third-party software such as
the MATLAB-like integrated development environment
[Spyder](https://www.spyder-ide.org):

```bash
conda install -c conda-forge spyder
```

Or alternatively, the web-based IDE
[JupyterLab](https://jupyterlab.readthedocs.io):

```bash
conda install -c conda-forge jupyterlab
```

If you do one of the above steps in the `pg` environment, then it will
automatically find pyGIMLi. But you may not want to install JupyterLab or
Spyder for every different environment. To use your existing JupyterLab
installation in the `base` environment with pyGIMLi in the `pg` environment,
follow these steps:

```bash
conda activate pg
conda install ipykernel
conda activate base
conda install -c conda-forge nb_conda_kernels
jupyter lab
```

## pyGIMLi on Google Colab

Even though still experimental, pyGIMLi can be run on Google Colab without any
installation on your own computer. Just create a new notebook and install the
pyGIMLi package via pip:

```bash
!pip install pygimli tetgen
```

Some preinstalled packages may pull in an incompatible NumPy version, so you
might have to uninstall them first:

```python
!pip uninstall -y numba tensorflow pytensor thinc
!pip install pygimli tetgen
```

## Staying up-to-date

Update your pyGIMLi installation from time to time to get the newest
functionality:

::::{tab-set}
:::{tab-item} pip
:sync: pm-pip

```bash
pip install -U pygimli
```

:::

:::{tab-item} conda
:sync: pm-conda

```bash
conda update -c gimli -c conda-forge pygimli
```

:::
::::

If something went wrong and you are running an old, no longer supported Python
version, consider a fresh install in a new clean environment.

## Development version

The conda packages follow our release rhythm. To work with the latest Python
code from `git` while still using the pre-built C++ core, first create an
environment with only `pgcore`:


::::{tab-set}
:::{tab-item} pip
:sync: pm-pip

```bash
python -m venv .venv --prompt pgcore
source .venv/bin/activate # Mac/Linux
.venv/Scripts/activate # windows
pip install pgcore
```

:::

:::{tab-item} conda
:sync: pm-conda

```bash
conda create -n pgcore -c gimli -c conda-forge pgcore
conda activate pgcore
```

:::
::::

Retrieve the source code with git:

```bash
git clone https://github.com/gimli-org/gimli
cd gimli
```

and install pyGIMLi as an editable (development) package:

::::{tab-set}

:::{tab-item} pip
:sync: pm-pip

```bash
pip install --no-build-isolation --no-deps -e .
```

:::

:::{tab-item} conda
:sync: pm-conda

```bash
conda develop .
```

:::

::::

Alternatively you could set the `PYTHONPATH` variable, but then you would have to take care of dependencies yourself. Later you can update the pyGIMLi code with:

```bash
git pull
```

Only if you need recent changes to the C++ core itself do you have to compile pyGIMLi with your system toolchain, as described in the
[build guide](https://www.pygimli.org/compilation.html#sec-build).
