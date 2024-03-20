Installation
============

.. raw:: html

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


On all platforms, we recommend install pyGIMLi via the conda package manager
and the cleanest way is to create a new environment.
Besides large installers like `Anaconda <https://anaconda.org/>`_
that come with many (great) packages, you can save space and install a 
light-weight installer like
`Miniconda <https://docs.anaconda.com/free/miniconda/miniconda-install/>`_
or `Miniforge <https://conda-forge.org/download/>`_.

.. note::

    After installing Anaconda or Miniconda, make sure that you have at least
    version 23.10.0 of conda (`conda --version`).

.. TODO: Temporarily disabled. PDF needs to update version automatically.
.. A **step-by-step guide for Windows users** can be found `here
.. <https://www.pygimli.org/_downloads/pygimli_win_install_guide.pdf>`_.

To avoid conflicts with other packages, we recommend to install pyGIMLi in a
separate environment. Here we call this environment `pg`, but you can give
it any name. Note that this environment has to be created only once.

Open a terminal (Linux & Mac) or the Anaconda Prompt (Windows) and type:

.. code-block:: bash

    conda create -n pg -c gimli -c conda-forge "pygimli>=1.5.0"

If you want to use pyGIMLi from the command line, you have to activate the
environment. You can put this line in your `~/.bashrc` file so that it is
activated automatically if you open a terminal.

.. code-block:: bash

    conda activate pg

After that you can use pyGIMLi with your text editor of choice and a terminal.

Usage with Spyder or JupyterLab
-------------------------------

Depending on your preferences, you can also install third-party software such as
the MATLAB-like integrated development environment (https://www.spyder-ide.org):

.. code-block:: bash

    conda install -c conda-forge spyder

Or alternatively, the web-based IDE JupyterLab (https://jupyterlab.readthedocs.io):

.. code-block:: bash

    conda install -c conda-forge jupyterlab

If you do one of the above steps in the `pg` environment, then it will
automatically find pyGIMLi. But you may not want to install JupyterLab or
Spyder for every different environment. To use your existing JupyterLab
installation in the `base` environment with pyGIMLi in the `pg` environment,
follow these steps:

.. code-block:: bash

    conda activate pg
    conda install ipykernel
    conda activate base
    conda install -c conda-forge nb_conda_kernels
    jupyter lab

pyGIMLi on Google Colab
-----------------------
pyGIMLi can be run on Google Colab without any installation on your own
computer. Just create a new Notebook and install pyGIMLi via pip:

.. code:: python

    !pip install pygimli tetgen

It turns out that there are some packages preinstalled that lead to some
incompatibl numpy version, so you might have to uninstall them first.

.. code:: python

    !pip uninstall -y numba tensorflow pytensor thinc
    !pip install pygimli tetgen

Staying up-to-date
------------------

Update your pyGIMLi installation from time to time, if want to have the newest
functionality. 

.. code-block:: bash

    conda update -c gimli -c conda-forge pygimli

If there something went wrong and you are running an old, not further
supported python version, consider a fresh install in a new clean environment.
The only drawback of using conda is that you are bound to the rhythm in which we
update the conda packages. In order to work with the latest Python codes you
should create an environment with the latest pyGIMLi C++ core only,

.. code-block:: bash

    conda create -n pgcore -c gimli -c conda-forge pgcore
    
retrieve the source code by git

.. code-block:: bash

    git clone https://github.com/gimli-org/gimli
    cd gimli

and install pygimli as a development package using conda

.. code-block:: bash

    conda develop .

or using pip

.. code-block:: bash

    pip install --no-build-isolation --no-deps -e .

Alternatively you could set the PYTHONPATH variable but you would have to care
for dependencies by yourself.

Later you can just update the pygimli code by

.. code-block:: bash

    git pull
    
Only if you need recent changes to the C++ core, you have to compile
pyGIMLi using your systems toolchain as described in 
https://www.pygimli.org/compilation.html#sec-build
