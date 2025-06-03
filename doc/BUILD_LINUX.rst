.. _sec:install_lin:

Building on Linux
-----------------

Compiling with a virtual environment
....................................

If you don't want to use the conda environment we encourage the use of a
virtual python environment.
Assuming you have a proper build toolchain and the required libraries
(see Installation on Ubuntu below) installed.

.. code-block:: bash

    mkdir -p gimli
    cd gimli
    git clone https://github.com/gimli-org/gimli.git
    cd gimli
    deactivate # in case there are other venv active
    python -m venv .venv --prompt=gimli
    . .venv/bin/activate
    python -m pip install -r dev_requirements.txt
    cd ..

We recommend an out of source build:

.. code-block:: bash

    mkdir -p build
    cd build
    cmake ../gimli
    make -j 4 gimli
    make pygimli J=4
    cd ..

There is no more need to change the `PATH` or `LD_LIBRARY_PATH`.
In fact, check to remove any prior changes to these environment variables from
older installations.

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$(pwd)/gimli
    python -c 'import pygimli as pg; pg.version()'
    python -c 'import pygimli as pg; print(pg.report())'


**All following build procedures are outdated and may not work without
tweaking**


Curl installer
..............
(Probably outdated.)

If you are not using Anaconda, you can build pyGIMLi from source in the current
directory via:

.. code-block:: bash

    curl -Ls install.pygimli.org | bash

This script accepts a few more options. For help see:

.. code-block:: bash

    curl -Ls install.pygimli.org | bash -s help

If something goes wrong please take a look at the error message. In most cases
there are missing or outdated packages. Please have a look at the prerequisites
tab or in the following Ubuntu section.

If the installation fails you can try the following instructions for manual
installation.


Example Installation on Ubuntu
..............................

Last try on Ubuntu 22.04.03 (23-11-14)

.. code-block:: bash

    sudo apt-get install build-essential g++ subversion git cmake \
                 python3-dev python3-matplotlib python3-numpy python3-pyqt5 \
                 python3-scipy libboost-all-dev libedit-dev \
                 libsuitesparse-dev libopenblas-openmp-dev libumfpack5 \
                 libcppunit-dev clang

    mkdir -p ~/src/gimli
    cd ~/src/gimli
    git clone https://github.com/gimli-org/gimli.git

    mkdir -p build
    cd build
    cmake ../gimli
    make -j 4 gimli
    make pygimli J=4

To use the gimli installation you need to set some environment variables (this
example assumes that the **src** directory resides in your home directory):

.. code-block:: bash

    GIMLI_INSTALLATION=$HOME/src/gimli
    export PYTHONPATH=$PYTHONPATH:$GIMLI_INSTALLATION/gimli
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GIMLI_INSTALLATION/build/lib
    export PATH=$PATH:$GIMLI_INSTALLATION/build/bin

You can do a quick test of the pygimli build and installation with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.__version__)'


Detailed Installation on Debian
...............................
(Probably outdated.)

Tested on Debian 12 (Bookworm).

In order to build pygimli (and gimli) Python 3, install
the required packages:

.. code-block:: bash

    sudo apt-get install wget subversion git cmake mercurial g++ \
        libboost-all-dev libblas-dev liblapack-dev libopenblas-dev \
        libsuitesparse-dev libedit-dev python3-dev \
        python3-numpy python3-matplotlib \
        python3-setuptools

Create a directory for your installation, e.g., $HOME/src/gimli

.. code-block:: bash

    mkdir -p ~/src/gimli
    cd ~/src/gimli

Checkout the current sources for libgimli:

.. code-block:: bash

    git clone https://github.com/gimli-org/gimli.git

We use `cmake <https://cmake.org>`_ for compilation. We recommend using a
build directory parallel to the gimli (trunk) path:

.. code-block:: bash

    mkdir -p build

The main directory structure should looks like this:

.. code-block:: bash

    gimli/gimli
    gimli/build

Change to the build path

.. code-block:: bash

    cd build

If you want to compile for Python 3.8, alternatively use:

.. code-block:: bash

    cmake ../gimli -DPYVERSION=3.8

If the output complains about missing dependencies, install these and repeat
the the last step. To build the library just run `make`.

.. code-block:: bash

    make

To speed up the build process using more CPUs, use the `-j` flag, e.g.:

.. code-block:: bash

    make -j 8

The libraries will be installed in **build/lib** and some test applications are
installed in build/bin. If you want to build the Python bindings, call:

.. code-block:: bash

    make pygimli

You might add J=8 (`make pygimli J=8`) for using 8 jobs in parallel to speed up
the build (adapt this to the number of real cores of the computer). The library
_pygimli_.so library will be copied into the source path
**../gimli/pygimli** in the subdirectory core.

To use the gimli installation you need to set some environment variables (this
example assumes that the **src** directory resides in your home directory):

.. code-block:: bash

    GIMLI_INSTALLATION=$HOME/src/gimli
    export PYTHONPATH=$PYTHONPATH:$GIMLI_INSTALLATION/gimli
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GIMLI_INSTALLATION/build/lib
    export PATH=$PATH:$GIMLI_INSTALLATION/build/bin

If you want to use the C++ command line applications, call

.. code-block:: bash

    make apps

Compiled binaries will be written to `build/bin`.

You can do a quick test of the pygimli build and installation with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.__version__)'

You can test your gimli build with:

.. code-block:: bash

    make check

Note that the test will be very silent if you don't have *cppunit* installed.

If you install pytest with

.. code-block:: bash

    sudo apt-get install python-pytest python3-pytest

then you can run the internal test suite with

.. code-block:: bash

    python -c "import pygimli; pygimli.test()"


Using Docker to build in Debian Bookworm (for advanced users only!)
...................................................................

If you want to use a Docker container to build (and possibly use) pyGIMLi, you
can use the Dockerfile found in the
`scripts/ <https://github.com/gimli-org/gimli/tree/master/core/scripts>`_
subdirectory named *Dockerfile_DebianBookworm*. Please refer to the file for
further instructions.


Troubleshooting
...............

If you experience runtime problems on starting pygimli like:

.. code-block:: bash

    ImportError: /usr/lib/libboost_python.so: undefined symbol: PyClass_Type

It may happen that CMake estimates the wrong libboost_python version by choosing py2 version instead of py3.
You can force cmake to select the correct version with:

.. code-block:: bash

    cmake ../gimli -DBoost_PYTHON_LIBRARY=/usr/lib64/libboost_python3.so

castXML
.......

castXML (https://github.com/CastXML/CastXML/) is needed to generate the code for the python bindings.
Some systems provide castxml binary so the build system should detect it if installed.
As fallback solution the build system tries to install castxml binaries or try to compile there own if the binaries don't work.
You can enforce the local binary installation with:

.. code-block:: bash

    cmake ../../src/castXML/ -DCASTXML_LOCAL=1
    make

or the local binary compilation with:

.. code-block:: bash

    cmake ../../src/castXML/ -DCASTXML_LOCALSRC=1
    make


If castXML build complains about missing clang or llvm command, please go into
$(GIMLISRC)/../thirdParty/build-XXX-XXX/castXML and try configure and build cmake manually

.. code-block:: bash

    CC=clang-3.6 CXX=clang++-3.6 cmake ../../src/castXML/
    make

If you build castXML manually you can provide this binary to cmake via

.. code-block:: bash

    cmake ../gimli -DCASTER_EXECUTABLE=$(PATH_TO_CASTXML)


Useful cmake settings
.....................

You can rebuild and update all local generated third party software by setting
the CLEAN environment variable:

.. code-block:: bash

    CLEAN=1 cmake ../gimli

Use alternative c++ compiler.

.. code-block:: bash

    CC=clang CXX=clang++ cmake ../gimli

Define alternative python version. On default the version of your active python
version will be chosen. You will need numpy and boost-python builds with your
desired python version.

.. code-block:: bash

    cmake ../gimli -DPYVERSION=3.6

Build the library with debug and profiling flags

.. code-block:: bash

    cmake ../gimli -DCMAKE_BUILD_TYPE=Debug

Build the library with gcc build.in sanity check

.. code-block:: bash

    cmake ../gimli -DCMAKE_BUILD_TYPE=Debug -DASAN=1


Useful make commands
....................

More verbose build output to view the complete command line:

.. code-block:: bash

    make VERBOSE=1
