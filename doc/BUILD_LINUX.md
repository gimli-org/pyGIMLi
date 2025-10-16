(build_lin)=
# Building on Linux

## Compiling with a virtual environment

If you don't want to use the conda environment we encourage the use of a
virtual python environment.
Assuming you have a proper build toolchain and the required libraries
(see Installation on Ubuntu below) installed.

First we need to create a root directory for our installation,
e.g., `$HOME/src/gimli` and get the source code:

```bash
    mkdir -p gimli
    cd gimli
    git clone https://github.com/gimli-org/gimli.git
    #optionally change to the development branch
    git checkout dev
    # create a virtual environment for pyGIMLi, this can be at any place with any name
    # if you want easy VScode support consider `gimli/.venv`
    deactivate # in case there is another venv active
    python -m venv venv-build --prompt=gimli-build
    # activate the venv
    source venv-build/bin/activate
    # update pip is allways a good idea
    python -m pip install -U pip
    # install gimli as editable with its dependencies into to venv
    pip install -e ./gimli[all]
```

We need to compile the C++ part of pyGIMLi, which is currently done with cmake and not
with the pip internal build system.
We recommend an out of source build:

```bash
    mkdir -p build
    cd build
    cmake ../gimli
    make -j 4 gimli
    make pygimli J=4
    cd ..
```

There is no more need to change the `PATH` or `LD_LIBRARY_PATH`.
In fact, check to remove any prior changes to these environment variables from
older installations.
If the build was successful he copy the compiled libraries back into the
source tree which is already editable known to the venv,
so you can test the pygimli build with:

```bash
    python -c 'import pygimli as pg; pg.version()'
    python -c 'import pygimli as pg; print(pg.Report())'
```

As long nothing changes in the C++ part of pygimli you can just update pyGIMLi
but just pulling the latest changes from git.
If you end the terminal session you can reactivate the venv with:

```bash
    source venv-build/bin/activate
```

## Example Installation on Ubuntu

Last try on Ubuntu 22.04.03 (23-11-14)

```bash
    sudo apt-get install build-essential g++ subversion git cmake \
                 python3-dev python3-matplotlib python3-numpy python3-pyqt5 \
                 python3-scipy libboost-all-dev libedit-dev \
                 libsuitesparse-dev libopenblas-openmp-dev libumfpack5 \
                 libomp-dev doxygen \
                 libcppunit-dev clang
```
Rest see above.

## Useful cmake settings

You can rebuild and update all local generated third party software by setting
the CLEAN environment variable:

```bash
CLEAN=1 cmake ../gimli
```

Use an alternative C++ compiler:
```bash
CC=clang CXX=clang++ cmake ../gimli
```

Build the library with debug and profiling flags:
```bash
cmake ../gimli -DCMAKE_BUILD_TYPE=Debug
```

Build the library with GCC address sanitizer check:
```bash
cmake ../gimli -DCMAKE_BUILD_TYPE=Debug -DASAN=1
```

## Useful make commands

Show more verbose build output:
```bash
make VERBOSE=1
```
