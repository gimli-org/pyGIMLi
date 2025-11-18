(build_mac)=
# Building on MAC OSX

The current working solution is based on `this discussion on GitHub
<https://github.com/gimli-org/gimli/discussions/603>`_.  Many thanks to Robin
Thibaut!

```bash

    # Install Homebrew
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

    # Install dependencies via brew
    brew install cmake
    brew install wget
    brew install mercurial

    # Install Miniforge
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh

    conda create -n pygimli_env python
    conda activate pygimli_env
    conda install -c conda-forge boost numpy scipy matplotlib openblas suitesparse

    # Clone the pygimli repository
    git clone https://github.com/gimli-org/gimli.git source
    mkdir build
    cd build

    PYTHON_EXEC=$(which python3)
    PYTHON_INC=$(python3 -c 'import sysconfig; print(sysconfig.get_path("include"))')
    PYTHON_LIB=$(python3-config --configdir)

    export CPLUS_INCLUDE_PATH=$PYTHON_INC

    cmake -DPYTHON_EXECUTABLE=$PYTHON_EXEC -DPYTHON_LIBRARY=$PYTHON_LIB -DPYTHON_INCLUDE_DIR=$PYTHON_INC ../source
    make -j 8
    make pygimli J=8
```

## Troubleshooting

If you encounter problems, you may have to specify some paths manually, e.g.:

```bash
    cmake -DPYTHON_EXECUTABLE=$PYTHON_EXEC -DPYTHON_LIBRARY=$PYTHON_LIB
          -DPYTHON_INCLUDE_DIR=$PYTHON_INC \
          -DUMFPACK_LIBRARIES=~/minforge3/base/lib/libumfpack.dylib \
          -DUMFPACK_INCLUDES=~/minforge3/base/include \
          -DCHOLMOD_LIBRARIES=~/minforge3/base/lib/libcholmod.dylib \
          -DCHOLMOD_INCLUDE_DIRS=~/minforge3/base/include \
          -DBLAS_openblas_LIBRARY=~/minforge3/base/lib/libopenblas.dylib \
          -DOpenBLAS_INCLUDE_DIR=~/minforge3/base/include \
          -DBoost_PYTHON_LIBRARY=~/minforge3/base/lib/libboost_python310.dylib \
          ../source
```
