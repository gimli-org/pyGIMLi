(build_win)=
# Building on Windows


**This is heavily outdated, new version comming soon.**

You need a Python installation with the modules numpy and matplotlib. There are
minimalistic installers and pre-packaged distributions. We recommend `Anaconda
<http://www.continuum.io/>`_ 64bit but `WinPython
<http://winpython.github.io/#releases>`_ will to the job too.

You also need a Linux-like command shell along with a gcc compiler.
Although there might be different solutions (Cygwin, Git Bash, MinGW/MSYS),
we only support the MSYS2 (Minimal System 2) hosted at http://www.msys2.org.
As computers and modern Windows (>=7) are 64bit we only test this.
Avoid installing into strange Windows folders, e.g., c:\ProgramData\mys64 is fine.

After installing MSYS, start the console once so it builds your personal home
directory where you find a .bashrc file, e.g. in

.. code-block::

    c:\ProgramData\mys64\home\YOUR_USERNAME

Edit .bashrc so that the Anaconda or WinPython installation path is added to your default
PATH.

e.g.:

.. code-block:: bash

    export ANACONDA=/c/ProgramData/Anaconda3
    export PATH=$PATH:$ANACONDA

This is necessary since gimli needs to know valid python installation and
version. Ideally the following one-liner will suffice to compile pyGIMLi in the
current directory.

**Note: The script will automatically take care of requirements and updates of MSYS2.
And also needs to modify/patch some of the llvm system files.**

.. code:: bash

    curl -Ls install.pygimli.org | bash

This script accepts a few more options. For help see

.. code:: bash

    curl -Ls install.pygimli.org | bash -s help

If everything runs fine, including some tests, the script will tell you some
additional PATH and PYTHONPATH settings for your .bashrc to use pygimli inside
the console or any IDE like Spyder (coming along with WinPython).

If something goes wrong, please take a look on the error message.

You can alse try the following instructions for manual installation.


## Manual installation

Make sure to have an updated msys2 environment. Run at least:

.. code-block:: bash

    pacman -Su
    pacman -Sy

to update your local package databases. See https://www.msys2.org/
for further instructions.

To get a complete working toolchain you need some packages installed.

.. code-block:: bash

    pacman -S make tar git subversion mercurial unzip wget patch

.. code-block:: bash

    pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran
    pacman -S mingw-w64-x86_64-openblas mingw-w64-x86_64-suitesparse
    pacman -S mingw-w64-x86_64-doxygen mingw-w64-x86_64-cppunit
    pacman -S mingw-w64-x86_64-llvm mingw-w64-x86_64-clang

The rest of the installation is like the linux way with some small differences.

Prepare the directory structure as described above:

The build is performed via cmake. While calling cmake *MSYS* users should tell
using the Unix makefile generator to find the correct gcc versions:

.. code-block:: bash

    cmake ../gimli -G 'Unix Makefiles'

If cmake complains about missing python stuff, make sure the Python interpreter
is in your execution path.


**Problems with cmake configuration**

If cmake can't install pygccxml or pyplusplus then you can provide those packages using pip from the anaconda distribution.
First make sure the needed scripts are in your path.

.. code-block:: bash

    export PATH=$PATH:$ANACONDA/Scripts

Then you can install those both packages in your user space

.. code-block:: bash

   pip install pygccxml --user
   pip install pyplusplus --user

If cmake complains about misssig numpy, python can't probably import numpy, which you can test:

.. code-block:: bash

    python -c 'import numpy'

Probably anaconda additional needs another path setting, don't ask me why

.. code-block:: bash

   export PATH=$PATH:$ANACONDA/Library/bin

Now python should be able to find numpy and cmake will work as supposed and you can continue the build process.


To build the library, just run

.. code-block:: bash

    make -j2

You might add the option -jN to use a number of N CPUs in parallel.
To build pygimli, run

.. code-block:: bash

    make pygimli J=2

You might add J=N to use a number of N CPUs in parallel.
Building pygimli takes some time and you can grab a coffee (or two).
If it finishes without errors you just need to set the environment:
(note that pygimli is still built in-place, in pygimli/core)

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/gimli/python
    export PATH=$PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

Compiled binaries will be written to `build/bin`.

You can test the pygimli build with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.version())'
