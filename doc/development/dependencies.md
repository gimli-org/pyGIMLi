.. _sec:install_pre:

Prerequisites
-------------

To build pyGIMLi from source, the following tools are required:

* subversion, git, mercurial, wget, tar
* cmake >= 2.8.10
* gcc >= 4.4
* >=Python-3.5 | >=Python-2.7
* numpy-dev
* >=matplotlib-3.0
* >=clang++-3.6.0 (3.7.0, 3.8.0)
* libclang-3.7-dev
* >=llvm-3.6.0 (3.7.0, 3.8.0)
* libz-dev
* python-setuptools

    tested on:

    * gentoo x86_64: gcc-4.4.5, gcc-4.5.3, gcc-4.5.4, gcc-4.9.2 gcc-5.3.0
    * Debian 3.2.46-1 x86_64: gcc-4.7.2
    * Ubuntu 16.04 LTS with gcc-5.4.0
    * Arch Linux gcc-5.2.0
    * CentOS
    * MinGW32: gcc-4.6.2/4, gcc-4.7.2, gcc-5.2.0
    * MinGW64: gcc-4.5.4, gcc-5.2.0, gcc-6.3.0, gcc-7.1.0

Optional Prerequisites
......................

These tools can be installed system-wide with your native package manager (i.e.
apt-get). If not found, the provided build scripts will automatically download
and compile the necessary components.

* libboost >=1.46 (python) (1.46, 1.48, 1.49, 1.51, 1.52, 1.53, 1.57)
* blas and lapack for suitesparse (system or auto via cmake)
* SuiteSparse (http://faculty.cse.tamu.edu/davis/suitesparse.html)
* cppunit
* procps
* triangle (http://www.cs.cmu.edu/~quake/triangle.html)

Prerequisites automatically installed by the installer
......................................................

These tools are required to create the Python bindings and are likely to be
outdated in your distribution and should therefore not be installed
system-wide. The build scripts will install them automatically.

* castxml
* pygccxml
* pyplusplus
