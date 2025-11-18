#!/bin/env bash

pacman --noconfirm -Suy

echo "#########################################################################"
echo "### Maybe you need to restart msys after update of system components! ###"
echo "#########################################################################"

pacman -S --needed --noconfirm \
        make \
        tar \
        git \
        unzip \
        wget \
        diffutils \
        patch 

# default cmake does not recognize WIN32 settings
pacman -r --noconfirm cmake

# essentials
pacman -S --needed --noconfirm \
        mingw-w64-x86_64-cmake \
        mingw-w64-x86_64-gcc \
        mingw-w64-x86_64-gcc-fortran \
        mingw-w64-x86_64-openblas \
        mingw-w64-x86_64-suitesparse \
        mingw-w64-x86_64-llvm \
        mingw-w64-x86_64-clang

pacman -S --noconfirm \
        mingw-w64-x86_64-cppunit \
        mingw-w64-x86_64-doxygen

# we need to replace these hard coded paths from the LLVM cmake config
#sed -i 's/C:\/repo\/mingw-w64-clang\/src\/build-x86_64/C:\/msys64\/mingw64/'  /mingw64/share/llvm/cmake/LLVMConfig.cmake
#sed -i 's/C:\/repo\/mingw-w64-clang\/src\/llvm-3.7.0.src\/cmake\/modules/C:\/msys64\/mingw64\/share\/llvm\/cmake/' /mingw64/share/llvm/cmake/LLVMConfig.cmake

#sed -i 's/FATAL_ERROR "The imported target/WARNING "The imported target/' /mingw64/share/llvm/cmake/LLVMExports.cmake
