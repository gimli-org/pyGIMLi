# - Try to find CHOLMOD
# Once done this will define
#
#  CHOLMOD_FOUND        - system has CHOLMOD
#  CHOLMOD_INCLUDE_DIRS - include directories for CHOLMOD
#  CHOLMOD_LIBRARIES    - libraries for CHOLMOD
#  HAVE_LIBCHOLMOD      - system has CHOLMOD

#=============================================================================
# Copyright (C) 2010-2011 Garth N. Wells, Anders Logg and Johannes Ring
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================
#
# Modified 2013 for libgimli by Carsten
#
#
message(STATUS "Checking for package 'CHOLMOD'")

# Find packages that CHOLMOD depends on
# set(CMAKE_LIBRARY_PATH ${BLAS_DIR}/lib $ENV{BLAS_DIR}/lib ${CMAKE_LIBRARY_PATH})
# set(CMAKE_LIBRARY_PATH ${LAPACK_DIR}/lib $ENV{LAPACK_DIR}/lib ${CMAKE_LIBRARY_PATH})
# set(CMAKE_LIBRARY_PATH ${EXTERNAL_DIR}/lib $ENV{EXTERNAL_DIR}/lib ${CMAKE_LIBRARY_PATH})

find_package(AMD QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(ParMETIS 4.0.2 QUIET)

find_path(SUITESPARSE_DIR SuiteSparse_demo.m
	${EXTERNAL_DIR}/SuiteSparse
	)


# Check for header file
find_path(CHOLMOD_INCLUDE_DIRS 
    cholmod.h
    HINTS 
        ../../_build_env/include # for conda
	      ${EXTERNAL_DIR}/include 
        $ENV{EXTERNAL_DIR}/include
	      ${SUITESPARSE_DIR}/CHOLMOD/include
    PATH_SUFFIXES 
        suitesparse 
        ufsparse
    DOC "Directory where the CHOLMOD header is located"
 )

# Check for CHOLMOD library
find_library(CHOLMOD_LIBRARY
    cholmod
    HINTS 
        ../../_build_env/lib # for conda
        ${EXTERNAL_DIR}/lib
        $ENV{EXTERNAL_DIR}/lib
    DOC 
        "The CHOLMOD library"
)

# Check for SUITESPARSECONFIG library
find_library(SUITESPARSE_LIBRARY suitesparseconfig
  HINTS ${EXTERNAL_DIR}/lib ${CCOLAMD_DIR}/lib $ENV{EXTERNAL_DIR}/lib $ENV{CCOLAMD_DIR}/lib
  DOC "The SUITESPARSECONFIG library"
  )

find_library(AMD_LIBRARY amd
  HINTS ${EXTERNAL_DIR}/lib ${AMD_DIR}/lib $ENV{EXTERNAL_DIR}/lib $ENV{AMD_DIR}/lib
  DOC "The AMD library"
  )

# Check for CAMD library
find_library(CAMD_LIBRARY camd
  HINTS ${EXTERNAL_DIR}/lib ${CAMD_DIR}/lib $ENV{EXTERNAL_DIR}/lib $ENV{CAMD_DIR}/lib
  DOC "The CAMD library"
  )

# Check for COLAMD library
find_library(COLAMD_LIBRARY colamd
  HINTS ${EXTERNAL_DIR}/lib ${COLAMD_DIR}/lib $ENV{EXTERNAL_DIR}/lib $ENV{COLAMD_DIR}/lib
  DOC "The COLAMD library"
  )

# Check for CCOLAMD library
find_library(CCOLAMD_LIBRARY ccolamd
  HINTS ${EXTERNAL_DIR}/lib ${CCOLAMD_DIR}/lib $ENV{EXTERNAL_DIR}/lib $ENV{CCOLAMD_DIR}/lib
  DOC "The CCOLAMD library"
  )


if (CHOLMOD_LIBRARY)
    set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARY})

    if (0)
      #needed? 
      # Collect libraries (order is important)
      if (AMD_FOUND)
          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${AMD_LIBRARIES})
      endif()
      if (CAMD_LIBRARY)
          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${CAMD_LIBRARY})
      endif()
      if (AMD_LIBRARY)
          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${AMD_LIBRARY})
      endif()
      if (COLAMD_LIBRARY)	
          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${COLAMD_LIBRARY})
      endif()
      if (CCOLAMD_LIBRARY)
          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${CCOLAMD_LIBRARY})
      endif()
      if (SUITESPARSE_LIBRARY)
          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${SUITESPARSE_LIBRARY})
      endif()

      # Don't link against system-wide blas when making conda package
      if (NOT ENV{CONDA_BUILD})
          message(STATUS "adding LAPACK_LIBRARIES to CHOLMOD PATH: ${LAPACK_LIBRARIES}")

          if (PARMETIS_FOUND)
              set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${PARMETIS_LIBRARIES})
          endif()
          if (LAPACK_FOUND)
              set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${LAPACK_LIBRARIES})
          endif()
          if (BLAS_FOUND)
              set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES})
          endif()
      endif()

      if (CONDA_BUILD)
          find_library(GFORTRAN_LIBRARY gfortran
              DOC "The gfortran library"
          )

          set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${GFORTRAN_LIBRARY})
      else()

          find_program(GFORTRAN_EXECUTABLE gfortran)
          if (GFORTRAN_EXECUTABLE)
              execute_process(COMMAND ${GFORTRAN_EXECUTABLE} -print-file-name=libgfortran.so
                              OUTPUT_VARIABLE GFORTRAN_LIBRARY
                              OUTPUT_STRIP_TRAILING_WHITESPACE)
              if (EXISTS "${GFORTRAN_LIBRARY}")
                  set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${GFORTRAN_LIBRARY})
              endif()
          endif(GFORTRAN_EXECUTABLE)
      endif()

      IF(WIN32)
    
      ELSE(WIN32)
      # On unix system, debug and release have the same name
          if (APPLE)
          else(APPLE)
              set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} rt)
          endif(APPLE)
      ENDIF(WIN32)
    endif()
endif ()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CHOLMOD
  "CHOLMOD could not be found. Be sure to set EXTERNAL_DIR."
  CHOLMOD_LIBRARIES CHOLMOD_INCLUDE_DIRS)

message(STATUS "CHOLMOD_LIBRARIES: ${CHOLMOD_LIBRARIES}")

mark_as_advanced(
  CHOLMOD_INCLUDE_DIRS
  CHOLMOD_LIBRARY
  CHOLMOD_LIBRARIES
  CAMD_LIBRARY
  AMD_LIBRARY
  COLAMD_LIBRARY
  CCOLAMD_LIBRARY
)
