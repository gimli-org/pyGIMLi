cmake_minimum_required(VERSION 3.15)

if (POLICY CMP0025) # recognize appleclang
    cmake_policy(SET CMP0025 NEW)
endif()
if (POLICY CMP0094)
    cmake_policy(SET CMP0094 NEW) #Python_FIND_STRATEGY LOCATION
endif()

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed.
            Please make a new directory (called a build directory)
            and run CMake from there. You may need to remove CMakeCache.txt." )
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/core/cmake")
include(GIMLImacros)

add_custom_target(pygimli)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Sets the configuration to build (Release, Debug, etc...)")
endif()

project(libgimli)
include(CheckIncludeFileCXX)

message(STATUS "CMAKE_SYSTEM:           ${CMAKE_SYSTEM}")
message(STATUS "CMAKE_SYSTEM_NAME:      ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_SYSTEM_VERSION:   ${CMAKE_SYSTEM_VERSION}")
message(STATUS "CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "CMAKE_SYSTEM_INFO_FILE: ${CMAKE_SYSTEM_INFO_FILE}")

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    execute_process (
        COMMAND bash -c "awk -F= '/^ID=/{print $2}' /etc/os-release |tr -d '\n' | tr -d '\"'"
            OUTPUT_VARIABLE OS_VAR
    )
    message(STATUS "OS: ${OS_VAR}")
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(OS_VAR "Windows")
    message(STATUS "OS: ${OS_VAR}")
elseif(CMAKE_SYSTEM_NAME MATCHES "MacOS")
    set(OS_VAR "MacOS")
    message(STATUS "MAC OS: ${OS_VAR}")
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(OS_VAR "Darwin")
    message(STATUS "Darwin: ${OS_VAR}")
else()
    message(FATAL "No OS detected!")
endif()

message(STATUS "UNIX:   " ${UNIX} " on all UNIX-like OS's, including Apple OS X and CygWin")
message(STATUS "WIN32:  " ${WIN32} " on Windows. Prior to 2.8.4 this included CygWin")
message(STATUS "APPLE:  " ${APPLE} " on Apple systems. Does not imply the system is Mac OS X.")
message(STATUS "MINGW:  " ${MINGW} " using the MinGW compiler in Windows")
message(STATUS "MSYS:   " ${MSYS} "  using the MSYS developer environment in Windows")
message(STATUS "CYGWIN: " ${CYGWIN} " using the CygWin version of cmake")
message(STATUS "MSVC:   " ${MSVC} " ")

message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")

#if (APPLE) // check if needed!
#    set(CMAKE_MACOSX_RPATH 1)
#    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#endif()

# Check if conda package is created
if (DEFINED ENV{CONDA_BUILD})
    message(STATUS "Conda package is being created.")
    message(STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
    set(CONDA_BUILD TRUE)
    set(Boost_INCLUDE_DIR "${CMAKE_PREFIX_PATH}/include")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-class-memaccess") # new since gcc.8 maybe need to be checked, not on apple
else()
    set(CONDA_BUILD FALSE)
endif()

message("***********************************************")
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message("***********************************************")

if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANGXX)

    # c++17 index access to std::share_ptr
    add_definitions(-std=c++17)
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-Os")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-comment")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-class-memaccess")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-vla")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-int-in-bool-context")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs") # pg build
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations") # loot of boost complaints about bind
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-overloaded-virtual")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_BIND_GLOBAL_PLACEHOLDERS")

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -ansi -pedantic -fno-omit-frame-pointer -Wall")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -std=c++11 -ansi -pedantic -fno-omit-frame-pointer -Wall")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp")

    ### debug settings
    # run with -DCMAKE_BUILD_TYPE=Debug
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -p -g -ansi -pedantic -fno-common -fno-omit-frame-pointer -Wall")

    if (WIN32 OR CONDA_BUILD)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-parentheses") # from boost since gcc.8 maybe need to be checked
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-attributes") # from gcc-8.2 + boost
    endif()

    # special section for clang compiler
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-overloaded-virtual")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Qunused-arguments")

        if (APPLE)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-warning-option") # on apple clang
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template") # on apple clang
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-infinite-recursion")
            # If it's libc++ and you're on <= 10.8, you need to compile with
            # clang++ -stdlib=libc++. If it's libstdc++ and you're on 10.9 or later,
            # you need to compile with clang++ -stdlib=libstdc++.
            add_compile_options(-stdlib=libstdc++)
            add_compile_options(-stdlib=libc++)
        endif()
    else() # if not clang but gcc

    endif()

    if (NOT WIN32 AND ASAN)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    endif()

    if (APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-infinite-recursion")
    endif()
elseif(MSVC)
	# MSVC complain a lot of possible unsecure std function
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

find_package(Git)
if (GIT_FOUND)
    message(STATUS ${GIT_FOUND})
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIMLI_GIT_TAG
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
        )
    message(STATUS "getting tag from git description: GIMLI_GIT_TAG=${GIMLI_GIT_TAG}")

    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIMLI_GIT_TAG_FULLVERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
        )
    message(STATUS "getting full version from git description: GIMLI_GIT_TAG_FULLVERSION=${GIMLI_GIT_TAG_FULLVERSION}")

    execute_process(COMMAND ${GIT_EXECUTABLE} describe --all --contains
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIMLI_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "getting branch name from git description: GIMLI_GIT_BRANCH=${GIMLI_GIT_BRANCH}")

    if (NOT GIMLI_GIT_BRANCH STREQUAL "master")
        set(LIBGIMLI_VERSION "${GIMLI_GIT_TAG_FULLVERSION}(${GIMLI_GIT_BRANCH})")
    else()
        set(LIBGIMLI_VERSION ${GIMLI_GIT_TAG_FULLVERSION})
    endif()
    message(STATUS "setting core version from git description to: ${LIBGIMLI_VERSION}")

    if (WHL_POST)
        set(LIBGIMLI_VERSION_WHL "${GIMLI_GIT_TAG}.post${WHL_POST}")
    else()
        set(LIBGIMLI_VERSION_WHL "${GIMLI_GIT_TAG}")
    endif()
else()
    set(LIBGIMLI_VERSION "untagged")
    set(LIBGIMLI_VERSION_WHL "untagged")
endif()

set(PACKAGE_NAME  \"${PROJECT_NAME}\" )
set(PACKAGE_VERSION  \"${LIBGIMLI_VERSION}\" )
set(PACKAGE_BUGREPORT  \"carsten@pygimli.org\")
set(PACKAGE_AUTHORS  \"carsten@pygimli.org thomas@pygimli.org florian@pygimli.org\")

################################################################################
# Check for libs and other packages we might use.
################################################################################
if (NOT THIRDPARTY_DIR)
    if (NOT ADDRESSMODEL)
        if ("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
            message(STATUS "Target is 64 bits")
            set (ADDRESSMODEL "64")
        else()
            message(STATUS "Target is 32 bits")
            set (ADDRESSMODEL "32")
        endif()
    endif()

    set(TARGETNAME "-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}-${ADDRESSMODEL}")
    set(ENV{TARGETNAME} ${TARGETNAME})
    set(THIRDPARTY_DIR ${PROJECT_SOURCE_DIR}/../thirdParty/)

    get_filename_component(EXTERNAL_DIR "${THIRDPARTY_DIR}/dist${TARGETNAME}" ABSOLUTE)
    get_filename_component(EXTERNAL_BUILD_DIR "${THIRDPARTY_DIR}/build${TARGETNAME}" ABSOLUTE)
    message(STATUS "ThirdParty set to: ${THIRDPARTY_DIR}")
    message(STATUS "External set to: ${EXTERNAL_DIR}")
endif()

if (NOT EXTERNAL_DIR)
set(EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external/)
message(STATUS "External set to: ${EXTERNAL_DIR}")
endif()

find_package(Threads REQUIRED)

################################################################################
#  Python
################################################################################
if (NOT PYVERSION)
    # Q&D python version guess .. while Python_FIND_STRATEGY LOCATION does not work in every case
    execute_process(COMMAND python --version
        OUTPUT_VARIABLE
            PYVERSION_GUESS
        OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_STRIP_TRAILING_WHITESPACE
    )
    STRING(REPLACE "Python " "" PYVERSION_GUESS "${PYVERSION_GUESS}")
    message(STATUS "No PYVERSION set, guessing from current python executable: ${PYVERSION_GUESS}")
    set(PYVERSION ${PYVERSION_GUESS})
endif ()

cprint(green "********** PYTHON settings ***********")
message(STATUS "Searching for python (${PYVERSION})" )
set(Python_FIND_STRATEGY LOCATION) # from cmake 3.15
set(Python_ADDITIONAL_VERSIONS=${PYVERSION})
set(Python_FIND_VIRTUALENV FIRST)

if (OS_VAR STREQUAL "almalinux")
    ## almalinux used for manylinux builds,
    ## EXACT fails since Python_Includes are not part of
    ## build virtual enviroment but on /opt
    find_package(Python ${PYVERSION} COMPONENTS Interpreter REQUIRED)
    find_package(Python ${PYVERSION} COMPONENTS Development.Module REQUIRED)
    find_package(Python ${PYVERSION} COMPONENTS NumPy REQUIRED)
elseif (WIN32)
    find_package(Python ${PYVERSION} COMPONENTS Interpreter REQUIRED)
    #py312 find python env, but Development fails, however include and lib seems okay
    find_package(Python ${PYVERSION} COMPONENTS Development)
    find_package(Python ${PYVERSION} COMPONENTS NumPy REQUIRED)
else()
    find_package(Python ${PYVERSION} EXACT COMPONENTS Interpreter REQUIRED)
    find_package(Python ${PYVERSION} EXACT COMPONENTS Development REQUIRED)
    find_package(Python ${PYVERSION} EXACT COMPONENTS NumPy REQUIRED)
endif()

set (PYTHON_FOUND TRUE CACHE BOOL "Is python available")

message(STATUS "Python_FOUND: ${Python_FOUND}")
message(STATUS "Python_Interpreter_FOUND: ${Python_Interpreter_FOUND}")
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")
message(STATUS "Python_INTERPRETER_ID: ${Python_INTERPRETER_ID}")
message(STATUS "Python_STDLIB: ${Python_STDLIB}")
message(STATUS "Python_STDARCH: ${Python_STDARCH}")
message(STATUS "Python_SITELIB: ${Python_SITELIB}")
message(STATUS "Python_SITEARCH: ${Python_SITEARCH}")
message(STATUS "Python_SOABI: ${Python_SOABI}")
message(STATUS "Python_INCLUDE_DIRS: ${Python_INCLUDE_DIRS}")
message(STATUS "Python_Development_FOUND: ${Python_Development_FOUND}")
message(STATUS "Python_Development.Embed_FOUND: ${Python_Development.Embed_FOUND}")
message(STATUS "Python_Development.Embed: $<TARGET_FILE:Python::Python>")
message(STATUS "Python_Development.Module_FOUND: ${Python_Development.Module_FOUND}")
message(STATUS "Python_Development.Module: $<TARGET_FILE:Python::Module>")
message(STATUS "Python_LIBRARIES: ${Python_LIBRARIES}")
message(STATUS "Python_LINK_OPTIONS: ${Python_LINK_OPTIONS}")
message(STATUS "Python_VERSION: ${Python_VERSION}")
message(STATUS "Python_NumPy: ${Python_NumPy_FOUND}, ${Python_NumPy_INCLUDE_DIRS}, ${Python_NumPy_VERSION}")


################################################################################
#  Boost
################################################################################
cprint(green "********** Boost.Python settings ***********")
#set(Boost_USE_STATIC_LIBS ON)

if (BUILD_BOOST)
    message(STATUS "Forced boost build.")
    build_package(Boost boost)
endif()

# weirdly: cannot run prior or some internal cache setting prevent manual search later on almalinux??
# if (NOT Boost_PYTHON_LIBRARY)
#     find_package(Boost COMPONENTS "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")

#     if (NOT Boost_FOUND)
#         unset(Boost_INCLUDE_DIR CACHE)
#         unset(Boost_INCLUDE_DIR)
#     endif()
# endif()
message(STATUS "Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")

#exit()
function(find_boost_root_from_distname_file)
    set (BOOST_DIST_HINT_FILE "${EXTERNAL_DIR}/.boost-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}.dist")
    message(STATUS "locking for ${BOOST_DIST_HINT_FILE}")

    if (EXISTS ${BOOST_DIST_HINT_FILE})
        message(STATUS "Found: ${BOOST_DIST_HINT_FILE}")
        FILE(READ "${BOOST_DIST_HINT_FILE}" BOOST_DIST_NAME)
        STRING(REGEX REPLACE "\n" "" BOOST_DIST_NAME "${BOOST_DIST_NAME}")
        get_filename_component(Boost_ROOT "${EXTERNAL_DIR}/${BOOST_DIST_NAME}" ABSOLUTE)
    else()
        message(STATUS "there is no ${BOOST_DIST_HINT_FILE}")
    endif()

    if (EXISTS ${Boost_ROOT})
        set(Boost_ROOT ${Boost_ROOT} CACHE FILEPATH "Boost root distribution")
        set(Boost_INCLUDE_DIR ${Boost_ROOT}/include CACHE FILEPATH "Boost inlude path")
        get_filename_component(Boost_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ABSOLUTE)
        message(STATUS "Boost_ROOT from hint file: ${Boost_ROOT}")
        message(STATUS "Boost_INCLUDE_DIR from hint file: ${Boost_INCLUDE_DIR}")
    else()
        message(STATUS "Boost_ROOT does not exist: ${Boost_ROOT}")
    endif()
endfunction(find_boost_root_from_distname_file)

message(STATUS "Trying to guess boost installation:")

if (NOT Boost_INCLUDE_DIR OR NOT Boost_ROOT)
    message(STATUS "Boost_INCLUDE_DIR not set")

    if (Boost_ROOT)
        message(STATUS "Found Boost_ROOT: ${Boost_ROOT}")
        get_filename_component(Boost_ROOT "${Boost_ROOT}" ABSOLUTE)
        set(Boost_INCLUDE_DIR ${Boost_ROOT}/include CACHE FILEPATH "Boost inlude path")
        message(STATUS "set Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
    else()
        message(STATUS "Boost_ROOT not set: building local")
        if (NOT Boost_INCLUDE_DIR)
            message(STATUS " .. building local without Boost_INCLUDE_DIR")
            find_boost_root_from_distname_file()
            find_or_build_package(Boost boost)
            if (NOT Boost_ROOT)
                find_boost_root_from_distname_file()
            endif()
        else()
            message(STATUS " .. building local with Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
            find_boost_root_from_distname_file()
            find_or_build_package_check(Boost boost Boost_INCLUDE_DIR False)
            find_boost_root_from_distname_file()
        endif()
    endif(Boost_ROOT)
endif(NOT Boost_INCLUDE_DIR OR NOT Boost_ROOT)

message(STATUS "Found Boost_ROOT: ${Boost_ROOT}")
message(STATUS "Found Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")
if (NOT Boost_PYTHON_LIBRARY)
    # Who sets this? It prevents successfull search directly after build
    unset(Boost_USE_DEBUG_RUNTIME)
    find_package(Boost COMPONENTS "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
    #find_or_build_package(Boost boost)
endif()

#exit()
mark_as_advanced(Boost_INCLUDE_DIR Boost_ROOT)

if (NOT Boost_PYTHON_LIBRARY)
    ################################################################################
    # Manuall check for boost python
    ################################################################################
    function(find_boost_python_manual bp_version)
        message(STATUS "----------- searching for libboost_python${bp_version}")
        find_package(Boost COMPONENTS "python${bp_version}")

        if (${Boost_PYTHON${bp_version}_FOUND})
            set(Boost_PYTHON_FOUND ${Boost_PYTHON${bp_version}_FOUND} PARENT_SCOPE)
            set(Boost_PYTHON_LIBRARY ${Boost_PYTHON${bp_version}_LIBRARY} PARENT_SCOPE)
            # set(Boost_PYTHON_FOUND ${Boost_PYTHON${bp_version}_FOUND} CACHE BOOL "boost python lib found")
            # set(Boost_PYTHON_LIBRARY ${Boost_PYTHON${bp_version}_LIBRARY_RELEASE} CACHE FILEPATH "boost python lib(find_boost_python_manual)")

            message(STATUS "----------- found: ${Boost_PYTHON_LIBRARY}")
        else()
            message(STATUS "----------- not found.")
        endif()
    endfunction()

    message(STATUS "boost-python not found by cmake ... trying to determine boost-python manually ...")
    if (WIN32)
        #set(bp_version "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}") # conda-forge
        set(bp_version "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}-mt-x64") # win64, boost-1.83
        #set(bp_version "python") # win64, boost-1.83
        string(TOUPPER ${bp_version} BP_VERSION)

        message(STATUS "... (win32) looking for ${bp_version} in Boost_ROOT=${Boost_ROOT}")
        find_package(Boost COMPONENTS "${bp_version}")

        ## local var evaluated before cache, and cache setting does not overwrite local.
        ## so we need to set booth here for repeated cmake calls
        if (Boost_PYTHON_FOUND)
            set(Boost_PYTHON_FOUND ${Boost_${bp_version}_FOUND} CACHE BOOL "boost python found")
            set(Boost_PYTHON_LIBRARY ${Boost_${BP_VERSION}_LIBRARY} CACHE FILEPATH "boost python lib")
        else()
            set(Boost_PYTHON_LIBRARY ${Boost_ROOT}/lib/libboost_${bp_version}.dll.a)
            message(STATUS "... (win32) looking for ${Boost_PYTHON_LIBRARY}")

            if (EXISTS ${Boost_PYTHON_LIBRARY})
                message(STATUS "... found.")
                set(Boost_PYTHON_FOUND TRUE CACHE BOOL "boost python found")
                set(Boost_PYTHON_FOUND TRUE)

                set(Boost_PYTHON_LIBRARY ${Boost_PYTHON_LIBRARY} CACHE FILEPATH "boost python lib")
            else()
                message(STATUS "Sorry but no boost-python found.")
            endif()
        endif()
    else()
        # the name for py3 boost-python library lacks probably py3 suffix, which is different for different OS
        #set (Boost_PYTHON_FOUND False)
        #set (Boost_PYTHON_LIBRARY "None")
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python37.so # gentoo >=boost-1.70)
            set(bp_version "${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            #/usr/lib/x86_64-linux-gnu/libboost_python-py34.so (debian)
            set(bp_version "-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python3-py37.so # weird debian)
            set(bp_version "${Python_VERSION_MAJOR}-py${Python_VERSION_MAJOR}${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python3.so.1.64.0 (arch style)
            set(bp_version "${Python_VERSION_MAJOR}")
            find_boost_python_manual(${bp_version})
        endif()
        if (NOT Boost_PYTHON_FOUND)
            # /usr/lib/libboost_python-3.4.so (gentoo <boost-1.70)
            set(bp_version "-${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}")
            find_boost_python_manual(${bp_version})
        endif()
    endif()

    #set(Boost_PYTHON_LIBRARY "/opt/homebrew/Caskroom/mambaforge/base/lib/libboost_python39.dylib")
    if (Boost_PYTHON_FOUND)
        message(STATUS "Boost_PYTHON_FOUND:" ${Boost_PYTHON_FOUND})
        message(STATUS "Boost_PYTHON_LIBRARY:" ${Boost_PYTHON_LIBRARY})
    else()
        message(FATAL_ERROR "Sorry but no boost-python found. Please search for "
                       "any libboost_python* on your system and "
                       "provide it to the authors so we can try to fix "
                       "this here.")
    endif()
endif()

if (NOT USE_BOOST_THREAD)
    set(Boost_THREAD_FOUND OFF)
    set(BOOST_BIND_FOUND OFF)
    set(USE_IPC OFF)
    set(THREADS "std")
else()
    status(ERROR "boost threads are not longer supported")
endif()

################################################################################
#   misc third party stuff
################################################################################
cprint(green "********** Misc settings ***********")

find_package(OpenMP)

# Find packages that CHOLMOD depends on
set(CMAKE_LIBRARY_PATH ${EXTERNAL_DIR}/lib $ENV{EXTERNAL_DIR}/lib ${CMAKE_LIBRARY_PATH})

if (J) # dummy to avoid error msg
endif()

set (UMFPACK_FOUND FALSE)

find_library(Triangle_LIBRARIES
            NAMES
                libtri.a
	            tri
                triangle
                libtriangle
            PATHS
                ${Triangle_PREFIX_PATH}
                ${EXTERNAL_DIR}
                ${PROJECT_SOURCE_DIR}/external
                ${PROJECT_BINARY_DIR}/external
            PATH_SUFFIXES
                lib
            )

message(STATUS "Triangle_LIBRARIES = ${Triangle_LIBRARIES}")

find_or_build_package(Triangle triangle LOCAL)
#find_package(Triangle)

if (CONDA_BUILD)
# (231115) there is a missing symlink to libopenblas.so in ${CMAKE_PREFIX_PATH}/lib
# so we need to search for libblas.so
    set(BLA_VENDOR Generic)
elseif (APPLE)
    ## check if openblas is used in runtime
    set(BLA_VENDOR OpenBLAS)
    #find_or_build_package(BLAS lapack)
else()
    set(BLA_VENDOR OpenBLAS)
endif()
#set(CMAKE_FIND_DEBUG_MODE TRUE)
find_package(BLAS REQUIRED)
set(CMAKE_FIND_DEBUG_MODE TRUE)
set(CMAKE_FIND_DEBUG_MODE FALSE)
# message(STATUS "BLAS_LINKER_FLAGS = ${BLAS_LINKER_FLAGS}")
# message(STATUS "BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

################################################################################
# Manuall check for openblas
################################################################################

#if (";${BLAS_LIBRARIES};" MATCHES "openblas")
message(STATUS "BLAS found:" ${BLAS_FOUND})
if (BLAS_FOUND)
    message(STATUS "OpenBLAS is used: ${BLAS_LIBRARIES}")

    if (MINGW)
        get_filename_component(MINGW_ROOT ${BLAS_LIBRARIES} PATH)
        get_filename_component(MINGW_ROOT ${MINGW_ROOT}/../ ABSOLUTE)
        message(STATUS "MINGW_ROOT: " ${MINGW_ROOT})

        FIND_PATH(OpenBLAS_INCLUDE_DIR
            NAMES
                cblas.h
            PATHS
                ${MINGW_ROOT}/include/OpenBLAS/
        )
    endif()

    if (NOT OpenBLAS_INCLUDE_DIR)
        message(STATUS "search for OpenBLAS cblas.h in prefered paths.")

        FIND_PATH(OpenBLAS_INCLUDE_DIR
            NAMES
                cblas.h
            PATHS
                ../../_build_env/include/ # for conda
                ${CMAKE_PREFIX_PATH}/include
                /usr/include/openblas/
                /usr/include/OpenBLAS/
                /mingw64/include/OpenBLAS/
            HINTS
                NO_DEFAULT_PATH
            )
        #set(CMAKE_FIND_DEBUG_MODE TRUE)
        #set(CMAKE_FIND_DEBUG_MODE FALSE)
        message(STATUS "OpenBLAS cblas.h: ${OpenBLAS_INCLUDE_DIR}")
    endif()

    if (NOT OpenBLAS_INCLUDE_DIR)
        message(STATUS "search for OpenBLAS cblas.h in default paths.")

        FIND_PATH(OpenBLAS_INCLUDE_DIR
            NAMES
                cblas.h
        )
        message(STATUS "OpenBLAS cblas.h: ${OpenBLAS_INCLUDE_DIR}")
    endif()

    if (OpenBLAS_INCLUDE_DIR)
        set(OPENBLAS_CBLAS_FOUND TRUE)
    else()
        set(OPENBLAS_CBLAS_FOUND FALSE)
    endif()

    set(OPENBLAS_FOUND TRUE)

else()
    set(OPENBLAS_FOUND FALSE)
    set(OPENBLAS_CBLAS_FOUND FALSE)
endif()

#exit()
#find_package(suitesparse)
#exit()

#find_or_build_package(CHOLMOD suitesparse)
find_package(UMFPACK REQUIRED)
find_package(CHOLMOD REQUIRED)

if (NOT NOCPPUNIT)
    find_package(CppUnit)
    if (CPPUNIT_FOUND)
        option (LIBGIMLI_BUILD_TESTS "Build unittests" ON)
    endif (CPPUNIT_FOUND)
endif()

if (NOT NOREADPROC)
    if(WIN32)
        set(READPROC_FOUND FALSE)
    else()
        find_package(readproc)
        message("READPROC_FOUND ${READPROC_FOUND}")
    endif()
else()
    set(READPROC_FOUND FALSE)
endif()

#exit()

if (NOT CASTER)
    set(CASTER "castxml")

    if (CMAKE_COMPILER_IS_GNUCC)
        string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${CMAKE_CXX_COMPILER_VERSION})
        list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
        list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)

        if (${GCC_MAJOR} EQUAL 5)
            message(STATUS "Found gcc version above 5: choosing castxml")
            set(CASTER "castxml")
        endif()

    endif()
endif()

find_program(CLANG_EXECUTABLE clang)
find_or_build_package(castxml castxml LOCAL)
set(CASTER_EXECUTABLE ${CASTXML_EXECUTABLE})

if (CASTER STREQUAL "gccxml")
    find_or_build_package(gccxml gccxml LOCAL)
    set(CASTER_EXECUTABLE ${GCCXML_EXECUTABLE})
else()
    if (NOT CASTER_EXECUTABLE)
        find_program(CASTXML castxml)

        if (NOT CASTXML OR CASTXML_LOCAL OR CASTXML_LOCALSRC)
            set(CASTER_EXECUTABLE ${CASTXML_EXECUTABLE})
        else()
            set(CASTER_EXECUTABLE ${CASTXML})
        endif()
    endif()
endif()

if(CLANG_EXECUTABLE)
    message(STATUS "clang++ found : ${CLANG_EXECUTABLE}")
else()
    message(FATAL_ERROR "Can't found program: clang")
endif()

find_or_build_package(castxml castxmlbin LOCAL)
set(CASTER_EXECUTABLE ${CASTXML_EXECUTABLE})

find_or_build_package(pygccxml pygccxml LOCAL)
find_or_build_package(pyplusplus pygccxml LOCAL)

if (CASTER_EXECUTABLE)
    set (CASTER_FOUND TRUE)
endif()

find_package(Doxygen)
find_package(Sphinx 3.1)

if (SPHINX_FOUND AND NOT SKIP_SPHINX)
    add_subdirectory(doc EXCLUDE_FROM_ALL)
    message(STATUS "Search recursive for related doc files from here: " ${CMAKE_SOURCE_DIR})

    set (RECDOCTYPES "*.rst" "*.md" "*.png" "*.svg" "*.bib" "paper/*.py")

    foreach (DTYPE ${RECDOCTYPES})
        file(GLOB_RECURSE DFILES
                LIST_DIRECTORIES false
                RELATIVE "${CMAKE_SOURCE_DIR}"
            "./doc/${DTYPE}")
        message(STATUS "Append DOC files: (${DTYPE}) " ${DFILES})
        list(APPEND DOCFILES ${DFILES})
    endforeach()

    set(RSTFILES ${DOCFILES} CACHE INTERNAL "RST source file that we need to copy")

    execute_process(COMMAND ${SPHINX_EXECUTABLE} --version
                    OUTPUT_VARIABLE SPHINX_VERSION
    )

    add_custom_target(removedoc
       COMMAND
       ${CMAKE_MAKE_PROGRAM} clean -f ${CMAKE_CURRENT_BINARY_DIR}/doc/
    )

endif()

message(STATUS "**********************************************************************")
message(STATUS "************************* Dependencies found *************************")
message(STATUS "**********************************************************************")
message(STATUS "THREADS              : ${THREADS} ${CMAKE_THREAD_LIBS_INIT} ${Boost_THREAD_LIBRARIES}")
message(STATUS "USE_OPENBLAS         : ${OPENBLAS_FOUND}, OPENBLAS_CBLAS_FOUND :${OPENBLAS_CBLAS_FOUND}")
message(STATUS "BLAS_LIBRARIES       : ${BLAS_LIBRARIES}")
message(STATUS "OpenBLAS_INCLUDE_DIR : ${OpenBLAS_INCLUDE_DIR}")
message(STATUS "CHOLMOD_LIBRARIES    : ${CHOLMOD_LIBRARIES}")
message(STATUS "UMFPACK_LIBRARIES    : ${UMFPACK_LIBRARIES}")
message(STATUS "TRIANGLE_FOUND       : ${TRIANGLE_FOUND} Triangle_LIBRARIES: ${Triangle_LIBRARIES}")
message(STATUS "Python_EXECUTABLE    : ${Python_EXECUTABLE}" )
message(STATUS "Python_Dev.Mod_FOUND : ${Python_Development.Module_FOUND}" )
message(STATUS "Python_Dev_FOUND     : ${Python_Development_FOUND}" )
message(STATUS "Python_INCLUDE_DIRS  : ${Python_INCLUDE_DIRS}")
message(STATUS "Boost_INCLUDE_DIR    : ${Boost_INCLUDE_DIR}")
message(STATUS "Boost_PYTHON_LIBRARY : ${Boost_PYTHON_LIBRARY}" )
message(STATUS "Python_NumPy         : ${Python_NumPy_FOUND} ${Python_NumPy_INCLUDE_DIRS} ver: ${Python_NumPy_VERSION}")
message(STATUS "CASTER_FOUND         : ${CASTER_FOUND} Caster: ${CASTER_EXECUTABLE}")
message(STATUS "PYGCCXML_FOUND       : ${PYGCCXML_FOUND} PYGCCXML: ${PYGCCXML_PATH}" )
message(STATUS "PYPLUSPLUS_FOUND     : ${PYPLUSPLUS_FOUND} PYPLUSPLUS_PATH: ${PYPLUSPLUS_PATH}" )
message(STATUS "")
message(STATUS "**********************************************************************")
message(STATUS "*************** Optional Dependencies found **************************")
message(STATUS "**********************************************************************")
message(STATUS "CPPUNIT             : ${CPPUNIT_FOUND} ${CPPUNIT_LIBRARIES}")
message(STATUS "Doxygen             : ${DOXYGEN_FOUND}")
message(STATUS "Sphinx              : ${SPHINX_FOUND} ${SPHINX_EXECUTABLE} ${SPHINX_VERSION}")
message(STATUS "**********************************************************************")

if (Python_EXECUTABLE
    AND TRIANGLE_FOUND
    AND (Python_Development.Module_FOUND OR Python_Development_FOUND)
    AND Boost_PYTHON_LIBRARY
    AND Python_NumPy_FOUND
    AND CASTER_FOUND
    AND PYGCCXML_FOUND
    AND PYPLUSPLUS_FOUND
)
    set (PYGIMLI 1)
    message(STATUS "")
    message(STATUS "pygimli can be build. run: make pygimli")
    message(STATUS "**********************************************************************")
    message(STATUS "")
else()
    #set (PYGIMLI 1)
    # temporary disable
    message (FATAL_ERROR "pygimli cannot be build due to some missing packages.
    Check above for python, castxml, pygccxml, pyplusplus, boost_python, numpy")
endif()

configure_file("${PROJECT_SOURCE_DIR}/core/config.cmake.h.in"
               "${PROJECT_BINARY_DIR}/config.cmake.h" )
add_definitions(-DHAVE_CONFIG_CMAKE_H)

################################################################################
# Define the configurable options
################################################################################
include(CMakeDependentOption)
cmake_dependent_option( GIMLI_BUILD_SHARED_LIBS_WITH_STATIC_DEPENDENCIES
"For WIN32 and APPLE where the dependency pack is used, prefer the static dependency libs over the shared/dynamic ones.  NOTE: On Windows you also need to be mindful of which C/C++ runtime setting has been used to compile the various
components - they must all match or you will get crashes, heap corruption and/or
other issues." FALSE "WIN32 OR APPLE" FALSE)

################################################################################
#
################################################################################
set (LIBRARY_INSTALL_DIR lib)
set (INCLUDE_INSTALL_DIR include/gimli/)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
if (WIN32)
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
else()
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# set up install sub-directories
if (CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS "${CMAKE_INSTALL_PREFIX}/lib64")
    set( GIMLI_LIB_INSTALL_DIR lib64 )
elseif (CMAKE_SIZEOF_VOID_P EQUAL 4 AND EXISTS "${CMAKE_INSTALL_PREFIX}/lib32")
    set( GIMLI_LIB_INSTALL_DIR lib32 )
else()
    set( GIMLI_LIB_INSTALL_DIR lib )
endif()

set(GIMLI_VER_INSTALL_SUBDIR "/${CMAKE_PROJECT_NAME}-${GIMLI_VERSION_MAJOR}" )
set(GIMLI_MODULE_INSTALL_DIR "${GIMLI_LIB_INSTALL_DIR}/${CMAKE_PROJECT_NAME}-${GIMLI_VERSION_MAJOR}.${GIMLI_VERSION_MINOR}" )
set(GIMLI_INCLUDE_INSTALL_DIR "include${GIMLI_VER_INSTALL_SUBDIR}" )
#set( GIMLI_PYTHON_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/pygimli)

################################################################################
# Add main header locations (for everything we build)
################################################################################
include_directories(${libgimli_SOURCE_DIR}/core/src)
include_directories("${PROJECT_BINARY_DIR}")

################################################################################
# add virtual env for whl package test
################################################################################
if (NOT TESTVENV)
    set(TESTVENV "${CMAKE_BINARY_DIR}/wheeltest_venv")
endif()

set(PYTHON_VIRTUAL_ENV ${TESTVENV})
if (NOT EXISTS "${PYTHON_VIRTUAL_ENV}")
    message(STATUS "Creating virtual python environment for testing [${PYTHON_VIRTUAL_ENV}]")
    execute_process(COMMAND "${Python_EXECUTABLE}" -m venv "${PYTHON_VIRTUAL_ENV}")
else ()
    message(STATUS "Using virtual python environment for testing [${PYTHON_VIRTUAL_ENV}]")
endif ()

# update the environment with VIRTUAL_ENV variable (mimic the activate script) and search again
set(ENV{VIRTUAL_ENV} "${PYTHON_VIRTUAL_ENV}")
set(Python3_FIND_VIRTUALENV ONLY)
unset(Python3_EXECUTABLE)
find_package(Python3 COMPONENTS Interpreter)

message(STATUS "Found python executable for building: ${Python_EXECUTABLE}")
message(STATUS "Found python executable for testing (venv): ${Python3_EXECUTABLE}")

set(WHEELHOUSE ${CMAKE_BINARY_DIR}/wheelhouse)

################################################################################
# define test and checks
################################################################################
# if (LIBGIMLI_BUILD_TESTS)
#      enable_testing()
# endif(LIBGIMLI_BUILD_TESTS)

################################################################################
# descend into subdirs
################################################################################
add_subdirectory(core/src)
add_subdirectory(core/tests EXCLUDE_FROM_ALL)

add_subdirectory(core/python EXCLUDE_FROM_ALL)
add_subdirectory(core/pgcore)
add_subdirectory(pygimli)

add_custom_target(check)
if (CPPUNIT_FOUND)
    add_dependencies(check gtest)
endif (CPPUNIT_FOUND)
add_dependencies(check pgtest)

add_custom_target(checkall DEPENDS doc)
