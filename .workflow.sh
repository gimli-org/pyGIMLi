#!/usr/bin/env bash
# General workflow for continuous integration script
# Assumes sources being up to date

# leave any virtual environment
deactivate 2>/dev/null || true

function BLUE(){
    echo -e '\033[0;34;49m'
}
function GREEN(){
    echo -e '\033[0;32m'
}
function YELLOW(){
    echo -e '\033[0;33;49m'
}
function RED(){
    echo -e '\033[0;31;49m'
}
function NCOL(){
    echo -e '\033[0m'
}

if [[ $- == *i* ]]; then
    # run with . *.sh or source *.sh
    echo "Running in interactive mode. Aliases will work"
else
    # run with bash *.sh
    echo "Running on NON interactive mode. Aliases will not work. Setting them now."
    shopt -s expand_aliases ## else aliases will be ignored in this bash session

    # make the script exit on every fail
    set -e

    alias return='exit'
fi

if [ ! -z "$VIRTUAL_ENV" ]; then
    RED
    echo "active virtual environment $VIRTUAL_ENV : leave first"
    echo "call deactivate"
    NCOL
    return
fi

function use_venv(){
    local venv_path=$1
    if [ -d $venv_path ]; then
        GREEN
        echo "Activating virtual environment $venv_path"
        NCOL
        if [ $OS == "Windows" ] || [ $OS == "Windows_NT" ]; then
            # windows
            source $venv_path/Scripts/activate
        else
            # linux / mac
            source $venv_path/bin/activate
        fi
    else
        GREEN
        echo "Virtual environment $venv_path does not exist. Creating it."
        NCOL
        new_venv $venv_path
    fi
}


function new_venv(){
    local venv_path=$1
    GREEN
    echo "Creating fresh virtual environment $venv_path"
    NCOL
    $BASEPYTHON --version
    deactivate 2>/dev/null || true

    rm -rf $venv_path
    $BASEPYTHON -m venv $venv_path
    use_venv $venv_path

    echo "Updating pip ..."
    python -m pip install --upgrade pip
    pip install uv
}


function testReport(){
    if [ -n "$VIRTUAL_ENV" ]; then
        GREEN
        echo "Testing build in active virtual environment: $VIRTUAL_ENV"
        echo "Python: $(python -V 2>&1) ($(python -c 'import sys; print(sys.executable)'))"
        NCOL
    else
        YELLOW
        echo "No active virtual environment."
        echo "Python: $(python -V 2>&1) ($(python -c 'import sys; print(sys.executable)'))"
        NCOL
    fi

    python -c 'import pygimli; print(pygimli.version())'
    python -c 'import pygimli; print(pygimli.Report())'
}


function clean(){
    GREEN
    echo "*** Cleaning ...                 ***"
    NCOL

    pushd $PROJECT_ROOT
        echo "clear pygimli cache"
        rm -rf ~/.cache/pygimli
        echo "remove SOURCE_BUILDS:"
        rm -rf $PROJECT_DIST
        rm -rf $PROJECT_SRC/build
        rm -rf $PROJECT_SRC/dist
        rm -rf $PROJECT_SRC/*.egg-info
        echo "remove BUILD_DIR: $BUILD_DIR"
        rm -rf $BUILD_DIR
        echo "remove VENV_BUILD: $VENV_BUILD"
        rm -rf $VENV_BUILD
        echo "remove VENV_TEST: $VENV_TEST"
        rm -rf $VENV_TEST
        echo "remove VENV_DOC: $VENV_DOC"
        rm -rf $VENV_DOC
        echo "remove VENV_PYGIMLI: $VENV_PYGIMLI"
        rm -rf $VENV_PYGIMLI
    popd
}


function build_pre(){
    GREEN
    echo "*** build_pre (Prepare building environment) ***"
    NCOL

    ## clean previous build, dist, and doc artifacts
    clean

    pushd $PROJECT_ROOT
        new_venv $VENV_BUILD

        echo "uv pip install -e $PROJECT_SRC/[build]"
        uv pip install -e $PROJECT_SRC/[build]
        rm -rf $BUILD_DIR
        mkdir -p $BUILD_DIR
    popd
}


function build(){
    GREEN
    echo "*** build (Building now ...) ***"
    NCOL

    if [ ! -d $BUILD_DIR ]; then
        YELLOW
        echo "BUILD_DIR: $BUILD_DIR not found. Preparing now (build_pre)."
        NCOL
        build_pre
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_BUILD

        pushd $BUILD_DIR
            if [ "$OS" == "Windows" ] || [ "$OS" == "Windows_NT" ]; then

                ### windows MSYS2 configuration
                cmake -G "Unix Makefiles" $PROJECT_SRC
            elif [ "$OS" == "MacOS" ] || [ "$(uname -s)" == "Darwin" ]; then

                ### MacOS configuration
                echo "MacOS build with custom openblas and umfpack/cholmod paths. Expecting \$CMAKE_PREFIX to be set."
                cmake \
                    -DNOREADPROC=1 \
                    -DBLAS_openblas_LIBRARY=$CMAKE_PREFIX/lib/libopenblas.dylib \
                    -DUMFPACK_LIBRARIES=$CMAKE_PREFIX/lib/libumfpack.dylib \
                    -DUMFPACK_INCLUDES=$CMAKE_PREFIX/include/suitesparse \
                    -DCHOLMOD_LIBRARIES=$CMAKE_PREFIX/lib/libcholmod.dylib \
                    -DCHOLMOD_INCLUDE_DIRS=$CMAKE_PREFIX/include/suitesparse \
                    -DOpenBLAS_INCLUDE_DIR=$CMAKE_PREFIX/include \
                    $PROJECT_SRC
            else

                ### Linux configuration
                cmake $PROJECT_SRC
            fi

            make -j $GIMLI_NUM_THREADS
            make pygimli J=$GIMLI_NUM_THREADS

            build_whls
    popd
}


function build_whls(){
    GREEN
    echo "*** build_wheels (Building wheels now ...) ***"
    NCOL

    if [ "$OS" == "Windows" ] || [ "$OS" == "Windows_NT" ]; then
        LIBGIMLI=$(ls $BUILD_DIR/bin/libgimli.dll | head -n 1)
        echo "Checking for libgimli in $BUILD_DIR/bin: $LIBGIMLI"
    elif [ "$OS" == "MacOS" ] || [ "$(uname -s)" == "Darwin" ]; then
        LIBGIMLI=$(ls $BUILD_DIR/lib/libgimli.dylib | head -n 1)
        echo "Checking for libgimli in $BUILD_DIR/lib: $LIBGIMLI"
    else
        LIBGIMLI=$(ls $BUILD_DIR/lib/libgimli.so | head -n 1)
        echo "Checking for libgimli in $BUILD_DIR/lib: $LIBGIMLI"
    fi

    if [ -z $LIBGIMLI ]; then
        YELLOW
        echo "libgimli not found in $BUILD_DIR/lib. Building now (build)."
        NCOL
        build
    else
        GREEN
        echo "libgimli found: $LIBGIMLI."
        NCOL
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_BUILD

        pushd $BUILD_DIR

        # create pgcore wheel
            make whlpgcoreCopyLibs

            pushd $BUILD_DIR/core/pgcore
                WHEELHOUSE=$BUILD_DIR/wheelhouse
                mkdir -p $WHEELHOUSE
                mkdir -p $BUILD_DIR/dist/

                python -m pip wheel --wheel-dir=$WHEELHOUSE .

                WHLFILE=$(ls $WHEELHOUSE/pgcore*.whl | head -n 1)

                if [ "$OS" == "MacOS" ] || [ "$(uname -s)" == "Darwin" ]; then
                    BLUE
                    echo "Repairing pgcore whl for MacOS ($WHLFILE)"
                    NCOL
                    delocate-wheel -v $WHLFILE

                elif [ "$OS" == "Windows" ] || [ "$OS" == "Windows_NT" ]; then
                    BLUE
                    echo "Repairing pgcore whl for Windows ($WHLFILE)"
                    NCOL
                    delvewheel repair $WHLFILE --add-path $BUILD_DIR/bin/

                else
                    BLUE
                    echo "Repairing pgcore whl for Linux ($WHLFILE)"
                    NCOL
                    echo "Build pgcore whl for Linux"
                fi
            popd

            GREEN
            echo "Copying pgcore whl ($WHLFILE) to build dist $BUILD_DIR/dist/"
            cp $WHLFILE $BUILD_DIR/dist/
            NCOL
        popd

        ### create pygimli wheel
        pushd $PROJECT_SRC
            python -m build --wheel --no-isolation --outdir $BUILD_DIR/dist/
        popd
    popd
}


function build_post(){
    GREEN
    echo "*** build_post (Testing build) ***"
    NCOL

    PGCORE_WHL=$(ls $BUILD_DIR/dist/pgcore*.whl | head -n 1)
    PG_WHL=$(ls $BUILD_DIR/dist/pygimli*.whl | head -n 1)

    if [ -z $PGCORE_WHL ] || [ -z $PG_WHL ]; then
        YELLOW
        echo "pgcore or pygimli whl's not found in $BUILD_DIR/dist. Building now (build_whls)"
        NCOL
        build_whls
    else
        GREEN
        echo "Whl's found in build dist:"
        echo -e "\t pgcore: $PGCORE_WHL"
        echo -e "\t pygimli: $PG_WHL"
        NCOL
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_BUILD

        # special case for windows .. pgcore install to ensuse mingw runtime libs are found
        if [ "$OS" == "Windows" ] || [ "$OS" == "Windows_NT" ]; then
            # windows MSYS2
            uv pip install $BUILD_DIR/dist/pgcore*.whl
        fi

        testReport

        mkdir -p $PROJECT_DIST

        GREEN
        echo "Copying built whl files to project dist: $PROJECT_DIST"
        NCOL
        cp $BUILD_DIR/dist/pgcore*.whl $PROJECT_DIST/
        cp $BUILD_DIR/dist/pygimli*.whl $PROJECT_DIST/
    popd
}


function test_pre(){
    GREEN
    echo "*** test_pre (Prepare testing environment) ***"
    NCOL

    if [ ! -f $PROJECT_DIST/pgcore*.whl ]; then
        YELLOW
        echo "pgcore wheel not found in project dist. Building first. (build_post)"
        NCOL
        build_post
    fi

    pushd $PROJECT_ROOT
        new_venv $VENV_TEST
        # not needed to install pgcore in editable after build for linux

        # special case for windows .. pgcore install to ensuse mingw runtime libs are found
        if [ "$OS" == "Windows" ] || [ "$OS" == "Windows_NT" ]; then
            # windows MSYS2
            uv pip install $PROJECT_DIST/pgcore*.whl
        fi
        uv pip install -e $PROJECT_SRC/[test]
    popd
}


function test(){
    GREEN
    echo "*** test (Testing now ...) ***"
    NCOL

    pushd $PROJECT_ROOT
        if [ ! -d $VENV_TEST ]; then
            test_pre
        else
            use_venv $VENV_TEST
        fi

        testReport
        python -c 'import pygimli; pygimli.test(show=False, abort=True)'
    popd
}


function doc_pre(){
    GREEN
    echo "*** doc_pre (Preparing documentation ...) ***"
    NCOL

    if [ ! -f $PROJECT_DIST/pgcore*.whl ]; then
        YELLOW
        echo "pgcore wheel not found in project dist. Building first. (build_post)"
        NCOL
        build_post
    fi

    pushd $PROJECT_ROOT
        new_venv $VENV_DOC

        # TODO find a way to install the whl file with optional deps
        uv pip install $PROJECT_SRC/[doc]
        uv pip uninstall pygimli ## remove rudimentary pygimli in doc venv
        uv pip install --force-reinstall $PROJECT_DIST/pgcore*.whl
        uv pip install $PROJECT_DIST/pygimli*.whl

        testReport
    popd
}

function doc(){
    GREEN
    echo "*** doc (Creating documentation) ***"
    NCOL

    use_venv $VENV_DOC
    if python -c 'import sphinx' &>/dev/null; then
        GREEN
        echo "sphinx is installed"
        NCOL
    else
        YELLOW
        echo "sphinx is NOT installed (calling doc_pre)"
        NCOL
        doc_pre
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_DOC

        pushd $BUILD_DIR
            #touch CMakeCache.txt # to renew search for sphinx
            cmake $PROJECT_SRC

            python -c 'import pygimli as pg; pg.version()'

            #make clean-gallery # should not be necessary
            if [ -x "$(command -v xvfb-run)" ]; then
                # xvfb is necessary for headless display of pyvista plots
                echo "xvfb-run available: using it to build docs"
                xvfb-run make doc
            else
                echo "xvfb-run not available: building docs without it"
                make doc
            fi
        popd
    popd
}

function doc_post(){
    GREEN
    echo "*** doc_post (Deploying html) ***"
    NCOL

    if [ ! -f $BUILD_DIR/doc/_build/html/index.html ]; then
        YELLOW
        echo "Documentation html not found in dist. Building first. (doc)"
        NCOL
        doc
    fi

    pushd $BUILD_DIR/doc/_build
        cp -r html $PROJECT_DIST/html
        tar -czvf $PROJECT_DIST/html.tgz html
    popd

    #rsync -aP $DISTPATH/html user@pygimli.org:DEV_HTML_PATH --delete

    # SNAP_PATH=`date +"%Y%m%d"`
    # mkdir -p ~/snapshots/$SNAP_PATH
    # rm -f ~/snapshots/latest
    # ln -s ~/snapshots/$SNAP_PATH ~/snapshots/latest
    # cp $PROJECT_DIST/*.whl ~/snapshots/latest/
    # source ~/snapshots/venv-oscar-py310/bin/activate
    # python -m pip install --force ~/snapshots/latest/pgcore*
    # python -m pip install --force ~/snapshots/latest/pygimli*
}

function install(){
    GREEN
    echo "*** Install (Creating editable installation in build venv) ***"
    NCOL

    if [ ! -f $PROJECT_DIST/pgcore*.whl ]; then
        YELLOW
        echo "pgcore wheel not found in project dist. Building first.(build_post)"
        NCOL
        build_post
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_PYGIMLI

        uv pip install $PROJECT_DIST/pgcore*.whl
        uv pip install -e $PROJECT_SRC/[opt]

        testReport
    popd

    GREEN
    echo "Editable installation created in venv: $VENV_PYGIMLI"
    echo ""
    echo "To use it, call: "
    echo ""
    echo "source $VENV_PYGIMLI/bin/activate #(linux/macos)"
    echo "source $VENV_PYGIMLI/Scripts/activate #(windows)"
    NCOL
}

function help(){
    echo ""
    echo "run: ${BASH_SOURCE[0]} TARGET"
    echo ""
    echo "TARGETS:"
    echo "    help       show this help"
    echo "    clean      remove build artifacts for PYVERSION"
    echo "    build_pre  prepare build environment venv"
    echo "    build      [build_pre] build project"
    echo "    build_whls [build] build whl files for pgcore and pygimli"
    echo "    build_post [build] test build"
    echo "    test_pre   [build] prepare test environment"
    echo "    test       [test_pre] run tests"
    echo "    doc_pre    [build] prepare documentation environment"
    echo "    doc        [doc_pre] build documentation"
    echo "    doc_post   [doc] deploy documentation"
    echo "    install    [build] Create default editable installation in venv"
    echo "    all        [clean build test doc]"
    echo ""
    echo "ENVIRONMENT variables:"
    echo "    BASEPYTHON   base python interpreter. Default system python3."
    echo "    PYVERSION    python version (e.g. 3.11, 3.14t) if no BASEPYTHON is given"
    echo "    SOURCE_DIR   source directory (default: top-level directory of the project)"
    echo ""
    echo "Examples:"
    echo "    bash ${BASH_SOURCE[0]} clean build test doc"
    echo "    PYVERSION=3.12 bash ${BASH_SOURCE[0]} all"
    echo ""
}
function all(){
    clean
    build_post
    test
    doc
}

JOB_NUM=0

if [ ! -z $GITHUB_ACTIONS ]; then
    GREEN
    echo "GITHUB action runner on WORKSPACE=$WORKSPACE at $RUNNER_NAME"
    echo "RUNNER_ARCH=$RUNNER_ARCH"
    echo "RUNNER_TEMP=$RUNNER_TEMP"
    echo "GITHUB_REF_NAME=$GITHUB_REF_NAME"
    echo "GITHUB_JOB=$GITHUB_JOB"
    echo "GITHUB_REPOSITORY=$GITHUB_REPOSITORY"
    echo "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"
    NCOL

    WORKSPACE=$GITHUB_WORKSPACE
    JOB_NUM=$GITHUB_RUN_NUMBER
    OS=$RUNNER_OS
    SYSTEM=$(uname -s)
elif [ -z $WORKSPACE ]; then
    WORKSPACE=$(realpath $(pwd))
    if [ -z $OS ]; then
        uname_s=$(uname -s 2>/dev/null || echo)
        case "$uname_s" in
            MINGW*|MSYS*|CYGWIN*|Windows_NT)
                OS=Windows
                ;;
            Darwin*)
                OS=MacOS
                ;;
            *)
                OS=Linux
                ;;
        esac
    fi
    SYSTEM=$(uname -s)
    GREEN
    echo "Local Build (no Jenkins) on WORKSPACE=$WORKSPACE"
    NCOL
else
    GREEN
    echo "Unknown CI Build on WORKSPACE=$WORKSPACE"
    NCOL
    WORKSPACE=$(realpath $(pwd))
fi

if [ -z $SOURCE_DIR ]; then
    # derive SOURCE_DIR from location of this script (top-level project dir)
    SCRIPT_REALPATH=$(realpath "${BASH_SOURCE[0]}")
    SCRIPT_DIR=$(dirname "$SCRIPT_REALPATH")
    CAND="$SCRIPT_DIR"
    # climb up until we find a project marker (.git, pyproject.toml)
    while [ "$CAND" != "/" ] && [ ! -e "$CAND/.git" ] && [ ! -e "$CAND/pyproject.toml" ]; do
        CAND=$(dirname "$CAND")
    done
    if [ "$CAND" = "/" ]; then
        # fallback: use the immediate directory name of the script
        RED
        echo "Could not find project root marker (.git or pyproject.toml). Using script directory name as SOURCE_DIR."
        NCOL
    else
        SOURCE_DIR=$(basename "$CAND")
    fi
    echo "Using SOURCE_DIR=$SOURCE_DIR"
else
    echo "Using SOURCE_DIR=$SOURCE_DIR (forced by env setting SOURCE_DIR)"
fi

function abspath() {
    local p="$1"
    local abs

    # Prefer realpath if available
    if command -v realpath >/dev/null 2>&1; then
        abs=$(realpath -m "$p" 2>/dev/null) || abs=$(realpath "$p" 2>/dev/null)
    fi

    # Fallback to python if realpath not present or failed
    if [ -z "$abs" ]; then
        RED
        echo "Could not determine absolute path for $p"
        return
        NCOL
    fi

    # Normalize empty result
    abs=${abs:-$p}
    # Convert MSYS (/c/...) paths to Windows style (c:/...)
    if [[ "$abs" =~ ^/([A-Za-z])/(.*) ]]; then
        local drive="${BASH_REMATCH[1],,}"
        local rest="${BASH_REMATCH[2]}"
        abs="${drive}:/${rest}"
    fi

    # Normalize backslashes to forward slashes (pip accepts forward slashes on Windows)
    abs="${abs//\\//}"

    echo "$abs"
}

PROJECT_ROOT=$WORKSPACE
PROJECT_SRC=$(abspath $PROJECT_ROOT/$SOURCE_DIR)

if [ -z $BASEPYTHON ]; then
    if [ -z $PYVERSION ]; then

        if [ ! -x "$(command -v python3)" ]; then
            RED
            echo "python3 not found in PATH. Please install python3 or set BASEPYTHON to a valid python interpreter."
            echo "e.g., BASEPYTHON=../../miniconda3/python bash ${BASH_SOURCE[0]}"
            NCOL
            return
        else
            echo "Using system python3 as BASEPYTHON"
            BASEPYTHON=python3
        fi
    else
        echo "building for python: $PYVERSION (forced by setting PYVERSION)"
        BASEPYTHON=python$PYVERSION
    fi
fi

PYVERSION=$($BASEPYTHON -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')

if [ "$OS" == "MacOS" ] || [ "$(uname -s)" == "Darwin" ]; then
    GIMLI_NUM_THREADS=$((`sysctl -n hw.ncpu` - 2))
    alias realpath='grealpath'
else
    GIMLI_NUM_THREADS=$((`nproc --all` - 2))
fi

OPENBLAS_CORETYPE="ARMV8"
DISPLAY=':99.0'
PYVISTA_OFF_SCREEN=True

BUILD_DIR=$(realpath $WORKSPACE/build-py$PYVERSION)
VENV_BUILD=$(realpath $WORKSPACE/venv-build-py$PYVERSION)
VENV_TEST=$(realpath $WORKSPACE/venv-test-py$PYVERSION)
VENV_DOC=$(realpath $WORKSPACE/venv-doc-py$PYVERSION)
VENV_PYGIMLI=$(realpath $WORKSPACE/venv-pygimli-py$PYVERSION)
PROJECT_DIST=$(realpath $WORKSPACE/dist-py$PYVERSION)

echo "JOB_NUM=$JOB_NUM"
echo "WORKSPACE=$WORKSPACE"
echo "PROJECT_SRC=$PROJECT_SRC"
echo "OS=$OS"
echo "SYSTEM=$SYSTEM"
echo "NUM_THREADS=$GIMLI_NUM_THREADS"
echo "PYTHON=$BASEPYTHON"
echo "PYVERSION=$PYVERSION"
echo "VENV_BUILD=$VENV_BUILD"
echo "BUILD_DIR=$BUILD_DIR"
echo "VENV_TEST=$VENV_TEST"
echo "VENV_DOC=$VENV_DOC"
echo "VENV_PYGIMLI=$VENV_PYGIMLI"
echo "PROJECT_DIST=$PROJECT_DIST"

echo -e "\nStarting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")

[ $# -lt 1 ] && echo "No workflow target specified" && help

for arg in $@
do
    case $arg in
    clean)
        clean;;
    build_pre)
        build_pre;;
    build)
        build;;
    build_whls)
        build_whls;;
    build_post)
        build_post;;
    test_pre)
        test_pre;;
    test)
        test;;
    doc_pre)
        doc_pre;;
    doc)
        doc;;
    doc_post)
        doc_post;;
    install)
        install;;
    all)
        all;;
    help)
        help
        return;;
    *)
        echo "Don't know what to do."
        help;;
    esac
done

end=$(date +"%s")
echo "Ending automatic build #$BUILD_NUMBER".
diff=$(($end-$start))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

#deactivate 2>/dev/null || true