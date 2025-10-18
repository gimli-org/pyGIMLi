#!/usr/bin/env bash
# General workflow for continuous integration script
# Assumes sources being up to date

# leave any virtual environment
deactivate 2>/dev/null || true

if [ ! -z "$VIRTUAL_ENV" ]; then
    RED
    echo "active virtual environment $VIRTUAL_ENV : leave first"
    NCOL
    return
fi

# bash -xe ./gimli/core/scripts/jenkins-nf.sh
if [[ $- == *i* ]]; then
    # run with . *.sh or source *.sh
    echo "Running in interactive mode. Aliases will work"
else
    # run with bash *.sh
    echo "Running on NON interactive mode. Aliases will not work. Setting them now."
    shopt -s expand_aliases ## else aliases will be ignored in this bash session

    # make the script exit on every fail
    set -e
    alias python='python3'
    alias return='exit'
fi

function GREEN(){
    echo -e '\033[0;32m'
}
function BLUE(){
    echo -e '\033[0;34;49m'
}
function RED(){
    echo -e '\033[0;31;49m'
}
function NCOL(){
    echo -e '\033[0m'
}

function use_venv(){
    local venv_path=$1
    if [ -d $venv_path ]; then
        GREEN
        echo "*** Activating virtual environment $venv_path ***"
        NCOL
        source $venv_path/bin/activate
    else
        GREEN
        echo "*** Virtual environment $venv_path does not exist. Creating it. ***"
        NCOL
        new_venv $venv_path
    fi
}
function new_venv(){
    local venv_path=$1
    GREEN
    echo "*** Creating fresh virtual environment $venv_path ***"
    NCOL
    python --version
    if [ -x deactivate ]; then
        echo "leave old virtual environment"
        deactivate
    fi
    rm -rf $venv_path
    python -m venv $venv_path
    use_venv $venv_path

    echo "Updating pip ..."
    python -m pip install --upgrade pip
    pip install uv
}

function clean(){
    GREEN
    echo "*** Cleaning ...                 ***"
    NCOL

    pushd $PROJECT_ROOT
        echo "clear pygimli cache"
        rm -rf ~/.cache/pygimli
        echo "remove SOURCE_BUILDS:"
        rm -rf $PROJECT_SRC/build
        rm -rf $PROJECT_SRC/dist
        rm -rf $PROJECT_SRC/*.egg-info
        echo "remove BUILD_DIR: $BUILD_DIR"
        rm -rf $BUILD_DIR
        echo "remove BUILD_VENV: $BUILD_VENV"
        rm -rf $BUILD_VENV
        echo "remove TEST_VENV: $TEST_VENV"
        rm -rf $TEST_VENV
        echo "remove DOC_VENV: $DOC_VENV"
        rm -rf $DOC_VENV
    popd
}

function build_pre(){
    GREEN
    echo "*** Prepare building environment ***"
    NCOL

    pushd $PROJECT_ROOT
        new_venv $BUILD_VENV
        echo "pip install -e $PROJECT_SRC/[build]"
        uv pip install -e $PROJECT_SRC/[build]
        rm -rf $BUILD_DIR
        mkdir -p $BUILD_DIR
    popd
}

function build(){
    GREEN
    echo "*** Building now ... ***"
    NCOL

    pushd $PROJECT_ROOT

        if [ ! -d $BUILD_DIR ]; then
            build_pre
        fi
        use_venv $BUILD_VENV
        pushd $BUILD_DIR

            cmake $PROJECT_SRC
            make -j $GIMLI_NUM_THREADS
            make pygimli J=$GIMLI_NUM_THREADS
            make whlTest
        popd
    popd
}

function build_post(){
    GREEN
    echo "*** Testing build ***"
    NCOL
    pushd $PROJECT_ROOT
        use_venv $BUILD_VENV
        python -c 'import pygimli; print(pygimli.version())'
        python -c 'import pygimli; print(pygimli.Report())'
    popd
}

function test_pre(){
    GREEN
    echo "*** Prepare testing environment ***"
    NCOL

    pushd $PROJECT_ROOT
        new_venv $TEST_VENV
        if [ ! -d $BUILD_DIR ]; then
            build_pre
            build
            build_post
        fi
        uv pip install -e $PROJECT_SRC/[test]
    popd
}
function test(){
    GREEN
    echo "*** Testing now ... ***"
    NCOL

    pushd $PROJECT_ROOT
        if [ ! -d $TEST_VENV ]; then
            test_pre
        else
            use_venv $TEST_VENV
        fi

        python -c 'import pygimli; print(pygimli.version())'
        python -c 'import pygimli; print(pygimli.Report())'
        python -c 'import pygimli; pygimli.test(show=False, abort=True)'
    popd
}

function doc_pre(){
    GREEN
    echo "*** Preparing documentation ... ***"
    NCOL

    pushd $PROJECT_ROOT
        new_venv $DOC_VENV
        if [ ! -d $BUILD_DIR ]; then
            build_pre
            build
            build_post
        fi
        uv pip install $PROJECT_SRC/[doc]
        uv pip install pygimli ## remove rudimentary pygimli in doc venv
        uv pip install $BUILD_DIR/wheelhouse/pgcore*manylinux*.whl
        uv pip install $BUILD_DIR/wheelhouse/pygimli*.whl
        python -c 'import pygimli; print(pygimli.version())'
        python -c 'import pygimli; print(pygimli.Report())'
    popd
}

function doc(){
    GREEN
    echo "*** Creating documentation ***"
    NCOL

    pushd $PROJECT_ROOT
        use_venv $DOC_VENV
        if [ ! -d $BUILD_DIR ]; then
            doc_pre
        fi
        pushd $BUILD_DIR
            touch CMakeCache.txt # to renew search for sphinx
            cmake $PROJECT_SRC
            make clean-gallery
            if [ -x xvfb-run > 2/dev/null || true ]; then
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

function help(){
    echo ""
    echo run: ${BASH_SOURCE[0]} "help | clean | build_pre | build | build_post | test_pre | test |
    | doc_pre | doc | all"
    echo ""
}
function all(){
    clean
    build
    test
    doc
}

# Show system information
lsb_release -d
uname -a
env

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
elif [ -z $WORKSPACE ]; then
    WORKSPACE=$(realpath $(pwd))
    GREEN
    echo "Local Build (no Jenkins) on WORKSPACE=$WORKSPACE"
    NCOL
else
    GREEN
    echo "CI Build on WORKSPACE=$WORKSPACE"
    NCOL
    WORKSPACE=$(realpath $(pwd))
fi

if [ -z $SOURCE_DIR ]; then
    SOURCE_DIR=gimli
    echo "Using SOURCE_DIR=$SOURCE_DIR"
else
    echo "Using SOURCE_DIR=$SOURCE_DIR (forced by env setting SOURCE_DIR)"
fi

PROJECT_ROOT=$WORKSPACE
PROJECT_SRC=$PROJECT_ROOT/$SOURCE_DIR

python --version

if [ -z $PYVERSION ]; then
    PYVERSION=$(python -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')
    echo "building for python: $PYVERSION"
else
    echo "building for python: $PYVERSION (forced by setting PYVERSION)"
fi

GIMLI_NUM_THREADS=$((`nproc --all` - 2))

OPENBLAS_CORETYPE="ARMV8"
DISPLAY=':99.0'
PYVISTA_OFF_SCREEN=True

BUILD_VENV=$(realpath $WORKSPACE/venv-build-py$PYVERSION)
BUILD_DIR=$(realpath $WORKSPACE/build-py$PYVERSION)
TEST_VENV=$(realpath $WORKSPACE/venv-test-py$PYVERSION)
DOC_VENV=$(realpath $WORKSPACE/venv-doc-py$PYVERSION)


echo "WORKSPACE=$WORKSPACE"
echo "JOB_NUM=$JOB_NUM"
echo "PROJECT_SRC=$PROJECT_SRC"
echo "OS=$OS"
echo "NUM_THREADS=$GIMLI_NUM_THREADS"
echo "BUILD_VENV=$BUILD_VENV"
echo "BUILD_DIR=$BUILD_DIR"
echo "TEST_VENV=$TEST_VENV"
echo "DOC_VENV=$DOC_VENV"

echo "Starting automatic build #$BUILD_NUMBER on" `date`
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