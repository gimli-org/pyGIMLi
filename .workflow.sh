#!/usr/bin/env bash
# General workflow for continuous integration script
# Assumes sources being up to date

# leave any virtual environment
deactivate 2>/dev/null || true

function green(){
    # green message
    echo -e '\033[0;32m'$1'\033[0m'
}

function blue(){
    # blue message
    echo -e '\033[0;34;49m'$1'\033[0m'
}

function yellow(){
    # yellow message
    echo -e '\033[0;33;49m'$1'\033[0m'
}

function red(){
    # red message
    echo -e '\033[0;31;49m'$1'\033[0m'
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
    red "Active virtual environment $VIRTUAL_ENV - leave first"
    red "Please call 'deactivate'"
    return
fi

function use_venv(){
    local venv_path=$1
    if [ -d $venv_path ]; then
        green "Activating virtual environment $venv_path"

        if [ $OS == "Windows" ] || [ $OS == "Windows_NT" ]; then
            # windows
            source $venv_path/Scripts/activate
        else
            # linux / mac
            source $venv_path/bin/activate
        fi
    else
        green "Virtual environment $venv_path does not exist. Creating it."
        new_venv $venv_path
    fi
}


function new_venv(){
    local venv_path=$1
    green "Creating fresh virtual environment $venv_path"
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
        green "Testing build in active virtual environment: $VIRTUAL_ENV"
        green "Python: $(python -V 2>&1) ($(python -c 'import sys; print(sys.executable)'))"
    else
        yellow "No active virtual environment."
        yellow "Python: $(python -V 2>&1) ($(python -c 'import sys; print(sys.executable)'))"
    fi

    python -c 'import pygimli; print(pygimli.version())'
    python -c 'import pygimli; print(pygimli.Report())'
}


function clean(){
    green "*** clean  (Remove old build artifacts) ***"

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
    green "*** build_pre (Prepare building environment) ***"
    ## clean previous build, dist, and doc artifacts
    clean

    pushd $PROJECT_ROOT
        new_venv $VENV_BUILD

        # this prevents pgcore being in dependencies since it can't be resolved
        # for new platforms or new python versions. Maybe add manual installation
        # of build prerequisites here if necessary.

        echo "Check if pgcore can be installed from PyPI"
        if uv pip install pgcore --dry-run 2>/dev/null; then
            green "pgcore is available on PyPI"
            echo "uv pip install -e $PROJECT_SRC/[build]"
            uv pip install -e $PROJECT_SRC/[build]
        else
            yellow "pgcore is not available for the $PYVERSION platform ${OS} on PyPI."
            yellow "Installing build dependencies manually"

            uv pip install "setuptools>=75.8.2"
            uv pip install "numpy>=2.1.3"
            uv pip install "pygccxml"
            uv pip install "pyplusplus"
            uv pip install build twine wheel
            uv pip install "auditwheel; sys_platform == 'linux'"
            uv pip install "delvewheel; sys_platform == 'win32'"
            uv pip install "delocate; sys_platform == 'darwin'"
            uv pip install "scooby"
        fi

        rm -rf $BUILD_DIR
        mkdir -p $BUILD_DIR
    popd
}


function build(){
    green "*** build (Building now ...) ***"

    if [ ! -d $BUILD_DIR ]; then
        yellow "BUILD_DIR: $BUILD_DIR not found. Preparing now (build_pre)."
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
    green "*** build_wheels (Building wheels now ...) ***"

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
        yellow "libgimli not found in $BUILD_DIR/lib. Building now (build)."
        build
    else
        green "libgimli found: $LIBGIMLI."
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_BUILD

        ### create pygimli wheel
        pushd $PROJECT_SRC
            python -m build --wheel --no-isolation --outdir $BUILD_DIR/dist/
        popd

        pushd $BUILD_DIR

        # create pgcore wheel
            make whlpgcoreCopyLibs

            pushd $BUILD_DIR/core/pgcore
                WHEELHOUSE=$BUILD_DIR/wheelhouse
                mkdir -p $WHEELHOUSE
                mkdir -p $BUILD_DIR/dist/

                python -m pip wheel --wheel-dir=$WHEELHOUSE .

                WHLFILE=$(ls $WHEELHOUSE/pgcore*.whl | head -n 1)

                if [ ! -z "$AUDITWHEEL_POLICY" ] && [ ! -z "$AUDITWHEEL_PLAT" ] ; then
                    blue "Repairing pgcore whl for $AUDITWHEEL_POLICY ($WHLFILE)"
                    auditwheel repair $WHLFILE -w $BUILD_DIR/dist/

                else
                    if [ "$OS" == "MacOS" ] || [ "$(uname -s)" == "Darwin" ]; then
                        blue "Repairing pgcore whl for $OS ($WHLFILE)"
                        delocate-wheel -v $WHLFILE

                    elif [ "$OS" == "Windows" ] || [ "$OS" == "Windows_NT" ]; then
                        blue "Repairing pgcore whl for $OS ($WHLFILE)"
                        delvewheel repair $WHLFILE --add-path $BUILD_DIR/bin/
                    elif [ ! -z "$AUDITWHEEL_POLICY" ] && [ ! -z "$AUDITWHEEL_PLAT" ] ; then
                        yellow "Unknown OS ($OS) for repairing pgcore whl."
                    fi

                    green "Copying pgcore whl ($WHLFILE) to build dist $BUILD_DIR/dist/"
                    cp $WHLFILE $BUILD_DIR/dist/
                fi
            popd

        popd

    popd
}


function install_WHL_E(){
    opt=$1
    green "*** install pygimli $opt from whl files (editable) ***"

    pushd $PROJECT_ROOT
        uv pip install $PROJECT_DIST/pgcore*.whl
        uv pip install --editable $PROJECT_SRC$opt
    popd
    testReport
}


function install_WHL(){
    opt=$1
    green "*** install pygimli $opt from whl files (non editable) ***"

    pushd $PROJECT_ROOT
        uv pip uninstall pygimli pgcore
        uv pip install --force-reinstall $PROJECT_DIST/pgcore*.whl
        WHLFILE=$(ls $PROJECT_DIST/pygimli*.whl | head -n 1)
        uv pip install "$WHLFILE$opt"
    popd
    testReport
}


function build_post(){
    green "*** build_post (Testing build) ***"

    PGCORE_WHL=$(ls $BUILD_DIR/dist/pgcore*.whl | head -n 1)
    PG_WHL=$(ls $BUILD_DIR/dist/pygimli*.whl | head -n 1)

    if [ -z $PGCORE_WHL ] || [ -z $PG_WHL ]; then
        yellow "pgcore or pygimli whl's not found in $BUILD_DIR/dist. Building now (build_whls)"
        build_whls
    else
        green "Whl's found in build dist:"
        green "\t pgcore: $PGCORE_WHL"
        green "\t pygimli: $PG_WHL"
    fi

    pushd $PROJECT_ROOT
        mkdir -p $PROJECT_DIST

        green "Copying built whl files to project dist: $PROJECT_DIST"

        cp $BUILD_DIR/dist/pgcore*.whl $PROJECT_DIST/
        cp $BUILD_DIR/dist/pygimli*.whl $PROJECT_DIST/

        if [ -d $PROJECT_ROOT/dist-manylinux ]; then
            blue "Copying built whl files to manylinux dist: $PROJECT_ROOT/dist-manylinux"
            cp $BUILD_DIR/dist/pgcore*.whl $PROJECT_ROOT/dist-manylinux/
        fi

        # test editable whl install
        use_venv $VENV_BUILD
        install_WHL_E [build]

        # test whl install
        use_venv $VENV_BUILD'_WHL'
        install_WHL
        deactivate
        rm -rf $VENV_BUILD'_WHL'

    popd
}


function need_build_post(){
    if [ ! -f $PROJECT_DIST/pgcore*.whl ]; then
        yellow "pgcore wheel not found in project dist. Building first. (build_post)"
        build_post
    fi
}


function test_pre(){
    green "*** test_pre (Prepare testing environment) ***"
    need_build_post

    pushd $PROJECT_ROOT
        new_venv $VENV_TEST
        install_WHL_E [test]
    popd
}


function test(){
    green "*** test (Testing now ...) ***"

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
    green "*** doc_pre (Preparing documentation ...) ***"
    need_build_post

    pushd $PROJECT_ROOT
        new_venv $VENV_DOC
        install_WHL_E [doc]
    popd
}


function doc(){
    green "*** doc (Creating documentation) ***"

    use_venv $VENV_DOC
    if python -c 'import sphinx' &>/dev/null; then
        green "sphinx is installed"
    else
        yellow "sphinx is NOT installed (calling doc_pre)"
        doc_pre
    fi

    pushd $PROJECT_ROOT
        use_venv $VENV_DOC

        pushd $BUILD_DIR
            #touch CMakeCache.txt # to renew search for sphinx
            cmake $PROJECT_SRC

            if [ ! -z "$SKIP_GALLERY" ]; then
                NG='-NG'
            fi

            #make clean-gallery # should not be necessary
            if [ -x "$(command -v xvfb-run)" ]; then
                # xvfb is necessary for headless display of pyvista plots
                echo "xvfb-run available: using it to build docs"
                xvfb-run make doc$NG
            else
                echo "xvfb-run not available: building docs without it"
                make doc$NG
            fi
        popd
    popd
}


function doc_post(){
    green "*** doc_post (Deploying html) ***"

    if [ ! -f $BUILD_DIR/doc/_build/html/index.html ]; then
        yellow "Documentation html not found in dist. Building first. (doc)"
        doc
    fi

    use_venv $VENV_DOC
    VERSION=$(python -c 'import pygimli; print(pygimli.__version__.split("(")[0].strip())')
    green "Deploying documentation for pygimli version: $VERSION"

    if [ "$RUNNER_NAME" == "pgserver" ]; then
        # Should only run on pgserver
        green "On pgserver - deploying documentation to /var/www/html"
        rsync -avP --delete $BUILD_DIR/doc/_build/html/ /var/www/html
    else
        yellow "Not on pgserver - skipping documentation deployment"
    fi

    mkdir -p $ARCHIVE_DIR
    pushd $BUILD_DIR/doc/_build
        tar -czf $ARCHIVE_DIR/html-$VERSION.tgz html
    popd
    green "Archived documentation to $ARCHIVE_DIR/html-$VERSION.tgz"
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
    green "*** Install (Creating editable installation in build venv) ***"
    need_build_post

    pushd $PROJECT_ROOT
        use_venv $VENV_PYGIMLI
        install_WHL_E [opt]
        testReport
    popd

    green "Editable installation created in venv: $VENV_PYGIMLI"
    green ""
    green "To use it, call: "
    green ""
    green "source $VENV_PYGIMLI/bin/activate #(linux/macos)"
    green "source $VENV_PYGIMLI/Scripts/activate #(windows)"
}


function twine(){
    WHL_FILE=$1
    green "*** twine upload whl file: $WHL_FILE ***"
    python -m twine upload --repository pypi $WHL_FILE
}


function deploy(){
    target_whl=$1
    if [ -z $target_whl ]; then
        red "give target_whl=\"pygimli\" or \"pgcore\" as argument"
        return
    fi
    green "*** Deploy $target_whl file ***"
    need_build_post
    use_venv $VENV_BUILD

    TARGET_WHL=$(ls $PROJECT_DIST/$target_whl*.whl | head -n 1)

    if [ -z $TARGET_WHL ]; then
        red "Target whl file $target_whl not found in $PROJECT_DIST"
        return
    else
        blue "FOUND: $TARGET_WHL"
    fi
    echo "Do you wish to deploy ? (can't be undone)"
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) twine $TARGET_WHL; break;;
            No ) return;;
        esac
    done
    return
}


function help(){
    echo ""
    echo "run: ${BASH_SOURCE[0]} TARGET"
    echo ""
    echo "TARGETS:"
    echo "    help           show this help"
    echo "    clean          remove build artifacts for PYVERSION"
    echo "    build_pre      prepare build environment venv"
    echo "    build          [build_pre] build project"
    echo "    build_whls     [build] build whl files for pgcore and pygimli"
    echo "    build_post     [build] test build"
    echo "    test_pre       [build] prepare test environment"
    echo "    test           [test_pre] run tests"
    echo "    doc_pre        [build] prepare documentation environment"
    echo "    doc            [doc_pre] build documentation"
    echo "    doc_post       [doc] deploy documentation"
    echo "    install        [build] Create default editable installation in venv"
    echo "    deploy         [build_post] deploy wheel files"
    echo "    all            [clean build test doc]"
    echo ""
    echo "ENVIRONMENT variables:"
    echo "    BASEPYTHON   base python interpreter. Default system python3."
    echo "    PYVERSION    python version (e.g. 3.11, 3.14t) if no BASEPYTHON is given"
    echo "    SOURCE_DIR   source directory (default: top-level directory of the project)"
    echo "    SKIP_GALLERY Skip building the gallery in the documentation"
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
    green "GITHUB action runner on WORKSPACE=$WORKSPACE at $RUNNER_NAME"
    green "RUNNER_ARCH=$RUNNER_ARCH"
    green "RUNNER_TEMP=$RUNNER_TEMP"
    green "GITHUB_REF_NAME=$GITHUB_REF_NAME"
    green "GITHUB_JOB=$GITHUB_JOB"
    green "GITHUB_REPOSITORY=$GITHUB_REPOSITORY"
    green "GITHUB_WORKSPACE=$GITHUB_WORKSPACE"

    WORKSPACE=$GITHUB_WORKSPACE
    JOB_NUM=$GITHUB_RUN_NUMBER
    OS=$RUNNER_OS
    SYSTEM=$(uname -s)

elif [ -z $WORKSPACE ]; then
    WORKSPACE=$(realpath $(pwd))
    if [ -z $OS ]; then
        uname_s=$(uname -s 2>/dev/null || echo)
        blue "Determining OS from uname: $uname_s"
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
    green "Local Build (no Jenkins) on WORKSPACE=$WORKSPACE"
else
    green "Unknown CI Build on WORKSPACE=$WORKSPACE"
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
        red "Could not find project root marker (.git or pyproject.toml). Using script directory name as SOURCE_DIR."
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
        red "Could not determine absolute path for $p"
        return
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
            red "python3 not found in PATH. Please install python3 or set BASEPYTHON to a valid python interpreter."
            red "e.g., BASEPYTHON=../../miniconda3/python bash ${BASH_SOURCE[0]}"
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
ARCHIVE_DIR=$(realpath $WORKSPACE/archive)

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
    deploy)
        deploy $2;;
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
