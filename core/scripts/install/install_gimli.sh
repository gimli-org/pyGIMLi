#!/usr/bin/env bash

if [[ $0 == *.sh ]]; then
    SCRIPT_REPO=$(dirname $0)
    GET="cat"
    echo "Using local scripts from: $SCRIPT_REPO"
else
    SCRIPT_REPO='-Ls https://raw.githubusercontent.com/gimli-org/gimli/dev/core/scripts/install'
    GET="curl"
    echo "Fetching scripts from: $SCRIPT_REPO"
fi

# Default values
GIMLI_ROOT=$(pwd)/gimli
GIMLI_SOURCE_DIR=$GIMLI_ROOT/gimli
GIMLI_BUILD_DIR=$GIMLI_ROOT/build
PYGIMLI_SOURCE_DIR=$GIMLI_SOURCE_DIR

PARALLEL_BUILD=2
PYTHON_MAJOR=3
UPDATE_ONLY=0
RUN_TESTS=0
BRANCH=''
CLEAN=0

# Utility functions
expand_tilde()
# Better handling of ~ character in path argument
{
    case "$1" in
    (\~)        echo "$HOME";;
    (\~/*)      echo "$HOME/${1#\~/}";;
    (\~[^/]*/*) local user=$(eval echo ${1%%/*})
                echo "$user/${1#*/}";;
    (\~[^/]*)   eval echo ${1};;
    (*)         echo "$1";;
    esac
}

function boxprint()
# Print arguments in colored box
{
  local s="$*"
  tput setaf 3
  echo " -${s//?/-}-
| ${s//?/ } |
| $(tput setaf 4)$s$(tput setaf 3) |
| ${s//?/ } |
 -${s//?/-}-"
  tput sgr 0
}

function is_virtual_env()
# Returns 0 if in virtual environment, 1 otherwise
# See PEP 0668: https://peps.python.org/pep-0668/
{
    case "$(python -c "import sys; \\
        print(sys.base_prefix != sys.prefix or hasattr(sys, 'real_prefix'))")" in
        "True") return 0 ;;
        *)      return 1 ;;
    esac
}

boxprint "Geophysical Inversion and Modelling Library (www.gimli.org)"

help(){
    echo "====================================================================="
    echo "Command line options:"
    echo "---------------------------------------------------------------------"
    echo " add these option with out leading - at the end of the command line"
    echo " e.g.    "
    echo " curl -Ls install.pygimli.org | bash -s py=3 path=./gimli j=8"
    echo "---------------------------------------------------------------------"
    echo "h|help"
    echo "      This help"
    echo "path=GIMLI_ROOT "
    echo "      Path for the gimli root dir. Default path=$GIMLI_ROOT"
    echo "py=PYTHON_MAJOR"
    echo "      Build for Python 2 or 3. Default py=$PYTHON_MAJOR"
    echo "j=PARALLEL_BUILD"
    echo "      Number of parallel compiler jobs. Default j=$PARALLEL_BUILD"
    echo "u|update"
    echo "      Just update your gimli installation. "
    echo "      The build path will not be removed in the first."
    echo "      This may work or may not work .. please use at own risk"
    echo "b|branch=branch"
    echo "      Checkout with a given git branch name. Default=''"
    echo "t|test"
    echo "      Run full testsite. Note 3d tests might fail without tetgen installation Default=False''"
    exit
}

for i in "$@"; do
    case "$i" in
        *gimli_root=*|*path=*)
            GIMLI_ROOT=$(expand_tilde `echo $i | sed 's/[-a-zA-Z0-9]*=//'`)
        ;;
        *j=*)
            PARALLEL_BUILD=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
        ;;
        *py=*)
            PYTHON_MAJOR=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
        ;;
        *u*|*update*)
            UPDATE_ONLY=1
        ;;
        *b*|*branch*)
            BRANCH=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
        ;;
        *C*|*clean*)
            CLEAN=1
        ;;
        *t*|*test*)
            RUN_TESTS=1
        ;;
        *h*|*help*)
            help
        ;;
    esac
done

# echo "=========================================="
if [ $(uname) == "Darwin" ]; then
    SYSTEM='mac'
    echo "Determining system ... DARWIN system found."
    # need to be tested
    CMAKE_GENERATOR='Unix Makefiles'
elif [ $(uname -o) == "Msys" ]; then
    if [ $(uname -m) == "x86_64" ]; then
        echo "Determining system ... Msys WIN64 system found."
        SYSTEM='win64'
        # just run this for the initial call .. after a mingw-shell restart these setting is automatic
        export PATH=$PATH:/mingw64/bin
    else
        echo "Determining system ... Msys WIN32 system found."
        echo "Warning!! WIN32 system not longer supported."
        SYSTEM='win32'
        # just run this for the initial call .. after a mingw-shell restart these setting is automatic
        export PATH=$PATH:/mingw32/bin
    fi
    # check for winpython
    "$GET" $SCRIPT_REPO/install_$SYSTEM'_winpython.sh' | bash
    [ -f $GIMLI_ROOT/.bash_hint_python ] && source $GIMLI_ROOT/.bash_hint_python

    # for mingw build
    CMAKE_GENERATOR='Unix Makefiles'

elif [ $(uname -o) == "GNU/Linux" ]; then
    echo "Determining system ... LINUX system found."
    SYSTEM='linux'
    CMAKE_GENERATOR='Unix Makefiles'
fi
echo "------------------------------------------"

if [ $(uname) == "Darwin" ]; then
    # mac lacks readlink -f
    export GIMLI_ROOT=$(python -c "import os; print(os.path.realpath('"${GIMLI_ROOT//\'/\\\'}"'))")
else
    # python way creates wrong leading /c/ for mysys
    export GIMLI_ROOT=$(readlink -f $GIMLI_ROOT)
fi

export UPDATE_ONLY=$UPDATE_ONLY
export RUN_TESTS=$RUN_TESTS
PYTHON_MAJOR=`python -c 'import sys; print(sys.version_info.major)'`
PYTHON_MINOR=`python -c 'import sys; print(sys.version_info.minor)'`
PYVERSION=$PYTHON_MAJOR'.'$PYTHON_MINOR

echo "Installing at: "$GIMLI_ROOT
echo "Build for Python="$PYVERSION
echo "Parallelize with j="$PARALLEL_BUILD
echo "Update only: " $UPDATE_ONLY
echo "Run tests: " $RUN_TESTS
echo "Branch: " $BRANCH
echo "CLEAN build: " $CLEAN

getGIMLI(){
    mkdir -p $GIMLI_ROOT
    pushd $GIMLI_ROOT

        if [ -d $GIMLI_SOURCE_DIR ]; then
            pushd $GIMLI_SOURCE_DIR
                git pull
            popd
        else
            git clone https://github.com/gimli-org/gimli.git
        fi

        if [ -n "$BRANCH" ]; then
            pushd $GIMLI_SOURCE_DIR
                echo "Switching to branch: " $BRANCH
                git checkout $BRANCH
            popd
        fi

        chmod +x $PYGIMLI_SOURCE_DIR/apps/*

        if is_virtual_env; then
            echo "Installing requirements into virtual enviroment"
            pip3 install -r $PYGIMLI_SOURCE_DIR/dev_requirements.txt
        else
            echo "Installing requirements into user space"
            pip3 install -r $PYGIMLI_SOURCE_DIR/dev_requirements.txt --user
        fi

    popd
}

buildGIMLI(){
    pushd $GIMLI_ROOT

        [ $UPDATE_ONLY -eq 0 ] && rm -rf $GIMLI_BUILD_DIR

        mkdir -p $GIMLI_BUILD_DIR

        pushd $GIMLI_BUILD_DIR
            echo "Installing for python=$PYVERSION"
            cmake -G "$CMAKE_GENERATOR" $GIMLI_SOURCE_DIR -DPYVERSION=$PYVERSION          

            make -j$PARALLEL_BUILD && make pygimli J=$PARALLEL_BUILD
        popd
    popd
}

testGIMLI(){
    echo ""
    export PYTHONPATH=$PYGIMLI_SOURCE_DIR
    python -c 'import pygimli as pg; print("pygimli version:", pg.__version__)'

    if [ $RUN_TESTS -eq 1 ]; then
        if [ -x "$(command -v pytest)" ]; then
            python -c 'import pygimli as pg; pg.test()'
        else
            echo "no pytest found"
        fi
    fi
    echo "export PYTHONPATH=$PYGIMLI_SOURCE_DIR" > $GIMLI_ROOT/.bash_hint_pygimli
    echo "export PATH=$PYGIMLI_SOURCE_DIR/apps:\$PATH" >> $GIMLI_ROOT/.bash_hint_pygimli
}

echo -e "\n========================================================================="
echo "Installing system prerequisites for:" $SYSTEM

"$GET" $SCRIPT_REPO/install_$SYSTEM'_prereqs.sh' | bash


if [ "$GET" == "curl" ]; then
    echo -e "\n========================================================================="
    echo "Get GIMLi sources"
    echo "-------------------------------------------------------------------------"
    getGIMLI
else
    echo -e "\n========================================================================="
    echo "local call .. skipp fetching sources"
    echo "-------------------------------------------------------------------------"
fi

echo -e "\n========================================================================="
echo "BUILD: pyGIMLi for" $SYSTEM
echo "-------------------------------------------------------------------------"
buildGIMLI

echo -e "\n========================================================================="
echo "TEST: pyGIMLi installation                                                "
echo "-------------------------------------------------------------------------"
testGIMLI

echo -e "\n========================================================================="
echo "Finialize:                                                               "
echo "-------------------------------------------------------------------------"
echo "Set the following setting to use pygimli, either locally"
echo "per session or permanently in your $HOME/.bashrc"
echo "-------------------------------------------------------------------------"
echo ""
[ -f $GIMLI_ROOT/.bash_hint_python ] && cat $GIMLI_ROOT/.bash_hint_python
[ -f $GIMLI_ROOT/.bash_hint_pygimli ] && cat $GIMLI_ROOT/.bash_hint_pygimli
echo ""

#https://raw.githubusercontent.com/gimli-org/gimli/master/scripts/install/install**
