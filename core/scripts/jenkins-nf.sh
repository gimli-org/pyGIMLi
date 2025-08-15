#!/bin/bash
# make the script exit on every fail
set -e

# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called for the newfea branch from the jenkins workspace, i.e.,
# bash -xe ./gimli/core/scripts/jenkins-nf.sh
if [[ $- == *i* ]]; then
    # run with . *.sh
    echo "Running on interactive mode. Aliases will work"
else
    # run with bash *.sh
    echo "Running on NON interactive mode. Aliases will not work. Setting then now."
    shopt -s expand_aliases ## else aliases will be ignored in this bash session
    alias python='python3'
    alias return='exit'
fi

GREEN(){
    echo -e '\033[0;32m'
}
BLUE(){
    echo -e '\033[0;34;49m'
}
NCOL(){
    echo -e '\033[0m'
}
function clean(){
    GREEN
    echo "************************************"
    echo "*** Cleaning ...                 ***"
    echo "************************************"
    NCOL

    pushd $PROJECT_ROOT
        echo "remove PROJECT_BLD: $PROJECT_BLD"
        rm -rf $PROJECT_BLD
        echo "remove BUILD_ENV: $BUILD_ENV"
        rm -rf $BUILD_ENV
        echo "remove TEST_ENV: $TEST_ENV"
        rm -rf $TEST_ENV
    popd
}
function build(){
    GREEN
    echo "************************************"
    echo "*** Building ...                 ***"
    echo "************************************"
    NCOL

    pushd $PROJECT_ROOT

        python -m venv $BUILD_ENV
        source $BUILD_ENV/bin/activate

        echo "Updating pip ..."
        python -m pip install --upgrade pip

        echo "Installing or update build venv ..."
        python -m pip install --quiet -r $PROJECT_SRC/core/scripts/requirements.txt

        mkdir -p $PROJECT_BLD
        pushd $PROJECT_BLD

            if [ $core_update -ge 1 ]; then
                echo "# Core changes detected. #"
                rm -rf CMakeCache.txt
            else
                echo "# No core changes detected. #"
            fi

            cmake $PROJECT_SRC -DTESTVENV=$TEST_ENV
            make -j $GIMLI_NUM_THREADS
            make pygimli J=$GIMLI_NUM_THREADS
            make whlTest

        popd
    popd
}
function test(){
    GREEN
    echo "************************************"
    echo "*** Testing ...                  ***"
    echo "************************************"
    NCOL
    source $TEST_ENV/bin/activate
}
function pylint(){
    source $TEST_ENV/bin/activate
    rm -f pylint.log

    python -m pylint --rcfile=$PROJECT_SRC/.pylintrc \
            --output-format=parseable \
            --msg-template='{path}:{line}: [{msg_id}, {obj}] {msg} ({symbol})' \
            $PROJECT_SRC/$PROJECT | tee pylint.log

    PYLINT_SCORE=$(grep "Your code has been rated at" pylint.log | tr -d -c '0-9/. ' | tr -s ' ' | cut -d '/' -f1 | tr -d "[:blank:]")
    echo "Pylint score: $PYLINT_SCORE/10"
}
function docclean(){
    GREEN
    echo "************************************"
    echo "*** Cleaning old documentation   ***"
    echo "************************************"
    NCOL
}
function doc(){
    GREEN
    echo "************************************"
    echo "*** Creating documentation       ***"
    echo "************************************"
    NCOL
    source $TEST_ENV/bin/activate
}
function deploy(){
    GREEN
    echo "************************************"
    echo "*** Deploying ...                ***"
    echo "************************************"
    NCOL
    PGCOREWHL=`ls -dtr1 $PROJECT_BLD/wheelhouse/pgcore* | tail -1`
    PGWHL=`ls -dtr1 $PROJECT_BLD/wheelhouse/pygimli* | tail -1`

    BLUE
    echo "copy final whl $PGCOREWHL $WHEELHOUSE"
    echo "copy final whl $PGWHL $WHEELHOUSE"
    NCOL

    mkdir -p $WHEELHOUSE
    cp $PGCOREWHL $WHEELHOUSE
    cp $PGWHL $WHEELHOUSE

}
function skip(){
    GREEN
    echo "Skipping continuous integration build and tests."
    NCOL
}
function help(){
    echo ""
    echo run: ${BASH_SOURCE[0]} "build|test|doc|docclean|depclean|all"
    echo ""
}
function all(){
    build
    test
    pylint
    docclean
    doc
    deploy
}
start=$(date +"%s")
echo "Starting automatic build #$BUILD_NUMBER on" `date`

if [ -z $WORKSPACE ]; then
    GREEN
    echo "Local Build (no Jenkins)"
    NCOL
    WORKSPACE=$(pwd)
else
    GREEN
    echo "Jenkins Build"
    NCOL
    echo "JOB_NAME=$JOB_NAME"
    echo "JOB_BASE_NAME=$JOB_BASE_NAME"
    echo "JENKINS_HOME=$JENKINS_HOME"
    echo "WORKSPACE=$WORKSPACE"
    echo "BUILD_TAG=$BUILD_TAG"
fi

PROJECT=gimli

AGENTS_ROOT=$WORKSPACE/..
PROJECT_ROOT=$(realpath $WORKSPACE)
PROJECT_SRC=$PROJECT_ROOT/$PROJECT
WHEELHOUSE=$(realpath $AGENTS_ROOT/wheelhouse) # path for resulting whls

# Show system information
lsb_release -d
uname -a
python --version

if [ -z $PYVERSION ]; then
    PYVERSION=$(python -c 'import sys; print(f"{sys.version_info.major}{sys.version_info.minor}")')
    echo "building for python: $PYVERSION"
else
    echo "building for python: $PYVERSION (forced by setting PYVERSION)"
fi

PYPLATFORM=cp$PYVERSION-cp$PYVERSION
PROJECT_BLD=$PROJECT_ROOT/build_$PYPLATFORM
BUILD_ENV=$PROJECT_ROOT/build_$PYPLATFORM'_venv'
TEST_ENV=$PROJECT_ROOT/test_$PYPLATFORM'_venv'

echo "AGENTS_ROOT=$AGENTS_ROOT"
echo "PROJECT_ROOT=$PROJECT_ROOT"
echo "PROJECT_SRC=$PROJECT_SRC"
echo "PROJECT_BLD=$PROJECT_BLD"
echo "BUILD_ENV=$BUILD_ENV"
echo "TEST_ENV=$TEST_ENV"
echo "WHEELHOUSE=$WHEELHOUSE"

echo "Starting automatic build #$BUILD_NUMBER on" `date` "ROOT: $PROJECT_ROOT"
start=$(date +"%s")

export GIMLI_NUM_THREADS=$((`nproc --all` - 4))

# Check if core was changed
CORE_NEEDS_UPDATE=$(git --git-dir=$PROJECT_SRC/.git diff --name-only $GIT_PREVIOUS_COMMIT $GIT_COMMIT | grep -c core/src || true)
echo "Core update: $CORE_NEEDS_UPDATE"

# Show last change to repo in build log
echo `git --git-dir $PROJECT_SRC/.git log -1 --pretty="Last change by %cn (%h): %B"`
LAST_COMMIT_MSG=`git --git-dir $PROJECT_SRC/.git log -1`

# Find if there is a [CI cmd] in the last commit message
if [[ "$LAST_COMMIT_MSG" == *"[CI"* ]]; then
    CI_CMD=`echo -e $LAST_COMMIT_MSG | sed 's/.*\[CI \([^]]*\)\].*/\1/g'`
    echo "Custom CI command from git command message run before all:" $CI_CMD
    #$CI_CMD
fi

[ $# -lt 1 ] && help

for arg in $@
do
    case $arg in
    build)
        build;;
    test)
        test;;
    doc)
        doc;;
    pylint)
        pylint;;
    docclean)
        docclean;;
    depclean)
        depclean;;
    deploy)
        deploy;;
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
