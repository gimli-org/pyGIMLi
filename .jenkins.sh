# Current bug in openblas core detection: https://stackoverflow.com/a/66057286
export OPENBLAS_CORETYPE="ARMV8"

# This script is for continuous integration using Jenkins (http://jenkins-ci.org/)
# It is called from the parent directory, i.e. bash -xe trunk/.jenkins.sh
rm -rf /var/lib/jenkins/.cache/pygimli # We need a better cleaning process

echo "Starting automatic build #$BUILD_NUMBER on" `date`
start=$(date +"%s")

# Show last change to repo in build log
echo `git --git-dir trunk/.git log -1 --pretty="Last change by %cn (%h): %B"`

# link python to python3
ln -sf /usr/bin/python3 python
export PATH=`pwd`:$PATH

# Make polytools available (needed by mt.createCube for example)
export PATH="/var/lib/jenkins/workspace/pyBERT/build/bin/":$PATH

# Show system information
lsb_release -d
uname -a
gcc --version
cmake --version
python --version
python -c "import numpy; print(numpy.__version__)"
python -c "import matplotlib; print(matplotlib.__version__)"

# Check if core was changed
core_update=$(git --git-dir=trunk/.git diff --name-only $GIT_PREVIOUS_COMMIT $GIT_COMMIT | grep -c core/src || true)

# Set this to 1 if you want clean build (also of dependencies)
export CLEAN=0

export GIMLI_NUM_THREADS=4

################
#  Main build  #
################

# just do this if something is wrong with the thirdparty sources
#rm -rf thirdParty/src
#rm -rf build # Uncomment for clean build (expensive, but necessary sometimes)
rm -f build/build_tests.html # remove old test report
#rm -f build/CMakeCache.txt # clean old cache

mkdir -p build
cd build

if [ ! -f CMakeCache.txt ]; then
    # Always rebuild core when Cmake cache does not exist
    core_update=2
fi

if [[ $core_update -ge 1 ]]; then
  echo "# Core changes detected. #"
  touch CMakeCache.txt

  cmake ../trunk -DPYTHON_EXECUTABLE=/usr/bin/python3

  make -j 8 gimli
  make pygimli J=4
  #make gtest # Deactivate temporarily
else
  echo "# No core changes detected. #"
fi

#############################
#  Testing & documentation  #
#############################

# Setup for 3D visualizations
# ------------------------------------------------
export DISPLAY=:99.0
export PYVISTA_OFF_SCREEN=true
# ------------------------------------------------

# Test pygimli
export PYTHONPATH=`pwd`/../trunk:$PYTHONPATH
OMP_THREAD_LIMIT=4 python -c "import pygimli; print(pygimli.Report()); pygimli.test(show=False, abort=True, htmlreport=\"build_tests.html\", devTests=True)"

#Build documentation
#make clean-gallery # TMP for fast Jenkins builds
make doc # = doxygen, sphinxapi, sphinxpdf, sphinxhtml
end=$(date +"%s")
echo "Ending automatic build #$BUILD_NUMBER".
diff=$(($end-$start))
echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."

# If this script fails, a mail with the log file will be send to mail@pygimli.org.
# If it succeeds, the documentation will be uploaded to pygimli.org.
# In any case, a badge icon of the current status and the log file will be uploaded.
