#!/bin/bash
# Script for compiling the code with DEBUG flag in order to count the number of warnings that occur during the compile process
# Run this script in the upper most folder of the PICLas repository.
#
# INPUTARGS: supply "-D" type variables to use cmake instead of ccmake by directly defining the compiler flags
# ./tools/test_max_warnings.sh -DPICLAS_EQNSYSNAME=drift_diffusion -DPICLAS_TIMEDISCMETHOD=Explicit-FV -DLIBS_USE_PETSC=ON
INPUTARGS="$*"

# 0. Get path to the script and location of the CMakeLists.txt file
SCRIPTDIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo -e "$SCRIPTDIR"
CMAKELISTDIR=$(echo "${SCRIPTDIR%/*}")
echo $CMAKELISTDIR

# 1. Create sub-directory
rm -rf build_test_max_warnings
mkdir build_test_max_warnings
cd build_test_max_warnings

# 2. Run cmake
if [[ -n "$INPUTARGS" ]]; then
  cmake tee -DCMAKE_BUILD_TYPE=DEBUG -DPICLAS_BUILD_POSTI=ON -DPOSTI_BUILD_SUPERB=ON $INPUTARGS $CMAKELISTDIR
else
  ccmake tee -DCMAKE_BUILD_TYPE=DEBUG -DPICLAS_BUILD_POSTI=ON -DPOSTI_BUILD_SUPERB=ON $INPUTARGS $CMAKELISTDIR
fi

# 3. Run make
if [[ -f "build.ninja" ]]; then
  ninja -j0 2>&1 | tee compile_output.txt
else
  make -j 2>&1 | tee compile_output.txt
fi

# 4. Extract the number of warning
WARNINGS=$(grep "Warning:" compile_output.txt | wc -l)
echo -e " "
echo -e " "
echo -e " "
echo -e "Found $WARNINGS warnings!"