#!/bin/bash
# Script for compiling the code with DEBUG flag in order to count the number of warnings that occur during the compile process
# Run this script in the upper most folder of the PICLas repository.

# 1. Run script
rm -r build_test_max_warnings || true
mkdir build_test_max_warnings
cd build_test_max_warnings
ccmake tee -DCMAKE_BUILD_TYPE=DEBUG -DPICLAS_BUILD_POSTI=ON -DPOSTI_BUILD_SUPERB=ON ..
make -j 2>&1 | tee compile_output.txt
WARNINGS=$(grep "Warning:" compile_output.txt | wc -l)
echo -e " "
echo -e " "
echo -e " "
echo -e "Found $WARNINGS warnings!"
