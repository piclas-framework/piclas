#!/bin/bash
mkdir build_test_warnings
cd build_test_warnings
ccmake tee -DCMAKE_BUILD_TYPE=DEBUG ..
make -j 2>&1 | tee compile_output.txt
WARNINGS=$(grep "Warning:" compile_output.txt | wc -l)
echo -e " "
echo -e " "
echo -e " "
echo -e "Found $WARNINGS warnings!"
