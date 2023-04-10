# Unit tests

Unit tests are used to test individual key units of the source code. Currently these key routines include:

* Functionality of Read-In tools.
* Functionality of matrix inversion via external library

## Integration of unit test with CTest

These unit tests are integrated into the **PICLas** build process using the [CTest](https://cmake.org/Wiki/CMake/Testing_With_CTest) tool. Usually CTest will be run every time you
build **PICLas** and give you an overview on the exit status of each test that looks something like this:

~~~~
Test project /home/piclas/build_
    Start 1: ReadInTools
1/2 Test #1: ReadInTools ......................   Passed    0.12 sec
    Start 2: MatrixInverse
2/2 Test #2: MatrixInverse ....................   Passed    0.12 sec

100% tests passed, 0 tests failed out of 2

Total Test time (real) =   0.24 sec
~~~~

To manually run the tests after a build use the CTest command

~~~~
ctest
~~~~

in your build directory. The manual page of CTest can give you an overview of all available options.

If you don't want to run the test after each build there is a CMake option called PICLAS_UNITTESTS that can be used to turn the tests on and off.
This is an advanced option that CCMake will only show if you enter the advanced mode by pressing the t key. 

## Implementation of unit tests

All unit tests are implemented in FORTRAN and can be found in the subdirectory unitTests in your **PICLas** directory alongside a seperate CMakeLists.txt and some binary input and reference files.

### CMakeLists.txt

The CMakeLists.txt defines a custom function called add_unit_test which can be used in the CMakeLists.txt to add a single test to the CTest tool. The syntax is

~~~~
add_unit_test(NAME SOURCEFILE.F90)
~~~~

All tests are defined using this function. At the end of the CMakeLists.txt a custom target all_tests is defined which includes all unit tests and will run the ctest command after it has been build.

The whole CMakeLists.txt content is included in the main CMakeLists.txt if the option PICLAS_UNITTESTS is set to ON (default) by CMake.

### General unit test structure

The general structure of the unit tests is the same in all cases. They are implemented as FORTRAN programs. The unit test will call a function or subroutine from the **PICLas** framework with input either set in the program itself or read from a binary file.
The output of this call will then be compared to some precomputed reference results (also stored as binary files) with a certain tolerance to account for differences in e.g. compiler versions and system architecture. If the results are within the given tolerance,
the test will be passed, otherwise it will fail by returning a value other than 0.

The programs usually also contain a command line option that can be uses to generate the reference solution from a code version that is known to work correctly. 

Have a look at the source code of one of the already implemented unit tests if you want to have a more detailed idea about how to implement your own tests.

### Generation of reference mesh data

Some of the unit tests require parts of the mesh data structure to be able to call the functions to be tested. For this purpose, a curved single element is created and all the mesh data stored as a binary file called ``UnittestElementData.bin``. This binary file can then be read during runtime
by the unit test programs. 

To generate the curved single element mesh, run **HOPR** with the parameter file provided in the ``unitTest`` subdirectory of **PICLas**. To generate the binary file, run **PICLas** with the following command line argument and the parameter file 
provided in the ``unitTest`` subdirectory:


~~~~
piclas --generateUnittestReferenceData parameter.ini
~~~~
