#=================================================================================================================================
# Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
# This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
# For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
#
# FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
#=================================================================================================================================

#add_unit_test(<target> <sources>...)
#
# Adds test executable, <target>, built from <sources> and
# adds the test so that CTest will run it.
# Both the executable and the test will be named <target>.
#
function(add_unit_test target)
    # Add executable and libraries
    ADD_EXECUTABLE(${target} ${ARGN})
    TARGET_LINK_LIBRARIES(${target} libpiclasshared)
    # Set properties of target
    SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${PICLAS_COMPILE_FLAGS})
    SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
    SET_TARGET_PROPERTIES(${target}
                          PROPERTIES
                          ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/lib"
                          LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/lib"
                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/unitTests/bin")
    # Add this test to the global list of tests
    SET(UNIT_TEST_TARGETS ${UNIT_TEST_TARGETS} ${target} PARENT_SCOPE)
    # Add the test
    ADD_TEST(NAME ${target}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/unitTests
    COMMAND  ${CMAKE_BINARY_DIR}/unitTests/bin/${target})
endfunction()

# Enable testing in general
ENABLE_TESTING()
# ... and make it parallel
INCLUDE(ProcessorCount)
PROCESSORCOUNT(N)
SET(ENV{CTEST_PARALLEL_LEVEL} "${N}")
MESSAGE(STATUS "Executing unit tests with $ENV{CTEST_PARALLEL_LEVEL} processes")

# Add the unit tests
add_unit_test(ReadInTools   ./unitTests/ReadInTools.f90)
add_unit_test(MatrixInverse ./unitTests/MatrixInverse.f90)
add_unit_test(MortarBasis   ./unitTests/MortarBasis.f90)
#add_unit_test(NodesAndWeights ./unitTests/NodesAndWeights.f90)
#add_unit_test(Vandermonde ./unitTests/Vandermonde.f90)
#add_unit_test(DerivativeMatrix ./unitTests/DerivativeMatrix.f90)
#add_unit_test(ChangeBasis ./unitTests/ChangeBasis.f90)
#add_unit_test(SurfInt ./unitTests/SurfInt.f90)
#add_unit_test(ProlongToFace ./unitTests/ProlongToFace.f90)

# Perform unit tests when the last test has been compiled
ADD_CUSTOM_TARGET( all_tests ALL
                   DEPENDS ${UNIT_TEST_TARGETS})
ADD_CUSTOM_COMMAND(TARGET all_tests
                   COMMENT "Run unit tests"
                   POST_BUILD COMMAND ctest --output-on-failure -j${N}
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
