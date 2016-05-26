CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
IF ("${CMAKE_VERSION}" VERSION_GREATER 3.0.0)
  cmake_policy (SET CMP0046 OLD)
ENDIF()

# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHORNET)
# Set options that need to be enabled BEFORE ANY PROJECT OR LANGUAGE COMMAND
IF (ONHORNET GREATER -1)
  SET(CMAKE_C_COMPILER       cc)
  SET(CMAKE_CXX_COMPILER     CC)
  SET(CMAKE_Fortran_COMPILER ftn)
ENDIF()

PROJECT(Boltzplatz)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
ENABLE_LANGUAGE(Fortran)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# =========================================================================
# Options
# =========================================================================
#include(ExternalProject)
OPTION(BOLTZPLATZ_MPI "Compile SINGLE or MPI version" ON)
SET(BOLTZPLATZ_EQNSYSNAME "maxwell" CACHE STRING "Used equation system")
SET_PROPERTY(CACHE BOLTZPLATZ_EQNSYSNAME PROPERTY STRINGS maxwell maxwell_pois electrostatic electrostatic_pois poisson)
OPTION(BOLTZPLATZ_PARTICLES "Particles" ON)
SET(BOLTZPLATZ_TIMEDISCMETHOD "RK4" CACHE STRING "Time integration method")
SET_PROPERTY(CACHE BOLTZPLATZ_TIMEDISCMETHOD PROPERTY STRINGS Euler-Explicit RK3 RK4 RK14 ParticleEuler-FieldRK4 DSMC RESERVOIR ParticleEuler-FieldEulerImplicit ParticleExplicit-FieldImplicitO3 ParticleExplicit-FieldImplicitO4 Ascher3 ParticleImplicit IPEF3 IPEF4 ImplicitO3 ImplicitO4 FieldRelaxation FieldRelaxation-AdaptiveTimeStep LD LD-DSMC)
OPTION(BOLTZPLATZ_GAUSS_NODES "Gauss nodes" ON)
SET(BOLTZPLATZ_POLYNOMIAL_DEGREE "N" CACHE STRING "N = compile with variable polynomial degree. Use any other number to compile with fixed polynomial degree.")
OPTION(BOLTZPLATZ_HLRS "enable the settings for the HLRS" OFF)
OPTION(BOLTZPLATZ_HDF5F90 "enable old HDF5 interface" OFF)
OPTION(BOLTZPLATZ_MKL "use MKL libraries instead of BLAS/LAPACK" OFF)
OPTION(BOLTZPLATZ_CODE_ANALYZE "Code Analyze" OFF)
OPTION(BOLTZPLATZ_CTAGS "Use Ctags" OFF)
ADD_DEFINITIONS(-DPP_nElems=NELEMZ)
OPTION(BOLTZPLATZ_HDG "enable Hybridized Discontinuous Galerkin (HDG) depending on choice of TIMEDISCMETHOD" OFF)

# =========================================================================
# Machine environements and specific library paths
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHLRS)
IF (ONHLRS GREATER -1)
  ADD_DEFINITIONS(-DLUSTRE -DHLRS)
  MESSAGE(STATUS "Compiling on HLRS system")
  SET(BOLTZPLATZ_HLRS ON)
ENDIF()
#STRING(COMPARE EQUAL ${CMAKEHOST} "juqueen" ONJUQUEEN)
#IF (ONJUQUEEN)
  #ADD_DEFINITIONS(-DGPFS)
#ENDIF()

# =========================================================================
# USERBLOCK + Preproc_flags
# =========================================================================

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    get_cmake_property(_vars CACHE_VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

SET(configuration ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake)
FILE(WRITE ${configuration} "")
GETLISTOFVARSSTARTINGWITH("BOLTZPLATZ_" boltzplatzVars)
FOREACH (_var IN LISTS boltzplatzVars)
   GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
   SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF")) 
   IF(${boolian})
      FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
   ELSE()
      IF(${_var})
         STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
      ENDIF()
      FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
   ENDIF()
ENDFOREACH()

# =========================================================================
# General preproc DEFINES
# =========================================================================
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# Equation system
# =========================================================================

IF(BOLTZPLATZ_EQNSYSNAME STREQUAL "maxwell")
  #SET(BOLTZPLATZ_C_PREPROC ${F2P_C_PREPROC} -Dmaxwell)
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell)
  SET(BOLTZPLATZ_HDG OFF)
ELSEIF(BOLTZPLATZ_EQNSYSNAME STREQUAL "maxwell_pois")
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell_pois -DPOIS)
  SET(BOLTZPLATZ_HDG OFF)
ELSEIF(BOLTZPLATZ_EQNSYSNAME STREQUAL "electrostatic_pois")
  ADD_DEFINITIONS(-DPP_nVar=4 -Delectostratic -DPOIS)
  SET(BOLTZPLATZ_HDG OFF)
ELSEIF(BOLTZPLATZ_EQNSYSNAME STREQUAL "electrostatic")
  ADD_DEFINITIONS(-DPP_nVar=4 -Delectrostatic)
  SET(BOLTZPLATZ_HDG OFF)
ELSEIF(BOLTZPLATZ_EQNSYSNAME STREQUAL "poisson")
  ADD_DEFINITIONS(-DPP_nVar=1 -Dpoisson)
  SET(BOLTZPLATZ_HDG ON)
ELSEIF(BOLTZPLATZ_EQNSYSNAME STREQUAL "magnetostatic")
  ADD_DEFINITIONS(-DPP_nVar=3 -Dmagnetostatic)
  SET(BOLTZPLATZ_HDG ON)
ELSEIF(BOLTZPLATZ_EQNSYSNAME STREQUAL "magnetostatic_poisson")
  ADD_DEFINITIONS(-DPP_nVar=4 -Dmagnetostatic_poisson)
  SET(BOLTZPLATZ_HDG ON)
ELSE()
   MESSAGE(ERROR "Unknown equation system.")
ENDIF()

# =========================================================================
# PARTICLES
# =========================================================================

IF(BOLTZPLATZ_PARTICLES)
  ADD_DEFINITIONS(-DPARTICLES)
ENDIF(BOLTZPLATZ_PARTICLES)

# =========================================================================
# CODE_ANALYZE
# =========================================================================

IF(BOLTZPLATZ_CODE_ANALYZE)
  ADD_DEFINITIONS(-DCODE_ANALYZE)
ENDIF(BOLTZPLATZ_CODE_ANALYZE)

# =========================================================================
# Timedisc method
# =========================================================================
IF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "Euler-Explicit")
   IF(BOLTZPLATZ_HDG)
      ADD_DEFINITIONS(-DPP_HDG -DPP_TimeDiscMethod=500)
   ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=0)
   ENDIF()
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RK3")
   IF(BOLTZPLATZ_HDG)
      ADD_DEFINITIONS(-DPP_HDG -DPP_TimeDiscMethod=501 -DLSERK)
   ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=1 -DLSERK)
   ENDIF()
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RK4")
   IF(BOLTZPLATZ_HDG)
      ADD_DEFINITIONS(-DPP_HDG -DPP_TimeDiscMethod=502 -DLSERK)
   ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=2 -DLSERK)
   ENDIF()
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RK14")
   IF(BOLTZPLATZ_HDG)
      ADD_DEFINITIONS(-DPP_HDG -DPP_TimeDiscMethod=506 -DLSERK)
   ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=6 -DLSERK)
   ENDIF()
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ParticleEuler-FieldRK4")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=5)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "DSMC")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=4)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RESERVOIR")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=42)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ParticleEuler-FieldEulerImplicit")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=100 -DIMEX)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ParticleExplicit-FieldImplicitO3")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=101 -DIMEX)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ParticleExplicit-FieldImplicitO4")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=102 -DIMEX)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "Ascher3")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=103 -DIMEX)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ParticleImplicit")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=110 -DIMPA)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "IPEF3")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=111 -DIMPA)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "IPEF4")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=112 -DIMPA)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ImplicitO3")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=121 -DIMPA)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "ImplicitO4")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=122 -DIMPA)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "FieldRelaxation")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=200)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "FieldRelaxation-AdaptiveTimeStep")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=201)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "LD")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=1000)
ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "LD-DSMC")
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=1001)
ELSE()
   MESSAGE(ERROR "Unknown timedisc method system.")
ENDIF()

# =========================================================================
# HDG
# =========================================================================

IF(BOLTZPLATZ_HDG)
   ADD_DEFINITIONS(-DPP_HDG)
   IF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "Euler-Explicit")
      MESSAGE(STATUS "use HDG with Euler-Explicit")
   ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RK3")
      MESSAGE(STATUS "use HDG with RK3")
   ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RK4")
      MESSAGE(STATUS "use HDG with RK4")
   ELSEIF(BOLTZPLATZ_TIMEDISCMETHOD STREQUAL "RK14")
      MESSAGE(STATUS "use HDG with RK14")
   ELSE()
      MESSAGE(ERROR ": Unknown timedisc method system for HDG.")
   ENDIF()
ENDIF(BOLTZPLATZ_HDG)

# =========================================================================
# Node type
# =========================================================================
IF(BOLTZPLATZ_GAUSS_NODES)
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSE()
  ADD_DEFINITIONS(-DPP_NodeType=2)
ENDIF(BOLTZPLATZ_GAUSS_NODES)

# =========================================================================
# Polynomial degree
# =========================================================================
ADD_DEFINITIONS(-DPP_N=${BOLTZPLATZ_POLYNOMIAL_DEGREE})

# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
   SET (CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: None Debug Release Profile."
      FORCE)
   IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile SANI)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
   ENDIF()
endif (NOT CMAKE_BUILD_TYPE)

# =========================================================================
# boltzplatz sources & targets
# =========================================================================

# include boltzplatz.h 
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/boltzplatz.h DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

# collect all F90 files needed for boltzplatzlib
FILE(GLOB_RECURSE boltzplatzF90 ./src/analyze/*.f90
                                ./src/dg/*.f90
                                ./src/equations/${BOLTZPLATZ_EQNSYSNAME}/*.f90
                                ./src/filter/*.f90
                                ./src/globals/*.f90
                                ./src/interpolation/basis.f90
                                ./src/interpolation/interpolation.f90
                                ./src/interpolation/interpolation_vars.f90
                                ./src/interpolation/changeBasis.f90
                                ./src/interpolation/prolongtoface.f90
                                ./src/io_hdf5/*.f90
                                ./src/linearsolver/*.f90
                                ./src/loadbalance/*.f90
                                ./src/mesh/*.f90
                                ./src/mpi/*.f90
                                ./src/output/*.f90
                                ./src/pml/*.f90
                                ./src/precond/*.f90
                                ./src/readintools/*.f90
                                ./src/recordpoints/*.f90
                                ./src/restart/*.f90
                                ./src/timedisc/*.f90
                                ./src/utils/*.f90 )


FILE(GLOB_RECURSE hdgF90 ./src/hdg/*.f90)
FILE(GLOB mainF90 ./src/boltzplatz.f90)
FILE(GLOB_RECURSE particlesF90 ./src/particles/*.f90
                               ./src/interpolation/eval_xyz.f90)
FILE(GLOB noparticlesF90 ./src/particles/analyze/particle_analyze_vars.f90 
                         ./src/particles/particle_mesh/particle_surfaces_vars.f90)
LIST(APPEND boltzplatzF90 ${mainF90})
IF (BOLTZPLATZ_PARTICLES)
   LIST(APPEND boltzplatzF90 ${particlesF90})
 ELSE()
   LIST(APPEND boltzplatzF90 ${noparticlesF90})
ENDIF(BOLTZPLATZ_PARTICLES)
IF(BOLTZPLATZ_HDG)
  LIST(APPEND boltzplatzF90 ${hdgF90})
ENDIF(BOLTZPLATZ_HDG)

# define boltzplatz library 
ADD_LIBRARY(boltzplatzlibF90 OBJECT ${boltzplatzF90})
SET_TARGET_PROPERTIES(boltzplatzlibF90 PROPERTIES LINKER_LANGUAGE Fortran)

#define boltzplatzlib
ADD_LIBRARY(boltzplatzlib STATIC $<TARGET_OBJECTS:boltzplatzlibF90>)
SET_TARGET_PROPERTIES(boltzplatzlib PROPERTIES LINKER_LANGUAGE Fortran)
TARGET_LINK_LIBRARIES(boltzplatzlib ${LINK_LIBS})


# define boltzplatz executable
ADD_EXECUTABLE(boltzplatz ./src/boltzplatz.f90)
TARGET_LINK_LIBRARIES(boltzplatz boltzplatzlib) 
SET_TARGET_PROPERTIES(boltzplatz PROPERTIES LINKER_LANGUAGE Fortran)
INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/boltzplatz DESTINATION bin)
INSTALL(FILES  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)

# =========================================================================
# FLAGS
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see boltzplatz.h)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# =========================================================================
# COMPILER STUFF
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  # set compiler
  #IF (BOLTZPLATZ_MPI)
  #   SET (CMAKE_Fortran_COMPILER mpif90)
  #ENDIF()
  #environment specific variables
  IF(ONHLRS GREATER -1)
    SET(BOLTZPLATZ_INSTRUCTION "-march=core-avx2")
  ELSE()
    SET(BOLTZPLATZ_INSTRUCTION "-march=native -finline-functions")
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  SET (CMAKE_Fortran_FLAGS_SANI    "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall -fsanitize=address,undefined -fno-omit-frame-pointer -Wc-binding-type")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(boltzplatzlibF90 PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
  SET_TARGET_PROPERTIES(boltzplatz PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   # ifort
   #environment specific variables
   IF(ONHORNET GREATER -1)
     SET(BOLTZPLATZ_INSTRUCTION "-xCORE-AVX2")
   ELSE()
      SET(BOLTZPLATZ_INSTRUCTION "-xHost")
   ENDIF()
   IF (BOLTZPLATZ_MPI AND ONHORNET LESS 0)
      SET(CMAKE_Fortran_COMPILER mpif90)
   ENDIF()
   # set Flags
   SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -assume bscc -shared-intel -DINTEL")
   SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 ${BOLTZPLATZ_INSTRUCTION} -vec-report0 -no-prec-div ")
   # -ipo is not working
   # -fp-model precise -> identical results 85s
   # -fp-model strict  -> identical results 96s
   # -rec-sqrt -fma -no-fma
   # additional -fp-model precise -fp-model exept -fp-model strict
   SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -p -O2 ${BOLTZPLATZ_INSTRUCTION} -vec-report0 -no-prec-div")
   #SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays -init=zero")
   SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O2 -fpe0 -traceback -warn all -fstack-protector -assume protect_parens -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -ftrapuv -init=snan -init=arrays -init=zero")
   # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(boltzplatzlibF90 PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
   SET_TARGET_PROPERTIES(boltzplatz    PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
ELSE()
   MESSAGE(ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(BOLTZPLATZ_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  LIST(GET MPI_INCLUDE_PATH 0 MPI_INCLUDE_PATH)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${MPI_INCLUDE_PATH}/../lib)
  #INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  MESSAGE(STATUS "MPI include dir: " ${MPI_INCLUDE_PATH})
  ADD_DEFINITIONS(-DMPI)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_INCLUDE_PATH}")
  SET (LINK_LIBS ${LINK_LIBS} ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(boltzplatzlib ${MPI_Fortran_LIBRARIES})
  TARGET_LINK_LIBRARIES(boltzplatz       ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LIBRARIES})
  MESSAGE(STATUS "Building boltzplatz with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
ENDIF()

# =========================================================================
# HDF5 library
# =========================================================================
#HDF5
SET(BOLTZPLATZ_HDF5_FOUND FALSE)
#IF(ONHLRS GREATER -1)
  #SET(HDF5_USE_STATIC_LIBRARIES TRUE)
#ENDIF()
FIND_PACKAGE (HDF5 NAMES hdf5 COMPONENTS Fortran static QUIET)
IF(HDF5_FOUND)
  LIST(GET HDF5_INCLUDE_DIR 0 HDF5_INCLUDE_DIR)
  MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR_FORTRAN})
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN})
  SET(CMAKE_Fortran_FLAGS "${CMAKE_FORTRAN_FLAGS} -I${HDF5_INCLUDE_DIR_FORTRAN}")
  TARGET_LINK_LIBRARIES(boltzplatz ${HDF5_FORTRAN_STATIC_LIBRARY} ${HDF5_LIBRARIES})
  MESSAGE(STATUS "Building boltzplatz with parallel HDF5 (HyperSlab)")
  ADD_DEFINITIONS(-DH5TOOLSDIR=\"${HDF5_INCLUDE_DIR}/../bin/\")
  SET(BOLTZPLATZ_HDF5_FOUND TRUE)
ELSE()
  SET(BOLTZPLATZ_HDF5_DIR $ENV{HDF5_DIR} CACHE PATH "Directory of HDF5")
  IF(EXISTS "${BOLTZPLATZ_HDF5_DIR}/include/")
    INCLUDE_DIRECTORIES (${BOLTZPLATZ_HDF5_DIR}/include)
    MESSAGE(STATUS "HDF5 include dir: " ${BOLTZPLATZ_HDF5_DIR}/include)
    SET(BOLTZPLATZ_HDF5_FOUND TRUE)
  ENDIF()
  IF(EXISTS "${BOLTZPLATZ_HDF5_DIR}/include/static")
    INCLUDE_DIRECTORIES (${BOLTZPLATZ_HDF5_DIR}/include/static)
    MESSAGE(STATUS "HDF5 include dir: " ${BOLTZPLATZ_HDF5_DIR}/include/static)
    SET(BOLTZPLATZ_HDF5_FOUND TRUE)
  ENDIF()
  IF(NOT BOLTZPLATZ_HDF5_FOUND)
    MESSAGE(ERROR "No hdf5 installation found!")
  ENDIF()
  FIND_LIBRARY(BOLTZPLATZ_HDF51 NAMES libhdf5_fortran.a hdf5_fortran PATHS  ${BOLTZPLATZ_HDF5_DIR}/lib)
  TARGET_LINK_LIBRARIES(boltzplatz ${BOLTZPLATZ_HDF51})
  IF(EXISTS "${BOLTZPLATZ_HDF5_DIR}/lib/libhdf5_hl_f90cstub.a")
    FIND_LIBRARY(BOLTZPLATZ_HDF52 NAMES libhdf5_f90cstub.a hdf5_f90cstub PATHS ${BOLTZPLATZ_HDF5_DIR}/lib)
    TARGET_LINK_LIBRARIES(boltzplatz ${BOLTZPLATZ_HDF52})
  ENDIF()
  FIND_LIBRARY(BOLTZPLATZ_HDF53 NAMES libhdf5.a hdf5 PATHS ${BOLTZPLATZ_HDF5_DIR}/lib)
  TARGET_LINK_LIBRARIES(boltzplatz ${BOLTZPLATZ_HDF53} -ldl -lz)
  MESSAGE(STATUS "Building boltzplatz with parallel HDF5 (HyperSlab)")
  IF (ONHORNET GREATER -1)
     # customized location of the tools on hazelhen
    ADD_DEFINITIONS(-DH5TOOLSDIR=\"${BOLTZPLATZ_HDF5_DIR}/../../../../hdf5/1.8.14/bin/\")
  ELSE()
    ADD_DEFINITIONS(-DH5TOOLSDIR=\"${BOLTZPLATZ_HDF5_DIR}/bin/\")
  ENDIF()
ENDIF()
IF(NOT BOLTZPLATZ_HDF5_FOUND)
   MESSAGE(ERROR "No hdf5 installation found!")
ENDIF()
IF (BOLTZPLATZ_HDF5F90)
   ADD_DEFINITIONS(-DHDF5_F90)
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
IF (NOT BOLTZPLATZ_MKL)
  # use Lapack/Blas for GNU
  MESSAGE(STATUS "use BLAS/Lapack library")
  FIND_PACKAGE(LAPACK)
  IF (LAPACK_FOUND)
    TARGET_LINK_LIBRARIES(boltzplatz    ${LAPACK_LIBRARIES}) 
    TARGET_LINK_LIBRARIES(boltzplatzlib ${LAPACK_LIBRARIES})
  ENDIF()
ELSE()
  MESSAGE(STATUS "use MKL library")
  LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  FIND_PACKAGE(MKL REQUIRED)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
  LINK_DIRECTORIES(${MKL_LIBRARY})
  TARGET_LINK_LIBRARIES(boltzplatz ${MKL_LIBRARY})
ENDIF()


# =========================================================================
# Userblock
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET boltzplatzlib PRE_BUILD COMMAND touch ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatzlib PRE_BUILD COMMAND rm ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( CMAKE \)]}\" > ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cat ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( GIT REVISIONS \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && git rev-parse --abbrev-ref HEAD >>  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt && git log --oneline -1  >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && echo "master: " >>  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt && git log master --oneline -1  >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( GIT DIFF \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && git diff -p origin/master >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( GIT URL \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && git config --get remote.origin.url >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)

# =========================================================================
# Install
# =========================================================================
INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/boltzplatz DESTINATION bin)
INSTALL(FILES  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)

# =========================================================================
# final message
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: BOLZPLATZ EXECUTABLE GENERATED!")

# =========================================================================
# Ctags
# =========================================================================
IF(BOLTZPLATZ_CTAGS)
  ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ctags -R -f ctags.txt src)
  ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
  ADD_DEPENDENCIES(boltzplatz tags)
ENDIF(BOLTZPLATZ_CTAGS)
