CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
IF ("${CMAKE_VERSION}" VERSION_GREATER 3.0.0)
  cmake_policy (SET CMP0046 OLD)
ENDIF()
PROJECT(Boltzplatz Fortran)

# todo: MPI is not working properbly, numerical recepies are still missing

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
include(ExternalProject)
SET(BOLTZPLATZ_EQNSYSNAME "maxwell" CACHE STRING "Used equation system")
OPTION(BOLTZPLATZ_MPI "Compile SINGLE or MPI version" ON)
OPTION(BOLTZPLATZ_HDF5F90 "enable old HDF5 interface" OFF)
OPTION(ENABLE_PARTICLES "Particles" OFF)
SET(BOLTZPLATZ_POLYNOMIAL_DEGREE "N" CACHE STRING "N = compile with variable polynomial degree. Use any other number to compile with fixed polynomial degree.")
OPTION(BOLTZPLATZ_MKL "use MKL libraries instead of BLAS/LAPACK" OFF)
ADD_DEFINITIONS(-DPP_nElems=NELEMZ)

# =========================================================================
# Machine environements and specific library paths
# =========================================================================
SET (HDF5DOWNLOAD http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${POSTI_HDF5VERSION}/src/hdf5-${POSTI_HDF5VERSION}.tar.gz )
#SET (TECIODOWNLOAD ${CMAKE_CURRENT_SOURCE_DIR}/share/tecio-${POSTI_TECIOVERSION}.tar.gz)
SET (TECIODOWNLOAD git@129.69.43.151:libs/TECPLOT.git)

SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHORNET)
IF (ONHORNET GREATER -1)
  ADD_DEFINITIONS(-DLUSTRE -DHLRS)
  MESSAGE(STATUS "Compiling on Hazelhen")
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
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DBASEDIR=\"\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"\"")


# =========================================================================
# Equation system
# =========================================================================
IF(BOLTZPLATZ_WHICH_EQNSYS STREQUAL "maxwell")
  #SET(BOLTZPLATZ_C_PREPROC ${F2P_C_PREPROC} -Dmaxwell)
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell)
ELSEIF(BOLTZPLATZ_EQSYSNAME STREQUAL "maxwell_pois")
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell -DPOIS)
ELSEIF(BOLTZPLATZ_EQSYSNAME STREQUAL "electrostatic_pois")
  ADD_DEFINITIONS(-DPP_nVar=4 -Delectostratic -DPOIS)
ELSEIF(BOLTZPLATZ_EQSYSNAME STREQUAL "electrostatic")
  ADD_DEFINITIONS(-DPP_nVar=4 -Delectrostatic)
ELSE()
   MESSAGE(ERROR "Unknown equation system.")
ENDIF()


#STRING(FIND "${BOLTZPLATZ_EQSYSNAME}" "_pois" ISPOIS)
#IF (ISPOIS GREATER -1)
#  SET(BOLTZPLATZ_C_PREPROC ${BOLTZPLATZ_C_PREPROC} -DPOIS)
#  ADD_DEFINITIONS(-DPOIS)
#ENDIF()

# =========================================================================
# PARTICLES
# =========================================================================

IF(ENABLE_PARTICLES)
  #SET(F2P_C_PREPROC ${F2P_C_PREPROC} -DPARTICLE)
  ADD_DEFINITIONS(-DPARTICLE)
ENDIF(ENABLE_PARTICLES)

# =========================================================================
# Node type
# =========================================================================
ADD_DEFINITIONS(-DPP_NodeType=1)


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
                                ./src/equations/${BOLTZPLATZ_EQSYSNAME}/*.f90
                                ./src/filter/*.f90
                                ./src/globals/*.f90
                                ./src/interpolation/*.f90
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


FILE(GLOB mainF90 ./src/*.f90)
FILE(GLOB particlesF90 ./src/particles/*.f90)
LIST(APPEND boltzplatzF90 ${mainF90})
IF (ENABLE_PARTICLES)
   LIST(APPEND boltzplatzF90 ${particlesF90})
ENDIF()


# define boltzplatz library 
ADD_LIBRARY(boltzplatzlibF90 OBJECT ${boltzplatzF90})
SET_TARGET_PROPERTIES(boltzplatzlibF90 PROPERTIES LINKER_LANGUAGE Fortran)
  
#define boltzplatzlib
ADD_LIBRARY(boltzplatzlib STATIC $<TARGET_OBJECTS:boltzplatzlibF90>)
SET_TARGET_PROPERTIES(boltzplatzlib PROPERTIES LINKER_LANGUAGE Fortran)

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
  #environment specific variables
  IF(ONHORNET GREATER -1)
    SET(BOLTZPLATZ_INSTRUCTION "-march=core-avx2")
  ELSE()
    SET(BOLTZPLATZ_INSTRUCTION "-march=native -finline-functions")
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(boltzplatzlibF90 PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
  SET_TARGET_PROPERTIES(boltzplatz PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   # ifort
   #environment specific variables
   IF(ONHORNET)
     SET(BOLTZPLATZ_INSTRUCTION "-xCORE-AVX2")
   ELSE()
     SET(BOLTZPLATZ_INSTRUCTION "-xHost")
   ENDIF()
   # set Flags
   SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -DINTEL")
   SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 ${BOLTZPLATZ_INSTRUCTION} -vec-report0 -no-prec-div")
   SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -p -O2 ${BOLTZPLATZ_INSTRUCTION} -vec-report0 -no-prec-div")
   SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit")
   # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(boltzplatzlibF90 PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
   SET_TARGET_PROPERTIES(boltzplatz    PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
ELSE()
   MESSAGE(ERROR "Unknown compiler")
ENDIF()


# =========================================================================
# HDF5 library
# =========================================================================
#HDF5
SET(BOLTZPLATZ_HDF5_FOUND FALSE)
FIND_PACKAGE (HDF5 NAMES hdf5 COMPONENTS Fortran static QUIET)
IF(HDF5_FOUND)
   LIST(GET HDF5_INCLUDE_DIR 0 HDF5_INCLUDE_DIR)
   MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR_FORTRAN})
   INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN})
   SET(CMAKE_Fortran_FLAGS "${CMAKE_FORTRAN_FLAGS} -I${HDF5_INCLUDE_DIR_FORTRAN}")
   TARGET_LINK_LIBRARIES(boltzplatz ${HDF5_FORTRAN_STATIC_LIBRARY} ${HDF5_LIBRARIES})
   MESSAGE(STATUS "Building boltzplatz with parallel HDF5 (HyperSlab)")
   SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DH5TOOLSDIR=\"\\\"${HDF5_INCLUDE_DIR}/../bin/\\\"\"")
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
   FIND_LIBRARY(BOLTZPLATZ_HDF51 hdf5_fortran  ${BOLTZPLATZ_HDF5_DIR}/lib)
   IF(EXISTS "${BOLTZPLATZ_HDF5_DIR}/lib/libhdf5_hl_f90cstub.a")
     FIND_LIBRARY(BOLTZPLATZ_HDF52 hdf5_f90cstub ${BOLTZPLATZ_HDF5_DIR}/lib)
   ENDIF()
   FIND_LIBRARY(BOLTZPLATZ_HDF53 hdf5          ${BOLTZPLATZ_HDF5_DIR}/lib)
   TARGET_LINK_LIBRARIES(boltzplatz ${BOLTZPLATZ_HDF51} ${BOLTZPLATZ_HDF52} ${BOLTZPLATZ_HDF53} -ldl -lz)
   MESSAGE(STATUS "Building boltzplatz with parallel HDF5 (HyperSlab)")
   SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DH5TOOLSDIR=\"\\\"${BOLTZPLATZ_HDF5_DIR}/bin/\\\"\"")
ENDIF()
IF(NOT BOLTZPLATZ_HDF5_FOUND)
   MESSAGE(ERROR "No hdf5 installation found!")
ENDIF()
IF (BOLTZPLATZ_HDF5F90)
   ADD_DEFINITIONS(-DHDF5_F90)
ENDIF()

# =========================================================================
# TECPLOT library
# =========================================================================
SET(TECDIR "/opt/tecio-2013" CACHE STRING "Directory of tecplot-library")
FIND_LIBRARY(TECLIB tecio ${TECDIR})

# =========================================================================
# Tecplot
# =========================================================================
SET(BOLTZPLATZ_TECIOVERSION "2013" CACHE STRING "TECIO version tag")
SET(BOLTZPLATZ_TECIODIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${BOLTZPLATZ_COMPILERNAME}/tecio-${BOLTZPLATZ_TECIOVERSION})
FIND_PACKAGE(TECIO QUIET)
IF (NOT TECIO_FOUND)
  OPTION(BOLTZPLATZ_BUILDTECIO "Compile and build tecio library" ON)
ELSE()
  OPTION(BOLTZPLATZ_BUILDTECIO "Compile and build tecio library" OFF)
ENDIF()
IF (BOLTZPLATZ_BUILDTECIO)
  IF (NOT EXISTS "${BOLTZPLATZ_TECIODIR}/src/tecio/tecio-${BOLTZPLATZ_TECIOVERSION}/libtecio.a")
    ExternalProject_Add(tecio
      GIT_REPOSITORY ${TECIODOWNLOAD}
      PREFIX ${BOLTZPLATZ_TECIODIR}
      PATCH_COMMAND tar -xvJf tecio-${BOLTZPLATZ_TECIOVERSION}.tar.bz2
      CONFIGURE_COMMAND ""
      BUILD_IN_SOURCE 1
      BUILD_COMMAND cd tecio-${BOLTZPLATZ_TECIOVERSION} && ./Runmake linuxg27x64.24 -tecio 
      INSTALL_COMMAND "" 
    )
  ENDIF()
  SET(TECIO_LIBRARY ${BOLTZPLATZ_TECIODIR}/src/tecio/tecio-${BOLTZPLATZ_TECIOVERSION}/libtecio.a)
  SET(TECIO_LIBRARIES ${BOLTZPLATZ_TECIODIR}/src/tecio/tecio-${BOLTZPLATZ_TECIOVERSION}/libtecio.a -lstdc++)
ENDIF()
#TARGET_LINK_LIBRARIES(boltzplatz ${TECLIB})
#TARGET_LINK_LIBRARIES(boltzplatz stdc++)



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
  LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
  FIND_PACKAGE(MKL REQUIRED)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
  LINK_DIRECTORIES(${MKL_LIBRARY})
  TARGET_LINK_LIBRARIES(boltzplatz ${MKL_LIBRARY})
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(BOLTZPLATZ_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  LIST(GET MPI_INCLUDE_PATH 0 MPI_INCLUDE_PATH)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${MPI_INCLUDE_PATH}/../lib)
  MESSAGE(STATUS "MPI include dir: " ${MPI_INCLUDE_PATH})
  ADD_DEFINITIONS(-DMPI)
  TARGET_LINK_LIBRARIES(boltzplatzlib ${MPI_Fortran_LIBRARIES})
  TARGET_LINK_LIBRARIES(boltzplatz       ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LIBRARIES})
  MESSAGE(STATUS "Building boltzplatz with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
ENDIF()

# =========================================================================
# Userblock
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET boltzplatzlib PRE_BUILD COMMAND touch ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatzlib PRE_BUILD COMMAND rm ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( CMAKE \)]}\" > ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cat ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.svn)
  ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( SVNREVISIONS \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && svn st -q -v src Makefile Makefile.defs CMakeLists.txt >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( SVNDIFF \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && svn diff src Makefile Makefile.defs CMakeLists.txt >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( SVNURL \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && svn info | grep "^URL:" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ELSE()
  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git/svn/.metadata)
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( GIT-SVN REVISIONS \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && git svn log --oneline -1 | cut -d '|' -f1 >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND echo \"{[\( GIT-SVN DIFF \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && git diff -p trunk >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ELSE()
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND
      echo \"{[\( GIT REVISIONS \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET buildeos POST_BUILD COMMAND
      cd ${CMAKE_CURRENT_SOURCE_DIR} && git rev-parse --abbrev-ref HEAD >>  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY2}/userblock.txt &&
      git log --oneline -1  >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY2}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET buildeos POST_BUILD COMMAND
      cd ${CMAKE_CURRENT_SOURCE_DIR} && echo "master: " >>  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY2}/userblock.txt &&
      git log origin/master --oneline -1  >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY2}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND
      echo \"{[\( GIT DIFF \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
    ADD_CUSTOM_COMMAND(TARGET boltzplatz POST_BUILD COMMAND
      cd ${CMAKE_CURRENT_SOURCE_DIR} && git diff -p origin/master >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
  ENDIF()
ENDIF()

# =========================================================================
# Install
# =========================================================================
INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/boltzplatz DESTINATION bin)
INSTALL(FILES  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)
