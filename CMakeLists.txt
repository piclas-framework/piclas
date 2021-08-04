# =========================================================================
# CMake configuration
# =========================================================================
# NOTE: use minimum CMake version required by tools/libs (Paraview, HDF5, FFTW)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1)
# Policies (see https://cmake.org/cmake/help/v3.13/manual/cmake-policies.7.html)
# VERSION_GREATER_EQUAL is incompatible with CMAKE versions < 3.7
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.6.0"))
  SET(GITSHALLOW "GIT_SHALLOW ON")
ENDIF()
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.12.0"))
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()
# Required to link visulib against visuReader
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.13.0"))
  CMAKE_POLICY(SET CMP0079 NEW)
ENDIF()

# =========================================================================
# Project definition
# =========================================================================
PROJECT(PICLas Fortran C CXX)

# Folder for custom CMake source files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)

# =========================================================================
# Performance mode
# =========================================================================
OPTION(PICLAS_PERFORMANCE "Enable performance optimizations (e.g. PGO for GNU)" OFF)

# =========================================================================
# Machine environment
# =========================================================================
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeListsMachine.txt)

# =========================================================================
# Store the current commit information
# =========================================================================
EXECUTE_PROCESS(COMMAND git rev-parse HEAD OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
SET(commit ${CMAKE_CURRENT_SOURCE_DIR}/src/commit.h)
MESSAGE(STATUS "Current git commit ${GIT_COMMIT} will be written to ${commit}")
FILE(WRITE ${commit} "! Define the current commit hash. The default must remain empty, i.e., ''. Do not commit the changed file!\n#define GIT_CURRENT_COMMIT  ''\n")
EXECUTE_PROCESS(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/cmake/setCommitHash.sh ${commit})

# =========================================================================
# Check IPO support:
# =========================================================================
# we need to have languages enabled and compilers defined for this
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.9.0"))
  CMAKE_POLICY(SET CMP0069 NEW)
  INCLUDE(CheckIPOSupported)
  check_ipo_supported(RESULT HASIPO OUTPUT error)
ELSE()
  SET(HASIPO FALSE)
ENDIF()

# =========================================================================
# Output paths
# =========================================================================
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Set custom install dir (needs to be done after project has been defined!)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug, Release, Profile, Sanitize (only GNU), Nitro (only GNU)." FORCE)
   IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile Sanitize Nitro)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
   ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
     SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
   ENDIF()
ENDIF (NOT CMAKE_BUILD_TYPE)

STRING(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE_LC)
IF (BUILD_TYPE_LC MATCHES "debug" OR BUILD_TYPE_LC MATCHES "sanitize")
  ADD_DEFINITIONS(-DUSE_DEBUG=1)
ELSE()
  ADD_DEFINITIONS(-DUSE_DEBUG=0)
  IF (HASIPO)
    # enable IPO globally (IPO branding: Intel => IPO, GNU => LTO, PGI => IPA)
    IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      # Do not use the standard CMake LTO option for GNU (-flto -fno-fat-lto-objects), as it does not allow speed-up during linking
      SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
      # Static libraries require either fat LTO objects (increases compilation time) or the use of linker plugins (per default enabled); the jobserver option reduces linking time
      # More information at: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html#Optimize-Options
      SET(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -flto=jobserver -fuse-linker-plugin")
      # Use the GCC wrapper for the complete toolchain to provide the path to the linker plugin (this might be the problem with using the CMAKE IPO)
      # SET(CMAKE_AR  "gcc-ar")
      # set(CMAKE_Fortran_COMPILER_AR "${CMAKE_AR}")
      # set(CMAKE_Fortran_COMPILER_RANLIB "${CMAKE_RANLIB}")
      # SET(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
      # SET(CMAKE_Fortran_ARCHIVE_FINISH TRUE)
      # SET(CMAKE_NM  "gcc-nm")
      # SET(CMAKE_RANLIB  "gcc-ranlib")
    ELSE()
      SET(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)  # enable IPO globally
    ENDIF()
  ENDIF()
ENDIF()

# =========================================================================
# Location of binary and filenames
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see piclas.h)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# add basedir for location of corresponding userblock-file
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# echo COMPILE_DEFINITIONS
# =========================================================================
GET_DIRECTORY_PROPERTY(comp_defs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
SET(comp_defs_str "COMPILE_DEFINITIONS = ")
FOREACH (def ${comp_defs})
   SET(comp_defs_str "${comp_defs_str} -D${def}")
ENDFOREACH()
ADD_CUSTOM_TARGET(preproc_defines COMMAND echo ${comp_defs_str})

# =========================================================================
# USERBLOCK + Preproc_flags
# =========================================================================

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    GET_CMAKE_PROPERTY(_vars CACHE_VARIABLES)
    STRING(REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    SET(${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# A function to get all user defined variables with a specified prefix
SET(configuration ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake)
FILE(WRITE ${configuration} "")
GETLISTOFVARSSTARTINGWITH("PICLAS_" piclasVars)
GETLISTOFVARSSTARTINGWITH("HDF5_"   piclasVars2)
GETLISTOFVARSSTARTINGWITH("CMAKE_"  piclasVars3)
GETLISTOFVARSSTARTINGWITH("POSTI_"  piclasVars4)
FOREACH (_var IN LISTS piclasVars piclasVars2 piclasVars3 piclasVars4)
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
  # Skip empty variables
  IF(NOT ${${_var}} STREQUAL "")
    SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF"))
    IF(${boolian})
      FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
    ELSE()
      IF(${_var})
        STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
      ENDIF()
      FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

# =========================================================================
# ADD LIBRARIES
# =========================================================================
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/CMakeListsLib.txt)

# =========================================================================
# PICLas
# =========================================================================
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/CMakeLists.txt)

# =========================================================================
# Posti
# =========================================================================
OPTION(PICLAS_BUILD_POSTI "Build POSTI toolset" ON)
IF(PICLAS_BUILD_POSTI)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/posti/CMakeLists.txt)
ELSE()
  GETLISTOFVARSSTARTINGWITH("POSTI_" postiVars)
  FOREACH (_var IN LISTS postiVars)
    UNSET(${_var} CACHE)
  ENDFOREACH()
ENDIF()

# =========================================================================
# Userblock
# =========================================================================
#MESSAGE(STATUS "Running: generateuserblock.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}  ${CMAKE_CACHEFILE_DIR} ${CMAKE_CACHE_MAJOR_VERSION}.${CMAKE_CACHE_MINOR_VERSION}.${CMAKE_CACHE_PATCH_VERSION} ${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals_vars.f90")
MESSAGE(STATUS "Running: generateuserblock.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_VERSION} ${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals_vars.f90")
#ADD_CUSTOM_COMMAND(TARGET libpiclasstatic PRE_BUILD COMMAND
  #sh ${CMAKE_CURRENT_SOURCE_DIR}/tools/userblock/generateuserblock.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}  ${CMAKE_CACHEFILE_DIR} "${CMAKE_CACHE_MAJOR_VERSION}.${CMAKE_CACHE_MINOR_VERSION}.${CMAKE_CACHE_PATCH_VERSION}" "${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals_vars.f90")
ADD_CUSTOM_COMMAND(TARGET libpiclasstatic PRE_BUILD COMMAND
  sh ${CMAKE_CURRENT_SOURCE_DIR}/tools/userblock/generateuserblock.sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_VERSION} "${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals_vars.f90")

# =========================================================================
# generate .piclas in $HOME
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET piclas     POST_BUILD COMMAND echo \"export PICLAS_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}\" > $ENV{HOME}/.piclas)
ADD_CUSTOM_COMMAND(TARGET piclas     POST_BUILD COMMAND echo \"alias piclas='${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/piclas'\" >> $ENV{HOME}/.piclas)

# =========================================================================
# Install
# =========================================================================
INSTALL(FILES    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)
