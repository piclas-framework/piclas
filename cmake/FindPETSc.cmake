#=================================================================================================================================
# Copyright (c) preCICE
# https://github.com/precice/precice
# GNU LESSER GENERAL PUBLIC LICENSE Version 3
#==================================================================================================================================
# This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
#==================================================================================================================================
#
# FindPETSc
# ---------
#
# Locates the PETSc library using pkg-config module PETSc
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following IMPORTED target:
#
#  PETSc::PETSc        - the PETSc library
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
#  PETSC_FOUND          - if false, do not try to link to PETSc
#  PETSC_LIBRARIES      - a list of the full paths to all libraries
#  PETSC_INCLUDE_DIRS   - a list of all include directories
#  PETSC_VERSION        - the full version of PETSc MAJOR.MINOR.PATCH
#  PETSC_VERSION_MAJOR  - the MAJOR part of PETSC_VERSION
#  PETSC_VERSION_MINOR  - the MINOR part of PETSC_VERSION
#  PETSC_VERSION_PATCH  - the PATCH part of PETSC_VERSION
#
# Author: Frédéric Simonis @fsimonis
# Modifications for PICLas: Paul Nizenkov

cmake_policy(VERSION 3.10)

# Generate a argument for cmake pkg-config call (Note: variable uses PETSc spelling due to FIND_PACKAGE call in CMakelistsLib.txt)
# Fall-back option is available if PkgConfig is not available
if(PETSc_FIND_QUIETLY)
  find_package(PkgConfig QUIET)
else()
  find_package(PkgConfig)
endif()

if(PKG_CONFIG_FOUND)
  # pkg-config will look at /usr by default, add additional search paths here
  IF(DEFINED ENV{PETSC_DIR})
    SET(PETSC_DIR "$ENV{PETSC_DIR}/lib/pkgconfig")
    # Add our PETSc to pkg-config search path
    SET(ENV{PKG_CONFIG_PATH} ${PETSC_DIR})
  ENDIF()
  # Build the pkg-config version spec
  set(_pkg_version_spec "")
  if(DEFINED PETSC_FIND_VERSION)
    if(PETSC_FIND_VERSION_EXACT)
      set(_pkg_version_spec "=${PETSC_FIND_VERSION}")
    else()
      set(_pkg_version_spec ">=${PETSC_FIND_VERSION}")
    endif()
  endif()

  # Allow system flags
  set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
  set(ENV{PKG_CONFIG_ALLOW_SYSTEM_LIBS} 1)

  # Use pkg-config to find PETSc
  set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "YES")

  if(PETSc_FIND_QUIETLY)
    pkg_check_modules(PETSC QUIET IMPORTED_TARGET GLOBAL "petsc${_pkg_version_spec}")
  else()
    pkg_check_modules(PETSC       IMPORTED_TARGET GLOBAL "petsc${_pkg_version_spec}")
  endif()

  unset(_pkg_version_spec)

  # Extract version parts from the version information
  if(PETSC_VERSION)
    set(_petsc_versions "")
    string(REGEX MATCHALL "[0-9]+" _petsc_versions ${PETSC_VERSION})
    list(GET _petsc_versions 0 _petsc_version_major)
    list(GET _petsc_versions 1 _petsc_version_minor)
    list(GET _petsc_versions 2 _petsc_version_patch)

    set(PETSC_VERSION ${PETSC_VERSION} CACHE STRING "Full version of PETSc")
    set(PETSC_VERSION_MAJOR ${_petsc_version_major} CACHE INTERNAL "Major version of PETSc")
    set(PETSC_VERSION_MINOR ${_petsc_version_minor} CACHE INTERNAL "Minor version of PETSc")
    set(PETSC_VERSION_PATCH ${_petsc_version_patch} CACHE INTERNAL "Patch version of PETSc")

    unset(_petsc_versions)
    unset(_petsc_version_major)
    unset(_petsc_version_minor)
    unset(_petsc_version_patch)
  endif()
else()
  # Fall-back if pkg-config is not available
  IF(NOT DEFINED PETSC_DIR)
    SET(PETSC_DIR $ENV{PETSC_DIR})
  ENDIF()

  # Set PETSC_FOUND to true if the path is defined
  if(DEFINED PETSC_DIR AND NOT PETSC_DIR STREQUAL "")
    SET(PETSC_FOUND TRUE)
  endif()

  # Check the include folder, the libraries and variables, get the version
  IF(PETSC_FOUND)
    IF(EXISTS "${PETSC_DIR}/include")
      SET(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include")
    ELSE()
      SET(PETSC_FOUND FALSE)
      MESSAGE(SEND_ERROR "Error: PETSc include not found in ${PETSC_DIR}")
    ENDIF()

    IF(EXISTS "${PETSC_DIR}/lib/libpetsc.so")
      SET(PETSC_LIBRARIES "${PETSC_DIR}/lib/libpetsc.so")
      SET(PETSC_LINK_LIBRARIES "${PETSC_DIR}/lib/libpetsc.so")
    ELSEIf(EXISTS "${PETSC_DIR}/lib/libpetsc.a")
      SET(PETSC_LIBRARIES "${PETSC_DIR}/lib/libpetsc.a")
      SET(PETSC_LINK_LIBRARIES "${PETSC_DIR}/lib/libpetsc.a")
    ELSE()
      SET(PETSC_FOUND FALSE)
      MESSAGE(SEND_ERROR "Error: PETSc library not found in ${PETSC_DIR}")
    ENDIF()

    IF(EXISTS ${PETSC_DIR}/conf/petscvariables)
      FILE(STRINGS ${PETSC_DIR}/conf/petscvariables PETSC_VARIABLES NEWLINE_CONSUME)
    ELSEIf(EXISTS ${PETSC_DIR}/lib/petsc/conf/petscvariables)
      FILE(STRINGS ${PETSC_DIR}/lib/petsc/conf/petscvariables PETSC_VARIABLES NEWLINE_CONSUME)
    ELSE()
      SET(PETSC_FOUND FALSE)
      MESSAGE(SEND_ERROR "Error: PETSc variables not found in ${PETSC_DIR}")
    ENDIF()

    # Find "^#define PETSC_VERSION_MAJOR" and get only the numbers and remove trailing line breaks
    SET(PETSC_VERSION_FILE "${PETSC_INCLUDE_DIRS}/petscversion.h")
    IF(EXISTS "${PETSC_VERSION_FILE}")
      EXECUTE_PROCESS(COMMAND cat "${PETSC_VERSION_FILE}" COMMAND grep "^#define PETSC_VERSION_MAJOR" COMMAND grep -o "[[:digit:]]*" COMMAND tr -d '\n' OUTPUT_VARIABLE PETSC_VERSION_MAJOR)
      EXECUTE_PROCESS(COMMAND cat "${PETSC_VERSION_FILE}" COMMAND grep "^#define PETSC_VERSION_MINOR" COMMAND grep -o "[[:digit:]]*" COMMAND tr -d '\n' OUTPUT_VARIABLE PETSC_VERSION_MINOR)
      EXECUTE_PROCESS(COMMAND cat "${PETSC_VERSION_FILE}" COMMAND grep "^#define PETSC_VERSION_SUBMINOR" COMMAND grep -o "[[:digit:]]*" COMMAND tr -d '\n' OUTPUT_VARIABLE PETSC_VERSION_PATCH)
      set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_PATCH}" CACHE STRING "Full version of PETSc")
    ELSE()
      SET(PETSC_FOUND FALSE)
      MESSAGE(SEND_ERROR "Error: PETSc version file not found in ${PETSC_DIR}")
    ENDIF()
  ENDIF()
endif()
unset(_petsc_quiet_arg)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc
  REQUIRED_VARS PETSC_FOUND PETSC_INCLUDE_DIRS PETSC_LIBRARIES PETSC_LINK_LIBRARIES
  VERSION_VAR PETSC_VERSION
  )

mark_as_advanced(PETSC_INCLUDE_DIRS PETSC_LIBRARIES PETSC_VERSION_MAJOR PETSC_VERSION_MINOR PETSC_VERSION_PATCH VERSION_VAR PETSC_VERSION)