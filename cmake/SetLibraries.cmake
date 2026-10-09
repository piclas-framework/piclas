# =========================================================================
# Set download locations depending on git origin
# =========================================================================
SET(LIBS_DLPATH "https://piclas.boltzplatz.eu/piclas/")
# Origin pointing to PICLas gitlab
IF("${GIT_ORIGIN}" MATCHES "piclas.boltzplatz.eu" AND "${GIT_ORIGIN}" MATCHES "^git@")
  SET(LIBS_DLPATH "git@piclas.boltzplatz.eu:piclas/")
ENDIF()
MESSAGE (STATUS "Setting download locations depending on git origin [LIBS_DLPATH] ${LIBS_DLPATH}")

# Unset leftover variables from previous runs
UNSET(linkedlibs CACHE)

# =========================================================================
# MPI
# =========================================================================
# Try to find system MPI
SET(MPI_DETERMINE_LIBRARY_VERSION TRUE)
FIND_PACKAGE(MPI QUIET)
IF (MPI_FOUND)
  MESSAGE (STATUS "[MPI] found in system libraries")
  OPTION(LIBS_USE_MPI "Compile SINGLE or MPI version" ON)
ELSE()
  MESSAGE (STATUS "[MPI] not found in system libraries")
  OPTION(LIBS_USE_MPI "Compile SINGLE or MPI version" OFF)
ENDIF()

IF(LIBS_USE_MPI)
  # If library is specifically requested, it is required
  SET(MPI_DETERMINE_LIBRARY_VERSION TRUE)
  FIND_PACKAGE(MPI REQUIRED)

  IF (NOT MPI_Fortran_NO_INTERROGATE)
    FOREACH(DIR ${MPI_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    FOREACH(DIR ${MPI_Fortran_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    LIST(APPEND linkedlibs ${MPI_Fortran_LIBRARIES})
  ENDIF()

  MARK_AS_ADVANCED(FORCE MPI_LIBRARY MPI_EXTRA_LIBRARY)
  REMOVE_DEFINITIONS(-DLIBS_MPICH_FIX_SHM_INTERFACE)

  # Detect MPI implementation and version since it changes some MPI definitions
  IF(MPI_C_LIBRARY_VERSION_STRING MATCHES ".*CRAY MPICH.*" AND MPI_C_VERSION_MAJOR VERSION_EQUAL "3")
    SET(LIBS_MPI_NAME "Cray MPICH")
    STRING(REGEX MATCH "([0-9]+)\\.([0-9]+)" MPI_C_LIBRARY_VERSION ${MPI_C_LIBRARY_VERSION_STRING})
    # Cray MPICH in combination with GNU has problems with calling the same MPI routine
    # with different arguments in the same compilation unit
    IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      SET(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    ENDIF()
  ELSEIF(MPI_C_LIBRARY_VERSION_STRING MATCHES ".*MPICH.*" AND MPI_C_VERSION_MAJOR VERSION_GREATER_EQUAL "3")
    SET(LIBS_MPI_NAME "MPICH")
    STRING(REGEX MATCH "([0-9]+)\\.([0-9]+)" MPI_C_LIBRARY_VERSION ${MPI_C_LIBRARY_VERSION_STRING})
    # Missing interface added in 4.2, see https://github.com/pmodels/mpich/pull/6727
    IF(${MPI_C_LIBRARY_VERSION} VERSION_LESS_EQUAL "4.1")
      ADD_COMPILE_DEFINITIONS(LIBS_MPICH_FIX_SHM_INTERFACE=1)
    ENDIF()
    # MPICH in combination with GNU has problems with calling the same MPI routine
    # with different arguments in the same compilation unit
    IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      SET(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    ENDIF()
  ELSEIF(MPI_C_LIBRARY_VERSION_STRING MATCHES ".*Open MPI.*" AND MPI_C_VERSION_MAJOR VERSION_EQUAL "3")
    SET(LIBS_MPI_NAME "OpenMPI")
    STRING(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" MPI_C_LIBRARY_VERSION ${MPI_C_LIBRARY_VERSION_STRING})
  ELSEIF(MPI_C_LIBRARY_VERSION_STRING MATCHES ".*HPE MPT.*" AND MPI_C_VERSION_MAJOR VERSION_EQUAL "3")
    SET(LIBS_MPI_NAME "HPE MPT")
    STRING(REGEX MATCH "([0-9]+)\\.([0-9]+)" MPI_C_LIBRARY_VERSION ${MPI_C_LIBRARY_VERSION_STRING})
    # HPE MPT is deprecated and should not be supported any more
    # MESSAGE(FATAL_ERROR "HPE MPT not supported any more")
  ELSEIF(MPI_C_LIBRARY_VERSION_STRING MATCHES ".*Intel.*" AND MPI_C_VERSION_MAJOR VERSION_EQUAL "3")
    SET(LIBS_MPI_NAME "Intel MPI")
    STRING(REGEX MATCH "([0-9]+)\\.([0-9]+)" MPI_C_LIBRARY_VERSION ${MPI_C_LIBRARY_VERSION_STRING})
  ELSE()
    MESSAGE(FATAL_ERROR "Cannot detect supported MPI type or version. Valid options are Cray MPICH, IntelMPI, MPICH, HPE MPT, and OpenMPI supporting MPI version 3.x")
    # HPE MPT is deprecated and should not be supported any more
    # MESSAGE(FATAL_ERROR "Cannot detect supported MPI type or version. Valid options are Cray MPICH, IntelMPI, MPICH, and OpenMPI supporting MPI version 3.x")
  ENDIF()

  MESSAGE(STATUS "Compiling with [${LIBS_MPI_NAME}] (${BoldBlue}v${MPI_C_LIBRARY_VERSION}${ColourReset})")
  ADD_COMPILE_DEFINITIONS(USE_MPI=1)

  # LUMI needs even more help here
  IF("${CMAKE_FQDN_HOST}" MATCHES ".can")
    SET(MPI_C_COMPILER       cc)
    SET(MPI_CXX_COMPILER     CC)
    SET(MPI_Fortran_COMPILER ftn)
  ENDIF()
ELSE()
  ADD_COMPILE_DEFINITIONS(USE_MPI=0)
ENDIF()


# =========================================================================
# Add the libraries
# =========================================================================
# Set directory to compile external libraries
IF(LIBS_USE_MPI)
  SET(LIBS_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}-MPI)
ELSE()
  SET(LIBS_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID})
ENDIF()
MARK_AS_ADVANCED(FORCE LIBS_EXTERNAL_LIB_DIR)


# =========================================================================
# HDF5 library
# =========================================================================
# Try to find system HDF5 using CMake
SET(LIBS_HDF5_CMAKE TRUE)

# Set preferences for HDF5 library
SET(HDF5_USE_STATIC_LIBRARIES TRUE)
IF (LIBS_USE_MPI)
  SET(HDF5_PREFER_PARALLEL TRUE)
  FIND_PROGRAM(HDF5_COMPILER h5pcc)
ELSE()
  SET(HDF5_PREFER_PARALLEL FALSE)
  FIND_PROGRAM(HDF5_COMPILER h5cc)
ENDIF()

# When using the configure version, CMake takes the directory of the first HDF5 compiler found
# > h5cc  - serial   version
# > h5pcc - parallel version
# > Thus, we need to prepend the PATH to ensure we are picking the correct one first
IF(NOT "${HDF5_COMPILER}" STREQUAL "" AND NOT "${HDF5_COMPILER}" STREQUAL "HDF5_COMPILER-NOTFOUND")
  SET(ORIGINAL_PATH_ENV "$ENV{PATH}")
  GET_FILENAME_COMPONENT(HDF5_PARENT_DIR ${HDF5_COMPILER} DIRECTORY)
  SET(ENV{PATH} "${HDF5_PARENT_DIR}:$ENV{PATH}")
ENDIF()

IF (NOT LIBS_BUILD_HDF5)
  # ParaView requires the HL libs but we cannot change the search later
  SET(HDF5_COMPONENTS C Fortran)
  FIND_PACKAGE(HDF5 QUIET COMPONENTS ${HDF5_COMPONENTS})

  # Could not find the static version, look for the shared library
  IF(NOT HDF5_FOUND)
    UNSET(HDF5_USE_STATIC_LIBRARIES)
    # ParaView requires the HL libs but we cannot change the search later
    FIND_PACKAGE(HDF5 QUIET COMPONENTS ${HDF5_COMPONENTS})
  ENDIF()

  IF (HDF5_FOUND)
    MESSAGE (STATUS "[HDF5] found in system libraries [${HDF5_DIR}]")
    SET(LIBS_BUILD_HDF5 OFF CACHE BOOL "Compile and build HDF5 library")
  ELSE()
    MESSAGE (STATUS "[HDF5] not found in system libraries")
    SET(LIBS_BUILD_HDF5 ON  CACHE BOOL "Compile and build HDF5 library")
  ENDIF()
ENDIF()

# Use system HDF5
IF(NOT LIBS_BUILD_HDF5)
  # Unset leftover paths from old CMake runs
  UNSET(HDF5_VERSION CACHE)
  UNSET(HDF5_DEFINITIONS)
  UNSET(HDF5_LIBRARIES)
  UNSET(HDF5_INCLUDE_DIR_FORTRAN)
  UNSET(HDF5_INCLUDE_DIR)
  UNSET(HDF5_DIFF_EXECUTABLE)

  # If library is specifically requested, it is required
  # > ParaView requires the HL libs but we cannot change the search later
  FIND_PACKAGE(HDF5 REQUIRED COMPONENTS ${HDF5_COMPONENTS})

  # If fortran module files cannot be found in the HDF5_INCLUDE_DIR set by FIND_PACKAGE(HDF5), obtain the correct path from the target properties (supposedly HDF5_INCLUDE_DIR_FORTRAN)
  # > NOTE: Depending on HDF5 config (flag HDF5_INSTALL_MOD_FORTRAN) and version, mod-files can be located in include/ or mod/, or subdirectories
  FIND_FILE(PATH_MODFILES h5a.mod ${HDF5_INCLUDE_DIR})
  IF(NOT PATH_MODFILES)
    IF(HDF5_USE_STATIC_LIBRARIES)
      GET_TARGET_PROPERTY(HDF5_INCLUDE_DIR hdf5_fortran-static INTERFACE_INCLUDE_DIRECTORIES)
    ELSE()
      GET_TARGET_PROPERTY(HDF5_INCLUDE_DIR hdf5_fortran-shared INTERFACE_INCLUDE_DIRECTORIES)
    ENDIF()
  ENDIF()

  # Check if HDF5 is parallel
  # > HDF5_IS_PARALLEL is set by FIND_PACKAGE(HDF5)
  IF(LIBS_USE_MPI)
    IF(NOT HDF5_IS_PARALLEL)
      MESSAGE(FATAL_ERROR "HDF5 is not built with parallel support. Please install a parallel version of HDF5 or build it yourself.")
    ENDIF()

    # If HDF5 Fortran library is not set, get it from the Fortran target
    IF("${HDF5_Fortran_LIBRARY_hdf5_fortran}" STREQUAL "")
      GET_PROPERTY(HDF5_Fortran_LIBRARY_hdf5_fortran TARGET hdf5::hdf5_fortran PROPERTY LOCATION)
    ENDIF()

    IF(NOT "${HDF5_Fortran_LIBRARY_hdf5_fortran}" STREQUAL "")
      IF(APPLE)
        EXECUTE_PROCESS(COMMAND nm -gU      ${HDF5_Fortran_LIBRARY_hdf5_fortran} COMMAND grep mpio_f08 OUTPUT_VARIABLE HDF5_USES_MPIF08 RESULT_VARIABLE GREP_RESULT OUTPUT_STRIP_TRAILING_WHITESPACE)
      ELSE()
        EXECUTE_PROCESS(COMMAND readelf -Ws ${HDF5_Fortran_LIBRARY_hdf5_fortran} COMMAND grep mpio_f08 OUTPUT_VARIABLE HDF5_USES_MPIF08 RESULT_VARIABLE GREP_RESULT OUTPUT_STRIP_TRAILING_WHITESPACE)
      ENDIF()
    # Cray might still not provide anything, set mpi_f08 to false
    ELSE()
      SET(GREP_RESULT 1)
    ENDIF()

    IF(GREP_RESULT EQUAL 0)
      SET(HDF5_MPI_VERSION "[mpi_f08]")
      SET(HDF5_HAS_MPIF08 TRUE)
    ELSE()
      SET(HDF5_MPI_VERSION "[mpi]")
      SET(HDF5_HAS_MPIF08 FALSE)
    ENDIF()
  ENDIF()

  # If HDF5 Fortran library is not set, get it from the Fortran target
  UNSET(GREP_RESULT)
  IF("${HDF5_C_LIBRARY_hdf5_c}" STREQUAL "")
    GET_PROPERTY(HDF5_C_LIBRARY_hdf5_c TARGET hdf5::hdf5 PROPERTY LOCATION)
  ENDIF()

  # Set build status to system
  SET(HDF5_BUILD_STATUS "system")
ELSE()
  MESSAGE(STATUS "Setting [HDF5] to self-build")
  # Origin pointing to Github
  IF("${GIT_ORIGIN}" STREQUAL "")
    # Use https://github.com/HDFGroup/hdf5.git when using piclas without git
    SET (HDF5DOWNLOAD "https://github.com/HDFGroup/hdf5.git")
  ELSEIF("${GIT_ORIGIN}" MATCHES ".github.com")
    # Use https://github.com/HDFGroup/hdf5.git when using piclas from github.com
    SET (HDF5DOWNLOAD "https://github.com/HDFGroup/hdf5.git")
  ELSE()
    IF("${GIT_ORIGIN}" MATCHES "https://gitlab-ci-token")
      # Use https://piclas.boltzplatz.eu/piclas/hdf5.git when gitlab runner needs to download libs
      STRING(REGEX REPLACE "/piclas.git" "/hdf5.git" HDF5DOWNLOAD ${GIT_ORIGIN})
    ELSE()
      # Use https://piclas.boltzplatz.eu/piclas/hdf5.git when using piclas from piclas.boltzplatz.eu
      SET (HDF5DOWNLOAD ${LIBS_DLPATH}hdf5.git )
    ENDIF()
  ENDIF()
  SET(HDF5_DOWNLOAD ${HDF5DOWNLOAD} CACHE STRING "HDF5 Download-link")
  MESSAGE(STATUS "Setting [HDF5] download link: ${HDF5DOWNLOAD}")
  MARK_AS_ADVANCED(FORCE HDF5_DOWNLOAD)

  # Set HDF5 tag / version
  SET(HDF5_STR "2.2.0")
  SET(HDF5_TAG "${HDF5_STR}" CACHE STRING "HDF5 version tag")
  MARK_AS_ADVANCED(FORCE HDF5_TAG)
  MESSAGE(STATUS "Setting [HDF5] download tag: ${BoldBlue}${HDF5_TAG}${ColourReset}")

  # Set HDF5 build dir
  SET(LIBS_HDF5_DIR ${LIBS_EXTERNAL_LIB_DIR}/HDF5)

  # Check if HDF5 was already built
  UNSET(HDF5_FOUND)
  UNSET(HDF5_VERSION)
  UNSET(HDF5_INCLUDE_DIR)
  UNSET(HDF5_LIBRARIES)
  UNSET(HDF5_Fortran_LIBRARIES)
  FIND_PACKAGE(HDF5 ${HDF5_STR} QUIET COMPONENTS C Fortran HDF5_PREFER_PARALLEL=${LIBS_USE_MPI} PATHS ${LIBS_HDF5_DIR} NO_DEFAULT_PATH)

  # CMake does not correctly pick-up HDF5 parallel support, thus set it manually
  SET(HDF5_IS_PARALLEL ${LIBS_USE_MPI})
  IF(HDF5_IS_PARALLEL)
    SET(HDF5_MPI_VERSION "[mpi_f08]")
    SET(HDF5_HAS_MPIF08 TRUE)
  ENDIF()

  IF(HDF5_FOUND)
    # If re-running CMake, it might wrongly pick-up the system HDF5
    IF(NOT EXISTS ${LIBS_HDF5_DIR}/lib/libhdf5.a)
      UNSET(HDF5_FOUND)
      SET(HDF5_VERSION     ${HDF5_STR})
    ENDIF()

    # CMake might fail to set the HDF5 paths
    IF(HDF5_FOUND AND "${HDF5_LIBRARIES}" STREQUAL "")
      SET(HDF5_INCLUDE_DIR       ${LIBS_HDF5_DIR}/include)
      # WARNING: The order of the following libraries matters! They need to be listed from the most dependent to the least dependent.
    SET(HDF5_LIBRARIES         ${LIBS_HDF5_DIR}/lib/libhdf5_hl.a ${LIBS_HDF5_DIR}/lib/libhdf5.a ${LIBS_HDF5_DIR}/lib/libhdf5_tools.a)
    SET(HDF5_Fortran_LIBRARIES ${LIBS_HDF5_DIR}/lib/libhdf5_fortran.a ${LIBS_HDF5_DIR}/lib/libhdf5_f90cstub.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl_fortran.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl_f90cstub.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl.a ${LIBS_HDF5_DIR}/lib/libhdf5.a ${LIBS_HDF5_DIR}/lib/libhdf5_tools.a)
    ENDIF()
  ENDIF()

  # Check again if HDF5 was found
  IF(NOT HDF5_FOUND)
    # Set parallel build with maximum number of threads
    INCLUDE(ProcessorCount)
    PROCESSORCOUNT(N)

    # Find zlib. We require the dynamically linked version, since on some systems
    # (e.g. Ubuntu 26.04) zlib is built without PIC.
    # ZLIB_USE_STATIC_LIBS = OFF by default
    FIND_PACKAGE(ZLIB REQUIRED)

    # Let CMake take care of download, configure and build
    EXTERNALPROJECT_ADD(HDF5
      GIT_REPOSITORY     ${HDF5_DOWNLOAD}
      GIT_TAG            ${HDF5_TAG}
      GIT_PROGRESS       TRUE
      GIT_SHALLOW        ON
      PREFIX             ${LIBS_HDF5_DIR}
      UPDATE_COMMAND     ""
      # HDF5 explicitely needs "make" to configure
      CMAKE_GENERATOR    "Unix Makefiles"
      BUILD_COMMAND      make -j${N}
      # Set the CMake arguments for HDF5
      CMAKE_ARGS         -DCMAKE_BUILD_TYPE=None -DCMAKE_INSTALL_PREFIX=${LIBS_HDF5_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DHDF5_INSTALL_LIB_DIR=lib -DHDF5_INSTALL_CMAKE_DIR=lib/cmake/hdf5 -DCMAKE_POLICY_DEFAULT_CMP0175=OLD -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DZLIB_INCLUDE_DIR=${ZLIB_INCLUDE_DIRS} -DZLIB_LIBRARY=${ZLIB_LIBRARIES} -DHDF5_ENABLE_SZIP_SUPPORT=OFF -DHDF5_ENABLE_PARALLEL=${LIBS_USE_MPI}
      # Set the build byproducts
      # WARNING: The order of the following libraries matters! They need to be listed from the most dependent to the least dependent.
      BUILD_BYPRODUCTS ${LIBS_HDF5_DIR}/lib/libhdf5_fortran.a ${LIBS_HDF5_DIR}/lib/libhdf5_f90cstub.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl_fortran.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl_f90cstub.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl.a ${LIBS_HDF5_DIR}/lib/libhdf5.a ${LIBS_HDF5_DIR}/lib/libhdf5_tools.a ${LIBS_HDF5_DIR}/bin/h5diff
    )

    # Add CMake HDF5 to the list of self-built externals
    LIST(PREPEND SELFBUILTEXTERNALS HDF5)

    # Set HDF5 version and MPI support
    SET(HDF5_VERSION ${HDF5_STR})

    # Set HDF5 paths
    # > NOTE: For self-built HDF5, we use a specific version, of which we know the installation directory of the fortran module files
    SET(HDF5_INCLUDE_DIR       ${LIBS_HDF5_DIR}/include)
    SET(HDF5_DIFF_EXECUTABLE   ${LIBS_HDF5_DIR}/bin/h5diff)
    # WARNING: The order of the following libraries matters! They need to be listed from the most dependent to the least dependent.
    SET(HDF5_LIBRARIES         ${LIBS_HDF5_DIR}/lib/libhdf5_hl.a ${LIBS_HDF5_DIR}/lib/libhdf5.a ${LIBS_HDF5_DIR}/lib/libhdf5_tools.a)
    SET(HDF5_Fortran_LIBRARIES ${LIBS_HDF5_DIR}/lib/libhdf5_fortran.a ${LIBS_HDF5_DIR}/lib/libhdf5_f90cstub.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl_fortran.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl_f90cstub.a ${LIBS_HDF5_DIR}/lib/libhdf5_hl.a ${LIBS_HDF5_DIR}/lib/libhdf5.a ${LIBS_HDF5_DIR}/lib/libhdf5_tools.a)
  ENDIF()

  # Set build status to self-built
  SET(HDF5_BUILD_STATUS "self-built")
ENDIF()

# HDF5 1.14 references build directory
# > https://github.com/HDFGroup/hdf5/issues/2422
IF(HDF5_VERSION VERSION_EQUAL "1.14")
  LIST(FILTER HDF5_INCLUDE_DIR EXCLUDE REGEX "src/H5FDsubfiling")
ENDIF()

# Actually add the HDF5 paths (system/self-built) to the linking paths, including the library containing dlopen/dlclose (usually -ldl on UNIX machines)
# > INFO: We could also use the HDF5::HDF5/hdf5::hdf5/hdf5::hdf5_fortran targets here but they are not set before compiling self-built HDF5
INCLUDE_DIRECTORIES(BEFORE ${HDF5_INCLUDE_DIR})
LIST(PREPEND linkedlibs ${HDF5_Fortran_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS})
IF(${HDF5_IS_PARALLEL})
  MESSAGE(STATUS "Compiling with ${HDF5_BUILD_STATUS} [HDF5] (${BoldBlue}v${HDF5_VERSION}${ColourReset}) with parallel support ${HDF5_MPI_VERSION}")
ELSE()
  MESSAGE(STATUS "Compiling with ${HDF5_BUILD_STATUS} [HDF5] (${BoldBlue}v${HDF5_VERSION}${ColourReset}) without parallel support")
ENDIF()

# Hide all the HDF5 libs paths
MARK_AS_ADVANCED(FORCE HDF5_DIR)
MARK_AS_ADVANCED(FORCE HDF5_C_INCLUDE_DIR)
MARK_AS_ADVANCED(FORCE HDF5_DIFF_EXECUTABLE)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_INCLUDE_DIR)
MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_dl)
MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_hdf5)
MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_m)
MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_sz)
MARK_AS_ADVANCED(FORCE HDF5_C_LIBRARY_z)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_dl)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_hdf5)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_hdf5_fortran)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_m)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_sz)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_z)
MARK_AS_ADVANCED(FORCE HDF5_hdf5_LIBRARY_hdf5)
MARK_AS_ADVANCED(FORCE HDF5_hdf5_LIBRARY_RELEASE)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_hdf5_fortran)
MARK_AS_ADVANCED(FORCE HDF5_Fortran_LIBRARY_hdf5_fortran_RELEASE)

# Restore the original PATH - only if it has been modified above, i.e. ORIGINAL_PATH_ENV is actually defined
IF(NOT "${HDF5_COMPILER}" STREQUAL "" AND NOT "${HDF5_COMPILER}" STREQUAL "HDF5_COMPILER-NOTFOUND")
  SET(ENV{PATH} "${ORIGINAL_PATH_ENV}")
ENDIF()

# =========================================================================
# Math library
# =========================================================================
# Try to find system LAPACK/OpenBLAS
IF (NOT LIBS_BUILD_MATH_LIB)
  FIND_PACKAGE(LAPACK QUIET)
ENDIF()

IF (LAPACK_FOUND)
  MESSAGE (STATUS "[BLAS/Lapack] found in system libraries")
  SET(LIBS_BUILD_MATH_LIB OFF CACHE BOOL "Compile and build math library")
ELSE()
  MESSAGE (STATUS "[BLAS/Lapack] not found in system libraries")
  SET(LIBS_BUILD_MATH_LIB ON  CACHE BOOL "Compile and build math library")
ENDIF()

# Use system LAPACK/MKL
IF(NOT LIBS_BUILD_MATH_LIB)
  # If library is specifically requested, it is required
  FIND_PACKAGE(LAPACK REQUIRED)
  IF (LAPACK_FOUND)
    LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
    MESSAGE(STATUS "Compiling with system [BLAS/Lapack]")
  ENDIF()

  # VDM inverse, replace lapack with analytical solution
  IF (CMAKE_BUILD_TYPE MATCHES "Debug" AND CMAKE_FQDN_HOST MATCHES "hawk\.hww\.hlrs\.de$")
    MESSAGE(STATUS "Compiling PICLas in debug mode on Hawk with system math lib. Setting VDM inverse to analytical solution")
    ADD_COMPILE_DEFINITIONS(VDM_ANALYTICAL)
  ENDIF()

# Build LAPACK/OpenBLAS in PICLAS
ELSE()
  # Offer LAPACK and OpenBLAS
  SET (LIBS_BUILD_MATH_LIB_VENDOR OpenBLAS CACHE STRING "Choose the type of math lib vendor, options are: LAPACK, OpenBLAS.")
  SET_PROPERTY(CACHE LIBS_BUILD_MATH_LIB_VENDOR PROPERTY STRINGS LAPACK OpenBLAS)

  # Build LAPACK
  IF (LIBS_BUILD_MATH_LIB_VENDOR MATCHES "LAPACK")
    # Origin pointing to Github
    IF("${GIT_ORIGIN}" STREQUAL "")
      # Use https://github.com/Reference-LAPACK/lapack.git when using piclas without git
      SET (MATHLIB_DOWNLOAD "https://github.com/Reference-LAPACK/lapack.git")
    ELSEIF("${GIT_ORIGIN}" MATCHES ".github.com")
      # Use https://github.com/Reference-LAPACK/lapack.git when using piclas from github.com
      SET (MATHLIB_DOWNLOAD "https://github.com/Reference-LAPACK/lapack.git")
    ELSE()
      IF("${GIT_ORIGIN}" MATCHES "https://gitlab-ci-token")
        # Use https://piclas.boltzplatz.eu/piclas/lapack.git when gitlab runner needs to download libs
        STRING(REGEX REPLACE "/piclas.git" "/lapack.git" MATHLIB_DOWNLOAD ${GIT_ORIGIN})
      ELSE()
        # Use https://piclas.boltzplatz.eu/piclas/lapack.git when using piclas from piclas.boltzplatz.eu
        SET (MATHLIB_DOWNLOAD ${LIBS_DLPATH}lapack.git )
      ENDIF()
    ENDIF()
    MESSAGE(STATUS "Downloading from ${MATHLIB_DOWNLOAD}")
    SET (MATH_LIB_DOWNLOAD ${MATHLIB_DOWNLOAD} CACHE STRING "LAPACK Download-link" FORCE)
    SET (MATH_LIB_TAG "v3.12.1")
    MARK_AS_ADVANCED(FORCE MATH_LIB_DOWNLOAD)
    MARK_AS_ADVANCED(FORCE MATH_LIB_TAG)
  # Build OpenBLAS
  ELSEIF (LIBS_BUILD_MATH_LIB_VENDOR MATCHES "OpenBLAS")
    IF("${GIT_ORIGIN}" STREQUAL "")
      # Use https://github.com/xianyi/OpenBLAS when using piclas without git
      SET (MATHLIB_DOWNLOAD "https://github.com/OpenMathLib/OpenBLAS.git")
    ELSEIF("${GIT_ORIGIN}" MATCHES ".github.com")
      # Use https://github.com/xianyi/OpenBLAS when using piclas from github.com
      SET (MATHLIB_DOWNLOAD "https://github.com/OpenMathLib/OpenBLAS.git")
    ELSE()
      IF("${GIT_ORIGIN}" MATCHES "https://gitlab-ci-token")
        # Use https://piclas.boltzplatz.eu/piclas/OpenBLAS.git when gitlab runner needs to download libs
        STRING(REGEX REPLACE "/piclas.git" "/OpenBLAS.git" MATHLIB_DOWNLOAD ${GIT_ORIGIN})
      ELSE()
        # Use https://piclas.boltzplatz.eu/piclas/OpenBLAS.git when using piclas from piclas.boltzplatz.eu
        SET (MATHLIB_DOWNLOAD ${LIBS_DLPATH}OpenBLAS.git )
      ENDIF()
    ENDIF()
    MESSAGE(STATUS "Downloading from ${MATHLIB_DOWNLOAD}")
    SET (MATH_LIB_DOWNLOAD ${MATHLIB_DOWNLOAD} CACHE STRING "OpenBLAS Download-link" FORCE)
    SET (MATH_LIB_TAG "v0.3.29")
    MARK_AS_ADVANCED(FORCE MATH_LIB_DOWNLOAD)
    MARK_AS_ADVANCED(FORCE MATH_LIB_TAG)
  # Unknown math lib vendor
  ELSE()
    MESSAGE(FATAL_ERROR "Unknown math lib vendor")
  ENDIF()

  # Set math libs build dir
  SET(LIBS_MATH_DIR  ${LIBS_EXTERNAL_LIB_DIR}/${LIBS_BUILD_MATH_LIB_VENDOR})

  IF (LIBS_BUILD_MATH_LIB_VENDOR MATCHES "LAPACK")
    # Check if math lib was already built
    IF (NOT EXISTS "${LIBS_MATH_DIR}/lib/liblapack.so")
      # Let CMake take care of download, configure and build
      EXTERNALPROJECT_ADD(${LIBS_BUILD_MATH_LIB_VENDOR}
        GIT_REPOSITORY ${MATH_LIB_DOWNLOAD}
        GIT_TAG ${MATH_LIB_TAG}
        GIT_PROGRESS TRUE
        GIT_SHALLOW ON
        PREFIX ${LIBS_MATH_DIR}
        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_PREFIX=${LIBS_MATH_DIR} -DBLAS++=OFF -DLAPACK++=OFF -DBUILD_SHARED_LIBS=ON -DCBLAS=OFF -DLAPACKE=OFF -DBUILD_TESTING=OFF
        BUILD_BYPRODUCTS ${LIBS_MATH_DIR}/lib/liblapack.so ${LIBS_MATH_DIR}/lib/libblas.so
      )

      LIST(APPEND SELFBUILTEXTERNALS ${LIBS_BUILD_MATH_LIB_VENDOR})
    ENDIF()
  ELSEIF (LIBS_BUILD_MATH_LIB_VENDOR MATCHES "OpenBLAS")
    # Check if math lib was already built
    IF (NOT EXISTS "${LIBS_MATH_DIR}/libopenblas.so")
      # Let CMake take care of download, configure and build
      EXTERNALPROJECT_ADD(${LIBS_BUILD_MATH_LIB_VENDOR}
        GIT_REPOSITORY ${MATH_LIB_DOWNLOAD}
        GIT_TAG ${MATH_LIB_TAG}
        GIT_PROGRESS TRUE
        GIT_SHALLOW ON
        PREFIX ${LIBS_MATH_DIR}
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_BYPRODUCTS ${LIBS_MATH_DIR}/src/${LIBS_BUILD_MATH_LIB_VENDOR}/libopenblas.so
        BUILD_IN_SOURCE TRUE
        INSTALL_COMMAND ""
      )

      LIST(APPEND SELFBUILTEXTERNALS ${LIBS_BUILD_MATH_LIB_VENDOR})
    ENDIF()
  ENDIF()

  IF (LIBS_BUILD_MATH_LIB_VENDOR MATCHES "LAPACK")
    # Set math lib paths
    UNSET(MATH_LIB_LIBRARIES)
    SET(MATH_LIB_LIBRARIES              ${LIBS_MATH_DIR}/lib)

    UNSET(LAPACK_LIBRARY)
    UNSET(BLAS_LIBRARY)
    UNSET(LAPACK_LIBRARIES)

    SET(LAPACK_LIBRARY                  ${MATH_LIB_LIBRARIES}/liblapack.so)
    SET(BLAS_LIBRARY                    ${MATH_LIB_LIBRARIES}/libblas.so)
    SET(LAPACK_LIBRARIES                ${LAPACK_LIBRARY}${BLAS_LIBRARY})

    # Actually add the math lib paths to the linking paths
    INCLUDE_DIRECTORIES (${MATH_LIB_LIBRARIES})
    LIST(APPEND linkedlibs ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
    MESSAGE(STATUS "Compiling with self-built [LAPACK]")
  ELSEIF (LIBS_BUILD_MATH_LIB_VENDOR MATCHES "OpenBLAS")
    # Set math lib paths
    SET(MATH_LIB_LIBRARIES              ${LIBS_MATH_DIR}/src/${LIBS_BUILD_MATH_LIB_VENDOR})

    UNSET(LAPACK_LIBRARY)
    UNSET(LAPACK_LIBRARIES)

    SET(LAPACK_LIBRARY                  ${MATH_LIB_LIBRARIES}/libopenblas.so)
    SET(LAPACK_LIBRARIES                ${LAPACK_LIBRARY}${BLAS_LIBRARY})

    # Actually add the math lib paths to the linking paths
    INCLUDE_DIRECTORIES (${MATH_LIB_LIBRARIES})
    LIST(APPEND linkedlibs ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
    MESSAGE(STATUS "Compiling with self-built [OpenBLAS]")
  ENDIF()
ENDIF()


# =========================================================================
# HOPR pre-processor
# =========================================================================
IF(LIBS_BUILD_HOPR)
  SET (HOPRDOWNLOAD "https://github.com/hopr-framework/hopr.git")
  SET (HOPR_DOWNLOAD ${HOPRDOWNLOAD} CACHE STRING "HOPR Download-link")
  MESSAGE (STATUS "HOPR download link: ${HOPRDOWNLOAD}")
  MARK_AS_ADVANCED(FORCE HOPR_DOWNLOAD)

  IF (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/bin/hopr")
    SET(HOPR_FOUND FALSE)
    MESSAGE(STATUS "HOPR executable not found in: ${CMAKE_CURRENT_BINARY_DIR}. HOPR will be built and/or installed.")
  ELSE()
    SET(HOPR_FOUND TRUE)
    MESSAGE(STATUS "HOPR has already been built.")
  ENDIF()
  # Set HOPR build dir
  SET(LIBS_HOPR_DIR  ${LIBS_EXTERNAL_LIB_DIR}/hopr/build)
  MESSAGE(STATUS "HOPR directory: ${LIBS_HOPR_DIR}")
  # Check if HOPR was already built and the executable is in the piclas/build/bin folder
  IF (NOT HOPR_FOUND)
    SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)
    # Let CMake take care of download, configure and build
    EXTERNALPROJECT_ADD(HOPR
    GIT_REPOSITORY ${HOPR_DOWNLOAD}
    GIT_TAG ${HOPR_TAG}
    GIT_PROGRESS FALSE
    GIT_SHALLOW ON
    PREFIX ${LIBS_HOPR_DIR}
    # Avoids rebuilding during PICLas recompilation if HOPR has already been built in share folder
    UPDATE_DISCONNECTED true
    # Avoid output
    LOG_CONFIGURE ON
    LOG_INSTALL ON
    LOG_BUILD ON
    LOG_OUTPUT_ON_FAILURE ON
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR} -DLIBS_BUILD_HDF5=OFF -DLIBS_BUILD_CGNS=ON -DLIBS_USE_MPI=OFF -DCMAKE_BUILD_TYPE=Release
    # BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/bin/hopr
    )
    # HDF5 has to be built before HOPR
    IF(LIBS_BUILD_HDF5)
      ExternalProject_Add_StepDependencies(HOPR build HDF5)
    ENDIF()
    # Make sure that the hopr binary is in the piclas/build/bin/ folder after piclas was build
    ExternalProject_Add_StepDependencies(HOPR install piclas)
    SET(HOPR_FOUND TRUE)
  ENDIF()
ENDIF()

# Download the HOPR executable directly and place it in the bin folder
IF(LIBS_DOWNLOAD_HOPR)
  MESSAGE(STATUS "HOPR executable download from: ${HOPR_DOWNLOAD_URL}")
  IF(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/bin/hopr")
    MESSAGE(STATUS "HOPR executable not found in: ${CMAKE_CURRENT_BINARY_DIR}/bin/hopr. HOPR will be downloaded/extracted.")
    # Check if .zip file has already been downloaded
    FILE(GLOB HOPR_ZIP_FILE ${CMAKE_CURRENT_BINARY_DIR}/hopr-*.zip)
    IF(NOT EXISTS "${HOPR_ZIP_FILE}")
      EXECUTE_PROCESS(COMMAND curl -L -O ${HOPR_DOWNLOAD_URL})
      FILE(GLOB HOPR_ZIP_FILE ${CMAKE_CURRENT_BINARY_DIR}/hopr-*.zip)
    ENDIF()
    # Check if AppImage has already been extracted
    FILE(GLOB HOPR_APP_FILE ${CMAKE_CURRENT_BINARY_DIR}/hopr-*.AppImage)
    IF(NOT EXISTS "${HOPR_APP_FILE}")
      IF(NOT EXISTS "${HOPR_ZIP_FILE}")
        MESSAGE(WARNING "FAILED to download hopr .zip file from ${HOPR_DOWNLOAD_URL}")
      ELSE()
        EXECUTE_PROCESS(COMMAND unzip ${HOPR_ZIP_FILE})
        FILE(GLOB HOPR_APP_FILE ${CMAKE_CURRENT_BINARY_DIR}/hopr-*.AppImage)
      ENDIF()
    ENDIF()
    # Check if AppImage file exists and create symbolic link into ./bin
    IF(NOT EXISTS "${HOPR_APP_FILE}")
      MESSAGE(WARNING "FAILED to extract hopr .zip file")
    ELSE()
      EXECUTE_PROCESS(COMMAND ln -s ${HOPR_APP_FILE} bin/hopr)
    ENDIF()
  ELSE()
    MESSAGE(STATUS "HOPR executable already exists under ${CMAKE_CURRENT_BINARY_DIR}/bin/hopr")
  ENDIF()
ENDIF()


# # =========================================================================
# # PAPI library
# # =========================================================================
# OPTION(LIBS_USE_PAPI "Use PAPI library to perform performance measurements (e.g. flop counts)." OFF)
# IF(LIBS_USE_PAPI)
#   # If library is specifically requested, it is required
#   FIND_PACKAGE(PAPI REQUIRED)
#   ADD_COMPILE_DEFINITIONS(PAPI)
#   LIST(APPEND linkedlibs ${PAPI_LIBRARIES})
#   INCLUDE_DIRECTORIES(${PAPI_INCLUDE_DIRS})
#   MESSAGE(STATUS "Compiling with [PAPI] benchmark support.")
# ENDIF()


# # =========================================================================
# # OPENMP library
# # =========================================================================
# # Try to find system OpenMP
# IF (NOT LIBS_USE_OPENMP)
#   FIND_PACKAGE(OpenMP QUIET)
# ENDIF()
# IF (OpenMP_FOUND)
#   MESSAGE (STATUS "[OpenMP] found in system libraries")
#   OPTION(LIBS_USE_OPENMP "Enable OpenMP" ON)
# ELSE()
#   MESSAGE (STATUS "[OpenMP] not found in system libraries")
#   OPTION(LIBS_USE_OPENMP "Enable OpenMP" OFF)
# ENDIF()
#
# IF(LIBS_USE_OPENMP)
#   IF ("${CMAKE_VERSION}" VERSION_LESS 3.1.0)
#     MESSAGE(WARNING "For finding OpenMP Fortran flags at least CMake version 3.1.0 is required. Please specify flags manually or use newer CMake version.")
#   ENDIF()
#   # If library is specifically requested, it is required
#   FIND_PACKAGE(OpenMP REQUIRED)
#   SET (CMAKE_Fortran_FLAGS_DEBUG          "${CMAKE_Fortran_FLAGS_DEBUG}          ${OpenMP_Fortran_FLAGS}")
#   SET (CMAKE_Fortran_FLAGS_RELEASE        "${CMAKE_Fortran_FLAGS_RELEASE}        ${OpenMP_Fortran_FLAGS}")
#   SET (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} ${OpenMP_Fortran_FLAGS}")
#   SET (CMAKE_CXX_FLAGS_DEBUG              "${CMAKE_CXX_FLAGS_DEBUG}              ${OpenMP_CXX_FLAGS}")
#   SET (CMAKE_CXX_FLAGS_RELEASE            "${CMAKE_CXX_FLAGS_RELEASE}            ${OpenMP_CXX_FLAGS}")
#   SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}     ${OpenMP_CXX_FLAGS}")
#   SET (CMAKE_EXE_LINKER_FLAGS             "${CMAKE_EXE_LINKER_FLAGS}             ${OpenMP_EXE_LINKER_FLAGS}")
#   ADD_COMPILE_DEFINITIONS(USE_OPENMP=1)
#   MESSAGE(STATUS "Compiling with [OpenMP] (${BoldBlue}v${OpenMP_Fortran_VERSION}${ColourReset})")
# ELSE()
#   ADD_DEFINITIONS(-DUSE_OPENMP=0)
# #  ENDIF()

#=========================================================================
# PETSc
#=========================================================================
OPTION(LIBS_USE_PETSC "Use the PETSc library" OFF)

# Look quietly for PETSc
FIND_PACKAGE(PETSc QUIET)
# Set the default value whether to build PETSc or use the system library
IF(PETSC_FOUND)
  CMAKE_DEPENDENT_OPTION(LIBS_BUILD_PETSC "Compile and build the PETSc library" OFF LIBS_USE_PETSC OFF)
ELSE()
  CMAKE_DEPENDENT_OPTION(LIBS_BUILD_PETSC "Compile and build the PETSc library" ON LIBS_USE_PETSC OFF)
ENDIF()

IF(LIBS_USE_PETSC)
  IF (LIBS_BUILD_PETSC)
    SET(LIBS_BUILD_PETSC_VERSION "3.25.1" CACHE STRING "PETSc self-built version tag")
    MARK_AS_ADVANCED(CLEAR LIBS_BUILD_PETSC_VERSION)
  ELSE()
    UNSET(LIBS_BUILD_PETSC_VERSION CACHE)
    MARK_AS_ADVANCED(FORCE LIBS_BUILD_PETSC_VERSION)
  ENDIF()

  IF (LIBS_BUILD_PETSC)
    # Set PETSc build dir
    SET(LIBS_PETSC_DIR      ${LIBS_EXTERNAL_LIB_DIR}/PETSc/build)
    SET(PETSC_BUILD_DIR     ${LIBS_EXTERNAL_LIB_DIR}/PETSc/linux-gnu-c)
    SET(PETSC_LIB_LIBRARIES ${LIBS_EXTERNAL_LIB_DIR}/PETSc/linux-gnu-c/lib)
    SET(PETSC_INCLUDE_DIRS  ${LIBS_EXTERNAL_LIB_DIR}/PETSc/linux-gnu-c/include)

    # Determine if PETSc should be built with debugging or optimization flags
    IF(CMAKE_BUILD_TYPE MATCHES "Debug" OR CMAKE_BUILD_TYPE MATCHES "Sanitize")
      SET(PETSC_COMPILER_FLAGS --with-debugging=1)
      SET(PETSC_OPTIMIZATION "")
      MESSAGE (STATUS "Compiling [PETSc] with debugging flags (${PETSC_COMPILER_FLAGS})!")
    ELSE()
      SET(PETSC_COMPILER_FLAGS --with-debugging=0)
      SET(PETSC_OPTIMIZATION "-O3 ${PICLAS_INSTRUCTION}")
      MESSAGE (STATUS "Compiling [PETSc] with optimization flags (${PETSC_COMPILER_FLAGS} ${PETSC_OPTIMIZATION})!")
    ENDIF()

    # Cofigure problem for cmake >=4
    # https://stackoverflow.com/questions/79534856/cannot-build-cmake-project-because-compatibility-with-cmake-3-5-has-been-remo
    # TODO: Remove this in the future or sight upper boundary if the issue is fixed
    IF(${CMAKE_VERSION} VERSION_LESS "4.0.0")
      SET(PETSC_CMAKEPOLICY "")
    ELSE()
      SET(PETSC_CMAKEPOLICY "3.5")
    ENDIF()

    # Fixes for GCC 15 and 16: Tested with PETSc version 3.22.5
    IF(${LIBS_BUILD_PETSC_VERSION} VERSION_LESS "3.24.0")
      IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        IF(${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "15.0.0")
          SET(HYPRE_COMPILER_FLAGS )
          SET(SCALAPACK_COMPILER_FLAGS )
        ELSE()
          # Fix "Error running make; make install on HYPRE" for 3.22.5 with GCC 15 and 16
          # GCC 15 changed the default C language standard from -std=gnu17 to -std=gnu23, and C23 added bool as a proper keyword — meaning older code that tries
          # to define bool via typedef now fails. This is exactly what HYPRE's older source does.
          # Set --download-hypre-configure-arguments=CFLAGS=-std=gnu17       for configure
          #     --download-hypre-cmake-arguments=-DCMAKE_C_FLAGS=-std=gnu17  for CMAKE
          SET(HYPRE_COMPILER_FLAGS -std=gnu17)
          # Fix "Error running make; make install on ScaLAPACK" for 3.22.5 with GCC 15 and 16: https://github.com/Reference-ScaLAPACK/scalapack/issues/129
          # "GCC-15 upped the default for -std= to gnu23 (more or less c23) from gnu18, which resulted in this error. The solution I found was to use -std=gnu90."
          # Set --download-scalapack-cmake-arguments=-DCMAKE_C_FLAGS=-std=gnu90
          SET(SCALAPACK_COMPILER_FLAGS -std=gnu90)
        ENDIF()
      ENDIF()
    ELSE()
      # Fix MPI Error: Name ‘mpi_comm_dup_fn’ at (1) is an ambiguous reference: PETSc internally uses use mpi (the old MPI Fortran module),
      # and when piclas uses mpi_f08, the compiler sees two conflicting sets of interfaces for the same MPI symbols (e.g. MPI_COMM_DUP_FN) causing the ambiguous reference error.
      SET(PETSC_MPI_COMPILER_FLAGS --with-mpi-ftn-module=mpi_f08)
    ENDIF()

    # Settings
    # --with-mpi-f90module-visibility=0       "With 0, mpi.mod will not be visible in use code (via petscsys.mod) - so mpi_f08 can now be used" (https://petsc.org/main/changes/315/)

    # Add PETSc project
    EXTERNALPROJECT_ADD(PETSc
        GIT_REPOSITORY "https://gitlab.com/petsc/petsc.git"
        GIT_TAG "v${LIBS_BUILD_PETSC_VERSION}"
        GIT_PROGRESS TRUE
        GIT_SHALLOW ON
        PREFIX ${LIBS_PETSC_DIR}
        INSTALL_DIR ${LIBS_EXTERNAL_LIB_DIR}/PETSc
        BUILD_IN_SOURCE TRUE
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND
          CMAKE_POLICY_VERSION_MINIMUM=${PETSC_CMAKEPOLICY}
          PETSC_DIR=${LIBS_PETSC_DIR}/src/PETSc
          ${LIBS_PETSC_DIR}/src/PETSc/configure
          --prefix=${PETSC_BUILD_DIR}
          ${PETSC_COMPILER_FLAGS}
          COPTFLAGS=${PETSC_OPTIMIZATION}
          CXXOPTFLAGS=${PETSC_OPTIMIZATION}
          FOPTFLAGS=${PETSC_OPTIMIZATION}
          ${PETSC_MPI_COMPILER_FLAGS}
          --with-shared-libraries=1
          --with-mpi-f90module-visibility=0
          --with-bison=0
          --download-hypre
          --download-hypre-configure-arguments=CFLAGS=${HYPRE_COMPILER_FLAGS} # -std=gnu17 for non-CMake HYPRE builds
          --download-mumps
          --download-scalapack
          --download-scalapack-cmake-arguments=-DCMAKE_C_FLAGS=${SCALAPACK_COMPILER_FLAGS} # -std=gnu90 for CMake ScaLAPACK builds
          --download-metis
          --download-parmetis     # requires metis
        # BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -j4
        BUILD_BYPRODUCTS ${PETSC_LIB_LIBRARIES}/libpetsc.so
    )

    # Actually add the PETSc lib paths to the linking paths
    INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES (${PETSC_LIB_LIBRARIES})
    LIST(APPEND linkedlibs ${PETSC_LIB_LIBRARIES}/libpetsc.so)
    LIST(APPEND SELFBUILTEXTERNALS PETSc)

    # If HDF5 is built in PIClas, it must occur before the PETSc compilation
    IF(LIBS_BUILD_HDF5)
      IF (NOT EXISTS "${LIBS_HDF5_DIR}/lib/libhdf5.a")
        ADD_DEPENDENCIES(PETSc HDF5)
      ENDIF()
    ENDIF()

    ADD_COMPILE_DEFINITIONS(USE_PETSC=1)
    MESSAGE(STATUS "Compiling with self-built [PETSc] (${BoldBlue}v${LIBS_BUILD_PETSC_VERSION}${ColourReset})")
  ELSE()
    IF(PETSC_FOUND)
      # Check if PETSc version needs FIX317
      IF (${PETSC_VERSION} VERSION_LESS 3.18)
        ADD_COMPILE_DEFINITIONS(USE_PETSC_FIX317=1)
      ELSE()
        ADD_COMPILE_DEFINITIONS(USE_PETSC_FIX317=0)
      ENDIF()

      INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})
      LIST(APPEND linkedlibs ${PETSC_LINK_LIBRARIES})

      ADD_COMPILE_DEFINITIONS(USE_PETSC=1)
      MESSAGE(STATUS "Compiling with system [PETSc] (${BoldBlue}v${PETSC_VERSION}${ColourReset}) [${PETSC_LINK_LIBRARIES}]")
    ELSE()
      MESSAGE(FATAL_ERROR "PETSc not found! Consider building PETSc with LIBS_BUILD_PETSC = ON.")
    ENDIF()
  ENDIF()

  # Do not allow PETSc + INT8 as this is not implemented yet
  IF(PICLAS_INTKIND8)
    MESSAGE(FATAL_ERROR "PETSc + INT8 is not implemented yet!")
  ENDIF()
ELSE()
  ADD_COMPILE_DEFINITIONS(USE_PETSC=0)
ENDIF(LIBS_USE_PETSC)
