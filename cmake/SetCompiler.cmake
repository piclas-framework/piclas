# TODO: Is this still required? Where does it need to go?
SET(FORCE_VDM_ANALYTICAL OFF)

# =========================================================================
# After settings specific compilers, enable named languages for cmake
# =========================================================================
ENABLE_LANGUAGE(Fortran C CXX)
INCLUDE(GNUInstallDirs)
MARK_AS_ADVANCED(FORCE C_PATH CXX_PATH Fortran_PATH)

# =========================================================================
# Set machine-specific definitions and settings
# =========================================================================
# HLRS Vulcan
IF (CMAKE_FQDN_HOST MATCHES "^cl[0-9]fr")
  MESSAGE(STATUS "Compiling on Vulcan for AMD EPYC Genoa")
  # Overwrite compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    SET(PICLAS_INSTRUCTION "-march=znver4 -mtune=znver4" CACHE STRING "Compiler optimization options")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(PICLAS_INSTRUCTION "-xCORE-AVX2" CACHE STRING "Compiler optimization options")
  ENDIF()
  # Set LUSTRE definition to account for filesystem and MPI implementation
  ADD_COMPILE_DEFINITIONS(LUSTRE)

# SuperMUC
ELSEIF (CMAKE_FQDN_HOST MATCHES "sng\.lrz\.de$")
  MESSAGE(STATUS "Compiling on SuperMUC")
  # Overwrite compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    SET(PICLAS_INSTRUCTION "-march=skylake-avx512 -mtune=skylake-avx512" CACHE STRING "Compiler optimization options")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(PICLAS_INSTRUCTION "-xSKYLAKE-AVX512" CACHE STRING "Compiler optimization options")
    # Explicitely enable usage of AVX512 registers
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-zmm-usage=high")
  ENDIF()
  # Set LUSTRE definition to account for filesystem and MPI implementation
  ADD_COMPILE_DEFINITIONS(LUSTRE)

# LUMI
ELSEIF(CMAKE_FQDN_HOST MATCHES "\.can$")
  MESSAGE(STATUS "Compiling on LUMI")
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(PICLAS_INSTRUCTION "-march=znver3 -mtune=znver3" CACHE STRING "Compiler optimization options")
  ELSE()
    MESSAGE(FATAL_ERROR "LUMI currently only supported using the GNU Compiler Collection (GCC). Please load/swap the following modules: LUMI PrgEnv-gnu cray-hdf5-parallel")
  ENDIF()

# IAG Prandtl
ELSEIF(CMAKE_FQDN_HOST MATCHES "^(prandtl|grafik.*)\.iag\.uni\-stuttgart\.de")
  MESSAGE(STATUS "Compiling on ${CMAKE_HOSTNAME}")
  SET(PICLAS_INSTRUCTION "-march=native -mtune=native" CACHE STRING "Compiler optimization options")
  # Set LUSTRE definition to account for filesystem
  ADD_COMPILE_DEFINITIONS(LUSTRE)

ELSE()
  MESSAGE(STATUS "Compiling on a generic machine [${CMAKE_HOSTNAME}]")
  # Set compiler target architecture
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Flang" OR CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    SET(PICLAS_INSTRUCTION "-march=native -mtune=native" CACHE STRING "Compiler optimization options")
  ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    SET(PICLAS_INSTRUCTION "-xHost" CACHE STRING "Compiler optimization options")
  ENDIF()
ENDIF()

MESSAGE(STATUS "Compiling Nitro/Release/Profile with [${CMAKE_Fortran_COMPILER_ID}] (v${CMAKE_Fortran_COMPILER_VERSION}) fortran compiler using PICLAS_INSTRUCTION [${PICLAS_INSTRUCTION}] instructions.")

# =========================================================================
# CHECK SUPPORT FOR VARIOUS FORTRAN (2003,2008) FEATURES
# =========================================================================
INCLUDE(CheckFortranSourceCompiles)
CHECK_FORTRAN_SOURCE_COMPILES(
"MODULE F03MOD
IMPLICIT NONE
PRIVATE
CONTAINS
! Check for MOVE_ALLOC feature
SUBROUTINE F03_MOVE_ALLOC()
REAL,ALLOCATABLE::a(:,:),b(:,:)
ALLOCATE(a(3,3))
CALL MOVE_ALLOC(a,b)
END SUBROUTINE F03_MOVE_ALLOC
END MODULE F03MOD

PROGRAM F03PROG
END"
Fortran2003Check
SRC_EXT F90)

IF(NOT Fortran2003Check)
  MESSAGE(FATAL_ERROR "Failed to compile basic Fortran2003 programm! Please ensure your compiler is up-to-date!")
ENDIF()

# =========================================================================
# CHECK IF GCC CONTAINS THE PARSE_ASSOCIATE BUG (GC C13.1-13.2)
# =========================================================================
INCLUDE(CheckFortranSourceCompiles)
CHECK_FORTRAN_SOURCE_COMPILES(
"MODULE GCC13MOD
IMPLICIT NONE
PRIVATE
CONTAINS
! Check for MOVE_ALLOC feature
SUBROUTINE GCC_PARSE_ASSOCIATE()
REAL::x1(3),x2(3)
ASSOCIATE(v1 => x1, v2 => x2)
v1 = 0
v2 = 1
v1 = 0 + v1
END ASSOCIATE
END SUBROUTINE GCC_PARSE_ASSOCIATE
END MODULE GCC13MOD

PROGRAM GCC13PROG
END"
GCC13Check
SRC_EXT F90)

IF(NOT GCC13Check)
  MESSAGE(WARNING "The requested compiler ${CMAKE_Fortran_COMPILER_ID} (v${CMAKE_Fortran_COMPILER_VERSION}) contains a bug in parse_associate. Ensure ASSOCIATE is only use with explicit array bounds or use a different compiler version. Please see the upstream issue, https://gcc.gnu.org/bugzilla/show_bug.cgi?id=109948")
ENDIF()

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# CMake can always request position independent code
# SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

# GNU Compiler Collection (GCC)
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set Flags (disable lto type warnings due to false positives with MATMUL, which is a known bug)
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -finit-real=snan -finit-integer=snan -Wno-lto-type-mismatch -lstdc++ -DGNU")
    IF(PICLAS_EXTRAE)
      SET (CMAKE_Fortran_FLAGS       "${CMAKE_Fortran_FLAGS} -finstrument-functions")
    ENDIF(PICLAS_EXTRAE)
    # Allow type mistach for argument types in MPI(CH) calls
    IF(MPI_C_LIBRARY_VERSION_STRING MATCHES ".*MPICH.*")
      SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    ENDIF()
  ENDIF()
  # initialize all variables as signalling NaNs to force the user to correctly initialize these data types
  SET (CMAKE_Fortran_FLAGS_NITRO       "${CMAKE_Fortran_FLAGS}  -Ofast ${PICLAS_INSTRUCTION}")
  SET (CMAKE_Fortran_FLAGS_RELEASE     "${CMAKE_Fortran_FLAGS}     -O3 ${PICLAS_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_PROFILE     "${CMAKE_Fortran_FLAGS} -pg -O3 ${PICLAS_INSTRUCTION} -finline-functions -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_NOGPROFILE  "${CMAKE_Fortran_FLAGS} -g1 -O3 ${PICLAS_INSTRUCTION} -fno-inline -fno-optimize-sibling-calls -fstack-arrays")
  SET (CMAKE_Fortran_FLAGS_DEBUG       "${CMAKE_Fortran_FLAGS} -g  -Og -ggdb3 -ffpe-trap=invalid,zero,overflow -fbounds-check -fbacktrace -Wall")
  SET (CMAKE_Fortran_FLAGS_SANITIZE    "${CMAKE_Fortran_FLAGS} -g  -Og -ggdb3 -ffpe-trap=invalid,zero,overflow,denorm -fbounds-check -fbacktrace -Wall -fsanitize=address,undefined,leak -fno-omit-frame-pointer -Wc-binding-type -Wuninitialized -pedantic")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (PICLAS_COMPILE_FLAGS "-xf95-cpp-input")
  ELSE()
    # Trailing white space required in case variable is unset!
    SET (PICLAS_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} ")
  ENDIF()

# AMD Optimized LLVM/CLANG
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -std=f2008 -lstdc++ -DFLANG")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE   "${CMAKE_Fortran_FLAGS}     -O3 ${PICLAS_INSTRUCTION} -finline-functions ")
  SET (CMAKE_Fortran_FLAGS_PROFILE   "${CMAKE_Fortran_FLAGS} -pg -O3 ${PICLAS_INSTRUCTION} -finline-functions ")
  SET (CMAKE_Fortran_FLAGS_DEBUG     "${CMAKE_Fortran_FLAGS} -g  -O0 -ggdb3 -ffpe-trap=invalid,zero,overflow -fbounds-check -finit-real=snan -fbacktrace -Wall")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (PICLAS_COMPILE_FLAGS "-xf95-cpp-input")
  ELSE()
    # Trailing white space required in case variable is unset!
    SET (PICLAS_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} ")
  ENDIF()

# Intel Compiler
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -lstdc++ -DINTEL")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE   "${CMAKE_Fortran_FLAGS}    -O3 ${PICLAS_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_PROFILE   "${CMAKE_Fortran_FLAGS} -p -O3 ${PICLAS_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG     "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
  # Compile flags depend on the generator
  IF(NOT "${CMAKE_GENERATOR}" MATCHES "Ninja")
    # add flags only for compiling not linking!
    SET (PICLAS_COMPILE_FLAGS "-fpp -allow nofpp_comments -assume bscc")
  ELSE()
    SET (PICLAS_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} -allow nofpp_comments -assume bscc")
  ENDIF()

# Cray Compiler
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  IF (NOT DEFINED C_FLAGS_INITIALIZED )
    SET (C_FLAGS_INITIALIZED "yes" CACHE INTERNAL "Flag if compiler flags are already initialized" )
    SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -ffree -s real64 -s integer64 -em -lstdc++ -hfp0 -DCRAY")
  ENDIF()
  SET (CMAKE_Fortran_FLAGS_RELEASE   "${CMAKE_Fortran_FLAGS} -s integer32  -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_PROFILE   "${CMAKE_Fortran_FLAGS} -s integer32  -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG     "${CMAKE_Fortran_FLAGS} -s integer32 -g -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (PICLAS_COMPILE_FLAGS "${NINJA_COLOR_DIAGNOSTICS} -F")
ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# Profile-Guided Optimization (PGO)
# =========================================================================
CMAKE_DEPENDENT_OPTION(USE_PGO "Enable profile-guided optimization (Only GNU Compiler supported)" OFF
                               "PICLAS_PERFORMANCE" OFF)
IF (USE_PGO)
  IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fprofile-use")
    SET(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -fprofile-generate")
  ELSE()
    MESSAGE(SEND_ERROR "Profile-guided optimization (PGO) currently only supported for GNU compiler. Either set USE_PGO=OFF or use the GNU compiler." )
  ENDIF()
ENDIF()

# =========================================================================
# SAVE CURRENT COMPILER FLAGS TO THE CACHE
# =========================================================================
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_RELEASE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_PROFILE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_NOGPROFILE)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_DEBUG)
MARK_AS_ADVANCED(FORCE CMAKE_Fortran_FLAGS_SANITIZE)
SET(CMAKE_Fortran_FLAGS            "${CMAKE_Fortran_FLAGS}"            CACHE STRING "Default compiler flags"  FORCE)
SET(CMAKE_Fortran_FLAGS_RELEASE    "${CMAKE_Fortran_FLAGS_RELEASE}"    CACHE STRING "Release compiler flags"  FORCE)
SET(CMAKE_Fortran_FLAGS_PROFILE    "${CMAKE_Fortran_FLAGS_PROFILE}"    CACHE STRING "Profile compiler flags"  FORCE)
SET(CMAKE_Fortran_FLAGS_NOGPROFILE "${CMAKE_Fortran_FLAGS_NOGPROFILE}" CACHE STRING "NoGProfile compiler flags" FORCE)
SET(CMAKE_Fortran_FLAGS_DEBUG      "${CMAKE_Fortran_FLAGS_DEBUG}"      CACHE STRING "Debug compiler flags"    FORCE)
SET(CMAKE_Fortran_FLAGS_SANITIZE   "${CMAKE_Fortran_FLAGS_SANITIZE}"   CACHE STRING "Sanitize compiler flags" FORCE)