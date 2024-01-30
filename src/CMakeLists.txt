# =========================================================================
# PICLas Options
# =========================================================================

# Build directory list
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
SUBDIRLIST(EQNSYSTEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/equations/)

OPTION(PICLAS_PARTICLES "Particles" ON)
# =========================================================================
# Time discretization
# =========================================================================
SET(PICLAS_TIMEDISCMETHOD "RK4" CACHE STRING "Time integration method")
SET_PROPERTY(CACHE PICLAS_TIMEDISCMETHOD PROPERTY STRINGS Euler-Explicit
                                                          Leapfrog
                                                          Boris-Leapfrog
                                                          Higuera-Cary
                                                          RK3
                                                          RK4
                                                          RK14
                                                          DSMC
                                                          FP-Flow
                                                          BGK-Flow
                                                          CN
                                                          ImplicitO3
                                                          ImplicitO4
                                                          Ascher3
                                                          ROS2
                                                          ROS3
                                                          ROS4
                                                          ROS46
                                                          ROS6
                                                          Radiation)

# =========================================================================
# Shared memory region splitting
# =========================================================================
# Get environment variable, e.g. PICLAS_SPLIT_TYPE=OMPI_COMM_TYPE_CORE
SET(PICLAS_SPLIT_TYPE "$ENV{PICLAS_SPLIT_TYPE}")
IF("${PICLAS_SPLIT_TYPE}" STREQUAL "")
  SET(PICLAS_SPLIT_TYPE "MPI_COMM_TYPE_SHARED")
ELSE()
  #IF("${PICLAS_SPLIT_TYPE}" STREQUAL "${PICLAS_SHARED_MEMORY}")
  MESSAGE(STATUS "Using user-defined environment variable [PICLAS_SPLIT_TYPE = ${PICLAS_SPLIT_TYPE}] for shared memory communicator splitting. Setting [PICLAS_SHARED_MEMORY = ${PICLAS_SPLIT_TYPE}]")
  #ENDIF()
ENDIF()
SET(PICLAS_SHARED_MEMORY "${PICLAS_SPLIT_TYPE}" CACHE STRING "Split type for creating new communicators based on colors and keys (requires MPI 3 or higher). Options with the prefix OMPI_ are specific to Open MPI.")
SET_PROPERTY(CACHE PICLAS_SHARED_MEMORY PROPERTY STRINGS MPI_COMM_TYPE_SHARED
                                                         OMPI_COMM_TYPE_CORE
                                                         PICLAS_COMM_TYPE_NODE)

ADD_DEFINITIONS(-DSharedMemoryMethod=${PICLAS_SHARED_MEMORY})
IF(PICLAS_SHARED_MEMORY STREQUAL "MPI_COMM_TYPE_SHARED")
  UNSET(PICLAS_SHARED_MEMORY_CORES CACHE)
  ADD_DEFINITIONS(-DCORE_SPLIT=0)
  MESSAGE(STATUS "Shared memory split type for subcommunicators set to node-level")
ELSEIF(PICLAS_SHARED_MEMORY STREQUAL "OMPI_COMM_TYPE_CORE")
  UNSET(PICLAS_SHARED_MEMORY_CORES CACHE)
  ADD_DEFINITIONS(-DCORE_SPLIT=1)
  MESSAGE(STATUS "Shared memory split type for subcommunicators set to core-level")
ELSEIF(PICLAS_SHARED_MEMORY STREQUAL "PICLAS_COMM_TYPE_NODE")
  SET(PICLAS_SHARED_MEMORY_CORES "2" CACHE STRING "Number of cores per node when setting PICLAS_SHARED_MEMORY=PICLAS_COMM_TYPE_NODE. All cores must be on the same physical node!")
  ADD_DEFINITIONS(-DCORE_SPLIT=${PICLAS_SHARED_MEMORY_CORES})
  MESSAGE(STATUS "Shared memory split type for subcommunicators set to sub-node-level with user-specific value [PICLAS_SHARED_MEMORY_CORES = ${PICLAS_SHARED_MEMORY_CORES}] cores per node")
ENDIF()

# =========================================================================
# MISC
# =========================================================================
CMAKE_DEPENDENT_OPTION(PICLAS_LOADBALANCE "Use load balance time measurement" ON "LIBS_USE_MPI" OFF)

OPTION(PICLAS_INTKIND8 "Enable particle numbers larger than 2147483647 (on most systems, huge(4-byte-integer)))" OFF)
IF(PICLAS_INTKIND8)
   ADD_DEFINITIONS(-DINTKIND8)
ENDIF()

OPTION(PICLAS_IMPLICIT_ANALYZE "IMPLICIT_ANALYZE" OFF)
IF(PICLAS_IMPLICIT_ANALYZE)
  ADD_DEFINITIONS(-DIMPLICIT_ANALYZE)
ENDIF(PICLAS_IMPLICIT_ANALYZE)

# =========================================================================
# Polynomial degree
# =========================================================================
SET(PICLAS_POLYNOMIAL_DEGREE "N" CACHE STRING "N: compile with variable polynomial degree. 1,2,3,...: compile with fixed polynomial degree for performance (untested).")
ADD_DEFINITIONS(-DPP_N=${PICLAS_POLYNOMIAL_DEGREE})

# =========================================================================
# Node type
# =========================================================================
SET(PICLAS_NODETYPE "GAUSS" CACHE STRING "Set of nodes used for representing the solution (GAUSS/GAUSS-LOBATTO).")
SET_PROPERTY(CACHE PICLAS_NODETYPE PROPERTY STRINGS "GAUSS" "GAUSS-LOBATTO")
IF(PICLAS_NODETYPE STREQUAL "GAUSS")
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSEIF(PICLAS_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_DEFINITIONS(-DPP_NodeType=2)
ELSE()
  MESSAGE(SEND_ERROR "Unknown node type.")
ENDIF()
ADD_DEFINITIONS(-DPP_nElems=NELEMZ)

# =========================================================================
# Equation system
# =========================================================================
SET(PICLAS_EQNSYSNAME "maxwell" CACHE STRING "Used equation system")
SET_PROPERTY(CACHE PICLAS_EQNSYSNAME PROPERTY STRINGS ${EQNSYSTEMS})

# =========================================================================
# Analytical Vandermonde
# =========================================================================
OPTION(PICLAS_VDM_ANALYTICAL "Alternative to matrix inversion: Compute inverse Vandermonde directly (required for inverse of Bezier and Legendre Vdm)" OFF)
# VDM inverse for Lagrange polynomials, replace lapack with analytical solution, for Bezier switch to LU inverse calculation
IF (FORCE_VDM_ANALYTICAL)
  SET(PICLAS_VDM_ANALYTICAL ON CACHE BOOL "Auto-Activated: Alternative to matrix inversion: Compute inverse Vandermonde directly (required for inverse of Bezier and Legendre Vdm)" FORCE)
ENDIF()
IF (PICLAS_VDM_ANALYTICAL)
  ADD_DEFINITIONS(-DVDM_ANALYTICAL)
ENDIF()

# =========================================================================
# Equation system
# =========================================================================
IF(PICLAS_EQNSYSNAME STREQUAL "maxwell")
  #SET(PICLAS_C_PREPROC ${F2P_C_PREPROC} -Dmaxwell)
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell)
  SET(PICLAS_HDG OFF)
ELSEIF(PICLAS_EQNSYSNAME STREQUAL "poisson")
  ADD_DEFINITIONS(-DPP_nVar=1 -Dpoisson)
  SET(PICLAS_HDG ON)
ELSEIF(PICLAS_EQNSYSNAME STREQUAL "maxwell_pois")
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell_pois -DPOIS)
  SET(PICLAS_HDG OFF)
  MESSAGE(SEND_ERROR "Equation system is currently not supported.")
ELSEIF(PICLAS_EQNSYSNAME STREQUAL "magnetostatic")
  ADD_DEFINITIONS(-DPP_nVar=3 -Dmagnetostatic)
  SET(PICLAS_HDG ON)
  MESSAGE(SEND_ERROR "Equation system is currently not supported.")
ELSEIF(PICLAS_EQNSYSNAME STREQUAL "magnetostatic_poisson")
  ADD_DEFINITIONS(-DPP_nVar=4 -Dmagnetostatic_poisson)
  SET(PICLAS_HDG ON)
  MESSAGE(SEND_ERROR "Equation system is currently not supported.")
ELSE()
   MESSAGE(SEND_ERROR "Unknown equation system.")
ENDIF()

# =========================================================================
# PARTICLES
# =========================================================================
IF(PICLAS_PARTICLES)
  ADD_DEFINITIONS(-DPARTICLES)
ELSE()
  IF(PICLAS_TIMEDISCMETHOD STREQUAL "DSMC")
    MESSAGE(SEND_ERROR "DSMC cannot be used without particles.")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "BGK-Flow")
    MESSAGE(SEND_ERROR "BGK-Flow cannot be used without particles.")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "FP-Flow")
    MESSAGE(SEND_ERROR "FP-Flow cannot be used without particles.")
  ENDIF()
ENDIF(PICLAS_PARTICLES)

# =========================================================================
# CODE_ANALYZE and EXTRAE instrumentation
# =========================================================================
IF(PICLAS_CODE_ANALYZE)
  ADD_DEFINITIONS(-DCODE_ANALYZE)
ENDIF(PICLAS_CODE_ANALYZE)

IF(PICLAS_EXTRAE)
  MESSAGE(STATUS "Activated extrae code tracing instrumentation using: $ENV{EXTRAE_HOME}")
  ADD_DEFINITIONS(-DEXTRAE)
ENDIF(PICLAS_EXTRAE)

# =========================================================================
# Timedisc method
# =========================================================================
SET(PICLAS_IMPLICIT OFF)
IF(PICLAS_TIMEDISCMETHOD STREQUAL "Euler-Explicit")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoisson.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=500)
  ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=0)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Higuera-Cary")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoissonByHigueraCary.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=507)
  ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=0)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Boris-Leapfrog")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoissonByBorisLeapfrog.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=508)
  ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=0)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Leapfrog")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoisson.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=509)
  ELSE()
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=0)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "RK3")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoissonByLSERK.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=501 -DLSERK)
  ELSE()
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByLSERK.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=1 -DLSERK)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "RK4")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoissonByLSERK.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=502 -DLSERK)
  ELSE()
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByLSERK.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=2 -DLSERK)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "RK14")
  IF(PICLAS_HDG)
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepPoissonByLSERK.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=506 -DLSERK)
  ELSE()
      FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByLSERK.f90)
      ADD_DEFINITIONS(-DPP_TimeDiscMethod=6 -DLSERK)
  ENDIF()
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "DSMC")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStep_DSMC.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=4)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "CN")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByImplicitRK.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=120 -DIMPA)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ImplicitO3")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByImplicitRK.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=121 -DIMPA)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ImplicitO4")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByImplicitRK.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=122 -DIMPA)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Ascher3")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByImplicitRK.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=123 -DIMPA)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS2")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByRosenbrock.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=130 -DROS)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS3")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByRosenbrock.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=131 -DROS)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS4")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByRosenbrock.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=132 -DROS)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS46")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByRosenbrock.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=133 -DROS)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS6")
  FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStepByRosenbrock.f90)
  ADD_DEFINITIONS(-DPP_TimeDiscMethod=134 -DROS)
  SET(PICLAS_IMPLICIT ON)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "FP-Flow")
   FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStep_FPFlow.f90)
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=300)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "BGK-Flow")
   FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStep_BGK.f90)
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=400)
ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Radiation")
   FILE(GLOB timediscF90 ${CMAKE_CURRENT_SOURCE_DIR}/src/timedisc/timedisc_TimeStep_Radiation.f90)
   ADD_DEFINITIONS(-DPP_TimeDiscMethod=600)
ELSE()
   MESSAGE(SEND_ERROR "Unknown time discretization method.")
ENDIF()

# =========================================================================
# HDG
# =========================================================================

IF(PICLAS_HDG)
  ADD_DEFINITIONS(-DUSE_HDG=1)
  IF(PICLAS_TIMEDISCMETHOD STREQUAL "Euler-Explicit")
      MESSAGE(STATUS "use HDG with Euler-Explicit")
    ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Higuera-Cary")
      MESSAGE(STATUS "use HDG with Higuera-Cary")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Boris-Leapfrog")
      MESSAGE(STATUS "use HDG with Boris-Leapfrog")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Leapfrog")
      MESSAGE(STATUS "use HDG with Leapfrog")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "RK3")
      MESSAGE(STATUS "use HDG with RK3")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "RK4")
      MESSAGE(STATUS "use HDG with RK4")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "RK14")
      MESSAGE(STATUS "use HDG with RK14")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "CN")
      MESSAGE(STATUS "use HDG with CN")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ImplicitO3")
     MESSAGE(STATUS "use HDG with ImplicitO3")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ImplicitO4")
     MESSAGE(STATUS "use HDG with ImplicitO4")
   ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "Ascher3")
     MESSAGE(STATUS "use HDG with Ascher3")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS2")
     MESSAGE(STATUS "use HDG with ROS2")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS4")
     MESSAGE(STATUS "use HDG with ROS4")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS46")
     MESSAGE(STATUS "use HDG with ROS46")
  ELSEIF(PICLAS_TIMEDISCMETHOD STREQUAL "ROS6")
     MESSAGE(STATUS "use HDG with ROS6")
  ELSE()
      MESSAGE(SEND_ERROR "Unknown time discretization method for HDG.")
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DUSE_HDG=0)
ENDIF(PICLAS_HDG)

# ========================================================================
# LOADBALANCE
# =========================================================================
IF(LIBS_USE_MPI)
  IF(PICLAS_LOADBALANCE)
    MESSAGE(STATUS "Building PICLas with LoadBalance")
    ADD_DEFINITIONS(-DUSE_LOADBALANCE=1)
  ELSE()
    ADD_DEFINITIONS(-DUSE_LOADBALANCE=0)
  ENDIF()
ELSE()
  ADD_DEFINITIONS(-DUSE_LOADBALANCE=0)
ENDIF()

# =========================================================================
# Coffee Lib
# =========================================================================
OPTION(CMAKE_ANS_COFFEE "Enable coffee lib" OFF)
MARK_AS_ADVANCED(CMAKE_ANS_COFFEE)
IF(CMAKE_ANS_COFFEE)
  ADD_DEFINITIONS(-DUSE_COFFEE=1)
ELSE()
  ADD_DEFINITIONS(-DUSE_COFFEE=0)
ENDIF()

# =========================================================================
# Unit test at end of compilation
# =========================================================================
# Check is unit tests can be activated
IF(BUILD_TYPE_LC MATCHES "sanitize" AND LIBS_USE_MPI)
  SET(ACTIVATE_UNIT_TEST OFF)
  # Check if unit tests were previously activated
  IF(PICLAS_UNITTESTS AND LIBS_MPI_NAME MATCHES "OpenMPI")
    MESSAGE(WARNING "Unit tests automatically disabled (PICLAS_UNITTESTS=OFF) with\n CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\n LIBS_USE_MPI = ${LIBS_USE_MPI}\nas the tests will always fail due to memory errors returned from the MPI library.")
    SET(PICLAS_UNITTESTS OFF)
  ENDIF()
ELSE()
  SET(ACTIVATE_UNIT_TEST ON)
ENDIF()

# Enable unit testing only when ACTIVATE_UNIT_TEST=ON
CMAKE_DEPENDENT_OPTION(PICLAS_UNITTESTS "Enable unit tests after build" ON "ACTIVATE_UNIT_TEST" OFF)

IF(PICLAS_UNITTESTS)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/unitTests/CMakeLists.txt)
ENDIF()

# =========================================================================
# PICLas sources & targets
# =========================================================================
FUNCTION(add_exec target)
  ADD_EXECUTABLE(       ${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} libpiclasshared)
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${PICLAS_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

FUNCTION(add_lib_shared target)
  ADD_LIBRARY(          ${target} SHARED ${ARGN})
  TARGET_LINK_LIBRARIES(${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.o libpiclasstatic userblocklib stacksizelib memusagelib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES PREFIX ""
                                             COMPILE_FLAGS ${PICLAS_COMPILE_FLAGS}
                                             LINKER_LANGUAGE Fortran)
  INSTALL(TARGETS ${target} LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                                                GROUP_READ GROUP_EXECUTE
                                                WORLD_READ WORLD_EXECUTE
                            DESTINATION lib)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

# include piclas.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

# collect all F90 files needed for libpiclasstatic
FILE(GLOB_RECURSE piclasF90 ./src/init/*.f90
                            ./src/readIMD/*.f90
                            ./src/analyze/*.f90
                            ./src/filter/*.f90
                            ./src/globals/*.f90
                            ./src/mortar/*.f90
                            ./src/interpolation/basis.f90
                            ./src/interpolation/interpolation.f90
                            ./src/interpolation/interpolation_vars.f90
                            ./src/interpolation/changeBasis.f90
                            ./src/io_hdf5/*.f90
                            ./src/loadbalance/*.f90
                            ./src/mesh/*.f90
                            ./src/interfaces/*.f90
                            ./src/dielectric/*.f90
                            ./src/mpi/*.f90
                            ./src/output/*.f90
                            ./src/readintools/*.f90
                            ./src/recordpoints/*.f90
                            ./src/restart/*.f90
                            ./src/utils/*.f90
                            ./unitTests/unittest.f90
                            ./unitTests/unittest_vars.f90)

LIST(APPEND timediscF90 ${timediscF90}
                        ./src/timedisc/timedisc.f90
                        ./src/timedisc/timedisc_init.f90
                        ./src/timedisc/timedisc_vars.f90)

LIST(APPEND piclasF90 ${timediscF90})

FILE(GLOB_RECURSE dgF90 ./src/dg/*.f90
                        ./src/interpolation/prolongtoface.f90
                        ./src/pml/*.f90
                        ./src/equations/${PICLAS_EQNSYSNAME}/*.f90)

FILE(GLOB_RECURSE hdgF90 ./src/hdg/*.f90
                         ./src/dg/dg.f90
                         ./src/dg/dg_vars.f90
                         ./src/interpolation/prolongtoface.f90
                         ./src/equations/${PICLAS_EQNSYSNAME}/equation.f90
                         ./src/equations/${PICLAS_EQNSYSNAME}/equation_tools.f90
                         ./src/equations/${PICLAS_EQNSYSNAME}/equation_vars.f90
                         ./src/equations/${PICLAS_EQNSYSNAME}/timeavg.f90
                         ./src/equations/${PICLAS_EQNSYSNAME}/timeavg_vars.f90
                         ./src/equations/${PICLAS_EQNSYSNAME}/getboundaryflux.f90)


FILE(GLOB_RECURSE implicitF90 ./src/linearsolver/linearsolver.f90
                              ./src/linearsolver/linearsolver_vars.f90
                              ./src/linearsolver/particlesolver.f90
                              ./src/linearsolver/predictor.f90
                              ./src/linearsolver/linearoperator/linearoperator.f90
                              ./src/linearsolver/newton.f90)

IF(PICLAS_EQNSYSNAME STREQUAL "maxwell_pois")
  FILE(GLOB_RECURSE precondF90 ./src/precond/precond.f90
                               ./src/precond/precond_vars.f90
                               ./src/precond/applypreconditioner.f90
                               ./src/linearsolver/linearoperator/sparse_ilu.f90
                               ./src/precond/csr/*.f90
                               ./src/precond/ilu/*.f90
                               ./src/precond/jacdg/*.f90
                               ./src/precond/jac_fd/*.f90
                               ./src/precond/jac_ex/${PICLAS_EQNSYSNAME}/*.f90)
ELSE()
  FILE(GLOB_RECURSE precondF90 ./src/precond/precond.f90
                               ./src/precond/precond_vars.f90
                               ./src/precond/applypreconditioner.f90
                               ./src/linearsolver/linearoperator/sparse_ilu.f90
                               ./src/precond/csr/*.f90
                               ./src/precond/ilu/*.f90
                               ./src/precond/jacdg/*.f90
                               ./src/precond/jac_fd/*.f90
                               ./src/precond/jac_ex/maxwell/*.f90)
ENDIF()

#FILE(GLOB mainF90 ./src/piclas.f90
                  #./src/piclaslib.f90)
FILE(GLOB_RECURSE particlesF90 ./src/particles/*.f90
                               ./src/posti/superB/*.f90
                               ./src/interpolation/eval_xyz.f90
                               ./src/radiation/ray_tracing/*.f90
                               ./src/radiation/radiative_transfer/radtrans_vars.f90
                               ./src/radiation/radiative_transfer/tracking/*.f90
                               ./src/radiation/radiation_solver/radiation_vars.f90)

FILE(GLOB noparticlesF90 ./src/particles/analyze/particle_analyze_vars.f90
                         ./src/particles/particle_mesh/particle_surfaces_vars.f90
                         ./src/particles/particle_mesh/particle_mesh_vars.f90
                         ./src/particles/particle_vars.f90
                         ./src/particles/emission/particle_emission_vars.f90
                         ./src/particles/emission/particle_surface_flux_vars.f90
                         ./src/posti/superB/superB_vars.f90
                         ./src/particles/boundary/particle_boundary_vars.f90)

FILE(GLOB_RECURSE extraeF90 ./src/extrae/*.f90)

FILE(GLOB_RECURSE radiationF90 ./src/radiation/*.f90)

IF(PICLAS_EXTRAE)
  LIST(APPEND piclasF90 ${extraeF90})
ENDIF(PICLAS_EXTRAE)

IF(PICLAS_HDG)
  LIST(APPEND piclasF90 ${hdgF90})
  IF(PICLAS_IMPLICIT)
    LIST(APPEND piclasF90 ${implicitF90})
  ENDIF(PICLAS_IMPLICIT)
ELSE()
  LIST(FILTER particlesF90 EXCLUDE REGEX ".*particle_br_electron_fluid.f90$")
  LIST(FILTER particlesF90 EXCLUDE REGEX ".*superb.f90$")
  LIST(APPEND piclasF90 ${dgF90})
  IF(PICLAS_IMPLICIT)
    LIST(APPEND piclasF90 ${implicitF90})
    LIST(APPEND piclasF90 ${precondF90})
  ENDIF(PICLAS_IMPLICIT)
ENDIF(PICLAS_HDG)

#LIST(APPEND piclasF90 ${mainF90})
IF (PICLAS_PARTICLES)
  LIST(APPEND piclasF90 ${particlesF90})
ELSE()
  LIST(APPEND piclasF90 ${noparticlesF90})
ENDIF(PICLAS_PARTICLES)
IF (PICLAS_TIMEDISCMETHOD STREQUAL "Radiation")
  LIST(APPEND piclasF90 ${radiationF90})
ENDIF()

# define libs
ADD_LIBRARY(          libpiclasstaticF90 OBJECT ${piclasF90} )
SET_TARGET_PROPERTIES(libpiclasstaticF90 PROPERTIES COMPILE_FLAGS ${PICLAS_COMPILE_FLAGS})
ADD_LIBRARY(          libpiclasstatic STATIC $<TARGET_OBJECTS:libpiclasstaticF90>)
SET_TARGET_PROPERTIES(libpiclasstatic PROPERTIES PREFIX ""
                                      OUTPUT_NAME "libpiclasf90")
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(libpiclasstaticF90 ${SELFBUILTEXTERNALS})
ENDIF()

ADD_LIBRARY(userblocklib STATIC ./src/output/read_userblock.c)
SET_TARGET_PROPERTIES(userblocklib PROPERTIES LANGUAGE C LINKER_LANGUAGE C
                                              COMPILE_FLAGS "-fPIC")

ADD_LIBRARY(stacksizelib STATIC ./src/globals/stacksize.c)
SET_TARGET_PROPERTIES(stacksizelib PROPERTIES LANGUAGE C LINKER_LANGUAGE C
                                              COMPILE_FLAGS "-fPIC")

ADD_LIBRARY(memusagelib STATIC ./src/globals/processmemusage.cxx)
SET_TARGET_PROPERTIES(memusagelib PROPERTIES LANGUAGE CXX LINKER_LANGUAGE CXX
                                             COMPILE_FLAGS "-fPIC")

# we link the statically built libs
ADD_LIB_SHARED(libpiclasshared ./src/piclaslib.f90)
SET_TARGET_PROPERTIES(libpiclasshared PROPERTIES OUTPUT_NAME "libpiclas")
ADD_DEPENDENCIES(libpiclasshared libpiclasstatic userblocklib stacksizelib memusagelib ${INTERNALLIBS})
ADD_EXEC(piclas ./src/piclas.f90)
ADD_DEPENDENCIES(piclas libpiclasstatic userblocklib stacksizelib memusagelib ${INTERNALLIBS})
# Make sure that the script is executed even though cmake is not explicitly executed after the commit is changed in git
ADD_DEPENDENCIES(piclas UpdateGitCommitHash)

# special compiler flags to avoid warnings for empty functions of testcase.f90 and overwrite of intrinsic abort in globals.f90
IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals.f90 PROPERTIES COMPILE_FLAGS "-Wno-intrinsic-shadow")
ENDIF()

# =========================================================================
# EXTRAE instrumentation
# =========================================================================
IF(PICLAS_EXTRAE)
  # MPI+OpenMP
  #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lompitracef -L$ENV{EXTRAE_HOME}/lib")
  # only MPI
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmpitracef -L$ENV{EXTRAE_HOME}/lib")
  # change binary name by adding the suffix '_extrae'
  SET_TARGET_PROPERTIES(piclas PROPERTIES OUTPUT_NAME "piclas_extrae")
ENDIF(PICLAS_EXTRAE)

# =========================================================================
# Shared memory debug flags
# =========================================================================
SET(PICLAS_DEBUG_MEMORY OFF CACHE BOOL "Output shared memory window sizes. 'USE MOD_MPI_Shared' cannot be used together with the 'ONLY:' statement in this case!")
IF(PICLAS_DEBUG_MEMORY)
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDEBUG_MEMORY")
ENDIF()

# =========================================================================
# MPI_WAIT measurement and output
# =========================================================================
CMAKE_DEPENDENT_OPTION(PICLAS_MEASURE_MPI_WAIT "Output measured and accumulated wall time during MPI_WAIT for each process." OFF "LIBS_USE_MPI" OFF)
IF(LIBS_USE_MPI)
  IF(PICLAS_MEASURE_MPI_WAIT)
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMEASURE_MPI_WAIT")
    MESSAGE(STATUS "Activated measurement of MPI_WAIT [-DMEASURE_MPI_WAIT]")
  ENDIF()
ENDIF()

# =========================================================================
# Ctags
# =========================================================================
OPTION(PICLAS_CTAGS "Use Ctags" ON)
IF(PICLAS_CTAGS)
  MARK_AS_ADVANCED(FORCE CTAGS_PATH)
  FIND_PROGRAM(CTAGS_PATH ctags-exuberant)
  IF(CTAGS_PATH)
    SET(             EQNS_WO_ACTUAL ${EQNSYSTEMS})
    LIST(REMOVE_ITEM EQNS_WO_ACTUAL ${PICLAS_EQNSYSNAME})
    SET(EXLUDES "")
    FOREACH(item ${EQNS_WO_ACTUAL})
       SET(EXLUDES "--exclude=src/equations/${item}/*" ${EXLUDES})
    ENDFOREACH()

    ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} src && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} --append src/posti)
    ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
    MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
    ADD_DEPENDENCIES(piclas tags)
  ENDIF(CTAGS_PATH)
ENDIF(PICLAS_CTAGS)
