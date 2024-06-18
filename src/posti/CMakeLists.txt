# =========================================================================
# POSTI
# =========================================================================

CMAKE_DEPENDENT_OPTION(POSTI_BUILD_SUPERB     "Build the superB tool for creating external magnetic field configurations"              ON "PICLAS_BUILD_POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_PICLAS2VTK "Build the piclas to VTK converter"                                                      ON "PICLAS_BUILD_POSTI" OFF)
CMAKE_DEPENDENT_OPTION(POSTI_BUILD_DMD        "Build the DMD tool"                                                                     OFF "PICLAS_BUILD_POSTI" OFF)

# Add CMakeLists

IF(POSTI_BUILD_SUPERB)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/posti/superB/CMakeLists.txt)
  ADD_DEFINITIONS(-DUSE_SUPER_B=1)
ELSE()
  ADD_DEFINITIONS(-DUSE_SUPER_B=0)
ENDIF()

IF(POSTI_BUILD_PICLAS2VTK)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/posti/piclas2vtk/CMakeLists.txt)
ENDIF()

IF(POSTI_BUILD_DMD)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/posti/dmd/CMakeLists.txt)
ENDIF()
