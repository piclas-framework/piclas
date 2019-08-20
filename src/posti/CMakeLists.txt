# =========================================================================
# POSTI
# =========================================================================

CMAKE_DEPENDENT_OPTION(POSTI_BUILD_VISU "Build the visu tool" ON "PICLAS_BUILDPOSTI" OFF)


# Add CMakeLists
IF(POSTI_BUILD_VISU)
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/src/posti/visu/CMakeLists.txt)
ENDIF()

