# =========================================================================
# piclas2vtk
# =========================================================================
SET(TOOL_NAME "piclas2vtk")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/posti/${TOOL_NAME}")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/${TOOL_NAME}.f90
                          ${TOOL_DIR}/${TOOL_NAME}_vars.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${PICLAS_COMPILE_FLAGS}
                                                    COMPILE_OPTIONS "-fPIC"
                                                    CMAKE_POSITION_INDEPENDENT_CODE ON)
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libpiclasshared)

# When using Ninja, CMake tries to parallelize the build resulting in missing module files. Add a dummy target to ensure synchronization
ADD_CUSTOM_TARGET(${TOOL_NAME}target)
ADD_DEPENDENCIES(${TOOL_NAME}target libpiclasshared)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 ${TOOL_NAME}target)

# Standalone Posti executable
add_exec(             ${TOOL_NAME} ${TOOL_DIR}/${TOOL_NAME}.f90)
TARGET_LINK_LIBRARIES(${TOOL_NAME} ${TOOL_NAME}lib libpiclasshared)
ADD_DEPENDENCIES(     ${TOOL_NAME} ${TOOL_NAME}lib libpiclasshared)
