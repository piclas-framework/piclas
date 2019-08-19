# =========================================================================
# VISU
# =========================================================================
SET(TOOL_NAME "ParaviewReader")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/visu")

MESSAGE(STATUS "VISU: building ${TOOL_NAME}")
INCLUDE_DIRECTORIES( ${TOOL_DIR})

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/CMakeLists.txt
                          ${TOOL_DIR}/converttovisu.f90
                          ${TOOL_DIR}/mappings.f90
                          ${TOOL_DIR}/paraviewReader/
                          ${TOOL_DIR}/pluginTypes_visu.h
                          ${TOOL_DIR}/plugin_visu.h
                          ${TOOL_DIR}/readstate.f90
                          ${TOOL_DIR}/visu_Cwrapper.f90
                          ${TOOL_DIR}/visu.f90
                          ${TOOL_DIR}/visumesh.f90
                          ${TOOL_DIR}/visu_vars.f90)


ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${PICLAS_COMPILE_FLAGS})
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libpiclasshared)

# Posti visu plugin for Paraview
FIND_PACKAGE(ParaView)
IF (ParaView_FOUND)
   INCLUDE(${PARAVIEW_USE_FILE})
   MESSAGE(STATUS ${PARAVIEW_USE_FILE})

   SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   INCLUDE_DIRECTORIES(${TOOL_DIR}/paraviewReader/)
   ADD_PARAVIEW_PLUGIN(visuReader "1.0"
     SERVER_MANAGER_XML     ${TOOL_DIR}/paraviewReader/visuReader.xml
     SERVER_MANAGER_SOURCES ${TOOL_DIR}/paraviewReader/visuReader.cpp
     REQUIRED_ON_SERVER
     DEPENDS               "${TOOL_DIR}/paraviewReader/visuReader.h"
     )

   ADD_DEPENDENCIES(     visuReader ${TOOL_NAME}lib libpiclasshared)
   TARGET_LINK_LIBRARIES(visuReader PUBLIC ${TOOL_NAME}lib libpiclasshared)
ENDIF()
