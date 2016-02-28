#
# Find the native TECIO includes and library
#
# TECIO_LIBRARIES   - List of fully qualified libraries to link against when using TECIO.
# TECIO_FOUND       - Do not attempt to use TECIO if "no" or undefined.

find_library(TECIO_LIBRARY tecio
  HINTS $ENV{TECIO_DIR}/lib
        $ENV{TECIO_ROOT}/lib
        /usr/local/lib
        /usr/lib
)

set(TECIO_FOUND "NO")
if(TECIO_LIBRARY)
  set( TECIO_LIBRARIES ${TECIO_LIBRARY} -lstdc++)
  set( TECIO_FOUND "YES" )
endif()

if(TECIO_FIND_REQUIRED AND NOT TECIO_FOUND)
  message(SEND_ERROR "Unable to find the requested TECIO libraries.")
endif()

# handle the QUIETLY and REQUIRED arguments and set TECIO_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TECIO DEFAULT_MSG TECIO_LIBRARY)


mark_as_advanced(
  TECIO_LIBRARY
)
