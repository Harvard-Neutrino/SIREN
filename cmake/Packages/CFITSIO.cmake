################################################################################
# Module to find cfitsio                                                       #
#                                                                              #
# This module defines:                                                         #
#                                                                              #
#   CFITSIO_FOUND                                                              #
#   CFITSIO_VERSION                                                            #
#   CFITSIO_LIBRARIES                                                          #
#   CFITSIO_INCLUDE_DIR                                                        #
#   CFITSIO_LIB_DIR                                                            #
#   CFITSIO_CPPFLAGS                                                           #
#   CFITSIO_LDFLAGS                                                            #
################################################################################

SET (CFITSIO_FIND_QUIETLY TRUE)
SET (CFITSIO_FIND_REQUIRED TRUE)

IF (NOT CFITSIO_FOUND)

  # Search user environment for headers, then default paths; extract version
  FIND_PATH (CFITSIO_INCLUDE_DIR fitsio.h
    PATHS $ENV{CFITSIOROOT}/include
    NO_DEFAULT_PATH)
  FIND_PATH (CFITSIO_INCLUDE_DIR fitsio.h)
  if(CFITSIO_INCLUDE_DIR)
    GET_FILENAME_COMPONENT (CFITSIOROOT ${CFITSIO_INCLUDE_DIR} PATH)
  else()
    FIND_PATH (CFITSIO_INCLUDE_DIR cfitsio/fitsio.h
      PATHS $ENV{CFITSIOROOT}/include
      )
    SET(CFITSIO_INCLUDE_DIR "${CFITSIO_INCLUDE_DIR}/cfitsio" CACHE PATH "Path to cfitsio headers" FORCE)
  endif()

  SET (CFITSIO_VERSION 0)
  IF (CFITSIO_INCLUDE_DIR)
    FILE (READ "${CFITSIO_INCLUDE_DIR}/fitsio.h" _cfitsio_VERSION)
    STRING (REGEX REPLACE ".*define CFITSIO_VERSION ([0-9]+\\.[0-9]+).*" "\\1"
            CFITSIO_VERSION "${_cfitsio_VERSION}")
  ENDIF (CFITSIO_INCLUDE_DIR)

  # Search user environment for libraries, then default paths
  FIND_LIBRARY (CFITSIO_LIBRARIES NAMES cfitsio
    PATHS $ENV{CFITSIOROOT}/lib
    NO_DEFAULT_PATH)
  FIND_LIBRARY (CFITSIO_LIBRARIES NAMES cfitsio)
  GET_FILENAME_COMPONENT (CFITSIO_LIB_DIR ${CFITSIO_LIBRARIES} PATH)

  # Set CFITSIO_FOUND and error out if cfitsio is not found
  INCLUDE (FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS (CFITSIO
    DEFAULT_MSG CFITSIO_LIBRARIES CFITSIO_INCLUDE_DIR)
  ADD_DEFINITIONS ("-I${CFITSIO_INCLUDE_DIR}")

  IF (CFITSIO_FOUND)
    # Set flags and print a status message
    MESSAGE (STATUS "CFITSIO version ${CFITSIO_VERSION} found:")

    SET (CFITSIO_CPPFLAGS "-I${CFITSIO_INCLUDE_DIR}")
    SET (CFITSIO_LDFLAGS "${CFITSIO_LIBRARIES}")
    
    MESSAGE (STATUS "  * includes: ${CFITSIO_INCLUDE_DIR}")
    MESSAGE (STATUS "  * libs:     ${CFITSIO_LIBRARIES}")
  ENDIF (CFITSIO_FOUND)

ENDIF (NOT CFITSIO_FOUND)

