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

include(FindPackageHandleStandardArgs)

set(CFITSIO_FIND_QUIETLY TRUE)
set(CFITSIO_FIND_REQUIRED TRUE)

if (NOT CFITSIO_FOUND)
    # Manually parse CPLUS_INCLUDE_PATH to add paths to search
    if (DEFINED ENV{CPLUS_INCLUDE_PATH})
        string(REPLACE ":" ";" CPLUS_INCLUDE_PATH_LIST "$ENV{CPLUS_INCLUDE_PATH}")
    else()
        set(CPLUS_INCLUDE_PATH_LIST "")
    endif()

    # Search user environment for headers, then default paths; extract version
    find_path(CFITSIO_INCLUDE_DIR
        NAMES fitsio.h
        PATHS
        $ENV{CFITSIOROOT}/include
        ${CPLUS_INCLUDE_PATH_LIST}
        ${CMAKE_INSTALL_PREFIX}/include
        ${CMAKE_PREFIX_PATH}/include
        /usr/include
        /usr/local/include
        /opt/local/include
        NO_DEFAULT_PATH
    )

    if (NOT CFITSIO_INCLUDE_DIR)
        find_path(CFITSIO_INCLUDE_DIR cfitsio/fitsio.h
            PATHS $ENV{CFITSIOROOT}/include)
    endif()

    if (CFITSIO_INCLUDE_DIR AND EXISTS "${CFITSIO_INCLUDE_DIR}/fitsio.h")
        get_filename_component(CFITSIOROOT ${CFITSIO_INCLUDE_DIR} PATH)
        set(CFITSIO_VERSION 0)
        file(STRINGS "${CFITSIO_INCLUDE_DIR}/fitsio.h" _cfitsio_VERSION REGEX "#define CFITSIO_VERSION[ \t]+([0-9]+\.[0-9]+)")
        string(REGEX REPLACE ".*#define CFITSIO_VERSION[ \t]+([0-9]+\.[0-9]+).*" "\\1" CFITSIO_VERSION "${_cfitsio_VERSION}")
    else()
        set(CFITSIO_INCLUDE_DIR "CFITSIO_INCLUDE_DIR-NOTFOUND")
    endif()

    if (DEFINED ENV{LD_LIBRARY_PATH})
        string(REPLACE ":" ";" LIBRARY_PATH_LIST "$ENV{LD_LIBRARY_PATH}")
    else()
        set(LIBRARY_PATH_LIST "")
    endif()

    # Search user environment for libraries, then default paths
    find_library(CFITSIO_LIBRARIES
        NAMES cfitsio
        PATHS
        $ENV{CFITSIOROOT}/lib
        ${LIBRARY_PATH_LIST}
        ${CMAKE_INSTALL_PREFIX}/lib
        ${CMAKE_PREFIX_PATH}/lib
        /usr/lib
        /usr/local/lib
        /opt/local/lib
        NO_DEFAULT_PATH
    )

    if (CFITSIO_LIBRARIES)
        get_filename_component(CFITSIO_LIB_DIR ${CFITSIO_LIBRARIES} PATH)
    else()
        set(CFITSIO_LIBRARIES "CFITSIO_LIBRARIES-NOTFOUND")
    endif()

    # Set CFITSIO_FOUND and error out if cfitsio is not found
    find_package_handle_standard_args(CFITSIO
        REQUIRED_VARS CFITSIO_LIBRARIES CFITSIO_INCLUDE_DIR
        VERSION_VAR CFITSIO_VERSION
    )

    if (CFITSIO_FOUND)
        # Set flags and print a status message
        message(STATUS "CFITSIO version ${CFITSIO_VERSION} found:")

        set(CFITSIO_CPPFLAGS "-I${CFITSIO_INCLUDE_DIR}")
        set(CFITSIO_LDFLAGS "${CFITSIO_LIBRARIES}")

        message(STATUS "  * includes: ${CFITSIO_INCLUDE_DIR}")
        message(STATUS "  * libs:     ${CFITSIO_LIBRARIES}")
    else()
        message(WARNING "CFITSIO not found. Please ensure CFITSIO is installed and the environment variables CFITSIOROOT, CPLUS_INCLUDE_PATH, LIBRARY_PATH, and LD_LIBRARY_PATH are set correctly.")
    endif()
endif()

