# MARLEYConfig.cmake
# This script is intended to be included by a CMakeLists.txt file.
# It locates and configures MARLEY for the current project using standard system locations
# or environment variables such as CMAKE_PREFIX_PATH and MARLEY_ROOT.

# Specify the names of the MARLEY targets
set(MARLEY_TARGETS MARLEY)

# Search for the MARLEY installation prefix
# If MARLEY_ROOT is set, prefer it, otherwise use CMAKE_PREFIX_PATH or default system locations
if(DEFINED ENV{MARLEY_ROOT})
    set(MARLEY_PREFIX "$ENV{MARLEY_ROOT}")
    message(STATUS "Using MARLEY_ROOT from environment: ${MARLEY_PREFIX}")
else()
    find_path(MARLEY_PREFIX
        NAMES include/marley/Particle.hh
        PATHS ${CMAKE_PREFIX_PATH} /usr /usr/local /opt/local /opt/homebrew
        DOC "Root directory for MARLEY installation"
    )
    if(NOT MARLEY_PREFIX)
        message(FATAL_ERROR "MARLEY not found. Please set MARLEY_ROOT or ensure MARLEY is installed in a standard location.")
    endif()
endif()

# Locate the MARLEY config file
set(MARLEY_CONFIG_DIR "${MARLEY_PREFIX}/share/marley/cmake")

if (EXISTS "${MARLEY_CONFIG_DIR}/marleyConfig.cmake")
    message(STATUS "MARLEY configuration file found at: ${MARLEY_CONFIG_DIR}")

    # Include the MARLEY config file
    include("${MARLEY_CONFIG_DIR}/marleyConfig.cmake")
else()
    message(FATAL_ERROR "MARLEY configuration file not found. Expected location: ${MARLEY_CONFIG_DIR}")
endif()

message(STATUS "MARLEY_PREFIX: ${MARLEY_PREFIX}")

# Locate the MARLEY library
find_library(MARLEY_LIBRARY
    NAMES libMARLEY libMARLEY.so libMARLEY${CMAKE_SHARED_LIBRARY_SUFFIX}
    PATHS "${MARLEY_PREFIX}/lib64" "${MARLEY_PREFIX}/lib"
    DOC "MARLEY library"
)

if (MARLEY_LIBRARY)
    message(STATUS "Found MARLEY library at: ${MARLEY_LIBRARY}")
else()
    message(FATAL_ERROR "MARLEY library not found in ${MARLEY_PREFIX}")
endif()

# Set MARLEY include directories
set(MARLEY_INCLUDE_DIR "${MARLEY_PREFIX}/include/marley")

if (EXISTS "${MARLEY_INCLUDE_DIR}")
    message(STATUS "MARLEY include directory found at: ${MARLEY_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "MARLEY include directory not found at: ${MARLEY_INCLUDE_DIR}")
endif()

# Create imported MARLEY targets if they do not already exist
foreach(target ${MARLEY_TARGETS})
    if (NOT TARGET ${target})
        add_library(${target} UNKNOWN IMPORTED)
        message(STATUS "Created imported target for ${target}")
    endif()
endforeach()

# Set properties for the MARLEY targets
set_target_properties(MARLEY PROPERTIES
    IMPORTED_LOCATION "${MARLEY_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MARLEY_INCLUDE_DIR}"
)

