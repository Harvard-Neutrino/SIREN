# cmake/Packages/MARLEY.cmake
#
# Optional MARLEY hookup.
#
# Usage (top-level CMakeLists.txt):
#   option(SIREN_WITH_MARLEY "Enable MARLEY support if available" ON)
#   option(SIREN_REQUIRE_MARLEY "Fail configure if MARLEY requested but not found" OFF)
#   include(cmake/Packages/MARLEY.cmake)
#   if(MARLEY_FOUND)
#     target_link_libraries(<your_target> PUBLIC MARLEY)
#     target_compile_definitions(<your_target> PUBLIC SIREN_HAS_MARLEY=1)
#   endif()

set(MARLEY_FOUND FALSE)

if(DEFINED SIREN_WITH_MARLEY AND NOT SIREN_WITH_MARLEY)
  message(STATUS "MARLEY support disabled (SIREN_WITH_MARLEY=OFF)")
  return()
endif()

set(_marley_ok TRUE)

# -----------------------------
# 1) Find installation prefix
# -----------------------------
unset(MARLEY_PREFIX CACHE)
unset(MARLEY_PREFIX)

if(DEFINED ENV{MARLEY_ROOT} AND NOT "$ENV{MARLEY_ROOT}" STREQUAL "")
  set(MARLEY_PREFIX "$ENV{MARLEY_ROOT}")
  message(STATUS "Using MARLEY_ROOT from environment: ${MARLEY_PREFIX}")
else()
  find_path(MARLEY_PREFIX
    NAMES include/marley/Particle.hh
    PATHS ${CMAKE_PREFIX_PATH} /usr /usr/local /opt/local /opt/homebrew
    DOC "Root directory for MARLEY installation"
  )
endif()

if(NOT MARLEY_PREFIX)
  set(_marley_ok FALSE)
  message(STATUS "MARLEY not found (no MARLEY_ROOT and not in CMAKE_PREFIX_PATH).")
endif()

# -----------------------------
# 2) Locate marleyConfig.cmake (optional; include if present)
# -----------------------------
if(_marley_ok)
  find_path(MARLEY_CONFIG_DIR
    NAMES marleyConfig.cmake
    HINTS "${MARLEY_PREFIX}"
    PATH_SUFFIXES
      share/marley/cmake
      lib/cmake/marley
      lib64/cmake/marley
      cmake
    DOC "Directory containing marleyConfig.cmake"
  )

  if(MARLEY_CONFIG_DIR AND EXISTS "${MARLEY_CONFIG_DIR}/marleyConfig.cmake")
    message(STATUS "MARLEY configuration file found at: ${MARLEY_CONFIG_DIR}/marleyConfig.cmake")
    include("${MARLEY_CONFIG_DIR}/marleyConfig.cmake")
  else()
    # Not fatal: config might be absent even if lib/headers exist.
    message(STATUS "MARLEY config not found (looked for marleyConfig.cmake under ${MARLEY_PREFIX}); will try manual import.")
  endif()
endif()

# -----------------------------
# 3) Find library + headers
# -----------------------------
if(_marley_ok)
  find_library(MARLEY_LIBRARY
    NAMES MARLEY
    PATHS "${MARLEY_PREFIX}/lib64" "${MARLEY_PREFIX}/lib"
    DOC "MARLEY library"
  )

  if(NOT MARLEY_LIBRARY)
    set(_marley_ok FALSE)
    message(STATUS "MARLEY library not found under ${MARLEY_PREFIX}/lib{,64}.")
  endif()
endif()

# NOTE: Your current file sets include dir to "${prefix}/include/marley".
# That only works if you include headers as <Particle.hh>.
# If you include as <marley/Particle.hh> (typical), the include dir should be "${prefix}/include".
if(_marley_ok)
  set(MARLEY_INCLUDE_DIR "${MARLEY_PREFIX}/include")
  if(NOT EXISTS "${MARLEY_INCLUDE_DIR}/marley/Particle.hh")
    set(_marley_ok FALSE)
    message(STATUS "MARLEY headers not found at: ${MARLEY_INCLUDE_DIR}/marley/Particle.hh")
  endif()
endif()

# -----------------------------
# 4) Create/normalize imported target "MARLEY"
# -----------------------------
if(_marley_ok)
  if(NOT TARGET MARLEY)
    add_library(MARLEY UNKNOWN IMPORTED)
    message(STATUS "Created imported target: MARLEY")
  endif()

  set_target_properties(MARLEY PROPERTIES
    IMPORTED_LOCATION "${MARLEY_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MARLEY_INCLUDE_DIR}"
  )

  set(MARLEY_FOUND TRUE)
  message(STATUS "MARLEY enabled: prefix=${MARLEY_PREFIX}")
endif()

# -----------------------------
# 5) Optional hard requirement behavior
# -----------------------------
if(NOT MARLEY_FOUND)
  if(DEFINED SIREN_REQUIRE_MARLEY AND SIREN_REQUIRE_MARLEY)
    message(FATAL_ERROR "SIREN_REQUIRE_MARLEY=ON but MARLEY could not be found.")
  endif()
endif()

