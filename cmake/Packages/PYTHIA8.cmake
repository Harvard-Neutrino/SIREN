# cmake/Packages/PYTHIA8.cmake
#
# Optional Pythia8 (+ LHAPDF) hookup, modeled on MARLEY.cmake.
#
# Usage (top-level CMakeLists.txt):
#   option(SIREN_WITH_PYTHIA8 "Enable Pythia8 support if available" OFF)
#   option(SIREN_REQUIRE_PYTHIA8 "Fail configure if Pythia8 requested but not found" OFF)
#   include(cmake/Packages/PYTHIA8.cmake)
#   if(PYTHIA8_FOUND)
#     target_link_libraries(<your_target> PUBLIC PYTHIA8)
#     target_compile_definitions(<your_target> PUBLIC SIREN_HAS_PYTHIA8=1)
#   endif()

set(PYTHIA8_FOUND FALSE)

if(DEFINED SIREN_WITH_PYTHIA8 AND NOT SIREN_WITH_PYTHIA8)
  message(STATUS "Pythia8 support disabled (SIREN_WITH_PYTHIA8=OFF)")
  return()
endif()

set(_pythia8_ok TRUE)

# -----------------------------
# 1) Find Pythia8 installation prefix
# -----------------------------
unset(PYTHIA8_PREFIX CACHE)
unset(PYTHIA8_PREFIX)

if(DEFINED ENV{PYTHIA8_ROOT} AND NOT "$ENV{PYTHIA8_ROOT}" STREQUAL "")
  set(PYTHIA8_PREFIX "$ENV{PYTHIA8_ROOT}")
  message(STATUS "Using PYTHIA8_ROOT from environment: ${PYTHIA8_PREFIX}")
elseif(DEFINED PYTHIA8_DIR AND NOT "${PYTHIA8_DIR}" STREQUAL "")
  set(PYTHIA8_PREFIX "${PYTHIA8_DIR}")
  message(STATUS "Using PYTHIA8_DIR from cache/cmdline: ${PYTHIA8_PREFIX}")
else()
  find_path(PYTHIA8_PREFIX
    NAMES include/Pythia8/Pythia.h
    PATHS ${CMAKE_PREFIX_PATH} /usr /usr/local /opt/local /opt/homebrew
    DOC "Root directory for Pythia8 installation"
  )
endif()

if(NOT PYTHIA8_PREFIX)
  set(_pythia8_ok FALSE)
  message(STATUS "Pythia8 not found (no PYTHIA8_ROOT/PYTHIA8_DIR and not in CMAKE_PREFIX_PATH).")
endif()

# -----------------------------
# 2) Headers
# -----------------------------
if(_pythia8_ok)
  set(PYTHIA8_INCLUDE_DIR "${PYTHIA8_PREFIX}/include")
  if(NOT EXISTS "${PYTHIA8_INCLUDE_DIR}/Pythia8/Pythia.h")
    set(_pythia8_ok FALSE)
    message(STATUS "Pythia8 headers not found at: ${PYTHIA8_INCLUDE_DIR}/Pythia8/Pythia.h")
  endif()
endif()

# -----------------------------
# 3) Library
# -----------------------------
if(_pythia8_ok)
  find_library(PYTHIA8_LIBRARY
    NAMES pythia8 Pythia8
    PATHS "${PYTHIA8_PREFIX}/lib64" "${PYTHIA8_PREFIX}/lib"
    NO_DEFAULT_PATH
    DOC "Pythia8 library"
  )

  if(NOT PYTHIA8_LIBRARY)
    set(_pythia8_ok FALSE)
    message(STATUS "Pythia8 library not found under ${PYTHIA8_PREFIX}/lib{,64}.")
  endif()
endif()

# -----------------------------
# 4) Optional LHAPDF (PythiaDISCrossSection uses LHAPDF6 PDF set)
# -----------------------------
set(_lhapdf_ok TRUE)
unset(LHAPDF_PREFIX CACHE)
unset(LHAPDF_PREFIX)

if(_pythia8_ok)
  if(DEFINED ENV{LHAPDF_ROOT} AND NOT "$ENV{LHAPDF_ROOT}" STREQUAL "")
    set(LHAPDF_PREFIX "$ENV{LHAPDF_ROOT}")
    message(STATUS "Using LHAPDF_ROOT from environment: ${LHAPDF_PREFIX}")
  elseif(DEFINED LHAPDF_DIR AND NOT "${LHAPDF_DIR}" STREQUAL "")
    set(LHAPDF_PREFIX "${LHAPDF_DIR}")
    message(STATUS "Using LHAPDF_DIR from cache/cmdline: ${LHAPDF_PREFIX}")
  else()
    find_path(LHAPDF_PREFIX
      NAMES include/LHAPDF/LHAPDF.h
      PATHS ${CMAKE_PREFIX_PATH} /usr /usr/local /opt/local /opt/homebrew
      DOC "Root directory for LHAPDF installation"
    )
  endif()

  if(NOT LHAPDF_PREFIX)
    set(_lhapdf_ok FALSE)
    message(STATUS "LHAPDF not found (no LHAPDF_ROOT/LHAPDF_DIR and not in CMAKE_PREFIX_PATH); required for Pythia8 DIS support.")
  else()
    set(LHAPDF_INCLUDE_DIR "${LHAPDF_PREFIX}/include")
    find_library(LHAPDF_LIBRARY
      NAMES LHAPDF
      PATHS "${LHAPDF_PREFIX}/lib64" "${LHAPDF_PREFIX}/lib"
      NO_DEFAULT_PATH
      DOC "LHAPDF library"
    )
    if(NOT LHAPDF_LIBRARY)
      set(_lhapdf_ok FALSE)
      message(STATUS "LHAPDF library not found under ${LHAPDF_PREFIX}/lib{,64}.")
    endif()
  endif()

  if(NOT _lhapdf_ok)
    set(_pythia8_ok FALSE)
  endif()
endif()

# -----------------------------
# 5) Create imported target "PYTHIA8" (carries Pythia8 + LHAPDF link/include)
# -----------------------------
if(_pythia8_ok)
  if(NOT TARGET PYTHIA8)
    add_library(PYTHIA8 UNKNOWN IMPORTED)
    message(STATUS "Created imported target: PYTHIA8")
  endif()

  set_target_properties(PYTHIA8 PROPERTIES
    IMPORTED_LOCATION "${PYTHIA8_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INCLUDE_DIR};${LHAPDF_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${LHAPDF_LIBRARY};dl;pthread"
  )

  set(PYTHIA8_FOUND TRUE)
  message(STATUS "Pythia8 enabled: prefix=${PYTHIA8_PREFIX} (LHAPDF prefix=${LHAPDF_PREFIX})")
endif()

# -----------------------------
# 6) Optional hard requirement behavior
# -----------------------------
if(NOT PYTHIA8_FOUND)
  if(DEFINED SIREN_REQUIRE_PYTHIA8 AND SIREN_REQUIRE_PYTHIA8)
    message(FATAL_ERROR "SIREN_REQUIRE_PYTHIA8=ON but Pythia8 could not be found.")
  endif()
endif()
