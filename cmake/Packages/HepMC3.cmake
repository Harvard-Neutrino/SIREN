# cmake/Packages/HepMC3.cmake
#
# Optional HepMC3 hookup (ROOT-free). Prefers the upstream HepMC3Config.cmake,
# falling back to a manual find. Sets HEPMC3_FOUND and ensures a usable
# HepMC3::HepMC3 imported target.
#
# Usage (top-level CMakeLists.txt):
#   option(SIREN_WITH_HEPMC3 "Enable HepMC3 output support if available" ON)
#   option(SIREN_REQUIRE_HEPMC3 "Fail configure if HepMC3 requested but not found" OFF)
#   include(HepMC3)
#   if(HEPMC3_FOUND)
#     target_link_libraries(<target> PUBLIC HepMC3::HepMC3)
#     target_compile_definitions(<target> PUBLIC SIREN_HAS_HEPMC3=1)
#   endif()

set(HEPMC3_FOUND FALSE)

if(DEFINED SIREN_WITH_HEPMC3 AND NOT SIREN_WITH_HEPMC3)
  message(STATUS "HepMC3 support disabled (SIREN_WITH_HEPMC3=OFF)")
  return()
endif()

# 1) Preferred: upstream package config (defines HepMC3::HepMC3).
find_package(HepMC3 CONFIG QUIET)
if(HepMC3_FOUND AND TARGET HepMC3::HepMC3)
  set(HEPMC3_FOUND TRUE)
  message(STATUS "Found HepMC3 via config: ${HepMC3_DIR}")
  return()
endif()

# 2) Fallback: manual find + UNKNOWN IMPORTED target.
find_path(HEPMC3_INCLUDE_DIR
  NAMES HepMC3/GenEvent.h
  HINTS $ENV{HEPMC3_ROOT}
  PATHS ${CMAKE_PREFIX_PATH} /usr /usr/local /opt/local /opt/homebrew
  PATH_SUFFIXES include
  DOC "HepMC3 include directory")

find_library(HEPMC3_LIBRARY
  NAMES HepMC3
  HINTS $ENV{HEPMC3_ROOT}
  PATHS ${CMAKE_PREFIX_PATH} /usr /usr/local /opt/local /opt/homebrew
  PATH_SUFFIXES lib lib64
  DOC "HepMC3 core library")

if(HEPMC3_INCLUDE_DIR AND HEPMC3_LIBRARY)
  if(NOT TARGET HepMC3::HepMC3)
    add_library(HepMC3::HepMC3 UNKNOWN IMPORTED)
    set_target_properties(HepMC3::HepMC3 PROPERTIES
      IMPORTED_LOCATION "${HEPMC3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${HEPMC3_INCLUDE_DIR}")
  endif()
  set(HEPMC3_FOUND TRUE)
  message(STATUS "Found HepMC3 (manual): ${HEPMC3_LIBRARY}")
endif()

if(NOT HEPMC3_FOUND)
  if(SIREN_REQUIRE_HEPMC3)
    message(FATAL_ERROR "HepMC3 requested (SIREN_REQUIRE_HEPMC3=ON) but not found")
  else()
    message(STATUS "HepMC3 not found; siren.io builds without HepMC3 support")
  endif()
endif()
