# cmake/Packages/ACHILLES.cmake
#
# Optional Achilles hookup, for AchillesTrialGenerator. Achilles must expose the
# geometry-driver interface (EventGen::InjectRay, GenerateSingleEvent,
# VertexEnvelope, TotalXSec), as on AchillesGen/Achilles branch
# feature/nugeometry-interface.
#
#   -DSIREN_WITH_ACHILLES=ON -DCMAKE_PREFIX_PATH=<achilles install prefix>
#
# Put the Achilles install first on CMAKE_PREFIX_PATH so that its own fmt,
# spdlog and yaml-cpp are found rather than other installations. Achilles'
# HepMC3 must be the same version as SIREN's; both are loaded as
# @rpath/libHepMC3.4.dylib and share one image.

set(ACHILLES_FOUND FALSE)

if(NOT SIREN_WITH_ACHILLES)
  message(STATUS "Achilles support disabled (SIREN_WITH_ACHILLES=OFF)")
  return()
endif()

find_package(Achilles CONFIG QUIET)

if(Achilles_FOUND AND TARGET achilles::Achilles)
  add_library(ACHILLES INTERFACE IMPORTED GLOBAL)
  target_link_libraries(ACHILLES INTERFACE achilles::Achilles)
  target_compile_definitions(ACHILLES INTERFACE SIREN_HAS_ACHILLES)
  set(ACHILLES_FOUND TRUE)
  message(STATUS "Achilles found: ${Achilles_DIR}")
elseif(SIREN_REQUIRE_ACHILLES)
  message(FATAL_ERROR "SIREN_REQUIRE_ACHILLES=ON but Achilles (achilles::Achilles) was not found; "
                      "add its install prefix to CMAKE_PREFIX_PATH")
else()
  message(STATUS "Achilles not found; AchillesTrialGenerator is not built")
endif()
