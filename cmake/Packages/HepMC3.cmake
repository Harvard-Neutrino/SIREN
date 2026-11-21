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

# The NuHepMC writer uses the Vector*Attribute types (VectorDoubleAttribute etc.)
# that HepMC3 3.2.x does not provide; 3.3 is the floor. Older versions are
# treated exactly like not-found.
set(HEPMC3_MINIMUM_VERSION 3.3)

if(DEFINED SIREN_WITH_HEPMC3 AND NOT SIREN_WITH_HEPMC3)
  message(STATUS "HepMC3 support disabled (SIREN_WITH_HEPMC3=OFF)")
  return()
endif()

# 1) Preferred: upstream package config (defines HepMC3::HepMC3). Requesting the
# minimum version makes config mode reject an older HepMC3 (HepMC3_FOUND stays
# false) before we ever adopt the target.
find_package(HepMC3 ${HEPMC3_MINIMUM_VERSION} CONFIG QUIET)
if(HepMC3_FOUND AND TARGET HepMC3::HepMC3)
  set(HEPMC3_FOUND TRUE)
  message(STATUS "Found HepMC3 via config: ${HepMC3_DIR} (version ${HepMC3_VERSION})")
endif()

# 2) Fallback: manual find + UNKNOWN IMPORTED target.
if(NOT HEPMC3_FOUND)
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
    # The manual path has no CMake package version to check, so read the
    # version macro straight out of the discovered headers. HepMC3/Version.h
    # defines HEPMC3_VERSION_CODE as 1000000*X + 1000*Y + Z (e.g. 3003001 for
    # 3.3.1); compare against the same floor used by the config-mode request.
    set(_hepmc3_manual_version_header "${HEPMC3_INCLUDE_DIR}/HepMC3/Version.h")
    set(HEPMC3_MANUAL_VERSION_OK FALSE)
    if(EXISTS "${_hepmc3_manual_version_header}")
      file(STRINGS "${_hepmc3_manual_version_header}" _hepmc3_version_code_line
        REGEX "^#define[ \t]+HEPMC3_VERSION_CODE[ \t]+[0-9]+")
      if(_hepmc3_version_code_line)
        string(REGEX REPLACE
          "^#define[ \t]+HEPMC3_VERSION_CODE[ \t]+([0-9]+).*" "\\1"
          _hepmc3_version_code "${_hepmc3_version_code_line}")
        # Decode 1000000*X + 1000*Y + Z back into X.Y.Z for VERSION_LESS.
        math(EXPR _hepmc3_ver_major "${_hepmc3_version_code} / 1000000")
        math(EXPR _hepmc3_ver_minor "(${_hepmc3_version_code} / 1000) % 1000")
        math(EXPR _hepmc3_ver_patch "${_hepmc3_version_code} % 1000")
        set(_hepmc3_manual_version "${_hepmc3_ver_major}.${_hepmc3_ver_minor}.${_hepmc3_ver_patch}")
        if(NOT _hepmc3_manual_version VERSION_LESS HEPMC3_MINIMUM_VERSION)
          set(HEPMC3_MANUAL_VERSION_OK TRUE)
        endif()
      endif()
    endif()

    if(HEPMC3_MANUAL_VERSION_OK)
      if(NOT TARGET HepMC3::HepMC3)
        add_library(HepMC3::HepMC3 UNKNOWN IMPORTED)
        set_target_properties(HepMC3::HepMC3 PROPERTIES
          IMPORTED_LOCATION "${HEPMC3_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${HEPMC3_INCLUDE_DIR}")
      endif()
      set(HEPMC3_FOUND TRUE)
      set(HepMC3_VERSION "${_hepmc3_manual_version}")
      message(STATUS "Found HepMC3 (manual): ${HEPMC3_LIBRARY} (version ${HepMC3_VERSION})")
    elseif(_hepmc3_manual_version)
      message(STATUS
        "Found HepMC3 ${_hepmc3_manual_version} (manual) at ${HEPMC3_LIBRARY}, "
        "but >= ${HEPMC3_MINIMUM_VERSION} is required for NuHepMC output; "
        "treating as not found")
    else()
      message(STATUS
        "Found a HepMC3 install (manual) at ${HEPMC3_LIBRARY} but could not "
        "determine its version from ${_hepmc3_manual_version_header}; "
        "treating as not found")
    endif()
  endif()
endif()

# 2b) Explicit version floor. Config mode should already have rejected an old
# HepMC3 above, but a config that does not honor the version request would still
# set HepMC3_VERSION; reject it here so we degrade exactly like not-found. The
# manual fallback resolves and checks its own version just above, so this only
# re-guards the config-mode path.
if(HEPMC3_FOUND AND DEFINED HepMC3_VERSION AND HepMC3_VERSION
    AND HepMC3_VERSION VERSION_LESS HEPMC3_MINIMUM_VERSION)
  message(STATUS
    "Found HepMC3 ${HepMC3_VERSION}, but >= ${HEPMC3_MINIMUM_VERSION} is "
    "required for NuHepMC output; treating as not found")
  set(HEPMC3_FOUND FALSE)
endif()

if(NOT HEPMC3_FOUND)
  if(SIREN_REQUIRE_HEPMC3)
    message(FATAL_ERROR "HepMC3 requested (SIREN_REQUIRE_HEPMC3=ON) but not found")
  else()
    message(STATUS "HepMC3 not found; siren.io builds without HepMC3 support")
  endif()
  return()
endif()

# 3) Probe optional zlib compression support (for gzip output). WriterGZ/ReaderGZ
#    only compile+link when HepMC3 was built with compression AND the gate macros
#    are defined; the HepMC3 config does not export them, so probe here and, only
#    if the probe links, let the io target define the macros itself. A failed
#    probe simply disables gzip -- it is never a configure error.
# bxzstr (HepMC3's compression backend) is header-only: the inflate/deflate calls
# are emitted into the CONSUMER'S translation unit and libHepMC3 exports none of
# them, so the probe (and later the io target) must link zlib itself. Without zlib
# the probe would fail to link and gzip would be (correctly) disabled -- but if zlib
# is present we must link it or the real build link-fails.
find_package(ZLIB QUIET)
if(ZLIB_FOUND)
  set(_hepmc3_gz_src "${CMAKE_CURRENT_BINARY_DIR}/hepmc3_gz_probe.cxx")
  file(WRITE "${_hepmc3_gz_src}"
  "#ifndef HEPMC3_USE_COMPRESSION\n#define HEPMC3_USE_COMPRESSION 1\n#endif\n"
  "#ifndef HEPMC3_Z_SUPPORT\n#define HEPMC3_Z_SUPPORT 1\n#endif\n"
  "#include <memory>\n"
  "#include <HepMC3/WriterGZ.h>\n#include <HepMC3/WriterAscii.h>\n"
  "#include <HepMC3/ReaderGZ.h>\n#include <HepMC3/ReaderAscii.h>\n"
  "#include <HepMC3/GenRunInfo.h>\n"
  "int main(){ HepMC3::WriterGZ<HepMC3::WriterAscii> w(\"probe.hepmc3.gz\", std::make_shared<HepMC3::GenRunInfo>());\n"
  "  HepMC3::ReaderGZ<HepMC3::ReaderAscii> r(\"probe.hepmc3.gz\"); return 0; }\n")
  try_compile(HEPMC3_HAS_COMPRESSION
    "${CMAKE_CURRENT_BINARY_DIR}/hepmc3_gz_probe"
    "${_hepmc3_gz_src}"
    LINK_LIBRARIES HepMC3::HepMC3 ZLIB::ZLIB
    CMAKE_FLAGS "-DCMAKE_CXX_STANDARD=17")
else()
  set(HEPMC3_HAS_COMPRESSION FALSE)
endif()
if(HEPMC3_HAS_COMPRESSION)
  message(STATUS "HepMC3 compression (gzip) available")
else()
  message(STATUS "HepMC3 compression not available; gzip output disabled")
endif()
