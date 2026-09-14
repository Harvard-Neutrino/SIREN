cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
cmake_policy(VERSION 3.20)

if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_POLICY_DEFAULT_CMP0042 NEW)
    cmake_policy(SET CMP0042 NEW)
    set(CMAKE_POLICY_DEFAULT_CMP0068 NEW)
    cmake_policy(SET CMP0068 NEW)
    set(MACOSX TRUE)
endif()

# point to appropriate cmake file locations
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_CURRENT_LIST_DIR}/cmake/Packages)

message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

set(CIBUILDWHEEL "$ENV{CIBUILDWHEEL}")
message(STATUS "CIBUILDWHEEL: ${CIBUILDWHEEL}")
if(CIBUILDWHEEL STREQUAL "1")
    set(CI_INSTALL_PREFIX "$ENV{CI_INSTALL_PREFIX}")
    message(STATUS "CI_INSTALL_PREFIX: ${CI_INSTALL_PREFIX}")
endif()


# parse pyproject.toml for the version
include(siren_pyproject)

project(siren VERSION ${SIREN_VERSION} LANGUAGES CXX)

option(SIREN_PYTHON_PACKAGE
    "Build and install the SIREN Python package"
    ON  # default ON, or OFF if you prefer
)
option(SIREN_WITH_MARLEY "Enable MARLEY support if available" ON)
option(SIREN_REQUIRE_MARLEY "Fail configure if MARLEY requested but not found" OFF)
option(SIREN_WITH_PYTHIA8 "Enable Pythia8 support if available" OFF)
option(SIREN_REQUIRE_PYTHIA8 "Fail configure if Pythia8 requested but not found" OFF)
option(SIREN_WITH_HEPMC3 "Enable HepMC3 output support if available" ON)
option(SIREN_REQUIRE_HEPMC3 "Fail configure if HepMC3 requested but not found" OFF)

# include cmake dependencies
include(GNUInstallDirs)
include(ExternalProject)

SET(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_C_STANDARD 99)

# set the build type and appropriate flags
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

function(apply_siren_compile_options target)
    # Check if the target is an INTERFACE target
    get_target_property(target_type ${target} TYPE)
    get_target_property(target_imported ${target} IMPORTED)

    if("${target_type}" STREQUAL "INTERFACE_LIBRARY")
        # Apply compile options to INTERFACE target
        target_compile_options(${target} INTERFACE
            -O2
            -Wall
            -fPIC
            $<$<CONFIG:Debug>:-g>
            $<$<CONFIG:Debug>:-O0>
            $<$<CONFIG:Release>:-O2>
            $<$<AND:$<COMPILE_LANG_AND_ID:CXX,GNU>,$<CONFIG:Release>>:-s>
            $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-stdlib=libc++>
        )
    elseif(target_imported)
        # Apply compile options to non-INTERFACE target (PRIVATE or PUBLIC)
        target_compile_options(${target} INTERFACE
            -O2
            -Wall
            -fPIC
            $<$<CONFIG:Debug>:-g>
            $<$<CONFIG:Debug>:-O0>
            $<$<CONFIG:Release>:-O2>
            $<$<AND:$<COMPILE_LANG_AND_ID:CXX,GNU>,$<CONFIG:Release>>:-s>
            $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-stdlib=libc++>
        )
    else()
        # Apply compile options to non-INTERFACE target (PRIVATE or PUBLIC)
        target_compile_options(${target} PRIVATE
            -O2
            -Wall
            -fPIC
            $<$<CONFIG:Debug>:-g>
            $<$<CONFIG:Debug>:-O0>
            $<$<CONFIG:Release>:-O2>
            $<$<AND:$<COMPILE_LANG_AND_ID:CXX,GNU>,$<CONFIG:Release>>:-s>
            $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-stdlib=libc++>
        )
    endif()
endfunction()

# Helper: strip local symbols from a shared object after build, in non-Debug
# configurations. Safe for both single- and multi-config generators:
#   * Single-config (Make, Ninja): gates on CMAKE_BUILD_TYPE at configure time,
#     because CMAKE_BUILD_TYPE is fixed.
#   * Multi-config (Xcode, Visual Studio, Ninja Multi-Config): gates on
#     $<CONFIG:Debug> at build time, because CMAKE_BUILD_TYPE is empty and
#     the actual config is only known at build time. When the config is
#     Debug, every generator expression evaluates to an empty string, the
#     COMMAND tokens are all empty, and CMake skips the rule.
# Skips silently if CMAKE_STRIP is not available on the host toolchain.
# CMake escapes the ELF origin token when writing install rules.
if(APPLE)
    set(SIREN_RPATH_ORIGIN "@loader_path")
else()
    set(SIREN_RPATH_ORIGIN "$ORIGIN")
endif()

function(siren_strip_module target)
    if(NOT CMAKE_STRIP)
        return()
    endif()
    get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if(_is_multi_config)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND $<$<NOT:$<CONFIG:Debug>>:${CMAKE_STRIP}>
                    $<$<NOT:$<CONFIG:Debug>>:-x>
                    $<$<NOT:$<CONFIG:Debug>>:$<TARGET_FILE:${target}>>
            COMMAND_EXPAND_LISTS
            COMMENT "Stripping local symbols from ${target} (release configs only)")
    elseif(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_STRIP} -x $<TARGET_FILE:${target}>
            COMMENT "Stripping local symbols from ${target}")
    endif()
endfunction()


# override install locations when building python extensions
if(DEFINED SKBUILD_PLATLIB_DIR)
    cmake_path(RELATIVE_PATH SKBUILD_HEADERS_DIR BASE_DIRECTORY ${SKBUILD_PLATLIB_DIR} OUTPUT_VARIABLE CMAKE_INSTALL_INCLUDEDIR)
    cmake_path(RELATIVE_PATH SKBUILD_PLATLIB_DIR BASE_DIRECTORY ${SKBUILD_PLATLIB_DIR} OUTPUT_VARIABLE CMAKE_INSTALL_LIBDIR)
    message(STATUS "Setting include dir to: ${CMAKE_INSTALL_INCLUDEDIR}")
    message(STATUS "Setting lib dir to: ${CMAKE_INSTALL_LIBDIR}")
    set(SIREN_INSTALL_INCLUDEDIR ${SKBUILD_HEADERS_DIR})
    set(SIREN_INSTALL_LIBDIR ${SKBUILD_PLATLIB_DIR})
else()
    set(SIREN_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
    set(SIREN_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()

# load python and pybind11
include(Python)
include(pybind11)

# load project dependencies
include(rk)
include(crundec3)
if(TARGET rk_static)
    apply_siren_compile_options(rk_static)
endif()
if(TARGET rk_shared)
    apply_siren_compile_options(rk_shared)
endif()
if(TARGET crundec_static)
    apply_siren_compile_options(crundec_static)
endif()
include(cereal)
if(TARGET cereal)
    apply_siren_compile_options(cereal)
endif()
set(DELABELLA_BUILD_STATIC ON CACHE BOOL "Build delabella static library" FORCE)
set(DELABELLA_BUILD_SHARED OFF CACHE BOOL "Build delabella shared library" FORCE)
include(delabella)
if(TARGET delabella_static)
    apply_siren_compile_options(delabella_static)
endif()
include(CFITSIO)
if(TARGET CFITSIO::cfitsio)
    apply_siren_compile_options(CFITSIO::cfitsio)
endif()
include(photospline)
if(TARGET photospline)
    apply_siren_compile_options(photospline)
    target_link_libraries(photospline PRIVATE CFITSIO::cfitsio)
endif()
include(googletest)
if(TARGET gtest)
    apply_siren_compile_options(gtest)
endif()
if(TARGET gtest_main)
    apply_siren_compile_options(gtest_main)
endif()
if(TARGET gmock)
    apply_siren_compile_options(gmock)
endif()
include(NamedType)
if(TARGET NamedType)
    apply_siren_compile_options(NamedType)
endif()
include(MARLEY)
if(TARGET MARLEY)
    apply_siren_compile_options(MARLEY)
endif()
include(PYTHIA8)
if(TARGET PYTHIA8)
    apply_siren_compile_options(PYTHIA8)
endif()
include(HepMC3)

# load macros for googletest
include(testing)

# add individual source folders
add_subdirectory(projects/utilities)
add_subdirectory(projects/serialization)
add_subdirectory(projects/math)
add_subdirectory(projects/dataclasses)
add_subdirectory(projects/geometry)
add_subdirectory(projects/detector)
add_subdirectory(projects/interactions)
add_subdirectory(projects/distributions)
add_subdirectory(projects/injection)
add_subdirectory(projects/io)

# define the target library
add_library(SIREN SHARED)
set_property(TARGET SIREN PROPERTY POSITION_INDEPENDENT_CODE ON)
apply_siren_compile_options(SIREN)
siren_strip_module(SIREN)

if(DEFINED MACOSX AND MACOSX)
if(CMAKE_VERSION VERSION_LESS 3.13)
    target_link_libraries(SIREN PUBLIC "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")
else()
    target_link_options(SIREN PUBLIC "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")
endif()
target_link_libraries(SIREN
    PRIVATE
    PUBLIC
        photospline
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_geometry
        SIREN_detector
        SIREN_interactions
        SIREN_distributions
        SIREN_injection
        SIREN_io
)
else()
target_link_libraries(SIREN
    PRIVATE
    PUBLIC
        photospline
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_geometry
        SIREN_detector
        SIREN_interactions
        SIREN_distributions
        SIREN_injection
        SIREN_io
)
endif()

if(PACKAGE_SHARED_DEPS)
    message(FATAL_ERROR
        "PACKAGE_SHARED_DEPS cannot produce a relocatable wheel. Use delocate-wheel, "
        "auditwheel repair, or delvewheel repair as in the cibuildwheel configuration.")
endif()

# Extension modules resolve Python through their host interpreter. Only C++
# executables embed Python; linking libpython into the wheel can load a second
# interpreter when wheel repair bundles that library.
target_link_libraries(SIREN INTERFACE
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:pybind11::embed>")

# Native development installation. Wheel builds select only PythonWheel.
install(TARGETS SIREN
    SIREN_utilities SIREN_serialization SIREN_math SIREN_dataclasses
    SIREN_geometry SIREN_detector SIREN_interactions SIREN_distributions
    SIREN_injection SIREN_io
    EXPORT ${PROJECT_NAME}Config
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Both wheel entry points install the targets, allowing CMake to replace build
# RPATHs. The core library and photospline share a package-local directory.
set(SIREN_PYTHON_MODULES utilities math dataclasses geometry detector
    interactions distributions injection hepmc3)
set(SIREN_WHEEL_LIBRARIES SIREN photospline)
if(TARGET spglam)
    list(APPEND SIREN_WHEEL_LIBRARIES spglam)
endif()
set_target_properties(${SIREN_WHEEL_LIBRARIES} PROPERTIES
    BUILD_WITH_INSTALL_RPATH FALSE
    INSTALL_RPATH "${SIREN_RPATH_ORIGIN}"
    INSTALL_RPATH_USE_LINK_PATH FALSE)
if(APPLE)
    set_target_properties(${SIREN_WHEEL_LIBRARIES} PROPERTIES
        INSTALL_NAME_DIR "@rpath"
        BUILD_WITH_INSTALL_NAME_DIR TRUE)
endif()
foreach(module IN LISTS SIREN_PYTHON_MODULES)
    set_target_properties(${module} PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH "${SIREN_RPATH_ORIGIN}/../siren.libs"
        INSTALL_RPATH_USE_LINK_PATH FALSE)
endforeach()
install(TARGETS ${SIREN_WHEEL_LIBRARIES}
    LIBRARY DESTINATION siren.libs COMPONENT PythonWheel EXCLUDE_FROM_ALL
    RUNTIME DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
install(TARGETS ${SIREN_PYTHON_MODULES}
    LIBRARY DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL
    RUNTIME DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/
    DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL
    PATTERN "__pycache__" EXCLUDE PATTERN "*.pyc" EXCLUDE)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/resources
    DESTINATION siren COMPONENT PythonWheel EXCLUDE_FROM_ALL
    PATTERN "__pycache__" EXCLUDE PATTERN "*.pyc" EXCLUDE)

include(siren_python_package)


# Export targets for use in downstream CMake projects
# -----------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
# Make importable from build directory
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT ${PROJECT_NAME}Config FILE "${PROJECT_NAME}Config.cmake"
       NAMESPACE ${PROJECT_NAME}::)

# Make importable from install location
set(_config_dir share/${PROJECT_NAME}/cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${_config_dir}
)
install(EXPORT ${PROJECT_NAME}Config
  FILE "${PROJECT_NAME}Config.cmake"
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION "${_config_dir}"
)

MESSAGE("")
MESSAGE("Done!")
MESSAGE("Run `make` then `make install`")
