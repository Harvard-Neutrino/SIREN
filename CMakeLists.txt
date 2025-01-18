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
include(pyproject)

project(siren VERSION ${SIREN_VERSION} LANGUAGES CXX)

# include cmake dependencies
include(GNUInstallDirs)
include(ExternalProject)

SET(CMAKE_CXX_STANDARD 14)
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
if(TARGET rk_static)
    apply_siren_compile_options(rk_static)
endif()
if(TARGET rk_shared)
    apply_siren_compile_options(rk_shared)
endif()
include(cereal)
if(TARGET cereal)
    apply_siren_compile_options(cereal)
endif()
include(delabella)
if(TARGET delabella_static)
    apply_siren_compile_options(delabella_static)
endif()
if(TARGET delabella_shared)
    apply_siren_compile_options(delabella_shared)
endif()
include(CFITSIO)
if(TARGET CFITSIO)
    apply_siren_compile_options(CFITSIO)
endif()
include(photospline)
if(TARGET photospline)
    apply_siren_compile_options(photospline)
    target_link_libraries(photospline PUBLIC CFITSIO)
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

# define the target library
add_library(SIREN SHARED)
set_property(TARGET SIREN PROPERTY POSITION_INDEPENDENT_CODE ON)
apply_siren_compile_options(SIREN)

if(DEFINED MACOSX AND MACOSX)
if(CMAKE_VERSION VERSION_LESS 3.13)
    target_link_libraries(SIREN PUBLIC "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")
else()
    target_link_options(SIREN PUBLIC "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")
endif()
target_link_libraries(SIREN
    PRIVATE
    PUBLIC
        delabella_shared
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
)
else()
target_link_libraries(SIREN
    PRIVATE
    PUBLIC
        photospline
        delabella_shared
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_geometry
        SIREN_detector
        SIREN_interactions
        SIREN_distributions
        SIREN_injection
)
endif()

# define the install path normally or for python package
if(DEFINED SKBUILD_PLATLIB_DIR)
    set_target_properties(SIREN PROPERTIES
            BUILD_WITH_INSTALL_RPATH FALSE
            LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN")
    if(DEFINED CIBUILDWHEEL AND CIBUILDWHEEL)
        message(STATUS "Setting SIREN install lib dir to: ${CI_INSTALL_PREFIX}/lib")
        message(STATUS "Setting SIREN install include dir to: ${CI_INSTALL_PREFIX}/include")
        install(TARGETS SIREN
            delabella_shared
            SIREN_utilities
            SIREN_serialization
            SIREN_math
            SIREN_dataclasses
            SIREN_geometry
            SIREN_detector
            SIREN_interactions
            SIREN_distributions
            SIREN_injection
            EXPORT ${PROJECT_NAME}Config
            LIBRARY DESTINATION "${CI_INSTALL_PREFIX}/lib"
            PUBLIC_HEADER DESTINATION "${CI_INSTALL_PREFIX}/include")
    else()
        install(TARGETS SIREN
            delabella_shared
            SIREN_utilities
            SIREN_serialization
            SIREN_math
            SIREN_dataclasses
            SIREN_geometry
            SIREN_detector
            SIREN_interactions
            SIREN_distributions
            SIREN_injection
            EXPORT ${PROJECT_NAME}Config
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren.libs
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    endif()
else()
    install(TARGETS SIREN
        delabella_shared
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_geometry
        SIREN_detector
        SIREN_interactions
        SIREN_distributions
        SIREN_injection
        EXPORT ${PROJECT_NAME}Config
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# optionally package runtime dependencies
if((DEFINED SKBUILD_PLATLIB_DIR) AND (PACKAGE_SHARED_DEPS))
    install(CODE "set(SIREN_LIB_FILE \"${PROJECT_BINARY_DIR}/${CMAKE_SHARED_MODULE_PREFIX}SIREN${CMAKE_SHARED_MODULE_SUFFIX}\")")
    install(CODE "set(PYTHON_DEP_LIB_DESTINATION \"${SKBUILD_PLATLIB_DIR}/siren.libs/\")")
    install(CODE [[
        file(GET_RUNTIME_DEPENDENCIES
            LIBRARIES ${SIREN_LIB_FILE}
            RESOLVED_DEPENDENCIES_VAR _r_deps
            UNRESOLVED_DEPENDENCIES_VAR _u_deps
            )
        foreach(_file ${_r_deps})
            file(INSTALL
                DESTINATION "${PYTHON_DEP_LIB_DESTINATION}"
                TYPE SHARED_LIBRARY
                FOLLOW_SYMLINK_CHAIN
                FILES "${_file}"
                )
        endforeach()
        list(LENGTH _u_deps _u_length)
        if("${_u_length}" GREATER 0)
            message(WARNING "Unresolved dependencies detected!")
        endif()
        ]]
    )
endif()

# install the python extensions
if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS utilities
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS math
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS dataclasses
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS geometry
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS detector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS interactions
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS distributions
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(TARGETS injection
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/resources
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/siren)
endif()

# Stage the python package within the build directory
set(PACKAGE_STAGING_DIR "${CMAKE_BINARY_DIR}/python_staging")

find_package(Python3 COMPONENTS Interpreter REQUIRED)  # or Python

# Configure script that copies extra files into the package staging area
configure_file("${PROJECT_SOURCE_DIR}/package/configure_python_package.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/configure_python_package.cmake" @ONLY)

# Clear out the staging area
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean"
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PACKAGE_STAGING_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PACKAGE_STAGING_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PACKAGE_STAGING_DIR}/siren
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean
    COMMENT "Clear out python package staging area"
    VERBATIM
)

# Copy shared object libraries into the staging area
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_extensions"
    DEPENDS
      SIREN
      utilities
      math
      dataclasses
      geometry
      detector
      interactions
      distributions
      injection
      "${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:utilities>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:math>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:dataclasses>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:geometry>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:detector>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:interactions>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:distributions>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:injection>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:SIREN>"
        "${PACKAGE_STAGING_DIR}/siren/"
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_extensions
    COMMENT "Copying shared object libraries into python package staging directory"
    VERBATIM
)

# Copy python files into the staging area
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_python
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean
        ${CMAKE_SOURCE_DIR}/python/Injector.py
        ${CMAKE_SOURCE_DIR}/python/Weighter.py
        ${CMAKE_SOURCE_DIR}/python/__init__.py
        ${CMAKE_SOURCE_DIR}/python/_util.py
        ${CMAKE_SOURCE_DIR}/python/resources.py
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/python
        ${PACKAGE_STAGING_DIR}/${PROJECT_NAME}
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_python
    COMMENT "Copying contents of python/ into the package staging directory"
    VERBATIM
)

# Copy resources into the staging area
file(GLOB_RECURSE RESOURCES_FILES ${CMAKE_SOURCE_DIR}/resources/*)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_resources
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean
        ${RESOURCES_FILES}
    COMMAND ${CMAKE_COMMAND} -E make_directory
        ${PACKAGE_STAGING_DIR}/${PROJECT_NAME}/resources
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/resources
        ${PACKAGE_STAGING_DIR}/${PROJECT_NAME}/resources
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_resources
    COMMENT "Copying resources/ into the staging directory"
    VERBATIM
)

# Parse top-level pyproject.toml and generate new pyproject.toml for package staging area
add_custom_command(
    OUTPUT ${PACKAGE_STAGING_DIR}/pyproject.toml
    COMMAND
        ${Python3_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/cmake/parse_pyproject.py
        ${CMAKE_SOURCE_DIR}/pyproject.toml
        ${PACKAGE_STAGING_DIR}/pyproject.toml
    DEPENDS
        ${CMAKE_SOURCE_DIR}/pyproject.toml
        ${CMAKE_SOURCE_DIR}/cmake/parse_pyproject.py
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean
    COMMAND ${CMAKE_COMMAND} -E touch ${PACKAGE_STAGING_DIR}/pyproject.toml
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean
    COMMENT "Parsing top-level pyproject.toml and generating new pyproject.toml for package staging area"
    VERBATIM
)

# Configure other files for the python package
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/.python_package_configured
    COMMAND ${CMAKE_COMMAND} -P
        ${CMAKE_CURRENT_BINARY_DIR}/configure_python_package.cmake
    DEPENDS
        ${PACKAGE_STAGING_DIR}/pyproject.toml
        ${CMAKE_CURRENT_BINARY_DIR}/configure_python_package.cmake
        ${PROJECT_SOURCE_DIR}/package/configure_python_package.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_resources
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_python
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_extensions
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.python_package_configured
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.stamp_clean
    COMMENT "Configuring other files for python package"
    VERBATIM
)

# Build the wheel from the staged package
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/.build_wheel
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/dist_wheels
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/dist_wheels
    COMMAND ${Python3_EXECUTABLE} -m pip wheel
            --no-deps
            --wheel-dir ${CMAKE_CURRENT_BINARY_DIR}/dist_wheels
            ${PACKAGE_STAGING_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/.python_package_configured
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.build_wheel
    COMMENT "Building wheel from the staged package"
)

set(WHEELS_DIR "${CMAKE_CURRENT_BINARY_DIR}/dist_wheels")
install(
    CODE
    "
    file(GLOB WHEELS \"${WHEELS_DIR}/*.whl\")
    message(STATUS \"Installing wheels: \$\{WHEELS\}\")
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install --no-deps --force-reinstall \$\{WHEELS\} --prefix=${CMAKE_INSTALL_PREFIX}
        COMMAND_ECHO STDOUT
    )
    "
)

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
