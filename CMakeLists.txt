cmake_minimum_required(VERSION 3.3.2 FATAL_ERROR)
cmake_policy(VERSION 3.3.2)

if(${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
endif()

# point to appropriate cmake file locations
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_CURRENT_LIST_DIR}/cmake/Packages)

message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

set(CIBUILDWHEEL $ENV{CIBUILDWHEEL})
message(STATUS "CIBUILDWHEEL: ${CIBUILDWHEEL}")
if(${CIBUILDWHEEL})
    set(CI_INSTALL_PREFIX $ENV{CI_INSTALL_PREFIX})
    message(STATUS "CI_INSTALL_PREFIX: ${CI_INSTALL_PREFIX}")
endif()

# parse pyproject.toml for the version
include(pyproject)

project(leptoninjector VERSION ${LI_VERSION} LANGUAGES CXX)

# include cmake dependencies
include(GNUInstallDirs)
include(ExternalProject)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_C_STANDARD 99)

# set the build type and appropriate flags
option(CMAKE_BUILD_TYPE "" "Release")
set(_FLAGS "-O2 -Wall -fPIC")
set(_FLAGS_DEBUG "-g -O0 -Wall -fPIC")
set(_FLAGS_RELEASE "-O2 -Wall -fPIC -s")
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(_FLAGS "${_FLAGS} -stdlib=libc++")
    set(_FLAGS_DEBUG "${_FLAGS_DEBUG} -stdlib=libc++")
    set(_FLAGS_RELEASE "${_FLAGS_RELEASE} -stdlib=libc++")
endif()
set(CMAKE_CXX_FLAGS ${_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS_RELEASE ${_FLAGS_RELEASE})

# override install locations when building python extensions
if(DEFINED SKBUILD)
    cmake_path(RELATIVE_PATH SKBUILD_HEADERS_DIR BASE_DIRECTORY ${SKBUILD_PLATLIB_DIR} OUTPUT_VARIABLE CMAKE_INSTALL_INCLUDEDIR)
    cmake_path(RELATIVE_PATH SKBUILD_PLATLIB_DIR BASE_DIRECTORY ${SKBUILD_PLATLIB_DIR} OUTPUT_VARIABLE CMAKE_INSTALL_LIBDIR)
    message(STATUS "Setting include dir to: ${CMAKE_INSTALL_INCLUDEDIR}")
    message(STATUS "Setting lib dir to: ${CMAKE_INSTALL_LIBDIR}")
    set(LI_INSTALL_INCLUDEDIR ${SKBUILD_HEADERS_DIR})
    set(LI_INSTALL_LIBDIR ${SKBUILD_PLATLIB_DIR})
else()
    set(LI_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
    set(LI_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()

# load python and pybind11
include(Python)
include(pybind11)

# load project dependencies
include(rk)
include(cereal)
include(delabella)
include(CFITSIO)
include(photospline)
include(googletest)
include(NamedType)

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
add_library(LeptonInjector SHARED)
set_property(TARGET LeptonInjector PROPERTY POSITION_INDEPENDENT_CODE ON)

if(${MACOSX})
if(CMAKE_VERSION VERSION_LESS 3.13)
    target_link_libraries(LeptonInjector PUBLIC "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")
else()
    target_link_options(LeptonInjector PUBLIC "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")
endif()
target_link_libraries(LeptonInjector
    PRIVATE
    PUBLIC
        delabella_shared
        photospline
        LI_utilities
        LI_serialization
        LI_math
        LI_dataclasses
        LI_geometry
        LI_detector
        LI_interactions
        LI_distributions
        LI_injection
)
else()
target_link_libraries(LeptonInjector
    PRIVATE
    PUBLIC
        photospline
        delabella_shared
        LI_utilities
        LI_serialization
        LI_math
        LI_dataclasses
        LI_geometry
        LI_detector
        LI_interactions
        LI_distributions
        LI_injection
)
endif()

# define the install path normally or for python package
if(DEFINED SKBUILD)
    set_target_properties(LeptonInjector PROPERTIES
            BUILD_WITH_INSTALL_RPATH FALSE
            LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN")
    if(${CIBUILDWHEEL})
        message(STATUS "Setting LeptonInjector install lib dir to: ${CI_INSTALL_PREFIX}/lib")
        message(STATUS "Setting LeptonInjector install include dir to: ${CI_INSTALL_PREFIX}/include")
        install(TARGETS LeptonInjector
            delabella_shared
            LI_utilities
            LI_serialization
            LI_math
            LI_dataclasses
            LI_geometry
            LI_detector
            LI_interactions
            LI_distributions
            LI_injection
            EXPORT ${PROJECT_NAME}Config
            LIBRARY DESTINATION "${CI_INSTALL_PREFIX}/lib"
            PUBLIC_HEADER DESTINATION "${CI_INSTALL_PREFIX}/include")
    else()
        install(TARGETS LeptonInjector
            delabella_shared
            LI_utilities
            LI_serialization
            LI_math
            LI_dataclasses
            LI_geometry
            LI_detector
            LI_interactions
            LI_distributions
            LI_injection
            EXPORT ${PROJECT_NAME}Config
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector.libs
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    endif()
else()
    install(TARGETS LeptonInjector
        delabella_shared
        LI_utilities
        LI_serialization
        LI_math
        LI_dataclasses
        LI_geometry
        LI_detector
        LI_interactions
        LI_distributions
        LI_injection
        EXPORT ${PROJECT_NAME}Config
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# optionally package runtime dependencies
if((DEFINED SKBUILD) AND (PACKAGE_SHARED_DEPS))
    install(CODE "set(LeptonInjector_LIB_FILE \"${PROJECT_BINARY_DIR}/${CMAKE_SHARED_MODULE_PREFIX}LeptonInjector${CMAKE_SHARED_MODULE_SUFFIX}\")")
    install(CODE "set(PYTHON_DEP_LIB_DESTINATION \"${SKBUILD_PLATLIB_DIR}/leptoninjector.libs/\")")
    install(CODE [[
        file(GET_RUNTIME_DEPENDENCIES
            LIBRARIES ${LeptonInjector_LIB_FILE}
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
if(DEFINED SKBUILD)
    install(TARGETS utilities
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS math
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS dataclasses
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS geometry
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS detector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS interactions
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS distributions
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS injection
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/resources
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
endif()

# Export targets for use in downstream CMake projects
# -----------------------------------------------------------------------------
include(CMakePackageConfigHelpers)
# Make importable from build directory
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT ${PROJECT_NAME}Config FILE ${PROJECT_NAME}Config.cmake)

# Make importable from install location
set(_config_dir share/${PROJECT_NAME}/cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${_config_dir}
)
install(EXPORT ${PROJECT_NAME}Config
  DESTINATION ${_config_dir}
)

MESSAGE("")
MESSAGE("Done!")
MESSAGE("Run `make` then `make install`")
