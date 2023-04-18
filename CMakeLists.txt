cmake_minimum_required(VERSION 3.3.2 FATAL_ERROR)
cmake_policy(VERSION 3.3.2)

# point to appropriate cmake file locations
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_CURRENT_LIST_DIR}/cmake/Packages)

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
    set(LI_INSTALL_INCLUDEDIR ${SKBUILD_HEADERS_DIR})
    set(LI_INSTALL_LIBDIR ${SKBUILD_PLATLIB_DIR})
else()
    set(LI_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
    set(LI_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif()

# load python and pybind11
include(pybind11)
if(CMAKE_VERSION VERSION_LESS 3.12.0)
    if(Python_ROOT_DIR)
        MESSAGE(WARNING "Python_ROOT_DIR is set but will be ignored by this version of CMake; set PYTHON_EXECUTABLE instead")
	endif()
	find_package(PythonInterp)
	find_package(PythonLibs)
	set(PYTHON_FOUND PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
else()
    if(PYTHON_EXECUTABLE)
        MESSAGE(WARNING "PYTHON_EXECUTABLE is set but will be ignored by this version of CMake; set Python_ROOT_DIR instead")
	endif()
    find_package(Python COMPONENTS Interpreter Development)
    set(PYTHON_FOUND Python_Interpreter_FOUND AND Python_Development_FOUND)
    set(PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
	set(PYTHON_INCLUDE_DIRS "${Python_INCLUDE_DIRS}")
	set(PYTHON_LIBRARIES "${Python_LIBRARIES}")
	set(PYTHON_VERSION_MAJOR "${Python_VERSION_MAJOR}")
	set(PYTHON_VERSION_MINOR "${Python_VERSION_MINOR}")
endif()

# load project dependencies
include(rk)
include(cereal)
include(delabella)
include(photospline)
include(googletest)

# load macros for googletest
include(testing)

# add individual source folders
add_subdirectory(projects/utilities)
add_subdirectory(projects/serialization)
add_subdirectory(projects/math)
add_subdirectory(projects/dataclasses)
add_subdirectory(projects/geometry)
add_subdirectory(projects/detector)
add_subdirectory(projects/crosssections)
add_subdirectory(projects/distributions)
add_subdirectory(projects/injection)

# define the target library
add_library(LeptonInjector SHARED)
set_property(TARGET LeptonInjector PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(LeptonInjector
    LI_utilities
    LI_serialization
    LI_math
    LI_dataclasses
    LI_geometry
    LI_detector
    LI_crosssections
    LI_distributions
    LI_injection
    -Wl,--whole-archive
    ${PROJECT_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_STATIC_LIBRARY_SUFFIX}
    delabella
    -Wl,--no-whole-archive
)

# define the install path normally or for python package
if(DEFINED SKBUILD)
    set_target_properties(LeptonInjector PROPERTIES
            BUILD_WITH_INSTALL_RPATH FALSE
            LINK_FLAGS "-Wl,-rpath,$ORIGIN")
    install(TARGETS LeptonInjector
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector.libs
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
else()
    install(TARGETS LeptonInjector
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
    install(TARGETS crosssections
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS distributions
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(TARGETS injection
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
    install(FILES ${PROJECT_SOURCE_DIR}/python/__init__.py DESTINATION ${CMAKE_INSTALL_LIBDIR}/leptoninjector)
endif()

MESSAGE("")
MESSAGE("Done!")
MESSAGE("Run `make` then `make install`")
