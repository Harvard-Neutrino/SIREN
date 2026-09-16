# Source-wheel builds already run the backend; never start a nested wheel build.
if(NOT SIREN_PYTHON_PACKAGE OR DEFINED SKBUILD)
    return()
endif()

option(SIREN_WHEEL_BUILD_ISOLATION
    "Install wheel backend requirements in an isolated environment" ON)
set(SIREN_WHEEL_PIP_OPTIONS)
if(NOT SIREN_WHEEL_BUILD_ISOLATION)
    list(APPEND SIREN_WHEEL_PIP_OPTIONS --no-build-isolation)
endif()
set(SIREN_WHEEL_ENV)
if(APPLE)
    if(CMAKE_OSX_DEPLOYMENT_TARGET)
        list(APPEND SIREN_WHEEL_ENV
            "MACOSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}")
    endif()
    if(CMAKE_OSX_ARCHITECTURES)
        set(wheel_archflags "")
        foreach(architecture IN LISTS CMAKE_OSX_ARCHITECTURES)
            string(APPEND wheel_archflags " -arch ${architecture}")
        endforeach()
        string(STRIP "${wheel_archflags}" wheel_archflags)
        list(APPEND SIREN_WHEEL_ENV "ARCHFLAGS=${wheel_archflags}")
    endif()
endif()

# The driver installs the PythonWheel component into a fresh staging tree on
# every build, so CMake's install rules are the only definition of the wheel
# contents, and rebuilds the wheel only when that staged tree or the packaging
# inputs change. Installing the wheel is a separate, explicit pip step; no
# CMake install component invokes pip (see README, "Installing the wheel").
add_custom_target(python_package ALL
    COMMAND ${CMAKE_COMMAND} -E env ${SIREN_WHEEL_ENV}
        ${Python_EXECUTABLE} "${CMAKE_CURRENT_LIST_DIR}/build_wheel.py"
        --source "${PROJECT_SOURCE_DIR}" --build "${CMAKE_BINARY_DIR}"
        --library-dir "${SIREN_WHEEL_LIBRARY_DIR}"
        --config $<CONFIG> --cmake "${CMAKE_COMMAND}" ${SIREN_WHEEL_PIP_OPTIONS}
    DEPENDS ${SIREN_WHEEL_LIBRARIES} ${SIREN_PYTHON_MODULES}
    COMMENT "Building a native wheel from the CMake-installed package"
    VERBATIM)
