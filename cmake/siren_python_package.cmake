# Source-wheel builds already run the backend; never start a nested wheel build.
if(NOT SIREN_PYTHON_PACKAGE OR DEFINED SKBUILD)
    return()
endif()

set(WHEELS_DIR "${CMAKE_BINARY_DIR}/dist_wheels")
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
file(GLOB_RECURSE PYTHON_PACKAGE_FILES LIST_DIRECTORIES false CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/python/*" "${PROJECT_SOURCE_DIR}/resources/*")
# Match the PythonWheel directory-install exclusions. Local imports must not
# invalidate a wheel whose payload excludes Python bytecode.
list(FILTER PYTHON_PACKAGE_FILES EXCLUDE REGEX "/__pycache__(/|$)|\\.pyc$")
# A removed input must invalidate the wheel as well as an added/edited input.
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/python_package_files.txt"
    CONTENT "${PYTHON_PACKAGE_FILES}\n")
set(wheel_input_files ${PYTHON_PACKAGE_FILES}
    "${PROJECT_SOURCE_DIR}/package/CMakeLists.txt"
    "${PROJECT_SOURCE_DIR}/pyproject.toml" "${PROJECT_SOURCE_DIR}/README.md"
    "${PROJECT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_LIST_DIR}/build_wheel.py"
    "${CMAKE_CURRENT_LIST_DIR}/wheel_rpath.py")
foreach(target IN LISTS SIREN_WHEEL_LIBRARIES SIREN_PYTHON_MODULES)
    list(APPEND wheel_input_files "$<TARGET_FILE:${target}>")
endforeach()
list(JOIN wheel_input_files "\n" wheel_input_manifest)
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/wheel_inputs-$<CONFIG>.txt"
    CONTENT "${wheel_input_manifest}\n")

# Check contents on every invocation, including mtime-preserving restores.
# The driver only stages/packages when its recorded inputs actually change.
add_custom_target(python_package ALL
    COMMAND ${CMAKE_COMMAND} -E env ${SIREN_WHEEL_ENV}
        ${Python_EXECUTABLE} "${CMAKE_CURRENT_LIST_DIR}/build_wheel.py"
        --source "${PROJECT_SOURCE_DIR}" --build "${CMAKE_BINARY_DIR}"
        --library-dir "${SIREN_WHEEL_LIBRARY_DIR}"
        --inputs "${CMAKE_CURRENT_BINARY_DIR}/wheel_inputs-$<CONFIG>.txt"
        --config $<CONFIG> --cmake "${CMAKE_COMMAND}" ${SIREN_WHEEL_PIP_OPTIONS}
    DEPENDS ${SIREN_WHEEL_LIBRARIES} ${SIREN_PYTHON_MODULES} ${PYTHON_PACKAGE_FILES}
        "${CMAKE_CURRENT_BINARY_DIR}/python_package_files.txt"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/wheel_inputs-$<CONFIG>.txt"
        "${CMAKE_CURRENT_LIST_DIR}/build_wheel.py"
        "${PROJECT_SOURCE_DIR}/package/CMakeLists.txt"
        "${PROJECT_SOURCE_DIR}/pyproject.toml" "${PROJECT_SOURCE_DIR}/README.md"
        "${PROJECT_SOURCE_DIR}/LICENSE"
    COMMENT "Building a native wheel from the CMake-installed package"
    VERBATIM)

# Select replacement/isolation using the actual destination, including DESTDIR.
install(CODE "
    file(GLOB WHEELS \"${WHEELS_DIR}/*.whl\")
    list(LENGTH WHEELS WHEEL_COUNT)
    if(NOT WHEEL_COUNT EQUAL 1)
        message(FATAL_ERROR \"Build the python_package target before installing SIREN, or configure -DSIREN_PYTHON_PACKAGE=OFF for a native-only install\")
    endif()
    execute_process(
        COMMAND \"${Python_EXECUTABLE}\" \"${CMAKE_CURRENT_LIST_DIR}/install_wheel.py\"
            --prefix \"\$\{CMAKE_INSTALL_PREFIX\}\" \$\{WHEELS\}
        RESULT_VARIABLE WHEEL_INSTALL_RESULT
        COMMAND_ECHO STDOUT)
    if(NOT WHEEL_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR \"SIREN wheel installation failed\")
    endif()
" COMPONENT PythonPackage)
