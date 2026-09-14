# Source-wheel builds already run the backend; never start a nested wheel build.
if(NOT SIREN_PYTHON_PACKAGE OR DEFINED SKBUILD)
    return()
endif()

set(PACKAGE_STAGING_DIR "${CMAKE_BINARY_DIR}/python_staging")
set(WHEELS_DIR "${CMAKE_BINARY_DIR}/dist_wheels")
file(GLOB_RECURSE PYTHON_PACKAGE_FILES LIST_DIRECTORIES false CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/python/*" "${PROJECT_SOURCE_DIR}/resources/*")
# A removed input must invalidate the wheel as well as an added/edited input.
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/python_package_files.txt"
    CONTENT "${PYTHON_PACKAGE_FILES}\n")

add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/.build_wheel"
    COMMAND ${CMAKE_COMMAND} -E rm -rf "${PACKAGE_STAGING_DIR}" "${WHEELS_DIR}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${PACKAGE_STAGING_DIR}" "${WHEELS_DIR}"
    COMMAND ${CMAKE_COMMAND} --install "${CMAKE_BINARY_DIR}"
        --config $<CONFIG> --prefix "${PACKAGE_STAGING_DIR}" --component PythonWheel
    COMMAND ${CMAKE_COMMAND} -E copy
        "${PROJECT_SOURCE_DIR}/pyproject.toml" "${PROJECT_SOURCE_DIR}/README.md"
        "${PROJECT_SOURCE_DIR}/LICENSE" "${PACKAGE_STAGING_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy "${PROJECT_SOURCE_DIR}/package/CMakeLists.txt"
        "${PACKAGE_STAGING_DIR}/CMakeLists.txt"
    COMMAND ${Python_EXECUTABLE} -m pip wheel --no-deps
        --wheel-dir "${WHEELS_DIR}" "${PACKAGE_STAGING_DIR}"
    COMMAND ${CMAKE_COMMAND} -E touch "${CMAKE_CURRENT_BINARY_DIR}/.build_wheel"
    DEPENDS ${SIREN_WHEEL_LIBRARIES} ${SIREN_PYTHON_MODULES} ${PYTHON_PACKAGE_FILES}
        "${CMAKE_CURRENT_BINARY_DIR}/python_package_files.txt"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
        "${PROJECT_SOURCE_DIR}/package/CMakeLists.txt"
        "${PROJECT_SOURCE_DIR}/pyproject.toml" "${PROJECT_SOURCE_DIR}/README.md"
        "${PROJECT_SOURCE_DIR}/LICENSE"
    COMMENT "Building a native wheel from the CMake-installed package"
    VERBATIM)
add_custom_target(python_package ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/.build_wheel")

install(CODE "
    file(GLOB WHEELS \"${WHEELS_DIR}/*.whl\")
    list(LENGTH WHEELS WHEEL_COUNT)
    if(NOT WHEEL_COUNT EQUAL 1)
        message(FATAL_ERROR \"Build the python_package target before installing SIREN\")
    endif()
    execute_process(
        COMMAND \"${Python_EXECUTABLE}\" -m pip install --no-deps --force-reinstall
            \$\{WHEELS\} --prefix=\$\{CMAKE_INSTALL_PREFIX\}
        RESULT_VARIABLE WHEEL_INSTALL_RESULT
        COMMAND_ECHO STDOUT)
    if(NOT WHEEL_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR \"SIREN wheel installation failed\")
    endif()
")
