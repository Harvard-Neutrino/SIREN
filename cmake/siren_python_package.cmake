if(${SIREN_PYTHON_PACKAGE})

# Stage the python package within the build directory
set(PACKAGE_STAGING_DIR "${CMAKE_BINARY_DIR}/python_staging")

file(MAKE_DIRECTORY ${PACKAGE_STAGING_DIR})

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
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.pyproject
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
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/.python_package_configured
        ${PACKAGE_STAGING_DIR}/pyproject.toml
        ${CMAKE_CURRENT_BINARY_DIR}/configure_python_package.cmake
        ${PROJECT_SOURCE_DIR}/package/configure_python_package.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_resources
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_python
        ${CMAKE_CURRENT_BINARY_DIR}/.stamp_copy_extensions
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/.build_wheel
    COMMENT "Building wheel from the staged package"
)

add_custom_target(
    python_package ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/.build_wheel
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

endif()
