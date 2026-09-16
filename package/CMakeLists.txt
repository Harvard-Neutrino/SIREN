cmake_minimum_required(VERSION 3.20)
project(siren_prebuilt LANGUAGES NONE)

# Binaries were installed by the outer build, including their install RPATHs.
# Use the same backend/metadata as a source build to retain native wheel tags.
install(DIRECTORY siren DESTINATION . COMPONENT PythonWheel)
if(NOT DEFINED SIREN_WHEEL_LIBRARY_DIR)
    message(FATAL_ERROR "The outer wheel build must provide SIREN_WHEEL_LIBRARY_DIR")
endif()
install(DIRECTORY "${SIREN_WHEEL_LIBRARY_DIR}" DESTINATION . COMPONENT PythonWheel OPTIONAL)
