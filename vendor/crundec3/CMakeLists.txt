cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
cmake_policy(VERSION 3.5)

project(crundec VERSION 1.7 LANGUAGES CXX)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_C_STANDARD 99)

LIST(APPEND crundec_SOURCES
    ${PROJECT_SOURCE_DIR}/CRunDec3.1/CRunDec.cpp
)

LIST(APPEND crundec_HEADERS
    ${PROJECT_SOURCE_DIR}/CRunDec3.1/CRunDec.h
)

add_library(crundec_static STATIC ${crundec_SOURCES})
set_target_properties(crundec_static PROPERTIES EXPORT_NAME crundec)

target_include_directories(crundec_static PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS crundec_static
    EXPORT ${PROJECT_NAME}Config
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/CRunDec3.1/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.m"
)

### Export targets for use with the install tree
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
export(EXPORT ${PROJECT_NAME}Config FILE ${PROJECT_NAME}Config.cmake)

set(_config_dir share/${PROJECT_NAME}/cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${_config_dir}
)
install(EXPORT ${PROJECT_NAME}Config
  DESTINATION ${_config_dir}
)

message(STATUS "Added crundec3")
