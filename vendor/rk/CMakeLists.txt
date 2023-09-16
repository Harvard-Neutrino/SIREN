cmake_minimum_required(VERSION 3.3.2 FATAL_ERROR)
cmake_policy(VERSION 3.3.2)

project(rk VERSION 1.7 LANGUAGES CXX)

include(GNUInstallDirs)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_C_STANDARD 99)

LIST(APPEND rk_SOURCES
    ${PROJECT_SOURCE_DIR}/rk/rk.cc
    ${PROJECT_SOURCE_DIR}/rk/geom3.cc
)

LIST(APPEND rk_HEADERS
    ${PROJECT_SOURCE_DIR}/rk/rk.hh
    ${PROJECT_SOURCE_DIR}/rk/geom3.hh
    ${PROJECT_SOURCE_DIR}/rk/geom3_Matrix3x3.icc
    ${PROJECT_SOURCE_DIR}/rk/geom3_Point3.icc
    ${PROJECT_SOURCE_DIR}/rk/geom3_Rotation3.icc
    ${PROJECT_SOURCE_DIR}/rk/geom3_UnitVector3.icc
    ${PROJECT_SOURCE_DIR}/rk/geom3_Vector3.icc
    ${PROJECT_SOURCE_DIR}/rk/rk_Boost.icc
    ${PROJECT_SOURCE_DIR}/rk/rk_P4.icc
    ${PROJECT_SOURCE_DIR}/rk/rk_Point4.icc
    ${PROJECT_SOURCE_DIR}/rk/LT.hh
    ${PROJECT_SOURCE_DIR}/rk/rk_LT.icc
    ${PROJECT_SOURCE_DIR}/rk/rkIO.hh
)

add_library(rk_static STATIC ${rk_SOURCES})
add_library(rk_shared SHARED ${rk_SOURCES})
set_target_properties(rk_static PROPERTIES EXPORT_NAME rk)
set_target_properties(rk_shared PROPERTIES EXPORT_NAME rk)

target_include_directories(rk_static PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(rk_shared PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS rk_static
    EXPORT ${PROJECT_NAME}Config
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS rk_shared
    EXPORT ${PROJECT_NAME}Config
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/rk/"
    #EXPORT ${PROJECT_NAME}Config
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.hh"
    PATTERN "*.icc"
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
export(EXPORT ${PROJECT_NAME}Config FILE ${PROJECT_NAME}Config.cmake)

# Make importable from install location
set(_config_dir share/${PROJECT_NAME}/cmake)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${_config_dir}
)
install(EXPORT ${PROJECT_NAME}Config
  DESTINATION ${_config_dir}
)

message("Added rk!")
