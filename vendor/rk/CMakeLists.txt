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

add_library(rk STATIC ${rk_SOURCES})
add_library(rk_shared SHARED ${rk_SOURCES})

target_include_directories(rk PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(rk_shared PUBLIC ${PROJECT_SOURCE_DIR})

install(TARGETS rk_shared
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/rk/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.hh"
    PATTERN "*.icc"
)
message("Added rk!")
