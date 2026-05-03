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

message(STATUS "Added crundec3")
