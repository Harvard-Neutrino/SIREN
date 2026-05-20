
#core files
LIST (APPEND geometry_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Geometry.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Box.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Cone.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Trd.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Cylinder.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Sphere.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/ExtrPoly.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Polycone.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Polyhedra.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/GeometryMesh.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/MeshBuilder.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Placement.cxx
)
add_library(SIREN_geometry OBJECT ${geometry_SOURCES})
set_property(TARGET SIREN_geometry PROPERTY POSITION_INDEPENDENT_CODE ON)

# Woop's watertight ray-triangle intersection requires that each multiply rounds
# independently (no FMA contraction). noinline prevents cross-call-boundary FMA;
# this flag prevents FMA within the function body.
if(NOT MSVC)
    set_source_files_properties(
        ${PROJECT_SOURCE_DIR}/projects/geometry/private/GeometryMesh.cxx
        PROPERTIES COMPILE_OPTIONS "-ffp-contract=off"
    )
endif()
target_include_directories(SIREN_geometry PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/geometry/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/geometry/>
)

target_link_libraries(SIREN_geometry
    PUBLIC
        photospline
        SIREN_serialization
        SIREN_math
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
)
apply_siren_compile_options(SIREN_geometry)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/geometry/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

pybind11_add_module(geometry ${PROJECT_SOURCE_DIR}/projects/geometry/private/pybindings/geometry.cxx)
target_link_libraries(geometry PRIVATE SIREN)
if(DEFINED SKBUILD)
set_target_properties(geometry PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
