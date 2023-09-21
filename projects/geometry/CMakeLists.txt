
#core files
LIST (APPEND geometry_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Geometry.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Box.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Cylinder.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Sphere.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/ExtrPoly.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/GeometryMesh.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/MeshBuilder.cxx
    ${PROJECT_SOURCE_DIR}/projects/geometry/private/Placement.cxx
)
add_library(LI_geometry OBJECT ${geometry_SOURCES})
set_property(TARGET LI_geometry PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_geometry PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/geometry/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/LeptonInjector/geometry/>
)

target_link_libraries(LI_geometry
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
    PUBLIC
        photospline
        LI_serialization
        LI_math
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/geometry/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_Geometry ${PROJECT_SOURCE_DIR}/projects/geometry/private/test/Geometry_TEST.cxx)
package_add_test(UnitTest_ExtrPoly ${PROJECT_SOURCE_DIR}/projects/geometry/private/test/ExtrPoly_TEST.cxx)

pybind11_add_module(geometry ${PROJECT_SOURCE_DIR}/projects/geometry/private/pybindings/geometry.cxx)
target_link_libraries(geometry PRIVATE LeptonInjector)
if(DEFINED SKBUILD)
set_target_properties(geometry PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
