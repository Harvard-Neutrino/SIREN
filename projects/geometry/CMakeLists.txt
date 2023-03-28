
#core files
LIST (APPEND geometry_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Geometry.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Box.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Cylinder.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Sphere.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/ExtrPoly.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/GeometryMesh.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/MeshBuilder.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Placement.cxx
)
add_library(LI_geometry OBJECT ${geometry_SOURCES})
set_property(TARGET LI_geometry PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_geometry PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/geometry/public/
)
add_dependencies(LI_geometry rk)
add_dependencies(LI_geometry LI_serialization LI_math LI_geometry)


target_link_libraries(LI_geometry photospline)
target_link_libraries(LI_geometry ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_geometry PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_geometry PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_geometry PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_geometry PUBLIC ${CMAKE_SOURCE_DIR}/projects/geometry/public/)

install(TARGETS LI_geometry
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/geometry/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_Geometry ${CMAKE_SOURCE_DIR}/projects/geometry/private/test/Geometry_TEST.cxx)
package_add_test(UnitTest_ExtrPoly ${CMAKE_SOURCE_DIR}/projects/geometry/private/test/ExtrPoly_TEST.cxx)

pybind11_add_module(geometry ${CMAKE_SOURCE_DIR}/projects/geometry/private/pybindings/geometry.cxx)
target_link_libraries(geometry PRIVATE LI_geometry)
