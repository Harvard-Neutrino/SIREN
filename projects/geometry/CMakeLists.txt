
#core files
LIST (APPEND geometry_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Geometry.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Vector3D.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Matrix3D.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Polynomial.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Placement.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/Quaternion.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/EulerAngles.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/GeometryMesh.cxx
    ${CMAKE_SOURCE_DIR}/projects/geometry/private/MeshBuilder.cxx
)
add_library(LI_geometry SHARED ${geometry_SOURCES})
target_include_directories(LI_geometry PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/geometry/public/
)
add_dependencies(LI_geometry rk)
add_dependencies(LI_geometry LI_serialization)

target_link_libraries(LI_geometry ${HDF5_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(LI_geometry photospline)
target_link_libraries(LI_geometry ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_geometry PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_geometry PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)

install(TARGETS LI_geometry
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/geometry/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

