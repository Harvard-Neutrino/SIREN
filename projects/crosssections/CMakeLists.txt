
#core files
LIST (APPEND crosssections_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/CrossSection.cxx
)
add_library(LI_crosssections SHARED ${crosssections_SOURCES})
target_include_directories(LI_crosssections PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/crosssections/public/
)
add_dependencies(LI_crosssections rk)
add_dependencies(LI_crosssections LI_serialization LI_utilities LI_geometry LI_detector)

target_link_libraries(LI_crosssections ${HDF5_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(LI_crosssections photospline)
target_link_libraries(LI_crosssections ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_crosssections PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/geometry/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/detector/public/)

install(TARGETS LI_crosssections
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/crosssections/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

