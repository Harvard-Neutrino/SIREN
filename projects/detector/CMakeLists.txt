
#core files
LIST (APPEND detector_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/detector/private/AtomicMasses.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/DensityDist.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/EarthModelCalculator.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/EarthModel.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/MaterialModel.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/Path.cxx
)
add_library(LI_detector SHARED ${detector_SOURCES})
target_include_directories(LI_detector PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/detector/public/
)
add_dependencies(LI_detector rk)
add_dependencies(LI_detector LI_serialization LI_utilities LI_geometry)

target_link_libraries(LI_detector ${HDF5_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(LI_detector photospline)
target_link_libraries(LI_detector ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_detector PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/geometry/public/)

install(TARGETS LI_detector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/detector/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

