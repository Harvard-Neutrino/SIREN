
#core files
LIST (APPEND distributions_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/ColumnDepthPositionDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/CylinderVolumePositionDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/DecayRangePositionDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/Distributions.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/OrientedCylinderPositionDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/PrimaryDirectionDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/PrimaryEnergyDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/PrimaryInjector.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/PrimaryNeutrinoHelicityDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/TargetMomentumDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/distributions/private/VertexPositionDistribution.cxx
)
add_library(LI_distributions SHARED ${distributions_SOURCES})
target_include_directories(LI_distributions PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/distributions/public/
)
add_dependencies(LI_distributions rk)
add_dependencies(LI_distributions LI_serialization LI_utilities LI_math LI_geometry LI_detector LI_crosssections)

target_link_libraries(LI_distributions ${HDF5_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(LI_distributions photospline)
target_link_libraries(LI_distributions ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_distributions PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_distributions PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_distributions PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_distributions PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_distributions PUBLIC ${CMAKE_SOURCE_DIR}/projects/geometry/public/)
target_include_directories(LI_distributions PUBLIC ${CMAKE_SOURCE_DIR}/projects/detector/public/)
target_include_directories(LI_distributions PUBLIC ${CMAKE_SOURCE_DIR}/projects/crosssections/public/)

install(TARGETS LI_distributions
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/distributions/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

