
#core files
LIST (APPEND detector_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/detector/private/AtomicMasses.cxx

    ${CMAKE_SOURCE_DIR}/projects/detector/private/Axis1D.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/CartesianAxis1D.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/RadialAxis1D.cxx

    ${CMAKE_SOURCE_DIR}/projects/detector/private/Distribution1D.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/ConstantDistribution1D.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/PolynomialDistribution1D.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/ExponentialDistribution1D.cxx

    ${CMAKE_SOURCE_DIR}/projects/detector/private/DensityDistribution.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/EarthModel.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/MaterialModel.cxx
    ${CMAKE_SOURCE_DIR}/projects/detector/private/Path.cxx
)
add_library(LI_detector OBJECT ${detector_SOURCES})
set_property(TARGET LI_detector PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_detector PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/detector/public/
)
add_dependencies(LI_detector rk)
add_dependencies(LI_detector LI_serialization LI_utilities LI_math LI_dataclasses LI_geometry)

target_link_libraries(LI_detector photospline)
target_link_libraries(LI_detector ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_detector PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/dataclasses/public/)
target_include_directories(LI_detector PUBLIC ${CMAKE_SOURCE_DIR}/projects/geometry/public/)

install(TARGETS LI_detector
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/detector/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_Axis ${CMAKE_SOURCE_DIR}/projects/detector/private/test/Axis_TEST.cxx)
package_add_test(UnitTest_DensityDistribution ${CMAKE_SOURCE_DIR}/projects/detector/private/test/DensityDistribution_TEST.cxx)
package_add_test(UnitTest_Distribution1D ${CMAKE_SOURCE_DIR}/projects/detector/private/test/Distribution1D_TEST.cxx)
package_add_test(UnitTest_EarthModel ${CMAKE_SOURCE_DIR}/projects/detector/private/test/EarthModel_TEST.cxx)
package_add_test(UnitTest_MaterialModel ${CMAKE_SOURCE_DIR}/projects/detector/private/test/MaterialModel_TEST.cxx)
package_add_test(UnitTest_Path ${CMAKE_SOURCE_DIR}/projects/detector/private/test/Path_TEST.cxx)

pybind11_add_module(Detector ${CMAKE_SOURCE_DIR}/projects/detector/private/pybindings/detector.cxx)
target_link_libraries(Detector PRIVATE LI_detector)
