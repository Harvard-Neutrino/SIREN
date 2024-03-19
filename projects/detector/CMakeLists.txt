
#core files
LIST (APPEND detector_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/detector/private/AtomicMasses.cxx

    ${PROJECT_SOURCE_DIR}/projects/detector/private/Axis1D.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/CartesianAxis1D.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/RadialAxis1D.cxx

    ${PROJECT_SOURCE_DIR}/projects/detector/private/Distribution1D.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/ConstantDistribution1D.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/PolynomialDistribution1D.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/ExponentialDistribution1D.cxx

    ${PROJECT_SOURCE_DIR}/projects/detector/private/DensityDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/DetectorModel.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/MaterialModel.cxx
    ${PROJECT_SOURCE_DIR}/projects/detector/private/Path.cxx
)
add_library(SIREN_detector OBJECT ${detector_SOURCES})
set_property(TARGET SIREN_detector PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_detector PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/detector/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/detector/>
)

target_link_libraries(SIREN_detector
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
    PUBLIC
        photospline
        SIREN_serialization
        SIREN_utilities
        SIREN_math
        SIREN_dataclasses
        SIREN_geometry
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/detector/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_Axis ${PROJECT_SOURCE_DIR}/projects/detector/private/test/Axis_TEST.cxx)
package_add_test(UnitTest_DensityDistribution ${PROJECT_SOURCE_DIR}/projects/detector/private/test/DensityDistribution_TEST.cxx)
package_add_test(UnitTest_Distribution1D ${PROJECT_SOURCE_DIR}/projects/detector/private/test/Distribution1D_TEST.cxx)
package_add_test(UnitTest_DetectorModel ${PROJECT_SOURCE_DIR}/projects/detector/private/test/DetectorModel_TEST.cxx)
package_add_test(UnitTest_MaterialModel ${PROJECT_SOURCE_DIR}/projects/detector/private/test/MaterialModel_TEST.cxx)
package_add_test(UnitTest_Path ${PROJECT_SOURCE_DIR}/projects/detector/private/test/Path_TEST.cxx)

pybind11_add_module(detector ${PROJECT_SOURCE_DIR}/projects/detector/private/pybindings/detector.cxx)
target_link_libraries(detector PRIVATE SIREN)
if(DEFINED SKBUILD)
set_target_properties(detector PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
