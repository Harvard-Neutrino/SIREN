
#core files
LIST (APPEND crosssections_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/Decay.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/CrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/CrossSectionCollection.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/DipoleFromTable.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/NeutrissimoDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/DISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/HNLFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/ElasticScattering.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/DummyCrossSection.cxx
)
add_library(LI_crosssections OBJECT ${crosssections_SOURCES})
set_property(TARGET LI_crosssections PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_crosssections PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/crosssections/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/LeptonInjector/crosssections/>
)

target_link_libraries(LI_crosssections
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
    PUBLIC
        photospline
        LI_serialization
        LI_utilities
        LI_math
        LI_dataclasses
        LI_detector
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/crosssections/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_DISFromSpline ${PROJECT_SOURCE_DIR}/projects/crosssections/private/test/DISFromSpline_TEST.cxx)
package_add_test(UnitTest_DipoleFromTable ${PROJECT_SOURCE_DIR}/projects/crosssections/private/test/DipoleFromTable_TEST.cxx)
#package_add_test(UnitTest_ElasticScattering ${PROJECT_SOURCE_DIR}/projects/crosssections/private/test/ElasticScattering_TEST.cxx)

pybind11_add_module(crosssections ${PROJECT_SOURCE_DIR}/projects/crosssections/private/pybindings/crosssections.cxx)
target_link_libraries(crosssections PRIVATE LeptonInjector photospline rk_static)
pybind11_add_module(pyDarkNewsCrossSectionSerializer ${PROJECT_SOURCE_DIR}/projects/crosssections/private/pybindings/pyDarkNewsCrossSectionSerializer.cxx)
target_link_libraries(pyDarkNewsCrossSectionSerializer PRIVATE LeptonInjector photospline rk_static pybind11::embed)
if(DEFINED SKBUILD)
set_target_properties(crosssections PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
