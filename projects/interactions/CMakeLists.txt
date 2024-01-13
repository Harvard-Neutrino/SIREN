
#core files
LIST (APPEND interactions_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/Decay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/InteractionCollection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DipoleFromTable.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DarkNewsCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DarkNewsDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/NeutrissimoDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/ElasticScattering.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DummyCrossSection.cxx
)
add_library(LI_interactions OBJECT ${interactions_SOURCES})
set_property(TARGET LI_interactions PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_interactions PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/interactions/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/LeptonInjector/interactions/>
)

target_link_libraries(LI_interactions
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
        pybind11::embed
        pybind11::python_link_helper
    PUBLIC
        photospline
        LI_serialization
        LI_utilities
        LI_math
        LI_dataclasses
        LI_detector
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/interactions/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_DISFromSpline ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/DISFromSpline_TEST.cxx)
package_add_test(UnitTest_DipoleFromTable ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/DipoleFromTable_TEST.cxx)
#package_add_test(UnitTest_ElasticScattering ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/ElasticScattering_TEST.cxx)

pybind11_add_module(interactions ${PROJECT_SOURCE_DIR}/projects/interactions/private/pybindings/interactions.cxx)
target_link_libraries(interactions PRIVATE LeptonInjector photospline rk_static pybind11::embed)
pybind11_add_module(pyDarkNewsSerializer ${PROJECT_SOURCE_DIR}/projects/interactions/private/pybindings/pyDarkNewsSerializer.cxx)
target_link_libraries(pyDarkNewsSerializer PRIVATE LeptonInjector photospline rk_static pybind11::embed)
if(DEFINED SKBUILD)
set_target_properties(interactions PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
