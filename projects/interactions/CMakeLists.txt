
#core files
LIST (APPEND interactions_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/Decay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/InteractionCollection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DipoleFromTable.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DarkNewsCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DarkNewsDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/NeutrissimoDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/ElasticScattering.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DummyCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyDarkNewsCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyDarkNewsDecay.cxx
)
add_library(SIREN_interactions OBJECT ${interactions_SOURCES})
set_property(TARGET SIREN_interactions PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_interactions PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/interactions/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/interactions/>
    ${PYTHON_INCLUDE_DIRS}
)

target_link_libraries(SIREN_interactions
    INTERFACE
        siren_compile_options
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
        pybind11::embed
    PUBLIC
        photospline
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_detector
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
target_link_libraries(interactions PRIVATE SIREN photospline rk_static pybind11::embed)
#pybind11_add_module(pyDarkNewsSerializer ${PROJECT_SOURCE_DIR}/projects/interactions/private/pybindings/pyDarkNewsSerializer.cxx)
#target_link_libraries(pyDarkNewsSerializer PRIVATE SIREN photospline rk_static pybind11::embed)
if(DEFINED SKBUILD)
set_target_properties(interactions PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
