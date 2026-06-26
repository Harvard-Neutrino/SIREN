
#core files
LIST (APPEND interactions_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/Decay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/InteractionCollection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DarkNewsCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DarkNewsDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLDipoleDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLDipoleFromTable.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLDISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLDipoleDISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/HNLDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/ElectroweakDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/ElasticScattering.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DummyCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyDarkNewsCrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/pyDarkNewsDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CharmMesonDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CharmMesonDecay3Body.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DMesonELoss.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/QuarkDISFromSpline.cxx

)
if(TARGET MARLEY)
    LIST (APPEND interactions_SOURCES
        ${PROJECT_SOURCE_DIR}/projects/interactions/private/MarleyCrossSection.cxx
    )
endif()
if(TARGET PYTHIA8)
    LIST (APPEND interactions_SOURCES
        ${PROJECT_SOURCE_DIR}/projects/interactions/private/PythiaDISCrossSection.cxx
    )
endif()
add_library(SIREN_interactions OBJECT ${interactions_SOURCES})
set_property(TARGET SIREN_interactions PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(SIREN_interactions PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/interactions/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/interactions/>
    ${PYTHON_INCLUDE_DIRS}
)
find_package(GSL REQUIRED)
target_link_libraries(SIREN_interactions
    PRIVATE
        pybind11::embed
        GSL::gsl
        GSL::gslcblas
    PUBLIC
        photospline
        $<BUILD_INTERFACE:crundec_static>
        $<BUILD_INTERFACE:rk_static>
        $<$<TARGET_EXISTS:MARLEY>:MARLEY>
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_detector
)
apply_siren_compile_options(SIREN_interactions)

if(TARGET PYTHIA8)
    target_link_libraries(SIREN_interactions PUBLIC PYTHIA8)
    target_compile_definitions(SIREN_interactions PUBLIC SIREN_HAS_PYTHIA8)
endif()

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/interactions/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_DISFromSpline ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/DISFromSpline_TEST.cxx)
package_add_test(UnitTest_DipoleFromTable ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/DipoleFromTable_TEST.cxx)
package_add_test(UnitTest_HNLDISFromSpline ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/HNLDISFromSpline_TEST.cxx)
package_add_test(UnitTest_HNLDecay ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/HNLDecay_TEST.cxx)
#package_add_test(UnitTest_ElasticScattering ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/ElasticScattering_TEST.cxx)

pybind11_add_module(interactions ${PROJECT_SOURCE_DIR}/projects/interactions/private/pybindings/interactions.cxx)
target_link_libraries(interactions PRIVATE SIREN photospline rk_static crundec_static pybind11::embed GSL::gsl GSL::gslcblas)
if(TARGET PYTHIA8)
    target_link_libraries(interactions PRIVATE PYTHIA8)
    target_compile_definitions(interactions PRIVATE SIREN_HAS_PYTHIA8)
endif()
if(TARGET MARLEY)
    target_compile_definitions(interactions PRIVATE SIREN_HAS_MARLEY)
endif()
#pybind11_add_module(pyDarkNewsSerializer ${PROJECT_SOURCE_DIR}/projects/interactions/private/pybindings/pyDarkNewsSerializer.cxx)
#target_link_libraries(pyDarkNewsSerializer PRIVATE SIREN photospline rk_static pybind11::embed)
siren_strip_module(interactions)
if(DEFINED SKBUILD)
set_target_properties(interactions PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
package_add_test(UnitTest_CharmMesonDecay ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/CharmMesonDecay_TEST.cxx)

# Charm-DIS interaction-depth closure invariant (no Pythia8 needed: uses a mock
# that exercises the real base-class CrossSection::TotalCrossSectionAllFinalStates).
package_add_test(UnitTest_CharmDISClosure ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/CharmDISClosure_TEST.cxx)

# Real PythiaDISCrossSection charm-DIS closure/normalization regression test.
# Only built with Pythia8 support; needs charm-DIS spline files via env vars
# (SIREN_PYTHIA_TEST_DSDXDY, SIREN_PYTHIA_TEST_SIGMA), otherwise it SKIPs.
if(TARGET PYTHIA8)
    package_add_test(UnitTest_PythiaDISCharmClosure ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/PythiaDISCharmClosure_TEST.cxx)
    target_compile_definitions(UnitTest_PythiaDISCharmClosure PRIVATE SIREN_HAS_PYTHIA8)
    target_link_libraries(UnitTest_PythiaDISCharmClosure PYTHIA8)
endif()
