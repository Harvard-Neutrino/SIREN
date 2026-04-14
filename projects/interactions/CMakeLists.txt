
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
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CharmDISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/Hadronization.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CharmHadronization.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CharmMesonDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/CharmMesonDecay3Body.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/DMesonELoss.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/QuarkDISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/interactions/private/PythiaDISCrossSection.cxx

)
if(TARGET MARLEY)
    LIST (APPEND interactions_SOURCES
        ${PROJECT_SOURCE_DIR}/projects/interactions/private/MarleyCrossSection.cxx
    )
endif()
add_library(SIREN_interactions OBJECT ${interactions_SOURCES})
set_property(TARGET SIREN_interactions PROPERTY POSITION_INDEPENDENT_CODE ON)

# Pythia8 + LHAPDF for PythiaDISCrossSection
# Defaults target our cluster installs; override via -DPYTHIA8_DIR=... on cmake cmdline
set(PYTHIA8_DIR "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/miaochenjin/local/pythia8315"
    CACHE PATH "Pythia8 installation root")
set(LHAPDF_DIR "/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_8_x86_64"
    CACHE PATH "LHAPDF installation root")

target_include_directories(SIREN_interactions PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/interactions/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/interactions/>
    ${PYTHON_INCLUDE_DIRS}
    ${PYTHIA8_DIR}/include
    ${LHAPDF_DIR}/include
)

target_link_libraries(SIREN_interactions
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
        pybind11::embed
    PUBLIC
        photospline
        $<$<TARGET_EXISTS:MARLEY>:MARLEY>
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_detector
)
apply_siren_compile_options(SIREN_interactions)

target_link_directories(SIREN_interactions PUBLIC ${PYTHIA8_DIR}/lib ${LHAPDF_DIR}/lib)
target_link_libraries(SIREN_interactions PUBLIC pythia8 LHAPDF dl pthread)

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
target_link_directories(interactions PRIVATE ${PYTHIA8_DIR}/lib ${LHAPDF_DIR}/lib)
target_link_libraries(interactions PRIVATE SIREN photospline rk_static pybind11::embed pythia8 LHAPDF dl pthread)
if(TARGET MARLEY)
    target_compile_definitions(interactions PRIVATE SIREN_HAS_MARLEY)
endif()
#pybind11_add_module(pyDarkNewsSerializer ${PROJECT_SOURCE_DIR}/projects/interactions/private/pybindings/pyDarkNewsSerializer.cxx)
#target_link_libraries(pyDarkNewsSerializer PRIVATE SIREN photospline rk_static pybind11::embed)
if(DEFINED SKBUILD)
set_target_properties(interactions PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
package_add_test(UnitTest_CharmMesonDecay ${PROJECT_SOURCE_DIR}/projects/interactions/private/test/CharmMesonDecay_TEST.cxx)
