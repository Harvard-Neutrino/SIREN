
#core files
LIST (APPEND dataclasses_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/Particle.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/InteractionSignature.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/InteractionRecord.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/InteractionTree.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/DecaySignature.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/ParticleID.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/ParticleType.cxx
)
add_library(SIREN_dataclasses OBJECT ${dataclasses_SOURCES})
set_property(TARGET SIREN_dataclasses PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_dataclasses PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/dataclasses/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/dataclasses/>
)

target_link_libraries(SIREN_dataclasses
    INTERFACE
    siren_compile_options
    PUBLIC
    photospline
    SIREN_utilities
    SIREN_serialization
    SIREN_math
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/dataclasses/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

#package_add_test(UnitTest_InteractionSignature ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/test/InteractionSignature_TEST.cxx)
#package_add_test(UnitTest_InteractionRecord ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/test/InteractionRecord_TEST.cxx)
#package_add_test(UnitTest_InteractionTree ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/test/InteractionTree_TEST.cxx)

pybind11_add_module(dataclasses ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/pybindings/dataclasses.cxx)
target_link_libraries(dataclasses PRIVATE SIREN)
if(DEFINED SKBUILD)
set_target_properties(dataclasses PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
