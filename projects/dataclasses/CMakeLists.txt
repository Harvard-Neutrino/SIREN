
#core files
LIST (APPEND dataclasses_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/Particle.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/InteractionSignature.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/InteractionRecord.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/InteractionTree.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/DecaySignature.cxx
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/DecayRecord.cxx
)
add_library(LI_dataclasses OBJECT ${dataclasses_SOURCES})
set_property(TARGET LI_dataclasses PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_dataclasses PUBLIC
    ${PROJECT_SOURCE_DIR}/projects/dataclasses/public/
)
add_dependencies(LI_dataclasses LI_serialization LI_utilities LI_math)

target_link_libraries(LI_dataclasses photospline)

target_include_directories(LI_dataclasses PUBLIC ${PROJECT_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_dataclasses PUBLIC ${PROJECT_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_dataclasses PUBLIC ${PROJECT_SOURCE_DIR}/projects/math/public/)

install(TARGETS LI_dataclasses
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/dataclasses/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

pybind11_add_module(dataclasses ${PROJECT_SOURCE_DIR}/projects/dataclasses/private/pybindings/dataclasses.cxx)
target_link_libraries(dataclasses PRIVATE LeptonInjector)
if(DEFINED SKBUILD)
set_target_properties(dataclasses PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
