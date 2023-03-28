
#core files
LIST (APPEND dataclasses_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/Particle.cxx
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/InteractionSignature.cxx
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/InteractionRecord.cxx
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/InteractionTree.cxx
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/DecaySignature.cxx
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/DecayRecord.cxx
)
add_library(LI_dataclasses OBJECT ${dataclasses_SOURCES})
set_property(TARGET LI_dataclasses PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_dataclasses PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/dataclasses/public/
)
add_dependencies(LI_dataclasses rk)
add_dependencies(LI_dataclasses LI_serialization LI_utilities LI_math)

target_link_libraries(LI_dataclasses photospline)
target_link_libraries(LI_dataclasses ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_dataclasses PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_dataclasses PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_dataclasses PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_dataclasses PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)

install(TARGETS LI_dataclasses
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/dataclasses/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

pybind11_add_module(dataclasses ${CMAKE_SOURCE_DIR}/projects/dataclasses/private/pybindings/dataclasses.cxx)
target_link_libraries(dataclasses PRIVATE LeptonInjector)
