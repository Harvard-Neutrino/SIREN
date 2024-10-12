
#core files
LIST (APPEND utilities_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/utilities/private/Interpolator.cxx
    ${PROJECT_SOURCE_DIR}/projects/utilities/private/Random.cxx
    ${PROJECT_SOURCE_DIR}/projects/utilities/private/StringManipulation.cxx
)
add_library(SIREN_utilities OBJECT ${utilities_SOURCES})
set_property(TARGET SIREN_utilities PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_utilities PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/utilities/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/utilities/>
)

target_link_libraries(SIREN_utilities
    INTERFACE
        siren_compile_options
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/utilities/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

pybind11_add_module(utilities ${PROJECT_SOURCE_DIR}/projects/utilities/private/pybindings/utilities.cxx)
target_link_libraries(utilities PRIVATE SIREN)
if(DEFINED SKBUILD)
set_target_properties(utilities PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
