
#core files
LIST (APPEND utilities_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/utilities/private/Interpolator.cxx
    ${PROJECT_SOURCE_DIR}/projects/utilities/private/Random.cxx
)
add_library(LI_utilities OBJECT ${utilities_SOURCES})
set_property(TARGET LI_utilities PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_utilities PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/utilities/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/utilities/>
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
