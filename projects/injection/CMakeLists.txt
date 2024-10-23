
#core files
LIST (APPEND injection_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Process.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Injector.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/WeightingUtils.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Weighter.cxx
)
add_library(SIREN_injection OBJECT ${injection_SOURCES})
set_property(TARGET SIREN_injection PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_injection PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/injection/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/injection/>
)

target_link_libraries(SIREN_injection
    PUBLIC
        photospline
        SIREN_utilities
        SIREN_serialization
        SIREN_math
        SIREN_dataclasses
        SIREN_geometry
        SIREN_detector
        SIREN_interactions
        SIREN_distributions
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
        pybind11::embed
)
apply_siren_compile_options(SIREN_injection)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/injection/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

#package_add_test(UnitTest_Injector ${PROJECT_SOURCE_DIR}/projects/injection/private/test/Injector_TEST.cxx)
if(NOT ${CIBUILDWHEEL})
package_add_test(UnitTest_CCM_HNL ${PROJECT_SOURCE_DIR}/projects/injection/private/test/CCM_HNL_TEST.cxx)
target_link_libraries(UnitTest_CCM_HNL pybind11::embed)
endif()

pybind11_add_module(injection ${PROJECT_SOURCE_DIR}/projects/injection/private/pybindings/injection.cxx)
target_link_libraries(injection PRIVATE SIREN photospline pybind11::embed)
if(DEFINED SKBUILD)
set_target_properties(injection PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
