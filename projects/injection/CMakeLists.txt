
#core files
LIST (APPEND injection_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Process.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Injector.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/WeightingUtils.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Weighter.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/TwoBodyKinematics.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/GeometryVolume.cxx
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
package_add_test(UnitTest_Weighter ${PROJECT_SOURCE_DIR}/projects/injection/private/test/Weighter_TEST.cxx)
package_add_test(UnitTest_PhaseSpaceReview ${PROJECT_SOURCE_DIR}/projects/injection/private/test/PhaseSpaceReview_TEST.cxx)
if(NOT CIBUILDWHEEL)
    target_link_libraries(UnitTest_PhaseSpaceReview pybind11::headers)
endif()
# Skip under cibuildwheel: package_add_test is a no-op there, so the following
# target_link_libraries would reference a target that was never created.
if(NOT CIBUILDWHEEL)
package_add_test(UnitTest_CharmSerialization ${PROJECT_SOURCE_DIR}/projects/injection/private/test/CharmSerialization_TEST.cxx)
# The injection headers pull in Pybind11Trampoline.h, so this TU needs the
# pybind11 headers (header-only; no Python embedding required).
target_link_libraries(UnitTest_CharmSerialization pybind11::headers)
endif()
package_add_test(UnitTest_TwoBodyKinematics ${PROJECT_SOURCE_DIR}/projects/injection/private/test/TwoBodyKinematics_TEST.cxx)
if(NOT CIBUILDWHEEL)
package_add_test(UnitTest_CCM_HNL ${PROJECT_SOURCE_DIR}/projects/injection/private/test/CCM_HNL_TEST.cxx)
target_link_libraries(UnitTest_CCM_HNL pybind11::embed)
endif()

pybind11_add_module(injection ${PROJECT_SOURCE_DIR}/projects/injection/private/pybindings/injection.cxx)
target_link_libraries(injection PRIVATE SIREN photospline pybind11::embed)
siren_strip_module(injection)
if(DEFINED SKBUILD)
set_target_properties(injection PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
