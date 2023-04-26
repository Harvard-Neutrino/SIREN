
#core files
LIST (APPEND injection_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Process.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/InjectorBase.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/WeightingUtils.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/TreeWeighter.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/Weighter.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/ColumnDepthLeptonInjector.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/CylinderVolumeLeptonInjector.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/DecayRangeLeptonInjector.cxx
    ${PROJECT_SOURCE_DIR}/projects/injection/private/RangedLeptonInjector.cxx
)
add_library(LI_injection OBJECT ${injection_SOURCES})
set_property(TARGET LI_injection PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_injection PUBLIC
    ${PROJECT_SOURCE_DIR}/projects/injection/public/
)
add_dependencies(LI_injection rk)
add_dependencies(LI_injection LI_serialization LI_utilities LI_math LI_dataclasses LI_geometry LI_detector LI_crosssections LI_distributions)

target_link_libraries(LI_injection photospline)
target_link_libraries(LI_injection rk)

target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/extern/rk/include)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/dataclasses/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/geometry/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/detector/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/crosssections/public/)
target_include_directories(LI_injection PUBLIC ${PROJECT_SOURCE_DIR}/projects/distributions/public/)

install(TARGETS LI_injection
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/injection/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

#package_add_test(UnitTest_Injector ${PROJECT_SOURCE_DIR}/projects/injection/private/test/Injector_TEST.cxx)
package_add_test(UnitTest_CCM_HNL ${PROJECT_SOURCE_DIR}/projects/injection/private/test/CCM_HNL_TEST.cxx)

pybind11_add_module(injection ${PROJECT_SOURCE_DIR}/projects/injection/private/pybindings/injection.cxx)
target_link_libraries(injection PRIVATE LeptonInjector photospline)
if(DEFINED SKBUILD)
set_target_properties(injection PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
