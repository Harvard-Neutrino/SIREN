
#core files
LIST (APPEND injection_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/injection/private/InjectorBase.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/WeightingUtils.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/TreeWeighter.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/Weighter.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/ColumnDepthLeptonInjector.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/CylinderVolumeLeptonInjector.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/DecayRangeLeptonInjector.cxx
    ${CMAKE_SOURCE_DIR}/projects/injection/private/RangedLeptonInjector.cxx
)
add_library(LI_injection OBJECT ${injection_SOURCES})
set_property(TARGET LI_injection PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_injection PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/injection/public/
)
add_dependencies(LI_injection rk)
add_dependencies(LI_injection LI_serialization LI_utilities LI_math LI_dataclasses LI_geometry LI_detector LI_crosssections LI_distributions)

target_link_libraries(LI_injection photospline)
target_link_libraries(LI_injection ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_injection PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/dataclasses/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/geometry/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/detector/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/crosssections/public/)
target_include_directories(LI_injection PUBLIC ${CMAKE_SOURCE_DIR}/projects/distributions/public/)

install(TARGETS LI_injection
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/injection/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_Injector ${CMAKE_SOURCE_DIR}/projects/injection/private/test/Injector_TEST.cxx)
