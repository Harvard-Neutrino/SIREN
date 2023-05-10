
#core files
LIST (APPEND crosssections_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/Decay.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/CrossSection.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/CrossSectionCollection.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/DipoleFromTable.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/NeutrissimoDecay.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/DISFromSpline.cxx
    ${PROJECT_SOURCE_DIR}/projects/crosssections/private/ElasticScattering.cxx
)
add_library(LI_crosssections OBJECT ${crosssections_SOURCES})
set_property(TARGET LI_crosssections PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_crosssections PUBLIC
    ${PROJECT_SOURCE_DIR}/projects/crosssections/public/
)
add_dependencies(LI_crosssections rk)
add_dependencies(LI_crosssections LI_serialization LI_utilities LI_math LI_dataclasses LI_detector)

target_link_libraries(LI_crosssections photospline)
target_link_libraries(LI_crosssections rk)

target_include_directories(LI_crosssections PUBLIC ${PROJECT_SOURCE_DIR}/vendor/rk/include)
target_include_directories(LI_crosssections PUBLIC ${PROJECT_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_crosssections PUBLIC ${PROJECT_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_crosssections PUBLIC ${PROJECT_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_crosssections PUBLIC ${PROJECT_SOURCE_DIR}/projects/detector/public/)
target_include_directories(LI_crosssections PUBLIC ${PROJECT_SOURCE_DIR}/projects/dataclasses/public/)

install(TARGETS LI_crosssections
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/crosssections/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_DISFromSpline ${PROJECT_SOURCE_DIR}/projects/crosssections/private/test/DISFromSpline_TEST.cxx)
package_add_test(UnitTest_DipoleFromTable ${PROJECT_SOURCE_DIR}/projects/crosssections/private/test/DipoleFromTable_TEST.cxx)
#package_add_test(UnitTest_ElasticScattering ${PROJECT_SOURCE_DIR}/projects/crosssections/private/test/ElasticScattering_TEST.cxx)

pybind11_add_module(crosssections ${PROJECT_SOURCE_DIR}/projects/crosssections/private/pybindings/crosssections.cxx)
target_link_libraries(crosssections PRIVATE LeptonInjector photospline rk)
if(DEFINED SKBUILD)
set_target_properties(crosssections PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
