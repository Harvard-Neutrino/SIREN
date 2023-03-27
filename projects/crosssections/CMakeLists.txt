
#core files
LIST (APPEND crosssections_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/Decay.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/CrossSection.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/CrossSectionCollection.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/DipoleFromTable.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/NeutrissimoDecay.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/DISFromSpline.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/ElasticScattering.cxx
)
add_library(LI_crosssections OBJECT ${crosssections_SOURCES})
set_property(TARGET LI_crosssections PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_crosssections PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/crosssections/public/
)
add_dependencies(LI_crosssections rk)
add_dependencies(LI_crosssections LI_serialization LI_utilities LI_math LI_dataclasses LI_detector)

target_link_libraries(LI_crosssections photospline)
target_link_libraries(LI_crosssections ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_crosssections PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/detector/public/)
target_include_directories(LI_crosssections PUBLIC ${CMAKE_SOURCE_DIR}/projects/dataclasses/public/)

install(TARGETS LI_crosssections
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/crosssections/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_DISFromSpline ${CMAKE_SOURCE_DIR}/projects/crosssections/private/test/DISFromSpline_TEST.cxx)
package_add_test(UnitTest_DipoleFromTable ${CMAKE_SOURCE_DIR}/projects/crosssections/private/test/DipoleFromTable_TEST.cxx)
#package_add_test(UnitTest_ElasticScattering ${CMAKE_SOURCE_DIR}/projects/crosssections/private/test/ElasticScattering_TEST.cxx)

pybind11_add_module(CrossSections ${CMAKE_SOURCE_DIR}/projects/crosssections/private/pybindings/crosssections.cxx)
target_link_libraries(CrossSections PRIVATE LI_crosssections)
