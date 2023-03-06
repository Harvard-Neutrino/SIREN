
#core files
LIST (APPEND crosssections_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/CrossSection.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/CrossSectionCollection.cxx
    ${CMAKE_SOURCE_DIR}/projects/crosssections/private/DipoleFromTable.cxx
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

package_add_test(UnitTest_CrossSection ${CMAKE_SOURCE_DIR}/projects/crosssections/private/test/CrossSection_TEST.cxx)
