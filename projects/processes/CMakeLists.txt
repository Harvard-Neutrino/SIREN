
#core files
LIST (APPEND processes_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/processes/private/CrossSection.cxx
    ${CMAKE_SOURCE_DIR}/projects/processes/private/CrossSectionCollection.cxx
    ${CMAKE_SOURCE_DIR}/projects/processes/private/DipoleFromTable.cxx
    ${CMAKE_SOURCE_DIR}/projects/processes/private/DISFromSpline.cxx
    ${CMAKE_SOURCE_DIR}/projects/processes/private/ElasticScattering.cxx
)
add_library(LI_processes OBJECT ${processes_SOURCES})
set_property(TARGET LI_processes PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_processes PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/processes/public/
)
add_dependencies(LI_processes rk)
add_dependencies(LI_processes LI_serialization LI_utilities LI_math LI_dataclasses LI_detector)

target_link_libraries(LI_processes photospline)
target_link_libraries(LI_processes ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_processes PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_processes PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)
target_include_directories(LI_processes PUBLIC ${CMAKE_SOURCE_DIR}/projects/utilities/public/)
target_include_directories(LI_processes PUBLIC ${CMAKE_SOURCE_DIR}/projects/math/public/)
target_include_directories(LI_processes PUBLIC ${CMAKE_SOURCE_DIR}/projects/detector/public/)
target_include_directories(LI_processes PUBLIC ${CMAKE_SOURCE_DIR}/projects/dataclasses/public/)

install(TARGETS LI_processes
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/processes/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_DISFromSpline ${CMAKE_SOURCE_DIR}/projects/processes/private/test/DISFromSpline_TEST.cxx)
package_add_test(UnitTest_DipoleFromTable ${CMAKE_SOURCE_DIR}/projects/processes/private/test/DipoleFromTable_TEST.cxx)
package_add_test(UnitTest_ElasticScattering ${CMAKE_SOURCE_DIR}/projects/processes/private/test/ElasticScattering_TEST.cxx)
