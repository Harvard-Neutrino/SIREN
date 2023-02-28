
#core files
LIST (APPEND math_SOURCES
    ${CMAKE_SOURCE_DIR}/projects/math/private/Vector3D.cxx
    ${CMAKE_SOURCE_DIR}/projects/math/private/Matrix3D.cxx
    ${CMAKE_SOURCE_DIR}/projects/math/private/Polynomial.cxx
    ${CMAKE_SOURCE_DIR}/projects/math/private/Quaternion.cxx
    ${CMAKE_SOURCE_DIR}/projects/math/private/EulerAngles.cxx
)
add_library(LI_math OBJECT ${math_SOURCES})
set_property(TARGET LI_math PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_math PUBLIC
    ${CMAKE_SOURCE_DIR}/projects/math/public/
)
add_dependencies(LI_math rk)
add_dependencies(LI_math LI_serialization)

target_link_libraries(LI_math photospline)
target_link_libraries(LI_math ${CMAKE_BINARY_DIR}/extern/rk/lib/${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX})

target_include_directories(LI_math PUBLIC ${CMAKE_BINARY_DIR}/extern/rk/include)
target_include_directories(LI_math PUBLIC ${CMAKE_SOURCE_DIR}/projects/serialization/public/)

install(TARGETS LI_math
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/projects/math/public/"
    DESTINATION "include"
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

