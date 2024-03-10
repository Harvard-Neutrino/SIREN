
#core files
LIST (APPEND math_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/math/private/Vector3D.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Matrix3D.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Polynomial.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Quaternion.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/EulerAngles.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Interpolation.cxx
)
add_library(SIREN_math OBJECT ${math_SOURCES})
set_property(TARGET SIREN_math PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(SIREN_math PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/math/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SIREN/math/>
)

if(${MACOSX})
target_link_libraries(SIREN_math
    PUBLIC
        photospline
        delabella_shared
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
        SIREN_serialization
)
else()
target_link_libraries(SIREN_math
    PUBLIC
        photospline
        delabella_shared
    PRIVATE
        -Wl,--whole-archive
        $<BUILD_INTERFACE:rk_static>
        -Wl,--no-whole-archive
        SIREN_serialization
)
endif()

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/math/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_MathMethods ${PROJECT_SOURCE_DIR}/projects/math/private/test/MathMethods_TEST.cxx)
package_add_test(UnitTest_Polynomial ${PROJECT_SOURCE_DIR}/projects/math/private/test/Polynomial_TEST.cxx)
package_add_test(UnitTest_Quaternion ${PROJECT_SOURCE_DIR}/projects/math/private/test/Quaternion_TEST.cxx)
package_add_test(UnitTest_EulerAngles ${PROJECT_SOURCE_DIR}/projects/math/private/test/EulerAngles_TEST.cxx)
package_add_test(UnitTest_Vector3D ${PROJECT_SOURCE_DIR}/projects/math/private/test/Vector3D_TEST.cxx)
package_add_test(UnitTest_Transform ${PROJECT_SOURCE_DIR}/projects/math/private/test/Transform_TEST.cxx)
package_add_test(UnitTest_InterpolationOperators ${PROJECT_SOURCE_DIR}/projects/math/private/test/InterpolationOperators_TEST.cxx)
package_add_test(UnitTest_Indexing ${PROJECT_SOURCE_DIR}/projects/math/private/test/Indexing_TEST.cxx)
package_add_test(UnitTest_Interpolator ${PROJECT_SOURCE_DIR}/projects/math/private/test/Interpolator_TEST.cxx)

pybind11_add_module(math ${PROJECT_SOURCE_DIR}/projects/math/private/pybindings/math.cxx)
target_link_libraries(math PRIVATE SIREN)
if(DEFINED SKBUILD)
set_target_properties(math PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,\\\$ORIGIN/../siren.libs/")
endif()
