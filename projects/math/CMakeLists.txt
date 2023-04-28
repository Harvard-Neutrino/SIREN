
#core files
LIST (APPEND math_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/math/private/Vector3D.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Matrix3D.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Polynomial.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Quaternion.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/EulerAngles.cxx
    ${PROJECT_SOURCE_DIR}/projects/math/private/Interpolation.cxx
)
add_library(LI_math OBJECT ${math_SOURCES})
set_property(TARGET LI_math PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_math PUBLIC
    ${PROJECT_SOURCE_DIR}/projects/math/public/
)
add_dependencies(LI_math rk)
add_dependencies(LI_math LI_serialization)

if(${MACOSX})
target_link_libraries(LI_math PUBLIC
    photospline
    rk
    delabella
)
else()
target_link_libraries(LI_math PUBLIC
    photospline
    rk
    -Wl,--whole-archive
    delabella
    -Wl,--no-whole-archive
)
endif()

target_include_directories(LI_math PUBLIC ${PROJECT_SOURCE_DIR}/vendor/rk)
target_include_directories(LI_math PUBLIC ${PROJECT_SOURCE_DIR}/projects/serialization/public/)

install(TARGETS LI_math
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
target_link_libraries(math PRIVATE LeptonInjector)
if(DEFINED SKBUILD)
set_target_properties(math PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
