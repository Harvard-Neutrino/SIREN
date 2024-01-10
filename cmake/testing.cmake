mark_as_advanced(
	BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
	gmock_build_tests gtest_build_samples gtest_build_tests
	gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)


set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

macro(package_add_test TESTNAME)
	add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main LeptonInjector)
    add_dependencies(${TESTNAME} rk_static)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
	set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

ENABLE_TESTING()

##package_add_test(UnitTest_Medium Medium_TEST.cxx)
#package_add_test(UnitTest_Axis Axis_TEST.cxx)
#package_add_test(UnitTest_Distribution Distribution_TEST.cxx)
#package_add_test(UnitTest_EarthModel EarthModel_TEST.cxx)
#package_add_test(UnitTest_ExtrPoly ExtrPoly_TEST.cxx)
#package_add_test(UnitTest_Geometry Geometry_TEST.cxx)
#package_add_test(UnitTest_Vector3D Vector3D_TEST.cxx)
##package_add_test(UnitTest_Sector Sector_TEST.cxx)
#package_add_test(UnitTest_MaterialModel MaterialModel_TEST.cxx)
#package_add_test(UnitTest_MathMethods MathMethods_TEST.cxx)
#package_add_test(UnitTest_Path Path_TEST.cxx)
#package_add_test(UnitTest_Polynomial Polynomial_TEST.cxx)
#package_add_test(UnitTest_DensityDistribution DensityDistribution_TEST.cxx)
#
#package_add_test(UnitTest_Quaternion Quaternion_TEST.cxx)
#
#package_add_test(UnitTest_DipoleTable DipoleTable_TEST.cxx)
#package_add_test(UnitTest_CrossSection CrossSection_TEST.cxx)
#package_add_test(UnitTest_Injector Injector_TEST.cxx)
#package_add_test(UnitTest_ElasticScattering ElasticScattering_TEST.cxx)
##package_add_test(UnitTest_IceCube IceCube_TEST.cxx)
