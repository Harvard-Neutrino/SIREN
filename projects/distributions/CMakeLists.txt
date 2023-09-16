
#core files
LIST (APPEND distributions_SOURCES
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/Distributions.cxx

    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/direction/PrimaryDirectionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/direction/IsotropicDirection.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/direction/FixedDirection.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/direction/Cone.cxx

    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/energy/PrimaryEnergyDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/energy/ModifiedMoyalPlusExponentialEnergyDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/energy/TabulatedFluxDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/energy/PowerLaw.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/energy/Monoenergetic.cxx

    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/helicity/PrimaryNeutrinoHelicityDistribution.cxx

    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/type/PrimaryInjector.cxx

    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/RangeFunction.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/DecayRangeFunction.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/DepthFunction.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/LeptonDepthFunction.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/VertexPositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/PointSourcePositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/SecondaryPositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/RangePositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/DecayRangePositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/OrientedCylinderPositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/CylinderVolumePositionDistribution.cxx
    ${PROJECT_SOURCE_DIR}/projects/distributions/private/primary/vertex/ColumnDepthPositionDistribution.cxx

    ${PROJECT_SOURCE_DIR}/projects/distributions/private/target/momentum/TargetMomentumDistribution.cxx

)
add_library(LI_distributions OBJECT ${distributions_SOURCES})
set_property(TARGET LI_distributions PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(LI_distributions PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/projects/distributions/public/>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/LeptonInjector/distributions/>
)

target_link_libraries(LI_distributions
    PRIVATE
        $<BUILD_INTERFACE:rk_static>
    PUBLIC
        photospline
        LI_serialization
        LI_utilities
        LI_math
        LI_dataclasses
        LI_geometry
        LI_detector
        LI_crosssections
)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/projects/distributions/public/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.tcc"
)

package_add_test(UnitTest_PrimaryDirectionDistribution ${PROJECT_SOURCE_DIR}/projects/distributions/private/test/PrimaryDirectionDistribution_TEST.cxx)
package_add_test(UnitTest_PrimaryEnergyDistribution ${PROJECT_SOURCE_DIR}/projects/distributions/private/test/PrimaryEnergyDistribution_TEST.cxx)

pybind11_add_module(distributions ${PROJECT_SOURCE_DIR}/projects/distributions/private/pybindings/distributions.cxx)
target_link_libraries(distributions PRIVATE LeptonInjector)
if(DEFINED SKBUILD)
set_target_properties(distributions PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        LINK_FLAGS "-Wl,-rpath,$ORIGIN/../leptoninjector.libs/")
endif()
