#include <gtest/gtest.h>
#include "SIREN/distributions/primary/vertex/FixedTargetPositionDistribution.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Placement.h"

#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/primary/energy/Monoenergetic.h"
#include "SIREN/distributions/primary/energy/PowerLaw.h"
#include "SIREN/distributions/primary/direction/IsotropicDirection.h"

using namespace siren::distributions;

TEST(WeightableDistribution, OperatorLessSameType) {
    Monoenergetic a(10.0);
    Monoenergetic b(20.0);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
    EXPECT_FALSE(a < a);
}

TEST(WeightableDistribution, OperatorLessDifferentTypes) {
    Monoenergetic mono(10.0);
    PowerLaw power(2.0, 1.0, 100.0);

    bool result1 = mono < power;
    bool result2 = power < mono;
    EXPECT_NE(result1, result2);
}

TEST(WeightableDistribution, OperatorLessDifferentFamilies) {
    Monoenergetic mono(10.0);
    IsotropicDirection iso;

    bool result1 = mono < iso;
    bool result2 = iso < mono;
    EXPECT_NE(result1, result2);
}

TEST(WeightableDistribution, OperatorEqualSameType) {
    Monoenergetic a(10.0);
    Monoenergetic b(10.0);
    Monoenergetic c(20.0);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
}

TEST(WeightableDistribution, OperatorEqualDifferentTypes) {
    Monoenergetic mono(10.0);
    PowerLaw power(2.0, 1.0, 100.0);
    IsotropicDirection iso;
    EXPECT_FALSE(mono == power);
    EXPECT_FALSE(mono == iso);
}

TEST(WeightableDistribution, StrictWeakOrdering) {
    Monoenergetic a(10.0);
    Monoenergetic b(20.0);
    Monoenergetic c(30.0);
    EXPECT_FALSE(a < a);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
    EXPECT_TRUE(b < c);
    EXPECT_TRUE(a < c);
}

TEST(WeightableDistribution, FixedTargetIdentityIncludesBounds) {
    using namespace siren::geometry;
    Cylinder cylinder(Placement(), 0.05, 0.0, 0.1);
    auto fid = std::make_shared<Sphere>(Placement(), 1.0, 0.0);
    auto same_fid = std::make_shared<Sphere>(Placement(), 1.0, 0.0);
    auto other_fid = std::make_shared<Sphere>(Placement(), 2.0, 0.0);
    FixedTargetPositionDistribution a(cylinder, fid, 25.0);
    FixedTargetPositionDistribution same(cylinder, same_fid, 25.0);
    FixedTargetPositionDistribution length(cylinder, fid, 5.0);
    FixedTargetPositionDistribution volume(cylinder, other_fid, 25.0);
    FixedTargetPositionDistribution absent(cylinder, 25.0);
    EXPECT_TRUE(a == same);
    EXPECT_FALSE(a < same);
    EXPECT_FALSE(same < a);
    for(auto const * other : {&length, &volume, &absent}) {
        EXPECT_FALSE(a == *other);
        EXPECT_NE(a < *other, *other < a);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
