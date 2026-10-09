#include <gtest/gtest.h>

#include "SIREN/distributions/primary/direction/Cone.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/detector/CartesianAxisExponentialDensityDistribution.h"
#include "SIREN/detector/CartesianAxisPolynomialDensityDistribution.h"

TEST(PublicHeaderCoexistence, DistinctTypesAreDeclared) {
    EXPECT_GT(sizeof(siren::distributions::Cone), 0u);
    EXPECT_GT(sizeof(siren::geometry::Cone), 0u);
    EXPECT_GT(sizeof(siren::detector::CartesianAxisExponentialDensityDistribution), 0u);
    EXPECT_GT(sizeof(siren::detector::CartesianAxisPolynomialDensityDistribution), 0u);
}
