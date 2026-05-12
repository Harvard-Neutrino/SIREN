#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "DistributionTest.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/utilities/Random.h"

using namespace siren::distributions;
using siren::geometry::Cylinder;


TEST(CylinderVolumePositionDistribution, Constructor) {
    Cylinder cylinder(10.0, 1.0, 4.0);
    CylinderVolumePositionDistribution dist(cylinder);
    EXPECT_EQ(dist.Name(), "CylinderVolumePositionDistribution");
}


TEST(CylinderVolumePositionDistribution, GenerationProbabilityZeroOutsideShell) {
    double inner = 1.0;
    double outer = 10.0;
    double height = 4.0;
    CylinderVolumePositionDistribution dist(Cylinder(outer, inner, height));

    siren::dataclasses::InteractionRecord record;
    record.interaction_vertex = {0.5, 0, 0};        // r < inner
    EXPECT_EQ(dist.GenerationProbability(nullptr, nullptr, record), 0.0);
    record.interaction_vertex = {15.0, 0, 0};       // r > outer
    EXPECT_EQ(dist.GenerationProbability(nullptr, nullptr, record), 0.0);
    record.interaction_vertex = {5.0, 0, height};   // |z| > h/2
    EXPECT_EQ(dist.GenerationProbability(nullptr, nullptr, record), 0.0);
}


TEST(CylinderVolumePositionDistribution, GenerationProbabilityConstantInsideShell) {
    double inner = 1.0;
    double outer = 10.0;
    double height = 4.0;
    CylinderVolumePositionDistribution dist(Cylinder(outer, inner, height));

    double V = M_PI * (outer*outer - inner*inner) * height;
    double expected = 1.0 / V;

    siren::dataclasses::InteractionRecord record;
    for(double r : {2.0, 5.0, 9.0}) {
        for(double z : {-1.5, 0.0, 1.5}) {
            record.interaction_vertex = {r, 0, z};
            EXPECT_NEAR(dist.GenerationProbability(nullptr, nullptr, record), expected, expected * 1e-9);
        }
    }
}


TEST(CylinderVolumePositionDistribution, SamplePositionInsideShell) {
    double inner = 1.0;
    double outer = 10.0;
    double height = 4.0;
    CylinderVolumePositionDistribution dist(Cylinder(outer, inner, height));

    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    size_t N = 10000;
    for(size_t i = 0; i < N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        record.SetDirection({0, 0, 1});
        dist.Sample(rand, nullptr, nullptr, record);
        auto v = record.GetInteractionVertex();
        double r = std::sqrt(v[0]*v[0] + v[1]*v[1]);
        EXPECT_GE(r, inner);
        EXPECT_LE(r, outer);
        EXPECT_GE(v[2], -height/2.0);
        EXPECT_LE(v[2],  height/2.0);
    }
}


// Tier 2: uniform-volume sampling. r^2 should be uniform on [inner^2, outer^2].
TEST(CylinderVolumePositionDistribution, SamplingIsUniformInRadius) {
    double inner = 1.0;
    double outer = 10.0;
    double height = 4.0;
    CylinderVolumePositionDistribution dist(Cylinder(outer, inner, height));

    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    size_t N = 100000;
    size_t n_bins = 20;
    DistributionTest test(inner*inner, outer*outer, n_bins);

    for(size_t i = 0; i < N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        record.SetDirection({0, 0, 1});
        dist.Sample(rand, nullptr, nullptr, record);
        auto v = record.GetInteractionVertex();
        double r2 = v[0]*v[0] + v[1]*v[1];
        test.AddValue(r2);
    }

    std::vector<double> expect(n_bins, double(N) / double(n_bins));
    EXPECT_TRUE(test.TestContents(3, expect));
}


// Tier 2: z should be uniform on [-h/2, h/2].
TEST(CylinderVolumePositionDistribution, SamplingIsUniformInZ) {
    double inner = 1.0;
    double outer = 10.0;
    double height = 4.0;
    CylinderVolumePositionDistribution dist(Cylinder(outer, inner, height));

    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    size_t N = 100000;
    size_t n_bins = 20;
    DistributionTest test(-height/2.0, height/2.0, n_bins);

    for(size_t i = 0; i < N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        record.SetDirection({0, 0, 1});
        dist.Sample(rand, nullptr, nullptr, record);
        auto v = record.GetInteractionVertex();
        test.AddValue(v[2]);
    }

    std::vector<double> expect(n_bins, double(N) / double(n_bins));
    EXPECT_TRUE(test.TestContents(3, expect));
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
