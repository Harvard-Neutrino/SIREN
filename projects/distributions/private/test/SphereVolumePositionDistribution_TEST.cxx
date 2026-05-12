#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "DistributionTest.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/distributions/primary/vertex/SphereVolumePositionDistribution.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/utilities/Random.h"

using namespace siren::distributions;
using siren::geometry::Sphere;


TEST(SphereVolumePositionDistribution, Constructor) {
    Sphere sphere(10.0, 1.0);
    SphereVolumePositionDistribution dist(sphere);
    EXPECT_EQ(dist.Name(), "SphereVolumePositionDistribution");
}


TEST(SphereVolumePositionDistribution, GenerationProbabilityZeroOutsideShell) {
    double inner = 1.0;
    double outer = 10.0;
    SphereVolumePositionDistribution dist(Sphere(outer, inner));

    siren::dataclasses::InteractionRecord record;
    record.interaction_vertex = {0.5, 0, 0};   // r < inner
    EXPECT_EQ(dist.GenerationProbability(nullptr, nullptr, record), 0.0);
    record.interaction_vertex = {15.0, 0, 0};  // r > outer
    EXPECT_EQ(dist.GenerationProbability(nullptr, nullptr, record), 0.0);
}


TEST(SphereVolumePositionDistribution, GenerationProbabilityConstantInsideShell) {
    double inner = 1.0;
    double outer = 10.0;
    SphereVolumePositionDistribution dist(Sphere(outer, inner));

    double V = (4.0/3.0) * M_PI * (outer*outer*outer - inner*inner*inner);
    double expected = 1.0 / V;

    siren::dataclasses::InteractionRecord record;
    for(double r : {2.0, 5.0, 9.0}) {
        record.interaction_vertex = {r, 0, 0};
        EXPECT_NEAR(dist.GenerationProbability(nullptr, nullptr, record), expected, expected * 1e-9);
    }
}


TEST(SphereVolumePositionDistribution, SamplePositionInsideShell) {
    double inner = 1.0;
    double outer = 10.0;
    SphereVolumePositionDistribution dist(Sphere(outer, inner));

    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    size_t N = 10000;
    for(size_t i = 0; i < N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        record.SetDirection({0, 0, 1});
        dist.Sample(rand, nullptr, nullptr, record);
        auto v = record.GetInteractionVertex();
        double r = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        EXPECT_GE(r, inner);
        EXPECT_LE(r, outer);
    }
}


// Tier 2: uniform-volume sampling. r^3 should be uniform on [inner^3, outer^3].
TEST(SphereVolumePositionDistribution, SamplingIsUniformInVolume) {
    double inner = 1.0;
    double outer = 10.0;
    SphereVolumePositionDistribution dist(Sphere(outer, inner));

    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    size_t N = 100000;
    size_t n_bins = 20;
    DistributionTest test(inner*inner*inner, outer*outer*outer, n_bins);

    for(size_t i = 0; i < N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        record.SetDirection({0, 0, 1});
        dist.Sample(rand, nullptr, nullptr, record);
        auto v = record.GetInteractionVertex();
        double r = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        test.AddValue(r * r * r);
    }

    std::vector<double> expect(n_bins, double(N) / double(n_bins));
    EXPECT_TRUE(test.TestContents(3, expect));
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
