#include <cmath>
#include <memory>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/distributions/primary/vertex/ColumnDepthPositionDistribution.h"
#include "SIREN/distributions/primary/vertex/DepthFunction.h"
#include "SIREN/distributions/primary/vertex/PointSourcePositionDistribution.h"
#include "SIREN/distributions/primary/vertex/RangeFunction.h"
#include "SIREN/distributions/primary/vertex/RangePositionDistribution.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"

namespace {
using namespace siren;
using dataclasses::InteractionRecord;
using dataclasses::InteractionSignature;
using dataclasses::ParticleType;

constexpr double mass = 0.1;
constexpr double energy = 1.0;
constexpr double width = 2e-16;
constexpr double length = 25.0;
constexpr double radius = 3.0;

class ConstantDecay : public interactions::Decay {
public:
    mutable InteractionRecord observed;
    bool equal(Decay const & other) const override { return this == &other; }
    double TotalDecayWidthAllFinalStates(InteractionRecord const & record) const override {
        observed = record;
        return width;
    }
    double TotalDecayWidth(ParticleType) const override { return width; }
    double TotalDecayWidth(InteractionRecord const &) const override { return width; }
    double DifferentialDecayWidth(InteractionRecord const &) const override { return width; }
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &,
                         std::shared_ptr<utilities::SIREN_random>) const override {}
    std::vector<InteractionSignature> GetPossibleSignatures() const override { return {}; }
    std::vector<InteractionSignature> GetPossibleSignaturesFromParent(ParticleType) const override { return {}; }
    double FinalStateProbability(InteractionRecord const &) const override { return 1.0; }
    std::vector<std::string> DensityVariables() const override { return {}; }
};

class ZeroRange : public distributions::RangeFunction {
public:
    double operator()(ParticleType const &, double) const override { return 0.0; }
protected:
    bool equal(RangeFunction const &) const override { return true; }
    bool less(RangeFunction const &) const override { return false; }
};

class ZeroDepth : public distributions::DepthFunction {
public:
    double operator()(ParticleType const &, double) const override { return 0.0; }
protected:
    bool equal(DepthFunction const &) const override { return true; }
    bool less(DepthFunction const &) const override { return false; }
};

class PrimaryDecayPosition : public testing::TestWithParam<int> {};

TEST_P(PrimaryDecayPosition, PreservesKinematicsAndSamplesExponentialFlight) {
    std::shared_ptr<distributions::VertexPositionDistribution> position;
    if(GetParam() == 0) {
        position = std::make_shared<distributions::PointSourcePositionDistribution>(math::Vector3D(0, 0, 0), length);
    } else if(GetParam() == 1) {
        position = std::make_shared<distributions::RangePositionDistribution>(
            radius, length / 2, std::make_shared<ZeroRange>(), std::set<ParticleType>{});
    } else {
        position = std::make_shared<distributions::ColumnDepthPositionDistribution>(
            radius, length / 2, std::make_shared<ZeroDepth>());
    }
    auto detector = std::make_shared<detector::DetectorModel>();
    auto decay = std::make_shared<ConstantDecay>();
    auto collection = std::make_shared<interactions::InteractionCollection>(
        ParticleType::N4, std::vector<std::shared_ptr<interactions::Decay>>{decay});
    auto random = std::make_shared<utilities::SIREN_random>(17);
    double momentum = std::sqrt(energy * energy - mass * mass);
    double lambda = momentum / mass * utilities::Constants::hbarc / width;
    double normalization = -std::expm1(-length / lambda);
    double expected_mean = lambda - length * std::exp(-length / lambda) / normalization;
    double total_distance = 0.0;
    constexpr int samples = 4000;
    for(int i = 0; i < samples; ++i) {
        dataclasses::PrimaryDistributionRecord record(ParticleType::N4);
        record.SetMass(mass);
        record.SetEnergy(energy);
        record.SetDirection({0, 0, 1});
        record.SetHelicity(-0.5);
        position->Sample(random, detector, collection, record);
        ASSERT_EQ(decay->observed.signature.primary_type, ParticleType::N4);
        ASSERT_DOUBLE_EQ(decay->observed.primary_mass, mass);
        ASSERT_DOUBLE_EQ(decay->observed.primary_momentum[0], energy);
        ASSERT_NEAR(decay->observed.primary_momentum[3], momentum, 1e-14);
        ASSERT_DOUBLE_EQ(decay->observed.primary_helicity, -0.5);

        double distance = record.GetInteractionVertex()[2] - record.GetInitialPosition()[2];
        ASSERT_GT(distance, 0.0);
        ASSERT_LT(distance, length);
        total_distance += distance;

        InteractionRecord event;
        record.FinalizeAvailable(event);
        double density = std::exp(-distance / lambda) / (lambda * normalization);
        if(GetParam() != 0) density /= M_PI * radius * radius;
        EXPECT_NEAR(position->GenerationProbability(detector, collection, event), density, density * 1e-10);
    }
    EXPECT_NEAR(total_distance / samples, expected_mean, 0.5);
}

INSTANTIATE_TEST_CASE_P(AllPositionSamplers, PrimaryDecayPosition, testing::Values(0, 1, 2));
} // namespace
