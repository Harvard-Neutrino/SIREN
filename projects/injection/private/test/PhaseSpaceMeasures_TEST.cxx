#include <gtest/gtest.h>

#include "PhaseSpaceTestFixtures.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/injection/FinalStateProbability.h"
#include "SIREN/injection/Weighter.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>

namespace {

using siren::dataclasses::InteractionRecord;
using siren::injection::MultiChannelPhaseSpace;
using siren::injection::PhaseSpaceChannel;
using siren::injection::PhaseSpaceMeasure;
using siren::injection::PhaseSpaceTopology;
using siren::injection::test::ConstantChannel;
using siren::injection::test::TwoBodyDecayRecord;

InteractionRecord ScatteringRecord(
    double beam_energy,
    double target_mass,
    double outgoing_mass,
    double recoil_mass)
{
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 0.0;
    record.primary_momentum = {beam_energy, 0.0, 0.0, beam_energy};
    record.target_mass = target_mass;
    record.interaction_vertex = {0.0, 0.0, -100.0};
    record.secondary_masses = {outgoing_mass, recoil_mass};
    record.secondary_momenta.resize(2);
    return record;
}

InteractionRecord BoostedAsymmetricTwoBodyDecayRecord() {
    constexpr double parent_mass = 4.0;
    constexpr double daughter_0_mass = 0.4;
    constexpr double daughter_1_mass = 1.7;
    constexpr double gamma = 1.5;
    double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    constexpr double daughter_1_cos_rest = 0.35;

    double p_rest = siren::injection::TwoBodyRestMomentum(
        parent_mass, daughter_1_mass, daughter_0_mass);
    double sin_rest = std::sqrt(
        1.0 - daughter_1_cos_rest * daughter_1_cos_rest);
    double px_1_rest = p_rest * sin_rest;
    double pz_1_rest = p_rest * daughter_1_cos_rest;
    double E_1_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_1_mass, daughter_0_mass);
    double E_0_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_0_mass, daughter_1_mass);

    auto boost = [beta](double E_rest, double px_rest, double pz_rest) {
        constexpr double gamma = 1.5;
        return std::array<double, 4>{
            gamma * (E_rest + beta * pz_rest),
            px_rest,
            0.0,
            gamma * (pz_rest + beta * E_rest)};
    };

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = parent_mass;
    record.primary_momentum = {
        gamma * parent_mass, 0.0, 0.0, gamma * beta * parent_mass};
    record.secondary_masses = {daughter_0_mass, daughter_1_mass};
    record.secondary_momenta = {
        boost(E_0_rest, -px_1_rest, -pz_1_rest),
        boost(E_1_rest, px_1_rest, pz_1_rest)};
    return record;
}

double DecayLabJacobian(
    InteractionRecord const & record,
    int daughter_index,
    double expected_cos_rest)
{
    int other_index = 1 - daughter_index;
    double parent_p = record.primary_momentum[3];
    double parent_E = record.primary_momentum[0];
    double beta = parent_p / parent_E;
    double gamma = parent_E / record.primary_mass;
    double daughter_mass = record.secondary_masses[daughter_index];
    double other_mass = record.secondary_masses[other_index];
    double p_rest = siren::injection::TwoBodyRestMomentum(
        record.primary_mass, daughter_mass, other_mass);
    double E_rest = siren::injection::TwoBodyRestEnergy(
        record.primary_mass, daughter_mass, other_mass);
    auto const & momentum = record.secondary_momenta[daughter_index];
    double p_lab = std::sqrt(
        momentum[1] * momentum[1] +
        momentum[2] * momentum[2] +
        momentum[3] * momentum[3]);
    double cos_lab = momentum[3] / p_lab;

    auto solutions = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, E_rest, daughter_mass, cos_lab);
    double best_jacobian = 0.0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (auto const & solution : solutions) {
        if (!solution.valid) continue;
        double distance = std::abs(
            solution.cos_theta_rest - expected_cos_rest);
        if (distance < best_distance) {
            best_distance = distance;
            best_jacobian = solution.jacobian;
        }
    }
    return best_jacobian;
}

InteractionRecord AsymmetricThreeBodyDecayRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.secondary_masses = {0.5, 1.0, 1.5};
    record.secondary_momenta = {
        {std::sqrt(0.5 * 0.5 + 0.7 * 0.7), 0.7, 0.0, 0.0},
        {std::sqrt(1.0 * 1.0 + 0.2 * 0.2 + 0.5 * 0.5),
         -0.2, 0.5, 0.0},
        {std::sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5),
         -0.5, -0.5, 0.0}};
    record.primary_mass =
        record.secondary_momenta[0][0] +
        record.secondary_momenta[1][0] +
        record.secondary_momenta[2][0];
    record.primary_momentum = {record.primary_mass, 0.0, 0.0, 0.0};
    return record;
}

class CountingConventionChannel final : public PhaseSpaceChannel {
public:
    mutable int topology_calls = 0;
    mutable int measure_calls = 0;

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override
    {
        return 1.0;
    }

    std::string Name() const override { return "CountingConvention"; }
    PhaseSpaceTopology Topology() const override {
        ++topology_calls;
        return PhaseSpaceTopology::Decay2Body;
    }
    PhaseSpaceMeasure Measure() const override {
        ++measure_calls;
        return PhaseSpaceMeasure::SolidAngleRest();
    }
};

siren::dataclasses::InteractionSignature SignatureWithSecondaries(size_t count) {
    siren::dataclasses::InteractionSignature signature;
    signature.secondary_types.assign(
        count, siren::dataclasses::ParticleType::unknown);
    return signature;
}

class MixedSignatureDecay final : public siren::interactions::Decay {
public:
    MixedSignatureDecay()
        : signatures_{SignatureWithSecondaries(2), SignatureWithSecondaries(3)} {}

    bool equal(siren::interactions::Decay const & other) const override {
        return dynamic_cast<MixedSignatureDecay const *>(&other) != nullptr;
    }
    double TotalDecayWidthAllFinalStates(InteractionRecord const &) const override {
        return 1.0;
    }
    double TotalDecayWidth(siren::dataclasses::ParticleType) const override {
        return 1.0;
    }
    double TotalDecayWidth(InteractionRecord const &) const override {
        return 1.0;
    }
    double DifferentialDecayWidth(InteractionRecord const &) const override {
        return 1.0;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override {}
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return signatures_;
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParent(
        siren::dataclasses::ParticleType) const override {
        return signatures_;
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return 1.0;
    }
    std::vector<std::string> DensityVariables() const override {
        return {"cos_theta"};
    }
    PhaseSpaceMeasure MeasureForSignature(
        siren::dataclasses::InteractionSignature const & signature) const override {
        return signature.secondary_types.size() == 2
            ? PhaseSpaceMeasure::SolidAngleRest()
            : PhaseSpaceMeasure::HelicityAngles();
    }

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
};

class MixedSignatureCrossSection final
    : public siren::interactions::CrossSection {
public:
    explicit MixedSignatureCrossSection(
        PhaseSpaceMeasure measure = PhaseSpaceMeasure::MandelstamQ2())
        : signatures_{SignatureWithSecondaries(2), SignatureWithSecondaries(3)}
        , measure_(measure) {}

    bool equal(siren::interactions::CrossSection const & other) const override {
        return dynamic_cast<MixedSignatureCrossSection const *>(&other) != nullptr;
    }
    double TotalCrossSection(InteractionRecord const &) const override {
        return 1.0;
    }
    double DifferentialCrossSection(InteractionRecord const &) const override {
        return 1.0;
    }
    double InteractionThreshold(InteractionRecord const &) const override {
        return 0.0;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override {}
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        return {siren::dataclasses::ParticleType::unknown};
    }
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(
        siren::dataclasses::ParticleType) const override {
        return GetPossibleTargets();
    }
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        return {siren::dataclasses::ParticleType::unknown};
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return signatures_;
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParents(
        siren::dataclasses::ParticleType,
        siren::dataclasses::ParticleType) const override {
        return signatures_;
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return 1.0;
    }
    std::vector<std::string> DensityVariables() const override {
        return {};
    }
    PhaseSpaceMeasure MeasureForSignature(
        siren::dataclasses::InteractionSignature const &) const override {
        return measure_;
    }

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
    PhaseSpaceMeasure measure_;
};

MultiChannelPhaseSpace TwoChannelMixture(std::vector<double> weights) {
    return MultiChannelPhaseSpace(
        {std::make_shared<ConstantChannel>(), std::make_shared<ConstantChannel>()},
        std::move(weights));
}

TEST(MultiChannelWeights, RejectsSizeMismatch) {
    EXPECT_THROW(TwoChannelMixture({1.0}), siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, NormalizesPositiveWeightSum) {
    auto mixture = TwoChannelMixture({0.2, 0.2});
    EXPECT_EQ(mixture.Weights(), (std::vector<double>{0.5, 0.5}));
}

TEST(MultiChannelWeights, RejectsNegativeAndNonFiniteWeights) {
    EXPECT_THROW(TwoChannelMixture({-0.1, 1.1}),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(TwoChannelMixture({
        std::numeric_limits<double>::quiet_NaN(), 1.0}),
        siren::utilities::ConfigurationError);
    EXPECT_THROW(TwoChannelMixture({
        std::numeric_limits<double>::infinity(), 0.0}),
        siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, AcceptsNormalizedWeights) {
    auto mixture = TwoChannelMixture({0.25, 0.75});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_NO_THROW(mixture.Sample(random, nullptr, record));
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
}

TEST(MultiChannelConventions, DensityEvaluatesEachConventionOnce) {
    auto first = std::make_shared<CountingConventionChannel>();
    auto second = std::make_shared<CountingConventionChannel>();
    MultiChannelPhaseSpace mixture({first, second});
    InteractionRecord record;

    for (int i = 0; i < 100; ++i) {
        EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    }
    // One observation per channel at construction and for every evaluation.
    EXPECT_EQ(first->topology_calls, 101);
    EXPECT_EQ(first->measure_calls, 101);
    EXPECT_EQ(second->topology_calls, 101);
    EXPECT_EQ(second->measure_calls, 101);
}

TEST(MultiChannelConventions, DensityUsesChangedConvention) {
    auto channel = std::make_shared<ConstantChannel>(
        2.0, PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2());
    MultiChannelPhaseSpace mixture({channel});
    InteractionRecord record;

    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 2.0);
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::MandelstamQ2());

    channel->SetConvention(PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::FixedMassY());
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::FixedMassY());
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 2.0);
}

TEST(MultiChannelConventions, ScatteringLabAngleMixtureIsFatalNotConvertible) {
    auto first = std::make_shared<ConstantChannel>(
        1.0, PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2());
    auto second = std::make_shared<ConstantChannel>(
        1.0, PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2());
    MultiChannelPhaseSpace mixture({first, second});
    second->SetConvention(PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::SolidAngleLab());

    auto diagnostics = mixture.ValidateChannelsDetailed();
    ASSERT_EQ(diagnostics.size(), 1u);
    EXPECT_EQ(diagnostics[0].severity,
              MultiChannelPhaseSpace::ChannelDiagnostic::Severity::Fatal);
    EXPECT_NE(diagnostics[0].message.find("not convertible"),
              std::string::npos);

    InteractionRecord record = ScatteringRecord(10.0, 1.0, 0.1, 1.0);
    auto random = std::make_shared<siren::utilities::SIREN_random>(8675309);
    EXPECT_THROW(mixture.Density(nullptr, record),
                 siren::utilities::MeasureCompatibilityError);
    EXPECT_THROW(mixture.Sample(random, nullptr, record),
                 siren::utilities::MeasureCompatibilityError);
}

TEST(PhysicalAdapterSignature, PinsDecayTopologyAndMeasure) {
    auto decay = std::make_shared<MixedSignatureDecay>();
    auto two_body = SignatureWithSecondaries(2);
    auto three_body = SignatureWithSecondaries(3);

    siren::injection::PhysicalDecayChannel unpinned(decay);
    EXPECT_EQ(unpinned.Topology(), PhaseSpaceTopology::Unspecified);
    EXPECT_EQ(unpinned.Measure(), PhaseSpaceMeasure::Unspecified());

    siren::injection::PhysicalDecayChannel pinned_two(decay, two_body);
    EXPECT_EQ(pinned_two.Topology(), PhaseSpaceTopology::Decay2Body);
    EXPECT_EQ(pinned_two.Measure(), PhaseSpaceMeasure::SolidAngleRest());

    siren::injection::PhysicalDecayChannel pinned_three(decay, three_body);
    EXPECT_EQ(pinned_three.Topology(), PhaseSpaceTopology::Decay3Body);
    EXPECT_EQ(pinned_three.Measure(), PhaseSpaceMeasure::HelicityAngles());
}

TEST(PhysicalAdapterSignature, PinsCrossSectionTopologyAndMeasure) {
    auto cross_section = std::make_shared<MixedSignatureCrossSection>();
    auto two_body = SignatureWithSecondaries(2);
    auto three_body = SignatureWithSecondaries(3);

    siren::injection::PhysicalCrossSectionChannel unpinned(cross_section);
    EXPECT_EQ(unpinned.Topology(), PhaseSpaceTopology::Unspecified);
    EXPECT_EQ(unpinned.Measure(), PhaseSpaceMeasure::Unspecified());

    siren::injection::PhysicalCrossSectionChannel pinned_two(
        cross_section, two_body);
    EXPECT_EQ(pinned_two.Topology(), PhaseSpaceTopology::Scatter2to2);
    EXPECT_EQ(pinned_two.Measure(), PhaseSpaceMeasure::MandelstamQ2());

    siren::injection::PhysicalCrossSectionChannel pinned_three(
        cross_section, three_body);
    EXPECT_EQ(pinned_three.Topology(), PhaseSpaceTopology::Scatter2to3);
    EXPECT_EQ(pinned_three.Measure(), PhaseSpaceMeasure::MandelstamQ2());
}

TEST(InteractionMeasureDeclaration, UndeclaredMeasuresAreUnspecified) {
    // Density-variable names such as "Bjorken x" do not declare a measure;
    // only an override does.
    auto signature = SignatureWithSecondaries(2);
    siren::interactions::DummyCrossSection undeclared;
    EXPECT_EQ(undeclared.Measure(), PhaseSpaceMeasure::Unspecified());
    EXPECT_EQ(undeclared.MeasureForSignature(signature),
              PhaseSpaceMeasure::Unspecified());

    MixedSignatureCrossSection declared(PhaseSpaceMeasure::FixedMassYPhi());
    EXPECT_EQ(declared.MeasureForSignature(signature),
              PhaseSpaceMeasure::FixedMassYPhi());
}

TEST(AzimuthTaxonomy, PredicatesAndCompletionsAgree) {
    using siren::dataclasses::MeasureHasExplicitAzimuth;
    using siren::dataclasses::MeasureIntegratesAzimuth;
    using siren::dataclasses::MeasureWithExplicitAzimuth;

    std::pair<PhaseSpaceMeasure, PhaseSpaceMeasure> lifts[] = {
        {PhaseSpaceMeasure::MandelstamQ2(), PhaseSpaceMeasure::MandelstamQ2Phi()},
        {PhaseSpaceMeasure::FixedMassY(), PhaseSpaceMeasure::FixedMassYPhi()},
        {PhaseSpaceMeasure::BjorkenXY(), PhaseSpaceMeasure::BjorkenXYPhi()},
        {PhaseSpaceMeasure::MandelstamQ2Y(), PhaseSpaceMeasure::MandelstamQ2YPhi()},
    };
    for (auto const & [marginal, joint] : lifts) {
        EXPECT_TRUE(MeasureIntegratesAzimuth(marginal));
        EXPECT_FALSE(MeasureHasExplicitAzimuth(marginal));
        EXPECT_TRUE(MeasureHasExplicitAzimuth(joint));
        EXPECT_FALSE(MeasureIntegratesAzimuth(joint));
        EXPECT_EQ(MeasureWithExplicitAzimuth(marginal), joint);
        EXPECT_EQ(MeasureWithExplicitAzimuth(joint), joint);
    }
    EXPECT_TRUE(MeasureIntegratesAzimuth(PhaseSpaceMeasure::CosThetaRest()));
    EXPECT_EQ(MeasureWithExplicitAzimuth(PhaseSpaceMeasure::CosThetaRest()),
              PhaseSpaceMeasure::SolidAngleRest());
    for (auto const & solid_angle :
         {PhaseSpaceMeasure::SolidAngleRest(), PhaseSpaceMeasure::SolidAngleLab()}) {
        EXPECT_TRUE(MeasureHasExplicitAzimuth(solid_angle));
        EXPECT_FALSE(MeasureIntegratesAzimuth(solid_angle));
        EXPECT_EQ(MeasureWithExplicitAzimuth(solid_angle), solid_angle);
        EXPECT_FALSE(siren::dataclasses::PhaseSpaceDensityConvertible(
            PhaseSpaceTopology::Decay2Body, solid_angle, PhaseSpaceMeasure::CosThetaRest()));
    }
    EXPECT_FALSE(MeasureHasExplicitAzimuth(PhaseSpaceMeasure::Unspecified()));
    EXPECT_FALSE(MeasureIntegratesAzimuth(PhaseSpaceMeasure::Unspecified()));
}

TEST(WeightingConvention, ElectsOneCommonConventionForBothDirections) {
    siren::injection::PhaseSpaceConvention marginal{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2()};
    siren::injection::PhaseSpaceConvention joint{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2Phi()};

    EXPECT_EQ(
        siren::injection::ResolveCommonFinalStateConvention(marginal, joint),
        joint);
    EXPECT_EQ(
        siren::injection::ResolveCommonFinalStateConvention(joint, marginal),
        joint);

    siren::injection::PhaseSpaceConvention different_family{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::BjorkenXYPhi()};
    EXPECT_THROW(
        siren::injection::ResolveCommonFinalStateConvention(
            marginal, different_family),
        siren::utilities::MeasureCompatibilityError);

    siren::injection::PhaseSpaceConvention different_topology{
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure::SolidAngleRest()};
    EXPECT_THROW(
        siren::injection::ResolveCommonFinalStateConvention(
            joint, different_topology),
        siren::utilities::MeasureCompatibilityError);
}

TEST(PhysicalChannelAdapters, AcceptExplicitSignatureConventionOverride) {
    auto scatter_signature = SignatureWithSecondaries(2);
    auto cross_section = std::make_shared<MixedSignatureCrossSection>();
    siren::injection::PhaseSpaceConvention scatter_convention{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2Phi()};

    siren::injection::PhysicalCrossSectionChannel scatter(
        cross_section, scatter_signature, scatter_convention);
    EXPECT_EQ(scatter.Topology(), scatter_convention.topology);
    EXPECT_EQ(scatter.Measure(), scatter_convention.measure);

    auto decay_signature = SignatureWithSecondaries(2);
    auto decay = std::make_shared<MixedSignatureDecay>();
    siren::injection::PhaseSpaceConvention decay_convention{
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure::SolidAngleLab(1)};

    siren::injection::PhysicalDecayChannel decay_channel(
        decay, decay_signature, decay_convention);
    EXPECT_EQ(decay_channel.Topology(), decay_convention.topology);
    EXPECT_EQ(decay_channel.Measure(), decay_convention.measure);

    siren::injection::PhaseSpaceConvention wrong_topology{
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure::SolidAngleRest()};
    EXPECT_THROW(
        siren::injection::PhysicalCrossSectionChannel(
            cross_section, scatter_signature, wrong_topology),
        siren::utilities::ConfigurationError);

    siren::injection::PhaseSpaceConvention unspecified_measure{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::Unspecified()};
    EXPECT_THROW(
        siren::injection::PhysicalCrossSectionChannel(
            cross_section, scatter_signature, unspecified_measure),
        siren::utilities::ConfigurationError);
}

TEST(WeightingConvention, LiftsNaturalFixedMassYIntoJointProposalMeasure) {
    auto signature = SignatureWithSecondaries(2);
    auto cross_section = std::make_shared<MixedSignatureCrossSection>(
        PhaseSpaceMeasure::FixedMassY());
    auto interactions =
        std::make_shared<siren::interactions::InteractionCollection>(
            siren::dataclasses::ParticleType::unknown,
            std::vector<std::shared_ptr<siren::interactions::CrossSection>>{
                cross_section});

    InteractionRecord record;
    record.signature = signature;
    record.primary_momentum = {1.0, 0.0, 0.0, 1.0};
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    siren::injection::PhysicalProcess process(
        siren::dataclasses::ParticleType::unknown, interactions);
    process.SetWeightingMode(siren::dataclasses::VertexWeightingMode::Fixed());
    siren::injection::PhaseSpaceConvention joint_convention;
    joint_convention.topology = PhaseSpaceTopology::Scatter2to2;
    joint_convention.measure = PhaseSpaceMeasure::FixedMassYPhi();

    EXPECT_NEAR(
        siren::injection::ProcessFinalStateProbability(
            detector, process, record, joint_convention),
        1.0 / (2.0 * M_PI), 1e-14);
}

std::shared_ptr<siren::interactions::InteractionCollection> CollectionOf(
    std::shared_ptr<siren::interactions::CrossSection> cross_section) {
    return std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::unknown,
        std::vector<std::shared_ptr<siren::interactions::CrossSection>>{cross_section});
}

std::shared_ptr<MultiChannelPhaseSpace> ConstantMixture(
    PhaseSpaceTopology topology, PhaseSpaceMeasure measure) {
    auto mixture = std::make_shared<MultiChannelPhaseSpace>(
        std::vector<std::shared_ptr<PhaseSpaceChannel>>{std::make_shared<ConstantChannel>(1.0, topology, measure)});
    return mixture;
}

TEST(ProcessPhaseSpaceValidation, SetPhaseSpaceChecksOnlyTheMixture) {
    // Registration cannot know which process the proposal will be weighted
    // against, so it accepts any internally consistent mixture.
    auto signature = SignatureWithSecondaries(2);
    siren::injection::PhysicalProcess process(
        siren::dataclasses::ParticleType::unknown,
        CollectionOf(std::make_shared<MixedSignatureCrossSection>(
            PhaseSpaceMeasure::FixedMassYPhi())));
    EXPECT_NO_THROW(process.SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::FixedMassY())));
    EXPECT_NO_THROW(process.SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::SolidAngleRest())));

    auto first = std::make_shared<ConstantChannel>();
    auto second = std::make_shared<ConstantChannel>();
    auto inconsistent = std::make_shared<MultiChannelPhaseSpace>(
        std::vector<std::shared_ptr<PhaseSpaceChannel>>{first, second});
    second->SetConvention(PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2Phi());
    EXPECT_THROW(process.SetPhaseSpace(signature, inconsistent),
                 siren::utilities::MeasureCompatibilityError);
}

TEST(ProcessPhaseSpaceValidation, WeighterRejectsProposalWithoutCommonMeasure) {
    // A (Q2, phi) proposal cannot be compared with a Bjorken (x, y) density.
    auto signature = SignatureWithSecondaries(2);
    auto interactions = CollectionOf(
        std::make_shared<MixedSignatureCrossSection>(PhaseSpaceMeasure::BjorkenXY()));
    auto injection = std::make_shared<siren::injection::PrimaryInjectionProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    auto physical = std::make_shared<siren::injection::PhysicalProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    injection->SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2Phi()));
    EXPECT_THROW(siren::injection::PrimaryProcessWeighter(physical, injection, nullptr),
                 siren::utilities::MeasureCompatibilityError);

    // Registered on both processes, the proposal replaces the model density on
    // both sides, so its chart need not match the model's.
    physical->SetPhaseSpace(signature, injection->GetPhaseSpace(signature));
    EXPECT_NO_THROW(siren::injection::PrimaryProcessWeighter(physical, injection, nullptr));
}

TEST(ProcessPhaseSpaceValidation, WeighterLiftsPerCosThetaProposal) {
    // A per-cos(theta) proposal is compared with a solid-angle model density
    // by lifting the proposal's declared uniform azimuth.
    auto signature = SignatureWithSecondaries(2);
    auto interactions = std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::unknown,
        std::vector<std::shared_ptr<siren::interactions::Decay>>{
            std::make_shared<MixedSignatureDecay>()});
    auto injection = std::make_shared<siren::injection::PrimaryInjectionProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    auto physical = std::make_shared<siren::injection::PhysicalProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    injection->SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::CosThetaRest()));
    EXPECT_NO_THROW(siren::injection::PrimaryProcessWeighter(physical, injection, nullptr));
}

TEST(CommonMeasure, UnspecifiedMajorityCannotOutvoteSpecifiedChannel) {
    auto first = std::make_shared<ConstantChannel>(2.0);
    auto second = std::make_shared<ConstantChannel>(4.0);
    auto third = std::make_shared<ConstantChannel>(8.0);
    MultiChannelPhaseSpace mixture({first, second, third}, {0.25, 0.25, 0.5});
    first->SetConvention(PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::Unspecified());
    second->SetConvention(PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::Unspecified());

    auto diagnostics = mixture.ValidateChannelsDetailed();
    ASSERT_EQ(diagnostics.size(), 2u);
    for (auto const & diagnostic : diagnostics) {
        EXPECT_NE(diagnostic.message.find("Unspecified which is not convertible to SolidAngleRest"),
                  std::string::npos);
    }
    InteractionRecord record;
    EXPECT_THROW(mixture.Density(nullptr, record), siren::utilities::MeasureCompatibilityError);
    EXPECT_THROW(mixture.DensityBreakdown(nullptr, record), siren::utilities::MeasureCompatibilityError);
}

TEST(CommonMeasure, ZeroUnspecifiedDensityCannotBypassMixedMeasureRejection) {
    auto first = std::make_shared<ConstantChannel>(8.0);
    auto second = std::make_shared<ConstantChannel>(0.0);
    MultiChannelPhaseSpace mixture({first, second});
    second->SetConvention(PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::Unspecified());

    InteractionRecord record;
    EXPECT_THROW(mixture.Density(nullptr, record), siren::utilities::MeasureCompatibilityError);
}

TEST(CommonMeasure, AllUnspecifiedChannelsRemainUnspecified) {
    MultiChannelPhaseSpace mixture({
        std::make_shared<ConstantChannel>(
            2.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified())}, {0.25, 0.75});

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::Unspecified());
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 3.5);
}

TEST(DecayMeasureConversion, UsesLabMeasureDaughterIndexInBothDirections) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    PhaseSpaceMeasure lab_0 = PhaseSpaceMeasure::SolidAngleLab(0);
    PhaseSpaceMeasure lab_1 = PhaseSpaceMeasure::SolidAngleLab(1);
    EXPECT_NE(lab_0, lab_1);

    double jacobian_0 = DecayLabJacobian(record, 0, -0.35);
    double jacobian_1 = DecayLabJacobian(record, 1, 0.35);
    ASSERT_GT(jacobian_0, 0.0);
    ASSERT_GT(jacobian_1, 0.0);
    ASSERT_GT(std::abs(jacobian_1 - jacobian_0), 1e-3);

    MultiChannelPhaseSpace lab_to_rest({
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_1)}, {0.5, 0.5});
    EXPECT_EQ(
        lab_to_rest.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_NEAR(
        lab_to_rest.Density(nullptr, record),
        0.5 * 3.0 + 0.5 * 5.0 * jacobian_1, 1e-13);

    MultiChannelPhaseSpace rest_to_lab({
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body, lab_1),
        std::make_shared<ConstantChannel>(
            6.0, PhaseSpaceTopology::Decay2Body, lab_1),
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest())}, {0.25, 0.25, 0.5});
    EXPECT_EQ(rest_to_lab.CommonMeasure(), lab_1);
    EXPECT_NEAR(
        rest_to_lab.Density(nullptr, record),
        0.25 * 4.0 + 0.25 * 6.0 + 0.5 * 3.0 / jacobian_1,
        1e-13);

    MultiChannelPhaseSpace lab_daughter_conversion({
        std::make_shared<ConstantChannel>(
            7.0, PhaseSpaceTopology::Decay2Body, lab_0),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_1)}, {0.5, 0.5});
    EXPECT_EQ(lab_daughter_conversion.CommonMeasure(), lab_0);
    EXPECT_NEAR(
        lab_daughter_conversion.Density(nullptr, record),
        0.5 * 7.0 + 0.5 * 5.0 * jacobian_1 / jacobian_0,
        1e-13);
}

TEST(DecayMeasureConversion, RejectsMissingMomentaAndInvalidDaughterIndex) {
    auto make_mixture = [](PhaseSpaceMeasure lab_measure) {
        MultiChannelPhaseSpace mixture({
            std::make_shared<ConstantChannel>(
                3.0, PhaseSpaceTopology::Decay2Body,
                PhaseSpaceMeasure::SolidAngleRest()),
            std::make_shared<ConstantChannel>(
                5.0, PhaseSpaceTopology::Decay2Body, lab_measure)}, {0.5, 0.5});
        return mixture;
    };

    InteractionRecord missing_momentum =
        BoostedAsymmetricTwoBodyDecayRecord();
    missing_momentum.secondary_momenta.resize(1);
    EXPECT_THROW(
        make_mixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, missing_momentum),
        std::runtime_error);

    InteractionRecord missing_mass = BoostedAsymmetricTwoBodyDecayRecord();
    missing_mass.secondary_masses.resize(1);
    EXPECT_THROW(
        make_mixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, missing_mass),
        std::runtime_error);

    InteractionRecord complete = BoostedAsymmetricTwoBodyDecayRecord();
    EXPECT_THROW(
        make_mixture(PhaseSpaceMeasure::SolidAngleLab(2)).Density(
            nullptr, complete),
        std::runtime_error);
}

namespace {

MultiChannelPhaseSpace RestPlusLabMixture(PhaseSpaceMeasure lab_measure) {
    MultiChannelPhaseSpace mixture({
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_measure)}, {0.5, 0.5});
    return mixture;
}

} // anonymous namespace

TEST(DecayMeasureConversion, ThrowsWhenLabAngleOutsideAllowedCone) {
    // Heavy daughters on a fast parent are confined to a narrow forward
    // cone; a backward daughter direction is kinematically impossible
    // for these masses, and the conversion must fail loudly.
    constexpr double parent_mass = 1.0;
    constexpr double daughter_mass = 0.45;
    constexpr double beta = 0.9;
    const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    constexpr double p_backward = 0.3;
    const double E_backward = std::sqrt(
        p_backward * p_backward + daughter_mass * daughter_mass);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = parent_mass;
    record.primary_momentum = {
        gamma * parent_mass, 0.0, 0.0, gamma * beta * parent_mass};
    record.secondary_masses = {daughter_mass, daughter_mass};
    record.secondary_momenta = {
        {E_backward, 0.0, 0.0, -p_backward},
        {0.0, 0.0, 0.0, 0.0}};

    EXPECT_THROW(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(0)).Density(
            nullptr, record),
        std::runtime_error);
}

TEST(DecayMeasureConversion, ThrowsOnDegenerateDaughterMomentum) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    record.secondary_momenta[1] = {record.secondary_masses[1], 0.0, 0.0, 0.0};
    EXPECT_THROW(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, record),
        std::runtime_error);
}

TEST(DecayMeasureConversion, ThrowsOnSubThresholdRecordMasses) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    record.primary_mass = 0.99 * (record.secondary_masses[0] +
                                  record.secondary_masses[1]);
    EXPECT_THROW(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, record),
        std::runtime_error);
}

TEST(DecayMeasureConversion, ParentAtRestConvertsAsIdentity) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    record.primary_momentum = {record.primary_mass, 0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, record),
        0.5 * 3.0 + 0.5 * 5.0);
}

TEST(ThreeBodyMeasureConversion, DalitzCrossFactorizationHasUnitJacobian) {
    namespace J = siren::injection::phase_space_jacobian;
    InteractionRecord record = AsymmetricThreeBodyDecayRecord();
    PhaseSpaceMeasure common = PhaseSpaceMeasure::DalitzPair(0, 1, 2);
    PhaseSpaceMeasure alternate = PhaseSpaceMeasure::DalitzPair(1, 0, 2);
    ASSERT_NE(common, alternate);

    auto pair_mass_squared = [&record](int first, int second) {
        auto const & p1 = record.secondary_momenta[first];
        auto const & p2 = record.secondary_momenta[second];
        double E = p1[0] + p2[0];
        double px = p1[1] + p2[1];
        double py = p1[2] + p2[2];
        double pz = p1[3] + p2[3];
        return E * E - px * px - py * py - pz * pz;
    };

    constexpr double alternate_density = 7.0;
    double old_intermediate = J::Recursive2BodyDensityToDalitzDensity(
        alternate_density, record.primary_mass,
        record.secondary_masses[alternate.spectator],
        record.secondary_masses[alternate.pair_first],
        record.secondary_masses[alternate.pair_second],
        pair_mass_squared(alternate.pair_first, alternate.pair_second));
    double old_converted = J::DalitzDensityToRecursive2BodyDensity(
        old_intermediate, record.primary_mass,
        record.secondary_masses[common.spectator],
        record.secondary_masses[common.pair_first],
        record.secondary_masses[common.pair_second],
        pair_mass_squared(common.pair_first, common.pair_second));
    ASSERT_GT(std::abs(old_converted - alternate_density), 1e-3);

    for (PhaseSpaceTopology topology : {
             PhaseSpaceTopology::Decay3Body,
             PhaseSpaceTopology::Scatter2to3}) {
        SCOPED_TRACE(siren::dataclasses::PhaseSpaceTopologyName(topology));
        MultiChannelPhaseSpace mixture({
            std::make_shared<ConstantChannel>(2.0, topology, common),
            std::make_shared<ConstantChannel>(
                alternate_density, topology, alternate)}, {0.25, 0.75});

        EXPECT_EQ(mixture.CommonMeasure(), common);
        EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 5.75);
        auto contributions = mixture.DensityBreakdown(nullptr, record);
        ASSERT_EQ(contributions.size(), 2u);
        EXPECT_DOUBLE_EQ(contributions[0], 0.5);
        EXPECT_DOUBLE_EQ(contributions[1], 5.25);
    }
}

TEST(ScatteringMeasureConversion, UsesIncomingAndOutgoingCmMomenta) {
    namespace J = siren::injection::phase_space_jacobian;
    constexpr double s = 25.0;
    constexpr double m_beam = 0.5;
    constexpr double m_target = 1.0;
    constexpr double m_outgoing = 2.0;
    constexpr double m_recoil = 0.75;

    double p_in_sq = siren::injection::Kallen(
        s, m_beam * m_beam, m_target * m_target) / (4.0 * s);
    double p_out_sq = siren::injection::Kallen(
        s, m_outgoing * m_outgoing, m_recoil * m_recoil) / (4.0 * s);
    double expected = 2.0 * std::sqrt(p_in_sq * p_out_sq);

    EXPECT_DOUBLE_EQ(J::SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target, m_outgoing, m_recoil), expected);
    EXPECT_NE(J::SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target), expected);
}

TEST(ScatteringMeasureConversion, MultiChannelUsesInelasticJacobian) {
    MultiChannelPhaseSpace mixture({
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())}, {0.5, 0.5});

    InteractionRecord record;
    record.primary_mass = 0.5;
    record.target_mass = 1.0;
    record.primary_momentum = {12.0, 0.0, 0.0, 0.0};
    record.secondary_masses = {2.0, 0.75};
    double s = record.primary_mass * record.primary_mass
             + record.target_mass * record.target_mass
             + 2.0 * record.target_mass * record.primary_momentum[0];
    double jacobian = siren::injection::phase_space_jacobian::
        SolidAngleRestToMandelstamQ2AbsJacobian(
            s, record.primary_mass, record.target_mass,
            record.secondary_masses[0], record.secondary_masses[1]);

    double expected = 0.5 * 3.0 + 0.5 * 5.0 * jacobian / (2.0 * M_PI);
    EXPECT_NEAR(mixture.Density(nullptr, record), expected, 1e-14);
}


TEST(ScatteringMeasureConversion, MixedFixedMassYAndQ2HasNoInverseYInflation) {
    constexpr double target_mass = 0.020;
    constexpr double incident_energy = 0.300;
    constexpr double y_density = 3.0;
    constexpr double q2_density = 5.0;
    double jacobian = 2.0 * target_mass * incident_energy;

    MultiChannelPhaseSpace mixture({
        std::make_shared<ConstantChannel>(
            y_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::FixedMassY()),
        std::make_shared<ConstantChannel>(
            q2_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())}, {0.5, 0.5});

    InteractionRecord record;
    record.target_mass = target_mass;
    record.primary_momentum = {incident_energy, 0.0, 0.0, incident_energy};
    double expected = 0.5 * y_density / jacobian + 0.5 * q2_density;

    record.interaction_parameters["bjorken_y"] = 0.2;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), expected);
    record.interaction_parameters["bjorken_y"] = 0.8;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), expected);
}

TEST(ScatteringMeasureConversion, ExplicitAzimuthWinsAndLiftsMarginal) {
    constexpr double marginal_density = 6.0;
    constexpr double joint_density = 2.0;

    MultiChannelPhaseSpace mixture({
        std::make_shared<ConstantChannel>(
            marginal_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            marginal_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            joint_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2Phi())}, {0.25, 0.25, 0.5});

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(),
              PhaseSpaceMeasure::MandelstamQ2Phi());
    EXPECT_NEAR(
        mixture.Density(nullptr, record),
        0.5 * marginal_density / (2.0 * M_PI) + 0.5 * joint_density,
        1e-14);

    siren::injection::PhaseSpaceConvention common = mixture.CommonConvention();
    EXPECT_EQ(common.topology, PhaseSpaceTopology::Scatter2to2);
    EXPECT_EQ(common.measure, PhaseSpaceMeasure::MandelstamQ2Phi());
    EXPECT_DOUBLE_EQ(mixture.DensityIn(nullptr, record, common),
                     mixture.Density(nullptr, record));

    siren::injection::PhaseSpaceConvention wrong_topology = common;
    wrong_topology.topology = PhaseSpaceTopology::Decay2Body;
    EXPECT_THROW(
        mixture.DensityIn(nullptr, record, wrong_topology),
        siren::utilities::MeasureCompatibilityError);
}

TEST(ScatteringMeasureConversion, ExplicitAzimuthCannotBeMarginalizedPointwise) {
    InteractionRecord record;
    EXPECT_THROW(
        siren::injection::ConvertDensity(
            1.0,
            PhaseSpaceMeasure::FixedMassYPhi(),
            PhaseSpaceMeasure::FixedMassY(),
            PhaseSpaceTopology::Scatter2to2,
            record),
        siren::utilities::MeasureCompatibilityError);
    EXPECT_THROW(
        siren::injection::ConvertDensity(
            0.0,
            PhaseSpaceMeasure::FixedMassYPhi(),
            PhaseSpaceMeasure::FixedMassY(),
            PhaseSpaceTopology::Scatter2to2,
            record),
        siren::utilities::MeasureCompatibilityError);
}

TEST(ScatteringMeasureConversion, BjorkenXYConvertsOnlyToQ2Y) {
    InteractionRecord record;
    record.target_mass = 0.938;
    record.primary_momentum = {5.0, 0.0, 0.0, 5.0};
    record.interaction_parameters["bjorken_y"] = 0.4;

    double converted = siren::injection::ConvertDensity(
        3.0,
        PhaseSpaceMeasure::BjorkenXY(),
        PhaseSpaceMeasure::MandelstamQ2Y(),
        PhaseSpaceTopology::Scatter2to2,
        record);
    EXPECT_GT(converted, 0.0);
    EXPECT_THROW(
        siren::injection::ConvertDensity(
            3.0,
            PhaseSpaceMeasure::BjorkenXY(),
            PhaseSpaceMeasure::MandelstamQ2(),
            PhaseSpaceTopology::Scatter2to2,
            record),
        siren::utilities::MeasureCompatibilityError);
}

TEST(ScatteringMeasureConversion, RejectsDecayStyleLabBoost) {
    EXPECT_THROW(
        MultiChannelPhaseSpace({
            std::make_shared<ConstantChannel>(
                1.0, PhaseSpaceTopology::Scatter2to2,
                PhaseSpaceMeasure::SolidAngleRest()),
            std::make_shared<ConstantChannel>(
                1.0, PhaseSpaceTopology::Scatter2to2,
                PhaseSpaceMeasure::SolidAngleLab())}),
        siren::utilities::MeasureCompatibilityError);
}

TEST(CoreMixture, InvalidComponentCannotBeMaskedByPositivePeer) {
    for (double bad : {-1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        MultiChannelPhaseSpace mixture({std::make_shared<ConstantChannel>(bad),
                                       std::make_shared<ConstantChannel>(10.0)});
        EXPECT_THROW(mixture.Density(nullptr, TwoBodyDecayRecord()),
                     siren::utilities::WeightCalculationError);
    }
    EXPECT_THROW(MultiChannelPhaseSpace({nullptr}), siren::utilities::ConfigurationError);
    EXPECT_THROW(MultiChannelPhaseSpace({std::make_shared<ConstantChannel>(),
        std::make_shared<ConstantChannel>()}, {1e308, 1e308}), siren::utilities::ConfigurationError);
}

TEST(CoreMeasure, RestCosThetaLiftsUniformAzimuthButCannotMarginalize) {
    auto integrated = PhaseSpaceMeasure::CosThetaRest();
    auto joint = PhaseSpaceMeasure::SolidAngleRest();
    auto topology = PhaseSpaceTopology::Decay2Body;
    auto record = TwoBodyDecayRecord();
    EXPECT_NEAR(siren::injection::ConvertDensity(0.5, integrated, joint, topology, record),
                1.0 / (4.0 * M_PI), 1e-16);
    EXPECT_TRUE(siren::injection::PhaseSpaceDensityConvertible(topology, integrated, joint));
    EXPECT_FALSE(siren::injection::PhaseSpaceDensityConvertible(topology, joint, integrated));
    EXPECT_THROW(siren::injection::ConvertDensity(0.5, joint, integrated, topology, record),
                 siren::utilities::MeasureCompatibilityError);
}

} // namespace
