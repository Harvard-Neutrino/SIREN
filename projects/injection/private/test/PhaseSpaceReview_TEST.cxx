#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

#include "../DetectorDirectedChannelUtils.h"
#include "../InteractionRecordUtils.h"

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

// Reconstruct a reference overlap cap to identify the part of a known lens
// that a simple cone-intersection sampler would miss.
std::pair<siren::math::Vector3D, double> LegacyOverlapCap(
    double theta_kin,
    double theta_bound,
    double axis_sep,
    siren::math::Vector3D const & kin_axis,
    siren::math::Vector3D const & bound_axis)
{
    double sin_kin = std::sin(theta_kin);
    double cos_kin = std::cos(theta_kin);
    double sin_sep = std::sin(axis_sep);
    double cos_sep = std::cos(axis_sep);
    double cos_phi = (std::cos(theta_bound) - cos_kin * cos_sep)
                   / (sin_kin * sin_sep);
    cos_phi = std::clamp(cos_phi, -1.0, 1.0);
    double sin_phi = std::sqrt(1.0 - cos_phi * cos_phi);

    siren::math::Vector3D in_plane = bound_axis - kin_axis * cos_sep;
    in_plane.normalize();
    siren::math::Vector3D out_plane =
        siren::math::vector_product(kin_axis, in_plane);
    out_plane.normalize();

    siren::math::Vector3D p1 = kin_axis * cos_kin
        + in_plane * (sin_kin * cos_phi)
        + out_plane * (sin_kin * sin_phi);
    siren::math::Vector3D p2 = kin_axis * cos_kin
        + in_plane * (sin_kin * cos_phi)
        - out_plane * (sin_kin * sin_phi);
    siren::math::Vector3D midpoint = kin_axis * cos_kin
        + in_plane * (sin_kin * cos_phi);
    midpoint.normalize();

    double theta_deep = axis_sep - theta_bound;
    siren::math::Vector3D deepest = theta_deep > 1e-15
        ? kin_axis * std::cos(theta_deep)
            + in_plane * std::sin(theta_deep)
        : kin_axis;
    siren::math::Vector3D center = midpoint + deepest;
    center.normalize();
    double cos_min = std::min({
        siren::math::scalar_product(center, p1),
        siren::math::scalar_product(center, p2),
        siren::math::scalar_product(center, midpoint),
        siren::math::scalar_product(center, deepest)});
    double half_angle = std::min(
        1.2 * std::acos(std::clamp(cos_min, -1.0, 1.0)), M_PI);
    return {center, std::cos(half_angle)};
}

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

InteractionRecord TwoBodyDecayRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.5};
    record.secondary_momenta = {
        {11.0, 12.0, 13.0, 14.0},
        {21.0, 22.0, 23.0, 24.0}};
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

InteractionRecord ThresholdThreeBodyScatteringRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 0.1;
    record.primary_momentum = {0.1, 0.0, 0.0, 0.0};
    record.target_mass = 0.1;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.2, 0.1, 0.1};
    record.secondary_momenta = {
        {11.0, 12.0, 13.0, 14.0},
        {21.0, 22.0, 23.0, 24.0},
        {31.0, 32.0, 33.0, 34.0}};
    return record;
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


class ConstantChannel final : public PhaseSpaceChannel {
public:
    explicit ConstantChannel(
        double density = 1.0,
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure measure = PhaseSpaceMeasure::SolidAngleRest())
        : density_(density), topology_(topology), measure_(measure) {}

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override
    {
        return density_;
    }

    std::string Name() const override { return "Constant"; }
    PhaseSpaceTopology Topology() const override { return topology_; }
    PhaseSpaceMeasure Measure() const override { return measure_; }

private:
    double density_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;
};

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

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
};

class MixedSignatureCrossSection final
    : public siren::interactions::CrossSection {
public:
    MixedSignatureCrossSection()
        : signatures_{SignatureWithSecondaries(2), SignatureWithSecondaries(3)} {}

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
        return {"q2"};
    }

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
};

MultiChannelPhaseSpace TwoChannelMixture(std::vector<double> weights) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(),
        std::make_shared<ConstantChannel>()};
    mixture.weights = std::move(weights);
    return mixture;
}

TEST(MultiChannelWeights, RejectsSizeMismatch) {
    auto mixture = TwoChannelMixture({1.0});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_THROW(mixture.Sample(random, nullptr, record),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(mixture.Density(nullptr, record),
                 siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, RejectsUnnormalizedWeights) {
    auto mixture = TwoChannelMixture({0.2, 0.2});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_THROW(mixture.Sample(random, nullptr, record),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(mixture.Density(nullptr, record),
                 siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, RejectsNegativeAndNonFiniteWeights) {
    InteractionRecord record;

    EXPECT_THROW(TwoChannelMixture({-0.1, 1.1}).Density(nullptr, record),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(TwoChannelMixture({
        std::numeric_limits<double>::quiet_NaN(), 1.0}).Density(nullptr, record),
        siren::utilities::ConfigurationError);
    EXPECT_THROW(TwoChannelMixture({
        std::numeric_limits<double>::infinity(), 0.0}).Density(nullptr, record),
        siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, AcceptsNormalizedWeights) {
    auto mixture = TwoChannelMixture({0.25, 0.75});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_NO_THROW(mixture.Sample(random, nullptr, record));
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
}

TEST(MultiChannelConventions, CompatibleMixtureCachesSuccessfulValidation) {
    auto first = std::make_shared<CountingConventionChannel>();
    auto second = std::make_shared<CountingConventionChannel>();
    MultiChannelPhaseSpace mixture;
    mixture.channels = {first, second};
    mixture.weights = {0.5, 0.5};
    InteractionRecord record;

    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    int first_topology_calls = first->topology_calls;
    int first_measure_calls = first->measure_calls;
    int second_topology_calls = second->topology_calls;
    int second_measure_calls = second->measure_calls;

    // The change-detection fingerprint reads each channel's conventions by
    // value (pointer identity alone can collide when a freed channel's
    // address is reused), so every evaluation costs exactly one Topology and
    // one Measure probe per channel. A cached validation adds no second
    // probe; only a rebuild does.
    for (int i = 0; i < 100; ++i) {
        EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    }
    EXPECT_EQ(first->topology_calls, first_topology_calls + 100);
    EXPECT_EQ(first->measure_calls, first_measure_calls + 100);
    EXPECT_EQ(second->topology_calls, second_topology_calls + 100);
    EXPECT_EQ(second->measure_calls, second_measure_calls + 100);

    // Public channel replacement changes the fingerprint and rebuilds the
    // cached conventions on the next evaluation: the fingerprint probe plus
    // the rebuild probe, and nothing more.
    auto replacement = std::make_shared<CountingConventionChannel>();
    mixture.channels[1] = replacement;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    EXPECT_EQ(first->topology_calls, first_topology_calls + 102);
    EXPECT_EQ(first->measure_calls, first_measure_calls + 102);
    EXPECT_EQ(replacement->topology_calls, 2);
    EXPECT_EQ(replacement->measure_calls, 2);
}

TEST(MultiChannelConventions, ClearedAndReallocatedChannelRebuildsCache) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {std::make_shared<ConstantChannel>(
        2.0, PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2())};
    mixture.weights = {1.0};
    InteractionRecord record;

    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 2.0);
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::MandelstamQ2());

    // Destroy-then-reallocate channel swap: the replacement may be allocated
    // at the freed channel's address, so a pointer-only fingerprint can
    // serve the stale measures. The value-folded fingerprint must rebuild
    // regardless of where the allocator places the new channel.
    mixture.channels.clear();
    mixture.channels.push_back(std::make_shared<ConstantChannel>(
        4.0, PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::FixedMassY()));

    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::FixedMassY());
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 4.0);
}

TEST(MultiChannelConventions, ScatteringLabAngleMixtureIsFatalNotConvertible) {
    // ConvertDensity implements the rest<->lab boost only for Decay2Body, so
    // a Scatter2to2 mixture containing a SolidAngleLab channel must be
    // reported as a fatal incompatibility by validation, not promised as an
    // auto-conversion that then throws at evaluation time.
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleLab())};
    mixture.weights = {0.5, 0.5};

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

TEST(CommonMeasure, UnspecifiedMajorityCannotOutvoteSpecifiedChannel) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            2.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            8.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest())};
    mixture.weights = {0.25, 0.25, 0.5};

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_THROW(mixture.Density(nullptr, record), std::runtime_error);
    EXPECT_THROW(mixture.DensityBreakdown(nullptr, record), std::runtime_error);
}

TEST(CommonMeasure, ZeroUnspecifiedDensityCannotBypassMixedMeasureRejection) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            8.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            0.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_THROW(mixture.Density(nullptr, record), std::runtime_error);
}

TEST(CommonMeasure, AllUnspecifiedChannelsRemainUnspecified) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            2.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified())};
    mixture.weights = {0.25, 0.75};

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

    MultiChannelPhaseSpace lab_to_rest;
    lab_to_rest.channels = {
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_1)};
    lab_to_rest.weights = {0.5, 0.5};
    EXPECT_EQ(
        lab_to_rest.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_NEAR(
        lab_to_rest.Density(nullptr, record),
        0.5 * 3.0 + 0.5 * 5.0 * jacobian_1, 1e-13);

    MultiChannelPhaseSpace rest_to_lab;
    rest_to_lab.channels = {
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body, lab_1),
        std::make_shared<ConstantChannel>(
            6.0, PhaseSpaceTopology::Decay2Body, lab_1),
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest())};
    rest_to_lab.weights = {0.25, 0.25, 0.5};
    EXPECT_EQ(rest_to_lab.CommonMeasure(), lab_1);
    EXPECT_NEAR(
        rest_to_lab.Density(nullptr, record),
        0.25 * 4.0 + 0.25 * 6.0 + 0.5 * 3.0 / jacobian_1,
        1e-13);

    MultiChannelPhaseSpace lab_daughter_conversion;
    lab_daughter_conversion.channels = {
        std::make_shared<ConstantChannel>(
            7.0, PhaseSpaceTopology::Decay2Body, lab_0),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_1)};
    lab_daughter_conversion.weights = {0.5, 0.5};
    EXPECT_EQ(lab_daughter_conversion.CommonMeasure(), lab_0);
    EXPECT_NEAR(
        lab_daughter_conversion.Density(nullptr, record),
        0.5 * 7.0 + 0.5 * 5.0 * jacobian_1 / jacobian_0,
        1e-13);
}

TEST(DecayMeasureConversion, RejectsMissingMomentaAndInvalidDaughterIndex) {
    auto make_mixture = [](PhaseSpaceMeasure lab_measure) {
        MultiChannelPhaseSpace mixture;
        mixture.channels = {
            std::make_shared<ConstantChannel>(
                3.0, PhaseSpaceTopology::Decay2Body,
                PhaseSpaceMeasure::SolidAngleRest()),
            std::make_shared<ConstantChannel>(
                5.0, PhaseSpaceTopology::Decay2Body, lab_measure)};
        mixture.weights = {0.5, 0.5};
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
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_measure)};
    mixture.weights = {0.5, 0.5};
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
        MultiChannelPhaseSpace mixture;
        mixture.channels = {
            std::make_shared<ConstantChannel>(2.0, topology, common),
            std::make_shared<ConstantChannel>(
                alternate_density, topology, alternate)};
        mixture.weights = {0.25, 0.75};

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
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())};
    mixture.weights = {0.5, 0.5};

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

TEST(ScatteringMeasureConversion, FixedMassYUsesTwoMEWithoutYFactor) {
    namespace J = siren::injection::phase_space_jacobian;
    constexpr double target_mass = 0.020;
    constexpr double incident_energy = 0.300;
    constexpr double y_density = 3.0;
    constexpr double q2_density = 5.0;
    double jacobian = 2.0 * target_mass * incident_energy;

    EXPECT_DOUBLE_EQ(J::FixedMassYToMandelstamQ2AbsJacobian(
                         target_mass, incident_energy),
                     jacobian);
    EXPECT_DOUBLE_EQ(J::FixedMassYDensityToMandelstamQ2Density(
                         y_density, target_mass, incident_energy),
                     y_density / jacobian);
    EXPECT_DOUBLE_EQ(J::MandelstamQ2DensityToFixedMassYDensity(
                         q2_density, target_mass, incident_energy),
                     q2_density * jacobian);

    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirectedScatteringChannel energy_loss_channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::BjorkenY,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    siren::injection::DetectorDirectedScatteringChannel recoil_channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::RecoilY,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    EXPECT_EQ(energy_loss_channel.Measure(), PhaseSpaceMeasure::FixedMassYPhi());
    EXPECT_EQ(recoil_channel.Measure(), PhaseSpaceMeasure::FixedMassYPhi());
}

TEST(ScatteringMeasureConversion, MixedFixedMassYAndQ2HasNoInverseYInflation) {
    constexpr double target_mass = 0.020;
    constexpr double incident_energy = 0.300;
    constexpr double y_density = 3.0;
    constexpr double q2_density = 5.0;
    double jacobian = 2.0 * target_mass * incident_energy;

    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            y_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::FixedMassY()),
        std::make_shared<ConstantChannel>(
            q2_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    record.target_mass = target_mass;
    record.primary_momentum = {incident_energy, 0.0, 0.0, incident_energy};
    double expected = 0.5 * y_density / jacobian + 0.5 * q2_density;

    record.interaction_parameters["bjorken_y"] = 0.2;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), expected);
    record.interaction_parameters["bjorken_y"] = 0.8;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), expected);
}

TEST(ScatteringMeasureConversion, RejectsDecayStyleLabBoost) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleLab())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    EXPECT_THROW(mixture.Density(nullptr, record), std::runtime_error);
}

TEST(NestedKleissPittau, ConvertsInnerCreditsToOuterMeasure) {
    namespace J = siren::injection::phase_space_jacobian;
    auto inner = std::make_shared<MultiChannelPhaseSpace>();
    inner->channels = {
        std::make_shared<ConstantChannel>(
            2.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            8.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())};
    inner->weights = {0.5, 0.5};

    MultiChannelPhaseSpace outer;
    outer.channels = {
        std::make_shared<siren::injection::NestedMixtureChannel>(inner),
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            7.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest())};
    outer.weights = {0.4, 0.3, 0.3};

    InteractionRecord record;
    record.primary_mass = 0.5;
    record.target_mass = 1.0;
    record.primary_momentum = {12.0, 0.0, 0.0, 0.0};
    record.secondary_masses = {2.0, 0.75};
    double s = record.primary_mass * record.primary_mass
             + record.target_mass * record.target_mass
             + 2.0 * record.target_mass * record.primary_momentum[0];
    double q2_to_rest = J::MandelstamQ2PhiDensityToSolidAngleRestDensity(
        1.0 / (2.0 * M_PI), s, record.primary_mass, record.target_mass,
        record.secondary_masses[0], record.secondary_masses[1]);
    ASSERT_GT(std::abs(q2_to_rest - 1.0), 0.1);

    double g_outer = 0.4 * (5.0 * q2_to_rest) + 0.3 * 3.0 + 0.3 * 7.0;
    EXPECT_NEAR(outer.Density(nullptr, record), g_outer, 1e-14);

    constexpr double event_weight = 1.5;
    constexpr double w2 = event_weight * event_weight;
    outer.Accumulate(nullptr, record, event_weight, false, true);

    ASSERT_EQ(outer.kp_count_, 1);
    ASSERT_EQ(outer.kp_accumulator_.size(), 3u);
    EXPECT_NEAR(
        outer.kp_accumulator_[0], w2 * 5.0 * q2_to_rest / g_outer,
        1e-14);
    EXPECT_NEAR(outer.kp_accumulator_[1], w2 * 3.0 / g_outer, 1e-14);
    EXPECT_NEAR(outer.kp_accumulator_[2], w2 * 7.0 / g_outer, 1e-14);

    ASSERT_EQ(inner->kp_count_, 1);
    ASSERT_EQ(inner->kp_accumulator_.size(), 2u);
    EXPECT_NEAR(
        inner->kp_accumulator_[0], w2 * 2.0 * q2_to_rest / g_outer,
        1e-14);
    EXPECT_NEAR(
        inner->kp_accumulator_[1], w2 * 8.0 * q2_to_rest / g_outer,
        1e-14);
}

TEST(DirectedScatteringCone, SampleAndDensityHaveIdenticalSupport) {
    auto target = siren::geometry::Box(10.0, 1.0, 1.0).create();
    siren::injection::DetectorDirectedScatteringChannel channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    auto random = std::make_shared<siren::utilities::SIREN_random>(20260709);
    InteractionRecord record = ScatteringRecord(20.0, 1.0, 0.05, 1.0);

    int misses = 0;
    int zero_density = 0;
    for (int i = 0; i < 2000; ++i) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        siren::math::Vector3D direction(
            record.secondary_momenta[0][1],
            record.secondary_momenta[0][2],
            record.secondary_momenta[0][3]);
        direction.normalize();
        siren::math::Vector3D vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);
        if (!siren::injection::detail::DirectionHitsTarget(
                *target, vertex, direction)) {
            ++misses;
        }
        if (channel.Density(nullptr, record) <= 0.0) ++zero_density;
    }

    EXPECT_GT(misses, 0);
    EXPECT_EQ(zero_density, 0);
}

TEST(DirectedScatteringFallback, BackwardRecoilSightlineUsesUniformCM) {
    constexpr double beam_energy = 20.0;
    constexpr double target_mass = 1.0;
    constexpr double outgoing_mass = 0.05;
    constexpr double recoil_mass = 1.0;
    constexpr double total_energy = beam_energy + target_mass;
    constexpr double s = total_energy * total_energy
                       - beam_energy * beam_energy;

    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, -100.0)),
        5.0, 0.0).create();
    siren::injection::DetectorDirectedScatteringChannel channel(
        target, 1,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);

    InteractionRecord record = ScatteringRecord(
        beam_energy, target_mass, outgoing_mass, recoil_mass);
    record.interaction_vertex = {0.0, 0.0, 0.0};
    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        total_energy, 0.0, 0.0, beam_energy,
        std::sqrt(s), recoil_mass, outgoing_mass,
        siren::math::Vector3D(0.0, 0.0, 0.0), *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Disjoint);
    EXPECT_FALSE(channel.DirectingActive(record));

    double p_in_cm_sq = siren::injection::Kallen(
        s, 0.0, target_mass * target_mass) / (4.0 * s);
    double p_out_cm_sq = siren::injection::Kallen(
        s, recoil_mass * recoil_mass, outgoing_mass * outgoing_mass)
        / (4.0 * s);
    double expected_density =
        1.0 / (8.0 * M_PI * std::sqrt(p_in_cm_sq * p_out_cm_sq));

    auto random = std::make_shared<siren::utilities::SIREN_random>(424242);
    for (int i = 0; i < 1000; ++i) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        EXPECT_GT(record.secondary_momenta[1][3], 0.0);
        EXPECT_NEAR(channel.Density(nullptr, record), expected_density,
                    1e-12 * expected_density);
    }
}

TEST(DirectedScatteringBranches, DensityIncludesRootSelectionProbability) {
    constexpr double E1 = 1.0;
    constexpr double m1 = 0.0;
    constexpr double m2 = 0.1;
    constexpr double m3 = 0.1;
    constexpr double m4 = 0.1;
    constexpr double cos_theta = 0.328;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    siren::math::Vector3D direction(sin_theta, 0.0, cos_theta);

    double total_e = E1 + m2;
    double s = total_e * total_e - E1 * E1;
    double C = E1 * cos_theta;
    double A = 0.5 * (s + m3 * m3 - m4 * m4);
    double a = total_e * total_e - C * C;
    double b = -2.0 * A * C;
    double c = total_e * total_e * m3 * m3 - A * A;
    double discriminant = b * b - 4.0 * a * c;
    ASSERT_GT(discriminant, 0.0);
    double p3 = (-b + std::sqrt(discriminant)) / (2.0 * a);
    double E3 = std::sqrt(p3 * p3 + m3 * m3);

    siren::geometry::Placement placement(direction * 100.0);
    auto target = siren::geometry::Sphere(placement, 5.0, 0.0).create();
    siren::injection::DetectorDirectedScatteringChannel channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = m1;
    record.primary_momentum = {E1, 0.0, 0.0, E1};
    record.target_mass = m2;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {m3, m4};
    record.secondary_momenta = {
        {E3, p3 * sin_theta, 0.0, p3 * cos_theta},
        {total_e - E3, -p3 * sin_theta, 0.0, E1 - p3 * cos_theta}};

    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        total_e, 0.0, 0.0, E1, std::sqrt(s), m3, m4,
        siren::math::Vector3D(0.0, 0.0, 0.0), *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Overlap);
    double angular_density = 1.0 / geometry.omega_eff;
    double root_s = std::sqrt(s);
    double p_rest = siren::injection::TwoBodyRestMomentum(root_s, m3, m4);
    double E_rest = siren::injection::TwoBodyRestEnergy(root_s, m3, m4);
    auto solutions = siren::injection::SolveLabAngle(
        E1 / total_e, total_e / root_s,
        p_rest, E_rest, m3, cos_theta);
    double root_jacobian = 0.0;
    for (auto const & solution : solutions) {
        if (solution.valid &&
            std::abs(solution.p_lab - p3) < 1e-10) {
            root_jacobian = solution.jacobian;
        }
    }
    ASSERT_GT(root_jacobian, 0.0);
    double q2_jacobian = siren::injection::phase_space_jacobian::
        SolidAngleRestToMandelstamQ2AbsJacobian(
            s, m1, m2, m3, m4);
    double expected = 0.5 * angular_density
                    * root_jacobian / q2_jacobian;

    EXPECT_NEAR(channel.Density(nullptr, record), expected,
                1e-10 * std::max(1.0, expected));
}

TEST(DirectedScatteringQ2Mapping, DownscatteringRetainsNegativeQ2Support) {
    constexpr double m1 = 0.02;
    constexpr double m2 = 0.01;
    constexpr double m3 = 0.01;
    constexpr double m4 = 0.02;
    constexpr double E1 = 0.04;
    const double p1 = std::sqrt(E1 * E1 - m1 * m1);
    const double s = m1 * m1 + m2 * m2 + 2.0 * m2 * E1;
    const double root_s = std::sqrt(s);
    const double E1cm = (s + m1 * m1 - m2 * m2) / (2.0 * root_s);
    const double E3cm = (s + m3 * m3 - m4 * m4) / (2.0 * root_s);
    const double p1cm = std::sqrt(E1cm * E1cm - m1 * m1);
    const double p3cm = std::sqrt(E3cm * E3cm - m3 * m3);
    const double base = m1 * m1 + m3 * m3 - 2.0 * E1cm * E3cm;
    const double q2min = -(base + 2.0 * p1cm * p3cm);
    const double q2max = -(base - 2.0 * p1cm * p3cm);
    ASSERT_LT(q2min, 0.0);
    ASSERT_GT(q2max, 0.0);

    auto target = siren::geometry::Box(1.0, 1.0, 1.0).create();
    siren::injection::DetectorDirectedScatteringChannel channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        siren::injection::DetectorDirectedScatteringChannel::Q2Mode::Tabulated,
        0.0, {-0.001, 0.001}, {0.0, 1.0});

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = m1;
    record.primary_momentum = {E1, 0.0, 0.0, p1};
    record.target_mass = m2;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {m3, m4};
    record.secondary_momenta.resize(2);

    const double expected_density =
        1.0 / (2.0 * M_PI * (q2max - q2min));
    auto random = std::make_shared<siren::utilities::SIREN_random>(577215);
    int negative_q2 = 0;
    for (int i = 0; i < 512; ++i) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        double q2 = record.interaction_parameters.at("Q2");
        EXPECT_GE(q2, q2min - 1e-15);
        EXPECT_LE(q2, q2max + 1e-15);
        if (q2 < -1e-8) ++negative_q2;
        EXPECT_NEAR(channel.Density(nullptr, record), expected_density,
                    1e-11 * expected_density);
    }

    EXPECT_GT(negative_q2, 50);
}

TEST(PropagatorMapping, RejectsWindowContainingThePole) {
    constexpr double m_med2 = 0.0025;
    EXPECT_THROW(
        siren::injection::PropagatorMapping(m_med2, -0.007, 0.01),
        std::invalid_argument);
    EXPECT_THROW(
        siren::injection::PropagatorMapping(m_med2, -m_med2, 0.01),
        std::invalid_argument);
    EXPECT_NO_THROW(
        siren::injection::PropagatorMapping(m_med2, -0.002, 0.01));
    EXPECT_NO_THROW(
        siren::injection::PropagatorMapping(m_med2, -0.01, -0.004));
}

TEST(PropagatorMapping, BelowPoleWindowRemainsNormalized) {
    siren::injection::PropagatorMapping map(0.0025, -0.01, -0.004);
    for (double r : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        double x = map.Forward(r);
        EXPECT_GE(x, -0.01 - 1e-15);
        EXPECT_LE(x, -0.004 + 1e-15);
        EXPECT_NEAR(map.Inverse(x), r, 1e-9);
    }
    constexpr int n = 2000;
    constexpr double lo = -0.01;
    constexpr double hi = -0.004;
    constexpr double h = (hi - lo) / n;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        double weight = (i == 0 || i == n) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
        sum += weight * map.Density(lo + i * h);
    }
    EXPECT_NEAR(sum * h / 3.0, 1.0, 1e-6);
}

TEST(DirectedScatteringQ2Mapping, PropagatorPoleInsideWindowIsExcluded) {
    // Pins the propagator proposal to the window above the t-channel pole
    // when the pole sits inside the kinematic Q^2 range.
    constexpr double m1 = 1.2;
    constexpr double m2 = 1.0;
    constexpr double m3 = 1.0;
    constexpr double m4 = 1.1;
    constexpr double E1 = 5.0;
    constexpr double mediator_mass = 0.05;
    constexpr double m_med2 = mediator_mass * mediator_mass;
    const double p1 = std::sqrt(E1 * E1 - m1 * m1);
    const double s = m1 * m1 + m2 * m2 + 2.0 * m2 * E1;
    const double root_s = std::sqrt(s);
    const double E1cm = (s + m1 * m1 - m2 * m2) / (2.0 * root_s);
    const double E3cm = (s + m3 * m3 - m4 * m4) / (2.0 * root_s);
    const double p1cm = std::sqrt(E1cm * E1cm - m1 * m1);
    const double p3cm = std::sqrt(E3cm * E3cm - m3 * m3);
    const double base = m1 * m1 + m3 * m3 - 2.0 * E1cm * E3cm;
    const double q2min = -(base + 2.0 * p1cm * p3cm);
    const double q2max = -(base - 2.0 * p1cm * p3cm);
    ASSERT_LT(q2min, -m_med2);   // the pole sits inside the raw window
    ASSERT_GT(q2max, 0.0);

    // Mirror of the channel's restriction: proposal support above the pole.
    const double q2_floor = -m_med2 + 1e-9 * (q2max + m_med2);
    siren::injection::PropagatorMapping reference(
        m_med2, std::max(q2min, q2_floor), q2max);

    auto target = siren::geometry::Box(1.0, 1.0, 1.0).create();
    siren::injection::DetectorDirectedScatteringChannel channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        siren::injection::DetectorDirectedScatteringChannel::Q2Mode::Propagator,
        mediator_mass);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = m1;
    record.primary_momentum = {E1, 0.0, 0.0, p1};
    record.target_mass = m2;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {m3, m4};
    record.secondary_momenta.resize(2);

    auto random = std::make_shared<siren::utilities::SIREN_random>(141421);
    std::set<double> distinct;
    for (int i = 0; i < 512; ++i) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        double q2 = record.interaction_parameters.at("Q2");
        EXPECT_GT(q2, -m_med2);
        EXPECT_LE(q2, q2max + 1e-12);
        distinct.insert(q2);
        double density = channel.Density(nullptr, record);
        ASSERT_TRUE(std::isfinite(density));
        ASSERT_GT(density, 0.0);
        // The density is steep near the pole; the tolerance covers the
        // Q^2 reconstruction roundoff between the stored parameter and
        // the record's momenta.
        double expected_density = reference.Density(q2) / (2.0 * M_PI);
        EXPECT_NEAR(density, expected_density,
                    1e-3 * expected_density);
    }
    // Draws must spread over the window, not sit on one endpoint.
    EXPECT_GT(distinct.size(), 100u);
}

TEST(KinematicInjectionFailure, ScatteringFailuresAreRetryable) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirectedScatteringChannel geometry_channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    siren::injection::DetectorDirectedScatteringChannel mapping_channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        siren::injection::DetectorDirectedScatteringChannel::Q2Mode::Propagator,
        0.1);

    InteractionRecord record = ScatteringRecord(0.1, 1.0, 1.0, 1.0);
    auto momenta_before = record.secondary_momenta;
    auto random = std::make_shared<siren::utilities::SIREN_random>(173205);

    EXPECT_THROW(geometry_channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
    EXPECT_THROW(mapping_channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
}

TEST(KinematicInjectionFailure, ThreeBodyThresholdFailuresAreRetryable) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirected3BodyChannel direct_channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        PhaseSpaceTopology::Scatter2to3);
    siren::injection::DetectorDirected3BodyChannel recursive_channel(
        target, 0, 1, 2, 1,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        PhaseSpaceTopology::Scatter2to3);

    InteractionRecord record = ThresholdThreeBodyScatteringRecord();
    auto momenta_before = record.secondary_momenta;
    auto random = std::make_shared<siren::utilities::SIREN_random>(223607);

    EXPECT_DOUBLE_EQ(direct_channel.Density(nullptr, record), 0.0);
    EXPECT_THROW(direct_channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);

    EXPECT_DOUBLE_EQ(recursive_channel.Density(nullptr, record), 0.0);
    EXPECT_THROW(recursive_channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
}

TEST(KinematicInjectionFailure, IsotropicTwoBodyRejectsSubThresholdDecay) {
    siren::injection::Isotropic2BodyChannel channel(0);
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_mass = 0.9;
    record.primary_momentum = {0.9, 0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.5};
    auto momenta_before = record.secondary_momenta;
    auto random = std::make_shared<siren::utilities::SIREN_random>(264575);

    EXPECT_DOUBLE_EQ(channel.Density(nullptr, record), 0.0);
    EXPECT_THROW(channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
}

TEST(KinematicInjectionFailure, IsotropicTwoBodyAcceptsExactThreshold) {
    siren::injection::Isotropic2BodyChannel channel(1);
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_mass = 1.0;
    record.primary_momentum = {1.0, 0.0, 0.0, 0.0};
    record.secondary_masses = {0.4, 0.6};
    auto random = std::make_shared<siren::utilities::SIREN_random>(271828);

    EXPECT_NO_THROW(channel.Sample(random, nullptr, record));
    EXPECT_NEAR(channel.Density(nullptr, record), 1.0 / (4.0 * M_PI), 1e-15);
    for (size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(record.secondary_momenta[i][0],
                    record.secondary_masses[i], 1e-15);
        EXPECT_DOUBLE_EQ(record.secondary_momenta[i][1], 0.0);
        EXPECT_DOUBLE_EQ(record.secondary_momenta[i][2], 0.0);
        EXPECT_DOUBLE_EQ(record.secondary_momenta[i][3], 0.0);
    }
}

TEST(KinematicInjectionFailure, DirectedStepsRejectSubThresholdParent) {
    // The same acceptance rule as the isotropic channel: parent_mass below
    // the daughter-mass sum has no two-body split, so the shared directed
    // samplers must fail the attempt instead of propagating NaN kinematics,
    // and the density side must report zero support.
    constexpr double parent_energy = 2.0;
    constexpr double parent_mass = 0.9;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    double parent_momentum = std::sqrt(
        parent_energy * parent_energy - parent_mass * parent_mass);
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        1.0, 0.0).create();
    siren::math::Vector3D vertex(0.0, 0.0, 0.0);
    auto random = std::make_shared<siren::utilities::SIREN_random>(577215);
    double target_volume = 4.0 * M_PI / 3.0;

    EXPECT_THROW(
        siren::injection::detail::SampleDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume,
            random),
        siren::utilities::InjectionFailure);
    EXPECT_DOUBLE_EQ(
        siren::injection::detail::DensityDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            1.0, 0.0, 0.0, 1.0,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume),
        0.0);

    siren::injection::detail::AngularSectorBin bin{
        0.0, 1.0, 0.0, siren::injection::detail::kTwoPi};
    EXPECT_THROW(
        siren::injection::detail::SampleAngularSectorStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, bin, random),
        siren::utilities::InjectionFailure);
    EXPECT_DOUBLE_EQ(
        siren::injection::detail::DensityAngularSectorStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            1.0, 0.0, 0.0, 1.0,
            vertex, *target, bin),
        0.0);
}


TEST(ChannelValidation, TwoBodyConstructorsRejectInvalidDaughterIndices) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    for (int index : {-1, 2}) {
        EXPECT_THROW(
            (void)siren::injection::DetectorDirected2BodyChannel(
                target, index,
                siren::injection::DetectorDirected2BodyChannel::Mode::Cone),
            std::runtime_error);
        EXPECT_THROW(
            (void)siren::injection::DetectorDirectedAngularSectorChannel(
                target, 0.0, 1.0, 0.0, M_PI, index),
            std::runtime_error);
        EXPECT_THROW(
            (void)siren::injection::Isotropic2BodyChannel(index),
            std::runtime_error);
    }
}

TEST(ChannelValidation, MalformedTwoBodyRecordsFailBeforeWriting) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    auto directed =
        std::make_shared<siren::injection::DetectorDirected2BodyChannel>(
            target, 0,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    auto sector =
        std::make_shared<siren::injection::DetectorDirectedAngularSectorChannel>(
            target, 0.0, 1.0, 0.0, M_PI, 0);
    auto isotropic =
        std::make_shared<siren::injection::Isotropic2BodyChannel>(0);
    std::vector<std::shared_ptr<PhaseSpaceChannel>> channels{
        directed, sector, isotropic};
    auto random = std::make_shared<siren::utilities::SIREN_random>(360555);

    for (auto const & channel : channels) {
        InteractionRecord missing_mass = TwoBodyDecayRecord();
        missing_mass.secondary_masses.resize(1);
        auto mass_momenta_before = missing_mass.secondary_momenta;
        EXPECT_THROW(channel->Sample(random, nullptr, missing_mass),
                     std::runtime_error);
        EXPECT_EQ(missing_mass.secondary_momenta, mass_momenta_before);
        EXPECT_DOUBLE_EQ(channel->Density(nullptr, missing_mass), 0.0);

        InteractionRecord missing_momentum = TwoBodyDecayRecord();
        missing_momentum.secondary_momenta.resize(1);
        auto momenta_before = missing_momentum.secondary_momenta;
        EXPECT_THROW(channel->Sample(random, nullptr, missing_momentum),
                     std::runtime_error);
        EXPECT_EQ(missing_momentum.secondary_momenta, momenta_before);
        EXPECT_DOUBLE_EQ(channel->Density(nullptr, missing_momentum), 0.0);
    }

    InteractionRecord malformed = TwoBodyDecayRecord();
    malformed.secondary_masses.resize(1);
    EXPECT_FALSE(directed->DirectingActive(malformed));
    EXPECT_FALSE(sector->DirectingActive(malformed));
}

TEST(ChannelValidation, RecursiveThreeBodyRolesMustBePermutation) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    using Channel = siren::injection::DetectorDirected3BodyChannel;
    using Mode = siren::injection::DetectorDirected2BodyChannel::Mode;

    EXPECT_NO_THROW(
        (void)Channel(
            target, 0, 1, 2, 1,
            Channel::InvariantMassMode::Uniform,
            0.0, 0.0, 0.8, 0.0, Mode::Cone));
    EXPECT_THROW(
        (void)Channel(
            target, 0, 0, 2, 2,
            Channel::InvariantMassMode::Uniform,
            0.0, 0.0, 0.8, 0.0, Mode::Cone),
        std::runtime_error);
    EXPECT_THROW(
        (void)Channel(
            target, 0, 1, 3, 1,
            Channel::InvariantMassMode::Uniform,
            0.0, 0.0, 0.8, 0.0, Mode::Cone),
        std::runtime_error);
    EXPECT_THROW(
        (void)Channel(
            target, 0, 1, 2, 0,
            Channel::InvariantMassMode::Uniform,
            0.0, 0.0, 0.8, 0.0, Mode::Cone),
        std::runtime_error);
}

TEST(ChannelValidation, UnknownInvariantMassModeIsRejected) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    using Channel = siren::injection::DetectorDirected3BodyChannel;
    EXPECT_THROW(
        (void)Channel(
            target, 0, static_cast<Channel::InvariantMassMode>(999),
            0.0, 0.0, 0.8, 0.0,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone),
        std::invalid_argument);
}

TEST(ThreeBodyRoles, UnifiedRolesDefineBothFactorizationMeasures) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    using Channel = siren::injection::DetectorDirected3BodyChannel;
    using Mode = siren::injection::DetectorDirected2BodyChannel::Mode;

    Channel direct(
        target, 2, Channel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0, Mode::Cone);
    Channel recursive(
        target, 2, 0, 1, 1, Channel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0, Mode::Cone);

    for (auto const * channel : {&direct, &recursive}) {
        auto measure = channel->Measure();
        EXPECT_EQ(measure.type, PhaseSpaceMeasure::Type::Recursive2Body);
        EXPECT_EQ(measure.spectator, 2);
        EXPECT_EQ(measure.pair_first, 0);
        EXPECT_EQ(measure.pair_second, 1);
    }
}

TEST(SharedInteractionRecordUtils, ReadWriteAndValidationUseOneLayout) {
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_momentum = {5.0, 1.0, 2.0, 3.0};
    record.interaction_vertex = {4.0, 5.0, 6.0};

    ASSERT_TRUE(siren::injection::detail::HasSecondaryStorage(record, 2));
    auto primary = siren::injection::detail::ReadPrimary(record);
    EXPECT_DOUBLE_EQ(primary.e, 5.0);
    EXPECT_DOUBLE_EQ(primary.p.GetX(), 1.0);
    EXPECT_DOUBLE_EQ(primary.p.GetY(), 2.0);
    EXPECT_DOUBLE_EQ(primary.p.GetZ(), 3.0);
    auto vertex = siren::injection::detail::ReadVertex(record);
    EXPECT_DOUBLE_EQ(vertex.GetX(), 4.0);
    EXPECT_DOUBLE_EQ(vertex.GetY(), 5.0);
    EXPECT_DOUBLE_EQ(vertex.GetZ(), 6.0);

    siren::injection::detail::WriteSecondary(record, 1, {
        7.0, siren::math::Vector3D(8.0, 9.0, 10.0)});
    auto secondary = siren::injection::detail::ReadSecondary(record, 1);
    EXPECT_DOUBLE_EQ(secondary.e, 7.0);
    EXPECT_DOUBLE_EQ(secondary.p.GetX(), 8.0);
    EXPECT_DOUBLE_EQ(secondary.p.GetY(), 9.0);
    EXPECT_DOUBLE_EQ(secondary.p.GetZ(), 10.0);

    record.secondary_momenta.pop_back();
    EXPECT_FALSE(siren::injection::detail::HasSecondaryStorage(record, 2));
    EXPECT_THROW(
        siren::injection::detail::RequireSecondaryStorage(
            record, 2, "test"),
        std::runtime_error);
}

TEST(DirectedRejectionExhaustion, SharedStepRejectsInjectionAttempt) {
    // beta == 1 with a finite supplied parent mass deliberately makes the
    // inverse solver reject every direction. This provides a deterministic
    // exhaustion fixture for the same post-classification numerical-failure
    // path that can occur near a real kinematic boundary.
    constexpr double parent_energy = 10.0;
    constexpr double parent_momentum = 10.0;
    constexpr double parent_mass = 2.0;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        1.0, 0.0).create();
    siren::math::Vector3D vertex(0.0, 0.0, 0.0);
    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        parent_energy, 0.0, 0.0, parent_momentum,
        parent_mass, daughter_mass, other_mass, vertex, *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::BoundInKin);

    auto random = std::make_shared<siren::utilities::SIREN_random>(331663);
    double target_volume = 4.0 * M_PI / 3.0;
    EXPECT_THROW(
        siren::injection::detail::SampleDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
            random),
        siren::utilities::InjectionFailure);
}

TEST(DirectedRejectionExhaustion, AngularSectorRejectsInjectionAttempt) {
    constexpr double parent_energy = 10.0;
    constexpr double parent_momentum = 10.0;
    constexpr double parent_mass = 2.0;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        1.0, 0.0).create();
    siren::math::Vector3D vertex(0.0, 0.0, 0.0);
    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        parent_energy, 0.0, 0.0, parent_momentum,
        parent_mass, daughter_mass, other_mass, vertex, *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::BoundInKin);

    siren::injection::detail::AngularSectorBin bin{
        0.0, 1.0, 0.0, siren::injection::detail::kTwoPi};
    auto random = std::make_shared<siren::utilities::SIREN_random>(346410);
    EXPECT_THROW(
        siren::injection::detail::SampleAngularSectorStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, bin, random),
        siren::utilities::InjectionFailure);
}

TEST(TabulatedMappingSupport, ZeroOverlapWindowHasZeroDensity) {
    siren::injection::TabulatedMapping map(
        {1e-4, 1.0}, {0.0, 1.0}, 0.0, 5e-5);

    EXPECT_FALSE(map.HasSupport());
    EXPECT_DOUBLE_EQ(map.Density(2.5e-5), 0.0);
    EXPECT_DOUBLE_EQ(map.Inverse(2.5e-5), 0.0);
    EXPECT_THROW(map.Forward(0.5), std::runtime_error);
}

TEST(TabulatedMappingSupport, ImmutableTableIsSharedAcrossEventWindows) {
    std::shared_ptr<siren::injection::TabulatedMappingTable const> table =
        std::make_shared<siren::injection::TabulatedMappingTable>(
            std::vector<double>{0.0, 1.0, 2.0, 3.0},
            std::vector<double>{0.0, 0.2, 0.7, 1.0});
    ASSERT_EQ(table.use_count(), 1);

    siren::injection::TabulatedMapping low_window(table, 0.0, 2.0);
    siren::injection::TabulatedMapping high_window(table, 1.0, 3.0);
    EXPECT_EQ(table.use_count(), 3);
    EXPECT_EQ(low_window.table.get(), table.get());
    EXPECT_EQ(high_window.table.get(), table.get());
    EXPECT_GT(low_window.Density(0.5), 0.0);
    EXPECT_GT(high_window.Density(2.5), 0.0);
    EXPECT_DOUBLE_EQ(low_window.Density(2.5), 0.0);
    EXPECT_DOUBLE_EQ(high_window.Density(0.5), 0.0);
}

TEST(TabulatedMappingSupport, RejectsUnsortedAndDuplicateNodes) {
    EXPECT_THROW(
        (void)siren::injection::TabulatedMapping(
            {0.0, 2.0, 1.0}, {0.0, 0.5, 1.0}, 0.0, 2.0),
        std::runtime_error);
    EXPECT_THROW(
        (void)siren::injection::TabulatedMapping(
            {0.0, 1.0, 1.0}, {0.0, 0.5, 1.0}, 0.0, 1.0),
        std::runtime_error);
    EXPECT_THROW(
        (void)siren::injection::TabulatedMapping(
            {0.0, std::numeric_limits<double>::quiet_NaN(), 2.0},
            {0.0, 0.5, 1.0}, 0.0, 2.0),
        std::runtime_error);
}

TEST(TabulatedMappingSupport, ValidatesCdfWithoutRejectingPlateaus) {
    EXPECT_THROW(
        (void)siren::injection::TabulatedMapping(
            {0.0, 1.0, 2.0}, {0.0, 1.0, 0.5}, 0.0, 2.0),
        std::runtime_error);
    EXPECT_THROW(
        (void)siren::injection::TabulatedMapping(
            {0.0, 1.0, 2.0},
            {0.0, std::numeric_limits<double>::infinity(), 1.0},
            0.0, 2.0),
        std::runtime_error);

    siren::injection::TabulatedMapping plateau(
        {0.0, 1.0, 2.0}, {0.0, 0.5, 0.5}, 0.0, 2.0);
    EXPECT_TRUE(plateau.HasSupport());
    EXPECT_DOUBLE_EQ(plateau.Density(1.5), 0.0);
}

TEST(TabulatedMappingSupport, ThreeBodyDensityIsEvaluableOffTable) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirected3BodyChannel uniform_channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    siren::injection::DetectorDirected3BodyChannel tabulated_channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Tabulated,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        PhaseSpaceTopology::Decay3Body,
        {10.0, 20.0}, {0.0, 1.0});

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.2, 0.2};
    record.secondary_momenta.resize(3);
    auto random = std::make_shared<siren::utilities::SIREN_random>(264575);
    uniform_channel.Sample(random, nullptr, record);

    double density = -1.0;
    EXPECT_NO_THROW(density = tabulated_channel.Density(nullptr, record));
    EXPECT_DOUBLE_EQ(density, 0.0);

    auto momenta_before = record.secondary_momenta;
    EXPECT_THROW(tabulated_channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
}

TEST(TabulatedMappingSupport, ScatteringDensityIsEvaluableOffTable) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirectedScatteringChannel geometry_channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    siren::injection::DetectorDirectedScatteringChannel tabulated_channel(
        target, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
        siren::injection::DetectorDirectedScatteringChannel::Q2Mode::Tabulated,
        0.0, {100.0, 200.0}, {0.0, 1.0});

    InteractionRecord record = ScatteringRecord(10.0, 1.0, 0.1, 1.0);
    auto random = std::make_shared<siren::utilities::SIREN_random>(282843);
    geometry_channel.Sample(random, nullptr, record);

    double density = -1.0;
    EXPECT_NO_THROW(density = tabulated_channel.Density(nullptr, record));
    EXPECT_DOUBLE_EQ(density, 0.0);

    auto momenta_before = record.secondary_momenta;
    EXPECT_THROW(tabulated_channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
}

TEST(AngularSectorPhi, ConstructorCanonicalizesOneTurnRanges) {
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        5.0, 0.0).create();

    EXPECT_NO_THROW(siren::injection::DetectorDirectedAngularSectorChannel(
        target, 0.0, 1.0, -M_PI / 4.0, M_PI / 4.0, 0));
    EXPECT_NO_THROW(siren::injection::DetectorDirectedAngularSectorChannel(
        target, 0.0, 1.0, -M_PI, M_PI, 0));
    EXPECT_THROW(siren::injection::DetectorDirectedAngularSectorChannel(
        target, 0.0, 1.0, 0.0, 2.0 * M_PI + 1e-6, 0),
        std::runtime_error);
    EXPECT_THROW(siren::injection::DetectorDirectedAngularSectorChannel(
        target, 0.0, 1.0,
        std::numeric_limits<double>::quiet_NaN(), M_PI, 0),
        std::runtime_error);
}

TEST(AngularSectorPhi, CenteredWrappedRangeSampleDensityClosure) {
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        5.0, 0.0).create();
    siren::injection::DetectorDirectedAngularSectorChannel channel(
        target, 0.0, 1.0, -M_PI / 4.0, M_PI / 4.0, 0);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.5};
    record.secondary_momenta.resize(2);

    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        2.0, 0.0, 0.0, 0.0, 2.0, 0.5, 0.5,
        siren::math::Vector3D(0.0, 0.0, 0.0), *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Rest);
    ASSERT_FALSE(geometry.inside_geometry);
    double expected_density = 1.0 /
        ((1.0 - std::cos(geometry.theta_bound)) * (M_PI / 2.0));

    siren::math::Vector3D perp1, perp2;
    siren::injection::detail::SectorPerpFrame(
        geometry.to_center, perp1, perp2);
    auto random = std::make_shared<siren::utilities::SIREN_random>(8675309);
    int negative_phi = 0;
    int positive_phi = 0;

    for (int i = 0; i < 1000; ++i) {
        channel.Sample(random, nullptr, record);
        siren::math::Vector3D direction(
            record.secondary_momenta[0][1],
            record.secondary_momenta[0][2],
            record.secondary_momenta[0][3]);
        direction.normalize();
        double cos_to = siren::math::scalar_product(
            direction, geometry.to_center);
        siren::math::Vector3D perpendicular =
            direction - geometry.to_center * cos_to;
        double phi = std::atan2(
            siren::math::scalar_product(perpendicular, perp2),
            siren::math::scalar_product(perpendicular, perp1));
        if (phi < 0.0) ++negative_phi;
        if (phi > 0.0) ++positive_phi;

        EXPECT_NEAR(channel.Density(nullptr, record), expected_density,
                    1e-12 * expected_density);
    }

    EXPECT_GT(negative_phi, 400);
    EXPECT_GT(positive_phi, 400);
}

TEST(AngularSectorRest, AabbCenterProducesOnShellBackToBackMomenta) {
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 0.0)),
        2.0, 1.0).create();
    siren::math::Vector3D vertex(0.0, 0.0, 0.0);
    ASSERT_FALSE(target->IsInside(vertex));

    siren::injection::DetectorDirectedAngularSectorChannel channel(
        target, 0.0, 1.0, 0.0, 2.0 * M_PI, 0);
    InteractionRecord record = TwoBodyDecayRecord();

    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        record.primary_momentum[0],
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3],
        record.primary_mass,
        record.secondary_masses[0],
        record.secondary_masses[1],
        vertex, *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Rest);
    ASSERT_FALSE(geometry.inside_geometry);
    EXPECT_DOUBLE_EQ(geometry.theta_bound, M_PI);
    EXPECT_NEAR(geometry.to_center.magnitude(), 1.0, 1e-15);
    ASSERT_TRUE(channel.DirectingActive(record));

    const double expected_p = siren::injection::TwoBodyRestMomentum(
        record.primary_mass,
        record.secondary_masses[0],
        record.secondary_masses[1]);
    auto random = std::make_shared<siren::utilities::SIREN_random>(141421);
    for (int sample_index = 0; sample_index < 256; ++sample_index) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        for (int daughter = 0; daughter < 2; ++daughter) {
            auto const & momentum = record.secondary_momenta[daughter];
            double p2 = momentum[1] * momentum[1]
                      + momentum[2] * momentum[2]
                      + momentum[3] * momentum[3];
            EXPECT_NEAR(std::sqrt(p2), expected_p, 1e-14);
            EXPECT_NEAR(momentum[0] * momentum[0] - p2,
                        record.secondary_masses[daughter]
                            * record.secondary_masses[daughter],
                        1e-14);
        }
        for (int component = 0; component < 4; ++component) {
            EXPECT_NEAR(
                record.secondary_momenta[0][component]
                    + record.secondary_momenta[1][component],
                record.primary_momentum[component], 1e-14);
        }
        EXPECT_NEAR(channel.Density(nullptr, record),
                    1.0 / siren::injection::detail::kFourPi, 1e-15);
    }
}

TEST(AngularSectorCaching, SiblingMixtureClassifiesDirectedRegimeOnce) {
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        5.0, 0.0).create();
    MultiChannelPhaseSpace mixture;
    for (int sector = 0; sector < 4; ++sector) {
        mixture.channels.push_back(
            std::make_shared<
                siren::injection::DetectorDirectedAngularSectorChannel>(
                    target, 0.0, 1.0,
                    0.5 * M_PI * sector,
                    0.5 * M_PI * (sector + 1), 0));
    }
    mixture.weights.assign(4, 0.25);

    InteractionRecord record = TwoBodyDecayRecord();
    auto random = std::make_shared<siren::utilities::SIREN_random>(17320508);
    mixture.channels.front()->Sample(random, nullptr, record);

    siren::injection::detail::ResetDirectedRegimeCacheForTesting();
    mixture.Accumulate(nullptr, record, 1.0, true, false);
    EXPECT_EQ(
        siren::injection::detail::DirectedRegimeCacheMissesForTesting(), 1u);
    EXPECT_EQ(mixture.kp_count_, 1);
}

TEST(DirectedGeometryVolume, UnresolvedThinCompositeDoesNotBecomeFullAabb) {
    constexpr double outer_width = 10.0;
    constexpr double inner_width = 9.99999;
    auto outer = siren::geometry::Box(
        outer_width, outer_width, outer_width).create();
    auto inner = siren::geometry::Box(
        inner_width, inner_width, inner_width).create();
    auto shell = siren::geometry::BooleanGeometry(
        siren::geometry::BooleanOperation::SUBTRACTION,
        outer, inner).create();
    EXPECT_TRUE(std::isnan(siren::injection::ExactGeometryVolume(*shell)));
    EXPECT_THROW(
        siren::injection::DetectorDirected2BodyChannel(
            shell, 0,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume),
        std::runtime_error);
    EXPECT_THROW(
        siren::injection::DetectorDirected3BodyChannel(
            shell, 0,
            siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
            0.0, 0.0, 0.8, 0.0,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume),
        std::runtime_error);
    EXPECT_THROW(
        siren::injection::DetectorDirected3BodyChannel(
            shell, 0, 1, 2, 1,
            siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
            0.0, 0.0, 0.8, 0.0,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume),
        std::runtime_error);
    EXPECT_THROW(
        siren::injection::DetectorDirectedScatteringChannel(
            shell, 0,
            siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume),
        std::runtime_error);
}

TEST(DirectedGeometryVolume, AnalyticVolumesIncludeAngularCuts) {
    constexpr double cylinder_outer = 4.0;
    constexpr double cylinder_inner = 1.0;
    constexpr double cylinder_height = 5.0;
    constexpr double cylinder_delta_phi = M_PI / 3.0;
    siren::geometry::Cylinder cylinder(
        cylinder_outer, cylinder_inner, cylinder_height,
        0.2, cylinder_delta_phi);
    double expected_cylinder = 0.5 * cylinder_delta_phi
        * (cylinder_outer * cylinder_outer - cylinder_inner * cylinder_inner)
        * cylinder_height;
    EXPECT_NEAR(siren::injection::ExactGeometryVolume(cylinder),
                expected_cylinder, 1e-14 * expected_cylinder);

    constexpr double sphere_outer = 4.0;
    constexpr double sphere_inner = 1.0;
    constexpr double sphere_delta_phi = M_PI / 2.0;
    constexpr double sphere_start_theta = 0.4;
    constexpr double sphere_delta_theta = 0.7;
    siren::geometry::Sphere sphere(
        sphere_outer, sphere_inner,
        0.3, sphere_delta_phi,
        sphere_start_theta, sphere_delta_theta);
    double expected_sphere =
        (sphere_outer * sphere_outer * sphere_outer
         - sphere_inner * sphere_inner * sphere_inner) / 3.0
        * sphere_delta_phi
        * (std::cos(sphere_start_theta)
           - std::cos(sphere_start_theta + sphere_delta_theta));
    EXPECT_NEAR(siren::injection::ExactGeometryVolume(sphere),
                expected_sphere, 1e-14 * expected_sphere);
}

TEST(DirectedGeometryVolume, ExactCallerVolumeBypassesUnresolvedEstimate) {
    constexpr double outer_width = 10.0;
    constexpr double inner_width = 9.99999;
    auto shell = siren::geometry::BooleanGeometry(
        siren::geometry::BooleanOperation::SUBTRACTION,
        siren::geometry::Box(
            outer_width, outer_width, outer_width).create(),
        siren::geometry::Box(
            inner_width, inner_width, inner_width).create()).create();
    double exact_volume = outer_width * outer_width * outer_width
                        - inner_width * inner_width * inner_width;

    EXPECT_NO_THROW(siren::injection::DetectorDirected2BodyChannel(
        shell, 0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Volume,
        exact_volume));
    EXPECT_NO_THROW(siren::injection::DetectorDirected3BodyChannel(
        shell, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology::Decay3Body, {}, {}, exact_volume));
    EXPECT_NO_THROW(siren::injection::DetectorDirected3BodyChannel(
        shell, 0, 1, 2, 1,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology::Decay3Body, {}, {}, exact_volume));
    EXPECT_NO_THROW(siren::injection::DetectorDirectedScatteringChannel(
        shell, 0,
        siren::injection::DetectorDirectedScatteringChannel::Variable::Q2,
        siren::injection::DetectorDirected2BodyChannel::Mode::Volume,
        siren::injection::DetectorDirectedScatteringChannel::Q2Mode::Geometry,
        0.0, {}, {}, exact_volume));
}

TEST(DirectedOverlapSampling, CoversLensExcludedByLegacyCap) {
    constexpr double beta = 0.9;
    constexpr double theta_kin = 0.4;
    constexpr double theta_bound = 1.15;
    constexpr double axis_sep = 0.854;
    constexpr double daughter_mass = 1.0;
    constexpr double other_mass = 0.0;

    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    double cos_critical = std::cos(theta_kin);
    double rest_energy = daughter_mass * gamma
        * std::sqrt(1.0 - beta * beta * cos_critical * cos_critical);
    double rest_momentum = std::sqrt(
        rest_energy * rest_energy - daughter_mass * daughter_mass);
    double parent_mass = rest_energy + rest_momentum;
    double parent_energy = gamma * parent_mass;
    double parent_momentum = beta * parent_energy;

    siren::math::Vector3D parent_axis(0.0, 0.0, 1.0);
    siren::math::Vector3D bound_axis(
        std::sin(axis_sep), 0.0, std::cos(axis_sep));
    constexpr double center_distance = 10.0;
    double bounding_radius = center_distance * std::sin(theta_bound);
    double sphere_radius = bounding_radius / std::sqrt(3.0);
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(bound_axis * center_distance),
        sphere_radius, 0.0).create();
    double target_volume = (4.0 / 3.0) * M_PI
        * sphere_radius * sphere_radius * sphere_radius;
    siren::math::Vector3D origin(0.0, 0.0, 0.0);

    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        parent_energy, 0.0, 0.0, parent_momentum,
        parent_mass, daughter_mass, other_mass, origin, *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Overlap);
    ASSERT_NEAR(geometry.theta_kin, theta_kin, 1e-12);
    ASSERT_NEAR(geometry.theta_bound, theta_bound, 1e-12);

    auto [legacy_center, legacy_cos_half] = LegacyOverlapCap(
        geometry.theta_kin, geometry.theta_bound, geometry.axis_sep,
        parent_axis, geometry.to_center);
    auto random = std::make_shared<siren::utilities::SIREN_random>(314159);
    int outside_legacy_cap = 0;

    for (int i = 0; i < 1000; ++i) {
        auto sample = siren::injection::detail::SampleDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            origin, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
            random);

        EXPECT_GE(siren::math::scalar_product(sample.lab_dir, parent_axis),
                  std::cos(geometry.theta_kin) - 1e-12);
        EXPECT_GE(siren::math::scalar_product(sample.lab_dir, geometry.to_center),
                  std::cos(geometry.theta_bound) - 1e-12);
        if (siren::math::scalar_product(sample.lab_dir, legacy_center)
            < legacy_cos_half - 1e-12) {
            ++outside_legacy_cap;
        }

        double density = siren::injection::detail::DensityDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            sample.E_lab,
            sample.p_lab * sample.lab_dir.GetX(),
            sample.p_lab * sample.lab_dir.GetY(),
            sample.p_lab * sample.lab_dir.GetZ(),
            origin, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
        EXPECT_NEAR(density, sample.rest_density,
                    1e-10 * std::max(1.0, sample.rest_density));
    }

    // The reference cap misses about 27% of this lens, so a threshold well
    // below that number makes this a deterministic assertion.
    EXPECT_GT(outside_legacy_cap, 100);
}

TEST(DirectedRestSampling, BoundingSphereInteriorUsesFullSphereDensity) {
    auto target = siren::geometry::Box(2.0, 2.0, 2.0).create();
    siren::math::Vector3D vertex(1.2, 0.0, 0.0);
    constexpr double parent_mass = 2.0;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    constexpr double target_volume = 8.0;

    ASSERT_FALSE(target->IsInside(vertex));
    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        parent_mass, 0.0, 0.0, 0.0,
        parent_mass, daughter_mass, other_mass, vertex, *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Rest);
    ASSERT_FALSE(geometry.inside_geometry);
    EXPECT_DOUBLE_EQ(geometry.theta_bound, M_PI);
    EXPECT_DOUBLE_EQ(geometry.omega_bound,
                     siren::injection::detail::kFourPi);

    auto random = std::make_shared<siren::utilities::SIREN_random>(271828);
    int away_from_detector = 0;
    for (int i = 0; i < 1000; ++i) {
        auto sample = siren::injection::detail::SampleDirectedStep(
            parent_mass, 0.0, 0.0, 0.0,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
            random);

        if (siren::math::scalar_product(sample.lab_dir, geometry.to_center) < 0.0) {
            ++away_from_detector;
        }
        EXPECT_NEAR(sample.rest_density,
                    1.0 / siren::injection::detail::kFourPi, 1e-15);

        double density = siren::injection::detail::DensityDirectedStep(
            parent_mass, 0.0, 0.0, 0.0,
            parent_mass, daughter_mass, other_mass,
            sample.E_lab,
            sample.p_lab * sample.lab_dir.GetX(),
            sample.p_lab * sample.lab_dir.GetY(),
            sample.p_lab * sample.lab_dir.GetZ(),
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
        EXPECT_NEAR(density, sample.rest_density, 1e-15);
    }

    // A full-sphere sampler should populate both hemispheres.
    EXPECT_GT(away_from_detector, 400);
    EXPECT_LT(away_from_detector, 600);
}

TEST(DirectedRestSampling, OutsideTargetConeUsesOnlyInitializedRestBranch) {
    siren::math::Vector3D center(0.0, 0.0, 10.0);
    constexpr double radius = 1.0;
    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(center), radius, 0.0).create();
    siren::math::Vector3D vertex(0.0, 0.0, 0.0);
    constexpr double parent_mass = 2.0;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    constexpr double target_volume = (4.0 / 3.0) * M_PI
        * radius * radius * radius;

    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        parent_mass, 0.0, 0.0, 0.0,
        parent_mass, daughter_mass, other_mass, vertex, *target);
    ASSERT_EQ(geometry.regime,
              siren::injection::detail::DirectedRegime::Rest);
    ASSERT_FALSE(geometry.inside_geometry);
    ASSERT_GT(geometry.omega_bound, 0.0);
    ASSERT_LT(geometry.omega_bound,
              siren::injection::detail::kFourPi);
    const double expected_density = 1.0 / geometry.omega_bound;

    auto random = std::make_shared<siren::utilities::SIREN_random>(161803);
    for (int i = 0; i < 256; ++i) {
        auto sample = siren::injection::detail::SampleDirectedStep(
            parent_mass, 0.0, 0.0, 0.0,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
            random);

        EXPECT_TRUE(std::isfinite(sample.rest_density));
        EXPECT_DOUBLE_EQ(sample.rest_density, expected_density);
        EXPECT_NEAR(sample.p_lab,
                    siren::injection::TwoBodyRestMomentum(
                        parent_mass, daughter_mass, other_mass),
                    1e-15);

        double density = siren::injection::detail::DensityDirectedStep(
            parent_mass, 0.0, 0.0, 0.0,
            parent_mass, daughter_mass, other_mass,
            sample.E_lab,
            sample.p_lab * sample.lab_dir.GetX(),
            sample.p_lab * sample.lab_dir.GetY(),
            sample.p_lab * sample.lab_dir.GetZ(),
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
        EXPECT_DOUBLE_EQ(density, expected_density);
    }
}

TEST(DirectedThresholdKinematics, CriticalCosineClampPreventsNaNRecord) {
    // This mass combination reproducibly exposes a cancellation: E_rest
    // rounds just below m_daughter, so the unbounded critical-cosine formula
    // evaluates slightly above one at threshold.
    constexpr double daughter_mass = 0.001;
    constexpr double other_mass = 300.0;
    constexpr double parent_mass = daughter_mass + other_mass;
    constexpr double gamma = 2.0;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    const double parent_energy = gamma * parent_mass;
    const double parent_momentum = beta * parent_energy;
    const double p_rest = siren::injection::TwoBodyRestMomentum(
        parent_mass, daughter_mass, other_mass);
    const double E_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_mass, other_mass);

    double cos_critical = siren::injection::CriticalCosTheta(
        beta, gamma, p_rest, E_rest, daughter_mass);
    ASSERT_TRUE(std::isfinite(cos_critical));
    EXPECT_GE(cos_critical, -1.0);
    EXPECT_LE(cos_critical, 1.0);

    auto target = siren::geometry::Sphere(
        siren::geometry::Placement(
            siren::math::Vector3D(0.0, 0.0, 100.0)),
        1.0, 0.0).create();
    siren::math::Vector3D vertex(0.0, 0.0, 0.0);
    auto geometry = siren::injection::detail::ClassifyDirectedRegime(
        parent_energy, 0.0, 0.0, parent_momentum,
        parent_mass, daughter_mass, other_mass, vertex, *target);
    EXPECT_TRUE(std::isfinite(geometry.theta_kin));
    EXPECT_TRUE(std::isfinite(geometry.omega_kin));
    EXPECT_GE(geometry.theta_kin, 0.0);
    EXPECT_LE(geometry.theta_kin, M_PI);

    siren::injection::DetectorDirected2BodyChannel channel(
        target, 0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = parent_mass;
    record.primary_momentum = {
        parent_energy, 0.0, 0.0, parent_momentum};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {daughter_mass, other_mass};
    record.secondary_momenta.resize(2);

    auto random = std::make_shared<siren::utilities::SIREN_random>(57721);
    for (int i = 0; i < 100; ++i) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        for (auto const & momentum : record.secondary_momenta) {
            for (double component : momentum) {
                EXPECT_TRUE(std::isfinite(component));
            }
        }
        double density = channel.Density(nullptr, record);
        EXPECT_TRUE(std::isfinite(density));
        EXPECT_GT(density, 0.0);
    }
}

TEST(TwoBodyLabAngle, RejectsNaNAndClampsRoundoffAtAngularBoundary) {
    constexpr double parent_mass = 2.0;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    constexpr double beta = 0.5;
    const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    const double p_rest = siren::injection::TwoBodyRestMomentum(
        parent_mass, daughter_mass, other_mass);
    const double E_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_mass, other_mass);

    auto invalid = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, E_rest, daughter_mass,
        std::numeric_limits<double>::quiet_NaN());
    for (auto const & solution : invalid) {
        EXPECT_FALSE(solution.valid);
        EXPECT_TRUE(std::isfinite(solution.cos_theta_rest));
        EXPECT_TRUE(std::isfinite(solution.p_lab));
        EXPECT_TRUE(std::isfinite(solution.jacobian));
    }

    auto rounded = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, E_rest, daughter_mass,
        std::nextafter(1.0, 2.0));
    EXPECT_TRUE(rounded[0].valid || rounded[1].valid);
    for (auto const & solution : rounded) {
        if (!solution.valid) continue;
        EXPECT_TRUE(std::isfinite(solution.cos_theta_rest));
        EXPECT_TRUE(std::isfinite(solution.p_lab));
        EXPECT_TRUE(std::isfinite(solution.jacobian));
    }
}

TEST(SharedKinematics, StableBreakupMomentumBacksInjectionWrapper) {
    constexpr double mass_a = 0.001;
    constexpr double mass_b = 300.0;
    const double parent_mass = mass_a + mass_b + 1e-10;
    double canonical = siren::math::TwoBodyRestMomentum(
        parent_mass, mass_a, mass_b);

    EXPECT_TRUE(std::isfinite(canonical));
    EXPECT_GT(canonical, 0.0);
    EXPECT_DOUBLE_EQ(
        siren::injection::TwoBodyRestMomentum(
            parent_mass, mass_a, mass_b),
        canonical);
    EXPECT_DOUBLE_EQ(
        siren::injection::Kallen(7.0, 2.0, 1.0),
        siren::math::Kallen(7.0, 2.0, 1.0));
    EXPECT_DOUBLE_EQ(
        siren::math::TwoBodyRestMomentum(
            mass_a + mass_b, mass_a, mass_b),
        0.0);
}

TEST(SharedInterpolation, BracketUsesRightBinAtInteriorNodes) {
    std::vector<double> grid{0.0, 1.0, 3.0, 10.0};
    EXPECT_EQ(
        siren::math::InterpolationBracket(grid, -1.0),
        (std::pair<std::size_t, std::size_t>{0, 1}));
    EXPECT_EQ(
        siren::math::InterpolationBracket(grid, 1.0),
        (std::pair<std::size_t, std::size_t>{1, 2}));
    EXPECT_EQ(
        siren::math::InterpolationBracket(grid, 2.0),
        (std::pair<std::size_t, std::size_t>{1, 2}));
    EXPECT_EQ(
        siren::math::InterpolationBracket(grid, 20.0),
        (std::pair<std::size_t, std::size_t>{2, 3}));
}

TEST(SharedLorentzBoost, RkAdapterMatchesAnalyticBoost) {
    constexpr double parent_mass = 2.0;
    constexpr double gamma = 1.5;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    const double parent_energy = gamma * parent_mass;
    const double parent_pz = gamma * beta * parent_mass;
    constexpr double daughter_mass = 0.4;
    constexpr double rest_px = 0.3;
    constexpr double rest_py = -0.2;
    constexpr double rest_pz = 0.7;
    const double rest_energy = std::sqrt(
        daughter_mass * daughter_mass
        + rest_px * rest_px + rest_py * rest_py + rest_pz * rest_pz);

    auto lab = siren::injection::detail::BoostRestFrameToLab(
        parent_energy, 0.0, 0.0, parent_pz,
        rest_energy, rest_px, rest_py, rest_pz);
    EXPECT_NEAR(lab[0], gamma * (rest_energy + beta * rest_pz), 1e-14);
    EXPECT_NEAR(lab[1], rest_px, 1e-14);
    EXPECT_NEAR(lab[2], rest_py, 1e-14);
    EXPECT_NEAR(lab[3], gamma * (rest_pz + beta * rest_energy), 1e-14);
    EXPECT_NEAR(
        lab[0] * lab[0]
        - lab[1] * lab[1] - lab[2] * lab[2] - lab[3] * lab[3],
        daughter_mass * daughter_mass, 1e-13);
}

TEST(SharedPerpendicularFrame, IsOrthonormalForAllAxisBranches) {
    for (siren::math::Vector3D axis : {
             siren::math::Vector3D(0.0, 0.0, 1.0),
             siren::math::Vector3D(1.0, 0.0, 0.0),
             siren::math::Vector3D(0.3, -0.4, 0.5).normalized()}) {
        siren::math::Vector3D first, second;
        siren::injection::detail::SectorPerpFrame(axis, first, second);
        EXPECT_NEAR(first.magnitude(), 1.0, 1e-15);
        EXPECT_NEAR(second.magnitude(), 1.0, 1e-15);
        EXPECT_NEAR(siren::math::scalar_product(axis, first), 0.0, 1e-15);
        EXPECT_NEAR(siren::math::scalar_product(axis, second), 0.0, 1e-15);
        EXPECT_NEAR(siren::math::scalar_product(first, second), 0.0, 1e-15);
    }
}

TEST(ThreeBodyPowerLaw, AcceptsMasslessPairAtOffset) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirected3BodyChannel channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.0, 0.0};
    record.secondary_momenta.resize(3);
    auto random = std::make_shared<siren::utilities::SIREN_random>(161803);

    EXPECT_NO_THROW(channel.Sample(random, nullptr, record));
    ASSERT_GT(record.interaction_parameters.at("phase_space_s_X"), 0.0);
    for (auto const & momentum : record.secondary_momenta) {
        for (double component : momentum) EXPECT_TRUE(std::isfinite(component));
    }
    double density = channel.Density(nullptr, record);
    EXPECT_TRUE(std::isfinite(density));
    EXPECT_GT(density, 0.0);
}

TEST(PowerLawMapping, NuOneUsesNormalizedLogarithmicForm) {
    constexpr double offset = 2.0;
    constexpr double x_min = 1.0;
    constexpr double x_max = 100.0;
    siren::injection::PowerLawMapping map(
        1.0, offset, offset + x_min, offset + x_max);

    for (double r : {0.0, 0.1, 0.5, 0.9, 1.0}) {
        double s = map.Forward(r);
        EXPECT_TRUE(std::isfinite(s));
        EXPECT_NEAR(map.Inverse(s), r, 2e-14);
        EXPECT_TRUE(std::isfinite(map.Density(s)));
        EXPECT_GT(map.Density(s), 0.0);
    }
    EXPECT_NEAR(map.Forward(0.5), offset + 10.0, 2e-14);
    EXPECT_NEAR(map.Density(offset + 10.0),
                1.0 / (10.0 * std::log(x_max / x_min)), 1e-15);
    EXPECT_DOUBLE_EQ(map.Density(offset + 0.5), 0.0);
    EXPECT_DOUBLE_EQ(map.Density(offset + 101.0), 0.0);

    // Integrate in logarithmically spaced bins, matching the logarithmic
    // shape of the nu=1 branch's density.
    double integral = 0.0;
    constexpr int bins = 1000;
    double log_range = std::log(x_max / x_min);
    for (int i = 0; i < bins; ++i) {
        double x_lo = x_min * std::exp(log_range * i / bins);
        double x_hi = x_min * std::exp(log_range * (i + 1) / bins);
        double x_mid = std::sqrt(x_lo * x_hi);
        integral += map.Density(offset + x_mid) * (x_hi - x_lo);
    }
    EXPECT_NEAR(integral, 1.0, 1e-6);
}

TEST(PowerLawMapping, RejectsInvalidAndNonNormalizableRanges) {
    EXPECT_THROW(
        (void)siren::injection::PowerLawMapping(0.8, 0.1, 0.0, 1.0),
        std::runtime_error);
    EXPECT_THROW(
        (void)siren::injection::PowerLawMapping(1.0, 0.0, 0.0, 1.0),
        std::runtime_error);
    EXPECT_THROW(
        (void)siren::injection::PowerLawMapping(1.2, 0.0, 0.0, 1.0),
        std::runtime_error);
}

TEST(AnalyticMappingSupport, DensityIsZeroOutsideDeclaredInterval) {
    auto expect_bounded_support = [](
        siren::injection::Mapping1D const & mapping,
        double lower,
        double upper)
    {
        EXPECT_GT(mapping.Density(lower), 0.0);
        EXPECT_GT(mapping.Density(upper), 0.0);
        EXPECT_DOUBLE_EQ(
            mapping.Density(std::nextafter(
                lower, -std::numeric_limits<double>::infinity())),
            0.0);
        EXPECT_DOUBLE_EQ(
            mapping.Density(std::nextafter(
                upper, std::numeric_limits<double>::infinity())),
            0.0);
        EXPECT_DOUBLE_EQ(
            mapping.Density(std::numeric_limits<double>::quiet_NaN()),
            0.0);
    };

    siren::injection::BreitWignerMapping breit_wigner(2.0, 0.2, 1.0, 9.0);
    siren::injection::PowerLawMapping power_law(0.8, 0.0, 1.0, 4.0);
    siren::injection::PropagatorMapping propagator(0.5, 0.0, 4.0);
    siren::injection::UniformMapping uniform(0.0, 4.0);
    siren::injection::LogMapping logarithmic(1.0, 4.0);
    siren::injection::ExponentialMapping exponential(2.0, 0.0, 4.0);
    siren::injection::GaussianMapping gaussian(0.0, 1.0, -2.0, 2.0);

    expect_bounded_support(breit_wigner, 1.0, 9.0);
    expect_bounded_support(power_law, 1.0, 4.0);
    expect_bounded_support(propagator, 0.0, 4.0);
    expect_bounded_support(uniform, 0.0, 4.0);
    expect_bounded_support(logarithmic, 1.0, 4.0);
    expect_bounded_support(exponential, 0.0, 4.0);
    expect_bounded_support(gaussian, -2.0, 2.0);
}

TEST(ThreeBodyPowerLaw, NuOneSamplesAndReportsFiniteDensity) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirected3BodyChannel channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw,
        0.0, 0.0, 1.0, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.2, 0.3};
    record.secondary_momenta.resize(3);
    auto random = std::make_shared<siren::utilities::SIREN_random>(316228);

    for (int i = 0; i < 100; ++i) {
        ASSERT_NO_THROW(channel.Sample(random, nullptr, record));
        double s_X = record.interaction_parameters.at("phase_space_s_X");
        EXPECT_GT(s_X, 0.0);
        EXPECT_TRUE(std::isfinite(s_X));
        for (auto const & momentum : record.secondary_momenta) {
            for (double component : momentum) {
                EXPECT_TRUE(std::isfinite(component));
            }
        }
        double density = channel.Density(nullptr, record);
        EXPECT_TRUE(std::isfinite(density));
        EXPECT_GT(density, 0.0);
    }
}

TEST(ThreeBodyPowerLaw, RejectsRangeBelowOffsetBeforeWritingMomenta) {
    auto target = siren::geometry::Box(10.0, 10.0, 10.0).create();
    siren::injection::DetectorDirected3BodyChannel valid_channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw,
        0.0, 0.0, 0.8, 0.0,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);
    siren::injection::DetectorDirected3BodyChannel invalid_channel(
        target, 0,
        siren::injection::DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw,
        0.0, 0.0, 0.8, 0.1,
        siren::injection::DetectorDirected2BodyChannel::Mode::Cone);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.0, 0.0};
    record.secondary_momenta.resize(3);
    auto random = std::make_shared<siren::utilities::SIREN_random>(141421);
    valid_channel.Sample(random, nullptr, record);
    auto momenta_before = record.secondary_momenta;

    EXPECT_DOUBLE_EQ(invalid_channel.Density(nullptr, record), 0.0);
    EXPECT_THROW(invalid_channel.Sample(random, nullptr, record),
                 std::runtime_error);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
    for (auto const & momentum : record.secondary_momenta) {
        for (double component : momentum) EXPECT_TRUE(std::isfinite(component));
    }
}

} // namespace
