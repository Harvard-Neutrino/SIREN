#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/ConstantDensityDistribution.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/injection/WeightingUtils.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"

#include "../DetectorDirectedChannelUtils.h"

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

class CompletePrimaryDistribution final
    : public siren::distributions::VertexPositionDistribution {
private:
    std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        siren::dataclasses::PrimaryDistributionRecord &) const override
    {
        return {siren::math::Vector3D(0.0, 0.0, 0.0),
                siren::math::Vector3D(0.0, 0.0, 0.0)};
    }

public:
    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        siren::dataclasses::PrimaryDistributionRecord & record) const override
    {
        record.SetMass(0.0);
        record.SetFourMomentum({1.0, 0.0, 0.0, 1.0});
        record.SetInitialPosition({0.0, 0.0, 0.0});
        record.SetInteractionVertex({0.0, 0.0, 0.0});
        record.SetHelicity(0.0);
        record.SetInitialTime(0.0);
        record.SetInteractionTime(0.0);
    }

    double GenerationProbability(
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        InteractionRecord const &) const override
    {
        return 1.0;
    }

    std::string Name() const override { return "CompletePrimary"; }

    std::shared_ptr<siren::distributions::PrimaryInjectionDistribution>
    clone() const override
    {
        return std::make_shared<CompletePrimaryDistribution>();
    }

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        InteractionRecord const &) const override
    {
        return {siren::math::Vector3D(0.0, 0.0, 0.0),
                siren::math::Vector3D(0.0, 0.0, 0.0)};
    }

protected:
    bool equal(
        siren::distributions::WeightableDistribution const & other) const override
    {
        return dynamic_cast<CompletePrimaryDistribution const *>(&other) != nullptr;
    }

    bool less(
        siren::distributions::WeightableDistribution const &) const override
    {
        return false;
    }
};

siren::dataclasses::InteractionSignature SignatureWithSecondaries(size_t count) {
    siren::dataclasses::InteractionSignature signature;
    signature.secondary_types.assign(
        count, siren::dataclasses::ParticleType::unknown);
    return signature;
}

class MixedSignatureCrossSection final
    : public siren::interactions::CrossSection {
public:
    explicit MixedSignatureCrossSection(
        std::vector<std::string> density_variables = {"q2"})
        : signatures_{SignatureWithSecondaries(2), SignatureWithSecondaries(3)}
        , density_variables_(std::move(density_variables)) {}

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
        return density_variables_;
    }

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
    std::vector<std::string> density_variables_;
};

namespace {

} // anonymous namespace

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

        double expected_density = reference.Density(q2) / (2.0 * M_PI);
        EXPECT_NEAR(density, expected_density,
                    1e-3 * expected_density);
    }

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

TEST(KinematicInjectionFailure, DirectedStepsRejectSubThresholdParent) {

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

    try {
        siren::injection::detail::SampleDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Volume,
            random);
        FAIL() << "Expected a forbidden two-body split";
    } catch (siren::utilities::InjectionFailure const & failure) {
        EXPECT_EQ(failure.reason(),
                  siren::utilities::FailureReason::KinematicallyForbidden);
    }
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
    try {
        siren::injection::detail::SampleAngularSectorStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, bin, random);
        FAIL() << "Expected a forbidden angular-sector split";
    } catch (siren::utilities::InjectionFailure const & failure) {
        EXPECT_EQ(failure.reason(),
                  siren::utilities::FailureReason::KinematicallyForbidden);
    }
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

TEST(DirectedRejectionExhaustion, SharedStepRejectsInjectionAttempt) {

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
    try {
        siren::injection::detail::SampleDirectedStep(
            parent_energy, 0.0, 0.0, parent_momentum,
            parent_mass, daughter_mass, other_mass,
            vertex, *target, target_volume,
            siren::injection::DetectorDirected2BodyChannel::Mode::Cone,
            random);
        FAIL() << "Expected directed rejection exhaustion";
    } catch (siren::utilities::InjectionFailure const & failure) {
        EXPECT_EQ(failure.reason(),
                  siren::utilities::FailureReason::SamplingFailure);
    }
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

TEST(DirectedGeometryVolume, VolumeEstimateAgreesWithAnalyticVolume) {

    siren::geometry::Sphere sphere(
        siren::geometry::Placement(siren::math::Vector3D(1.0, -2.0, 3.0)), 10.0, 0.0);
    auto first = siren::injection::EstimateGeometryVolume(sphere);
    auto second = siren::injection::EstimateGeometryVolume(sphere);
    EXPECT_EQ(first.volume, second.volume);
    double exact = siren::injection::ExactGeometryVolume(sphere);
    EXPECT_LT(std::abs(first.volume - exact), 5.0 * first.standard_error);
    EXPECT_LT(first.standard_error, 0.01 * exact);
    EXPECT_TRUE(first.resolved);
}

namespace {

class CorruptedCrossingsGeometry : public siren::geometry::Geometry {
    int mode_;
    mutable unsigned long calls_ = 0;
    siren::geometry::Box box_;
public:
    explicit CorruptedCrossingsGeometry(int mode)
        : Geometry("CorruptedCrossingsGeometry"), mode_(mode), box_(2, 2, 2) {}
    std::shared_ptr<siren::geometry::Geometry> create() const override {
        return std::make_shared<CorruptedCrossingsGeometry>(mode_);
    }
    siren::geometry::AABB GetBoundingBox() const override { return box_.GetBoundingBox(); }
    std::vector<Intersection> ComputeIntersections(
        siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const override {
        auto hits = box_.ComputeIntersections(position, direction);
        if (++calls_ % 97 == 0 && hits.size() == 2) {
            if (mode_ == 0) hits.back().distance = std::numeric_limits<double>::quiet_NaN();
            if (mode_ == 1) hits.back().distance = std::numeric_limits<double>::infinity();
            if (mode_ == 2) {
                auto extra = hits.back();
                extra.distance += 1.0;
                extra.entering = true;
                hits.push_back(extra);
            }
        }
        return hits;
    }
protected:
    bool equal(Geometry const &) const override { return false; }
    bool less(Geometry const &) const override { return false; }
    void print(std::ostream & os) const override { os << "CorruptedCrossingsGeometry"; }
};
} // namespace

TEST(DirectedGeometryVolume, CorruptedCrossingsCannotBeChecked) {

    for (int mode = 0; mode < 3; ++mode) {
        CorruptedCrossingsGeometry geometry(mode);
        auto estimate = siren::injection::EstimateGeometryVolume(geometry);
        EXPECT_FALSE(estimate.resolved) << "mode " << mode;
        EXPECT_FALSE(estimate.reason.empty()) << "mode " << mode;
        for (double volume : {1e-6, 1.0, 4.0, 8.0}) {
            try {
                siren::injection::ResolveDetectorDirectedVolume(geometry, true, volume);
                ADD_FAILURE() << "mode " << mode << " accepted volume " << volume;
            } catch (std::runtime_error const & error) {
                EXPECT_NE(std::string(error.what()).find("cannot be checked"), std::string::npos)
                    << error.what();
            }
        }
    }

    siren::geometry::BooleanGeometry same(
        siren::geometry::BooleanOperation::UNION,
        std::make_shared<siren::geometry::Box>(2, 2, 2), std::make_shared<siren::geometry::Box>(1, 1, 1));
    EXPECT_NEAR(siren::injection::ResolveDetectorDirectedVolume(same, true, 8.0), 8.0, 0.0);
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
