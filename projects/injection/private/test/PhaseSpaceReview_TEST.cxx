#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

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


namespace {


} // anonymous namespace


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


} // namespace
