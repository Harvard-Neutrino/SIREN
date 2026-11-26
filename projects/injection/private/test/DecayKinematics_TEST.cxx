#include <gtest/gtest.h>

#include "PhaseSpaceTestFixtures.h"
#include "SIREN/math/FourMomentum.h"
#include "SIREN/math/Kinematics.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/dataclasses/RecordKinematics.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>

namespace {

using siren::dataclasses::InteractionRecord;
using siren::injection::test::TwoBodyDecayRecord;

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

TEST(RecordKinematics, ReadWriteAndValidationUseOneLayout) {
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_momentum = {5.0, 1.0, 2.0, 3.0};
    record.interaction_vertex = {4.0, 5.0, 6.0};

    ASSERT_TRUE(siren::dataclasses::HasSecondaryStorage(record, 2));
    auto primary = siren::dataclasses::PrimaryMomentum(record);
    EXPECT_DOUBLE_EQ(primary.Energy(), 5.0);
    EXPECT_DOUBLE_EQ(primary.Momentum().GetX(), 1.0);
    EXPECT_DOUBLE_EQ(primary.Momentum().GetY(), 2.0);
    EXPECT_DOUBLE_EQ(primary.Momentum().GetZ(), 3.0);
    auto vertex = siren::dataclasses::InteractionVertex(record);
    EXPECT_DOUBLE_EQ(vertex.GetX(), 4.0);
    EXPECT_DOUBLE_EQ(vertex.GetY(), 5.0);
    EXPECT_DOUBLE_EQ(vertex.GetZ(), 6.0);

    siren::dataclasses::SetSecondaryMomentum(record, 1, siren::math::FourMomentum::FromComponents({7.0, 8.0, 9.0, 10.0}));
    auto secondary = siren::dataclasses::SecondaryMomentum(record, 1);
    EXPECT_DOUBLE_EQ(secondary.Energy(), 7.0);
    EXPECT_DOUBLE_EQ(secondary.Momentum().GetX(), 8.0);
    EXPECT_DOUBLE_EQ(secondary.Momentum().GetY(), 9.0);
    EXPECT_DOUBLE_EQ(secondary.Momentum().GetZ(), 10.0);

    record.secondary_momenta.pop_back();
    EXPECT_FALSE(siren::dataclasses::HasSecondaryStorage(record, 2));
    EXPECT_THROW(
        siren::dataclasses::RequireSecondaryStorage(
            record, 2, "test"),
        std::runtime_error);
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

TEST(SharedLorentzBoost, RkAdapterMatchesAnalyticBoost) {
    constexpr double parent_mass = 2.0;
    constexpr double gamma = 1.5;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    const double parent_pz = gamma * beta * parent_mass;
    constexpr double daughter_mass = 0.4;
    constexpr double rest_px = 0.3;
    constexpr double rest_py = -0.2;
    constexpr double rest_pz = 0.7;
    const double rest_energy = std::sqrt(
        daughter_mass * daughter_mass
        + rest_px * rest_px + rest_py * rest_py + rest_pz * rest_pz);

    auto frame = siren::math::FourMomentum::FromMass(
        parent_mass, siren::math::Vector3D(0.0, 0.0, parent_pz));
    auto lab = siren::math::FourMomentum::FromMass(
        daughter_mass, siren::math::Vector3D(rest_px, rest_py, rest_pz))
        .BoostToLab(frame).Components();
    EXPECT_NEAR(lab[0], gamma * (rest_energy + beta * rest_pz), 1e-14);
    EXPECT_NEAR(lab[1], rest_px, 1e-14);
    EXPECT_NEAR(lab[2], rest_py, 1e-14);
    EXPECT_NEAR(lab[3], gamma * (rest_pz + beta * rest_energy), 1e-14);
    EXPECT_NEAR(
        lab[0] * lab[0]
        - lab[1] * lab[1] - lab[2] * lab[2] - lab[3] * lab[3],
        daughter_mass * daughter_mass, 1e-13);
}

TEST(OnShellDecay, RoundedEnergyTolerancePreservesInputAndFailureType) {
    siren::injection::Isotropic2BodyChannel channel;
    auto random = std::make_shared<siren::utilities::SIREN_random>(824);
    const double shell_energy = std::hypot(2., 3.);
    for (double delta : {-1.9e-5, 1.9e-5}) {
        auto record = TwoBodyDecayRecord();
        record.primary_momentum = {shell_energy * (1. + delta), 0., 0., 3.};
        const auto original = record.primary_momentum;
        EXPECT_GT(channel.Density(nullptr, record), 0.);
        channel.Sample(random, nullptr, record);
        EXPECT_EQ(record.primary_momentum, original);
        for (std::size_t i = 0; i < 4; ++i) {
            double expected = i == 0 ? shell_energy : original[i];
            EXPECT_NEAR(record.secondary_momenta[0][i] + record.secondary_momenta[1][i],
                        expected, 2e-15 * shell_energy);
        }
    }
    for (double delta : {-2.1e-5, 2.1e-5}) {
        auto record = TwoBodyDecayRecord();
        record.primary_momentum = {shell_energy * (1. + delta), 0., 0., 3.};
        const auto primary = record.primary_momentum;
        const auto daughters = record.secondary_momenta;
        EXPECT_DOUBLE_EQ(channel.Density(nullptr, record), 0.);
        try {
            channel.Sample(random, nullptr, record);
            FAIL() << "Expected an invalid decay parent";
        } catch (siren::utilities::InjectionFailure const & error) {
            EXPECT_EQ(error.reason(), siren::utilities::FailureReason::KinematicallyForbidden);
        }
        EXPECT_EQ(record.primary_momentum, primary);
        EXPECT_EQ(record.secondary_momenta, daughters);
    }
}

} // namespace
