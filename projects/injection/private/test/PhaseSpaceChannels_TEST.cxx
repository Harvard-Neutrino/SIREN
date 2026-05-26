#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <memory>

using siren::dataclasses::InteractionRecord;
using siren::dataclasses::ParticleType;
using siren::geometry::Placement;
using siren::geometry::Sphere;
using siren::injection::DetectorDirected2BodyChannel;
using siren::injection::DetectorDirected3BodyChannel;
using siren::injection::DetectorDirectedScatteringChannel;
using siren::injection::Isotropic2BodyChannel;
using siren::injection::MultiChannelPhaseSpace;
using siren::injection::PhaseSpaceConvention;
using siren::injection::PhysicalDecayChannel;
using siren::math::Vector3D;

namespace {

double SpatialMagnitude(std::array<double, 4> const & p) {
    return std::sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
}

double MassSquared(std::array<double, 4> const & p) {
    return p[0] * p[0] - p[1] * p[1] - p[2] * p[2] - p[3] * p[3];
}

bool RayHits(
    std::shared_ptr<siren::geometry::Geometry const> const & target,
    std::array<double, 3> const & vertex,
    std::array<double, 4> const & p)
{
    double mag = SpatialMagnitude(p);
    if (mag <= 0.0) return false;
    Vector3D pos(vertex);
    Vector3D dir(p[1] / mag, p[2] / mag, p[3] / mag);
    auto hits = target->Intersections(pos, dir);
    for (auto const & hit : hits) {
        if (hit.entering && hit.distance > 0.0) return true;
    }
    return false;
}

class CustomTestChannel : public siren::injection::PhaseSpaceChannel {
public:
    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override {
        return 1.0;
    }

    std::string Name() const override {
        return "CustomTest";
    }

    PhaseSpaceConvention Convention() const override {
        return PhaseSpaceConvention::Custom;
    }
};

class MixedArityDecay : public siren::interactions::Decay {
public:
    bool equal(siren::interactions::Decay const &) const override {
        return false;
    }

    double TotalDecayWidthAllFinalStates(InteractionRecord const &) const override {
        return 0.0;
    }

    double TotalDecayWidth(siren::dataclasses::ParticleType) const override {
        return 0.0;
    }

    double TotalDecayWidth(InteractionRecord const &) const override {
        return 0.0;
    }

    double DifferentialDecayWidth(InteractionRecord const &) const override {
        return 0.0;
    }

    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override {}

    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return {TwoBodySignature(), ThreeBodySignature()};
    }

    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType) const override {
        return GetPossibleSignatures();
    }

    double FinalStateProbability(InteractionRecord const &) const override {
        return 1.0;
    }

    std::vector<std::string> DensityVariables() const override {
        return {"CosTheta"};
    }

    static siren::dataclasses::InteractionSignature TwoBodySignature() {
        siren::dataclasses::InteractionSignature sig;
        sig.primary_type = ParticleType::N4;
        sig.target_type = ParticleType::Decay;
        sig.secondary_types = {ParticleType::NuLight, ParticleType::Gamma};
        return sig;
    }

    static siren::dataclasses::InteractionSignature ThreeBodySignature() {
        siren::dataclasses::InteractionSignature sig;
        sig.primary_type = ParticleType::N4;
        sig.target_type = ParticleType::Decay;
        sig.secondary_types = {
            ParticleType::NuLight,
            ParticleType::MuMinus,
            ParticleType::MuPlus
        };
        return sig;
    }
};

void ExpectConserved(InteractionRecord const & record, double tol) {
    std::array<double, 4> sum = {0, 0, 0, 0};
    for (auto const & p : record.secondary_momenta) {
        for (int i = 0; i < 4; ++i) sum[i] += p[i];
    }
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(sum[i], record.primary_momentum[i], tol);
    }
}

} // namespace

TEST(PhaseSpaceChannels, DetectorDirected3BodySamplesConservedPointingEvent) {
    auto target = Sphere(Placement(Vector3D(0, 0, 100000)), 1.0, 0.0).create();
    auto random = std::make_shared<siren::utilities::SIREN_random>(12345);

    InteractionRecord record;
    record.signature.primary_type = ParticleType::N4;
    record.signature.secondary_types = {
        ParticleType::NuLight,
        ParticleType::MuMinus,
        ParticleType::MuPlus
    };
    record.primary_mass = 1.0;
    double E = 20.0;
    double pz = std::sqrt(E * E - record.primary_mass * record.primary_mass);
    record.primary_momentum = {E, 0.0, 0.0, pz};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {
        0.0,
        siren::utilities::Constants::muonMass,
        siren::utilities::Constants::muonMass
    };
    record.secondary_momenta.resize(3);

    DetectorDirected3BodyChannel channel(
        target,
        0,
        1,
        2,
        1,
        DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0,
        0.0,
        0.8,
        0.0,
        DetectorDirected2BodyChannel::Mode::Volume);

    channel.Sample(random, nullptr, record);

    ExpectConserved(record, 1e-8);
    for (size_t i = 0; i < record.secondary_momenta.size(); ++i) {
        EXPECT_NEAR(MassSquared(record.secondary_momenta[i]),
                    record.secondary_masses[i] * record.secondary_masses[i],
                    1e-8);
    }
    EXPECT_TRUE(RayHits(target, record.interaction_vertex, record.secondary_momenta[1]));
    EXPECT_GT(channel.Density(nullptr, record), 0.0);
    EXPECT_EQ(channel.Convention(), PhaseSpaceConvention::Recursive2Body);
}

TEST(PhaseSpaceChannels, DetectorDirectedScatteringSamplesConservedPointingEvent) {
    auto target = Sphere(Placement(Vector3D(0, 0, 100000)), 1.0, 0.0).create();
    auto random = std::make_shared<siren::utilities::SIREN_random>(67890);

    InteractionRecord record;
    record.signature.primary_type = ParticleType::NuMu;
    record.signature.target_type = ParticleType::PPlus;
    record.signature.secondary_types = {
        ParticleType::N4,
        ParticleType::PPlus
    };
    record.primary_mass = 0.0;
    double E = 20.0;
    record.primary_momentum = {E, 0.0, 0.0, E};
    record.target_mass = siren::utilities::Constants::protonMass;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.05, record.target_mass};
    record.secondary_momenta.resize(2);

    DetectorDirectedScatteringChannel channel(
        target,
        0,
        DetectorDirectedScatteringChannel::Variable::Q2,
        DetectorDirected2BodyChannel::Mode::Volume);

    channel.Sample(random, nullptr, record);

    std::array<double, 4> initial = {
        record.primary_momentum[0] + record.target_mass,
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3]
    };
    std::array<double, 4> final = {0, 0, 0, 0};
    for (auto const & p : record.secondary_momenta) {
        for (int i = 0; i < 4; ++i) final[i] += p[i];
    }
    for (int i = 0; i < 4; ++i) EXPECT_NEAR(final[i], initial[i], 1e-8);

    for (size_t i = 0; i < record.secondary_momenta.size(); ++i) {
        EXPECT_NEAR(MassSquared(record.secondary_momenta[i]),
                    record.secondary_masses[i] * record.secondary_masses[i],
                    1e-8);
    }
    EXPECT_TRUE(RayHits(target, record.interaction_vertex, record.secondary_momenta[0]));
    EXPECT_TRUE(record.interaction_parameters.count("Q2"));
    EXPECT_TRUE(record.interaction_parameters.count("bjorken_y"));
    EXPECT_GT(channel.Density(nullptr, record), 0.0);
    EXPECT_EQ(channel.Convention(), PhaseSpaceConvention::MandelstamST);
}

TEST(PhaseSpaceChannels, ConventionValidationWarnsOnNonCustomMismatch) {
    auto target = Sphere(Placement(Vector3D(0, 0, 100000)), 1.0, 0.0).create();

    MultiChannelPhaseSpace mc;
    mc.channels = {
        std::make_shared<Isotropic2BodyChannel>(0),
        std::make_shared<DetectorDirectedScatteringChannel>(
            target,
            0,
            DetectorDirectedScatteringChannel::Variable::RecoilY,
            DetectorDirected2BodyChannel::Mode::Volume)
    };
    mc.weights = {0.5, 0.5};

    auto diagnostics = mc.ValidateConventions();
    ASSERT_FALSE(diagnostics.empty());
    EXPECT_NE(diagnostics.front().find("RestFrameSolidAngle"), std::string::npos);
    EXPECT_EQ(mc.CommonConvention(), PhaseSpaceConvention::RestFrameSolidAngle);
}

TEST(PhaseSpaceChannels, ConventionValidationRejectsCustomMix) {
    MultiChannelPhaseSpace mc;
    mc.channels = {
        std::make_shared<Isotropic2BodyChannel>(0),
        std::make_shared<CustomTestChannel>()
    };
    mc.weights = {0.5, 0.5};

    EXPECT_THROW(mc.ValidateConventions(), std::runtime_error);
}

TEST(PhaseSpaceChannels, PhysicalDecayConventionFromModel) {
    auto decay = std::make_shared<MixedArityDecay>();

    PhysicalDecayChannel generic(decay);
    PhysicalDecayChannel two_body(decay, MixedArityDecay::TwoBodySignature());
    PhysicalDecayChannel three_body(decay, MixedArityDecay::ThreeBodySignature());

    // MixedArityDecay doesn't override Convention(); the base class
    // detects mixed arities and returns Custom for all constructors.
    EXPECT_EQ(generic.Convention(), PhaseSpaceConvention::Custom);
    EXPECT_EQ(two_body.Convention(), PhaseSpaceConvention::Custom);
    EXPECT_EQ(three_body.Convention(), PhaseSpaceConvention::Custom);

    // The explicit-convention constructor always honors the user's choice.
    PhysicalDecayChannel explicit_rest(decay, PhaseSpaceConvention::RestFrameSolidAngle);
    EXPECT_EQ(explicit_rest.Convention(), PhaseSpaceConvention::RestFrameSolidAngle);
}

TEST(PhaseSpaceJacobians, RestFrameAndLabSolidAngleIntegralsAgree) {
    double M = 1.0;
    double E = 3.0;
    double beta = std::sqrt(E * E - M * M) / E;
    double gamma = E / M;
    double m_a = 0.1;
    double m_b = 0.2;
    double p_rest = siren::injection::TwoBodyRestMomentum(M, m_a, m_b);
    double e_rest = siren::injection::TwoBodyRestEnergy(M, m_a, m_b);
    double f_rest = 1.0 / (4.0 * M_PI);

    int N = 4000;
    double integral = 0.0;
    for (int i = 0; i < N; ++i) {
        double cos_lab = -1.0 + (i + 0.5) * 2.0 / N;
        double f_lab = 0.0;
        for (int solution = 0; solution < 2; ++solution) {
            double lab_to_rest =
                siren::injection::phase_space_jacobian::LabToRestFrameSolidAngleJacobian(
                    beta, gamma, p_rest, e_rest, m_a, cos_lab, solution);
            f_lab += f_rest * lab_to_rest;
        }
        integral += f_lab * (2.0 / N) * (2.0 * M_PI);
    }

    EXPECT_NEAR(integral, 1.0, 5e-3);
}

TEST(PhaseSpaceJacobians, RecursiveAndDalitzIntegralsAgree) {
    double M = 2.0;
    double m_s = 0.2;
    double m_1 = 0.3;
    double m_2 = 0.4;
    double s_min = (m_1 + m_2) * (m_1 + m_2);
    double s_max = (M - m_s) * (M - m_s);
    double recursive_density = 1.0 / (2.0 * (s_max - s_min));

    int N = 2000;
    double integral = 0.0;
    for (int i = 0; i < N; ++i) {
        double s_pair = s_min + (i + 0.5) * (s_max - s_min) / N;
        double jac =
            siren::injection::phase_space_jacobian::Recursive2BodyToDalitzAbsJacobian(
                M, m_s, m_1, m_2, s_pair);
        double dalitz_density =
            siren::injection::phase_space_jacobian::Recursive2BodyDensityToDalitzDensity(
                recursive_density, M, m_s, m_1, m_2, s_pair);
        double dalitz_range = 2.0 * jac;
        integral += dalitz_density * dalitz_range * (s_max - s_min) / N;
    }

    EXPECT_NEAR(integral, 1.0, 5e-4);
    EXPECT_EQ(
        siren::injection::phase_space_jacobian::Recursive2BodyToHelicityAnglesJacobian(),
        1.0);
}

TEST(PhaseSpaceJacobians, BjorkenAndQ2YIntegralsAgree) {
    double M = 0.938;
    double E = 5.0;
    double x_min = 0.1;
    double x_max = 0.9;
    double y_min = 0.2;
    double y_max = 0.8;
    double bjorken_density = 1.0 / ((x_max - x_min) * (y_max - y_min));

    int N = 1000;
    double integral = 0.0;
    for (int i = 0; i < N; ++i) {
        double y = y_min + (i + 0.5) * (y_max - y_min) / N;
        double q2_min = siren::injection::phase_space_jacobian::Q2FromBjorkenXY(
            x_min, y, M, E);
        double q2_max = siren::injection::phase_space_jacobian::Q2FromBjorkenXY(
            x_max, y, M, E);
        double q2_density =
            siren::injection::phase_space_jacobian::BjorkenXYDensityToQ2YDensity(
                bjorken_density, y, M, E);
        integral += q2_density * (q2_max - q2_min) * (y_max - y_min) / N;
    }

    EXPECT_NEAR(integral, 1.0, 1e-10);
}

TEST(PhaseSpaceChannels, SolidAngleDensityWhenInsideVolume) {
    auto target = Sphere(Placement(Vector3D(0.0, 0.0, 0.0)), 10.0, 0.0).create();
    auto random = std::make_shared<siren::utilities::SIREN_random>(42);

    InteractionRecord record;
    record.signature.primary_type = ParticleType::N4;
    record.signature.secondary_types = {ParticleType::NuLight, ParticleType::Gamma};
    record.primary_mass = 1.0;
    double E = 20.0;
    double pz = std::sqrt(E * E - record.primary_mass * record.primary_mass);
    record.primary_momentum = {E, 0.0, 0.0, pz};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.0, 0.0};
    record.secondary_momenta.resize(2);

    DetectorDirected2BodyChannel channel(target, 0, DetectorDirected2BodyChannel::Mode::Volume);

    double sphere_volume = (4.0 / 3.0) * M_PI * 1000.0;
    channel.SetVolume(sphere_volume);

    channel.Sample(random, nullptr, record);

    double density = channel.Density(nullptr, record);
    EXPECT_GT(density, 0.0);
}

TEST(PhaseSpaceChannels, SampleVolumeWhenInsideVolume) {
    auto target = Sphere(Placement(Vector3D(0.0, 0.0, 0.0)), 10.0, 0.0).create();
    auto random = std::make_shared<siren::utilities::SIREN_random>(42);

    InteractionRecord record;
    record.signature.primary_type = ParticleType::N4;
    record.signature.secondary_types = {ParticleType::NuLight, ParticleType::Gamma};
    record.primary_mass = 1.0;
    double E = 20.0;
    double pz = std::sqrt(E * E - record.primary_mass * record.primary_mass);
    record.primary_momentum = {E, 0.0, 0.0, pz};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.0, 0.0};
    record.secondary_momenta.resize(2);

    DetectorDirected2BodyChannel channel(target, 0, DetectorDirected2BodyChannel::Mode::Volume);

    channel.Sample(random, nullptr, record);

    ExpectConserved(record, 1e-8);
    EXPECT_GT(channel.Density(nullptr, record), 0.0);
}

TEST(PhaseSpaceChannels, ThrowsEarlyOnExtremelySmallVolumeRatio) {
    // Extremely hollow sphere: R = 10.0, r = 9.9999
    auto target = Sphere(Placement(Vector3D(0.0, 0.0, 0.0)), 10.0, 9.9999).create();

    EXPECT_THROW(
        DetectorDirected2BodyChannel(target, 0, DetectorDirected2BodyChannel::Mode::Volume),
        std::runtime_error
    );

    EXPECT_NO_THROW(
        DetectorDirected2BodyChannel(target, 0, DetectorDirected2BodyChannel::Mode::Cone)
    );
}
