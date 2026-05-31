#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <memory>
#include <numeric>

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
using siren::injection::PhaseSpaceTopology;
using siren::injection::PhaseSpaceMeasure;
using siren::injection::PhysicalDecayChannel;
using siren::injection::TwoBodyRestMomentum;
using siren::injection::TwoBodyRestEnergy;
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

    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Unspecified;
    }

    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::Unspecified();
    }
};

// A channel that reports a caller-chosen topology and measure (constant
// density), for driving MultiChannelPhaseSpace::Density's measure-conversion
// path directly.
class ConfigurableChannel : public siren::injection::PhaseSpaceChannel {
public:
    ConfigurableChannel(PhaseSpaceTopology topo, PhaseSpaceMeasure meas)
        : topo_(topo), meas_(meas) {}
    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}
    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override { return 1.0; }
    std::string Name() const override { return "Configurable"; }
    PhaseSpaceTopology Topology() const override { return topo_; }
    PhaseSpaceMeasure Measure() const override { return meas_; }
private:
    PhaseSpaceTopology topo_;
    PhaseSpaceMeasure meas_;
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

// ---- PropagatorMapping (t-channel 1/(Q^2 + m^2)^2 importance map) ----

TEST(InvariantMassMapping, PropagatorRoundTripAndDensity) {
    using siren::injection::PropagatorMapping;
    double m2 = 0.2 * 0.2;          // mediator mass squared (m_V2 = 200 MeV)
    double x_min = 0.0, x_max = 1.0;
    PropagatorMapping map(m2, x_min, x_max);

    // Forward(0) -> x_min, Forward(1) -> x_max
    EXPECT_NEAR(map.Forward(0.0), x_min, 1e-9);
    EXPECT_NEAR(map.Forward(1.0), x_max, 1e-9);

    // Forward and Inverse are mutual inverses across the unit interval.
    for (double r = 0.05; r < 1.0; r += 0.05) {
        double x = map.Forward(r);
        EXPECT_GE(x, x_min - 1e-9);
        EXPECT_LE(x, x_max + 1e-9);
        EXPECT_NEAR(map.Inverse(x), r, 1e-9);
    }

    // Density matches the analytic CDF slope dr/dx = g(x), checked by
    // finite difference of Inverse (the CDF).
    for (double x = 0.05; x < 1.0; x += 0.1) {
        double h = 1e-6;
        double slope = (map.Inverse(x + h) - map.Inverse(x - h)) / (2.0 * h);
        EXPECT_NEAR(map.Density(x), slope, 1e-4 * std::max(1.0, map.Density(x)));
    }

    // Density integrates to 1 over [x_min, x_max] (midpoint rule).
    int N = 200000;
    double dx = (x_max - x_min) / N;
    double integral = 0.0;
    for (int i = 0; i < N; ++i) {
        double x = x_min + (i + 0.5) * dx;
        integral += map.Density(x) * dx;
    }
    EXPECT_NEAR(integral, 1.0, 1e-4);

    // The map is forward-peaked: the median Q^2 sits far below the range
    // midpoint (propagator concentrates weight near x_min).
    EXPECT_LT(map.Forward(0.5), 0.2 * (x_max - x_min));
}

// ---- DetectorDirectedScatteringChannel in Propagator Q^2 mode ----

TEST(PhaseSpaceChannels, ScatteringPropagatorModeSampleEqualsDensity) {
    // chi (m1) + Ar(rest) -> chi'(m3) + Ar(m4), t-channel mediator m_V2.
    auto target = Sphere(Placement(Vector3D(0, 0, 5.0)), 2.0, 0.0).create();
    auto random = std::make_shared<siren::utilities::SIREN_random>(2024);

    double m_chi = 0.008, m_chi_prime = 0.050, M_Ar = 37.215, m_V2 = 0.200;

    auto make_record = [&]() {
        InteractionRecord r;
        r.signature.primary_type = ParticleType::N4;     // stand-in for chi
        r.signature.target_type = ParticleType::PPlus;
        r.signature.secondary_types = {ParticleType::N4, ParticleType::PPlus};
        r.primary_mass = m_chi;
        double E = 0.5;
        double p = std::sqrt(E * E - m_chi * m_chi);
        r.primary_momentum = {E, 0.0, 0.0, p};
        r.target_mass = M_Ar;
        r.interaction_vertex = {0.0, 0.0, 0.0};
        r.secondary_masses = {m_chi_prime, M_Ar};
        r.secondary_momenta.resize(2);
        return r;
    };

    DetectorDirectedScatteringChannel channel(
        target, 0,
        DetectorDirectedScatteringChannel::Variable::Q2,
        DetectorDirected2BodyChannel::Mode::Volume,
        DetectorDirectedScatteringChannel::Q2Mode::Propagator,
        m_V2);

    // Closure: E_g[1] = 1 when sampling from g and weighting by 1/g over
    // the marginal Q^2 measure.  Equivalently, with a self-consistent
    // Sample==Density the running mean of (analytic physical density)/g
    // is stable.  Here we check the simplest invariant: every sampled
    // event conserves 4-momentum, stays on shell, and has Density>0 at
    // exactly the sampled point, with Q^2 inside the kinematic range.
    int n = 4000;
    int n_ok = 0;
    double sum_inv_g = 0.0;   // E_g[1/g] should approach the Q^2 range width
    double q2_lo = 1e30, q2_hi = -1e30;
    for (int i = 0; i < n; ++i) {
        InteractionRecord r = make_record();
        channel.Sample(random, nullptr, r);

        // 4-momentum conservation
        std::array<double, 4> initial = {
            r.primary_momentum[0] + r.target_mass,
            r.primary_momentum[1], r.primary_momentum[2], r.primary_momentum[3]};
        std::array<double, 4> final = {0, 0, 0, 0};
        for (auto const & p : r.secondary_momenta)
            for (int k = 0; k < 4; ++k) final[k] += p[k];
        for (int k = 0; k < 4; ++k)
            EXPECT_NEAR(final[k], initial[k], 1e-6);

        // on-shell secondaries
        EXPECT_NEAR(MassSquared(r.secondary_momenta[0]), m_chi_prime * m_chi_prime, 1e-6);
        EXPECT_NEAR(MassSquared(r.secondary_momenta[1]), M_Ar * M_Ar, 1e-4);

        double g = channel.Density(nullptr, r);
        EXPECT_GT(g, 0.0);
        if (g > 0.0) { sum_inv_g += 1.0 / g; ++n_ok; }

        double Q2 = r.interaction_parameters["Q2"];
        q2_lo = std::min(q2_lo, Q2);
        q2_hi = std::max(q2_hi, Q2);
    }
    ASSERT_GT(n_ok, n / 2);

    // E_g[1/g] estimates the integral of dQ^2 over the support = (q2max - q2min).
    // The sampled spread [q2_lo, q2_hi] should bracket a positive width and
    // E_g[1/g] should be within ~10% of that width (Monte-Carlo tolerance).
    double mean_inv_g = sum_inv_g / n_ok;
    double width = q2_hi - q2_lo;
    EXPECT_GT(width, 0.0);
    EXPECT_GT(mean_inv_g, 0.0);
    EXPECT_NEAR(mean_inv_g / width, 1.0, 0.15);

    // Propagator peaking: most samples land in the lowest 20% of the Q^2 span.
    int n_low = 0;
    for (int i = 0; i < n; ++i) {
        InteractionRecord r = make_record();
        channel.Sample(random, nullptr, r);
        double Q2 = r.interaction_parameters["Q2"];
        if (Q2 < q2_lo + 0.2 * width) ++n_low;
    }
    EXPECT_GT(n_low, n / 2);
}

TEST(PhaseSpaceChannels, ScatteringPropagatorModeRejectsBadConfig) {
    auto target = Sphere(Placement(Vector3D(0, 0, 5.0)), 2.0, 0.0).create();
    // Propagator mode with non-positive mediator mass must throw.
    EXPECT_THROW(
        DetectorDirectedScatteringChannel(
            target, 0,
            DetectorDirectedScatteringChannel::Variable::Q2,
            DetectorDirected2BodyChannel::Mode::Volume,
            DetectorDirectedScatteringChannel::Q2Mode::Propagator,
            0.0),
        std::runtime_error);
    // Propagator/Tabulated mode requires Variable::Q2.
    EXPECT_THROW(
        DetectorDirectedScatteringChannel(
            target, 0,
            DetectorDirectedScatteringChannel::Variable::RecoilY,
            DetectorDirected2BodyChannel::Mode::Volume,
            DetectorDirectedScatteringChannel::Q2Mode::Propagator,
            0.2),
        std::runtime_error);
}

TEST(PhaseSpaceChannels, TopologyMismatchDetected) {
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

    // Mixing Decay2Body with Scatter2to2 is a topology mismatch
    EXPECT_THROW(mc.CommonTopology(), std::runtime_error);
    auto diagnostics = mc.ValidateChannels();
    ASSERT_FALSE(diagnostics.empty());
    EXPECT_NE(diagnostics.front().find("Topology"), std::string::npos);
}

TEST(PhaseSpaceChannels, UnspecifiedTopologyMismatchDetected) {
    MultiChannelPhaseSpace mc;
    mc.channels = {
        std::make_shared<Isotropic2BodyChannel>(0),
        std::make_shared<CustomTestChannel>()
    };
    mc.weights = {0.5, 0.5};

    // Mixing Decay2Body with Unspecified is a topology mismatch
    EXPECT_THROW(mc.CommonTopology(), std::runtime_error);
}

TEST(ConvertDensity, ThrowsOnUnsupportedMeasurePair) {
    // Two Decay2Body channels whose measures are not mutually convertible
    // under that topology (SolidAngleRest <-> MandelstamQ2 is defined only
    // for Scatter2to2).  Density() must throw rather than silently return an
    // unconverted density (former gap G6: silent identity fallback).
    MultiChannelPhaseSpace mc;
    mc.channels = {
        std::make_shared<ConfigurableChannel>(
            PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConfigurableChannel>(
            PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::MandelstamQ2())
    };
    mc.weights = {0.5, 0.5};

    InteractionRecord record;
    EXPECT_THROW(mc.Density(nullptr, record), std::runtime_error);
}

TEST(ConvertDensity, SameMeasureMixtureDoesNotConvertOrThrow) {
    // When every channel shares the measure, no conversion is attempted, so
    // Density() simply sums the weighted channel densities.
    MultiChannelPhaseSpace mc;
    mc.channels = {
        std::make_shared<ConfigurableChannel>(
            PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConfigurableChannel>(
            PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::SolidAngleRest())
    };
    mc.weights = {0.5, 0.5};

    InteractionRecord record;
    double d = 0.0;
    EXPECT_NO_THROW(d = mc.Density(nullptr, record));
    EXPECT_NEAR(d, 1.0, 1e-12);
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

// ================================================================== //
//  Category 1: Jacobian invertibility                                 //
// ================================================================== //

TEST(JacobianInvertibility, RestLabRoundTrip) {
    // Five kinematic configurations:
    // (M, mA, mB, gamma) covering massless daughter, equal masses,
    // threshold, moderate boost, ultra-relativistic
    struct Config {
        double M, mA, mB, gamma;
    };
    Config configs[] = {
        {1.0,  0.0,   0.2,  3.0},     // massless daughter A
        {1.0,  0.3,   0.3,  5.0},     // equal masses
        {0.7 + 1e-6, 0.3, 0.4, 2.0},  // threshold: M barely > mA+mB
        {2.0,  0.5,   0.8,  1.5},     // moderate boost
        {10.0, 0.1,   0.2, 100.0},    // ultra-relativistic
    };

    for (auto const & c : configs) {
        double beta = std::sqrt(1.0 - 1.0 / (c.gamma * c.gamma));
        double p_rest = siren::injection::TwoBodyRestMomentum(c.M, c.mA, c.mB);
        double e_rest = siren::injection::TwoBodyRestEnergy(c.M, c.mA, c.mB);

        for (int i = 0; i < 10; ++i) {
            double cos_theta_rest = -1.0 + (i + 0.5) * 2.0 / 10.0;
            auto lab = siren::injection::BoostToLab(
                beta, c.gamma, p_rest, e_rest, cos_theta_rest);
            auto solutions = siren::injection::SolveLabAngle(
                beta, c.gamma, p_rest, e_rest, c.mA, lab.cos_theta_lab);

            bool found = false;
            for (int s = 0; s < 2; ++s) {
                if (solutions[s].valid &&
                    std::abs(solutions[s].cos_theta_rest - cos_theta_rest) < 1e-6) {
                    EXPECT_NEAR(solutions[s].cos_theta_rest, cos_theta_rest, 1e-8);
                    found = true;
                }
            }
            EXPECT_TRUE(found)
                << "RestLabRoundTrip failed for M=" << c.M
                << " mA=" << c.mA << " mB=" << c.mB
                << " gamma=" << c.gamma
                << " cos_theta_rest=" << cos_theta_rest;
        }
    }
}

TEST(JacobianInvertibility, JacobianProductIsUnity) {
    struct Config {
        double M, mA, mB, gamma;
    };
    Config configs[] = {
        {1.0,  0.0,   0.2,  3.0},
        {1.0,  0.3,   0.3,  5.0},
        {0.7 + 1e-6, 0.3, 0.4, 2.0},
        {2.0,  0.5,   0.8,  1.5},
        {10.0, 0.1,   0.2, 100.0},
    };

    for (auto const & c : configs) {
        double beta = std::sqrt(1.0 - 1.0 / (c.gamma * c.gamma));
        double p_rest = siren::injection::TwoBodyRestMomentum(c.M, c.mA, c.mB);
        double e_rest = siren::injection::TwoBodyRestEnergy(c.M, c.mA, c.mB);

        for (int i = 0; i < 10; ++i) {
            double cos_theta_rest = -1.0 + (i + 0.5) * 2.0 / 10.0;
            auto lab = siren::injection::BoostToLab(
                beta, c.gamma, p_rest, e_rest, cos_theta_rest);

            // Find which solution matches
            auto solutions = siren::injection::SolveLabAngle(
                beta, c.gamma, p_rest, e_rest, c.mA, lab.cos_theta_lab);
            for (int s = 0; s < 2; ++s) {
                if (!solutions[s].valid) continue;
                if (std::abs(solutions[s].cos_theta_rest - cos_theta_rest) > 1e-8)
                    continue;

                double rest_to_lab =
                    siren::injection::phase_space_jacobian::RestFrameToLabSolidAngleJacobian(
                        beta, c.gamma, p_rest, e_rest, c.mA, lab.cos_theta_lab, s);
                double lab_to_rest =
                    siren::injection::phase_space_jacobian::LabToRestFrameSolidAngleJacobian(
                        beta, c.gamma, p_rest, e_rest, c.mA, lab.cos_theta_lab, s);

                if (rest_to_lab > 0.0 && lab_to_rest > 0.0) {
                    EXPECT_NEAR(rest_to_lab * lab_to_rest, 1.0, 1e-10)
                        << "Product != 1 for M=" << c.M
                        << " gamma=" << c.gamma
                        << " cos_theta_rest=" << cos_theta_rest;
                }
            }
        }
    }
}

TEST(JacobianInvertibility, BjorkenQ2RoundTrip) {
    double M = 0.938;
    double E = 5.0;
    int N = 20;
    for (int i = 0; i < N; ++i) {
        double x = 0.1 + (i + 0.5) * 0.8 / N;
        for (int j = 0; j < N; ++j) {
            double y = 0.2 + (j + 0.5) * 0.6 / N;
            double q2 = siren::injection::phase_space_jacobian::Q2FromBjorkenXY(
                x, y, M, E);
            double x_back = siren::injection::phase_space_jacobian::BjorkenXFromQ2Y(
                q2, y, M, E);
            EXPECT_NEAR(x_back, x, 1e-12)
                << "BjorkenQ2RoundTrip failed at x=" << x << " y=" << y;
        }
    }
}

TEST(JacobianInvertibility, BjorkenQ2JacobianProductIsUnity) {
    double M = 0.938;
    double E = 5.0;
    int N = 20;
    for (int i = 0; i < N; ++i) {
        double y = 0.1 + (i + 0.5) * 0.8 / N;
        double j_fwd = siren::injection::phase_space_jacobian::BjorkenXYToQ2YAbsJacobian(
            y, M, E);
        double j_inv = siren::injection::phase_space_jacobian::Q2YToBjorkenXYAbsJacobian(
            y, M, E);
        EXPECT_NEAR(j_fwd * j_inv, 1.0, 1e-12)
            << "BjorkenQ2 Jacobian product != 1 at y=" << y;
    }
}

TEST(JacobianInvertibility, SolidAngleMandelstamJacobianProductIsUnity) {
    // Several (s, m_beam, m_target) configurations
    struct Config {
        double s, m_beam, m_target;
    };
    // s must exceed (m_beam + m_target)^2 for physical scattering
    Config configs[] = {
        {10.0, 0.0,   0.938},    // massless beam
        {5.0,  0.106, 0.938},    // muon on proton
        {100.0, 0.5,  1.0},      // heavy particles
        {2.0,  0.01,  0.5},      // light beam, moderate target
        {50.0, 0.0,   0.0},      // both massless
    };

    for (auto const & c : configs) {
        double j_fwd =
            siren::injection::phase_space_jacobian::SolidAngleRestToMandelstamQ2AbsJacobian(
                c.s, c.m_beam, c.m_target);
        double j_inv =
            siren::injection::phase_space_jacobian::MandelstamQ2ToSolidAngleRestAbsJacobian(
                c.s, c.m_beam, c.m_target);
        if (j_fwd > 0.0 && j_inv > 0.0) {
            EXPECT_NEAR(j_fwd * j_inv, 1.0, 1e-12)
                << "SolidAngle-Mandelstam product != 1 at s=" << c.s
                << " m_beam=" << c.m_beam << " m_target=" << c.m_target;
        }
    }
}

// ================================================================== //
//  Category 2: Normalization integrals                                //
// ================================================================== //

TEST(JacobianIntegrals, SolidAngleRestToMandelstamQ2IntegralAgreement) {
    // 2->2 scattering: m_beam=0.008, m_target=37.215, E_beam=1.0
    // In the lab frame the target is at rest.
    double m_beam = 0.008;
    double m_target = 37.215;
    double E_beam = 1.0;
    double s = m_beam * m_beam + m_target * m_target + 2.0 * m_target * E_beam;
    double p_CM_sq = siren::injection::Kallen(
        s, m_beam * m_beam, m_target * m_target) / (4.0 * s);
    ASSERT_GT(p_CM_sq, 0.0);
    double Q2_max = 4.0 * p_CM_sq;  // backward scattering: cos_theta = -1

    // Isotropic density in solid angle: rho_Omega = 1/(4*pi)
    // Integral of rho_Omega dOmega over full sphere = 1.
    //
    // Converting to Q2 (azimuth-integrated):
    //   rho(Q2) = rho_Omega * 2*pi / (2*p_CM^2) = 1/(4*pi) * 2*pi / (2*p_CM^2)
    //           = 1 / (4*p_CM^2)
    //
    // Integral of rho(Q2) dQ2 from 0 to Q2_max:
    //   = (1 / (4*p_CM^2)) * 4*p_CM^2 = 1
    double rho_omega = 1.0 / (4.0 * M_PI);

    // The converted density is constant (does not depend on Q2), so
    // the integral is simply rho_Q2 * Q2_max.
    double rho_Q2 =
        siren::injection::phase_space_jacobian::SolidAngleRestDensityToMandelstamQ2Density(
            rho_omega, s, m_beam, m_target);
    double integral = rho_Q2 * Q2_max;

    EXPECT_NEAR(integral, 1.0, 5e-3)
        << "SolidAngle->MandelstamQ2 integral = " << integral;
}

// ================================================================== //
//  Category 5: Boundary and edge cases                                //
// ================================================================== //

TEST(JacobianBoundary, ZeroBoostJacobianIsUnity) {
    double beta = 0.0;
    double gamma = 1.0;
    double M = 1.0;
    double mA = 0.3;
    double mB = 0.4;
    double p_rest = siren::injection::TwoBodyRestMomentum(M, mA, mB);
    double e_rest = siren::injection::TwoBodyRestEnergy(M, mA, mB);

    for (int i = 0; i < 20; ++i) {
        double cos_theta = -1.0 + (i + 0.5) * 2.0 / 20.0;
        double jac =
            siren::injection::phase_space_jacobian::RestFrameToLabSolidAngleJacobian(
                beta, gamma, p_rest, e_rest, mA, cos_theta, 0);
        EXPECT_NEAR(jac, 1.0, 1e-10)
            << "Zero-boost Jacobian != 1 at cos_theta=" << cos_theta;
    }
}

TEST(JacobianBoundary, MasslessDaughterNoCriticalAngle) {
    // Massless daughter: all lab angles accessible, critical cos = -1
    double M = 2.0;
    double mA = 0.0;
    double mB = 0.5;
    double p_rest = siren::injection::TwoBodyRestMomentum(M, mA, mB);
    double e_rest = siren::injection::TwoBodyRestEnergy(M, mA, mB);

    // Several boost values
    double gammas[] = {1.5, 3.0, 10.0, 100.0};
    for (double gamma : gammas) {
        double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
        double cos_crit = siren::injection::CriticalCosTheta(
            beta, gamma, p_rest, e_rest, mA);
        EXPECT_DOUBLE_EQ(cos_crit, -1.0)
            << "Massless daughter should have cos_crit=-1 at gamma=" << gamma;
    }
}

TEST(JacobianBoundary, ThresholdDecaySmallMomentum) {
    double mA = 0.3;
    double mB = 0.4;
    double M = mA + mB + 1e-6;  // barely above threshold
    double p_rest = siren::injection::TwoBodyRestMomentum(M, mA, mB);
    EXPECT_GT(p_rest, 0.0);
    EXPECT_LT(p_rest, 1e-2);  // should be very small

    double e_rest = siren::injection::TwoBodyRestEnergy(M, mA, mB);
    double gamma = 2.0;
    double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));

    // SolveLabAngle should not produce NaN even at threshold
    auto lab = siren::injection::BoostToLab(
        beta, gamma, p_rest, e_rest, 0.0);
    EXPECT_TRUE(std::isfinite(lab.cos_theta_lab));
    EXPECT_TRUE(std::isfinite(lab.p_lab));

    auto solutions = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, e_rest, mA, lab.cos_theta_lab);
    for (int s = 0; s < 2; ++s) {
        if (solutions[s].valid) {
            EXPECT_TRUE(std::isfinite(solutions[s].cos_theta_rest));
            EXPECT_TRUE(std::isfinite(solutions[s].p_lab));
            EXPECT_TRUE(std::isfinite(solutions[s].jacobian));
        }
    }
}

TEST(JacobianBoundary, CriticalAngleMergePoint) {
    // Setup where beta_daughter_rest < beta_parent (two solutions exist)
    // Use a heavy daughter with a fast parent
    double M = 1.0;
    double mA = 0.9;   // heavy daughter
    double mB = 0.05;
    double p_rest = siren::injection::TwoBodyRestMomentum(M, mA, mB);
    double e_rest = siren::injection::TwoBodyRestEnergy(M, mA, mB);
    double beta_daughter_rest = p_rest / e_rest;

    // Need beta_parent > beta_daughter_rest
    double gamma = 3.0;
    double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    ASSERT_GT(beta, beta_daughter_rest)
        << "Test requires beta_parent > beta_daughter_rest";

    double cos_crit = siren::injection::CriticalCosTheta(
        beta, gamma, p_rest, e_rest, mA);
    ASSERT_GT(cos_crit, -1.0) << "Expected a finite critical angle";

    // Approach the critical angle from the accessible side (slightly
    // larger cos_theta_lab = slightly smaller angle from beam axis).
    // At the exact boundary the discriminant may be numerically zero
    // or negative, so nudge inward by a small amount.
    // Verify SolveLabAngle works near the critical angle.
    // Approach from the accessible side with decreasing offsets and
    // verify the two solutions converge.
    double offsets[] = {1e-4, 1e-6, 1e-8, 1e-10};
    double prev_gap = 1e10;
    for (double eps : offsets) {
        double cos_near = cos_crit + eps;
        if (cos_near > 1.0) cos_near = 1.0;

        auto solutions = siren::injection::SolveLabAngle(
            beta, gamma, p_rest, e_rest, mA, cos_near);
        int n_valid = 0;
        for (int s = 0; s < 2; ++s) {
            if (solutions[s].valid) {
                EXPECT_TRUE(std::isfinite(solutions[s].cos_theta_rest))
                    << "NaN near critical angle, solution " << s
                    << " eps=" << eps;
                n_valid++;
            }
        }
        if (n_valid < 2) continue;

        // The gap between the two solutions should shrink as eps -> 0
        double gap = std::abs(
            solutions[0].cos_theta_rest - solutions[1].cos_theta_rest);
        EXPECT_LT(gap, prev_gap + 1e-12)
            << "Gap should decrease approaching critical angle";
        prev_gap = gap;
    }
}

TEST(JacobianBoundary, ForwardBackwardAngles) {
    double M = 2.0;
    double mA = 0.3;
    double mB = 0.5;
    double gamma = 4.0;
    double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    double p_rest = siren::injection::TwoBodyRestMomentum(M, mA, mB);
    double e_rest = siren::injection::TwoBodyRestEnergy(M, mA, mB);

    // Forward: cos_theta_rest = +1
    {
        auto lab = siren::injection::BoostToLab(
            beta, gamma, p_rest, e_rest, 1.0);
        EXPECT_TRUE(std::isfinite(lab.cos_theta_lab));
        EXPECT_TRUE(std::isfinite(lab.p_lab));

        auto solutions = siren::injection::SolveLabAngle(
            beta, gamma, p_rest, e_rest, mA, lab.cos_theta_lab);
        bool found_forward = false;
        for (int s = 0; s < 2; ++s) {
            if (solutions[s].valid &&
                std::abs(solutions[s].cos_theta_rest - 1.0) < 1e-8) {
                found_forward = true;
                EXPECT_NEAR(solutions[s].cos_theta_rest, 1.0, 1e-10);
            }
        }
        EXPECT_TRUE(found_forward) << "Could not recover cos_theta_rest = +1";
    }

    // Backward: cos_theta_rest = -1
    {
        auto lab = siren::injection::BoostToLab(
            beta, gamma, p_rest, e_rest, -1.0);
        EXPECT_TRUE(std::isfinite(lab.cos_theta_lab));
        EXPECT_TRUE(std::isfinite(lab.p_lab));

        auto solutions = siren::injection::SolveLabAngle(
            beta, gamma, p_rest, e_rest, mA, lab.cos_theta_lab);
        bool found_backward = false;
        for (int s = 0; s < 2; ++s) {
            if (solutions[s].valid &&
                std::abs(solutions[s].cos_theta_rest - (-1.0)) < 1e-8) {
                found_backward = true;
                EXPECT_NEAR(solutions[s].cos_theta_rest, -1.0, 1e-10);
            }
        }
        EXPECT_TRUE(found_backward) << "Could not recover cos_theta_rest = -1";
    }
}

// ================================================================== //
//  Category 6: Convertibility group table                             //
// ================================================================== //

TEST(ConvertibilityGroups, ExhaustiveTable) {
    using siren::dataclasses::PhaseSpaceTopology;
    using siren::dataclasses::PhaseSpaceMeasure;
    using siren::dataclasses::MeasureConvertibilityGroup;

    // Decay2Body
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::SolidAngleRest()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::SolidAngleLab()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::Recursive2Body()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::DalitzPair()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::HelicityAngles()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::MandelstamQ2()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::BjorkenXY()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay2Body,
              PhaseSpaceMeasure::Unspecified()), -1);

    // Decay3Body
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::Recursive2Body()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::DalitzPair()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::HelicityAngles()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::SolidAngleRest()), 1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::SolidAngleLab()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::MandelstamQ2()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::BjorkenXY()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Decay3Body,
              PhaseSpaceMeasure::Unspecified()), -1);

    // DecayNBody
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::SolidAngleRest()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::SolidAngleLab()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::Recursive2Body()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::DalitzPair()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::HelicityAngles()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::MandelstamQ2()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::BjorkenXY()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::DecayNBody,
              PhaseSpaceMeasure::Unspecified()), -1);

    // Scatter2to2
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::SolidAngleRest()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::SolidAngleLab()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::MandelstamQ2()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::BjorkenXY()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::Recursive2Body()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::DalitzPair()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::HelicityAngles()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to2,
              PhaseSpaceMeasure::Unspecified()), -1);

    // Scatter2to3
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::Recursive2Body()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::DalitzPair()), 0);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::SolidAngleRest()), 1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::SolidAngleLab()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::HelicityAngles()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::MandelstamQ2()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::BjorkenXY()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Scatter2to3,
              PhaseSpaceMeasure::Unspecified()), -1);

    // Unspecified topology: all measures give -1
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::SolidAngleRest()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::SolidAngleLab()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::Recursive2Body()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::DalitzPair()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::HelicityAngles()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::MandelstamQ2()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::BjorkenXY()), -1);
    EXPECT_EQ(MeasureConvertibilityGroup(PhaseSpaceTopology::Unspecified,
              PhaseSpaceMeasure::Unspecified()), -1);
}

TEST(ConvertibilityGroups, CompatibilityIsSymmetric) {
    using siren::dataclasses::PhaseSpaceTopology;
    using siren::dataclasses::PhaseSpaceMeasure;
    using siren::dataclasses::PhaseSpaceCompatible;

    // Enumerate all topology-measure pairs
    PhaseSpaceTopology topologies[] = {
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceTopology::Decay3Body,
        PhaseSpaceTopology::DecayNBody,
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceTopology::Scatter2to3,
        PhaseSpaceTopology::Unspecified,
    };
    PhaseSpaceMeasure measures[] = {
        PhaseSpaceMeasure::SolidAngleRest(),
        PhaseSpaceMeasure::SolidAngleLab(),
        PhaseSpaceMeasure::Recursive2Body(),
        PhaseSpaceMeasure::DalitzPair(),
        PhaseSpaceMeasure::HelicityAngles(),
        PhaseSpaceMeasure::MandelstamQ2(),
        PhaseSpaceMeasure::BjorkenXY(),
        PhaseSpaceMeasure::Unspecified(),
    };

    // For every pair (T, Ma) x (T, Mb), verify symmetry
    for (auto topo : topologies) {
        for (auto ma : measures) {
            for (auto mb : measures) {
                bool ab = PhaseSpaceCompatible(topo, ma, topo, mb);
                bool ba = PhaseSpaceCompatible(topo, mb, topo, ma);
                EXPECT_EQ(ab, ba)
                    << "Asymmetric compatibility for topology="
                    << static_cast<int>(topo)
                    << " ma=" << static_cast<int>(ma.type)
                    << " mb=" << static_cast<int>(mb.type);
            }
        }
    }
}

TEST(ConvertibilityGroups, LegacyConventionRoundTrip) {
    using siren::dataclasses::PhaseSpaceConvention;
    using siren::dataclasses::PhaseSpaceTopology;
    using siren::dataclasses::PhaseSpaceMeasure;
    using siren::dataclasses::TopologyFromConvention;
    using siren::dataclasses::MeasureFromConvention;

    // RestFrameSolidAngle + n=2 -> (Decay2Body, SolidAngleRest)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::RestFrameSolidAngle, 2),
              PhaseSpaceTopology::Decay2Body);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::RestFrameSolidAngle),
              PhaseSpaceMeasure::SolidAngleRest());

    // LabFrameSolidAngle + n=2 -> (Decay2Body, SolidAngleLab)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::LabFrameSolidAngle, 2),
              PhaseSpaceTopology::Decay2Body);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::LabFrameSolidAngle),
              PhaseSpaceMeasure::SolidAngleLab());

    // Recursive2Body + n=3 -> (Decay3Body, Recursive2Body)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::Recursive2Body, 3),
              PhaseSpaceTopology::Decay3Body);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::Recursive2Body),
              PhaseSpaceMeasure::Recursive2Body());

    // Dalitz + n=3 -> (Decay3Body, DalitzPair)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::Dalitz, 3),
              PhaseSpaceTopology::Decay3Body);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::Dalitz),
              PhaseSpaceMeasure::DalitzPair());

    // HelicityAngles + n=3 -> (Decay3Body, HelicityAngles)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::HelicityAngles, 3),
              PhaseSpaceTopology::Decay3Body);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::HelicityAngles),
              PhaseSpaceMeasure::HelicityAngles());

    // BjorkenXY -> (Scatter2to2, BjorkenXY)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::BjorkenXY, 2),
              PhaseSpaceTopology::Scatter2to2);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::BjorkenXY),
              PhaseSpaceMeasure::BjorkenXY());

    // MandelstamST -> (Scatter2to2, MandelstamQ2)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::MandelstamST, 2),
              PhaseSpaceTopology::Scatter2to2);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::MandelstamST),
              PhaseSpaceMeasure::MandelstamQ2());

    // Custom -> (Unspecified, Unspecified)
    EXPECT_EQ(TopologyFromConvention(PhaseSpaceConvention::Custom, 2),
              PhaseSpaceTopology::Unspecified);
    EXPECT_EQ(MeasureFromConvention(PhaseSpaceConvention::Custom),
              PhaseSpaceMeasure::Unspecified());
}

// ================================================================== //
//  Category 7: Known analytic values                                  //
// ================================================================== //

TEST(JacobianAnalytic, KnownRestToLabValue) {
    double M = 1.0;
    double mA = 0.1;
    double mB = 0.2;
    double E_parent = 3.0;
    double gamma = E_parent / M;
    double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));

    double p_rest = siren::injection::TwoBodyRestMomentum(M, mA, mB);
    double e_rest = siren::injection::TwoBodyRestEnergy(M, mA, mB);

    // At cos_theta_rest = 0 (theta_rest = pi/2):
    //   p_parallel_rest = 0
    //   p_perp = p_rest
    //   p_lab_parallel = gamma * beta * e_rest
    //   p_lab_perp = p_rest
    //   p_lab = sqrt(p_rest^2 + gamma^2 * beta^2 * e_rest^2)
    //   E_lab = gamma * e_rest  (since p_parallel = 0)
    double cos_theta_rest = 0.0;
    auto lab = siren::injection::BoostToLab(
        beta, gamma, p_rest, e_rest, cos_theta_rest);

    double p_lab_par_expected = gamma * beta * e_rest;
    double p_lab_expected = std::sqrt(
        p_rest * p_rest + gamma * gamma * beta * beta * e_rest * e_rest);
    double E_lab_expected = gamma * e_rest;

    EXPECT_NEAR(lab.p_lab, p_lab_expected, 1e-12);
    EXPECT_NEAR(lab.E_lab, E_lab_expected, 1e-12);
    EXPECT_NEAR(lab.cos_theta_lab, p_lab_par_expected / p_lab_expected, 1e-12);

    // Verify the Jacobian from the framework matches the analytic
    // calculation via the definitions.
    //
    // From SolveLabAngle, the Jacobian dOmega_lab/dOmega_rest is:
    //   gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab| / p_lab^2
    //
    // At cos_theta_rest=0:
    //   dp_factor = |p_lab - beta * E_lab * cos_theta_lab|
    double dp_factor = std::abs(lab.p_lab - beta * lab.E_lab * lab.cos_theta_lab);
    double jac_analytic = gamma * p_rest * dp_factor / (lab.p_lab * lab.p_lab);

    // Get the framework Jacobian
    auto solutions = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, e_rest, mA, lab.cos_theta_lab);

    bool found = false;
    for (int s = 0; s < 2; ++s) {
        if (!solutions[s].valid) continue;
        if (std::abs(solutions[s].cos_theta_rest - cos_theta_rest) > 1e-8)
            continue;
        EXPECT_NEAR(solutions[s].jacobian, jac_analytic, 1e-12)
            << "Jacobian mismatch: framework=" << solutions[s].jacobian
            << " analytic=" << jac_analytic;
        found = true;
    }
    EXPECT_TRUE(found) << "Could not find matching solution for cos_theta_rest=0";
}

TEST(JacobianIntegrals, NonTrivialAngularDistributionIntegralAgreement) {
    // Use f(cos_theta) = (3/4) * (1 + cos^2(theta)) which integrates
    // to 1 over the unit sphere: integral (3/4)(1+cos^2) dOmega =
    // (3/4) * 2*pi * integral_{-1}^{1} (1+cos^2) d(cos) = (3/4)*2*pi*(8/3) = 4*pi
    // So f_Omega = (3/16*pi) * (1 + cos^2(theta)) is a normalized density.
    //
    // Convert to Q2 via the Jacobian and numerically integrate:
    //   f(Q2) = f_Omega * 2*pi / (2*p_CM^2)
    // where cos_theta = 1 - Q2/(2*p_CM^2)
    //
    // The integral of f(Q2) dQ2 from 0 to Q2_max should equal 1.

    double m_beam = 0.1;
    double m_target = 0.938;
    double E_beam = 3.0;
    double s = m_beam * m_beam + m_target * m_target + 2.0 * m_target * E_beam;
    double p_CM_sq = siren::injection::Kallen(
        s, m_beam * m_beam, m_target * m_target) / (4.0 * s);
    ASSERT_GT(p_CM_sq, 0.0);
    double Q2_max = 4.0 * p_CM_sq;

    // Numerically integrate the converted density
    int N = 10000;
    double integral = 0.0;
    double dQ2 = Q2_max / N;
    for (int i = 0; i < N; ++i) {
        double Q2 = (i + 0.5) * dQ2;
        double cos_theta = 1.0 - Q2 / (2.0 * p_CM_sq);

        // Normalized angular density: (3/(16*pi)) * (1 + cos^2)
        double f_omega = (3.0 / (16.0 * M_PI)) * (1.0 + cos_theta * cos_theta);

        // Convert to Q2 density using the Jacobian
        double f_Q2 = siren::injection::phase_space_jacobian::
            SolidAngleRestDensityToMandelstamQ2Density(
                f_omega, s, m_beam, m_target);

        integral += f_Q2 * dQ2;
    }

    EXPECT_NEAR(integral, 1.0, 1e-3)
        << "Non-trivial angular distribution integral via Q2 = " << integral;
}

// ================================================================ //
//  Auto-conversion integration tests                                 //
// ================================================================ //
//
// These tests put channels with different measures into a single
// MultiChannelPhaseSpace and verify that Density() auto-converts
// correctly.

namespace {

// Helper: build a 2-body decay record for a boosted parent.
InteractionRecord Make2BodyRecord(
    double M, double m_A, double m_B, double E_parent)
{
    InteractionRecord record;
    record.signature.primary_type = ParticleType::N4;
    record.signature.target_type = ParticleType::Decay;
    record.signature.secondary_types = {ParticleType::NuLight, ParticleType::Gamma};
    record.primary_mass = M;
    double p_parent = std::sqrt(E_parent * E_parent - M * M);
    record.primary_momentum = {E_parent, 0.0, 0.0, p_parent};
    record.primary_initial_position = {0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {m_A, m_B};
    record.secondary_momenta = {{0, 0, 0, 0}, {0, 0, 0, 0}};
    record.secondary_helicities = {0, 0};
    return record;
}

// A test channel that wraps Isotropic2BodyChannel but reports
// SolidAngleLab and evaluates its density in lab-frame measure.
class LabFrameIsotropicChannel : public siren::injection::PhaseSpaceChannel {
public:
    LabFrameIsotropicChannel() : iso_(0) {}

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> dm,
        InteractionRecord & record) const override
    {
        iso_.Sample(random, dm, record);
    }

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const & record) const override
    {
        // Isotropic in rest frame means 1/(4*pi) per dOmega_rest.
        // Convert to density per dOmega_lab using the Jacobian.
        double rest_density = 1.0 / (4.0 * M_PI);

        double M = record.primary_mass;
        double m_A = record.secondary_masses[0];
        double m_B = record.secondary_masses[1];
        double p_rest = TwoBodyRestMomentum(M, m_A, m_B);
        double E_rest = TwoBodyRestEnergy(M, m_A, m_B);

        double E_parent = record.primary_momentum[0];
        double px = record.primary_momentum[1];
        double py = record.primary_momentum[2];
        double pz = record.primary_momentum[3];
        double p_parent = std::sqrt(px*px + py*py + pz*pz);
        if (p_parent < 1e-15) return rest_density;

        double beta = p_parent / E_parent;
        double gamma = E_parent / M;

        double px_A = record.secondary_momenta[0][1];
        double py_A = record.secondary_momenta[0][2];
        double pz_A = record.secondary_momenta[0][3];
        double p_A = std::sqrt(px_A*px_A + py_A*py_A + pz_A*pz_A);
        if (p_A < 1e-15) return rest_density;

        double cos_theta_lab = (px_A*px + py_A*py + pz_A*pz) / (p_A * p_parent);

        // Find the matching Lorentz boost solution
        double E_A_lab = record.secondary_momenta[0][0];
        double p_par_lab = p_A * cos_theta_lab;
        double p_par_rest = gamma * (p_par_lab - beta * E_A_lab);
        double cos_rest_actual = (p_rest > 0) ? p_par_rest / p_rest : 0.0;

        auto solutions = siren::injection::SolveLabAngle(
            beta, gamma, p_rest, E_rest, m_A, cos_theta_lab);

        // |dOmega_lab/dOmega_rest| for the correct solution
        double J = 0.0;
        double best_dist = 1e30;
        for (auto const & sol : solutions) {
            if (!sol.valid) continue;
            double dist = std::abs(sol.cos_theta_rest - cos_rest_actual);
            if (dist < best_dist) {
                best_dist = dist;
                J = sol.jacobian;
            }
        }
        if (J <= 0.0 || !std::isfinite(J)) return 0.0;

        // rest_density / J converts from per-dOmega_rest to per-dOmega_lab
        return rest_density / J;
    }

    std::string Name() const override { return "LabFrameIsotropic"; }
    PhaseSpaceTopology Topology() const override { return PhaseSpaceTopology::Decay2Body; }
    PhaseSpaceMeasure Measure() const override { return PhaseSpaceMeasure::SolidAngleLab(); }

private:
    Isotropic2BodyChannel iso_;
};

} // anonymous namespace


TEST(AutoConversion, Decay2BodyRestPlusLabDensityAgreesWithPureRest) {
    // Two channels that sample identically (isotropic 2-body) but
    // report densities in different measures.  The MultiChannelPhaseSpace
    // should auto-convert so the combined density equals the individual
    // density (since both channels are isotropic, the combined density
    // should be 1/(4*pi) in rest-frame measure regardless of weights).

    auto rest_channel = std::make_shared<Isotropic2BodyChannel>(0);
    auto lab_channel = std::make_shared<LabFrameIsotropicChannel>();

    // Verify they report different measures
    ASSERT_EQ(rest_channel->Measure(), PhaseSpaceMeasure::SolidAngleRest());
    ASSERT_EQ(lab_channel->Measure(), PhaseSpaceMeasure::SolidAngleLab());

    // Same topology
    ASSERT_EQ(rest_channel->Topology(), lab_channel->Topology());

    MultiChannelPhaseSpace mc;
    mc.channels = {rest_channel, lab_channel};
    mc.weights = {0.5, 0.5};

    // CommonMeasure should pick SolidAngleRest (higher priority)
    EXPECT_EQ(mc.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());

    // Validation should pass (same topology, compatible measures)
    auto diags = mc.ValidateChannels();
    // Should have an info diagnostic about auto-conversion, not an error
    for (auto const & d : diags) {
        EXPECT_EQ(d.find("incompatibility"), std::string::npos)
            << "Unexpected incompatibility: " << d;
    }

    // Sample many events and verify the combined density equals 1/(4*pi)
    // for each event (since both channels sample isotropically).
    double M = 1.0, m_A = 0.1, m_B = 0.2, E = 5.0;
    InteractionRecord record = Make2BodyRecord(M, m_A, m_B, E);
    auto random = std::make_shared<siren::utilities::SIREN_random>(42);

    int N = 500;
    double max_deviation = 0.0;
    for (int i = 0; i < N; ++i) {
        InteractionRecord r = record;
        mc.Sample(random, nullptr, r);
        double combined = mc.Density(nullptr, r);
        double expected = 1.0 / (4.0 * M_PI);
        double deviation = std::abs(combined - expected) / expected;
        if (deviation > max_deviation) max_deviation = deviation;
    }

    EXPECT_LT(max_deviation, 1e-8)
        << "Combined rest+lab density deviates from 1/(4*pi) by "
        << max_deviation * 100 << "%";
}

TEST(AutoConversion, Decay2BodyMixedWeightsClosure) {
    // Closure test: sample from the multichannel with unequal weights,
    // compute the multichannel density, compute the physical density
    // (1/(4*pi)), and verify the ratio is stable (constant).

    auto rest_channel = std::make_shared<Isotropic2BodyChannel>(0);
    auto lab_channel = std::make_shared<LabFrameIsotropicChannel>();

    MultiChannelPhaseSpace mc;
    mc.channels = {rest_channel, lab_channel};
    mc.weights = {0.1, 0.9};

    double M = 0.5, m_A = 0.05, m_B = 0.1, E = 3.0;
    InteractionRecord record = Make2BodyRecord(M, m_A, m_B, E);
    auto random = std::make_shared<siren::utilities::SIREN_random>(123);

    double physical = 1.0 / (4.0 * M_PI);
    int N = 2000;
    double sum_w = 0.0;
    int n_valid = 0;

    for (int i = 0; i < N; ++i) {
        InteractionRecord r = record;
        mc.Sample(random, nullptr, r);
        double gen = mc.Density(nullptr, r);
        if (gen <= 0 || !std::isfinite(gen)) continue;
        double w = physical / gen;
        sum_w += w;
        ++n_valid;
    }

    ASSERT_GT(n_valid, N / 2) << "Too many invalid samples";

    // The average weight should be close to 1.0 (unbiased estimator
    // of the integral of the physical density, which is 1.0).
    double mean_w = sum_w / n_valid;
    EXPECT_NEAR(mean_w, 1.0, 0.05)
        << "Closure test: mean weight = " << mean_w
        << " (expected 1.0 for unbiased sampling)";
}

TEST(AutoConversion, Scatter2to2Q2PlusSolidAngleClosure) {
    // Test the SolidAngleRest <-> MandelstamQ2 auto-conversion for
    // Scatter2to2 topology.  We create two channels that sample from
    // the same isotropic (in CM) distribution but report density in
    // different measures.

    // Simple 2->2 scattering: m_beam=0.1, m_target=1.0 at E_lab=2.0
    double m_beam = 0.1;
    double m_target = 1.0;
    double E_lab = 2.0;
    double s = m_beam*m_beam + m_target*m_target + 2.0*m_target*E_lab;
    double p_CM_sq = siren::injection::Kallen(s, m_beam*m_beam, m_target*m_target) / (4.0 * s);
    double Q2_max = 4.0 * p_CM_sq;

    // Channel A: reports density per dOmega_rest = 1/(4*pi)
    // Channel B: reports density per dQ2 = 1/(4*pi) * 2*pi / (2*p_CM^2)
    //            = 1 / (4 * p_CM^2)
    // These are the same physical distribution (isotropic in CM)
    // expressed in different measures.

    double density_rest = 1.0 / (4.0 * M_PI);
    double density_Q2 = density_rest * 2.0 * M_PI / (2.0 * p_CM_sq);

    // We can verify the densities integrate to 1:
    // rest: integral of 1/(4pi) dOmega = 1
    // Q2:   integral of density_Q2 dQ2 from 0 to Q2_max = density_Q2 * 4*p_CM^2 = 1
    EXPECT_NEAR(density_Q2 * Q2_max, 1.0, 1e-10);

    // Now verify the auto-conversion: if we ask MultiChannelPhaseSpace
    // for the density in the common measure (SolidAngleRest), the
    // Q2 density should be converted back to rest-frame and yield
    // 1/(4*pi).

    // We test this by computing the conversion manually and checking
    // it matches:
    namespace J = siren::injection::phase_space_jacobian;
    double converted_back = J::MandelstamQ2DensityToSolidAngleRestDensity(
        density_Q2, s, m_beam, m_target);

    EXPECT_NEAR(converted_back, density_rest, 1e-10)
        << "Manual Q2->rest conversion: got " << converted_back
        << ", expected " << density_rest;
}

TEST(AutoConversion, CrossFactorizationRecursive2Body) {
    // Two Recursive2Body measures with different spectator/pair
    // assignments must convert through Dalitz.
    //
    // P -> 0 + 1 + 2 with masses (m0, m1, m2).
    // Factorization A: spectator=0, pair=(1,2)
    // Factorization B: spectator=2, pair=(0,1)
    //
    // A uniform density in factorization A should convert to the
    // correct density in factorization B.

    namespace J = siren::injection::phase_space_jacobian;

    double M = 1.0;
    double m0 = 0.1;
    double m1 = 0.15;
    double m2 = 0.2;

    // Create a phase-space point in the parent rest frame.
    // Factorization A: pair = (1,2), sample pair invariant mass
    // and isotropic decay angles.
    double s12_min = (m1 + m2) * (m1 + m2);
    double s12_max = (M - m0) * (M - m0);
    double s12 = 0.5 * (s12_min + s12_max);
    double sqrt_s12 = std::sqrt(s12);

    // Momenta in the pair rest frame
    double p_pair = TwoBodyRestMomentum(sqrt_s12, m1, m2);
    double p_spectator = TwoBodyRestMomentum(M, m0, sqrt_s12);

    // Build daughter momenta in parent rest frame.
    // Spectator along +z, pair system along -z.
    double cos_sub = 0.3;
    double sin_sub = std::sqrt(1.0 - cos_sub * cos_sub);

    double E_pair = std::sqrt(p_spectator * p_spectator + s12);
    double E0 = std::sqrt(p_spectator * p_spectator + m0 * m0);

    std::array<double, 4> p0 = {E0, 0.0, 0.0, p_spectator};

    // Boost daughters from pair rest frame to parent rest frame.
    // The pair moves along -z with speed beta_pair = p_spectator / E_pair.
    // Boost formula for direction -z:
    //   E' = gamma * (E - beta * pz)
    //   pz' = gamma * (pz - beta * E)
    double beta_pair = p_spectator / E_pair;
    double gamma_pair = E_pair / sqrt_s12;

    double E1_prf = std::sqrt(p_pair * p_pair + m1 * m1);
    double E2_prf = std::sqrt(p_pair * p_pair + m2 * m2);

    double p1x_prf = p_pair * sin_sub;
    double p1z_prf = p_pair * cos_sub;

    double E1_lab = gamma_pair * (E1_prf - beta_pair * p1z_prf);
    double p1z_lab = gamma_pair * (p1z_prf - beta_pair * E1_prf);
    double E2_lab = gamma_pair * (E2_prf - beta_pair * (-p1z_prf));
    double p2z_lab = gamma_pair * ((-p1z_prf) - beta_pair * E2_prf);

    std::array<double, 4> p1 = {E1_lab, p1x_prf, 0.0, p1z_lab};
    std::array<double, 4> p2 = {E2_lab, -p1x_prf, 0.0, p2z_lab};

    // Verify 4-momentum conservation
    for (int i = 0; i < 4; ++i) {
        double sum = p0[i] + p1[i] + p2[i];
        double expected = (i == 0) ? M : 0.0;
        EXPECT_NEAR(sum, expected, 1e-10)
            << "4-momentum not conserved in component " << i;
    }

    // Build the InteractionRecord
    InteractionRecord record;
    record.primary_mass = M;
    record.primary_momentum = {M, 0.0, 0.0, 0.0};
    record.secondary_masses = {m0, m1, m2};
    record.secondary_momenta = {p0, p1, p2};

    // Compute the Jacobian for each factorization at this point
    double jac_A = J::Recursive2BodyToDalitzAbsJacobian(M, m0, m1, m2, s12);

    // For factorization B: spectator=2, pair=(0,1).
    // Compute s_01 from the momenta.
    double E01 = p0[0] + p1[0];
    double px01 = p0[1] + p1[1];
    double py01 = p0[2] + p1[2];
    double pz01 = p0[3] + p1[3];
    double s01 = E01*E01 - px01*px01 - py01*py01 - pz01*pz01;

    double jac_B = J::Recursive2BodyToDalitzAbsJacobian(M, m2, m0, m1, s01);

    // The conversion factor from A to B is jac_A / jac_B:
    // density_B = density_A * (jac_A^{-1}) * jac_B = density_A * jac_B / jac_A
    // (because Recursive->Dalitz divides by jac, and Dalitz->Recursive multiplies by jac)
    double density_A = 1.0;
    double expected_B = density_A / jac_A * jac_B;

    EXPECT_GT(jac_A, 0.0);
    EXPECT_GT(jac_B, 0.0);
    EXPECT_GT(expected_B, 0.0);

    // Now test through ConvertDensity
    PhaseSpaceMeasure meas_A = PhaseSpaceMeasure::Recursive2Body(0, 1, 2);
    PhaseSpaceMeasure meas_B = PhaseSpaceMeasure::Recursive2Body(2, 0, 1);

    EXPECT_NE(meas_A, meas_B);
    EXPECT_TRUE(PhaseSpaceCompatible(
        PhaseSpaceTopology::Decay3Body, meas_A,
        PhaseSpaceTopology::Decay3Body, meas_B));

    // Build a MultiChannelPhaseSpace with two mock channels
    // that report the same density but different measures.
    // We can't easily build mock channels here, so test the
    // ConvertDensity function directly via the public
    // PhaseSpaceJacobian functions and verify consistency.

    // Direct Jacobian ratio
    double ratio_direct = jac_B / jac_A;

    // Through the two-step conversion:
    // Step 1: Recursive2Body(0,1,2) -> Dalitz(0,1,2)
    double dalitz_d = J::Recursive2BodyDensityToDalitzDensity(
        density_A, M, m0, m1, m2, s12);
    // Step 2: Dalitz -> Recursive2Body(2,0,1)
    double converted_B = J::DalitzDensityToRecursive2BodyDensity(
        dalitz_d, M, m2, m0, m1, s01);

    EXPECT_NEAR(converted_B, expected_B, 1e-10)
        << "Cross-factorization conversion: got " << converted_B
        << ", expected " << expected_B
        << " (ratio jac_B/jac_A = " << ratio_direct << ")";

    // Also verify that if both factorizations happen to produce the
    // same pair invariant mass at a symmetric point, the conversion
    // gives the correct ratio. For the general asymmetric case above,
    // jac_A != jac_B, so the ratio should differ from 1.
    if (std::abs(m0 - m2) > 1e-6) {
        EXPECT_NE(jac_A, jac_B)
            << "Jacobians should differ for asymmetric masses";
    }
}
