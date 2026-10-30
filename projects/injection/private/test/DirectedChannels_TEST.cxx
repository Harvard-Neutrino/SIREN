#include <gtest/gtest.h>
#include "../DetectorDirectedChannelUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <limits>
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
using siren::injection::PhaseSpaceTopology;
using siren::injection::PhaseSpaceMeasure;
using siren::injection::PhysicalDecayChannel;
using siren::injection::TwoBodyRestMomentum;
using siren::injection::TwoBodyRestEnergy;
using siren::math::Vector3D;

TEST(DirectedDecayBranches, ScaledInverseWeightsAndUniformScattering) {
    using namespace siren::injection;
    using namespace siren::injection::detail;
    std::array<TwoBodyLabSolution, 2> solutions{};
    for (auto & s : solutions) s.valid = true;
    for (auto scale : {1e-290, 1.0, 1e290}) {
        solutions[0].jacobian = scale;
        solutions[1].jacobian = 100.0 * scale;
        EXPECT_NEAR(InverseJacobianFirstProbability(scale, 100*scale), 100.0/101, 1e-15);
        for (int branch = 0; branch < 2; ++branch) {
            EXPECT_NEAR(DirectedRestJacobian(solutions, branch,
                DirectedBranchSelection::InverseJacobianWeighted)/scale, 100.0/101, 1e-15);
            EXPECT_DOUBLE_EQ(DirectedRestJacobian(solutions, branch,
                DirectedBranchSelection::Uniform), solutions[branch].jacobian/2);
        }
    }
    solutions[0].jacobian = solutions[1].jacobian = 0;
    EXPECT_DOUBLE_EQ(InverseJacobianFirstProbability(0, 0), .5);
    EXPECT_DOUBLE_EQ(DirectedRestJacobian(solutions, 0,
        DirectedBranchSelection::InverseJacobianWeighted), 0);
}

TEST(DirectedDecayBranches, ForwardBoostDerivativeAndBranchBalance) {
    using namespace siren::injection;
    using namespace siren::injection::detail;
    double gamma = 2, beta = std::sqrt(1-1/(gamma*gamma));
    double p = std::sqrt(.015*.015-.010*.010), e = .015;

    auto lab_cos = [&](double u) {
        double z = gamma*(p*u+beta*e);
        return z/std::sqrt(z*z+p*p*(1-u*u));
    };
    for (double c : {.9999, .99, .9, .8}) {
        auto solutions = SolveLabAngle(beta, gamma, p, e, .010, c);
        ASSERT_TRUE(solutions[0].valid);
        ASSERT_TRUE(solutions[1].valid);
        double q = DirectedRestJacobian(solutions, 0,
            DirectedBranchSelection::InverseJacobianWeighted);
        double branch0 = InverseJacobianFirstProbability(
            solutions[0].jacobian, solutions[1].jacobian);
        for (int i = 0; i < 2; ++i) {
            double u = solutions[i].cos_theta_rest, h = 1e-6;
            double derivative = std::abs((lab_cos(u+h)-lab_cos(u-h))/(2*h));
            EXPECT_NEAR(derivative, solutions[i].jacobian, 1e-7*derivative);
            EXPECT_NEAR(q/derivative, i == 0 ? branch0 : 1-branch0, 1e-7);
        }
    }
}

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
    EXPECT_EQ(channel.Measure().type, PhaseSpaceMeasure::Type::Recursive2Body);
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
    EXPECT_EQ(channel.Measure().type, PhaseSpaceMeasure::Type::MandelstamQ2Phi);
}

TEST(InvariantMassMapping, PropagatorRoundTripAndDensity) {
    using siren::injection::PropagatorMapping;
    double m2 = 0.2 * 0.2;          // mediator mass squared (m_V2 = 200 MeV)
    double x_min = 0.0, x_max = 1.0;
    PropagatorMapping map(m2, x_min, x_max);

    EXPECT_NEAR(map.Forward(0.0), x_min, 1e-9);
    EXPECT_NEAR(map.Forward(1.0), x_max, 1e-9);

    for (double r = 0.05; r < 1.0; r += 0.05) {
        double x = map.Forward(r);
        EXPECT_GE(x, x_min - 1e-9);
        EXPECT_LE(x, x_max + 1e-9);
        EXPECT_NEAR(map.Inverse(x), r, 1e-9);
    }

    for (double x = 0.05; x < 1.0; x += 0.1) {
        double h = 1e-6;
        double slope = (map.Inverse(x + h) - map.Inverse(x - h)) / (2.0 * h);
        EXPECT_NEAR(map.Density(x), slope, 1e-4 * std::max(1.0, map.Density(x)));
    }

    int N = 200000;
    double dx = (x_max - x_min) / N;
    double integral = 0.0;
    for (int i = 0; i < N; ++i) {
        double x = x_min + (i + 0.5) * dx;
        integral += map.Density(x) * dx;
    }
    EXPECT_NEAR(integral, 1.0, 1e-4);

    EXPECT_LT(map.Forward(0.5), 0.2 * (x_max - x_min));
}

TEST(PhaseSpaceChannels, ScatteringPropagatorModeSampleEqualsDensity) {

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

    double E1_lab = 0.5;
    double s_ana = m_chi * m_chi + M_Ar * M_Ar + 2.0 * M_Ar * E1_lab;
    double root_s = std::sqrt(s_ana);
    double E1cm = (s_ana + m_chi * m_chi - M_Ar * M_Ar) / (2.0 * root_s);
    double E3cm = (s_ana + m_chi_prime * m_chi_prime - M_Ar * M_Ar) / (2.0 * root_s);
    double p1cm = std::sqrt(E1cm * E1cm - m_chi * m_chi);
    double p3cm = std::sqrt(E3cm * E3cm - m_chi_prime * m_chi_prime);
    double base = m_chi * m_chi + m_chi_prime * m_chi_prime - 2.0 * E1cm * E3cm;
    double q2_min_ana = -(base + 2.0 * p1cm * p3cm);   // cos* = +1
    double q2_max_ana = -(base - 2.0 * p1cm * p3cm);   // cos* = -1
    ASSERT_GT(q2_max_ana, q2_min_ana);

    int n = 4000;
    int n_ok = 0;
    double sum_inv_g = 0.0;
    for (int i = 0; i < n; ++i) {
        InteractionRecord r = make_record();
        channel.Sample(random, nullptr, r);

        std::array<double, 4> initial = {
            r.primary_momentum[0] + r.target_mass,
            r.primary_momentum[1], r.primary_momentum[2], r.primary_momentum[3]};
        std::array<double, 4> final = {0, 0, 0, 0};
        for (auto const & p : r.secondary_momenta)
            for (int k = 0; k < 4; ++k) final[k] += p[k];
        for (int k = 0; k < 4; ++k)
            EXPECT_NEAR(final[k], initial[k], 1e-6);

        EXPECT_NEAR(MassSquared(r.secondary_momenta[0]), m_chi_prime * m_chi_prime, 1e-6);
        EXPECT_NEAR(MassSquared(r.secondary_momenta[1]), M_Ar * M_Ar, 1e-4);

        double g = channel.Density(nullptr, r);
        EXPECT_GT(g, 0.0);
        if (g > 0.0) { sum_inv_g += 1.0 / g; ++n_ok; }

        double Q2 = r.interaction_parameters["Q2"];
        EXPECT_GE(Q2, q2_min_ana - 1e-9);
        EXPECT_LE(Q2, q2_max_ana + 1e-9);
    }
    ASSERT_GT(n_ok, n / 2);

    double mean_inv_g = sum_inv_g / n_ok;
    double width = q2_max_ana - q2_min_ana;
    EXPECT_NEAR(mean_inv_g / (2.0 * M_PI * width), 1.0, 0.15);

    int n_low = 0;
    for (int i = 0; i < n; ++i) {
        InteractionRecord r = make_record();
        channel.Sample(random, nullptr, r);
        double Q2 = r.interaction_parameters["Q2"];
        if (Q2 < q2_min_ana + 0.2 * width) ++n_low;
    }
    EXPECT_GT(n_low, n / 2);
}

TEST(PhaseSpaceChannels, ScatteringPropagatorModeRejectsBadConfig) {
    auto target = Sphere(Placement(Vector3D(0, 0, 5.0)), 2.0, 0.0).create();

    EXPECT_THROW(
        DetectorDirectedScatteringChannel(
            target, 0,
            DetectorDirectedScatteringChannel::Variable::Q2,
            DetectorDirected2BodyChannel::Mode::Volume,
            DetectorDirectedScatteringChannel::Q2Mode::Propagator,
            0.0),
        std::runtime_error);

    EXPECT_THROW(
        DetectorDirectedScatteringChannel(
            target, 0,
            DetectorDirectedScatteringChannel::Variable::RecoilY,
            DetectorDirected2BodyChannel::Mode::Volume,
            DetectorDirectedScatteringChannel::Q2Mode::Propagator,
            0.2),
        std::runtime_error);
}

TEST(MultiChannelDensity, DensityBreakdownSumsToDensity) {

    auto target = Sphere(Placement(Vector3D(0.0, 0.0, 0.0)), 10.0, 0.0).create();
    auto random = std::make_shared<siren::utilities::SIREN_random>(7);

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

    MultiChannelPhaseSpace mc;
    mc.channels = {
        std::make_shared<Isotropic2BodyChannel>(0),
        std::make_shared<DetectorDirected2BodyChannel>(
            target, 0, DetectorDirected2BodyChannel::Mode::Volume)
    };
    mc.weights = {0.3, 0.7};

    for (int s = 0; s < 200; ++s) {
        mc.Sample(random, nullptr, record);
        auto breakdown = mc.DensityBreakdown(nullptr, record);
        ASSERT_EQ(breakdown.size(), mc.channels.size());

        double g = mc.Density(nullptr, record);
        double sum = std::accumulate(breakdown.begin(), breakdown.end(), 0.0);
        EXPECT_NEAR(sum, g, 1e-12 * std::max(1.0, std::abs(g)));

        for (size_t i = 0; i < mc.channels.size(); ++i) {
            double gi = mc.channels[i]->Density(nullptr, record);
            EXPECT_NEAR(breakdown[i], mc.weights[i] * gi,
                        1e-12 * std::max(1.0, std::abs(breakdown[i])));
        }
    }
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

    auto target = Sphere(Placement(Vector3D(0.0, 0.0, 0.0)), 10.0, 9.9999).create();

    EXPECT_THROW(
        DetectorDirected2BodyChannel(target, 0, DetectorDirected2BodyChannel::Mode::Volume),
        std::runtime_error
    );

    EXPECT_NO_THROW(
        DetectorDirected2BodyChannel(target, 0, DetectorDirected2BodyChannel::Mode::Cone)
    );
}

namespace {

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

        double E_A_lab = record.secondary_momenta[0][0];
        double p_par_lab = p_A * cos_theta_lab;
        double p_par_rest = gamma * (p_par_lab - beta * E_A_lab);
        double cos_rest_actual = (p_rest > 0) ? p_par_rest / p_rest : 0.0;

        auto solutions = siren::injection::SolveLabAngle(
            beta, gamma, p_rest, E_rest, m_A, cos_theta_lab);

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

        return rest_density / J;
    }

    std::string Name() const override { return "LabFrameIsotropic"; }
    PhaseSpaceTopology Topology() const override { return PhaseSpaceTopology::Decay2Body; }
    PhaseSpaceMeasure Measure() const override { return PhaseSpaceMeasure::SolidAngleLab(); }

private:
    Isotropic2BodyChannel iso_;
};

} // anonymous namespace

namespace {

std::vector<std::shared_ptr<siren::injection::PhaseSpaceChannel>>
TwoCompatibleChannels() {
    return {
        std::make_shared<Isotropic2BodyChannel>(0),
        std::make_shared<Isotropic2BodyChannel>(0),
    };
}

} // namespace

