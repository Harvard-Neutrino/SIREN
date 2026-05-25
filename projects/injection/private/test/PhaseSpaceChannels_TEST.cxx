#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
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
}
