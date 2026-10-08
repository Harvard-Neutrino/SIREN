#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <array>
#include <cmath>
#include <memory>
#include <vector>

using siren::dataclasses::InteractionRecord;
using siren::dataclasses::ParticleType;
using siren::geometry::Box;
using siren::geometry::Geometry;
using siren::geometry::Placement;
using siren::injection::DetectorDirectedAngularSectorChannel;
using siren::injection::Isotropic2BodyChannel;
using siren::injection::MultiChannelPhaseSpace;
using siren::injection::PhaseSpaceChannel;
using siren::math::Vector3D;

namespace {

constexpr double kTwoPi = 2.0 * M_PI;

// V1 -> chi chi record (mirrors tests/python/test_directed_regimes.py).
InteractionRecord MakeV1Record(double m_chi, double E_V1) {
    InteractionRecord record;
    record.signature.primary_type = ParticleType::N4;
    record.signature.target_type = ParticleType::Decay;
    record.signature.secondary_types = {ParticleType::NuLight, ParticleType::Gamma};
    record.primary_mass = 0.017;
    double pz = std::sqrt(std::max(E_V1 * E_V1 - 0.017 * 0.017, 0.0));
    record.primary_momentum = {E_V1, 0.0, 0.0, pz};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {m_chi, m_chi};
    record.secondary_momenta.resize(2);
    record.secondary_helicities = {0, 0};
    return record;
}

std::shared_ptr<Geometry const> BoxAt(double x, double y, double z, double w) {
    return Box(Placement(Vector3D(x, y, z)), w, w, w).create();
}

// Closure E_g[f/g] for a 50/50 (isotropic, directed) mixture, with f the
// isotropic physical density.  -> 1 iff the directed channel's Density is the
// true density its Sample draws from (Contract C1).
double Closure(std::shared_ptr<PhaseSpaceChannel> directed,
               double m_chi, double E_V1, int N, unsigned seed) {
    auto iso = std::make_shared<Isotropic2BodyChannel>(0);
    auto mc = std::make_shared<MultiChannelPhaseSpace>();
    mc->channels = {iso, directed};
    mc->weights = {0.5, 0.5};
    auto rng = std::make_shared<siren::utilities::SIREN_random>(seed);

    double sum_w = 0.0;
    int n = 0;
    for (int i = 0; i < N; ++i) {
        InteractionRecord record = MakeV1Record(m_chi, E_V1);
        mc->Sample(rng, nullptr, record);
        double f = iso->Density(nullptr, record);
        double g = mc->Density(nullptr, record);
        if (g > 0.0 && std::isfinite(g) && std::isfinite(f)) {
            sum_w += f / g;
            ++n;
        }
    }
    return (n > 0) ? sum_w / n : std::nan("");
}

}  // namespace

// A single full-cone sector (u in [0,1], phi in [0,2pi]) covers the whole
// bounding cone; closure must hold.
TEST(DirectedAngularSector, FullConeClosure) {
    auto target = BoxAt(0.0, 0.0, 20.0, 4.0);
    auto sector = std::make_shared<DetectorDirectedAngularSectorChannel>(
        target, 0.0, 1.0, 0.0, kTwoPi, 0);
    double c = Closure(sector, 0.008, 0.018, 20000, 101);
    EXPECT_GT(c, 0.9);
    EXPECT_LT(c, 1.1);
}

// A 2x4 partition of sectors tiles the cone disjointly; the mixture closure
// must hold (the sectors sum to the full-cone density).
TEST(DirectedAngularSector, PartitionClosure) {
    auto target = BoxAt(0.0, 0.0, 20.0, 4.0);
    auto inner = std::make_shared<MultiChannelPhaseSpace>();
    int n_theta = 2, n_phi = 4;
    for (int it = 0; it < n_theta; ++it) {
        double u_lo = double(it) / n_theta;
        double u_hi = double(it + 1) / n_theta;
        for (int ip = 0; ip < n_phi; ++ip) {
            double phi_lo = kTwoPi * ip / n_phi;
            double phi_hi = kTwoPi * (ip + 1) / n_phi;
            inner->channels.push_back(
                std::make_shared<DetectorDirectedAngularSectorChannel>(
                    target, u_lo, u_hi, phi_lo, phi_hi, 0));
        }
    }
    inner->weights.assign(inner->channels.size(), 1.0 / inner->channels.size());

    // Wrap as one channel via the mixture's own density.
    struct Wrap : PhaseSpaceChannel {
        std::shared_ptr<MultiChannelPhaseSpace> mc;
        void Sample(std::shared_ptr<siren::utilities::SIREN_random> r,
                    std::shared_ptr<siren::detector::DetectorModel const> d,
                    InteractionRecord & rec) const override { mc->Sample(r, d, rec); }
        double Density(std::shared_ptr<siren::detector::DetectorModel const> d,
                       InteractionRecord const & rec) const override { return mc->Density(d, rec); }
        std::string Name() const override { return "Wrap"; }
        siren::injection::PhaseSpaceTopology Topology() const override { return mc->CommonTopology(); }
        siren::injection::PhaseSpaceMeasure Measure() const override { return mc->CommonMeasure(); }
    };
    auto wrap = std::make_shared<Wrap>();
    wrap->mc = inner;

    double c = Closure(wrap, 0.008, 0.018, 20000, 202);
    EXPECT_GT(c, 0.9);
    EXPECT_LT(c, 1.1);
}

// Sample == Density: a direction drawn by a sector lands in that sector's own
// bin, so the sector reports a positive density there; conservation holds.
TEST(DirectedAngularSector, SelfDensityPositiveAndConserved) {
    auto target = BoxAt(0.0, 0.0, 20.0, 4.0);
    auto sector = std::make_shared<DetectorDirectedAngularSectorChannel>(
        target, 0.0, 0.5, 0.0, M_PI, 0);  // inner polar half, first azimuth half
    auto rng = std::make_shared<siren::utilities::SIREN_random>(303);

    int n_pos = 0, N = 2000;
    for (int i = 0; i < N; ++i) {
        InteractionRecord record = MakeV1Record(0.008, 0.018);
        sector->Sample(rng, nullptr, record);
        // momentum conservation: secondaries sum to the primary
        std::array<double, 4> sum{0, 0, 0, 0};
        for (auto const & p : record.secondary_momenta)
            for (int k = 0; k < 4; ++k) sum[k] += p[k];
        for (int k = 0; k < 4; ++k)
            EXPECT_NEAR(sum[k], record.primary_momentum[k], 1e-6);
        if (sector->Density(nullptr, record) > 0.0) ++n_pos;
    }
    // The vast majority of self-samples are in-bin (boundary numerics aside).
    EXPECT_GT(n_pos, 0.99 * N);
}

// Constructor guards reject degenerate bins.
TEST(DirectedAngularSector, ConstructorValidation) {
    auto target = BoxAt(0.0, 0.0, 20.0, 4.0);
    EXPECT_THROW(DetectorDirectedAngularSectorChannel(target, 0.5, 0.5, 0.0, 1.0, 0),
                 std::runtime_error);   // u_lo == u_hi
    EXPECT_THROW(DetectorDirectedAngularSectorChannel(target, -0.1, 0.5, 0.0, 1.0, 0),
                 std::runtime_error);   // u_lo < 0
    EXPECT_THROW(DetectorDirectedAngularSectorChannel(target, 0.0, 1.0, 1.0, 1.0, 0),
                 std::runtime_error);   // phi_lo == phi_hi
}
