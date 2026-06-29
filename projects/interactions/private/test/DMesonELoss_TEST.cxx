// Unit tests for DMesonELoss: D -> {same D, Hadrons} energy-loss model. The
// inelasticity z is a truncated Gaussian that is also the advertised density,
// so Sample == Density (closure). See DMesonELoss.h for the physics.
#include <cmath>
#include <algorithm>
#include <vector>
#include <memory>
#include <stdexcept>

#include <gtest/gtest.h>

#include "SIREN/interactions/DMesonELoss.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Errors.h"

using namespace siren::utilities;
using namespace siren::interactions;
using namespace siren::dataclasses;

namespace {

// Build an unfinalized D0 -> {D0, Hadrons} record at a given lab energy along z.
InteractionRecord make_d0_record(const InteractionSignature & sig, double E_D) {
    double mD = Constants::D0Mass;
    double p_D = std::sqrt(E_D * E_D - mD * mD);
    InteractionRecord rec;
    rec.signature = sig;
    rec.primary_mass = mD;
    rec.primary_momentum = {E_D, 0, 0, p_D};
    rec.primary_helicity = 0;
    rec.target_mass = Constants::protonMass;
    rec.target_helicity = 0;
    rec.secondary_momenta = {{0, 0, 0, 0}, {0, 0, 0, 0}};
    rec.secondary_masses = {mD, 0.0};
    rec.secondary_helicities = {0, 0};
    return rec;
}

} // namespace

// The D meson always loses energy: E_out in (mD, E_D).
TEST(DMesonELoss, EnergyMonotonicallyDecreases) {
    DMesonELoss xs;
    auto sigs = xs.GetPossibleSignaturesFromParents(ParticleType::D0, ParticleType::PPlus);
    ASSERT_EQ(sigs.size(), 1u);
    auto sig = sigs[0];
    EXPECT_EQ(sig.secondary_types[0], ParticleType::D0);
    EXPECT_EQ(sig.secondary_types[1], ParticleType::Hadrons);

    double mD = Constants::D0Mass;
    auto rng = std::make_shared<SIREN_random>();

    for (double E_D : {10.0, 100.0, 1000.0}) {
        for (int i = 0; i < 2000; ++i) {
            InteractionRecord rec = make_d0_record(sig, E_D);
            CrossSectionDistributionRecord cdr(rec);
            xs.SampleFinalState(cdr, rng);
            cdr.Finalize(rec);

            double E_out = rec.secondary_momenta[0][0];
            EXPECT_LT(E_out, E_D);   // never gains energy (z >= z_min_ > 0)
            EXPECT_GE(E_out, mD);    // stays on/above mass shell (no z>z_max_ leakage)
            EXPECT_GT(E_out, 0.0);
        }
    }
}

// Outgoing D on-shell, energy conserved, both products collinear with +z.
TEST(DMesonELoss, FourMomentumAndMassInvariants) {
    DMesonELoss xs;
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::D0, ParticleType::PPlus)[0];
    double mD = Constants::D0Mass;
    double E_D = 100.0;
    double p_D = std::sqrt(E_D * E_D - mD * mD);
    auto rng = std::make_shared<SIREN_random>();

    for (int i = 0; i < 3000; ++i) {
        InteractionRecord rec = make_d0_record(sig, E_D);
        CrossSectionDistributionRecord cdr(rec);
        xs.SampleFinalState(cdr, rng);
        cdr.Finalize(rec);

        const auto & pdm = rec.secondary_momenta[0];   // outgoing D
        const auto & ph  = rec.secondary_momenta[1];   // hadron

        // outgoing D on its mass shell.
        double m2 = pdm[0]*pdm[0] - pdm[1]*pdm[1] - pdm[2]*pdm[2] - pdm[3]*pdm[3];
        EXPECT_NEAR(std::sqrt(std::max(0.0, m2)), mD, 1e-6);

        // energy conserved; hadron carries the loss and is ~massless.
        EXPECT_NEAR(pdm[0] + ph[0], E_D, 1e-6);
        EXPECT_GT(ph[0], 0.0);
        double mh2 = ph[0]*ph[0] - ph[1]*ph[1] - ph[2]*ph[2] - ph[3]*ph[3];
        EXPECT_NEAR(mh2, 0.0, 1e-6);

        // Both products collinear with +z. Momentum is NOT conserved by design
        // (massless hadron), so assert direction, not p-balance.
        EXPECT_NEAR(pdm[1], 0.0, 1e-9);
        EXPECT_NEAR(pdm[2], 0.0, 1e-9);
        EXPECT_GT(pdm[3], 0.0);
        EXPECT_NEAR(ph[1], 0.0, 1e-9);
        EXPECT_NEAR(ph[2], 0.0, 1e-9);
        EXPECT_GT(ph[3], 0.0);
        EXPECT_GT(p_D, 0.0);   // incoming +z direction
    }
}

// TotalCrossSection positive, increasing in E, and throws on bad primary.
TEST(DMesonELoss, TotalCrossSectionPositiveAndThrows) {
    DMesonELoss xs;
    // Positive and increasing over the (>1 PeV) validity range.
    double last = -1.0;
    for (double E : {1e6, 1e7, 1e8, 1e9}) {
        double s = xs.TotalCrossSection(ParticleType::D0, E);
        EXPECT_GT(s, 0.0);
        if (last >= 0.0) EXPECT_GT(s, last);
        last = s;
    }
    EXPECT_THROW(xs.TotalCrossSection(ParticleType::PiPlus, 1e7), std::runtime_error);
}

// Sub-threshold primary fails recoverably with InjectionFailure.
TEST(DMesonELoss, SubThresholdThrows) {
    DMesonELoss xs;
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::D0, ParticleType::PPlus)[0];
    auto rng = std::make_shared<SIREN_random>();
    double mD = Constants::D0Mass;

    // Energy exactly at the D mass (D at rest): no valid final state.
    {
        InteractionRecord rec;
        rec.signature = sig;
        rec.primary_mass = mD;
        rec.primary_momentum = {mD, 0, 0, 0};
        rec.secondary_momenta = {{0, 0, 0, 0}, {0, 0, 0, 0}};
        rec.secondary_masses = {mD, 0.0};
        rec.secondary_helicities = {0, 0};
        CrossSectionDistributionRecord cdr(rec);
        EXPECT_THROW(xs.SampleFinalState(cdr, rng), siren::utilities::InjectionFailure);
    }
    // Energy below the D mass.
    {
        double E_D = 0.5 * mD;
        InteractionRecord rec;
        rec.signature = sig;
        rec.primary_mass = mD;
        rec.primary_momentum = {E_D, 0, 0, 0};
        rec.secondary_momenta = {{0, 0, 0, 0}, {0, 0, 0, 0}};
        rec.secondary_masses = {mD, 0.0};
        rec.secondary_helicities = {0, 0};
        CrossSectionDistributionRecord cdr(rec);
        EXPECT_THROW(xs.SampleFinalState(cdr, rng), siren::utilities::InjectionFailure);
    }
}

// Empirical sampled-z distribution matches FinalStateProbability bin-by-bin.
TEST(DMesonELoss, ZDensityClosureAcrossEnergies) {
    DMesonELoss xs;
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::D0, ParticleType::PPlus)[0];
    auto rng = std::make_shared<SIREN_random>();
    const double mD = Constants::D0Mass;
    const double zlo = 0.001;     // == z_min_

    // At low E the kinematic cut z <= 1 - mD/E truncates the Gaussian below its
    // mean (z0=0.56) -- where a fixed-interval density would mis-normalize and
    // break closure. The density applies the same E-dependent cut, so Sample ==
    // Density at every energy. Span low to high E to exercise that.
    for (double E_D : {5.0, 10.0, 50.0, 200.0, 1000.0, 3000.0}) {
        const double z_kin = 1.0 - mD / E_D;          // kinematic upper limit
        const double zhi = std::min(1.0, z_kin);      // == min(z_max_, z_kin)
        ASSERT_GT(zhi, zlo) << "E_D=" << E_D;
        const double p_D = std::sqrt(E_D * E_D - mD * mD);

        const int NB = 20;
        const double bw = (zhi - zlo) / NB;
        std::vector<long> counts(NB, 0);
        const int N = 200000;
        for (int i = 0; i < N; ++i) {
            InteractionRecord rec = make_d0_record(sig, E_D);
            CrossSectionDistributionRecord cdr(rec);
            xs.SampleFinalState(cdr, rng);
            cdr.Finalize(rec);
            double z = 1.0 - rec.secondary_momenta[0][0] / E_D;
            // sampler must respect the support the density advertises.
            EXPECT_GE(z, zlo - 1e-9) << "E_D=" << E_D;
            EXPECT_LE(z, zhi + 1e-9) << "E_D=" << E_D;
            int b = (int)((z - zlo) / bw);
            if (b < 0) b = 0;
            if (b >= NB) b = NB - 1;
            counts[b]++;
        }

        auto fsp_at_z = [&](double z) -> double {
            double E_out = E_D * (1.0 - z);
            double p_out = std::sqrt(std::max(0.0, E_out * E_out - mD * mD));
            InteractionRecord r;
            r.signature = sig;
            r.primary_mass = mD;
            r.primary_momentum = {E_D, 0, 0, p_D};
            r.target_mass = Constants::protonMass;
            r.secondary_momenta = {{E_out, 0, 0, p_out}, {E_D - E_out, 0, 0, E_D - E_out}};
            r.secondary_masses = {mD, 0.0};
            return xs.FinalStateProbability(r);
        };

        // bin-by-bin closure within ~4 Poisson sigma.
        for (int b = 0; b < NB; ++b) {
            double zc = zlo + (b + 0.5) * bw;
            double pred = fsp_at_z(zc);
            EXPECT_GE(pred, 0.0);
            if (counts[b] < 50) continue;
            double emp = (double)counts[b] / (N * bw);
            double sg = std::sqrt((double)counts[b]) / (N * bw);
            EXPECT_NEAR(emp, pred, 4.0 * sg + 0.02 * pred) << "E_D=" << E_D << " z=" << zc;
        }

        // density integrates to 1 over the realized support [zlo, zhi].
        int M = 4000;
        double integral = 0.0, dz = (zhi - zlo) / M;
        for (int i = 0; i <= M; ++i) {
            double z = zlo + i * dz;
            double w = (i == 0 || i == M) ? 0.5 : 1.0;
            integral += w * fsp_at_z(z) * dz;
        }
        EXPECT_NEAR(integral, 1.0, 2e-2) << "E_D=" << E_D;

        // density vanishes ABOVE the kinematic cut (the fixed-interval bug left
        // nonzero density where the sampler produces nothing).
        double z_above = zhi + 0.4 * (1.0 - zhi);
        if (z_above < 1.0 - 1e-9) {
            EXPECT_EQ(fsp_at_z(z_above), 0.0) << "E_D=" << E_D << " z_above=" << z_above;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
