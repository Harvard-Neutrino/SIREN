/**
 * Unit test for DMesonELoss (charm-meson energy-loss "cross section").
 *
 * DMesonELoss models a D meson losing a fraction z of its energy to a hadronic
 * vertex. The signature is  D -> {same D, Hadrons}  on a PPlus target. The
 * inelasticity z is drawn from a Gaussian (z0=0.56, sigma=0.2) truncated to
 * [z_min_, z_max_] = [0.001, 0.999] via rejection; the same truncated Gaussian
 * is the density returned by DifferentialCrossSection / FinalStateProbability,
 * so Sample == Density (closure).
 *
 * Tests:
 *  1. EnergyMonotonicallyDecreases -- outgoing D energy is always in
 *     (D mass, primary energy); the D never gains energy.
 *  2. FourMomentumAndMassInvariants -- outgoing D is on its mass shell,
 *     energy is conserved (E0 == E_D + E_H), and the two outgoing 3-momenta
 *     are collinear with the incoming direction (momentum is NOT conserved by
 *     design -- the hadron is massless and collinear).
 *  3. TotalCrossSectionPositiveAndThrows -- TotalCrossSection is positive and
 *     increasing in E over the validity range, and throws on an unsupported
 *     primary.
 *  4. SubThresholdThrows -- SampleFinalState throws siren::utilities::
 *     InjectionFailure when the primary energy is at or below the D mass.
 *  5. ZDensityClosure -- the empirical sampled-z distribution matches
 *     FinalStateProbability/DifferentialCrossSection (the truncated Gaussian)
 *     bin-by-bin within Poisson error (closure of the realized sampler against
 *     the advertised density).
 */
#include <cmath>
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

// --- Test 1: the D meson always loses energy --------------------------------

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
            // The D meson never gains energy (z >= z_min_ > 0).
            EXPECT_LT(E_out, E_D);
            // It stays on (or above) its mass shell -- no z>z_max_ leakage.
            EXPECT_GE(E_out, mD);
            EXPECT_GT(E_out, 0.0);
        }
    }
}

// --- Test 2: on-shell, energy conservation, collinearity --------------------

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

        // (a) outgoing D is on its mass shell.
        double m2 = pdm[0]*pdm[0] - pdm[1]*pdm[1] - pdm[2]*pdm[2] - pdm[3]*pdm[3];
        EXPECT_NEAR(std::sqrt(std::max(0.0, m2)), mD, 1e-6);

        // (b) energy is conserved: E0 == E_D_out + E_H.
        EXPECT_NEAR(pdm[0] + ph[0], E_D, 1e-6);
        // hadron carries the lost energy and is (nearly) massless.
        EXPECT_GT(ph[0], 0.0);
        double mh2 = ph[0]*ph[0] - ph[1]*ph[1] - ph[2]*ph[2] - ph[3]*ph[3];
        EXPECT_NEAR(mh2, 0.0, 1e-6);

        // (c) collinearity: both outgoing 3-momenta lie along +z (the incoming
        //     direction). Momentum is NOT conserved by design (massless hadron),
        //     so we assert direction, not p-balance.
        EXPECT_NEAR(pdm[1], 0.0, 1e-9);
        EXPECT_NEAR(pdm[2], 0.0, 1e-9);
        EXPECT_GT(pdm[3], 0.0);
        EXPECT_NEAR(ph[1], 0.0, 1e-9);
        EXPECT_NEAR(ph[2], 0.0, 1e-9);
        EXPECT_GT(ph[3], 0.0);
        // direction of p_D matches the incoming +z direction.
        EXPECT_GT(p_D, 0.0);
    }
}

// --- Test 3: total cross section positivity / monotonicity / throws ---------

TEST(DMesonELoss, TotalCrossSectionPositiveAndThrows) {
    DMesonELoss xs;
    // Positive and increasing with energy over the (>1 PeV) validity range.
    double last = -1.0;
    for (double E : {1e6, 1e7, 1e8, 1e9}) {
        double s = xs.TotalCrossSection(ParticleType::D0, E);
        EXPECT_GT(s, 0.0);
        if (last >= 0.0) EXPECT_GT(s, last);
        last = s;
    }
    // Unsupported primary throws.
    EXPECT_THROW(xs.TotalCrossSection(ParticleType::PiPlus, 1e7), std::runtime_error);
}

// --- Test 4: sub-threshold primary fails recoverably ------------------------

TEST(DMesonELoss, SubThresholdThrows) {
    DMesonELoss xs;
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::D0, ParticleType::PPlus)[0];
    auto rng = std::make_shared<SIREN_random>();
    double mD = Constants::D0Mass;

    // Primary energy exactly at the D mass (D at rest): no valid final state.
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
    // Primary energy below the D mass.
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

// --- Test 5: sampled-z density matches FinalStateProbability ----------------

TEST(DMesonELoss, ZDensityClosure) {
    DMesonELoss xs;
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::D0, ParticleType::PPlus)[0];
    double E_D = 1000.0;   // well above threshold, so the truncation is the only cut
    auto rng = std::make_shared<SIREN_random>();

    // FinalStateProbability returns dxs/txs == normalized truncated Gaussian in
    // z over [z_min_, z_max_]. Histogram z = 1 - E_out/E_D and compare the
    // empirical pdf to FinalStateProbability evaluated at each bin center.
    const double zlo = 0.001, zhi = 0.999;
    const int NB = 20;
    const double bw = (zhi - zlo) / NB;
    std::vector<long> counts(NB, 0);
    const int N = 200000;

    for (int i = 0; i < N; ++i) {
        InteractionRecord rec = make_d0_record(sig, E_D);
        CrossSectionDistributionRecord cdr(rec);
        xs.SampleFinalState(cdr, rng);
        cdr.Finalize(rec);
        double E_out = rec.secondary_momenta[0][0];
        double z = 1.0 - E_out / E_D;
        // All accepted samples lie in [z_min_, z_max_].
        EXPECT_GE(z, zlo - 1e-9);
        EXPECT_LE(z, zhi + 1e-9);
        int b = (int)((z - zlo) / bw);
        if (b < 0) b = 0;
        if (b >= NB) b = NB - 1;
        counts[b]++;
    }

    // Build a record at a chosen z and read back the normalized z-density.
    auto fsp_at_z = [&](double z) -> double {
        double E_out = E_D * (1.0 - z);
        double mD = Constants::D0Mass;
        double p_out = std::sqrt(std::max(0.0, E_out * E_out - mD * mD));
        InteractionRecord r;
        r.signature = sig;
        r.primary_mass = mD;
        r.primary_momentum = {E_D, 0, 0, std::sqrt(E_D * E_D - mD * mD)};
        r.target_mass = Constants::protonMass;
        r.secondary_momenta = {{E_out, 0, 0, p_out}, {E_D - E_out, 0, 0, E_D - E_out}};
        r.secondary_masses = {mD, 0.0};
        // FinalStateProbability = DifferentialCrossSection / TotalCrossSection,
        // i.e. the normalized truncated Gaussian density in z.
        return xs.FinalStateProbability(r);
    };

    // (a) The density is non-negative everywhere and the empirical pdf closes
    //     bin-by-bin within ~4 Poisson sigma.
    for (int b = 0; b < NB; ++b) {
        double zc = zlo + (b + 0.5) * bw;
        double pred = fsp_at_z(zc);
        EXPECT_GE(pred, 0.0);
        if (counts[b] < 50) continue;
        double emp = (double)counts[b] / (N * bw);
        double sigma = std::sqrt((double)counts[b]) / (N * bw);
        EXPECT_NEAR(emp, pred, 4.0 * sigma + 0.02 * pred);
    }

    // (b) The density integrates to 1 over [z_min_, z_max_] (trapezoid on a
    //     fine grid) -- it is a properly normalized pdf in z.
    int M = 2000;
    double integral = 0.0;
    double dz = (zhi - zlo) / M;
    for (int i = 0; i <= M; ++i) {
        double z = zlo + i * dz;
        double w = (i == 0 || i == M) ? 0.5 : 1.0;
        integral += w * fsp_at_z(z) * dz;
    }
    EXPECT_NEAR(integral, 1.0, 2e-2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
