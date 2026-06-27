/**
 * Contract test for QuarkDISFromSpline (slow-rescaling charm DIS sampler).
 *
 * Invariant under test: SampleFinalState samples (xi, y) AND an independently-
 * sampled fragmentation z and a uniform azimuth phi that set the D-meson
 * momentum, but the advertised density (DensityVariables / FinalStateProbability /
 * DifferentialCrossSection) accounts for (xi, y) only. The omitted z/phi factors
 * cancel in the weight ratio ONLY in the standard unbiased configuration (the
 * same cross-section object supplies both the injection and physical densities
 * and no biased phase-space channel is installed on the D kinematics).
 *
 * These tests do NOT require a spline file: the no-arg ctor only calls
 * normalize_pdf() and compute_cdf() for the fragmentation pdf.
 *
 * Tests:
 *  1. ContractPinsTwoDensityVariables -- DensityVariables() must contain exactly
 *     {"Bjorken xi", "Bjorken y"}. This is an intentional TRIPWIRE: if a future
 *     change adds the fragmentation-z (or phi) factor to the density, this test
 *     MUST be updated in lockstep, which forces the closure implications to be
 *     reconsidered (the z/phi factors then no longer simply cancel).
 *  2. FragmentationPdfNormalized -- the fragmentation pdf sample_pdf(z) is a
 *     properly normalized density over (0.001, 0.999): its integral is 1. This
 *     documents that the z density EXISTS and is normalized, i.e. the
 *     alternative (carry z in the density) path is feasible.
 */
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/interactions/QuarkDISFromSpline.h"

using namespace siren::interactions;

// --- Test 1: density-variable contract (tripwire) ---------------------------

TEST(QuarkDISDensityContract, ContractPinsTwoDensityVariables) {
    // The no-arg ctor builds only the fragmentation pdf/cdf -- no spline needed.
    QuarkDISFromSpline xs;
    std::vector<std::string> vars = xs.DensityVariables();

    // Tripwire: the density covers exactly (xi, y). The independently-sampled
    // fragmentation z and azimuth phi are deliberately NOT in the density. If
    // this assertion ever needs to change, the closure argument (z/phi cancel in
    // the weight ratio) must be re-derived -- do not simply bump the count.
    ASSERT_EQ(vars.size(), 2u);
    EXPECT_EQ(vars[0], "Bjorken xi");
    EXPECT_EQ(vars[1], "Bjorken y");
}

// --- Test 2: fragmentation pdf is normalized --------------------------------

TEST(QuarkDISDensityContract, FragmentationPdfNormalized) {
    QuarkDISFromSpline xs;

    // sample_pdf(z) divides the unnormalized fragmentation integrand by
    // fragmentation_integral (the integral of that integrand over [0.001,0.999]),
    // so it must integrate to 1 over the same interval. Composite-trapezoid on a
    // fine grid; the integrand is smooth and bounded away from the z=0,1 poles.
    const double zlo = 0.001, zhi = 0.999;
    const int M = 20000;
    const double dz = (zhi - zlo) / M;
    double integral = 0.0;
    for (int i = 0; i <= M; ++i) {
        double z = zlo + i * dz;
        double w = (i == 0 || i == M) ? 0.5 : 1.0;
        double f = xs.sample_pdf(z);
        EXPECT_GE(f, 0.0);   // a pdf is non-negative on its support
        integral += w * f * dz;
    }
    EXPECT_NEAR(integral, 1.0, 1e-3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
