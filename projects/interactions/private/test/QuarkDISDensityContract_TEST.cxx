// Contract tests for QuarkDISFromSpline (slow-rescaling charm DIS sampler).
// SampleFinalState also draws fragmentation z and azimuth phi, but the density
// covers (xi, y) only; z/phi cancel in the unbiased weight ratio. See
// QuarkDISFromSpline.h. No spline file needed: the no-arg ctor only builds the
// fragmentation pdf/cdf.
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/interactions/QuarkDISFromSpline.h"

using namespace siren::interactions;

TEST(QuarkDISDensityContract, ContractPinsTwoDensityVariables) {
    QuarkDISFromSpline xs;
    std::vector<std::string> vars = xs.DensityVariables();

    // Tripwire: density covers exactly (xi, y); z and phi are deliberately NOT
    // in it. If this must change, re-derive the closure argument (z/phi cancel
    // in the weight ratio) -- do not simply bump the count.
    ASSERT_EQ(vars.size(), 2u);
    EXPECT_EQ(vars[0], "Bjorken xi");
    EXPECT_EQ(vars[1], "Bjorken y");
}

// The fragmentation pdf sample_pdf(z) integrates to 1 over [0.001, 0.999].
TEST(QuarkDISDensityContract, FragmentationPdfNormalized) {
    QuarkDISFromSpline xs;

    // Composite-trapezoid on a fine grid; integrand is smooth and bounded away
    // from the z=0,1 poles.
    const double zlo = 0.001, zhi = 0.999;
    const int M = 20000;
    const double dz = (zhi - zlo) / M;
    double integral = 0.0;
    for (int i = 0; i <= M; ++i) {
        double z = zlo + i * dz;
        double w = (i == 0 || i == M) ? 0.5 : 1.0;
        double f = xs.sample_pdf(z);
        EXPECT_GE(f, 0.0);   // pdf non-negative on its support
        integral += w * f * dz;
    }
    EXPECT_NEAR(integral, 1.0, 1e-3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
