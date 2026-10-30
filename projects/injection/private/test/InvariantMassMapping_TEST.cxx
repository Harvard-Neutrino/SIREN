#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/injection/InvariantMassMapping.h"

using namespace siren::injection;

static const double M_N4 = 0.140;     // HNL mass
static const double M_MU = 0.10566;   // muon mass
static const double M_PION = 0.13957;  // charged pion mass
static const double M_ELEC = 0.000511; // electron mass
static const double M_Z  = 91.1876;   // Z boson mass
static const double G_Z  = 2.4952;    // Z boson width

TEST(InvariantMassMapping, BreitWignerRoundTrip) {
    BreitWignerMapping bw(M_Z, G_Z, 10.0 * 10.0, 200.0 * 200.0);

    for (double r : {0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99}) {
        double s = bw.Forward(r);
        double r_back = bw.Inverse(s);
        EXPECT_NEAR(r_back, r, 1e-12);
    }
}

TEST(InvariantMassMapping, BreitWignerPeakNearResonance) {
    BreitWignerMapping bw(M_Z, G_Z, 10.0 * 10.0, 200.0 * 200.0);

    double density_on_peak = bw.Density(M_Z * M_Z);
    double density_off_peak = bw.Density(50.0 * 50.0);

    EXPECT_GT(density_on_peak, density_off_peak * 100);
}

TEST(InvariantMassMapping, BreitWignerDensityNormalized) {
    BreitWignerMapping bw(M_Z, G_Z, 80.0 * 80.0, 100.0 * 100.0);

    int N = 100000;
    double integral = 0.0;
    double ds = (bw.s_max - bw.s_min) / N;
    for (int i = 0; i < N; ++i) {
        double s = bw.s_min + (i + 0.5) * ds;
        integral += bw.Density(s) * ds;
    }
    EXPECT_NEAR(integral, 1.0, 1e-3);
}

TEST(InvariantMassMapping, PowerLawRoundTrip) {
    PowerLawMapping pl(0.8, 0.0, 0.01, 100.0);

    for (double r : {0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99}) {
        double s = pl.Forward(r);
        double r_back = pl.Inverse(s);
        EXPECT_NEAR(r_back, r, 1e-12);
    }
}

TEST(InvariantMassMapping, PowerLawDensityNormalized) {
    PowerLawMapping pl(0.8, 0.0, 0.01, 100.0);

    int N = 100000;
    double integral = 0.0;
    double ds = (pl.s_max - pl.s_min) / N;
    for (int i = 0; i < N; ++i) {
        double s = pl.s_min + (i + 0.5) * ds;
        integral += pl.Density(s) * ds;
    }
    EXPECT_NEAR(integral, 1.0, 1e-3);
}

TEST(InvariantMassMapping, UniformRoundTrip) {
    UniformMapping u(1.0, 10.0);

    for (double r : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        double s = u.Forward(r);
        double r_back = u.Inverse(s);
        EXPECT_NEAR(r_back, r, 1e-14);
    }
    EXPECT_NEAR(u.Density(5.0), 1.0 / 9.0, 1e-14);
}
