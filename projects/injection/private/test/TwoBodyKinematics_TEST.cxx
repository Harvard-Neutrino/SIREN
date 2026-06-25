#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/injection/InvariantMassMapping.h"

using namespace siren::injection;

// Physical constants (GeV)
static const double M_N4 = 0.140;     // HNL mass
static const double M_MU = 0.10566;   // muon mass
static const double M_PION = 0.13957;  // charged pion mass
static const double M_ELEC = 0.000511; // electron mass
static const double M_Z  = 91.1876;   // Z boson mass
static const double G_Z  = 2.4952;    // Z boson width

// ================================================================== //
//  TwoBodyKinematics tests                                            //
// ================================================================== //

TEST(TwoBodyKinematics, KallenFunction) {
    // Kallen(a,b,c) = a^2+b^2+c^2-2(ab+ac+bc)
    // For equal b=c: Kallen(a,b,b) = a^2 + 2b^2 - 2ab - 2ab - 2b^2 = a^2 - 4ab
    double a = 1.0, b = 0.1;
    double lam = Kallen(a, b, b);
    EXPECT_NEAR(lam, a * a - 4.0 * a * b, 1e-14);

    // Known case: Kallen(1, 0, 0) = 1
    EXPECT_NEAR(Kallen(1.0, 0.0, 0.0), 1.0, 1e-14);
}

TEST(TwoBodyKinematics, RestMomentum) {
    // N4 -> mu + pi
    double p = TwoBodyRestMomentum(M_N4, M_MU, M_PION);
    // Should be > 0 if M_N4 > M_MU + M_PION... but 0.140 < 0.105 + 0.140.
    // This decay is kinematically forbidden!
    // Use a heavier N4.
    double M_heavy = 0.300;
    p = TwoBodyRestMomentum(M_heavy, M_MU, M_PION);
    EXPECT_GT(p, 0);

    double E_mu = TwoBodyRestEnergy(M_heavy, M_MU, M_PION);
    EXPECT_NEAR(E_mu * E_mu - p * p, M_MU * M_MU, 1e-10);
}

TEST(TwoBodyKinematics, BoostAndInverseRoundTrip) {
    // Generate a rest-frame angle, boost to lab, then solve back.
    // Should recover the original angle.
    double M_parent = 0.300;
    double m_A = M_MU;
    double m_B = M_PION;

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    // Parent has lab energy 2 GeV (highly boosted)
    double E_parent_lab = 2.0;
    double p_parent_lab = std::sqrt(E_parent_lab * E_parent_lab - M_parent * M_parent);
    double beta = p_parent_lab / E_parent_lab;
    double gamma = E_parent_lab / M_parent;

    for (double cos_rest : {-0.9, -0.5, 0.0, 0.3, 0.7, 0.99}) {
        LabFrameResult lab = BoostToLab(beta, gamma, p_rest, E_rest, cos_rest);

        auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, m_A, lab.cos_theta_lab);

        bool found = false;
        for (auto & sol : solutions) {
            if (sol.valid && std::abs(sol.cos_theta_rest - cos_rest) < 1e-8) {
                found = true;
                EXPECT_NEAR(sol.p_lab, lab.p_lab, 1e-8);
            }
        }
        EXPECT_TRUE(found) << "Failed to recover cos_rest=" << cos_rest;
    }
}

TEST(TwoBodyKinematics, MasslessNoCriticalAngle) {
    // For massless daughters (like photons from ALP decay),
    // all lab angles should be accessible.
    double M_parent = 0.100;
    double m_A = 0.0;
    double m_B = 0.0;

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    double E_parent_lab = 1.0;
    double p_parent_lab = std::sqrt(E_parent_lab * E_parent_lab - M_parent * M_parent);
    double beta = p_parent_lab / E_parent_lab;
    double gamma = E_parent_lab / M_parent;

    double cos_crit = CriticalCosTheta(beta, gamma, p_rest, E_rest, m_A);
    EXPECT_EQ(cos_crit, -1.0);

    // Only one solution for massless daughters
    auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, m_A, 0.5);
    int n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);
    EXPECT_EQ(n_valid, 1);
}

TEST(TwoBodyKinematics, MassiveTwoSolutions) {
    // For massive daughters when beta_parent > beta_daughter_rest,
    // lab angles below the critical angle should have two solutions.
    double M_parent = 0.300;
    double m_A = M_MU;  // 0.106 GeV
    double m_B = M_PION;  // 0.140 GeV

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    // High boost: E_parent = 5 GeV
    double E_parent_lab = 5.0;
    double p_parent_lab = std::sqrt(E_parent_lab * E_parent_lab - M_parent * M_parent);
    double beta = p_parent_lab / E_parent_lab;
    double gamma = E_parent_lab / M_parent;

    double beta_rest = p_rest / E_rest;

    // Check if two-solution regime exists
    if (beta > beta_rest) {
        double cos_crit = CriticalCosTheta(beta, gamma, p_rest, E_rest, m_A);
        EXPECT_GT(cos_crit, -1.0);
        EXPECT_LT(cos_crit, 1.0);

        // At a forward angle (well inside critical), expect two solutions
        auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, m_A, 0.999);
        int n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);
        EXPECT_EQ(n_valid, 2);
    }
}

TEST(TwoBodyKinematics, JacobianConsistency) {
    // Verify Jacobian by numerical differentiation.
    // The Jacobian dOmega_lab/dOmega_rest should satisfy:
    //   d(cos_theta_lab)/d(cos_theta_rest) = J * (p_rest / p_lab)^2
    double M_parent = 0.300;
    double m_A = M_MU;
    double m_B = M_PION;

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    double E_parent_lab = 2.0;
    double p_parent_lab = std::sqrt(E_parent_lab * E_parent_lab - M_parent * M_parent);
    double beta = p_parent_lab / E_parent_lab;
    double gamma = E_parent_lab / M_parent;

    double cos_rest = 0.3;
    double eps = 1e-6;

    LabFrameResult lab1 = BoostToLab(beta, gamma, p_rest, E_rest, cos_rest - eps);
    LabFrameResult lab2 = BoostToLab(beta, gamma, p_rest, E_rest, cos_rest + eps);

    double numerical_dcoslab_dcosrest = (lab2.cos_theta_lab - lab1.cos_theta_lab) / (2 * eps);

    // Get the analytic Jacobian at this point
    LabFrameResult lab_mid = BoostToLab(beta, gamma, p_rest, E_rest, cos_rest);
    auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, m_A, lab_mid.cos_theta_lab);

    bool checked = false;
    for (auto & sol : solutions) {
        if (sol.valid && std::abs(sol.cos_theta_rest - cos_rest) < 1e-5) {
            // The Jacobian is dOmega_lab/dOmega_rest.
            // Since phi is unchanged, this equals d(cos_theta_lab)/d(cos_theta_rest).
            // But our Jacobian includes the p_lab^2/p_rest^2 factor for solid angle.
            // So: J_solid_angle = |d(cos_theta_lab)/d(cos_theta_rest)|
            //                   = jacobian * (p_rest / p_lab)^2  [no, that's wrong]
            //
            // Actually the solid angle Jacobian is:
            // dOmega_lab/dOmega_rest = |d(cos_theta_lab)/d(cos_theta_rest)| * |d(phi_lab)/d(phi_rest)|
            // For a boost along z, d(phi_lab)/d(phi_rest) = 1.
            // So dOmega_lab/dOmega_rest = |d(cos_theta_lab)/d(cos_theta_rest)|
            //
            // Wait, our Jacobian sol.jacobian is defined as |dOmega_lab/dOmega_rest|
            // = p_lab^2 / (gamma * p_rest * |p_lab - beta * E_lab * cos_theta_lab|)
            // This should equal |d(cos_theta_lab)/d(cos_theta_rest)| since phi is unchanged.
            EXPECT_NEAR(sol.jacobian, std::abs(numerical_dcoslab_dcosrest), 1e-4)
                << "Jacobian mismatch at cos_rest=" << cos_rest;
            checked = true;
        }
    }
    EXPECT_TRUE(checked);
}


// ================================================================== //
//  InvariantMassMapping tests                                         //
// ================================================================== //

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

    // Numerical integration of density over [s_min, s_max]
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
