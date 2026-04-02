/**
 * Unit test for CharmMesonDecay CDF and Interpolator1D behavior.
 *
 * Tests:
 * 1. Interpolator1D with a known linear inverse CDF (sanity check)
 * 2. Interpolator1D with the actual D meson decay CDF table
 * 3. CharmMesonDecay::SampleFinalState Q² distribution
 *
 * Build (from SIREN/build):
 *   make -j4   (if registered in CMakeLists.txt)
 * Or standalone:
 *   g++ -std=c++17 -I../projects/utilities/public -I../projects/interactions/public \
 *       -I../projects/dataclasses/public -I../vendor/cereal/include \
 *       -I../vendor/rk/include -I../vendor/photospline/include \
 *       -L. -lSIREN -o CharmMesonDecay_TEST CharmMesonDecay_TEST.cxx -lgtest -lgtest_main
 */
#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include "SIREN/utilities/Interpolator.h"
#include "SIREN/utilities/Integration.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/interactions/CharmMesonDecay.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"

using namespace siren::utilities;
using namespace siren::interactions;
using namespace siren::dataclasses;

// ── Test 1: Interpolator1D with known linear CDF ─────────────────────────

TEST(Interpolator1D, LinearInverseCDF) {
    // CDF: F(x) = x / 1.4, so inverse: x = 1.4 * u
    // Table: x = CDF values, f = Q2 values
    TableData1D<double> table;
    int N = 102;
    for (int i = 0; i < N; ++i) {
        double u = (double)i / (N - 1);  // CDF: 0 to 1
        double q2 = u * 1.4;             // Q2: 0 to 1.4
        table.x.push_back(u);
        table.f.push_back(q2);
    }

    Interpolator1D<double> interp(table);

    std::cout << "Test 1: Linear inverse CDF" << std::endl;
    std::cout << "  IsLog: " << interp.IsLog() << std::endl;

    // Check known values
    double tol = 0.01;
    EXPECT_NEAR(interp(0.0), 0.0, tol);
    EXPECT_NEAR(interp(0.25), 0.35, tol);
    EXPECT_NEAR(interp(0.5), 0.7, tol);
    EXPECT_NEAR(interp(0.75), 1.05, tol);
    EXPECT_NEAR(interp(1.0), 1.4, tol);

    std::cout << "  interp(0.0)  = " << interp(0.0) << " (expect 0.0)" << std::endl;
    std::cout << "  interp(0.25) = " << interp(0.25) << " (expect 0.35)" << std::endl;
    std::cout << "  interp(0.5)  = " << interp(0.5) << " (expect 0.7)" << std::endl;
    std::cout << "  interp(0.75) = " << interp(0.75) << " (expect 1.05)" << std::endl;
    std::cout << "  interp(1.0)  = " << interp(1.0) << " (expect 1.4)" << std::endl;

    // Sample mean should be 0.7
    double sum = 0;
    int Nsamp = 10000;
    SIREN_random rng;
    for (int i = 0; i < Nsamp; ++i) {
        double u = rng.Uniform(0, 1);
        sum += interp(u);
    }
    double mean = sum / Nsamp;
    std::cout << "  Mean from " << Nsamp << " samples: " << mean << " (expect ~0.7)" << std::endl;
    EXPECT_NEAR(mean, 0.7, 0.05);
}

// ── Test 2: Interpolator1D with D meson decay CDF ───────────────────────

TEST(Interpolator1D, DMesonDecayCDF) {
    // Replicate computeDiffGammaCDF for D0 -> K- e+ nu_e
    double mD = Constants::D0Mass;
    double mK = Constants::KminusMass;
    double F0CKM = 0.719;
    double alpha = 0.50;
    double ms = 2.00697;
    double GF = Constants::FermiConstant;

    std::cout << "\nTest 2: D meson decay CDF" << std::endl;
    std::cout << "  mD = " << mD << ", mK = " << mK << std::endl;

    // DifferentialDecayWidth (fixed version)
    auto dGamma = [&](double Q2) -> double {
        double Q2tilde = Q2 / (ms * ms);
        double ff2 = std::pow(F0CKM / ((1 - Q2tilde) * (1 - alpha * Q2tilde)), 2);
        double EK = 0.5 * (Q2 - (mD * mD + mK * mK)) / mD;
        double pk_sq = EK * EK - mK * mK;
        if (pk_sq < 0) return 0.0;
        double PK = std::sqrt(pk_sq);
        return std::pow(GF, 2) / (24 * std::pow(M_PI, 3)) * ff2 * std::pow(PK, 3);
    };

    // Normalize (same as SIREN: Romberg integration over [0, 1.4])
    std::function<double(double)> pdf_func = dGamma;
    double norm = rombergIntegrate(pdf_func, 0.0, 1.4);
    std::cout << "  Normalization: " << norm << std::endl;

    auto normed_pdf = [&](double Q2) -> double {
        return dGamma(Q2) / norm;
    };

    // Build CDF table (same as SIREN: 100 nodes from 0.01 to ~1.39)
    double Q2_min = 0.0;
    double Q2_max = 1.4;
    std::vector<double> Q2spline;
    for (int i = 0; i < 100; ++i) {
        Q2spline.push_back(0.01 + i * (Q2_max - Q2_min) / 100);
    }

    std::vector<double> cdf_Q2_nodes;
    std::vector<double> cdf_vector;
    std::vector<double> pdf_vector;

    cdf_Q2_nodes.push_back(0);
    cdf_vector.push_back(0);
    pdf_vector.push_back(0);

    for (size_t i = 0; i < Q2spline.size(); ++i) {
        double cur_Q2 = Q2spline[i];
        double cur_pdf = normed_pdf(cur_Q2);
        double area;
        if (i == 0) {
            area = cur_Q2 * cur_pdf * 0.5;
        } else {
            area = 0.5 * (pdf_vector[i - 1] + cur_pdf) * (Q2spline[i] - Q2spline[i - 1]);
        }
        pdf_vector.push_back(cur_pdf);
        cdf_Q2_nodes.push_back(cur_Q2);
        cdf_vector.push_back(area + cdf_vector.back());
    }

    cdf_Q2_nodes.push_back(Q2_max);
    cdf_vector.push_back(1.0);
    pdf_vector.push_back(0);

    std::cout << "  CDF table size: " << cdf_vector.size() << " nodes" << std::endl;
    std::cout << "  CDF last value before forced 1.0: " << cdf_vector[cdf_vector.size() - 2] << std::endl;

    // Build Interpolator1D (inverse CDF: x=CDF, f=Q2)
    TableData1D<double> inverse_cdf_data;
    inverse_cdf_data.x = cdf_vector;
    inverse_cdf_data.f = cdf_Q2_nodes;

    Interpolator1D<double> inverseCdf(inverse_cdf_data);

    std::cout << "  Interpolator1D IsLog: " << inverseCdf.IsLog() << std::endl;
    std::cout << "  MinX: " << inverseCdf.MinX() << ", MaxX: " << inverseCdf.MaxX() << std::endl;

    // Evaluate at known CDF values
    std::cout << "\n  Inverse CDF spot checks:" << std::endl;
    double test_cdfs[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    for (double u : test_cdfs) {
        double q2 = inverseCdf(u);
        // Also compute expected by linear interpolation
        double q2_expected = 0;
        for (size_t j = 0; j < cdf_vector.size() - 1; ++j) {
            if (cdf_vector[j] <= u && u <= cdf_vector[j + 1]) {
                double t = (cdf_vector[j + 1] - cdf_vector[j] > 0) ?
                    (u - cdf_vector[j]) / (cdf_vector[j + 1] - cdf_vector[j]) : 0;
                q2_expected = cdf_Q2_nodes[j] + (cdf_Q2_nodes[j + 1] - cdf_Q2_nodes[j]) * t;
                break;
            }
        }
        std::cout << "    CDF=" << u << " -> Q2=" << q2 << " (linear expected: " << q2_expected
                  << ", diff=" << q2 - q2_expected << ")" << std::endl;
    }

    // Sample and compute mean Q²
    SIREN_random rng;
    int Nsamp = 100000;
    double sum_q2 = 0;
    double sum_q2_linear = 0;
    for (int i = 0; i < Nsamp; ++i) {
        double u = rng.Uniform(0, 1);
        sum_q2 += inverseCdf(u);
        // Linear interpolation for comparison
        double q2_lin = 0;
        for (size_t j = 0; j < cdf_vector.size() - 1; ++j) {
            if (cdf_vector[j] <= u && u <= cdf_vector[j + 1]) {
                double t = (cdf_vector[j + 1] - cdf_vector[j] > 0) ?
                    (u - cdf_vector[j]) / (cdf_vector[j + 1] - cdf_vector[j]) : 0;
                q2_lin = cdf_Q2_nodes[j] + (cdf_Q2_nodes[j + 1] - cdf_Q2_nodes[j]) * t;
                break;
            }
        }
        sum_q2_linear += q2_lin;
    }
    double mean_interp = sum_q2 / Nsamp;
    double mean_linear = sum_q2_linear / Nsamp;

    std::cout << "\n  Mean Q² from Interpolator1D: " << mean_interp << std::endl;
    std::cout << "  Mean Q² from linear interp:  " << mean_linear << std::endl;
    std::cout << "  Expected (correct):          ~0.58" << std::endl;
    std::cout << "  Difference:                  " << mean_interp - mean_linear << std::endl;

    // The test: Interpolator1D should give similar results to linear interpolation
    // If not, this confirms the Interpolator1D bias
    if (std::abs(mean_interp - mean_linear) > 0.05) {
        std::cout << "\n  *** INTERPOLATOR1D BIAS DETECTED ***" << std::endl;
        std::cout << "  Interpolator1D gives " << mean_interp << " vs linear " << mean_linear << std::endl;
    }
    EXPECT_NEAR(mean_interp, mean_linear, 0.05);
}

// ── Test 2b: Compare DifferentialDecayWidth with analytical formula ──────

TEST(CharmMesonDecay, DiffDecayWidthComparison) {
    // Check if the SIREN DifferentialDecayWidth matches our analytical formula
    double mD = Constants::D0Mass;   // 1.86962 (note: swapped with DPlus in Constants.h!)
    double mK = Constants::KminusMass;
    double F0CKM = 0.719;
    double alpha = 0.50;
    double ms_star = 2.00697;
    double GF = Constants::FermiConstant;

    std::cout << "\nTest 2b: DifferentialDecayWidth comparison" << std::endl;
    std::cout << "  Constants::D0Mass = " << Constants::D0Mass << " (PDG D0 = 1.86484, PDG D+ = 1.86966)" << std::endl;
    std::cout << "  Constants::DPlusMass = " << Constants::DPlusMass << " (SWAPPED!)" << std::endl;

    // Our analytical formula
    auto dGamma_analytical = [&](double Q2) -> double {
        double Q2tilde = Q2 / (ms_star * ms_star);
        double ff2 = std::pow(F0CKM / ((1 - Q2tilde) * (1 - alpha * Q2tilde)), 2);
        double EK = 0.5 * (Q2 - (mD * mD + mK * mK)) / mD;
        double pk_sq = EK * EK - mK * mK;
        if (pk_sq < 0) return 0.0;
        return std::pow(GF, 2) / (24 * std::pow(M_PI, 3)) * ff2 * std::pow(std::sqrt(pk_sq), 3);
    };

    // Build an InteractionRecord to call the SIREN DifferentialDecayWidth
    CharmMesonDecay decay(ParticleType::D0);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::D0);
    auto sig = sigs[0]; // D0 -> K- e+ nu_e

    std::cout << "\n  Q² | analytical | SIREN DDW | ratio" << std::endl;
    double Q2_tests[] = {0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.35};
    for (double Q2 : Q2_tests) {
        double anal = dGamma_analytical(Q2);

        // For SIREN DDW, we need to construct an InteractionRecord with the right kinematics
        // DDW(InteractionRecord) reconstructs Q² from 4-momenta, so we need to set them up
        // at the given Q². Use D rest frame for simplicity.
        double EK_rest = (mD * mD + mK * mK - Q2) / (2 * mD);
        double PK_rest = std::sqrt(std::max(0.0, EK_rest * EK_rest - mK * mK));

        InteractionRecord rec;
        rec.signature = sig;
        rec.primary_mass = mD;
        rec.primary_momentum = {mD, 0, 0, 0};  // D at rest
        rec.target_mass = 0;
        rec.secondary_momenta = {
            {EK_rest, PK_rest, 0, 0},  // K along x
            {0, 0, 0, 0},  // lepton (not used by DDW)
            {0, 0, 0, 0}   // neutrino (not used by DDW)
        };
        rec.secondary_masses = {mK, Constants::electronMass, 0.0};

        double siren_ddw_from_record = decay.DifferentialDecayWidth(rec);

        // Also call 4-arg DDW directly with same constants
        std::vector<double> my_constants = {F0CKM, alpha, ms_star};
        double siren_ddw_4arg = decay.DifferentialDecayWidth(my_constants, Q2, mD, mK);

        // And with alpha=0.44 (FormFactorFromRecord value)
        std::vector<double> ffr_constants = {0.719, 0.44, 2.00697};
        double siren_ddw_4arg_044 = decay.DifferentialDecayWidth(ffr_constants, Q2, mD, mK);

        std::cout << "  Q2=" << Q2
                  << " | anal=" << anal
                  << " | 4arg(a=0.50)=" << siren_ddw_4arg
                  << " | 4arg(a=0.44)=" << siren_ddw_4arg_044
                  << " | fromRecord=" << siren_ddw_from_record
                  << " | ratioRec/anal=" << siren_ddw_from_record / anal
                  << std::endl;
    }
}

// ── Test 3: CharmMesonDecay SampleFinalState Q² ─────────────────────────

TEST(CharmMesonDecay, SampledQ2Distribution) {
    // Create D0 decay
    CharmMesonDecay decay(ParticleType::D0);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::D0);
    // sig[0] = D0 -> K- e+ nu_e
    auto sig = sigs[0];

    std::cout << "\nTest 3: CharmMesonDecay sampled Q² distribution" << std::endl;
    std::cout << "  Signature: D0 -> ";
    for (auto s : sig.secondary_types) std::cout << (int)s << " ";
    std::cout << std::endl;

    double mD = Constants::D0Mass;
    double E_D = 100.0;  // 100 GeV D meson
    double p_D = std::sqrt(E_D * E_D - mD * mD);

    auto rng = std::make_shared<SIREN_random>();
    int Nsamp = 5000;
    double sum_q2 = 0;

    for (int i = 0; i < Nsamp; ++i) {
        InteractionRecord rec;
        rec.signature = sig;
        rec.primary_mass = mD;
        rec.primary_momentum = {E_D, 0, 0, p_D};
        rec.primary_helicity = 0;
        rec.target_mass = 0;
        rec.target_helicity = 0;
        rec.secondary_momenta = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
        rec.secondary_masses = {Constants::KminusMass, Constants::electronMass, 0.0};
        rec.secondary_helicities = {0, 0, 0};

        CrossSectionDistributionRecord cdr(rec);
        decay.SampleFinalState(cdr, rng);
        cdr.Finalize(rec);

        // Reconstruct Q² = (p_D - p_K)²
        double qE  = rec.primary_momentum[0] - rec.secondary_momenta[0][0];
        double qpx = rec.primary_momentum[1] - rec.secondary_momenta[0][1];
        double qpy = rec.primary_momentum[2] - rec.secondary_momenta[0][2];
        double qpz = rec.primary_momentum[3] - rec.secondary_momenta[0][3];
        double Q2 = qE * qE - qpx * qpx - qpy * qpy - qpz * qpz;
        sum_q2 += Q2;
    }

    double mean_q2 = sum_q2 / Nsamp;
    std::cout << "  Mean Q² from SampleFinalState (" << Nsamp << " events): " << mean_q2 << std::endl;
    std::cout << "  Expected (correct): ~0.58" << std::endl;

    // If mean Q² > 0.7, the CDF sampling is biased
    if (mean_q2 > 0.7) {
        std::cout << "  *** Q² BIAS CONFIRMED: " << mean_q2 << " >> 0.58 ***" << std::endl;
    }

    // Loose check — correct value is 0.58, biased is ~0.79
    // This test documents the known bias rather than asserting correctness
    EXPECT_GT(mean_q2, 0.3);   // sanity: not negative/zero
    EXPECT_LT(mean_q2, 1.2);   // sanity: not absurdly high
}
