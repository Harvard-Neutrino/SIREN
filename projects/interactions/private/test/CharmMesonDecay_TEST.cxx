/**
 * Unit test for CharmMesonDecay CDF / Interpolator1D behavior and the
 * SampleFinalState <-> FinalStateProbability closure invariant.
 *
 * Tests:
 * 1. Interpolator1D with a known linear inverse CDF (sanity check)
 * 2. Interpolator1D with the actual D meson decay CDF table
 * 3. SampleFinalState q^2 distribution closes with FinalStateProbability
 * 4. TotalDecayWidthForFinalState throws on unsupported signatures
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
#include <stdexcept>

#include <gtest/gtest.h>

#include <rk/rk.hh>
#include <rk/geom3.hh>

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

// --- Test 1: Interpolator1D with known linear CDF ------------------------

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

    // Check known values
    double tol = 0.01;
    EXPECT_NEAR(interp(0.0), 0.0, tol);
    EXPECT_NEAR(interp(0.25), 0.35, tol);
    EXPECT_NEAR(interp(0.5), 0.7, tol);
    EXPECT_NEAR(interp(0.75), 1.05, tol);
    EXPECT_NEAR(interp(1.0), 1.4, tol);

    // Sample mean should be 0.7
    double sum = 0;
    int Nsamp = 10000;
    SIREN_random rng;
    for (int i = 0; i < Nsamp; ++i) {
        double u = rng.Uniform(0, 1);
        sum += interp(u);
    }
    double mean = sum / Nsamp;
    EXPECT_NEAR(mean, 0.7, 0.05);
}

// --- Test 2: Interpolator1D with D meson decay CDF -----------------------

TEST(Interpolator1D, DMesonDecayCDF) {
    // Build an inverse-CDF table locally from a single-pole form factor for
    // D0 -> K- e+ nu_e (self-contained Interpolator1D sanity check).
    double mD = Constants::D0Mass;
    double mK = Constants::KMinusMass;
    double F0CKM = 0.719;
    double alpha = 0.50;
    double ms = 2.00697;
    double GF = Constants::FermiConstant;

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

    // Build Interpolator1D (inverse CDF: x=CDF, f=Q2)
    TableData1D<double> inverse_cdf_data;
    inverse_cdf_data.x = cdf_vector;
    inverse_cdf_data.f = cdf_Q2_nodes;

    Interpolator1D<double> inverse_cdf(inverse_cdf_data);

    // Sample and compute mean Q^2 from the interpolator and from a manual
    // linear interpolation; the two must agree (this guards the interpolator).
    SIREN_random rng;
    int Nsamp = 100000;
    double sum_q2 = 0;
    double sum_q2_linear = 0;
    for (int i = 0; i < Nsamp; ++i) {
        double u = rng.Uniform(0, 1);
        sum_q2 += inverse_cdf(u);
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

    EXPECT_NEAR(mean_interp, mean_linear, 0.05);
}

// --- Test 3: SampleFinalState q^2 closes with FinalStateProbability -------

namespace {
// Reconstruct q^2 = (p_D - p_K)^2 from a finalized record.
double reconstruct_q2(const InteractionRecord & rec) {
    double qE  = rec.primary_momentum[0] - rec.secondary_momenta[0][0];
    double qpx = rec.primary_momentum[1] - rec.secondary_momenta[0][1];
    double qpy = rec.primary_momentum[2] - rec.secondary_momenta[0][2];
    double qpz = rec.primary_momentum[3] - rec.secondary_momenta[0][3];
    return qE * qE - qpx * qpx - qpy * qpy - qpz * qpz;
}

// Build a record at a given q^2 in the D rest frame with hadron mass mK so that
// FinalStateProbability can be evaluated on a single mixture component.
InteractionRecord make_record_at_q2(const InteractionSignature & sig,
                                    double mD, double mK, double ml, double q2) {
    double EK = (mD * mD + mK * mK - q2) / (2 * mD);
    double PK = std::sqrt(std::max(0.0, EK * EK - mK * mK));
    InteractionRecord rec;
    rec.signature = sig;
    rec.primary_mass = mD;
    rec.primary_momentum = {mD, 0, 0, 0};
    rec.target_mass = 0;
    rec.secondary_momenta = {{EK, PK, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    rec.secondary_masses = {mK, ml, 0.0};
    return rec;
}
} // namespace

TEST(CharmMesonDecay, SampledQ2Distribution) {
    // D0 -> K- e+ nu_e with kinematic K/K*(892) mixing.
    CharmMesonDecay decay(ParticleType::D0);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::D0);
    auto sig = sigs[0];

    double mD = Constants::D0Mass;
    double ml = Constants::electronMass;
    double mK = Constants::KMinusMass;
    double mKstar = Constants::KPrimePlusMass;

    double E_D = 100.0;  // boosted D meson (q^2 is frame-independent)
    double p_D = std::sqrt(E_D * E_D - mD * mD);

    auto rng = std::make_shared<SIREN_random>();
    int Nsamp = 30000;
    double sum_q2 = 0;
    double sum_q2_sq = 0;

    // Histogram q^2 separately for the K and K* sub-populations.
    const int NB = 16;
    double q2lo = 0.0;
    double q2hi = (mD - mK) * (mD - mK);  // widest support (K mass)
    double bw = (q2hi - q2lo) / NB;
    std::vector<long> countK(NB, 0), countKstar(NB, 0);

    for (int i = 0; i < Nsamp; ++i) {
        InteractionRecord rec;
        rec.signature = sig;
        rec.primary_mass = mD;
        rec.primary_momentum = {E_D, 0, 0, p_D};
        rec.primary_helicity = 0;
        rec.target_mass = 0;
        rec.target_helicity = 0;
        rec.secondary_momenta = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
        rec.secondary_masses = {mK, ml, 0.0};
        rec.secondary_helicities = {0, 0, 0};

        CrossSectionDistributionRecord cdr(rec);
        decay.SampleFinalState(cdr, rng);
        cdr.Finalize(rec);

        double q2 = reconstruct_q2(rec);
        sum_q2 += q2;
        sum_q2_sq += q2 * q2;

        double sampled_mK = rec.secondary_masses[0];
        int b = (int)((q2 - q2lo) / bw);
        if (b < 0) b = 0;
        if (b >= NB) b = NB - 1;
        if (std::abs(sampled_mK - mK) < std::abs(sampled_mK - mKstar)) countK[b]++;
        else countKstar[b]++;
    }

    double mean_sampled = sum_q2 / Nsamp;
    double var_sampled = sum_q2_sq / Nsamp - mean_sampled * mean_sampled;
    double stderr_mean = std::sqrt(var_sampled / Nsamp);

    // --- (1) Normalization: integral of FinalStateProbability over q^2,
    //         summed over the K and K* mixture, must equal 1. ---
    auto fsp_integral = [&](double comp_mK) -> double {
        std::function<double(double)> integrand = [&](double q2) -> double {
            InteractionRecord r = make_record_at_q2(sig, mD, comp_mK, ml, q2);
            return decay.FinalStateProbability(r);
        };
        double a = ml * ml;
        double b = (mD - comp_mK) * (mD - comp_mK);
        return rombergIntegrate(integrand, a, b);
    };
    double normK = fsp_integral(mK);
    double normKstar = fsp_integral(mKstar);
    double total_norm = normK + normKstar;
    EXPECT_NEAR(total_norm, 1.0, 1e-2);

    // --- (2) Closure of the mean: mean_density (in-test quadrature of
    //         FinalStateProbability) must match mean_sampled within MC error. ---
    auto fsp_q2_integral = [&](double comp_mK) -> double {
        std::function<double(double)> integrand = [&](double q2) -> double {
            InteractionRecord r = make_record_at_q2(sig, mD, comp_mK, ml, q2);
            return q2 * decay.FinalStateProbability(r);
        };
        double a = ml * ml;
        double b = (mD - comp_mK) * (mD - comp_mK);
        return rombergIntegrate(integrand, a, b);
    };
    double mean_density = (fsp_q2_integral(mK) + fsp_q2_integral(mKstar)) / total_norm;
    EXPECT_NEAR(mean_sampled, mean_density, 4.0 * stderr_mean);

    // --- (3) Bin-by-bin closure per sub-population. For each filled bin the
    //         empirical density (count/(N*bw), with the sub-population folded
    //         in via N = total) must match FinalStateProbability at the bin
    //         center within ~3 sigma Poisson error. ---
    for (int b = 0; b < NB; ++b) {
        double q2c = q2lo + (b + 0.5) * bw;
        // K sub-population
        if (countK[b] > 30 && q2c < (mD - mK) * (mD - mK)) {
            double emp = (double)countK[b] / (Nsamp * bw);
            double sigma = std::sqrt((double)countK[b]) / (Nsamp * bw);
            InteractionRecord r = make_record_at_q2(sig, mD, mK, ml, q2c);
            double pred = decay.FinalStateProbability(r);
            EXPECT_NEAR(emp, pred, 4.0 * sigma + 0.02 * pred);
        }
        // K* sub-population
        if (countKstar[b] > 30 && q2c < (mD - mKstar) * (mD - mKstar)) {
            double emp = (double)countKstar[b] / (Nsamp * bw);
            double sigma = std::sqrt((double)countKstar[b]) / (Nsamp * bw);
            InteractionRecord r = make_record_at_q2(sig, mD, mKstar, ml, q2c);
            double pred = decay.FinalStateProbability(r);
            EXPECT_NEAR(emp, pred, 4.0 * sigma + 0.02 * pred);
        }
    }
}

// --- Test 4: TotalDecayWidthForFinalState fails loudly on bad signatures --

TEST(CharmMesonDecay, UnsupportedSignaturesThrow) {
    CharmMesonDecay decay(ParticleType::D0);

    // Positive control: a valid D0 signature has positive width.
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::D0);
    InteractionRecord good;
    good.signature = sigs[0];
    double w = 0.0;
    EXPECT_NO_THROW(w = decay.TotalDecayWidthForFinalState(good));
    EXPECT_GT(w, 0.0);

    // Unsupported primary type.
    InteractionRecord bad_primary;
    bad_primary.signature.primary_type = ParticleType::PiPlus;
    bad_primary.signature.target_type = ParticleType::Decay;
    bad_primary.signature.secondary_types = {ParticleType::Hadrons};
    EXPECT_THROW(decay.TotalDecayWidthForFinalState(bad_primary), std::runtime_error);

    // Matched primary (D0) but a signature with no implemented mode.
    InteractionRecord bad_secondaries;
    bad_secondaries.signature.primary_type = ParticleType::D0;
    bad_secondaries.signature.target_type = ParticleType::Decay;
    bad_secondaries.signature.secondary_types = {ParticleType::PiPlus, ParticleType::PiMinus};
    EXPECT_THROW(decay.TotalDecayWidthForFinalState(bad_secondaries), std::runtime_error);
}

// --- Test 5: analytic angle-average matches a numeric quadrature oracle -----
//
// The weighting code (FinalStateProbability via SampledQ2Density) uses the
// closed-form CharmMesonDecay::VAWeightAngleAverage. The closed form must
// reproduce a high-resolution numeric quadrature of the identical clamped V-A
// weight, evaluated with the same rk::P4 boosts SampleFinalState uses.
namespace {
double numericVAWeightAngleAverage(double mD, double mK, double ml, double m23) {
    double mnu = 0.0;
    double p1Abs = 0.5 * std::sqrt((mD - mK - m23) * (mD + mK + m23)
                                 * (mD + mK - m23) * (mD - mK + m23)) / mD;
    double p23Abs = 0.5 * std::sqrt((m23 - ml - mnu) * (m23 + ml + mnu)
                                  * (m23 + ml - mnu) * (m23 - ml + mnu)) / m23;
    if (p1Abs <= 0.0 || p23Abs <= 0.0) return 0.0;
    double wtMEmax = std::min(std::pow(mD, 4) / 16.0,
                              mD * (mD - mK - ml) * (mD - mK - mnu) * (mD - ml - mnu));
    rk::P4 p4m23(geom3::Vector3(0.0, 0.0, -p1Abs), m23);
    rk::Boost boost = p4m23.labBoost();
    rk::P4 p4K(geom3::Vector3(0.0, 0.0, p1Abs), mK);
    const int N = 20000;            // even -> composite Simpson
    double h = 2.0 / N, sum = 0.0;
    for (int i = 0; i <= N; ++i) {
        double c = -1.0 + i * h;
        double s = std::sqrt(std::max(0.0, 1.0 - c * c));
        geom3::Vector3 dir(s, 0.0, c);
        rk::P4 p4l = rk::P4(p23Abs * dir, ml).boost(boost);
        rk::P4 p4nu = rk::P4(-p23Abs * dir, mnu).boost(boost);
        double w = mD * p4l.e() * p4nu.dot(p4K);
        if (w < 0.0) w = 0.0;
        if (w > wtMEmax) w = wtMEmax;
        double wgt = (i == 0 || i == N) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += wgt * w;
    }
    return 0.5 * (sum * h / 3.0);
}
} // namespace

TEST(CharmMesonDecay, VAWeightAngleAverageMatchesNumericReference) {
    CharmMesonDecay d0(ParticleType::D0), dp(ParticleType::DPlus);
    struct Case { CharmMesonDecay* dec; double mD; double mK; double ml; };
    std::vector<Case> cases = {
        {&d0, Constants::D0Mass,    Constants::KMinusMass,    Constants::electronMass},
        {&d0, Constants::D0Mass,    Constants::KMinusMass,    Constants::muonMass},
        {&d0, Constants::D0Mass,    Constants::KPrimePlusMass, Constants::electronMass},
        {&dp, Constants::DPlusMass, Constants::K0Mass,        Constants::electronMass},
        {&dp, Constants::DPlusMass, Constants::K0Mass,        Constants::muonMass},
        {&dp, Constants::DPlusMass, Constants::KPrimePlusMass, Constants::muonMass},
    };
    for (auto & cs : cases) {
        double m23Min = cs.ml;
        double m23Max = cs.mD - cs.mK;
        const int NG = 40;
        for (int g = 1; g < NG; ++g) {
            double m23 = m23Min + (m23Max - m23Min) * (double)g / NG;
            double ana = cs.dec->VAWeightAngleAverage(cs.mD, cs.mK, cs.ml, m23);
            double num = numericVAWeightAngleAverage(cs.mD, cs.mK, cs.ml, m23);
            // Tolerance is quadrature-limited (Simpson degrades to O(h^2) at the
            // clamp kinks); 1e-4 relative is far tighter than any real algebra bug.
            double tol = 1e-4 * std::abs(num) + 1e-12;
            EXPECT_NEAR(ana, num, tol)
                << "mD=" << cs.mD << " mK=" << cs.mK << " ml=" << cs.ml << " m23=" << m23;
        }
    }
}

// --- Test 6: lab decay length L = beta*gamma*c*tau (cascade separation) ------
//
// The multi-cascade search reconstructs the separation L between the production
// cascade and the charm-decay cascade; the hard regime is below ~10 m. L is set
// by the D species proper lifetime and the lab boost, so it must be physically
// correct across the analysis energy band (TeV-PeV). Decay::TotalDecayLength
// returns beta*gamma*(1/Gamma)*hbarc, and because the modeled branching ratios
// sum to 1 the width recovers the PDG lifetime exactly.
namespace {
constexpr double tau_D0_s = 410.1e-15;   // proper lifetimes hardcoded in
constexpr double tau_Dp_s = 1040.0e-15;  // CharmMesonDecay (seconds).
constexpr double tau_Ds_s = 504.0e-15;
constexpr double c_m_per_s  = 2.99792458e8;    // speed of light [m/s]
// Note: Decay::TotalDecayLength returns METERS (SIREN base length unit; hbarc
// carries the cm->m conversion), consistent with the Taupede reco frame.

InteractionRecord make_boosted_D(ParticleType d, double mD, double E_D) {
    double p = std::sqrt(E_D * E_D - mD * mD);
    InteractionRecord rec;
    rec.signature.primary_type = d;
    rec.signature.target_type = ParticleType::Decay;
    rec.primary_mass = mD;
    rec.primary_momentum = {E_D, 0.0, 0.0, p};   // boosted along +z
    return rec;
}
} // namespace

TEST(CharmMesonDecay, LabDecayLengthIsBetaGammaCTau) {
    struct Case { ParticleType d; double mD; double tau; };
    std::vector<Case> cases = {
        {ParticleType::D0,     Constants::D0Mass,     tau_D0_s},
        {ParticleType::DPlus,  Constants::DPlusMass,  tau_Dp_s},
        {ParticleType::DsPlus, Constants::DsPlusMass, tau_Ds_s},
    };
    for (auto const & cs : cases) {
        CharmMesonDecay decay(cs.d);
        double prev_L = -1.0, prev_E = -1.0;
        for (double E_D : {1.0e4, 1.0e5, 1.0e6}) {   // 10 TeV, 100 TeV, 1 PeV
            InteractionRecord rec = make_boosted_D(cs.d, cs.mD, E_D);
            double L_m = decay.TotalDecayLength(rec);                  // SIREN [m]
            double p = std::sqrt(E_D * E_D - cs.mD * cs.mD);
            double betagamma = p / cs.mD;
            double expected_m = betagamma * c_m_per_s * cs.tau;        // physics truth [m]
            // 0.5% absorbs the ~0.01% rounding of SIREN's hbar/hbarc constants.
            EXPECT_NEAR(L_m, expected_m, 5e-3 * expected_m)
                << "species=" << static_cast<int>(cs.d) << " E_D=" << E_D;
            // beta*gamma is linear in E for E >> m: L(10x E) ~ 10x L.
            if (prev_L > 0.0)
                EXPECT_NEAR(L_m / prev_L, E_D / prev_E, 1e-3 * (E_D / prev_E));
            prev_L = L_m; prev_E = E_D;
        }
    }
}

TEST(CharmMesonDecay, LabDecayLengthSpeciesOrdering) {
    // At fixed boost energy the L ordering follows the lifetimes: D+ > Ds > D0.
    double E_D = 1.0e5;   // 100 TeV
    CharmMesonDecay d0(ParticleType::D0), dp(ParticleType::DPlus), ds(ParticleType::DsPlus);
    double L_D0 = d0.TotalDecayLength(make_boosted_D(ParticleType::D0,     Constants::D0Mass,     E_D));
    double L_Dp = dp.TotalDecayLength(make_boosted_D(ParticleType::DPlus,  Constants::DPlusMass,  E_D));
    double L_Ds = ds.TotalDecayLength(make_boosted_D(ParticleType::DsPlus, Constants::DsPlusMass, E_D));
    EXPECT_GT(L_Dp, L_Ds);
    EXPECT_GT(L_Ds, L_D0);
    // At 100 TeV a D0 travels a few meters -- squarely in the Taupede regime.
    EXPECT_GT(L_D0, 1.0);      // > 1 m
    EXPECT_LT(L_D0, 50.0);     // < 50 m
}

TEST(CharmMesonDecay, FinalStateProbabilityThrowsOnEmptySignature) {
    // finalize() does not copy the signature; FinalStateProbability must reject
    // an empty-signature record loudly instead of indexing out of bounds (UB).
    CharmMesonDecay decay(ParticleType::D0);
    InteractionRecord rec;   // default: empty signature, no secondaries
    EXPECT_THROW(decay.FinalStateProbability(rec), std::runtime_error);
}
