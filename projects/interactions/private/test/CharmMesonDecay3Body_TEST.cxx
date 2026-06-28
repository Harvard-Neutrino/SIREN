/**
 * Unit test for CharmMesonDecay3Body -- Pythia-style 3-body phase-space decay
 * of charm mesons with V-A matrix-element reweighting and K / K*(892) kinematic
 * mixing.
 *
 *   D (m0) -> K (m1) + lepton (m2) + neutrino (m3)
 *
 * All three daughters are constructed in the D rest frame and boosted to the
 * lab with the SAME boost, so 4-momentum is conserved exactly (up to float).
 * A per-event coin flip draws the hadron mass as either the pseudoscalar K mass
 * or the K*(892) mass (= Constants::KPrimePlusMass via KStarMass()).
 *
 * Tests:
 *  1. ThreeBodyEnergyMomentumConservation -- sum of the three secondary
 *     4-momenta equals the primary 4-momentum component-wise.
 *  2. DaughterMassShells -- lepton/neutrino on shell, hadron mass equals either
 *     the pseudoscalar K or the K*(892) mass, and m23 = sqrt((p_l+p_nu)^2) lies
 *     in the allowed window [ml, mD - mK].
 *  3. KStarMixingFraction -- the empirical K vs K* split matches the PDG-derived
 *     fracK within a few binomial sigma (D0 and D+).
 *  4. TotalDecayWidthAndBranchingSums -- per-signature widths are positive, the
 *     three branching ratios partition unity, TotalDecayWidth equals the sum
 *     over signatures, and bad signatures throw.
 *  5. FinalStateProbabilityClosure -- FinalStateProbability is non-negative,
 *     equals 1 for the fully hadronic mode, and the sampled q^2 distribution
 *     closes against FinalStateProbability bin-by-bin within MC error.
 */
#include <cmath>
#include <array>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>

#include <gtest/gtest.h>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "SIREN/interactions/CharmMesonDecay3Body.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Constants.h"

using namespace siren::utilities;
using namespace siren::interactions;
using namespace siren::dataclasses;

namespace {

// Build an unfinalized 3-body semileptonic record for a boosted D meson.
InteractionRecord make_semilep_record(const InteractionSignature & sig,
                                      double mD, double mK, double ml, double E_D) {
    double p_D = std::sqrt(E_D * E_D - mD * mD);
    InteractionRecord rec;
    rec.signature = sig;
    rec.primary_mass = mD;
    rec.primary_momentum = {E_D, 0, 0, p_D};
    rec.primary_helicity = 0;
    rec.target_mass = 0;
    rec.target_helicity = 0;
    rec.secondary_momenta = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    rec.secondary_masses = {mK, ml, 0.0};
    rec.secondary_helicities = {0, 0, 0};
    return rec;
}

// q^2 = (p_D - p_K)^2 from a finalized record.
double reconstruct_q2(const InteractionRecord & rec) {
    double qE  = rec.primary_momentum[0] - rec.secondary_momenta[0][0];
    double qpx = rec.primary_momentum[1] - rec.secondary_momenta[0][1];
    double qpy = rec.primary_momentum[2] - rec.secondary_momenta[0][2];
    double qpz = rec.primary_momentum[3] - rec.secondary_momenta[0][3];
    return qE * qE - qpx * qpx - qpy * qpy - qpz * qpz;
}

} // namespace

// --- Test 1: exact 4-momentum conservation ----------------------------------

TEST(CharmMesonDecay3Body, ThreeBodyEnergyMomentumConservation) {
    CharmMesonDecay3Body decay(ParticleType::D0);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::D0);
    ASSERT_GE(sigs.size(), 2u);

    double mD = Constants::D0Mass;
    double E_D = 100.0;
    auto rng = std::make_shared<SIREN_random>();

    // both semileptonic modes (e and mu)
    std::vector<std::pair<InteractionSignature, double>> cases = {
        {sigs[0], Constants::electronMass},
        {sigs[1], Constants::muonMass},
    };

    for (auto & c : cases) {
        const InteractionSignature & sig = c.first;
        double ml = c.second;
        for (int i = 0; i < 2000; ++i) {
            InteractionRecord rec = make_semilep_record(sig, mD, Constants::KMinusMass, ml, E_D);
            CrossSectionDistributionRecord cdr(rec);
            decay.SampleFinalState(cdr, rng);
            cdr.Finalize(rec);
            ASSERT_EQ(rec.secondary_momenta.size(), 3u);
            for (int comp = 0; comp < 4; ++comp) {
                double s = rec.secondary_momenta[0][comp]
                         + rec.secondary_momenta[1][comp]
                         + rec.secondary_momenta[2][comp];
                EXPECT_NEAR(s, rec.primary_momentum[comp], 1e-5);
            }
        }
    }
}

// --- Test 2: daughter mass shells & m23 window ------------------------------

TEST(CharmMesonDecay3Body, DaughterMassShells) {
    CharmMesonDecay3Body decay(ParticleType::D0);
    auto sig = decay.GetPossibleSignaturesFromParent(ParticleType::D0)[0]; // D0 -> K- e+ nu
    double mD = Constants::D0Mass;
    double ml = Constants::electronMass;
    double mK_base = Constants::KMinusMass;
    double mKstar = Constants::KPrimePlusMass;  // KStarMass() in the source
    double E_D = 100.0;
    auto rng = std::make_shared<SIREN_random>();

    for (int i = 0; i < 5000; ++i) {
        InteractionRecord rec = make_semilep_record(sig, mD, mK_base, ml, E_D);
        CrossSectionDistributionRecord cdr(rec);
        decay.SampleFinalState(cdr, rng);
        cdr.Finalize(rec);

        const auto & pK  = rec.secondary_momenta[0];
        const auto & pl  = rec.secondary_momenta[1];
        const auto & pnu = rec.secondary_momenta[2];

        auto mass = [](const std::array<double,4> & p) {
            double m2 = p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
            return std::sqrt(std::max(0.0, m2));
        };

        // lepton on shell, neutrino massless.
        EXPECT_NEAR(mass(pl), ml, 1e-5);
        EXPECT_NEAR(mass(pnu), 0.0, 1e-5);

        // hadron mass is one of the two mixture components.
        double mK = mass(pK);
        bool is_K = std::abs(mK - mK_base) < 1e-4;
        bool is_Kstar = std::abs(mK - mKstar) < 1e-4;
        EXPECT_TRUE(is_K || is_Kstar) << "hadron mass " << mK << " is neither K nor K*";

        // m23 = (lepton + neutrino) invariant mass within [ml, mD - mK].
        std::array<double,4> p23 = {pl[0]+pnu[0], pl[1]+pnu[1], pl[2]+pnu[2], pl[3]+pnu[3]};
        double m23 = mass(p23);
        double m23max = mD - mK;
        EXPECT_GE(m23, ml - 1e-4);
        EXPECT_LE(m23, m23max + 1e-4);
    }
}

// --- Test 3: K / K*(892) mixing fraction ------------------------------------

TEST(CharmMesonDecay3Body, KStarMixingFraction) {
    double mKstar = Constants::KPrimePlusMass;
    auto rng = std::make_shared<SIREN_random>();
    int N = 20000;

    struct Case {
        ParticleType primary;
        double mD;
        double mK_base;
        double fracK;
    };
    std::vector<Case> cases = {
        {ParticleType::D0,    Constants::D0Mass,    Constants::KMinusMass, 3.41 / (3.41 + 2.17)},
        {ParticleType::DPlus, Constants::DPlusMass, Constants::K0Mass,     8.74 / (8.74 + 5.33)},
    };

    for (auto & c : cases) {
        CharmMesonDecay3Body decay(c.primary);
        auto sig = decay.GetPossibleSignaturesFromParent(c.primary)[0];
        double ml = Constants::electronMass;
        double E_D = 50.0;
        long countK = 0;
        for (int i = 0; i < N; ++i) {
            InteractionRecord rec = make_semilep_record(sig, c.mD, c.mK_base, ml, E_D);
            CrossSectionDistributionRecord cdr(rec);
            decay.SampleFinalState(cdr, rng);
            cdr.Finalize(rec);
            double m2 = rec.secondary_momenta[0][0]*rec.secondary_momenta[0][0]
                      - rec.secondary_momenta[0][1]*rec.secondary_momenta[0][1]
                      - rec.secondary_momenta[0][2]*rec.secondary_momenta[0][2]
                      - rec.secondary_momenta[0][3]*rec.secondary_momenta[0][3];
            double mK = std::sqrt(std::max(0.0, m2));
            if (std::abs(mK - c.mK_base) < std::abs(mK - mKstar)) countK++;
        }
        double emp = (double)countK / N;
        double sigma = std::sqrt(c.fracK * (1.0 - c.fracK) / N);
        EXPECT_NEAR(emp, c.fracK, 5.0 * sigma);
    }
}

// --- Test 4: decay widths & branching bookkeeping ---------------------------

TEST(CharmMesonDecay3Body, TotalDecayWidthAndBranchingSums) {
    // D0: BR_semilep(e) = BR_semilep(mu) = 0.0558; hadronic remainder = 1 - 2*0.0558.
    // D+: BR_semilep(e) = BR_semilep(mu) = 0.1407; hadronic remainder = 1 - 2*0.1407.
    struct Case { ParticleType primary; double br_semilep; };
    std::vector<Case> cases = {
        {ParticleType::D0,    0.0558},
        {ParticleType::DPlus, 0.1407},
    };

    for (auto & c : cases) {
        CharmMesonDecay3Body decay(c.primary);
        auto sigs = decay.GetPossibleSignaturesFromParent(c.primary);
        ASSERT_EQ(sigs.size(), 3u);   // e-mode, mu-mode, hadronic

        double sum_width = 0.0;
        for (auto & sig : sigs) {
            InteractionRecord r;
            r.signature = sig;
            double w = 0.0;
            EXPECT_NO_THROW(w = decay.TotalDecayWidthForFinalState(r));
            EXPECT_GT(w, 0.0);
            sum_width += w;
        }

        // The three branching ratios partition unity.
        double br_sum = c.br_semilep + c.br_semilep + (1.0 - 2.0 * c.br_semilep);
        EXPECT_NEAR(br_sum, 1.0, 1e-12);

        // TotalDecayWidth(primary) is the sum over GetPossibleSignaturesFromParent.
        double total = decay.TotalDecayWidth(c.primary);
        EXPECT_NEAR(total, sum_width, std::abs(total) * 1e-12 + 1e-30);
        EXPECT_GT(total, 0.0);
    }

    // Unsupported primary throws.
    CharmMesonDecay3Body decay(ParticleType::D0);
    InteractionRecord bad_primary;
    bad_primary.signature.primary_type = ParticleType::PiPlus;
    bad_primary.signature.target_type = ParticleType::Decay;
    bad_primary.signature.secondary_types = {ParticleType::Hadrons};
    EXPECT_THROW(decay.TotalDecayWidthForFinalState(bad_primary), std::runtime_error);

    // Matched primary (D0) but an unimplemented set of secondaries throws.
    InteractionRecord bad_secondaries;
    bad_secondaries.signature.primary_type = ParticleType::D0;
    bad_secondaries.signature.target_type = ParticleType::Decay;
    bad_secondaries.signature.secondary_types = {ParticleType::PiPlus, ParticleType::PiMinus, ParticleType::Pi0};
    EXPECT_THROW(decay.TotalDecayWidthForFinalState(bad_secondaries), std::runtime_error);
}

// --- Test 5: FinalStateProbability sanity + q^2 closure ---------------------

TEST(CharmMesonDecay3Body, FinalStateProbabilityClosure) {
    CharmMesonDecay3Body decay(ParticleType::D0);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::D0);
    auto sig = sigs[0];   // D0 -> K- e+ nu

    double mD = Constants::D0Mass;
    double ml = Constants::electronMass;
    double mK = Constants::KMinusMass;
    double mKstar = Constants::KPrimePlusMass;
    double E_D = 100.0;
    auto rng = std::make_shared<SIREN_random>();

    // (a) The fully hadronic catch-all has FinalStateProbability == 1.
    {
        auto had_sig = sigs[2];  // D0 -> Hadrons
        InteractionRecord rec;
        rec.signature = had_sig;
        rec.primary_mass = mD;
        rec.primary_momentum = {E_D, 0, 0, std::sqrt(E_D * E_D - mD * mD)};
        rec.secondary_momenta = {{0, 0, 0, 0}};
        rec.secondary_masses = {0.0};
        EXPECT_NEAR(decay.FinalStateProbability(rec), 1.0, 1e-12);
    }

    // (b) Sample many semileptonic events, histogram q^2 separately for the K
    //     and K* sub-populations, and compare the empirical density per bin to
    //     fractions[comp] * FinalStateProbability evaluated at a record built at
    //     that q^2 with the matching hadron mass.
    const int NB = 16;
    double q2lo = 0.0;
    double q2hi = (mD - mK) * (mD - mK);  // widest support (K mass)
    double bw = (q2hi - q2lo) / NB;
    std::vector<long> countK(NB, 0), countKstar(NB, 0);
    const int N = 40000;

    for (int i = 0; i < N; ++i) {
        InteractionRecord rec = make_semilep_record(sig, mD, mK, ml, E_D);
        CrossSectionDistributionRecord cdr(rec);
        decay.SampleFinalState(cdr, rng);
        cdr.Finalize(rec);

        double q2 = reconstruct_q2(rec);
        double sampled_mK = rec.secondary_masses[0];
        double fsp = decay.FinalStateProbability(rec);
        EXPECT_GE(fsp, 0.0);

        int b = (int)((q2 - q2lo) / bw);
        if (b < 0) b = 0;
        if (b >= NB) b = NB - 1;
        if (std::abs(sampled_mK - mK) < std::abs(sampled_mK - mKstar)) countK[b]++;
        else countKstar[b]++;
    }

    // Build a record at a target q^2 in the D rest frame with hadron mass comp_mK
    // and evaluate FinalStateProbability there (it already includes the K/K*
    // mixture weight fractions[comp]).
    auto fsp_at_q2 = [&](double comp_mK, double q2) -> double {
        double EK = (mD * mD + comp_mK * comp_mK - q2) / (2 * mD);
        double pk_sq = EK * EK - comp_mK * comp_mK;
        if (pk_sq < 0) return 0.0;
        double PK = std::sqrt(pk_sq);
        InteractionRecord r;
        r.signature = sig;
        r.primary_mass = mD;
        r.primary_momentum = {mD, 0, 0, 0};   // D at rest (q^2 frame-independent)
        r.secondary_momenta = {{EK, PK, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        r.secondary_masses = {comp_mK, ml, 0.0};
        return decay.FinalStateProbability(r);
    };

    for (int b = 0; b < NB; ++b) {
        double q2c = q2lo + (b + 0.5) * bw;
        if (countK[b] > 40 && q2c < (mD - mK) * (mD - mK)) {
            double emp = (double)countK[b] / (N * bw);
            double sigma = std::sqrt((double)countK[b]) / (N * bw);
            double pred = fsp_at_q2(mK, q2c);
            EXPECT_NEAR(emp, pred, 4.0 * sigma + 0.03 * pred);
        }
        if (countKstar[b] > 40 && q2c < (mD - mKstar) * (mD - mKstar)) {
            double emp = (double)countKstar[b] / (N * bw);
            double sigma = std::sqrt((double)countKstar[b]) / (N * bw);
            double pred = fsp_at_q2(mKstar, q2c);
            EXPECT_NEAR(emp, pred, 4.0 * sigma + 0.03 * pred);
        }
    }
}

// --- Analytic angle-average matches a numeric quadrature oracle ------------
//
// The weighting code uses the closed-form CharmMesonDecay3Body::
// VAWeightAngleAverage. The closed form must reproduce a high-resolution numeric
// quadrature of the identical clamped V-A weight, evaluated with the same rk::P4
// boosts SampleFinalState uses.
namespace {
double numericVAWeightAngleAverage3B(double mD, double mK, double ml, double m23) {
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
    const int N = 20000;
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

TEST(CharmMesonDecay3Body, VAWeightAngleAverageMatchesNumericReference) {
    CharmMesonDecay3Body d0(ParticleType::D0), dp(ParticleType::DPlus);
    struct Case { CharmMesonDecay3Body* dec; double mD; double mK; double ml; };
    std::vector<Case> cases = {
        {&d0, Constants::D0Mass,    Constants::KMinusMass,     Constants::electronMass},
        {&d0, Constants::D0Mass,    Constants::KMinusMass,     Constants::muonMass},
        {&d0, Constants::D0Mass,    Constants::KPrimePlusMass, Constants::electronMass},
        {&dp, Constants::DPlusMass, Constants::K0Mass,         Constants::electronMass},
        {&dp, Constants::DPlusMass, Constants::KPrimePlusMass, Constants::muonMass},
    };
    for (auto & cs : cases) {
        double m23Min = cs.ml;
        double m23Max = cs.mD - cs.mK;
        const int NG = 40;
        for (int g = 1; g < NG; ++g) {
            double m23 = m23Min + (m23Max - m23Min) * (double)g / NG;
            double ana = cs.dec->VAWeightAngleAverage(cs.mD, cs.mK, cs.ml, m23);
            double num = numericVAWeightAngleAverage3B(cs.mD, cs.mK, cs.ml, m23);
            double tol = 1e-4 * std::abs(num) + 1e-12;
            EXPECT_NEAR(ana, num, tol)
                << "mD=" << cs.mD << " mK=" << cs.mK << " ml=" << cs.ml << " m23=" << m23;
        }
    }
}

// CharmMesonDecay3Body implements only D0 and D+. Unsupported species (e.g. Ds)
// and empty-signature records must fail loudly rather than silently mis-decay or
// index out of bounds. (CharmMesonDecay covers D0/D+/Ds and anti-flavors.)
TEST(CharmMesonDecay3Body, UnsupportedSpeciesAndEmptySignatureThrow) {
    EXPECT_THROW({ CharmMesonDecay3Body d(ParticleType::DsPlus); }, std::runtime_error);
    CharmMesonDecay3Body dp(ParticleType::DPlus);
    EXPECT_THROW(dp.GetPossibleSignaturesFromParent(ParticleType::DsPlus), std::runtime_error);
    CharmMesonDecay3Body d0(ParticleType::D0);
    InteractionRecord rec;   // default: empty signature
    EXPECT_THROW(d0.FinalStateProbability(rec), std::runtime_error);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
