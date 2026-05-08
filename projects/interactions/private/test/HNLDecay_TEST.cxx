#include <cmath>
#include <set>
#include <vector>
#include <algorithm>

#include <gtest/gtest.h>

#include "SIREN/interactions/HNLDecay.h"
#include "SIREN/interactions/HNLDipoleDecay.h"
#include "SIREN/interactions/ElectroweakDecay.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Sampling.h"

using namespace siren::interactions;
using namespace siren::dataclasses;
using namespace siren::utilities;

// ---------------------------------------------------------------------------
// HNLDecay
// ---------------------------------------------------------------------------

TEST(HNLDecay, ConstructorVector) {
    std::vector<double> mixing = {0.0, 1e-3, 0.0};
    ASSERT_NO_THROW(HNLDecay decay(0.3, mixing, HNLDecay::Majorana));
}

TEST(HNLDecay, ConstructorScalar) {
    ASSERT_NO_THROW(HNLDecay decay(0.3, 1e-3, HNLDecay::Dirac));
}

TEST(HNLDecay, GetHNLMass) {
    HNLDecay decay(0.5, 1e-3, HNLDecay::Majorana);
    EXPECT_DOUBLE_EQ(decay.GetHNLMass(), 0.5);
}

TEST(HNLDecay, SignaturesNonEmpty) {
    HNLDecay decay(0.3, 1e-3, HNLDecay::Majorana);
    auto sigs = decay.GetPossibleSignatures();
    EXPECT_GT(sigs.size(), 0u);
}

TEST(HNLDecay, SignaturesNoDuplicates) {
    HNLDecay decay(0.3, 1e-3, HNLDecay::Majorana);
    auto sigs = decay.GetPossibleSignatures();
    // Check no duplicate signatures by comparing all pairs
    for (size_t i = 0; i < sigs.size(); ++i) {
        for (size_t j = i + 1; j < sigs.size(); ++j) {
            bool same = (sigs[i].primary_type == sigs[j].primary_type &&
                         sigs[i].target_type == sigs[j].target_type &&
                         sigs[i].secondary_types == sigs[j].secondary_types);
            EXPECT_FALSE(same) << "Duplicate signature at indices " << i << " and " << j;
        }
    }
}

TEST(HNLDecay, SignaturesFromParent) {
    HNLDecay decay(0.3, 1e-3, HNLDecay::Majorana);
    auto sigs_n4 = decay.GetPossibleSignaturesFromParent(ParticleType::N4);
    auto sigs_n4bar = decay.GetPossibleSignaturesFromParent(ParticleType::N4Bar);
    EXPECT_GT(sigs_n4.size(), 0u);
    EXPECT_GT(sigs_n4bar.size(), 0u);
    // All N4 signatures should have N4 as primary
    for (auto const & sig : sigs_n4) {
        EXPECT_EQ(sig.primary_type, ParticleType::N4);
    }
}

TEST(HNLDecay, TotalDecayWidthPositive) {
    HNLDecay decay(0.3, 1e-3, HNLDecay::Majorana);
    double width = decay.TotalDecayWidth(ParticleType::N4);
    EXPECT_GT(width, 0.0);
}

TEST(HNLDecay, HeavierHNLDecaysFaster) {
    HNLDecay light(0.2, 1e-3, HNLDecay::Majorana);
    HNLDecay heavy(0.5, 1e-3, HNLDecay::Majorana);
    double width_light = light.TotalDecayWidth(ParticleType::N4);
    double width_heavy = heavy.TotalDecayWidth(ParticleType::N4);
    EXPECT_GT(width_heavy, width_light);
}

TEST(HNLDecay, LargerMixingDecaysFaster) {
    HNLDecay small_mix(0.3, 1e-4, HNLDecay::Majorana);
    HNLDecay large_mix(0.3, 1e-3, HNLDecay::Majorana);
    double width_small = small_mix.TotalDecayWidth(ParticleType::N4);
    double width_large = large_mix.TotalDecayWidth(ParticleType::N4);
    EXPECT_GT(width_large, width_small);
}

TEST(HNLDecay, MajoranaVsDiracWidths) {
    // Majorana HNL has additional decay channels
    HNLDecay dirac(0.3, 1e-3, HNLDecay::Dirac);
    HNLDecay majorana(0.3, 1e-3, HNLDecay::Majorana);
    double width_dirac = dirac.TotalDecayWidth(ParticleType::N4);
    double width_majorana = majorana.TotalDecayWidth(ParticleType::N4);
    EXPECT_GT(width_majorana, width_dirac);
}

TEST(HNLDecay, BelowPionThresholdStillDecays) {
    // HNL lighter than pion should still decay via 3-body leptonic channels
    HNLDecay decay(0.05, 1e-3, HNLDecay::Majorana);
    double width = decay.TotalDecayWidth(ParticleType::N4);
    EXPECT_GT(width, 0.0);
}

// ---------------------------------------------------------------------------
// HNLDipoleDecay
// ---------------------------------------------------------------------------

TEST(HNLDipoleDecay, Constructor) {
    std::vector<double> coupling = {0.0, 1e-6, 0.0};
    ASSERT_NO_THROW(HNLDipoleDecay decay(0.3, coupling, HNLDipoleDecay::Dirac));
}

TEST(HNLDipoleDecay, SignaturesNonEmpty) {
    HNLDipoleDecay decay(0.3, 1e-6, HNLDipoleDecay::Dirac);
    auto sigs = decay.GetPossibleSignatures();
    EXPECT_GT(sigs.size(), 0u);
}

TEST(HNLDipoleDecay, TotalDecayWidthPositive) {
    HNLDipoleDecay decay(0.3, 1e-6, HNLDipoleDecay::Dirac);
    double width = decay.TotalDecayWidth(ParticleType::NuF4);
    EXPECT_GT(width, 0.0);
}

// ---------------------------------------------------------------------------
// ElectroweakDecay
// ---------------------------------------------------------------------------

TEST(ElectroweakDecay, ConstructorW) {
    std::set<ParticleType> primaries = {ParticleType::WPlus, ParticleType::WMinus};
    ASSERT_NO_THROW(ElectroweakDecay decay(primaries));
}

TEST(ElectroweakDecay, ConstructorZ) {
    std::set<ParticleType> primaries = {ParticleType::Z0};
    ASSERT_NO_THROW(ElectroweakDecay decay(primaries));
}

TEST(ElectroweakDecay, SignaturesNonEmpty) {
    std::set<ParticleType> primaries = {ParticleType::WPlus, ParticleType::WMinus, ParticleType::Z0};
    ElectroweakDecay decay(primaries);
    auto sigs = decay.GetPossibleSignatures();
    EXPECT_GT(sigs.size(), 0u);
}

TEST(ElectroweakDecay, WDecayWidthPositive) {
    std::set<ParticleType> primaries = {ParticleType::WPlus};
    ElectroweakDecay decay(primaries);
    double width = decay.TotalDecayWidth(ParticleType::WPlus);
    EXPECT_GT(width, 0.0);
}

TEST(ElectroweakDecay, ZDecayWidthPositive) {
    std::set<ParticleType> primaries = {ParticleType::Z0};
    ElectroweakDecay decay(primaries);
    double width = decay.TotalDecayWidth(ParticleType::Z0);
    EXPECT_GT(width, 0.0);
}

// ---------------------------------------------------------------------------
// SampleFinalState and momentum conservation
// ---------------------------------------------------------------------------

TEST(HNLDecay, SampleTwoBodyConservesMomentum) {
    // N4 (0.5 GeV) -> nu + pi0 (2-body decay)
    double hnl_mass = 0.5;
    HNLDecay decay(hnl_mass, std::vector<double>{0, 0, 1e-3}, HNLDecay::ChiralNature::Dirac);

    InteractionRecord event;
    event.signature.primary_type = ParticleType::N4;
    event.signature.target_type = ParticleType::Decay;
    event.signature.secondary_types = {ParticleType::NuTau, ParticleType::Pi0};
    event.primary_mass = hnl_mass;
    double energy = 5.0;
    event.primary_momentum = {energy, 0, 0, std::sqrt(energy*energy - hnl_mass*hnl_mass)};
    event.primary_helicity = 1.0;

    auto rand = std::make_shared<SIREN_random>();
    for (int i = 0; i < 100; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);

        // Sum secondary momenta
        double psum[4] = {0, 0, 0, 0};
        for (size_t s = 0; s < event.secondary_momenta.size(); ++s) {
            for (int j = 0; j < 4; ++j)
                psum[j] += event.secondary_momenta[s][j];
        }
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(psum[j], event.primary_momentum[j],
                        1e-6 * std::max(1.0, std::abs(event.primary_momentum[j])))
                << "4-momentum not conserved in component " << j;
        }
    }
}

TEST(HNLDecay, SampleThreeBodyConservesMomentum) {
    // N4 (0.5 GeV) -> nu + e- + e+ (3-body leptonic)
    double hnl_mass = 0.5;
    HNLDecay decay(hnl_mass, std::vector<double>{0, 0, 1e-3}, HNLDecay::ChiralNature::Dirac);

    InteractionRecord event;
    event.signature.primary_type = ParticleType::N4;
    event.signature.target_type = ParticleType::Decay;
    event.signature.secondary_types = {ParticleType::NuLight, ParticleType::EMinus, ParticleType::EPlus};
    event.primary_mass = hnl_mass;
    double energy = 5.0;
    event.primary_momentum = {energy, 0, 0, std::sqrt(energy*energy - hnl_mass*hnl_mass)};
    event.primary_helicity = 1.0;

    auto rand = std::make_shared<SIREN_random>();
    for (int i = 0; i < 50; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);

        double psum[4] = {0, 0, 0, 0};
        for (size_t s = 0; s < event.secondary_momenta.size(); ++s) {
            for (int j = 0; j < 4; ++j)
                psum[j] += event.secondary_momenta[s][j];
        }
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(psum[j], event.primary_momentum[j],
                        1e-6 * std::max(1.0, std::abs(event.primary_momentum[j])))
                << "4-momentum not conserved in component " << j;
        }
    }
}

TEST(HNLDecay, SampleDifferentialWidthPositive) {
    // DifferentialDecayWidth should be > 0 for a valid 2-body signature
    double hnl_mass = 0.5;
    HNLDecay decay(hnl_mass, std::vector<double>{0, 0, 1e-3}, HNLDecay::ChiralNature::Dirac);

    InteractionRecord event;
    event.signature.primary_type = ParticleType::N4;
    event.signature.target_type = ParticleType::Decay;
    event.signature.secondary_types = {ParticleType::NuTau, ParticleType::Pi0};
    event.primary_mass = hnl_mass;
    double energy = 5.0;
    event.primary_momentum = {energy, 0, 0, std::sqrt(energy*energy - hnl_mass*hnl_mass)};
    event.primary_helicity = 1.0;

    // Need to populate secondary momenta for DifferentialDecayWidth
    auto rand = std::make_shared<SIREN_random>();
    CrossSectionDistributionRecord xsec_record(event);
    decay.SampleFinalState(xsec_record, rand);
    xsec_record.Finalize(event);

    double dw = decay.DifferentialDecayWidth(event);
    EXPECT_GT(dw, 0.0);

    double fp = decay.FinalStateProbability(event);
    EXPECT_GT(fp, 0.0);
    EXPECT_LE(fp, 1.0);
}

// ---------------------------------------------------------------------------
// Sampling utility
// ---------------------------------------------------------------------------

TEST(Sampling, UniformProposal) {
    auto random = std::make_shared<SIREN_random>();
    auto proposal = [&]() -> std::pair<std::vector<double>, double> {
        return std::make_pair(std::vector<double>{random->Uniform(-5, 5)}, 1.0);
    };
    auto likelihood = [](std::vector<double> const & x) -> double {
        return std::exp(-x[0] * x[0] / 2.0);
    };
    auto sample = MetropolisHasting_Sample(proposal, likelihood, random);
    ASSERT_EQ(sample.size(), 1u);
    EXPECT_TRUE(std::isfinite(sample[0]));
}

TEST(Sampling, ZeroLikelihoodDoesNotCrash) {
    auto random = std::make_shared<SIREN_random>();
    auto proposal = [&]() -> std::pair<std::vector<double>, double> {
        return std::make_pair(std::vector<double>{random->Uniform(0, 1)}, 1.0);
    };
    auto likelihood = [](std::vector<double> const & x) -> double {
        return 0.0;
    };
    ASSERT_NO_THROW(MetropolisHasting_Sample(proposal, likelihood, random));
}

TEST(Sampling, NonUniformProposal) {
    // Proposal N(0,2), target N(2,0.5) -- samples should converge to target mean
    auto random = std::make_shared<SIREN_random>();
    auto proposal = [&]() -> std::pair<std::vector<double>, double> {
        double u1 = random->Uniform(0, 1);
        double u2 = random->Uniform(0, 1);
        double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2 * M_PI * u2);
        double x = 2.0 * z;
        double density = std::exp(-x * x / 8.0) / (2.0 * std::sqrt(2.0 * M_PI));
        return std::make_pair(std::vector<double>{x}, density);
    };
    auto target = [](std::vector<double> const & x) -> double {
        double dx = x[0] - 2.0;
        return std::exp(-dx * dx / 0.5) / (0.5 * std::sqrt(2.0 * M_PI));
    };

    double sum = 0;
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        auto sample = MetropolisHasting_Sample(proposal, target, random, 200);
        sum += sample[0];
    }
    double mean = sum / n;
    // Mean should be near 2 (target), not 0 (proposal)
    EXPECT_NEAR(mean, 2.0, 0.3);
}

// ---------------------------------------------------------------------------
// Physics regression tests
// ---------------------------------------------------------------------------

TEST(HNLDecay, MesonChannelsZeroBelowThreshold) {
    // Below pion mass, all hadronic two-body channels should have zero width.
    HNLDecay decay(0.05, 1e-3, HNLDecay::Majorana);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::N4);
    InteractionRecord record;
    record.primary_mass = 0.05;
    for (auto const & sig : sigs) {
        if (sig.secondary_types.size() != 2) continue;
        ParticleType s1 = sig.secondary_types[1];
        // Look for any charged or neutral meson final state
        bool is_meson = (s1 == ParticleType::PiMinus || s1 == ParticleType::PiPlus ||
                        s1 == ParticleType::Pi0 || s1 == ParticleType::KMinus ||
                        s1 == ParticleType::KPlus || s1 == ParticleType::Eta);
        if (!is_meson) continue;
        record.signature = sig;
        EXPECT_EQ(decay.TotalDecayWidthForFinalState(record), 0.0)
            << "Sub-threshold meson channel returned non-zero width for secondary " << static_cast<int>(s1);
    }
}

TEST(HNLDecay, AllPartialWidthsNonNegative) {
    HNLDecay decay(0.5, 1e-3, HNLDecay::Majorana);
    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::N4);
    InteractionRecord record;
    record.primary_mass = 0.5;
    for (auto const & sig : sigs) {
        record.signature = sig;
        double w = decay.TotalDecayWidthForFinalState(record);
        EXPECT_GE(w, 0.0);
        EXPECT_TRUE(std::isfinite(w));
    }
}

TEST(HNLDecay, MajoranaApproxTwiceDirac) {
    // For a pure leptonic mixing, total Majorana width is ~2x Dirac width
    // (Majorana adds the charge-conjugate channels)
    HNLDecay dirac(0.5, std::vector<double>{0, 0, 1e-3}, HNLDecay::Dirac);
    HNLDecay majorana(0.5, std::vector<double>{0, 0, 1e-3}, HNLDecay::Majorana);
    double wD = dirac.TotalDecayWidth(ParticleType::N4);
    double wM = majorana.TotalDecayWidth(ParticleType::N4);
    EXPECT_NEAR(wM / wD, 2.0, 0.1);
}

TEST(HNLDecay, ThreeBodyOnShell) {
    // Sample many 3-body decays and verify each secondary is on its mass shell
    double hnl_mass = 0.5;
    HNLDecay decay(hnl_mass, std::vector<double>{0, 0, 1e-3}, HNLDecay::Dirac);

    InteractionRecord event;
    event.signature.primary_type = ParticleType::N4;
    event.signature.target_type = ParticleType::Decay;
    event.signature.secondary_types = {ParticleType::NuLight, ParticleType::EMinus, ParticleType::EPlus};
    event.primary_mass = hnl_mass;
    double energy = 5.0;
    event.primary_momentum = {energy, 0, 0, std::sqrt(energy*energy - hnl_mass*hnl_mass)};
    event.primary_helicity = 1.0;

    auto rand = std::make_shared<SIREN_random>();
    double me = 0.000511; // electron mass GeV
    for (int i = 0; i < 50; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);
        ASSERT_EQ(event.secondary_momenta.size(), 3u);
        for (size_t s = 0; s < 3; ++s) {
            double E = event.secondary_momenta[s][0];
            double px = event.secondary_momenta[s][1];
            double py = event.secondary_momenta[s][2];
            double pz = event.secondary_momenta[s][3];
            double m2 = E*E - px*px - py*py - pz*pz;
            double expected = (s == 0) ? 0.0 : me*me;
            EXPECT_NEAR(m2, expected, 1e-4)
                << "Secondary " << s << " not on mass shell";
        }
    }
}

TEST(ElectroweakDecay, WHadronicToLeptonicRatio) {
    // SM prediction: sum_q Gamma(W -> q qbar) / sum_l Gamma(W -> l nu) ~= 3 * sum |V_qq|^2 ~= 2
    // (sum over kinematically accessible quarks: ud, us, cd, cs)
    std::set<ParticleType> primaries = {ParticleType::WPlus};
    ElectroweakDecay decay(primaries);

    auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::WPlus);
    double w_lep = 0, w_had = 0;
    InteractionRecord record;
    record.primary_mass = siren::utilities::Constants::wMass;
    for (auto const & sig : sigs) {
        record.signature = sig;
        double w = decay.TotalDecayWidthForFinalState(record);
        ParticleType s0 = sig.secondary_types[0];
        bool is_lep = (s0 == ParticleType::EPlus || s0 == ParticleType::MuPlus || s0 == ParticleType::TauPlus);
        if (is_lep) w_lep += w;
        else w_had += w;
    }
    EXPECT_GT(w_lep, 0.0);
    EXPECT_GT(w_had, 0.0);
    // Hadronic / leptonic ratio should be ~2 (within ~10% of SM)
    EXPECT_NEAR(w_had / w_lep, 2.0, 0.3);
}

TEST(ElectroweakDecay, ZHadronicToLeptonicColorFactor) {
    // For a single quark flavor, Gamma(Z->qqbar) should include factor 3 from color
    std::set<ParticleType> primaries = {ParticleType::Z0};
    ElectroweakDecay decay(primaries);

    InteractionRecord record_d;
    record_d.primary_mass = siren::utilities::Constants::zMass;
    record_d.signature.primary_type = ParticleType::Z0;
    record_d.signature.secondary_types = {ParticleType::d, ParticleType::dBar};

    InteractionRecord record_l;
    record_l.primary_mass = siren::utilities::Constants::zMass;
    record_l.signature.primary_type = ParticleType::Z0;
    record_l.signature.secondary_types = {ParticleType::EMinus, ParticleType::EPlus};

    double w_d = decay.TotalDecayWidthForFinalState(record_d);
    double w_l = decay.TotalDecayWidthForFinalState(record_l);
    EXPECT_GT(w_d, 0.0);
    EXPECT_GT(w_l, 0.0);
    // Width ratio depends on cV,cA; the color factor of 3 should make w_d > w_l
    EXPECT_GT(w_d, w_l);
}

TEST(ElectroweakDecay, ZMismatchedFinalStateZero) {
    // Non-pair final states (e.g. e- e-, e- mu+) should yield zero width
    std::set<ParticleType> primaries = {ParticleType::Z0};
    ElectroweakDecay decay(primaries);

    InteractionRecord record;
    record.primary_mass = siren::utilities::Constants::zMass;
    record.signature.primary_type = ParticleType::Z0;
    record.signature.secondary_types = {ParticleType::EMinus, ParticleType::MuPlus};
    EXPECT_EQ(decay.TotalDecayWidthForFinalState(record), 0.0);

    record.signature.secondary_types = {ParticleType::d, ParticleType::sBar};
    EXPECT_EQ(decay.TotalDecayWidthForFinalState(record), 0.0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
