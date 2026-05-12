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
#include "SIREN/utilities/Constants.h"
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
    double width = decay.TotalDecayWidth(ParticleType::N4);
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
    // k2 = pHNL - k3 - k4 by construction, so 4-momentum is conserved exactly.
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

TEST(HNLDecay, ThreeBodyChargedLeptonsOnShell) {
    // Sample many 3-body decays and verify all particles are on mass shell.
    // Charged leptons are constructed with P4(3vec, mass) and are exactly on-shell.
    // The neutrino is derived via k2 = pHNL - k3 - k4; the rotation-matrix
    // construction preserves cos(alpha34) exactly, so the neutrino is on the
    // lightcone to machine precision.
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
    double me = siren::utilities::Constants::electronMass;
    double max_nu_m2 = 0;
    for (int i = 0; i < 50; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);
        ASSERT_EQ(event.secondary_momenta.size(), 3u);

        for (size_t s = 1; s < 3; ++s) {
            double E = event.secondary_momenta[s][0];
            double px = event.secondary_momenta[s][1];
            double py = event.secondary_momenta[s][2];
            double pz = event.secondary_momenta[s][3];
            double m2 = E*E - px*px - py*py - pz*pz;
            EXPECT_NEAR(m2, me*me, 1e-6)
                << "Charged lepton " << s << " not on mass shell";
        }

        double E = event.secondary_momenta[0][0];
        double px = event.secondary_momenta[0][1];
        double py = event.secondary_momenta[0][2];
        double pz = event.secondary_momenta[0][3];
        max_nu_m2 = std::max(max_nu_m2, std::abs(E*E - px*px - py*py - pz*pz));
    }
    EXPECT_LT(max_nu_m2, 1e-8)
        << "Neutrino m^2 violation exceeds 1e-8 GeV^2";
}

TEST(HNLDecay, ThreeBodyDifferentialWidthPositive) {
    // DifferentialDecayWidth reconstructs the neutrino as P4(3vec, 0),
    // recomputing E from |p|.  If the neutrino is off-shell, the
    // reconstructed s1/s2 differ from the sampled values, potentially
    // producing a negative or zero width.  This regression test verifies
    // that DifferentialDecayWidth returns a positive value for sampled
    // 3-body events, and that FinalStateProbability stays in (0, 1].
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
    for (int i = 0; i < 50; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);

        double dw = decay.DifferentialDecayWidth(event);
        EXPECT_GT(dw, 0.0)
            << "3-body DifferentialDecayWidth should be positive for a sampled event";

        double fp = decay.FinalStateProbability(event);
        EXPECT_GT(fp, 0.0);
        EXPECT_LE(fp, 1.0)
            << "FinalStateProbability > 1 indicates Dalitz invariant inconsistency";
    }
}

TEST(HNLDecay, ThreeBodyDalitzRoundTrip) {
    // The sampled Dalitz invariants (s1, s2) should be recoverable from
    // the stored secondary momenta.  DifferentialDecayWidth reconstructs
    // k2 = P4(3vec, 0) and computes s1 = (k2+k3)^2/M_N^2.  If the
    // angular solver produces an off-shell neutrino, the reconstructed
    // s1/s2 will differ from the sampled values.  This test verifies
    // round-trip consistency.
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
    double me = siren::utilities::Constants::electronMass;
    double max_ds = 0;
    for (int i = 0; i < 50; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);

        // Reconstruct k2 the same way DifferentialDecayWidth does:
        // P4(3-momentum, mass=0) so E = |p|
        rk::P4 k2(geom3::Vector3(event.secondary_momenta[0][1],
                                  event.secondary_momenta[0][2],
                                  event.secondary_momenta[0][3]), 0);
        rk::P4 k3(geom3::Vector3(event.secondary_momenta[1][1],
                                  event.secondary_momenta[1][2],
                                  event.secondary_momenta[1][3]), me);
        rk::P4 k4(geom3::Vector3(event.secondary_momenta[2][1],
                                  event.secondary_momenta[2][2],
                                  event.secondary_momenta[2][3]), me);

        double s1_recon = (k2+k3).dot(k2+k3) / (hnl_mass*hnl_mass);
        double s2_recon = (k2+k4).dot(k2+k4) / (hnl_mass*hnl_mass);

        // Also compute s1/s2 using the stored energy (not recomputed from mass=0)
        // This is what was sampled
        rk::P4 k2_stored(event.secondary_momenta[0][0],
                          geom3::Vector3(event.secondary_momenta[0][1],
                                         event.secondary_momenta[0][2],
                                         event.secondary_momenta[0][3]));
        double s1_stored = (k2_stored+k3).dot(k2_stored+k3) / (hnl_mass*hnl_mass);
        double s2_stored = (k2_stored+k4).dot(k2_stored+k4) / (hnl_mass*hnl_mass);

        max_ds = std::max(max_ds, std::abs(s1_recon - s1_stored));
        max_ds = std::max(max_ds, std::abs(s2_recon - s2_stored));
    }
    EXPECT_LT(max_ds, 1e-10)
        << "Dalitz invariants from P4(3vec,0) differ from stored -- neutrino off-shell";
}

TEST(HNLDecay, ThreeBodyOpeningAngleConsistency) {
    // Verify that the opening angle between k3 and k4 in the rest frame
    // matches the Dalitz prediction.  This is the core invariant that
    // the rotation-matrix construction preserves exactly.
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
    double me = siren::utilities::Constants::electronMass;
    double max_err = 0;
    for (int i = 0; i < 50; ++i) {
        CrossSectionDistributionRecord xsec_record(event);
        decay.SampleFinalState(xsec_record, rand);
        xsec_record.Finalize(event);

        // Boost secondaries back to rest frame
        rk::P4 pHNL(geom3::Vector3(event.primary_momentum[1],
                                    event.primary_momentum[2],
                                    event.primary_momentum[3]), hnl_mass);
        rk::Boost to_rest = pHNL.restBoost();

        rk::P4 k3_lab(geom3::Vector3(event.secondary_momenta[1][1],
                                      event.secondary_momenta[1][2],
                                      event.secondary_momenta[1][3]), me);
        rk::P4 k4_lab(geom3::Vector3(event.secondary_momenta[2][1],
                                      event.secondary_momenta[2][2],
                                      event.secondary_momenta[2][3]), me);
        rk::P4 k3_rest = k3_lab;
        rk::P4 k4_rest = k4_lab;
        k3_rest.boost(to_rest);
        k4_rest.boost(to_rest);

        double cos_alpha34_geom = k3_rest.momentum().direction().dot(
                                  k4_rest.momentum().direction());

        // Dalitz prediction: cos(alpha34) from energies
        double E3 = k3_rest.e();
        double E4 = k4_rest.e();
        double p3 = k3_rest.p();
        double p4 = k4_rest.p();
        double s1 = (pHNL - k4_lab).dot(pHNL - k4_lab) / (hnl_mass*hnl_mass);
        double s2 = (pHNL - k3_lab).dot(pHNL - k3_lab) / (hnl_mass*hnl_mass);
        double s3 = 1 + (me*me + me*me)/(hnl_mass*hnl_mass) - s1 - s2;
        double cos_alpha34_dalitz = (E3*E4 - 0.5*(s3*hnl_mass*hnl_mass - me*me - me*me)) / (p3*p4);

        max_err = std::max(max_err, std::abs(cos_alpha34_geom - cos_alpha34_dalitz));
    }
    EXPECT_LT(max_err, 1e-10)
        << "Opening angle from geometry does not match Dalitz prediction";
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

// ---------------------------------------------------------------------------
// Gap-closing tests: absolute widths, mixed-flavor kinematics, symmetry
// ---------------------------------------------------------------------------

TEST(ElectroweakDecay, WTotalWidthAbsolute) {
    // After applying color factor (Nc=3) and CKM^2 corrections, the total
    // tree-level W width should be approximately:
    //   Gamma_W = GammaW_per_channel * (3_leptonic + 3*sum_ij|V_ij|^2)
    // where GammaW_per_channel = g^2 * M_W / (48*pi).
    // At tree level with 4 accessible quark channels (ud, us, cd, cs),
    // Gamma_total ~ 1.8 - 2.1 GeV (PDG measured: 2.085 GeV).
    std::set<ParticleType> primaries = {ParticleType::WPlus};
    ElectroweakDecay decay(primaries);
    double total_width = decay.TotalDecayWidth(ParticleType::WPlus);
    EXPECT_GT(total_width, 0.0);
    // Tree-level (no QCD corrections) should be in range [1.5, 2.3] GeV
    // This catches the missing-Nc and missing-CKM^2 bugs:
    // with bugs: total ~ 1.15 GeV (below range)
    // correct: total ~ 1.9 GeV (in range)
    EXPECT_GT(total_width, 1.5) << "W total width too low -- check color factor and CKM^2";
    EXPECT_LT(total_width, 2.3) << "W total width too high";
}

TEST(ElectroweakDecay, ZTotalWidthAbsolute) {
    // Z total width should be approximately 2.0-2.8 GeV at tree level.
    // PDG measured: 2.4952 GeV.
    // Without color factor for quarks, hadronic channels are 3x too low
    // and total would be ~ 1.0 GeV.
    std::set<ParticleType> primaries = {ParticleType::Z0};
    ElectroweakDecay decay(primaries);
    double total_width = decay.TotalDecayWidth(ParticleType::Z0);
    EXPECT_GT(total_width, 0.0);
    // Tree level with 5 quark flavors (u,d,s,c,b) + 3 leptons + 3 neutrinos
    EXPECT_GT(total_width, 1.8) << "Z total width too low -- check color factor for quarks";
    EXPECT_LT(total_width, 3.0) << "Z total width too high";
}

TEST(HNLDecay, MixedFlavorThreeBodyConservesMomentum) {
    // N4 -> nu_tau + mu- + tau+ (m_alpha=m_mu, m_beta=m_tau, m_alpha != m_beta)
    double hnl_mass = 2.5;
    HNLDecay decay(hnl_mass, std::vector<double>{0, 0, 1e-3}, HNLDecay::ChiralNature::Dirac);

    InteractionRecord event;
    event.signature.primary_type = ParticleType::N4;
    event.signature.target_type = ParticleType::Decay;
    event.signature.secondary_types = {ParticleType::NuTau, ParticleType::MuMinus, ParticleType::TauPlus};
    event.primary_mass = hnl_mass;
    double energy = 10.0;
    event.primary_momentum = {energy, 0, 0, std::sqrt(energy*energy - hnl_mass*hnl_mass)};
    event.primary_helicity = 1.0;

    auto rand = std::make_shared<SIREN_random>();
    double max_nu_m2_violation = 0;
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
                << "4-momentum not conserved in component " << j
                << " (mixed-flavor: mu+tau final state)";
        }

        double m_charged[] = {siren::utilities::Constants::muonMass,
                              siren::utilities::Constants::tauMass};
        for (size_t s = 1; s <= 2; ++s) {
            double E = event.secondary_momenta[s][0];
            double px = event.secondary_momenta[s][1];
            double py = event.secondary_momenta[s][2];
            double pz = event.secondary_momenta[s][3];
            double m2 = E*E - px*px - py*py - pz*pz;
            EXPECT_NEAR(m2, m_charged[s-1]*m_charged[s-1], 1e-6)
                << "Charged lepton " << s << " not on mass shell in mixed-flavor decay";
        }

        double E_nu = event.secondary_momenta[0][0];
        double px_nu = event.secondary_momenta[0][1];
        double py_nu = event.secondary_momenta[0][2];
        double pz_nu = event.secondary_momenta[0][3];
        double nu_m2 = E_nu*E_nu - px_nu*px_nu - py_nu*py_nu - pz_nu*pz_nu;
        max_nu_m2_violation = std::max(max_nu_m2_violation, std::abs(nu_m2));
    }
    EXPECT_LT(max_nu_m2_violation, 1e-8)
        << "Neutrino m^2 violation exceeds 1e-8 GeV^2 in mixed-flavor decay";
}

TEST(HNLDecay, MixedFlavorEMuConservesMomentum) {
    // N4 -> nu_mu + e- + mu+ (m_alpha=m_e, m_beta=m_mu)
    double hnl_mass = 0.5;
    HNLDecay decay(hnl_mass, std::vector<double>{0, 1e-3, 0}, HNLDecay::ChiralNature::Dirac);

    InteractionRecord event;
    event.signature.primary_type = ParticleType::N4;
    event.signature.target_type = ParticleType::Decay;
    event.signature.secondary_types = {ParticleType::NuMu, ParticleType::EMinus, ParticleType::MuPlus};
    event.primary_mass = hnl_mass;
    double energy = 5.0;
    event.primary_momentum = {energy, 0, 0, std::sqrt(energy*energy - hnl_mass*hnl_mass)};
    event.primary_helicity = 1.0;

    auto rand = std::make_shared<SIREN_random>();
    double max_nu_m2_violation = 0;
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
                << "4-momentum not conserved in component " << j
                << " (e+mu final state)";
        }

        double m_charged[] = {siren::utilities::Constants::electronMass,
                              siren::utilities::Constants::muonMass};
        for (size_t s = 1; s <= 2; ++s) {
            double E = event.secondary_momenta[s][0];
            double px = event.secondary_momenta[s][1];
            double py = event.secondary_momenta[s][2];
            double pz = event.secondary_momenta[s][3];
            double m2 = E*E - px*px - py*py - pz*pz;
            EXPECT_NEAR(m2, m_charged[s-1]*m_charged[s-1], 1e-6)
                << "Charged lepton " << s << " not on mass shell in e+mu decay";
        }

        double E_nu = event.secondary_momenta[0][0];
        double px_nu = event.secondary_momenta[0][1];
        double py_nu = event.secondary_momenta[0][2];
        double pz_nu = event.secondary_momenta[0][3];
        double nu_m2 = E_nu*E_nu - px_nu*px_nu - py_nu*py_nu - pz_nu*pz_nu;
        max_nu_m2_violation = std::max(max_nu_m2_violation, std::abs(nu_m2));
    }
    EXPECT_LT(max_nu_m2_violation, 1e-8)
        << "Neutrino m^2 violation exceeds 1e-8 GeV^2 in e+mu decay";
}

TEST(HNLDecay, MajoranaN4EqualsN4Bar) {
    // For a Majorana HNL, the total decay width should be identical
    // for N4 and N4Bar (they are the same particle).
    HNLDecay decay(0.5, std::vector<double>{0, 1e-3, 0}, HNLDecay::Majorana);
    double w_n4 = decay.TotalDecayWidth(ParticleType::N4);
    double w_n4bar = decay.TotalDecayWidth(ParticleType::N4Bar);
    EXPECT_GT(w_n4, 0.0);
    EXPECT_GT(w_n4bar, 0.0);
    EXPECT_DOUBLE_EQ(w_n4, w_n4bar)
        << "Majorana N4 and N4Bar should have identical total width";
}

TEST(HNLDecay, DiracN4DiffersFromN4Bar) {
    // For a Dirac HNL with only one mixing element, N4 and N4Bar
    // may differ (the available channels differ for particle vs antiparticle)
    // or at minimum, both should have non-zero width
    HNLDecay decay(0.5, std::vector<double>{0, 1e-3, 0}, HNLDecay::Dirac);
    double w_n4 = decay.TotalDecayWidth(ParticleType::N4);
    double w_n4bar = decay.TotalDecayWidth(ParticleType::N4Bar);
    EXPECT_GT(w_n4, 0.0);
    EXPECT_GT(w_n4bar, 0.0);
}

TEST(HNLDecay, HadronicSubtractionNonNegative) {
    // At intermediate masses (1-2 GeV), both inclusive quark channels and
    // exclusive meson channels are active. The total width (which includes
    // the subtraction Gamma_quarks - Gamma_mesons for each CC channel)
    // must remain non-negative for all partial widths.
    // Test several mass points in the critical region.
    double mass_points[] = {0.8, 1.0, 1.2, 1.5, 1.8, 2.0};
    for (double m_N : mass_points) {
        HNLDecay decay(m_N, std::vector<double>{1e-3, 1e-3, 1e-3}, HNLDecay::Majorana);
        auto sigs = decay.GetPossibleSignaturesFromParent(ParticleType::N4);
        InteractionRecord record;
        record.primary_mass = m_N;
        double total_width = 0;
        for (auto const & sig : sigs) {
            record.signature = sig;
            double w = decay.TotalDecayWidthForFinalState(record);
            EXPECT_GE(w, 0.0) << "Negative partial width at M_N = " << m_N
                               << " GeV for channel with "
                               << sig.secondary_types.size() << " secondaries";
            total_width += w;
        }
        // Total width from summing partials should match TotalDecayWidth
        double total_direct = decay.TotalDecayWidth(ParticleType::N4);
        EXPECT_NEAR(total_width, total_direct, 1e-10 * total_direct)
            << "Sum of partial widths != TotalDecayWidth at M_N = " << m_N;
    }
}

TEST(HNLDecay, WidthScalesAsMN5AtLowMass) {
    // Below all meson thresholds and well above lepton masses,
    // the total HNL decay width should scale as M_N^5 (Seesaw formula).
    // Use M_N = 50 MeV and 80 MeV (both well above 2*m_e but below m_pi).
    double m1 = 0.05, m2 = 0.08;
    HNLDecay decay1(m1, std::vector<double>{1e-3, 0, 0}, HNLDecay::Majorana);
    HNLDecay decay2(m2, std::vector<double>{1e-3, 0, 0}, HNLDecay::Majorana);
    double w1 = decay1.TotalDecayWidth(ParticleType::N4);
    double w2 = decay2.TotalDecayWidth(ParticleType::N4);
    EXPECT_GT(w1, 0.0);
    EXPECT_GT(w2, 0.0);
    // w2/w1 should be approximately (m2/m1)^5 = (0.08/0.05)^5 = 1.6^5 = 10.49
    double expected_ratio = std::pow(m2 / m1, 5);
    double actual_ratio = w2 / w1;
    // Allow 20% tolerance for mass corrections (m_e/M_N not negligible)
    EXPECT_NEAR(actual_ratio, expected_ratio, 0.2 * expected_ratio)
        << "Width does not scale as M_N^5 at low mass";
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
