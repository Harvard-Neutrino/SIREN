// Regression tests for the PythiaDISCrossSection charm-DIS interaction-depth
// closure bug. Guards that per-signature TotalCrossSection is partitioned by
// fragmentation fraction (sum == inclusive sigma, not 3x) and that the
// generation side (TotalCrossSectionAllFinalStates) equals the physical sum.
// See PythiaDISCrossSection.h. Requires SIREN_WITH_PYTHIA8 plus spline files in
// env vars SIREN_PYTHIA_TEST_DSDXDY / SIREN_PYTHIA_TEST_SIGMA; SKIPs if unset.

#include <cmath>
#include <cstdlib>
#include <vector>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"

#ifdef SIREN_HAS_PYTHIA8
#include "SIREN/interactions/PythiaDISCrossSection.h"
#endif

using namespace siren::interactions;
using namespace siren::dataclasses;

namespace {
double expected_ff(ParticleType d) {
    // Raw fractions / 0.98 to fold in the unmodeled Lambda_c and renormalize to
    // sum to 1.0. Mirrors PythiaDISCrossSection::FragmentationFraction.
    if(d==ParticleType::D0 || d==ParticleType::D0Bar) return 0.6 / 0.98;
    if(d==ParticleType::DPlus || d==ParticleType::DMinus) return 0.23 / 0.98;
    if(d==ParticleType::DsPlus || d==ParticleType::DsMinus) return 0.15 / 0.98;
    return 0.0;
}
} // namespace

#ifdef SIREN_HAS_PYTHIA8
TEST(PythiaDISCharmClosure, FragmentationPartitionAndClosure) {
    const char* dd = std::getenv("SIREN_PYTHIA_TEST_DSDXDY");
    const char* tt = std::getenv("SIREN_PYTHIA_TEST_SIGMA");
    if(!dd || !tt) {
        GTEST_SKIP() << "Set SIREN_PYTHIA_TEST_DSDXDY and SIREN_PYTHIA_TEST_SIGMA to run.";
    }

    std::vector<ParticleType> primaries = { ParticleType::NuMu };
    std::vector<ParticleType> targets   = { ParticleType::PPlus };
    PythiaDISCrossSection xs(
        std::string(dd), std::string(tt),
        /*interaction_type=*/1, /*target_mass=*/0.9382720813,
        /*minimum_Q2=*/1.0, primaries, targets,
        /*pythia_data_path=*/"", /*pdf_set=*/"LHAPDF6:CT18NLO", /*units=*/"cm");

    const double ff_sum = (0.6 + 0.23 + 0.15) / 0.98;   // renormalized -> 1.0

    for(double E : {30.0, 100.0, 300.0}) {
        double sigma_incl = xs.TotalCrossSection(ParticleType::NuMu, E);
        ASSERT_GT(sigma_incl, 0.0);

        std::vector<InteractionSignature> sigs =
            xs.GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::PPlus);
        ASSERT_EQ(sigs.size(), 3u);

        double sum = 0.0;
        std::vector<double> per_sig;
        for(auto const & sig : sigs) {
            InteractionRecord rec;
            rec.signature = sig;
            rec.primary_momentum[0] = E;
            double v = xs.TotalCrossSection(rec);
            sum += v;
            per_sig.push_back(v);

            ParticleType d = ParticleType::unknown;
            for(auto t : sig.secondary_types) if(isD(t)) { d = t; break; }
            ASSERT_NE(d, ParticleType::unknown);
            EXPECT_NEAR(v, sigma_incl * expected_ff(d), sigma_incl * 1e-9)
                << "signature D-type cross section is not fragmentation-partitioned at E=" << E;
        }

        // The three values must differ (0.6 vs 0.23 vs 0.15), not all == sigma_incl.
        EXPECT_GT(std::abs(per_sig[0] - per_sig[1]), sigma_incl * 1e-6);

        // Sum is the partitioned inclusive sigma, not 3x.
        EXPECT_NEAR(sum, sigma_incl * ff_sum, sigma_incl * 1e-9);
        EXPECT_LT(sum, sigma_incl * 1.5); // would be 3x without the fix

        // Generation side == physical side (sum).
        InteractionRecord rec;
        rec.signature.primary_type = ParticleType::NuMu;
        rec.signature.target_type = ParticleType::PPlus;
        rec.primary_momentum[0] = E;
        double gen = xs.TotalCrossSectionAllFinalStates(rec);
        EXPECT_NEAR(gen, sum, sigma_incl * 1e-9) << "generation/physical interaction depth mismatch";
    }
}

// With only a total spline (empty differential filename), FinalStateProbability
// returns a constant 1.0: the intractable Pythia final-state density cancels in
// the unbiased weight, so only the total cross section matters.
TEST(PythiaDISCharmClosure, ConstantFinalStateProbabilityWithoutDifferential) {
    const char* tt = std::getenv("SIREN_PYTHIA_TEST_SIGMA");
    if(!tt) {
        GTEST_SKIP() << "Set SIREN_PYTHIA_TEST_SIGMA to run.";
    }
    std::vector<ParticleType> primaries = { ParticleType::NuMu };
    std::vector<ParticleType> targets   = { ParticleType::PPlus };
    PythiaDISCrossSection xs(
        /*differential=*/std::string(""), std::string(tt),
        /*interaction_type=*/1, /*target_mass=*/0.9382720813,
        /*minimum_Q2=*/1.0, primaries, targets,
        /*pythia_data_path=*/"", /*pdf_set=*/"LHAPDF6:CT18NLO", /*units=*/"cm");

    EXPECT_GT(xs.TotalCrossSection(ParticleType::NuMu, 100.0), 0.0);

    // FSP is exactly 1.0: the no-differential branch returns before kinematics.
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::PPlus)[0];
    InteractionRecord rec;
    rec.signature = sig;
    rec.primary_momentum = {100.0, 0.0, 0.0, 100.0};
    rec.target_mass = 0.9382720813;
    rec.secondary_momenta = {{40.0, 1.0, 0.0, 39.0}, {0.0, 0.0, 0.0, 0.0}, {1.86, 0.0, 0.0, 1.0}};
    rec.secondary_masses = {0.105, 0.0, 1.86};
    EXPECT_EQ(xs.FinalStateProbability(rec), 1.0);
}

// Charm is forced only in CC; NC charm must be rejected at construction (not
// later inside SampleFinalState) with an actionable error.
TEST(PythiaDISCharmClosure, NeutralCurrentRejectedAtConstruction) {
    const char* tt = std::getenv("SIREN_PYTHIA_TEST_SIGMA");
    if(!tt) GTEST_SKIP() << "Set SIREN_PYTHIA_TEST_SIGMA to run.";
    std::vector<ParticleType> primaries = { ParticleType::NuMu };
    std::vector<ParticleType> targets   = { ParticleType::PPlus };
    EXPECT_THROW({
        PythiaDISCrossSection xs(std::string(""), std::string(tt),
            /*interaction_type=*/2, 0.9382720813, 1.0, primaries, targets,
            "", "LHAPDF6:CT18NLO", "cm");
    }, std::runtime_error);
    // Positive control: CC still constructs cleanly.
    EXPECT_NO_THROW({
        PythiaDISCrossSection xs(std::string(""), std::string(tt),
            /*interaction_type=*/1, 0.9382720813, 1.0, primaries, targets,
            "", "LHAPDF6:CT18NLO", "cm");
    });
}

// Species tripwire: pin the renormalized fractions, the raw sum (0.98), and
// Lambda_c -> 0 so modeling Lambda_c forces a test update rather than silently
// shifting the species mix.
TEST(PythiaDISCharmClosure, FragmentationFractionSpeciesTripwire) {
    const char* tt = std::getenv("SIREN_PYTHIA_TEST_SIGMA");
    if(!tt) GTEST_SKIP() << "Set SIREN_PYTHIA_TEST_SIGMA to run.";
    std::vector<ParticleType> primaries = { ParticleType::NuMu };
    std::vector<ParticleType> targets   = { ParticleType::PPlus };
    PythiaDISCrossSection xs(std::string(""), std::string(tt),
        1, 0.9382720813, 1.0, primaries, targets, "", "LHAPDF6:CT18NLO", "cm");

    EXPECT_NEAR(xs.FragmentationFraction(ParticleType::D0),     0.6 / 0.98, 1e-12);
    EXPECT_NEAR(xs.FragmentationFraction(ParticleType::DPlus),  0.23 / 0.98, 1e-12);
    EXPECT_NEAR(xs.FragmentationFraction(ParticleType::DsPlus), 0.15 / 0.98, 1e-12);
    // Renormalized fractions recover the full inclusive sigma.
    double sum = xs.FragmentationFraction(ParticleType::D0)
               + xs.FragmentationFraction(ParticleType::DPlus)
               + xs.FragmentationFraction(ParticleType::DsPlus);
    EXPECT_NEAR(sum, 1.0, 1e-12);
    // Raw fractions sum to 0.98; the 0.02 gap is the folded unmodeled baryon.
    EXPECT_NEAR(0.6 + 0.23 + 0.15, 0.98, 1e-12);
    // Lambda_c (PDG 4122) intentionally unmodeled -> FF 0. Modeling it must update this.
    EXPECT_EQ(xs.FragmentationFraction(static_cast<ParticleType>(4122)), 0.0);
}

// Closure + partition must hold across the analysis band (TeV-PeV), not only
// <=300 GeV. Uses a wide total spline (100 GeV - 1 PeV) via SIREN_PYTHIA_WIDE_SIGMA.
TEST(PythiaDISCharmClosure, FragmentationClosureAtAnalysisEnergies) {
    const char* tt = std::getenv("SIREN_PYTHIA_WIDE_SIGMA");
    if(!tt) GTEST_SKIP() << "Set SIREN_PYTHIA_WIDE_SIGMA (wide total spline, 100 GeV-1 PeV) to run.";
    std::vector<ParticleType> primaries = { ParticleType::NuMu };
    std::vector<ParticleType> targets   = { ParticleType::PPlus };
    PythiaDISCrossSection xs(std::string(""), std::string(tt),
        1, 0.9382720813, 1.0, primaries, targets, "", "LHAPDF6:CT18NLO", "cm");
    const double ff_sum = (0.6 + 0.23 + 0.15) / 0.98;

    for(double E : {1.0e4, 1.0e5, 1.0e6}) {   // 10 TeV, 100 TeV, 1 PeV
        double sigma_incl = 0.0;
        ASSERT_NO_THROW(sigma_incl = xs.TotalCrossSection(ParticleType::NuMu, E))
            << "total spline does not cover E=" << E << " GeV";
        ASSERT_GT(sigma_incl, 0.0);

        auto sigs = xs.GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::PPlus);
        ASSERT_EQ(sigs.size(), 3u);
        double sum = 0.0;
        for(auto const & sig : sigs) {
            InteractionRecord rec;
            rec.signature = sig;
            rec.primary_momentum[0] = E;
            double v = xs.TotalCrossSection(rec);
            sum += v;
            ParticleType d = ParticleType::unknown;
            for(auto t : sig.secondary_types) if(isD(t)) { d = t; break; }
            ASSERT_NE(d, ParticleType::unknown);
            EXPECT_NEAR(v, sigma_incl * expected_ff(d), sigma_incl * 1e-9) << "E=" << E;
        }
        EXPECT_NEAR(sum, sigma_incl * ff_sum, sigma_incl * 1e-9) << "E=" << E;

        InteractionRecord rec;
        rec.signature.primary_type = ParticleType::NuMu;
        rec.signature.target_type = ParticleType::PPlus;
        rec.primary_momentum[0] = E;
        double gen = xs.TotalCrossSectionAllFinalStates(rec);
        EXPECT_NEAR(gen, sum, sigma_incl * 1e-9) << "generation/physical mismatch at E=" << E;
    }
}

// Out-of-range differential evaluation must RAISE, not return 0 (a silent zero
// on a sampled event biases its weight).
TEST(PythiaDISCharmClosure, DifferentialOutOfRangeRaises) {
    const char* dd = std::getenv("SIREN_PYTHIA_TEST_DSDXDY");
    const char* tt = std::getenv("SIREN_PYTHIA_TEST_SIGMA");
    if(!dd || !tt) GTEST_SKIP() << "Set SIREN_PYTHIA_TEST_DSDXDY and SIREN_PYTHIA_TEST_SIGMA to run.";
    std::vector<ParticleType> primaries = { ParticleType::NuMu };
    std::vector<ParticleType> targets   = { ParticleType::PPlus };
    PythiaDISCrossSection xs(std::string(dd), std::string(tt),
        1, 0.9382720813, 1.0, primaries, targets, "", "LHAPDF6:CT18NLO", "cm");

    // In range -> finite positive density (Q2 computed from x,y).
    EXPECT_GT(xs.DifferentialCrossSection(100.0, 0.1, 0.5, 0.105), 0.0);
    // Energy outside the spline extent raises.
    EXPECT_THROW(xs.DifferentialCrossSection(1.0e12, 0.1, 0.5, 0.105), std::runtime_error);
    // (x, y) outside the spline grid raises (explicit Q2 bypasses the Q2 cut).
    EXPECT_THROW(xs.DifferentialCrossSection(100.0, 1.0e-8, 0.5, 0.105, 5.0), std::runtime_error);
}
#endif // SIREN_HAS_PYTHIA8

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
