// Regression test for the PythiaDISCrossSection charm-DIS interaction-depth
// closure / normalization bug.
//
// Requires a Pythia8-enabled build (SIREN_WITH_PYTHIA8=ON) AND total/differential
// charm-DIS spline files supplied via environment variables:
//
//   SIREN_PYTHIA_TEST_DSDXDY  -> differential (dsdxdy) FITS spline
//   SIREN_PYTHIA_TEST_SIGMA   -> total (sigma) FITS spline
//
// If the variables are unset the test SKIPs (it only needs the spline tables,
// not LHAPDF/Pythia at runtime: the total-cross-section path never calls Pythia).
//
// What it guards:
//   * Per-signature TotalCrossSection must be PARTITIONED by fragmentation
//     fraction (D0:0.6, D+:0.23, Ds:0.15), not equal to the full inclusive sigma.
//   * sum over signatures of TotalCrossSection == sigma_inclusive * 0.98 (NOT 3x).
//   * TotalCrossSectionAllFinalStates (generation side) == that sum (physical
//     side) -> closure holds.
//
// Before the fix (no fragmentation fraction in TotalCrossSection) every signature
// returns the full inclusive sigma, the sum is 3x too large, and the partition
// assertions FAIL. After the fix they PASS.

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
    // Renormalized to sum to 1.0 over the implemented D species (the unmodeled
    // Lambda_c fraction is folded in by dividing each raw value by 0.98). Mirrors
    // PythiaDISCrossSection::FragmentationFraction.
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

            // identify the D meson in this signature and check FF partition
            ParticleType d = ParticleType::unknown;
            for(auto t : sig.secondary_types) if(isD(t)) { d = t; break; }
            ASSERT_NE(d, ParticleType::unknown);
            EXPECT_NEAR(v, sigma_incl * expected_ff(d), sigma_incl * 1e-9)
                << "signature D-type cross section is not fragmentation-partitioned at E=" << E;
        }

        // The three values must differ (0.6 vs 0.23 vs 0.15), i.e. NOT all == sigma_incl.
        EXPECT_GT(std::abs(per_sig[0] - per_sig[1]), sigma_incl * 1e-6);

        // Sum is the inclusive sigma (partitioned), not 3x.
        EXPECT_NEAR(sum, sigma_incl * ff_sum, sigma_incl * 1e-9);
        EXPECT_LT(sum, sigma_incl * 1.5); // would be 3x without the fix

        // Generation side (TotalCrossSectionAllFinalStates) == physical side (sum).
        InteractionRecord rec;
        rec.signature.primary_type = ParticleType::NuMu;
        rec.signature.target_type = ParticleType::PPlus;
        rec.primary_momentum[0] = E;
        double gen = xs.TotalCrossSectionAllFinalStates(rec);
        EXPECT_NEAR(gen, sum, sigma_incl * 1e-9) << "generation/physical interaction depth mismatch";
    }
}

// The differential spline is optional. With only a total spline (empty
// differential filename), FinalStateProbability must return a constant 1.0 --
// the Pythia final-state density is intractable but cancels in the unbiased
// weight, so only the total cross section matters. Needs only the total spline.
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

    // Total still works.
    EXPECT_GT(xs.TotalCrossSection(ParticleType::NuMu, 100.0), 0.0);

    // FinalStateProbability is exactly 1.0 for any well-formed record, since the
    // no-differential branch returns before touching kinematics.
    auto sig = xs.GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::PPlus)[0];
    InteractionRecord rec;
    rec.signature = sig;
    rec.primary_momentum = {100.0, 0.0, 0.0, 100.0};
    rec.target_mass = 0.9382720813;
    rec.secondary_momenta = {{40.0, 1.0, 0.0, 39.0}, {0.0, 0.0, 0.0, 0.0}, {1.86, 0.0, 0.0, 1.0}};
    rec.secondary_masses = {0.105, 0.0, 1.86};
    EXPECT_EQ(xs.FinalStateProbability(rec), 1.0);
}

// PythiaDISCrossSection forces charm only in charged current. NC charm is not
// forced by Z exchange, so it must be rejected at construction with an
// actionable error rather than failing later inside SampleFinalState.
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

// Species/fragmentation tripwire. The three modeled D fractions are renormalized
// (each /0.98) so the partitioned signatures recover the inclusive charm cross
// section; the unmodeled Lambda_c (~0.02 of the D0/D+/Ds total) is folded in.
// Pin the raw sum (0.98) and Lambda_c -> 0 so adding Lambda_c forces a test
// update rather than silently shifting the species mix.
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
    // Renormalized fractions over the modeled species recover the full inclusive sigma.
    double sum = xs.FragmentationFraction(ParticleType::D0)
               + xs.FragmentationFraction(ParticleType::DPlus)
               + xs.FragmentationFraction(ParticleType::DsPlus);
    EXPECT_NEAR(sum, 1.0, 1e-12);
    // Raw fractions sum to 0.98; the 0.02 gap is the folded unmodeled-baryon fraction.
    EXPECT_NEAR(0.6 + 0.23 + 0.15, 0.98, 1e-12);
    // Lambda_c (PDG 4122) is intentionally unmodeled -> FF 0. Modeling it must update this.
    EXPECT_EQ(xs.FragmentationFraction(static_cast<ParticleType>(4122)), 0.0);
}

// Closure + fragmentation partition must hold across the ANALYSIS energy band
// (TeV-PeV), not only <=300 GeV. Uses a wide total spline (100 GeV - 1 PeV); a
// spline that silently threw or mis-partitioned at PeV would crash or bias the
// production run at exactly the analysis energies.
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

// Out-of-range differential evaluation must RAISE, not silently return 0 (a
// silent zero on a sampled event biases its weight). Needs a differential spline.
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
    // Energy outside the spline extent must raise.
    EXPECT_THROW(xs.DifferentialCrossSection(1.0e12, 0.1, 0.5, 0.105), std::runtime_error);
    // (x, y) outside the spline grid must raise (explicit Q2 bypasses the Q2 cut).
    EXPECT_THROW(xs.DifferentialCrossSection(100.0, 1.0e-8, 0.5, 0.105, 5.0), std::runtime_error);
}
#endif // SIREN_HAS_PYTHIA8

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
