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
#endif // SIREN_HAS_PYTHIA8

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
