// Regression test for the charm-DIS interaction-depth closure invariant.
//
// Background: a charm-DIS cross section registers three D-type final-state
// signatures (D0, D+, Ds) that share a single inclusive charm total cross
// section sigma(E). The interaction depth that sets the vertex distribution is
// computed two ways that MUST agree:
//
//   generation side : CrossSection::TotalCrossSectionAllFinalStates(record)
//                     (called by the ~10 vertex/position distributions)
//   physical side   : sum over GetPossibleSignaturesFromParents of
//                     TotalCrossSection(signature)            (Weighter.tcc)
//
// The base-class default TotalCrossSectionAllFinalStates SUMS TotalCrossSection
// over the signatures, so the two sides agree only if no subclass overrides
// TotalCrossSectionAllFinalStates to short-circuit the sum. Additionally, the
// inclusive sigma must be PARTITIONED across the D species by fragmentation
// fraction, otherwise the sum triple-counts charm production.
//
// This test does not need Pythia8: it uses a mock that reproduces the two
// behaviors of PythiaDISCrossSection (TotalCrossSection independent of meson
// type; three registered D-type signatures) and exercises the real base-class
// CrossSection::TotalCrossSectionAllFinalStates. It guards the invariant the
// PythiaDISCrossSection fix relies on. The companion gtest
// PythiaDISCharmClosure_TEST exercises the real class on a Pythia-enabled build.

#include <cmath>
#include <vector>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"

using namespace siren::interactions;
using namespace siren::dataclasses;

namespace {

// Mock reproducing PythiaDISCrossSection's relevant semantics.
//   apply_ff     : multiply by FragmentationFraction per signature (as in QuarkDISFromSpline)
//   override_afs : short-circuit TotalCrossSectionAllFinalStates to TotalCrossSection
//                  instead of the base per-signature sum
class MockCharmXS : public CrossSection {
public:
    bool apply_ff;
    bool override_afs;
    MockCharmXS(bool ff, bool ovr) : apply_ff(ff), override_afs(ovr) {}

    // D0:D+/-:Ds = 0.60:0.23:0.15 renormalized to sum to 1.0 (Lambda_c not
    // modeled, its fraction redistributed). Mirrors
    // QuarkDISFromSpline::FragmentationFraction.
    static double FragmentationFraction(ParticleType d) {
        if(d==ParticleType::D0 || d==ParticleType::D0Bar) return 0.6 / 0.98;
        if(d==ParticleType::DPlus || d==ParticleType::DMinus) return 0.23 / 0.98;
        if(d==ParticleType::DsPlus || d==ParticleType::DsMinus) return 0.15 / 0.98;
        return 0.0;
    }
    // Stand-in for the 1-D inclusive charm total spline (independent of meson type).
    static double sigma_inclusive(double E) { return 1.0e-38 * std::log10(E); }

    std::vector<InteractionSignature> sigs(ParticleType primary, ParticleType target) const {
        InteractionSignature base;
        base.primary_type = primary;
        base.target_type = target;
        base.secondary_types = { ParticleType::MuMinus, ParticleType::D0 };
        std::vector<InteractionSignature> out;
        for(ParticleType d : { ParticleType::D0, ParticleType::DPlus, ParticleType::DsPlus }) {
            InteractionSignature s = base;
            s.secondary_types[1] = d;
            out.push_back(s);
        }
        return out;
    }

    double TotalCrossSection(InteractionRecord const & r) const override {
        double s = sigma_inclusive(r.primary_momentum[0]);
        if(apply_ff) {
            for(auto t : r.signature.secondary_types)
                if(isD(t)) { s *= FragmentationFraction(t); break; }
        }
        return s;
    }
    double TotalCrossSectionAllFinalStates(InteractionRecord const & r) const override {
        if(override_afs) return TotalCrossSection(r);
        return CrossSection::TotalCrossSectionAllFinalStates(r);
    }
    std::vector<InteractionSignature> GetPossibleSignaturesFromParents(ParticleType p, ParticleType t) const override {
        return sigs(p, t);
    }
    std::vector<InteractionSignature> GetPossibleSignatures() const override { return sigs(ParticleType::NuMu, ParticleType::PPlus); }
    std::vector<ParticleType> GetPossibleTargets() const override { return { ParticleType::PPlus }; }
    std::vector<ParticleType> GetPossibleTargetsFromPrimary(ParticleType) const override { return { ParticleType::PPlus }; }
    std::vector<ParticleType> GetPossiblePrimaries() const override { return { ParticleType::NuMu }; }
    double DifferentialCrossSection(InteractionRecord const &) const override { return 0.0; }
    double InteractionThreshold(InteractionRecord const &) const override { return 0.0; }
    double FinalStateProbability(InteractionRecord const &) const override { return 0.0; }
    std::vector<std::string> DensityVariables() const override { return {}; }
    void SampleFinalState(CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override {}
    bool equal(CrossSection const &) const override { return false; }
};

// Replicates the Weighter.tcc physical-side loop (InteractionProbability,
// NormalizedPositionProbability) and WeightingUtils CrossSectionProbability:
// sum TotalCrossSection over GetPossibleSignaturesFromParents.
double phys_path(MockCharmXS const & xs, ParticleType primary, ParticleType target, double E) {
    InteractionRecord fake;
    fake.signature.primary_type = primary;
    fake.signature.target_type = target;
    fake.primary_momentum[0] = E;
    double total = 0.0;
    for(auto const & sig : xs.GetPossibleSignaturesFromParents(primary, target)) {
        fake.signature = sig;
        fake.primary_momentum[0] = E;
        total += xs.TotalCrossSection(fake);
    }
    return total;
}

double gen_path(MockCharmXS const & xs, ParticleType primary, ParticleType target, double E) {
    InteractionRecord rec;
    rec.signature.primary_type = primary;
    rec.signature.target_type = target;
    rec.primary_momentum[0] = E;
    return xs.TotalCrossSectionAllFinalStates(rec);
}

} // namespace

// Reproduces the f1751c6b bug: the override makes the generation side report 1x
// while the physical side reports 3x -> closure broken by a factor of 3.
TEST(CharmDISClosure, OverrideBreaksClosureByFactorThree) {
    MockCharmXS xs(/*ff=*/false, /*override=*/true);
    for(double E : {10.0, 100.0, 1000.0}) {
        double gen = gen_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double phys = phys_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        EXPECT_NEAR(gen, MockCharmXS::sigma_inclusive(E), 1e-50);   // 1x
        EXPECT_NEAR(phys, 3.0 * MockCharmXS::sigma_inclusive(E), 1e-50); // 3x
        EXPECT_NEAR(gen / phys, 1.0 / 3.0, 1e-9);
    }
}

// Reproduces 4b7baf4a (override removed, no FF): the two sides agree (closure
// restored) but both equal 3x the inclusive sigma -> charm production overcounted.
TEST(CharmDISClosure, OverrideRemovedClosesButOvercountsThreeX) {
    MockCharmXS xs(/*ff=*/false, /*override=*/false);
    for(double E : {10.0, 100.0, 1000.0}) {
        double gen = gen_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double phys = phys_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        EXPECT_NEAR(gen, phys, 1e-50);                                  // closure ok
        EXPECT_NEAR(gen, 3.0 * MockCharmXS::sigma_inclusive(E), 1e-50); // but 3x high
    }
}

// With the fragmentation fraction applied per signature in TotalCrossSection
// (and the FFs renormalized to sum to 1.0), the generation-side and physical-
// side inclusive charm cross sections must agree AND both equal the full
// inclusive sigma -- partitioned across the three D species, not triple-counted.
TEST(CharmDISClosure, FragmentationFractionRestoresPhysicalNormalization) {
    MockCharmXS xs(/*ff=*/true, /*override=*/false);
    // FF sum = (0.6 + 0.23 + 0.15) / 0.98 == 1.0 (renormalized), so the
    // partitioned total recovers the full inclusive sigma exactly.
    for(double E : {10.0, 100.0, 1000.0}) {
        double gen = gen_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double phys = phys_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double s = MockCharmXS::sigma_inclusive(E);
        EXPECT_NEAR(gen, phys, 1e-50);                 // closure ok
        EXPECT_NEAR(gen, s, 1e-48);                    // and physically normalized (== inclusive sigma)
    }
}

// The three implemented fragmentation fractions must partition the inclusive
// charm cross section exactly: they sum to 1.0 (Lambda_c fraction redistributed).
TEST(CharmDISClosure, FragmentationFractionsSumToOne) {
    double sum = MockCharmXS::FragmentationFraction(ParticleType::D0)
               + MockCharmXS::FragmentationFraction(ParticleType::DPlus)
               + MockCharmXS::FragmentationFraction(ParticleType::DsPlus);
    EXPECT_NEAR(sum, 1.0, 1e-12);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
