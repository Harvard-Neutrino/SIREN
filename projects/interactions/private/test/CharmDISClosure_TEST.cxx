// Regression test for the charm-DIS interaction-depth closure invariant: the
// generation side (TotalCrossSectionAllFinalStates) and physical side (sum of
// per-signature TotalCrossSection, Weighter.tcc) must agree, and the inclusive
// sigma must be partitioned across the three D species by fragmentation fraction
// (else the sum triple-counts). Pythia-free: a mock reproduces PythiaDIS's two
// relevant behaviors and exercises the real base-class TotalCrossSectionAllFinalStates.
// Companion PythiaDISCharmClosure_TEST exercises the real class.

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

    // D0:D+/-:Ds = 0.60:0.23:0.15 each /0.98 to sum to 1.0 (unmodeled Lambda_c
    // redistributed). Mirrors QuarkDISFromSpline::FragmentationFraction.
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

// f1751c6b bug: override makes the generation side report 1x while the physical
// side reports 3x -> closure broken by a factor of 3.
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

// 4b7baf4a (override removed, no FF): sides agree (closure restored) but both
// equal 3x the inclusive sigma -> charm production overcounted.
TEST(CharmDISClosure, OverrideRemovedClosesButOvercountsThreeX) {
    MockCharmXS xs(/*ff=*/false, /*override=*/false);
    for(double E : {10.0, 100.0, 1000.0}) {
        double gen = gen_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double phys = phys_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        EXPECT_NEAR(gen, phys, 1e-50);                                  // closure ok
        EXPECT_NEAR(gen, 3.0 * MockCharmXS::sigma_inclusive(E), 1e-50); // but 3x high
    }
}

// With FF applied per signature (FFs renormalized to sum to 1.0), both sides
// agree AND equal the full inclusive sigma -- partitioned, not triple-counted.
TEST(CharmDISClosure, FragmentationFractionRestoresPhysicalNormalization) {
    MockCharmXS xs(/*ff=*/true, /*override=*/false);
    for(double E : {10.0, 100.0, 1000.0}) {
        double gen = gen_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double phys = phys_path(xs, ParticleType::NuMu, ParticleType::PPlus, E);
        double s = MockCharmXS::sigma_inclusive(E);
        EXPECT_NEAR(gen, phys, 1e-50);                 // closure ok
        EXPECT_NEAR(gen, s, 1e-48);                    // and physically normalized (== inclusive sigma)
    }
}

// The three implemented fragmentation fractions sum to 1.0.
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
