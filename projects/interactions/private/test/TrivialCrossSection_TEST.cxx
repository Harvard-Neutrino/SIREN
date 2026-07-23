// TrivialCrossSection invariants: tabulated total cross section
// (interpolation, threshold, clamping), exact pass-through final-state
// kinematics, unit final-state probability for matching signatures, signature
// filtering by both primary and target, and equality-preserving polymorphic
// serialization.

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>

#include <gtest/gtest.h>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/TrivialCrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/utilities/Random.h"

using namespace siren::interactions;
using namespace siren::dataclasses;

namespace {

InteractionRecord MakeRecord(double energy,
        ParticleType primary = ParticleType::NuMu,
        ParticleType target = ParticleType::Nucleon) {
    InteractionRecord record;
    record.signature.primary_type = primary;
    record.signature.target_type = target;
    record.signature.secondary_types = {primary, target};
    record.primary_momentum = {energy, 0.0, 0.0, energy};
    record.primary_mass = 0.0;
    record.primary_helicity = -0.5;
    record.target_mass = 0.939;
    record.target_helicity = 0.5;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    return record;
}

std::shared_ptr<TrivialCrossSection> MakeTabulated() {
    std::map<ParticleType, std::pair<std::vector<double>, std::vector<double>>> tables;
    tables[ParticleType::NuMu] = {{1.0, 3.0}, {2e-38, 6e-38}};
    return std::make_shared<TrivialCrossSection>(tables,
            std::vector<ParticleType>{ParticleType::Nucleon});
}

} // namespace

TEST(TrivialCrossSection, ConstantTable) {
    TrivialCrossSection xs(1e-38,
            {ParticleType::NuMu, ParticleType::NuE},
            {ParticleType::Nucleon});
    EXPECT_DOUBLE_EQ(1e-38, xs.TotalCrossSection(ParticleType::NuMu, 0.5));
    EXPECT_DOUBLE_EQ(1e-38, xs.TotalCrossSection(ParticleType::NuE, 100.0));
    EXPECT_DOUBLE_EQ(0.0, xs.TotalCrossSection(ParticleType::NuTau, 1.0));
    EXPECT_DOUBLE_EQ(0.0, xs.InteractionThreshold(MakeRecord(1.0)));
}

TEST(TrivialCrossSection, InterpolationThresholdClamp) {
    std::shared_ptr<TrivialCrossSection> xs = MakeTabulated();
    // Knots reproduce exactly; midpoints interpolate linearly.
    EXPECT_DOUBLE_EQ(2e-38, xs->TotalCrossSection(ParticleType::NuMu, 1.0));
    EXPECT_DOUBLE_EQ(4e-38, xs->TotalCrossSection(ParticleType::NuMu, 2.0));
    // Below the first knot the cross section vanishes and the threshold
    // reports the first knot.
    EXPECT_DOUBLE_EQ(0.0, xs->TotalCrossSection(ParticleType::NuMu, 0.5));
    EXPECT_DOUBLE_EQ(1.0, xs->InteractionThreshold(MakeRecord(0.5)));
    // Above the last knot the value clamps.
    EXPECT_DOUBLE_EQ(6e-38, xs->TotalCrossSection(ParticleType::NuMu, 10.0));
    // Record-based lookup uses the primary energy.
    EXPECT_DOUBLE_EQ(4e-38, xs->TotalCrossSection(MakeRecord(2.0)));
}

TEST(TrivialCrossSection, InvalidConstruction) {
    EXPECT_THROW(TrivialCrossSection(0.0, {ParticleType::NuMu}, {ParticleType::Nucleon}),
            std::runtime_error);
    EXPECT_THROW(TrivialCrossSection(1e-38, {}, {ParticleType::Nucleon}),
            std::runtime_error);
    EXPECT_THROW(TrivialCrossSection(1e-38, {ParticleType::NuMu}, {}),
            std::runtime_error);
    std::map<ParticleType, std::pair<std::vector<double>, std::vector<double>>> bad_order;
    bad_order[ParticleType::NuMu] = {{3.0, 1.0}, {2e-38, 6e-38}};
    EXPECT_THROW(TrivialCrossSection(bad_order, {ParticleType::Nucleon}),
            std::runtime_error);
    std::map<ParticleType, std::pair<std::vector<double>, std::vector<double>>> bad_length;
    bad_length[ParticleType::NuMu] = {{1.0, 3.0}, {2e-38}};
    EXPECT_THROW(TrivialCrossSection(bad_length, {ParticleType::Nucleon}),
            std::runtime_error);
}

TEST(TrivialCrossSection, SignatureFiltering) {
    std::shared_ptr<TrivialCrossSection> xs = MakeTabulated();

    std::vector<InteractionSignature> signatures =
        xs->GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::Nucleon);
    ASSERT_EQ(1u, signatures.size());
    EXPECT_EQ(ParticleType::NuMu, signatures[0].primary_type);
    EXPECT_EQ(ParticleType::Nucleon, signatures[0].target_type);
    ASSERT_EQ(2u, signatures[0].secondary_types.size());
    EXPECT_EQ(ParticleType::NuMu, signatures[0].secondary_types[0]);
    EXPECT_EQ(ParticleType::Nucleon, signatures[0].secondary_types[1]);

    // Unconfigured primaries and targets produce no signatures.
    EXPECT_TRUE(xs->GetPossibleSignaturesFromParents(ParticleType::NuE, ParticleType::Nucleon).empty());
    EXPECT_TRUE(xs->GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::PPlus).empty());
    EXPECT_TRUE(xs->GetPossibleTargetsFromPrimary(ParticleType::NuE).empty());

    // The nucleus-target convention is expressible.
    TrivialCrossSection argon(1e-38, {ParticleType::NuMu}, {ParticleType::Ar40Nucleus});
    std::vector<InteractionSignature> argon_signatures =
        argon.GetPossibleSignaturesFromParents(ParticleType::NuMu, ParticleType::Ar40Nucleus);
    ASSERT_EQ(1u, argon_signatures.size());
    EXPECT_EQ(ParticleType::Ar40Nucleus, argon_signatures[0].target_type);
}

TEST(TrivialCrossSection, FinalStateProbability) {
    std::shared_ptr<TrivialCrossSection> xs = MakeTabulated();

    // Matching signature: exactly unit probability.
    EXPECT_DOUBLE_EQ(1.0, xs->FinalStateProbability(MakeRecord(2.0)));

    // A record with different secondaries has zero probability.
    InteractionRecord mismatched = MakeRecord(2.0);
    mismatched.signature.secondary_types = {ParticleType::MuMinus, ParticleType::Hadrons};
    EXPECT_DOUBLE_EQ(0.0, xs->FinalStateProbability(mismatched));

    // Below threshold the total vanishes, so the probability is zero.
    EXPECT_DOUBLE_EQ(0.0, xs->FinalStateProbability(MakeRecord(0.5)));
}

TEST(TrivialCrossSection, PassThroughKinematics) {
    std::shared_ptr<TrivialCrossSection> xs = MakeTabulated();
    std::shared_ptr<siren::utilities::SIREN_random> random =
        std::make_shared<siren::utilities::SIREN_random>(1234);

    InteractionRecord record = MakeRecord(2.0);
    record.primary_momentum = {2.0, 0.3, -0.4, std::sqrt(2.0 * 2.0 - 0.3 * 0.3 - 0.4 * 0.4)};

    CrossSectionDistributionRecord xsec_record(record);
    xs->SampleFinalState(xsec_record, random);

    InteractionRecord finalized = record;
    xsec_record.Finalize(finalized);

    ASSERT_EQ(2u, finalized.secondary_momenta.size());
    for(size_t component = 0; component < 4; ++component) {
        EXPECT_DOUBLE_EQ(record.primary_momentum[component],
                finalized.secondary_momenta[0][component]);
    }
    EXPECT_DOUBLE_EQ(record.target_mass, finalized.secondary_momenta[1][0]);
    EXPECT_DOUBLE_EQ(0.0, finalized.secondary_momenta[1][1]);
    EXPECT_DOUBLE_EQ(0.0, finalized.secondary_momenta[1][2]);
    EXPECT_DOUBLE_EQ(0.0, finalized.secondary_momenta[1][3]);
}

TEST(TrivialCrossSection, EqualAndSerialization) {
    std::shared_ptr<TrivialCrossSection> xs = MakeTabulated();
    std::shared_ptr<TrivialCrossSection> same = MakeTabulated();
    TrivialCrossSection constant(1e-38, {ParticleType::NuMu}, {ParticleType::Nucleon});

    EXPECT_TRUE(xs->equal(*same));
    EXPECT_FALSE(xs->equal(constant));

    std::stringstream ss;
    {
        cereal::JSONOutputArchive archive(ss);
        std::shared_ptr<CrossSection> base = xs;
        archive(base);
    }
    std::shared_ptr<CrossSection> loaded;
    {
        cereal::JSONInputArchive archive(ss);
        archive(loaded);
    }
    ASSERT_TRUE(loaded);
    EXPECT_TRUE(xs->equal(*loaded));
    EXPECT_DOUBLE_EQ(4e-38, loaded->TotalCrossSection(MakeRecord(2.0)));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
