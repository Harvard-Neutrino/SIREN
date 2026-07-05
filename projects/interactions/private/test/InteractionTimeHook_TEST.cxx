// Pins the post-selection vertex-time hooks Decay::SampleDecayTime and
// CrossSection::SampleInteractionTime: the base hooks are the identity, and an
// override applied through SetInteractionTime shifts record.interaction_time
// and back-syncs the daughter production times. The wiring here mirrors
// Injector::SampleCrossSection (SampleFinalState -> hook -> SetInteractionTime
// if it differs -> Finalize).

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Random.h"

using namespace siren::interactions;
using namespace siren::dataclasses;

namespace {

// Value the hooks shift the vertex time to, distinct from the flight time.
constexpr double kShiftedTime = 123.5;

// Populate each secondary with a valid four-momentum so Finalize can derive the
// remaining fields. A real SampleFinalState fills these in.
void FillSecondaries(CrossSectionDistributionRecord & record) {
    for(auto & secondary : record.GetSecondaryParticleRecords()) {
        secondary.SetFourMomentum({5.0, 0.0, 0.0, 4.0});
        secondary.SetHelicity(0.0);
    }
}

// Minimal CrossSection whose SampleInteractionTime overrides the vertex time.
class TimeShiftXS : public CrossSection {
public:
    double TotalCrossSection(InteractionRecord const &) const override { return 1.0; }
    double DifferentialCrossSection(InteractionRecord const &) const override { return 0.0; }
    double InteractionThreshold(InteractionRecord const &) const override { return 0.0; }
    double FinalStateProbability(InteractionRecord const &) const override { return 0.0; }
    std::vector<std::string> DensityVariables() const override { return {}; }
    std::vector<ParticleType> GetPossibleTargets() const override { return { ParticleType::PPlus }; }
    std::vector<ParticleType> GetPossibleTargetsFromPrimary(ParticleType) const override { return { ParticleType::PPlus }; }
    std::vector<ParticleType> GetPossiblePrimaries() const override { return { ParticleType::NuMu }; }
    std::vector<InteractionSignature> GetPossibleSignatures() const override { return {}; }
    std::vector<InteractionSignature> GetPossibleSignaturesFromParents(ParticleType, ParticleType) const override { return {}; }
    void SampleFinalState(CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random>) const override { FillSecondaries(record); }
    bool equal(CrossSection const &) const override { return false; }

    double SampleInteractionTime(CrossSectionDistributionRecord const &, std::shared_ptr<siren::utilities::SIREN_random>) const override {
        return kShiftedTime;
    }
};

// Minimal Decay whose SampleDecayTime overrides the vertex time.
class TimeShiftDecay : public Decay {
public:
    double TotalDecayWidthAllFinalStates(InteractionRecord const &) const override { return 1.0; }
    double TotalDecayWidth(ParticleType) const override { return 1.0; }
    double TotalDecayWidth(InteractionRecord const &) const override { return 1.0; }
    double DifferentialDecayWidth(InteractionRecord const &) const override { return 0.0; }
    double FinalStateProbability(InteractionRecord const &) const override { return 0.0; }
    std::vector<std::string> DensityVariables() const override { return {}; }
    std::vector<InteractionSignature> GetPossibleSignatures() const override { return {}; }
    std::vector<InteractionSignature> GetPossibleSignaturesFromParent(ParticleType) const override { return {}; }
    void SampleFinalState(CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random>) const override { FillSecondaries(record); }
    bool equal(Decay const &) const override { return false; }

    double SampleDecayTime(CrossSectionDistributionRecord const &, std::shared_ptr<siren::utilities::SIREN_random>) const override {
        return kShiftedTime;
    }
};

// Build a parent record with a fixed vertex time and one secondary.
InteractionRecord MakeVertexRecord() {
    InteractionRecord r;
    r.signature.primary_type = ParticleType::NuMu;
    r.signature.target_type = ParticleType::PPlus;
    r.signature.secondary_types = { ParticleType::MuMinus };
    r.primary_momentum = { 10.0, 0.0, 0.0, 10.0 };
    r.interaction_vertex = { 1.0, 2.0, 3.0 };
    r.interaction_time = 17.0;
    return r;
}

// Replicate the Injector's post-selection wiring: sample the final state, call
// the hook, apply the override if it differs, then finalize.
template <typename Interaction>
InteractionRecord DriveVertex(Interaction const & interaction,
                              std::shared_ptr<siren::utilities::SIREN_random> random,
                              double (*hook)(Interaction const &, CrossSectionDistributionRecord const &, std::shared_ptr<siren::utilities::SIREN_random>)) {
    InteractionRecord record = MakeVertexRecord();
    CrossSectionDistributionRecord xsec_record(record);
    interaction.SampleFinalState(xsec_record, random);
    double proposed_time = hook(interaction, xsec_record, random);
    if(proposed_time != xsec_record.GetInteractionTime()) {
        xsec_record.SetInteractionTime(proposed_time);
    }
    xsec_record.Finalize(record);
    return record;
}

double xs_hook(TimeShiftXS const & xs, CrossSectionDistributionRecord const & r, std::shared_ptr<siren::utilities::SIREN_random> rand) {
    return xs.SampleInteractionTime(r, rand);
}
double decay_hook(TimeShiftDecay const & d, CrossSectionDistributionRecord const & r, std::shared_ptr<siren::utilities::SIREN_random> rand) {
    return d.SampleDecayTime(r, rand);
}

} // namespace

// The base-class hooks are the identity: they hand back the flight-time value
// already on the record, so default physics sees no change.
TEST(InteractionTimeHook, DefaultXSHookIsIdentity) {
    TimeShiftXS xs;
    InteractionRecord record = MakeVertexRecord();
    CrossSectionDistributionRecord xsec_record(record);
    // Call the base method explicitly to exercise the identity default.
    double t = xs.CrossSection::SampleInteractionTime(xsec_record, nullptr);
    EXPECT_DOUBLE_EQ(t, xsec_record.GetInteractionTime());
}

TEST(InteractionTimeHook, DefaultDecayHookIsIdentity) {
    TimeShiftDecay d;
    InteractionRecord record = MakeVertexRecord();
    CrossSectionDistributionRecord xsec_record(record);
    double t = d.Decay::SampleDecayTime(xsec_record, nullptr);
    EXPECT_DOUBLE_EQ(t, xsec_record.GetInteractionTime());
}

// An overriding cross section shifts the vertex time, and the daughter picks it
// up as its production time (the back-sync through secondary_times).
TEST(InteractionTimeHook, CrossSectionOverrideShiftsVertexAndDaughterTime) {
    TimeShiftXS xs;
    auto random = std::make_shared<siren::utilities::SIREN_random>(1234);
    InteractionRecord record = DriveVertex<TimeShiftXS>(xs, random, xs_hook);
    EXPECT_DOUBLE_EQ(record.interaction_time, kShiftedTime);
    ASSERT_EQ(record.secondary_times.size(), 1u);
    EXPECT_DOUBLE_EQ(record.secondary_times.at(0), kShiftedTime);
}

TEST(InteractionTimeHook, DecayOverrideShiftsVertexAndDaughterTime) {
    TimeShiftDecay d;
    auto random = std::make_shared<siren::utilities::SIREN_random>(1234);
    InteractionRecord record = DriveVertex<TimeShiftDecay>(d, random, decay_hook);
    EXPECT_DOUBLE_EQ(record.interaction_time, kShiftedTime);
    ASSERT_EQ(record.secondary_times.size(), 1u);
    EXPECT_DOUBLE_EQ(record.secondary_times.at(0), kShiftedTime);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
