#include <gtest/gtest.h>

#include "PhaseSpaceTestFixtures.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/ConstantDensityDistribution.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/FinalStateProbability.h"
#include "SIREN/injection/Weighter.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>

namespace {

using siren::dataclasses::InteractionRecord;
using siren::injection::MultiChannelPhaseSpace;
using siren::injection::PhaseSpaceChannel;
using siren::injection::PhaseSpaceMeasure;
using siren::injection::PhaseSpaceTopology;
using siren::injection::test::ConstantChannel;

class CompletePrimaryDistribution final
    : public siren::distributions::VertexPositionDistribution {
private:
    std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        siren::dataclasses::PrimaryDistributionRecord &) const override
    {
        return {siren::math::Vector3D(0.0, 0.0, 0.0),
                siren::math::Vector3D(0.0, 0.0, 0.0)};
    }

public:
    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        siren::dataclasses::PrimaryDistributionRecord & record) const override
    {
        record.SetMass(0.0);
        record.SetFourMomentum({1.0, 0.0, 0.0, 1.0});
        record.SetInitialPosition({0.0, 0.0, 0.0});
        record.SetInteractionVertex({0.0, 0.0, 0.0});
        record.SetHelicity(0.0);
        record.SetInitialTime(0.0);
        record.SetInteractionTime(0.0);
    }

    double GenerationProbability(
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        InteractionRecord const &) const override
    {
        return 1.0;
    }

    std::string Name() const override { return "CompletePrimary"; }

    std::shared_ptr<siren::distributions::PrimaryInjectionDistribution>
    clone() const override
    {
        return std::make_shared<CompletePrimaryDistribution>();
    }

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(
        std::shared_ptr<siren::detector::DetectorModel const>,
        std::shared_ptr<siren::interactions::InteractionCollection const>,
        InteractionRecord const &) const override
    {
        return {siren::math::Vector3D(0.0, 0.0, 0.0),
                siren::math::Vector3D(0.0, 0.0, 0.0)};
    }

protected:
    bool equal(
        siren::distributions::WeightableDistribution const & other) const override
    {
        return dynamic_cast<CompletePrimaryDistribution const *>(&other) != nullptr;
    }

    bool less(
        siren::distributions::WeightableDistribution const &) const override
    {
        return false;
    }
};

class RetryableFailureInjector final : public siren::injection::Injector {
public:
    RetryableFailureInjector(
        unsigned int attempts,
        std::shared_ptr<siren::injection::PrimaryInjectionProcess> process,
        std::shared_ptr<siren::utilities::SIREN_random> random)
        : Injector(attempts, nullptr, std::move(process), std::move(random)) {}

    std::shared_ptr<siren::interactions::Interaction> SelectChannel(
        InteractionRecord &,
        std::shared_ptr<siren::interactions::InteractionCollection>)
        const override
    {
        throw siren::utilities::InjectionFailure(
            "test event has no kinematically allowed phase space");
    }
};

siren::dataclasses::InteractionSignature SharedDecaySignature() {
    siren::dataclasses::InteractionSignature signature;
    signature.primary_type = siren::dataclasses::ParticleType::NuMu;
    signature.target_type = siren::dataclasses::ParticleType::Decay;
    signature.secondary_types = {
        siren::dataclasses::ParticleType::NuE,
        siren::dataclasses::ParticleType::NuEBar};
    return signature;
}

siren::dataclasses::InteractionSignature SharedCrossSectionSignature() {
    siren::dataclasses::InteractionSignature signature;
    signature.primary_type = siren::dataclasses::ParticleType::NuMu;
    signature.target_type = siren::dataclasses::ParticleType::Nucleon;
    signature.secondary_types = {
        siren::dataclasses::ParticleType::NuMu,
        siren::dataclasses::ParticleType::Nucleon};
    return signature;
}

class TaggedDecay : public siren::interactions::Decay {
public:
    TaggedDecay(int tag, double rate, double density)
        : tag_(tag), rate_(rate), density_(density) {}

    bool equal(siren::interactions::Decay const & other) const override {
        auto const * tagged = dynamic_cast<TaggedDecay const *>(&other);
        return tagged && tagged->tag_ == tag_ && tagged->rate_ == rate_
            && tagged->density_ == density_;
    }
    double TotalDecayWidthAllFinalStates(InteractionRecord const &) const override {
        return rate_;
    }
    double TotalDecayWidth(siren::dataclasses::ParticleType) const override {
        return rate_;
    }
    double TotalDecayWidth(InteractionRecord const &) const override {
        return rate_;
    }
    double TotalDecayLengthAllFinalStates(InteractionRecord const &) const override {
        return siren::utilities::Constants::cm / rate_;
    }
    double TotalDecayLength(InteractionRecord const &) const override {
        return siren::utilities::Constants::cm / rate_;
    }
    double DifferentialDecayWidth(InteractionRecord const &) const override {
        return density_ * rate_;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord & record,
        std::shared_ptr<siren::utilities::SIREN_random>) const override
    {
        record.SetInteractionParameter("selected_process", tag_);
        for (std::size_t i = 0;
             i < record.GetSecondaryParticleRecords().size(); ++i) {
            auto & secondary = record.GetSecondaryParticleRecord(i);
            secondary.SetMass(0.0);
            secondary.SetFourMomentum(
                i == 0 ? std::array<double, 4>{1.0, 0.0, 0.0, 1.0}
                       : std::array<double, 4>{1.0, 0.0, 0.0, -1.0});
            secondary.SetHelicity(static_cast<double>(tag_));
        }
    }
    double SampleDecayTime(
        siren::dataclasses::CrossSectionDistributionRecord const &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override
    {
        return 100.0 + tag_;
    }
    std::vector<double> SecondaryMasses(
        std::vector<siren::dataclasses::ParticleType> const & types) const override
    {
        return std::vector<double>(types.size(), 0.0);
    }
    std::vector<double> SecondaryHelicities(
        InteractionRecord const & record) const override
    {
        return std::vector<double>(
            record.signature.secondary_types.size(),
            static_cast<double>(tag_));
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return {SharedDecaySignature()};
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParent(
        siren::dataclasses::ParticleType primary) const override
    {
        return primary == SharedDecaySignature().primary_type
            ? std::vector<siren::dataclasses::InteractionSignature>{
                  SharedDecaySignature()}
            : std::vector<siren::dataclasses::InteractionSignature>{};
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return density_;
    }
    std::vector<std::string> DensityVariables() const override {
        return {"cos_theta"};
    }
    siren::dataclasses::PhaseSpaceTopology Topology() const override {
        return siren::dataclasses::PhaseSpaceTopology::Decay2Body;
    }
    siren::dataclasses::PhaseSpaceMeasure Measure() const override {
        return siren::dataclasses::PhaseSpaceMeasure::SolidAngleRest();
    }
    int Tag() const { return tag_; }

private:
    int tag_;
    double rate_;
    double density_;
};

class TaggedCrossSection final : public siren::interactions::CrossSection {
public:
    TaggedCrossSection(int tag, double rate, double density)
        : tag_(tag), rate_(rate), density_(density) {}

    bool equal(siren::interactions::CrossSection const & other) const override {
        auto const * tagged = dynamic_cast<TaggedCrossSection const *>(&other);
        return tagged && tagged->tag_ == tag_ && tagged->rate_ == rate_
            && tagged->density_ == density_;
    }
    double TotalCrossSection(InteractionRecord const &) const override {
        return rate_;
    }
    double DifferentialCrossSection(InteractionRecord const &) const override {
        return density_ * rate_;
    }
    double InteractionThreshold(InteractionRecord const &) const override {
        return 0.0;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord & record,
        std::shared_ptr<siren::utilities::SIREN_random>) const override
    {
        record.SetInteractionParameter("selected_process", tag_);
        for (std::size_t i = 0;
             i < record.GetSecondaryParticleRecords().size(); ++i) {
            auto & secondary = record.GetSecondaryParticleRecord(i);
            secondary.SetMass(0.0);
            secondary.SetFourMomentum(
                i == 0 ? std::array<double, 4>{1.0, 0.0, 0.0, 1.0}
                       : std::array<double, 4>{1.0, 0.0, 0.0, -1.0});
            secondary.SetHelicity(static_cast<double>(tag_));
        }
    }
    double SampleInteractionTime(
        siren::dataclasses::CrossSectionDistributionRecord const &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override
    {
        return 100.0 + tag_;
    }
    std::vector<double> SecondaryMasses(
        std::vector<siren::dataclasses::ParticleType> const & types) const override
    {
        return std::vector<double>(types.size(), 0.0);
    }
    std::vector<double> SecondaryHelicities(
        InteractionRecord const & record) const override
    {
        return std::vector<double>(
            record.signature.secondary_types.size(),
            static_cast<double>(tag_));
    }
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        return {siren::dataclasses::ParticleType::Nucleon};
    }
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(
        siren::dataclasses::ParticleType primary) const override
    {
        return primary == SharedCrossSectionSignature().primary_type
            ? GetPossibleTargets()
            : std::vector<siren::dataclasses::ParticleType>{};
    }
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        return {siren::dataclasses::ParticleType::NuMu};
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return {SharedCrossSectionSignature()};
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParents(
        siren::dataclasses::ParticleType primary,
        siren::dataclasses::ParticleType target) const override
    {
        auto signature = SharedCrossSectionSignature();
        return primary == signature.primary_type && target == signature.target_type
            ? std::vector<siren::dataclasses::InteractionSignature>{signature}
            : std::vector<siren::dataclasses::InteractionSignature>{};
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return density_;
    }
    std::vector<std::string> DensityVariables() const override {
        return {"Q2"};
    }
    siren::dataclasses::PhaseSpaceTopology Topology() const override {
        return siren::dataclasses::PhaseSpaceTopology::Scatter2to2;
    }
    siren::dataclasses::PhaseSpaceMeasure Measure() const override {
        return siren::dataclasses::PhaseSpaceMeasure::MandelstamQ2();
    }
    int Tag() const { return tag_; }

private:
    int tag_;
    double rate_;
    double density_;
};

std::shared_ptr<siren::detector::DetectorModel> SelectionDetector() {
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    siren::detector::DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = siren::geometry::Sphere(100.0, 0.0).create();
    world.density = siren::detector::ConstantDensityDistribution(1.0).create();
    detector->AddSector(world);
    return detector;
}

InteractionRecord SelectionRecord() {
    InteractionRecord record;
    record.signature.primary_type = siren::dataclasses::ParticleType::NuMu;
    record.primary_mass = 0.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 2.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.interaction_time = 7.0;
    return record;
}

std::shared_ptr<siren::injection::PrimaryInjectionProcess> SelectionProcess(
    std::shared_ptr<siren::interactions::InteractionCollection> interactions)
{
    auto process = std::make_shared<siren::injection::PrimaryInjectionProcess>(
        siren::dataclasses::ParticleType::NuMu, std::move(interactions));
    process->AddPrimaryInjectionDistribution(
        std::make_shared<CompletePrimaryDistribution>());
    return process;
}

TEST(ConcreteInteractionSelection, SameSignatureDecaysKeepRateProbability) {
    auto detector = SelectionDetector();
    for (bool reverse_order : {false, true}) {
        auto slow = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
        auto fast = std::make_shared<TaggedDecay>(2, 3.0, 6.0);
        std::vector<std::shared_ptr<siren::interactions::Decay>> models =
            reverse_order
                ? std::vector<std::shared_ptr<siren::interactions::Decay>>{
                      fast, slow}
                : std::vector<std::shared_ptr<siren::interactions::Decay>>{
                      slow, fast};
        auto interactions =
            std::make_shared<siren::interactions::InteractionCollection>(
                siren::dataclasses::ParticleType::NuMu, models);
        auto random = std::make_shared<siren::utilities::SIREN_random>(8675309);
        siren::injection::Injector injector(
            1, detector, SelectionProcess(interactions), random);

        int slow_count = 0;
        for (int i = 0; i < 400; ++i) {
            InteractionRecord record = SelectionRecord();
            auto selected = injector.SelectChannel(record, interactions);
            ASSERT_TRUE(selected);
            slow_count += selected.get() == slow.get();
        }
        EXPECT_GT(slow_count, 60) << "reverse_order=" << reverse_order;
        EXPECT_LT(slow_count, 140) << "reverse_order=" << reverse_order;
    }
}

TEST(ConcreteInteractionSelection, PureDecaysNeedNoMaterialWorld) {
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    auto decay = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
    auto interactions =
        std::make_shared<siren::interactions::InteractionCollection>(
            siren::dataclasses::ParticleType::NuMu,
            std::vector<std::shared_ptr<siren::interactions::Decay>>{decay});
    auto random = std::make_shared<siren::utilities::SIREN_random>(831609);
    siren::injection::Injector injector(
        1, detector, SelectionProcess(interactions), random);

    for (double distance : {0.0, 1e6, 1e12, 1e18}) {
        auto record = SelectionRecord();
        record.interaction_vertex = {distance, distance, distance};
        auto selected = injector.SelectChannel(record, interactions);
        EXPECT_EQ(selected.get(), decay.get());
        EXPECT_EQ(record.signature, SharedDecaySignature());
    }
}

TEST(ConcreteInteractionSelection, OnlyCrossSectionsNavigateMaterial) {
    class CountingSphere : public siren::geometry::Sphere {
    public:
        CountingSphere() : Sphere(100.0, 0.0) {}
        mutable int intersection_calls = 0;
        std::vector<Intersection> ComputeIntersections(
            siren::math::Vector3D const & position,
            siren::math::Vector3D const & direction) const override
        {
            ++intersection_calls;
            return Sphere::ComputeIntersections(position, direction);
        }
    };

    auto geometry = std::make_shared<CountingSphere>();
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    siren::detector::DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = geometry;
    world.density = siren::detector::ConstantDensityDistribution(1.0).create();
    detector->AddSector(world);

    auto decay = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
    auto cross_section = std::make_shared<TaggedCrossSection>(2, 1.0, 2.0);
    for (bool include_scattering : {false, true}) {
        std::vector<std::shared_ptr<siren::interactions::CrossSection>> scattering;
        if (include_scattering) scattering.push_back(cross_section);
        auto interactions =
            std::make_shared<siren::interactions::InteractionCollection>(
                siren::dataclasses::ParticleType::NuMu, scattering,
                std::vector<std::shared_ptr<siren::interactions::Decay>>{decay});
        auto random = std::make_shared<siren::utilities::SIREN_random>(831609);
        siren::injection::Injector injector(
            1, detector, SelectionProcess(interactions), random);
        auto record = SelectionRecord();
        geometry->intersection_calls = 0;
        auto selected = injector.SelectChannel(record, interactions);
        ASSERT_TRUE(selected);
        if (include_scattering) {
            EXPECT_GT(geometry->intersection_calls, 0);
        } else {
            EXPECT_EQ(geometry->intersection_calls, 0);
            EXPECT_EQ(selected.get(), decay.get());
        }
    }
}

TEST(ConcreteInteractionSelection, SameSignatureCrossSectionsKeepRateProbability) {
    auto detector = SelectionDetector();
    for (bool reverse_order : {false, true}) {
        auto slow = std::make_shared<TaggedCrossSection>(1, 1.0, 2.0);
        auto fast = std::make_shared<TaggedCrossSection>(2, 3.0, 6.0);
        std::vector<std::shared_ptr<siren::interactions::CrossSection>> models =
            reverse_order
                ? std::vector<std::shared_ptr<siren::interactions::CrossSection>>{
                      fast, slow}
                : std::vector<std::shared_ptr<siren::interactions::CrossSection>>{
                      slow, fast};
        auto interactions =
            std::make_shared<siren::interactions::InteractionCollection>(
                siren::dataclasses::ParticleType::NuMu, models);
        auto random = std::make_shared<siren::utilities::SIREN_random>(314159);
        siren::injection::Injector injector(
            1, detector, SelectionProcess(interactions), random);

        int slow_count = 0;
        for (int i = 0; i < 400; ++i) {
            InteractionRecord record = SelectionRecord();
            auto selected = injector.SelectChannel(record, interactions);
            ASSERT_TRUE(selected);
            slow_count += selected.get() == slow.get();
        }
        EXPECT_GT(slow_count, 60) << "reverse_order=" << reverse_order;
        EXPECT_LT(slow_count, 140) << "reverse_order=" << reverse_order;
    }
}

std::pair<int, int> GenerateTaggedDecayCounts(bool use_phase_space) {
    constexpr int sample_count = 240;
    auto detector = SelectionDetector();
    auto slow = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
    auto fast = std::make_shared<TaggedDecay>(2, 3.0, 6.0);
    auto interactions =
        std::make_shared<siren::interactions::InteractionCollection>(
            siren::dataclasses::ParticleType::NuMu,
            std::vector<std::shared_ptr<siren::interactions::Decay>>{
                slow, fast});
    auto process = SelectionProcess(interactions);
    if (use_phase_space) {
        auto phase_space = std::make_shared<MultiChannelPhaseSpace>(
            std::vector<std::shared_ptr<PhaseSpaceChannel>>{std::make_shared<ConstantChannel>()});
        process->SetPhaseSpace(SharedDecaySignature(), phase_space);
    }
    auto random = std::make_shared<siren::utilities::SIREN_random>(271828);
    siren::injection::Injector injector(
        sample_count, detector, process, random);

    int slow_count = 0;
    int fast_count = 0;
    for (int i = 0; i < sample_count; ++i) {
        siren::dataclasses::InteractionTree tree = injector.GenerateEvent();
        EXPECT_EQ(tree.tree.size(), 1u);
        if (tree.tree.empty()) continue;
        InteractionRecord const & record = tree.tree.front()->record;
        int tag = use_phase_space
            ? static_cast<int>(record.secondary_helicities.at(0))
            : static_cast<int>(
                  record.interaction_parameters.at("selected_process"));
        EXPECT_DOUBLE_EQ(record.interaction_time, 100.0 + tag);
        EXPECT_EQ(record.secondary_times.size(), 2u);
        if (record.secondary_times.size() == 2u) {
            EXPECT_DOUBLE_EQ(record.secondary_times[0], record.interaction_time);
            EXPECT_DOUBLE_EQ(record.secondary_times[1], record.interaction_time);
        }
        slow_count += tag == 1;
        fast_count += tag == 2;
    }
    return {slow_count, fast_count};
}

TEST(ConcreteInteractionSelection, FallbackSamplesTheSelectedDecay) {
    auto [slow_count, fast_count] = GenerateTaggedDecayCounts(false);
    EXPECT_GT(slow_count, 30);
    EXPECT_LT(slow_count, 90);
    EXPECT_EQ(slow_count + fast_count, 240);
}

TEST(ConcreteInteractionSelection, PhaseSpaceUsesSelectedDecayMetadataAndTime) {
    auto [slow_count, fast_count] = GenerateTaggedDecayCounts(true);
    EXPECT_GT(slow_count, 30);
    EXPECT_LT(slow_count, 90);
    EXPECT_EQ(slow_count + fast_count, 240);
}

std::pair<int, int> GenerateTaggedCrossSectionCounts(bool use_phase_space) {
    constexpr int sample_count = 240;
    auto detector = SelectionDetector();
    auto slow = std::make_shared<TaggedCrossSection>(1, 1.0, 2.0);
    auto fast = std::make_shared<TaggedCrossSection>(2, 3.0, 6.0);
    auto interactions =
        std::make_shared<siren::interactions::InteractionCollection>(
            siren::dataclasses::ParticleType::NuMu,
            std::vector<std::shared_ptr<siren::interactions::CrossSection>>{
                slow, fast});
    auto process = SelectionProcess(interactions);
    if (use_phase_space) {
        auto phase_space = std::make_shared<MultiChannelPhaseSpace>(
            std::vector<std::shared_ptr<PhaseSpaceChannel>>{std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())});
        process->SetPhaseSpace(SharedCrossSectionSignature(), phase_space);
    }
    auto random = std::make_shared<siren::utilities::SIREN_random>(161803);
    siren::injection::Injector injector(
        sample_count, detector, process, random);

    int slow_count = 0;
    int fast_count = 0;
    for (int i = 0; i < sample_count; ++i) {
        siren::dataclasses::InteractionTree tree = injector.GenerateEvent();
        EXPECT_EQ(tree.tree.size(), 1u);
        if (tree.tree.empty()) continue;
        InteractionRecord const & record = tree.tree.front()->record;
        int tag = use_phase_space
            ? static_cast<int>(record.secondary_helicities.at(0))
            : static_cast<int>(
                  record.interaction_parameters.at("selected_process"));
        EXPECT_DOUBLE_EQ(record.interaction_time, 100.0 + tag);
        slow_count += tag == 1;
        fast_count += tag == 2;
    }
    return {slow_count, fast_count};
}

TEST(ConcreteInteractionSelection, CrossSectionFallbackUsesSelectedModel) {
    auto [slow_count, fast_count] = GenerateTaggedCrossSectionCounts(false);
    EXPECT_GT(slow_count, 30);
    EXPECT_LT(slow_count, 90);
    EXPECT_EQ(slow_count + fast_count, 240);
}

TEST(ConcreteInteractionSelection, CrossSectionPhaseSpaceUsesSelectedMetadata) {
    auto [slow_count, fast_count] = GenerateTaggedCrossSectionCounts(true);
    EXPECT_GT(slow_count, 30);
    EXPECT_LT(slow_count, 90);
    EXPECT_EQ(slow_count + fast_count, 240);
}

TEST(ConcreteInteractionSelection, InvalidSelectedInteractionFailsLoudly) {
    auto detector = SelectionDetector();
    auto model = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
    auto interactions =
        std::make_shared<siren::interactions::InteractionCollection>(
            siren::dataclasses::ParticleType::NuMu,
            std::vector<std::shared_ptr<siren::interactions::Decay>>{model});
    auto random = std::make_shared<siren::utilities::SIREN_random>(42);
    siren::injection::Injector injector(
        1, detector, SelectionProcess(interactions), random);
    InteractionRecord record = SelectionRecord();
    record.signature = SharedDecaySignature();

    EXPECT_THROW(
        injector.SampleMatchingFinalState(record, nullptr),
        siren::utilities::ConfigurationError);
    EXPECT_THROW(
        injector.SampleMatchingFinalState(
            record,
            std::make_shared<siren::interactions::Interaction>()),
        siren::utilities::ConfigurationError);
}

TEST(ConcreteInteractionSelection, FixedDensityIsRateConditionalAndOrderIndependent) {
    auto detector = SelectionDetector();
    for (bool reverse_order : {false, true}) {
        auto slow = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
        auto fast = std::make_shared<TaggedDecay>(2, 3.0, 6.0);
        std::vector<std::shared_ptr<siren::interactions::Decay>> models =
            reverse_order
                ? std::vector<std::shared_ptr<siren::interactions::Decay>>{
                      fast, slow}
                : std::vector<std::shared_ptr<siren::interactions::Decay>>{
                      slow, fast};
        auto interactions =
            std::make_shared<siren::interactions::InteractionCollection>(
                siren::dataclasses::ParticleType::NuMu, models);
        InteractionRecord record = SelectionRecord();
        record.signature = SharedDecaySignature();

        siren::injection::PhysicalProcess process(
            siren::dataclasses::ParticleType::NuMu, interactions);
        process.SetWeightingMode(siren::dataclasses::VertexWeightingMode::Fixed());
        double density = siren::injection::ProcessFinalStateProbability(
            detector, process, record,
            siren::injection::PhaseSpaceConvention{
                PhaseSpaceTopology::Decay2Body,
                PhaseSpaceMeasure::SolidAngleRest()});
        EXPECT_DOUBLE_EQ(density, 5.0);
        EXPECT_DOUBLE_EQ(
            siren::injection::FixedVertexChannelSelectionProbability(
                detector, interactions, record),
            1.0);
    }
}

TEST(KinematicInjectionFailure, InjectorCountsFailureAsAnAttempt) {
    auto process = std::make_shared<siren::injection::PrimaryInjectionProcess>(
        siren::dataclasses::ParticleType::unknown, nullptr);
    process->AddPrimaryInjectionDistribution(
        std::make_shared<CompletePrimaryDistribution>());
    auto random = std::make_shared<siren::utilities::SIREN_random>(244949);
    RetryableFailureInjector injector(2, process, random);

    siren::dataclasses::InteractionTree first;
    EXPECT_NO_THROW(first = injector.GenerateEvent());
    EXPECT_TRUE(first.tree.empty());
    EXPECT_EQ(injector.InjectionAttempts(), 1u);
    EXPECT_EQ(injector.InjectedEvents(), 0u);

    siren::dataclasses::InteractionTree second;
    EXPECT_NO_THROW(second = injector.GenerateEvent());
    EXPECT_TRUE(second.tree.empty());
    EXPECT_EQ(injector.InjectionAttempts(), 2u);
    EXPECT_EQ(injector.InjectedEvents(), 0u);

    // The attempt budget is still enforced after two rejected draws.
    EXPECT_THROW(injector.GenerateEvent(), std::runtime_error);
}

TEST(CoreAtRest, WidthRatiosRemainFiniteForTinyWidths) {
    class RestDecay : public TaggedDecay {
    public:
        using TaggedDecay::TaggedDecay;
        double TotalDecayLength(InteractionRecord const &) const override { return 0.0; }
    };
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    auto slow = std::make_shared<RestDecay>(1, 1e-300, 2.0);
    auto fast = std::make_shared<RestDecay>(2, 3e-300, 6.0);
    auto interactions = std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::NuMu,
        std::vector<std::shared_ptr<siren::interactions::Decay>>{slow, fast});
    auto random = std::make_shared<siren::utilities::SIREN_random>(4187);
    siren::injection::Injector injector(400, detector, SelectionProcess(interactions), random);
    auto record = SelectionRecord();
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    int slow_count = 0;
    for (int i = 0; i < 400; ++i)
        slow_count += injector.SelectChannel(record, interactions).get() == slow.get();
    EXPECT_GT(slow_count, 60);
    EXPECT_LT(slow_count, 140);
    EXPECT_DOUBLE_EQ(siren::injection::CrossSectionProbability(detector, interactions, record), 5.0);
}

TEST(CoreAtRest, MaterialSelectionFailsBeforeGeometryNavigation) {
    auto detector = SelectionDetector();
    auto xs = std::make_shared<TaggedCrossSection>(1, 1.0, 1.0);
    auto interactions = std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::NuMu,
        std::vector<std::shared_ptr<siren::interactions::CrossSection>>{xs});
    siren::injection::Injector injector(1, detector, SelectionProcess(interactions),
        std::make_shared<siren::utilities::SIREN_random>(1));
    auto record = SelectionRecord();
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    EXPECT_THROW(injector.SelectChannel(record, interactions), siren::utilities::ConfigurationError);
}

TEST(CoreProbability, RecordOverloadUsesConfiguredProposal) {
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    auto decay = std::make_shared<TaggedDecay>(1, 1.0, 2.0);
    auto interactions = std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::NuMu,
        std::vector<std::shared_ptr<siren::interactions::Decay>>{decay});
    auto process = SelectionProcess(interactions);
    process->SetPhaseSpace(SharedDecaySignature(), std::make_shared<MultiChannelPhaseSpace>(
        std::vector<std::shared_ptr<PhaseSpaceChannel>>{std::make_shared<ConstantChannel>(0.25)}));
    siren::injection::Injector injector(5, detector, process,
        std::make_shared<siren::utilities::SIREN_random>(1));
    auto record = SelectionRecord();
    record.signature = SharedDecaySignature();
    EXPECT_DOUBLE_EQ(injector.GenerationProbability(record, process), 0.25);
    EXPECT_DOUBLE_EQ(injector.GenerationProbability(record), 1.25);
}

} // namespace
