#include <gtest/gtest.h>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"
#include "SIREN/dataclasses/VertexWeightingMode.h"
#include "SIREN/distributions/primary/direction/IsotropicDirection.h"
#include "SIREN/distributions/primary/mass/PrimaryMass.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/DetectorDirected3BodyChannel.h"
#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"
#include "SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/injection/Process.h"
#include "SIREN/interactions/CharmMesonDecay.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using siren::dataclasses::InteractionRecord;
using siren::dataclasses::InteractionSignature;
using siren::dataclasses::ParticleType;
using siren::geometry::Geometry;
using siren::geometry::Placement;
using siren::geometry::Sphere;
using siren::injection::DetectorDirected2BodyChannel;
using siren::injection::DetectorDirected3BodyChannel;
using siren::injection::DetectorDirectedAngularSectorChannel;
using siren::injection::DetectorDirectedScatteringChannel;
using siren::injection::Isotropic2BodyChannel;
using siren::injection::MultiChannelPhaseSpace;
using siren::injection::NestedMixtureChannel;
using siren::injection::PhaseSpaceChannel;
using siren::injection::PhaseSpaceTopology;
using siren::injection::PhysicalCrossSectionChannel;
using siren::injection::PhysicalDecayChannel;
using siren::injection::PhysicalProcess;
using siren::injection::TabulatedMappingTable;
using siren::interactions::CharmMesonDecay;
using siren::interactions::DummyCrossSection;
using siren::interactions::InteractionCollection;
using siren::math::Vector3D;

using ChannelPtr = std::shared_ptr<PhaseSpaceChannel>;

template<class OutputArchive, class T>
std::string SaveToString(T const & value) {
    std::stringstream stream;
    {
        OutputArchive archive(stream);
        archive(::cereal::make_nvp("Value", value));
    }
    return stream.str();
}

template<class InputArchive, class T>
T LoadFromString(std::string const & bytes) {
    std::stringstream stream(bytes);
    T loaded;
    InputArchive archive(stream);
    archive(::cereal::make_nvp("Value", loaded));
    return loaded;
}

template<class OutputArchive, class InputArchive, class T>
T RoundTripEqual(T const & value) {
    std::string first = SaveToString<OutputArchive>(value);
    T loaded = LoadFromString<InputArchive, T>(first);
    std::string second = SaveToString<OutputArchive>(loaded);
    EXPECT_EQ(first, second);
    return loaded;
}

std::shared_ptr<Geometry const> MakeTarget() {
    return Sphere(Placement(Vector3D(0.0, 0.0, 0.0)), 10.0, 0.0).create();
}

InteractionSignature DummySignature(ParticleType primary) {
    InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = ParticleType::Nucleon;
    signature.secondary_types = {primary, ParticleType::Nucleon};
    return signature;
}

InteractionSignature CharmSignature(
    std::shared_ptr<CharmMesonDecay> const & decay)
{
    auto signatures = decay->GetPossibleSignaturesFromParent(ParticleType::D0);
    for (auto const & signature : signatures) {
        if (signature.secondary_types.size() == 3) return signature;
    }
    throw std::runtime_error("CharmMesonDecay has no D0 three-body signature");
}

std::vector<std::pair<std::string, ChannelPtr>> MakeConcreteChannels() {
    auto target = MakeTarget();
    auto decay = std::make_shared<CharmMesonDecay>(ParticleType::D0);
    auto cross_section = std::make_shared<DummyCrossSection>();

    auto nested_mixture = std::make_shared<MultiChannelPhaseSpace>(
        std::vector<ChannelPtr>{
            std::make_shared<Isotropic2BodyChannel>(0),
            std::make_shared<Isotropic2BodyChannel>(1)},
        std::vector<double>{0.4, 0.6});
    auto nested = std::make_shared<NestedMixtureChannel>(nested_mixture);
    nested->label = "inner";

    std::vector<std::pair<std::string, ChannelPtr>> channels;
    channels.emplace_back(
        "Isotropic2Body",
        std::make_shared<Isotropic2BodyChannel>(1));
    channels.emplace_back(
        "PhysicalDecay",
        std::make_shared<PhysicalDecayChannel>(
            decay, CharmSignature(decay)));
    channels.emplace_back(
        "PhysicalCrossSection",
        std::make_shared<PhysicalCrossSectionChannel>(
            cross_section, DummySignature(ParticleType::NuE)));
    channels.emplace_back(
        "DetectorDirected2Body",
        std::make_shared<DetectorDirected2BodyChannel>(
            target, 1, DetectorDirected2BodyChannel::Mode::Volume, 100.0));
    channels.emplace_back(
        "DetectorDirectedAngularSector",
        std::make_shared<DetectorDirectedAngularSectorChannel>(
            target, 0.1, 0.8, 0.2, 4.0, 1));
    channels.emplace_back(
        "DetectorDirected3Body",
        std::make_shared<DetectorDirected3BodyChannel>(
            target,
            0,
            DetectorDirected3BodyChannel::InvariantMassMode::Tabulated,
            0.0,
            0.0,
            0.8,
            0.0,
            DetectorDirected2BodyChannel::Mode::Volume,
            PhaseSpaceTopology::Decay3Body,
            std::vector<double>{0.0, 1.0, 2.0},
            std::vector<double>{0.0, 0.25, 1.0},
            100.0));
    channels.emplace_back(
        "DetectorDirectedScattering",
        std::make_shared<DetectorDirectedScatteringChannel>(
            target,
            0,
            DetectorDirectedScatteringChannel::Variable::Q2,
            DetectorDirected2BodyChannel::Mode::Volume,
            DetectorDirectedScatteringChannel::Q2Mode::Tabulated,
            0.0,
            std::vector<double>{0.0, 1.0, 2.0},
            std::vector<double>{0.0, 0.5, 1.0},
            100.0));
    channels.emplace_back("NestedMixture", nested);
    return channels;
}

InteractionRecord ThreeBodyProbeRecord() {
    InteractionRecord record;
    record.signature.primary_type = ParticleType::N4;
    record.signature.target_type = ParticleType::Decay;
    record.signature.secondary_types = {
        ParticleType::NuE, ParticleType::MuMinus, ParticleType::MuPlus};
    record.primary_mass = 5.0;
    record.primary_momentum = {5.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.7, 0.9};
    record.secondary_momenta.resize(3);
    return record;
}

template<class Channel>
void ExpectThreeBodyDensityPreserved(
    std::shared_ptr<Channel> const & channel,
    unsigned int seed)
{
    InteractionRecord record = ThreeBodyProbeRecord();
    auto random = std::make_shared<siren::utilities::SIREN_random>(seed);
    channel->Sample(random, nullptr, record);
    double before = channel->Density(nullptr, record);
    ASSERT_GT(before, 0.0);

    ChannelPtr base = channel;
    ChannelPtr loaded_base = RoundTripEqual<
        cereal::BinaryOutputArchive,
        cereal::BinaryInputArchive>(base);
    auto loaded = std::dynamic_pointer_cast<DetectorDirected3BodyChannel>(
        loaded_base);
    ASSERT_NE(loaded, nullptr);
    EXPECT_DOUBLE_EQ(before, loaded->Density(nullptr, record));
}

} // namespace

TEST(PhaseSpaceSerialization, ConcreteChannelsRoundTripInJsonAndBinary) {
    auto channels = MakeConcreteChannels();
    ASSERT_EQ(channels.size(), 8u);

    for (auto const & entry : channels) {
        SCOPED_TRACE(entry.first);
        ChannelPtr json_loaded = RoundTripEqual<
            cereal::JSONOutputArchive,
            cereal::JSONInputArchive>(entry.second);
        ChannelPtr binary_loaded = RoundTripEqual<
            cereal::BinaryOutputArchive,
            cereal::BinaryInputArchive>(entry.second);
        ASSERT_NE(json_loaded, nullptr);
        ASSERT_NE(binary_loaded, nullptr);
        EXPECT_EQ(json_loaded->Name(), entry.second->Name());
        EXPECT_EQ(binary_loaded->Name(), entry.second->Name());
    }
}

TEST(PhaseSpaceSerialization, MultiChannelStateRoundTrips) {
    auto mixture = std::make_shared<MultiChannelPhaseSpace>(
        std::vector<ChannelPtr>{
            std::make_shared<Isotropic2BodyChannel>(0),
            std::make_shared<Isotropic2BodyChannel>(1)},
        std::vector<double>{0.3, 0.7},
        true);
    mixture->kp_accumulator_ = {1.25, 2.5};
    mixture->kp_count_ = 17;
    mixture->kp_succ_select_ = {3.5, 4.5};
    mixture->kp_fail_select_ = {0.25, 0.75};

    auto json_loaded = RoundTripEqual<
        cereal::JSONOutputArchive,
        cereal::JSONInputArchive>(mixture);
    auto binary_loaded = RoundTripEqual<
        cereal::BinaryOutputArchive,
        cereal::BinaryInputArchive>(mixture);

    auto expect_state = [&mixture](
        std::shared_ptr<MultiChannelPhaseSpace> const & loaded) {
        ASSERT_NE(loaded, nullptr);
        EXPECT_EQ(loaded->channels.size(), mixture->channels.size());
        EXPECT_EQ(loaded->weights, mixture->weights);
        EXPECT_EQ(loaded->allow_incompatible_, mixture->allow_incompatible_);
        EXPECT_EQ(loaded->kp_accumulator_, mixture->kp_accumulator_);
        EXPECT_EQ(loaded->kp_count_, mixture->kp_count_);
        EXPECT_EQ(loaded->kp_succ_select_, mixture->kp_succ_select_);
        EXPECT_EQ(loaded->kp_fail_select_, mixture->kp_fail_select_);
    };
    expect_state(json_loaded);
    expect_state(binary_loaded);
}

TEST(PhaseSpaceSerialization, NestedMixtureRecurses) {
    auto inner = std::make_shared<MultiChannelPhaseSpace>(
        std::vector<ChannelPtr>{
            std::make_shared<Isotropic2BodyChannel>(0),
            std::make_shared<Isotropic2BodyChannel>(1)},
        std::vector<double>{0.2, 0.8});
    inner->kp_accumulator_ = {4.0, 9.0};
    inner->kp_count_ = 6;
    inner->kp_succ_select_ = {1.0, 2.0};
    inner->kp_fail_select_ = {3.0, 4.0};

    auto nested = std::make_shared<NestedMixtureChannel>(inner);
    nested->label = "directed-group";
    auto outer = std::make_shared<MultiChannelPhaseSpace>(
        std::vector<ChannelPtr>{
            nested, std::make_shared<Isotropic2BodyChannel>(0)},
        std::vector<double>{0.65, 0.35});

    auto loaded = RoundTripEqual<
        cereal::BinaryOutputArchive,
        cereal::BinaryInputArchive>(outer);
    ASSERT_NE(loaded, nullptr);
    ASSERT_EQ(loaded->channels.size(), 2u);
    auto loaded_nested = std::dynamic_pointer_cast<NestedMixtureChannel>(
        loaded->channels.at(0));
    ASSERT_NE(loaded_nested, nullptr);
    ASSERT_NE(loaded_nested->mixture, nullptr);
    EXPECT_EQ(loaded_nested->label, "directed-group");
    EXPECT_EQ(loaded_nested->mixture->weights, inner->weights);
    EXPECT_EQ(
        loaded_nested->mixture->kp_accumulator_, inner->kp_accumulator_);
    EXPECT_EQ(loaded_nested->mixture->kp_count_, inner->kp_count_);
    EXPECT_EQ(
        loaded_nested->mixture->kp_succ_select_, inner->kp_succ_select_);
    EXPECT_EQ(
        loaded_nested->mixture->kp_fail_select_, inner->kp_fail_select_);
}

TEST(PhaseSpaceSerialization, PhysicalProcessVersion2PreservesMap) {
    auto cross_section = std::make_shared<DummyCrossSection>();
    auto interactions = std::make_shared<InteractionCollection>(
        ParticleType::NuE,
        std::vector<std::shared_ptr<siren::interactions::CrossSection>>{
            cross_section});
    auto process = std::make_shared<PhysicalProcess>(
        ParticleType::NuE, interactions);

    InteractionSignature first = DummySignature(ParticleType::NuE);
    InteractionSignature second = DummySignature(ParticleType::NuMu);
    auto first_channel = std::make_shared<PhysicalCrossSectionChannel>(
        cross_section, first);
    auto second_channel = std::make_shared<PhysicalCrossSectionChannel>(
        cross_section, second);
    process->SetPhaseSpace(
        first,
        std::make_shared<MultiChannelPhaseSpace>(
            std::vector<ChannelPtr>{first_channel},
            std::vector<double>{1.0}));
    process->SetPhaseSpace(
        second,
        std::make_shared<MultiChannelPhaseSpace>(
            std::vector<ChannelPtr>{first_channel, second_channel},
            std::vector<double>{0.4, 0.6}));

    auto loaded = RoundTripEqual<
        cereal::JSONOutputArchive,
        cereal::JSONInputArchive>(process);
    ASSERT_NE(loaded, nullptr);
    ASSERT_EQ(loaded->GetPhaseSpaceMap().size(), 2u);
    ASSERT_TRUE(loaded->HasPhaseSpace(first));
    ASSERT_TRUE(loaded->HasPhaseSpace(second));
    EXPECT_EQ(loaded->GetPhaseSpace(first)->channels.size(), 1u);
    EXPECT_EQ(loaded->GetPhaseSpace(second)->channels.size(), 2u);
}

TEST(PhaseSpaceSerialization, PhysicalProcessVersion1FixtureLoadsEmptyMap) {
    // Version 1 predates PhaseSpaceMap.
    constexpr char fixture[] = R"JSON({
        "PhysicalProcess": {
            "polymorphic_id": 1073741824,
            "ptr_wrapper": {
                "id": 2147483649,
                "data": {
                    "cereal_class_version": 1,
                    "PhysicalDistributions": [
                        {
                            "polymorphic_id": 2147483649,
                            "polymorphic_name": "siren::distributions::PrimaryMass",
                            "ptr_wrapper": {
                                "id": 2147483650,
                                "data": {
                                    "cereal_class_version": 0,
                                    "PrimaryMass": 0.0,
                                    "value0": {
                                        "cereal_class_version": 0,
                                        "value0": {
                                            "cereal_class_version": 0
                                        }
                                    }
                                }
                            }
                        },
                        {
                            "polymorphic_id": 2147483650,
                            "polymorphic_name": "siren::distributions::IsotropicDirection",
                            "ptr_wrapper": {
                                "id": 2147483651,
                                "data": {
                                    "cereal_class_version": 0,
                                    "value0": {
                                        "cereal_class_version": 0,
                                        "value0": {
                                            "value0": {}
                                        }
                                    }
                                }
                            }
                        }
                    ],
                    "value0": {
                        "cereal_class_version": 0,
                        "PrimaryType": 14,
                        "Interactions": {
                            "polymorphic_id": 1073741824,
                            "ptr_wrapper": {
                                "id": 2147483652,
                                "data": {
                                    "cereal_class_version": 0,
                                    "PrimaryType": 14,
                                    "TargetTypes": [2000002112],
                                    "CrossSections": [{
                                        "polymorphic_id": 2147483651,
                                        "polymorphic_name": "siren::interactions::DummyCrossSection",
                                        "ptr_wrapper": {
                                            "id": 2147483653,
                                            "data": {
                                                "cereal_class_version": 0,
                                                "value0": {"cereal_class_version": 0}
                                            }
                                        }
                                    }],
                                    "Decays": []
                                }
                            }
                        }
                    },
                    "WeightingMode": {
                        "cereal_class_version": 0,
                        "ComputeInteractionProbability": false,
                        "ComputePositionProbability": false,
                        "BoundSource": 2
                    }
                }
            }
        }
    })JSON";

    std::stringstream stream(fixture);
    std::shared_ptr<PhysicalProcess> loaded;
    cereal::JSONInputArchive archive(stream);
    archive(::cereal::make_nvp("PhysicalProcess", loaded));

    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(loaded->GetPrimaryType(), ParticleType::NuMu);
    EXPECT_EQ(loaded->GetPhysicalDistributions().size(), 2u);
    ASSERT_NE(loaded->GetInteractions(), nullptr);
    EXPECT_EQ(loaded->GetInteractions()->GetCrossSections().size(), 1u);
    EXPECT_FALSE(loaded->HasAnyPhaseSpace());
    EXPECT_TRUE(loaded->GetPhaseSpaceMap().empty());
    EXPECT_EQ(
        loaded->GetWeightingMode(),
        siren::dataclasses::VertexWeightingMode::Fixed());
}

TEST(PhaseSpaceSerialization, PhysicalProcessPreservesDecayAliasing) {
    auto decay = std::make_shared<CharmMesonDecay>(ParticleType::D0);
    InteractionSignature signature = CharmSignature(decay);
    auto interactions = std::make_shared<InteractionCollection>(
        ParticleType::D0,
        std::vector<std::shared_ptr<siren::interactions::Decay>>{decay});
    auto process = std::make_shared<PhysicalProcess>(
        ParticleType::D0, interactions);
    auto physical_channel = std::make_shared<PhysicalDecayChannel>(
        decay, signature);
    process->SetPhaseSpace(
        signature,
        std::make_shared<MultiChannelPhaseSpace>(
            std::vector<ChannelPtr>{physical_channel},
            std::vector<double>{1.0}));

    auto loaded = RoundTripEqual<
        cereal::BinaryOutputArchive,
        cereal::BinaryInputArchive>(process);
    ASSERT_NE(loaded, nullptr);
    ASSERT_NE(loaded->GetInteractions(), nullptr);
    ASSERT_EQ(loaded->GetInteractions()->GetDecays().size(), 1u);
    auto loaded_mixture = loaded->GetPhaseSpace(signature);
    ASSERT_NE(loaded_mixture, nullptr);
    ASSERT_EQ(loaded_mixture->channels.size(), 1u);
    auto loaded_channel = std::dynamic_pointer_cast<PhysicalDecayChannel>(
        loaded_mixture->channels.front());
    ASSERT_NE(loaded_channel, nullptr);
    EXPECT_EQ(
        loaded->GetInteractions()->GetDecays().front().get(),
        loaded_channel->GetDecay().get());
}

TEST(PhaseSpaceSerialization, ThreeBodyFactorizationsRestoreDensity) {
    auto target = MakeTarget();
    auto direct = std::make_shared<DetectorDirected3BodyChannel>(
        target,
        0,
        DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0,
        0.0,
        0.8,
        0.0,
        DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology::Decay3Body,
        std::vector<double>{},
        std::vector<double>{},
        100.0);
    auto recursive = std::make_shared<DetectorDirected3BodyChannel>(
        target,
        0,
        1,
        2,
        1,
        DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        0.0,
        0.0,
        0.8,
        0.0,
        DetectorDirected2BodyChannel::Mode::Volume,
        PhaseSpaceTopology::Decay3Body,
        std::vector<double>{},
        std::vector<double>{},
        100.0);

    ExpectThreeBodyDensityPreserved(direct, 1201);
    ExpectThreeBodyDensityPreserved(recursive, 1202);
}

TEST(PhaseSpaceSerialization, TamperedTableIsRejected) {
    constexpr char fixture[] = R"JSON({
        "S": [0.0, 2.0, 1.0],
        "CDF": [0.0, 0.5, 1.0]
    })JSON";
    TabulatedMappingTable table(
        std::vector<double>{0.0, 1.0},
        std::vector<double>{0.0, 1.0});
    std::stringstream stream(fixture);
    cereal::JSONInputArchive archive(stream);
    EXPECT_THROW(table.load(archive, 0), std::runtime_error);
}

TEST(PhaseSpaceSerialization, TamperedMeasureIsRejected) {
    constexpr char fixture[] = R"JSON({
        "Type": 999,
        "Spectator": 0,
        "PairFirst": 1,
        "PairSecond": 2
    })JSON";
    siren::dataclasses::PhaseSpaceMeasure measure;
    std::stringstream stream(fixture);
    cereal::JSONInputArchive archive(stream);
    EXPECT_THROW(measure.load(archive, 0), std::runtime_error);
}
