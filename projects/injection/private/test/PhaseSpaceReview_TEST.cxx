#include <gtest/gtest.h>
#include "../LorentzBoostUtils.h"
#include "SIREN/math/Kinematics.h"

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
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/PhysicalChannelAdapters.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/injection/WeightingUtils.h"
#include "SIREN/injection/Weighter.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Random.h"

#include "../InteractionRecordUtils.h"

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

InteractionRecord ScatteringRecord(
    double beam_energy,
    double target_mass,
    double outgoing_mass,
    double recoil_mass)
{
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 0.0;
    record.primary_momentum = {beam_energy, 0.0, 0.0, beam_energy};
    record.target_mass = target_mass;
    record.interaction_vertex = {0.0, 0.0, -100.0};
    record.secondary_masses = {outgoing_mass, recoil_mass};
    record.secondary_momenta.resize(2);
    return record;
}

InteractionRecord TwoBodyDecayRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 2.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 0.0};
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.5};
    record.secondary_momenta = {
        {11.0, 12.0, 13.0, 14.0},
        {21.0, 22.0, 23.0, 24.0}};
    return record;
}

InteractionRecord BoostedAsymmetricTwoBodyDecayRecord() {
    constexpr double parent_mass = 4.0;
    constexpr double daughter_0_mass = 0.4;
    constexpr double daughter_1_mass = 1.7;
    constexpr double gamma = 1.5;
    double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    constexpr double daughter_1_cos_rest = 0.35;

    double p_rest = siren::injection::TwoBodyRestMomentum(
        parent_mass, daughter_1_mass, daughter_0_mass);
    double sin_rest = std::sqrt(
        1.0 - daughter_1_cos_rest * daughter_1_cos_rest);
    double px_1_rest = p_rest * sin_rest;
    double pz_1_rest = p_rest * daughter_1_cos_rest;
    double E_1_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_1_mass, daughter_0_mass);
    double E_0_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_0_mass, daughter_1_mass);

    auto boost = [beta](double E_rest, double px_rest, double pz_rest) {
        constexpr double gamma = 1.5;
        return std::array<double, 4>{
            gamma * (E_rest + beta * pz_rest),
            px_rest,
            0.0,
            gamma * (pz_rest + beta * E_rest)};
    };

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = parent_mass;
    record.primary_momentum = {
        gamma * parent_mass, 0.0, 0.0, gamma * beta * parent_mass};
    record.secondary_masses = {daughter_0_mass, daughter_1_mass};
    record.secondary_momenta = {
        boost(E_0_rest, -px_1_rest, -pz_1_rest),
        boost(E_1_rest, px_1_rest, pz_1_rest)};
    return record;
}

double DecayLabJacobian(
    InteractionRecord const & record,
    int daughter_index,
    double expected_cos_rest)
{
    int other_index = 1 - daughter_index;
    double parent_p = record.primary_momentum[3];
    double parent_E = record.primary_momentum[0];
    double beta = parent_p / parent_E;
    double gamma = parent_E / record.primary_mass;
    double daughter_mass = record.secondary_masses[daughter_index];
    double other_mass = record.secondary_masses[other_index];
    double p_rest = siren::injection::TwoBodyRestMomentum(
        record.primary_mass, daughter_mass, other_mass);
    double E_rest = siren::injection::TwoBodyRestEnergy(
        record.primary_mass, daughter_mass, other_mass);
    auto const & momentum = record.secondary_momenta[daughter_index];
    double p_lab = std::sqrt(
        momentum[1] * momentum[1] +
        momentum[2] * momentum[2] +
        momentum[3] * momentum[3]);
    double cos_lab = momentum[3] / p_lab;

    auto solutions = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, E_rest, daughter_mass, cos_lab);
    double best_jacobian = 0.0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (auto const & solution : solutions) {
        if (!solution.valid) continue;
        double distance = std::abs(
            solution.cos_theta_rest - expected_cos_rest);
        if (distance < best_distance) {
            best_distance = distance;
            best_jacobian = solution.jacobian;
        }
    }
    return best_jacobian;
}

InteractionRecord ThresholdThreeBodyScatteringRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = 0.1;
    record.primary_momentum = {0.1, 0.0, 0.0, 0.0};
    record.target_mass = 0.1;
    record.interaction_vertex = {0.0, 0.0, 0.0};
    record.secondary_masses = {0.2, 0.1, 0.1};
    record.secondary_momenta = {
        {11.0, 12.0, 13.0, 14.0},
        {21.0, 22.0, 23.0, 24.0},
        {31.0, 32.0, 33.0, 34.0}};
    return record;
}

InteractionRecord AsymmetricThreeBodyDecayRecord() {
    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.secondary_masses = {0.5, 1.0, 1.5};
    record.secondary_momenta = {
        {std::sqrt(0.5 * 0.5 + 0.7 * 0.7), 0.7, 0.0, 0.0},
        {std::sqrt(1.0 * 1.0 + 0.2 * 0.2 + 0.5 * 0.5),
         -0.2, 0.5, 0.0},
        {std::sqrt(1.5 * 1.5 + 0.5 * 0.5 + 0.5 * 0.5),
         -0.5, -0.5, 0.0}};
    record.primary_mass =
        record.secondary_momenta[0][0] +
        record.secondary_momenta[1][0] +
        record.secondary_momenta[2][0];
    record.primary_momentum = {record.primary_mass, 0.0, 0.0, 0.0};
    return record;
}

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

class ConstantChannel final : public PhaseSpaceChannel {
public:
    explicit ConstantChannel(
        double density = 1.0,
        PhaseSpaceTopology topology = PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure measure = PhaseSpaceMeasure::SolidAngleRest())
        : density_(density), topology_(topology), measure_(measure) {}

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override
    {
        return density_;
    }

    std::string Name() const override { return "Constant"; }
    PhaseSpaceTopology Topology() const override { return topology_; }
    PhaseSpaceMeasure Measure() const override { return measure_; }

private:
    double density_;
    PhaseSpaceTopology topology_;
    PhaseSpaceMeasure measure_;
};

class CountingConventionChannel final : public PhaseSpaceChannel {
public:
    mutable int topology_calls = 0;
    mutable int measure_calls = 0;

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random>,
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord &) const override {}

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const>,
        InteractionRecord const &) const override
    {
        return 1.0;
    }

    std::string Name() const override { return "CountingConvention"; }
    PhaseSpaceTopology Topology() const override {
        ++topology_calls;
        return PhaseSpaceTopology::Decay2Body;
    }
    PhaseSpaceMeasure Measure() const override {
        ++measure_calls;
        return PhaseSpaceMeasure::SolidAngleRest();
    }
};

siren::dataclasses::InteractionSignature SignatureWithSecondaries(size_t count) {
    siren::dataclasses::InteractionSignature signature;
    signature.secondary_types.assign(
        count, siren::dataclasses::ParticleType::unknown);
    return signature;
}

class MixedSignatureDecay final : public siren::interactions::Decay {
public:
    MixedSignatureDecay()
        : signatures_{SignatureWithSecondaries(2), SignatureWithSecondaries(3)} {}

    bool equal(siren::interactions::Decay const & other) const override {
        return dynamic_cast<MixedSignatureDecay const *>(&other) != nullptr;
    }
    double TotalDecayWidthAllFinalStates(InteractionRecord const &) const override {
        return 1.0;
    }
    double TotalDecayWidth(siren::dataclasses::ParticleType) const override {
        return 1.0;
    }
    double TotalDecayWidth(InteractionRecord const &) const override {
        return 1.0;
    }
    double DifferentialDecayWidth(InteractionRecord const &) const override {
        return 1.0;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override {}
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return signatures_;
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParent(
        siren::dataclasses::ParticleType) const override {
        return signatures_;
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return 1.0;
    }
    std::vector<std::string> DensityVariables() const override {
        return {"cos_theta"};
    }
    PhaseSpaceMeasure MeasureForSignature(
        siren::dataclasses::InteractionSignature const & signature) const override {
        return signature.secondary_types.size() == 2
            ? PhaseSpaceMeasure::SolidAngleRest()
            : PhaseSpaceMeasure::HelicityAngles();
    }

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
};

class MixedSignatureCrossSection final
    : public siren::interactions::CrossSection {
public:
    explicit MixedSignatureCrossSection(
        PhaseSpaceMeasure measure = PhaseSpaceMeasure::MandelstamQ2())
        : signatures_{SignatureWithSecondaries(2), SignatureWithSecondaries(3)}
        , measure_(measure) {}

    bool equal(siren::interactions::CrossSection const & other) const override {
        return dynamic_cast<MixedSignatureCrossSection const *>(&other) != nullptr;
    }
    double TotalCrossSection(InteractionRecord const &) const override {
        return 1.0;
    }
    double DifferentialCrossSection(InteractionRecord const &) const override {
        return 1.0;
    }
    double InteractionThreshold(InteractionRecord const &) const override {
        return 0.0;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override {}
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        return {siren::dataclasses::ParticleType::unknown};
    }
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(
        siren::dataclasses::ParticleType) const override {
        return GetPossibleTargets();
    }
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        return {siren::dataclasses::ParticleType::unknown};
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignatures() const override {
        return signatures_;
    }
    std::vector<siren::dataclasses::InteractionSignature>
    GetPossibleSignaturesFromParents(
        siren::dataclasses::ParticleType,
        siren::dataclasses::ParticleType) const override {
        return signatures_;
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return 1.0;
    }
    std::vector<std::string> DensityVariables() const override {
        return {};
    }
    PhaseSpaceMeasure MeasureForSignature(
        siren::dataclasses::InteractionSignature const &) const override {
        return measure_;
    }

private:
    std::vector<siren::dataclasses::InteractionSignature> signatures_;
    PhaseSpaceMeasure measure_;
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

MultiChannelPhaseSpace TwoChannelMixture(std::vector<double> weights) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(),
        std::make_shared<ConstantChannel>()};
    mixture.weights = std::move(weights);
    return mixture;
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
        auto phase_space = std::make_shared<MultiChannelPhaseSpace>();
        phase_space->channels = {std::make_shared<ConstantChannel>()};
        phase_space->weights = {1.0};
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
        auto phase_space = std::make_shared<MultiChannelPhaseSpace>();
        phase_space->channels = {std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())};
        phase_space->weights = {1.0};
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

        double density = siren::injection::SelectedFinalStateProbability(
            detector, interactions, record,
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

TEST(MultiChannelWeights, RejectsSizeMismatch) {
    auto mixture = TwoChannelMixture({1.0});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_THROW(mixture.Sample(random, nullptr, record),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(mixture.Density(nullptr, record),
                 siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, RejectsUnnormalizedWeights) {
    auto mixture = TwoChannelMixture({0.2, 0.2});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_THROW(mixture.Sample(random, nullptr, record),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(mixture.Density(nullptr, record),
                 siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, RejectsNegativeAndNonFiniteWeights) {
    InteractionRecord record;

    EXPECT_THROW(TwoChannelMixture({-0.1, 1.1}).Density(nullptr, record),
                 siren::utilities::ConfigurationError);
    EXPECT_THROW(TwoChannelMixture({
        std::numeric_limits<double>::quiet_NaN(), 1.0}).Density(nullptr, record),
        siren::utilities::ConfigurationError);
    EXPECT_THROW(TwoChannelMixture({
        std::numeric_limits<double>::infinity(), 0.0}).Density(nullptr, record),
        siren::utilities::ConfigurationError);
}

TEST(MultiChannelWeights, AcceptsNormalizedWeights) {
    auto mixture = TwoChannelMixture({0.25, 0.75});
    auto random = std::make_shared<siren::utilities::SIREN_random>(1);
    InteractionRecord record;

    EXPECT_NO_THROW(mixture.Sample(random, nullptr, record));
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
}

TEST(MultiChannelConventions, CompatibleMixtureCachesSuccessfulValidation) {
    auto first = std::make_shared<CountingConventionChannel>();
    auto second = std::make_shared<CountingConventionChannel>();
    MultiChannelPhaseSpace mixture;
    mixture.channels = {first, second};
    mixture.weights = {0.5, 0.5};
    InteractionRecord record;

    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    int first_topology_calls = first->topology_calls;
    int first_measure_calls = first->measure_calls;
    int second_topology_calls = second->topology_calls;
    int second_measure_calls = second->measure_calls;

    for (int i = 0; i < 100; ++i) {
        EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    }
    EXPECT_EQ(first->topology_calls, first_topology_calls + 100);
    EXPECT_EQ(first->measure_calls, first_measure_calls + 100);
    EXPECT_EQ(second->topology_calls, second_topology_calls + 100);
    EXPECT_EQ(second->measure_calls, second_measure_calls + 100);

    // Public channel replacement changes the fingerprint and rebuilds the
    // cached conventions on the next evaluation: the fingerprint probe plus
    // the rebuild probe, and nothing more.
    auto replacement = std::make_shared<CountingConventionChannel>();
    mixture.channels[1] = replacement;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 1.0);
    EXPECT_EQ(first->topology_calls, first_topology_calls + 102);
    EXPECT_EQ(first->measure_calls, first_measure_calls + 102);
    EXPECT_EQ(replacement->topology_calls, 2);
    EXPECT_EQ(replacement->measure_calls, 2);
}

TEST(MultiChannelConventions, ClearedAndReallocatedChannelRebuildsCache) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {std::make_shared<ConstantChannel>(
        2.0, PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2())};
    mixture.weights = {1.0};
    InteractionRecord record;

    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 2.0);
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::MandelstamQ2());

    mixture.channels.clear();
    mixture.channels.push_back(std::make_shared<ConstantChannel>(
        4.0, PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::FixedMassY()));

    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::FixedMassY());
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 4.0);
}

TEST(MultiChannelConventions, ScatteringLabAngleMixtureIsFatalNotConvertible) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleLab())};
    mixture.weights = {0.5, 0.5};

    auto diagnostics = mixture.ValidateChannelsDetailed();
    ASSERT_EQ(diagnostics.size(), 1u);
    EXPECT_EQ(diagnostics[0].severity,
              MultiChannelPhaseSpace::ChannelDiagnostic::Severity::Fatal);
    EXPECT_NE(diagnostics[0].message.find("not convertible"),
              std::string::npos);

    InteractionRecord record = ScatteringRecord(10.0, 1.0, 0.1, 1.0);
    auto random = std::make_shared<siren::utilities::SIREN_random>(8675309);
    EXPECT_THROW(mixture.Density(nullptr, record),
                 siren::utilities::MeasureCompatibilityError);
    EXPECT_THROW(mixture.Sample(random, nullptr, record),
                 siren::utilities::MeasureCompatibilityError);
}

TEST(PhysicalAdapterSignature, PinsDecayTopologyAndMeasure) {
    auto decay = std::make_shared<MixedSignatureDecay>();
    auto two_body = SignatureWithSecondaries(2);
    auto three_body = SignatureWithSecondaries(3);

    siren::injection::PhysicalDecayChannel unpinned(decay);
    EXPECT_EQ(unpinned.Topology(), PhaseSpaceTopology::Unspecified);
    EXPECT_EQ(unpinned.Measure(), PhaseSpaceMeasure::Unspecified());

    siren::injection::PhysicalDecayChannel pinned_two(decay, two_body);
    EXPECT_EQ(pinned_two.Topology(), PhaseSpaceTopology::Decay2Body);
    EXPECT_EQ(pinned_two.Measure(), PhaseSpaceMeasure::SolidAngleRest());

    siren::injection::PhysicalDecayChannel pinned_three(decay, three_body);
    EXPECT_EQ(pinned_three.Topology(), PhaseSpaceTopology::Decay3Body);
    EXPECT_EQ(pinned_three.Measure(), PhaseSpaceMeasure::HelicityAngles());
}

TEST(PhysicalAdapterSignature, PinsCrossSectionTopologyAndMeasure) {
    auto cross_section = std::make_shared<MixedSignatureCrossSection>();
    auto two_body = SignatureWithSecondaries(2);
    auto three_body = SignatureWithSecondaries(3);

    siren::injection::PhysicalCrossSectionChannel unpinned(cross_section);
    EXPECT_EQ(unpinned.Topology(), PhaseSpaceTopology::Unspecified);
    EXPECT_EQ(unpinned.Measure(), PhaseSpaceMeasure::Unspecified());

    siren::injection::PhysicalCrossSectionChannel pinned_two(
        cross_section, two_body);
    EXPECT_EQ(pinned_two.Topology(), PhaseSpaceTopology::Scatter2to2);
    EXPECT_EQ(pinned_two.Measure(), PhaseSpaceMeasure::MandelstamQ2());

    siren::injection::PhysicalCrossSectionChannel pinned_three(
        cross_section, three_body);
    EXPECT_EQ(pinned_three.Topology(), PhaseSpaceTopology::Scatter2to3);
    EXPECT_EQ(pinned_three.Measure(), PhaseSpaceMeasure::MandelstamQ2());
}

TEST(InteractionMeasureDeclaration, UndeclaredMeasuresAreUnspecified) {
    // Density-variable names such as "Bjorken x" do not declare a measure;
    // only an override does.
    auto signature = SignatureWithSecondaries(2);
    siren::interactions::DummyCrossSection undeclared;
    EXPECT_EQ(undeclared.Measure(), PhaseSpaceMeasure::Unspecified());
    EXPECT_EQ(undeclared.MeasureForSignature(signature),
              PhaseSpaceMeasure::Unspecified());

    MixedSignatureCrossSection declared(PhaseSpaceMeasure::FixedMassYPhi());
    EXPECT_EQ(declared.MeasureForSignature(signature),
              PhaseSpaceMeasure::FixedMassYPhi());
}

TEST(AzimuthTaxonomy, PredicatesAndCompletionsAgree) {
    using siren::dataclasses::MeasureHasExplicitAzimuth;
    using siren::dataclasses::MeasureIntegratesAzimuth;
    using siren::dataclasses::MeasureWithExplicitAzimuth;

    std::pair<PhaseSpaceMeasure, PhaseSpaceMeasure> lifts[] = {
        {PhaseSpaceMeasure::MandelstamQ2(), PhaseSpaceMeasure::MandelstamQ2Phi()},
        {PhaseSpaceMeasure::FixedMassY(), PhaseSpaceMeasure::FixedMassYPhi()},
        {PhaseSpaceMeasure::BjorkenXY(), PhaseSpaceMeasure::BjorkenXYPhi()},
        {PhaseSpaceMeasure::MandelstamQ2Y(), PhaseSpaceMeasure::MandelstamQ2YPhi()},
    };
    for (auto const & [marginal, joint] : lifts) {
        EXPECT_TRUE(MeasureIntegratesAzimuth(marginal));
        EXPECT_FALSE(MeasureHasExplicitAzimuth(marginal));
        EXPECT_TRUE(MeasureHasExplicitAzimuth(joint));
        EXPECT_FALSE(MeasureIntegratesAzimuth(joint));
        EXPECT_EQ(MeasureWithExplicitAzimuth(marginal), joint);
        EXPECT_EQ(MeasureWithExplicitAzimuth(joint), joint);
    }
    EXPECT_TRUE(MeasureHasExplicitAzimuth(PhaseSpaceMeasure::SolidAngleRest()));
    EXPECT_FALSE(MeasureIntegratesAzimuth(PhaseSpaceMeasure::SolidAngleRest()));
    EXPECT_EQ(MeasureWithExplicitAzimuth(PhaseSpaceMeasure::SolidAngleRest()),
              PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_FALSE(MeasureHasExplicitAzimuth(PhaseSpaceMeasure::Unspecified()));
    EXPECT_FALSE(MeasureIntegratesAzimuth(PhaseSpaceMeasure::Unspecified()));
}

TEST(WeightingConvention, ElectsOneCommonConventionForBothDirections) {
    siren::injection::PhaseSpaceConvention marginal{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2()};
    siren::injection::PhaseSpaceConvention joint{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2Phi()};

    EXPECT_EQ(
        siren::injection::ResolveCommonFinalStateConvention(marginal, joint),
        joint);
    EXPECT_EQ(
        siren::injection::ResolveCommonFinalStateConvention(joint, marginal),
        joint);

    siren::injection::PhaseSpaceConvention different_family{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::BjorkenXYPhi()};
    EXPECT_THROW(
        siren::injection::ResolveCommonFinalStateConvention(
            marginal, different_family),
        siren::utilities::MeasureCompatibilityError);

    siren::injection::PhaseSpaceConvention different_topology{
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure::SolidAngleRest()};
    EXPECT_THROW(
        siren::injection::ResolveCommonFinalStateConvention(
            joint, different_topology),
        siren::utilities::MeasureCompatibilityError);
}

TEST(PhysicalChannelAdapters, AcceptExplicitSignatureConventionOverride) {
    auto scatter_signature = SignatureWithSecondaries(2);
    auto cross_section = std::make_shared<MixedSignatureCrossSection>();
    siren::injection::PhaseSpaceConvention scatter_convention{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::MandelstamQ2Phi()};

    siren::injection::PhysicalCrossSectionChannel scatter(
        cross_section, scatter_signature, scatter_convention);
    EXPECT_EQ(scatter.Topology(), scatter_convention.topology);
    EXPECT_EQ(scatter.Measure(), scatter_convention.measure);

    auto decay_signature = SignatureWithSecondaries(2);
    auto decay = std::make_shared<MixedSignatureDecay>();
    siren::injection::PhaseSpaceConvention decay_convention{
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure::SolidAngleLab(1)};

    siren::injection::PhysicalDecayChannel decay_channel(
        decay, decay_signature, decay_convention);
    EXPECT_EQ(decay_channel.Topology(), decay_convention.topology);
    EXPECT_EQ(decay_channel.Measure(), decay_convention.measure);

    siren::injection::PhaseSpaceConvention wrong_topology{
        PhaseSpaceTopology::Decay2Body,
        PhaseSpaceMeasure::SolidAngleRest()};
    EXPECT_THROW(
        siren::injection::PhysicalCrossSectionChannel(
            cross_section, scatter_signature, wrong_topology),
        siren::utilities::ConfigurationError);

    siren::injection::PhaseSpaceConvention unspecified_measure{
        PhaseSpaceTopology::Scatter2to2,
        PhaseSpaceMeasure::Unspecified()};
    EXPECT_THROW(
        siren::injection::PhysicalCrossSectionChannel(
            cross_section, scatter_signature, unspecified_measure),
        siren::utilities::ConfigurationError);
}

TEST(WeightingConvention, LiftsNaturalFixedMassYIntoJointProposalMeasure) {
    auto signature = SignatureWithSecondaries(2);
    auto cross_section = std::make_shared<MixedSignatureCrossSection>(
        PhaseSpaceMeasure::FixedMassY());
    auto interactions =
        std::make_shared<siren::interactions::InteractionCollection>(
            siren::dataclasses::ParticleType::unknown,
            std::vector<std::shared_ptr<siren::interactions::CrossSection>>{
                cross_section});

    InteractionRecord record;
    record.signature = signature;
    siren::injection::PhaseSpaceConvention joint_convention;
    joint_convention.topology = PhaseSpaceTopology::Scatter2to2;
    joint_convention.measure = PhaseSpaceMeasure::FixedMassYPhi();

    EXPECT_NEAR(
        siren::injection::SelectedFinalStateProbability(
            nullptr, interactions, record, joint_convention),
        1.0 / (2.0 * M_PI), 1e-14);
}

std::shared_ptr<siren::interactions::InteractionCollection> CollectionOf(
    std::shared_ptr<siren::interactions::CrossSection> cross_section) {
    return std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::unknown,
        std::vector<std::shared_ptr<siren::interactions::CrossSection>>{cross_section});
}

std::shared_ptr<MultiChannelPhaseSpace> ConstantMixture(
    PhaseSpaceTopology topology, PhaseSpaceMeasure measure) {
    auto mixture = std::make_shared<MultiChannelPhaseSpace>();
    mixture->channels = {std::make_shared<ConstantChannel>(1.0, topology, measure)};
    mixture->weights = {1.0};
    return mixture;
}

TEST(ProcessPhaseSpaceValidation, SetPhaseSpaceChecksOnlyTheMixture) {
    // Registration cannot know which process the proposal will be weighted
    // against, so it accepts any internally consistent mixture.
    auto signature = SignatureWithSecondaries(2);
    siren::injection::PhysicalProcess process(
        siren::dataclasses::ParticleType::unknown,
        CollectionOf(std::make_shared<MixedSignatureCrossSection>(
            PhaseSpaceMeasure::FixedMassYPhi())));
    EXPECT_NO_THROW(process.SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::FixedMassY())));
    EXPECT_NO_THROW(process.SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::SolidAngleRest())));

    auto inconsistent = std::make_shared<MultiChannelPhaseSpace>();
    inconsistent->channels = {
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2Phi()),
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::SolidAngleRest())};
    inconsistent->weights = {0.5, 0.5};
    EXPECT_THROW(process.SetPhaseSpace(signature, inconsistent),
                 siren::utilities::MeasureCompatibilityError);
}

TEST(ProcessPhaseSpaceValidation, WeighterRejectsProposalWithoutCommonMeasure) {
    // A (Q2, phi) proposal cannot be compared with a Bjorken (x, y) density.
    auto signature = SignatureWithSecondaries(2);
    auto interactions = CollectionOf(
        std::make_shared<MixedSignatureCrossSection>(PhaseSpaceMeasure::BjorkenXY()));
    auto injection = std::make_shared<siren::injection::PrimaryInjectionProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    auto physical = std::make_shared<siren::injection::PhysicalProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    injection->SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Scatter2to2, PhaseSpaceMeasure::MandelstamQ2Phi()));
    EXPECT_THROW(siren::injection::PrimaryProcessWeighter(physical, injection, nullptr),
                 siren::utilities::MeasureCompatibilityError);

    // Registered on both processes, the proposal replaces the model density on
    // both sides, so its chart need not match the model's.
    physical->SetPhaseSpace(signature, injection->GetPhaseSpace(signature));
    EXPECT_NO_THROW(siren::injection::PrimaryProcessWeighter(physical, injection, nullptr));
}

TEST(ProcessPhaseSpaceValidation, WeighterLiftsPerCosThetaProposal) {
    // A per-cos(theta) proposal is compared with a solid-angle model density
    // by lifting the proposal's declared uniform azimuth.
    auto signature = SignatureWithSecondaries(2);
    auto interactions = std::make_shared<siren::interactions::InteractionCollection>(
        siren::dataclasses::ParticleType::unknown,
        std::vector<std::shared_ptr<siren::interactions::Decay>>{
            std::make_shared<MixedSignatureDecay>()});
    auto injection = std::make_shared<siren::injection::PrimaryInjectionProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    auto physical = std::make_shared<siren::injection::PhysicalProcess>(
        siren::dataclasses::ParticleType::unknown, interactions);
    injection->SetPhaseSpace(signature, ConstantMixture(
        PhaseSpaceTopology::Decay2Body, PhaseSpaceMeasure::CosThetaRest()));
    EXPECT_NO_THROW(siren::injection::PrimaryProcessWeighter(physical, injection, nullptr));
}

TEST(CommonMeasure, UnspecifiedMajorityCannotOutvoteSpecifiedChannel) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            2.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            8.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest())};
    mixture.weights = {0.25, 0.25, 0.5};

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_THROW(mixture.Density(nullptr, record), std::runtime_error);
    EXPECT_THROW(mixture.DensityBreakdown(nullptr, record), std::runtime_error);
}

TEST(CommonMeasure, ZeroUnspecifiedDensityCannotBypassMixedMeasureRejection) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            8.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            0.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_THROW(mixture.Density(nullptr, record), std::runtime_error);
}

TEST(CommonMeasure, AllUnspecifiedChannelsRemainUnspecified) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            2.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified()),
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::Unspecified())};
    mixture.weights = {0.25, 0.75};

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(), PhaseSpaceMeasure::Unspecified());
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 3.5);
}

TEST(DecayMeasureConversion, UsesLabMeasureDaughterIndexInBothDirections) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    PhaseSpaceMeasure lab_0 = PhaseSpaceMeasure::SolidAngleLab(0);
    PhaseSpaceMeasure lab_1 = PhaseSpaceMeasure::SolidAngleLab(1);
    EXPECT_NE(lab_0, lab_1);

    double jacobian_0 = DecayLabJacobian(record, 0, -0.35);
    double jacobian_1 = DecayLabJacobian(record, 1, 0.35);
    ASSERT_GT(jacobian_0, 0.0);
    ASSERT_GT(jacobian_1, 0.0);
    ASSERT_GT(std::abs(jacobian_1 - jacobian_0), 1e-3);

    MultiChannelPhaseSpace lab_to_rest;
    lab_to_rest.channels = {
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_1)};
    lab_to_rest.weights = {0.5, 0.5};
    EXPECT_EQ(
        lab_to_rest.CommonMeasure(), PhaseSpaceMeasure::SolidAngleRest());
    EXPECT_NEAR(
        lab_to_rest.Density(nullptr, record),
        0.5 * 3.0 + 0.5 * 5.0 * jacobian_1, 1e-13);

    MultiChannelPhaseSpace rest_to_lab;
    rest_to_lab.channels = {
        std::make_shared<ConstantChannel>(
            4.0, PhaseSpaceTopology::Decay2Body, lab_1),
        std::make_shared<ConstantChannel>(
            6.0, PhaseSpaceTopology::Decay2Body, lab_1),
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest())};
    rest_to_lab.weights = {0.25, 0.25, 0.5};
    EXPECT_EQ(rest_to_lab.CommonMeasure(), lab_1);
    EXPECT_NEAR(
        rest_to_lab.Density(nullptr, record),
        0.25 * 4.0 + 0.25 * 6.0 + 0.5 * 3.0 / jacobian_1,
        1e-13);

    MultiChannelPhaseSpace lab_daughter_conversion;
    lab_daughter_conversion.channels = {
        std::make_shared<ConstantChannel>(
            7.0, PhaseSpaceTopology::Decay2Body, lab_0),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_1)};
    lab_daughter_conversion.weights = {0.5, 0.5};
    EXPECT_EQ(lab_daughter_conversion.CommonMeasure(), lab_0);
    EXPECT_NEAR(
        lab_daughter_conversion.Density(nullptr, record),
        0.5 * 7.0 + 0.5 * 5.0 * jacobian_1 / jacobian_0,
        1e-13);
}

TEST(DecayMeasureConversion, RejectsMissingMomentaAndInvalidDaughterIndex) {
    auto make_mixture = [](PhaseSpaceMeasure lab_measure) {
        MultiChannelPhaseSpace mixture;
        mixture.channels = {
            std::make_shared<ConstantChannel>(
                3.0, PhaseSpaceTopology::Decay2Body,
                PhaseSpaceMeasure::SolidAngleRest()),
            std::make_shared<ConstantChannel>(
                5.0, PhaseSpaceTopology::Decay2Body, lab_measure)};
        mixture.weights = {0.5, 0.5};
        return mixture;
    };

    InteractionRecord missing_momentum =
        BoostedAsymmetricTwoBodyDecayRecord();
    missing_momentum.secondary_momenta.resize(1);
    EXPECT_THROW(
        make_mixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, missing_momentum),
        std::runtime_error);

    InteractionRecord missing_mass = BoostedAsymmetricTwoBodyDecayRecord();
    missing_mass.secondary_masses.resize(1);
    EXPECT_THROW(
        make_mixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, missing_mass),
        std::runtime_error);

    InteractionRecord complete = BoostedAsymmetricTwoBodyDecayRecord();
    EXPECT_THROW(
        make_mixture(PhaseSpaceMeasure::SolidAngleLab(2)).Density(
            nullptr, complete),
        std::runtime_error);
}

namespace {

MultiChannelPhaseSpace RestPlusLabMixture(PhaseSpaceMeasure lab_measure) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Decay2Body,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Decay2Body, lab_measure)};
    mixture.weights = {0.5, 0.5};
    return mixture;
}

} // anonymous namespace

TEST(DecayMeasureConversion, ThrowsWhenLabAngleOutsideAllowedCone) {
    // Heavy daughters on a fast parent are confined to a narrow forward
    // cone; a backward daughter direction is kinematically impossible
    // for these masses, and the conversion must fail loudly.
    constexpr double parent_mass = 1.0;
    constexpr double daughter_mass = 0.45;
    constexpr double beta = 0.9;
    const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    constexpr double p_backward = 0.3;
    const double E_backward = std::sqrt(
        p_backward * p_backward + daughter_mass * daughter_mass);

    InteractionRecord record;
    record.signature.secondary_types = {
        siren::dataclasses::ParticleType::unknown,
        siren::dataclasses::ParticleType::unknown};
    record.primary_mass = parent_mass;
    record.primary_momentum = {
        gamma * parent_mass, 0.0, 0.0, gamma * beta * parent_mass};
    record.secondary_masses = {daughter_mass, daughter_mass};
    record.secondary_momenta = {
        {E_backward, 0.0, 0.0, -p_backward},
        {0.0, 0.0, 0.0, 0.0}};

    EXPECT_THROW(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(0)).Density(
            nullptr, record),
        std::runtime_error);
}

TEST(DecayMeasureConversion, ThrowsOnDegenerateDaughterMomentum) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    record.secondary_momenta[1] = {record.secondary_masses[1], 0.0, 0.0, 0.0};
    EXPECT_THROW(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, record),
        std::runtime_error);
}

TEST(DecayMeasureConversion, ThrowsOnSubThresholdRecordMasses) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    record.primary_mass = 0.99 * (record.secondary_masses[0] +
                                  record.secondary_masses[1]);
    EXPECT_THROW(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, record),
        std::runtime_error);
}

TEST(DecayMeasureConversion, ParentAtRestConvertsAsIdentity) {
    InteractionRecord record = BoostedAsymmetricTwoBodyDecayRecord();
    record.primary_momentum = {record.primary_mass, 0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(
        RestPlusLabMixture(PhaseSpaceMeasure::SolidAngleLab(1)).Density(
            nullptr, record),
        0.5 * 3.0 + 0.5 * 5.0);
}

TEST(ThreeBodyMeasureConversion, DalitzCrossFactorizationHasUnitJacobian) {
    namespace J = siren::injection::phase_space_jacobian;
    InteractionRecord record = AsymmetricThreeBodyDecayRecord();
    PhaseSpaceMeasure common = PhaseSpaceMeasure::DalitzPair(0, 1, 2);
    PhaseSpaceMeasure alternate = PhaseSpaceMeasure::DalitzPair(1, 0, 2);
    ASSERT_NE(common, alternate);

    auto pair_mass_squared = [&record](int first, int second) {
        auto const & p1 = record.secondary_momenta[first];
        auto const & p2 = record.secondary_momenta[second];
        double E = p1[0] + p2[0];
        double px = p1[1] + p2[1];
        double py = p1[2] + p2[2];
        double pz = p1[3] + p2[3];
        return E * E - px * px - py * py - pz * pz;
    };

    constexpr double alternate_density = 7.0;
    double old_intermediate = J::Recursive2BodyDensityToDalitzDensity(
        alternate_density, record.primary_mass,
        record.secondary_masses[alternate.spectator],
        record.secondary_masses[alternate.pair_first],
        record.secondary_masses[alternate.pair_second],
        pair_mass_squared(alternate.pair_first, alternate.pair_second));
    double old_converted = J::DalitzDensityToRecursive2BodyDensity(
        old_intermediate, record.primary_mass,
        record.secondary_masses[common.spectator],
        record.secondary_masses[common.pair_first],
        record.secondary_masses[common.pair_second],
        pair_mass_squared(common.pair_first, common.pair_second));
    ASSERT_GT(std::abs(old_converted - alternate_density), 1e-3);

    for (PhaseSpaceTopology topology : {
             PhaseSpaceTopology::Decay3Body,
             PhaseSpaceTopology::Scatter2to3}) {
        SCOPED_TRACE(siren::dataclasses::PhaseSpaceTopologyName(topology));
        MultiChannelPhaseSpace mixture;
        mixture.channels = {
            std::make_shared<ConstantChannel>(2.0, topology, common),
            std::make_shared<ConstantChannel>(
                alternate_density, topology, alternate)};
        mixture.weights = {0.25, 0.75};

        EXPECT_EQ(mixture.CommonMeasure(), common);
        EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), 5.75);
        auto contributions = mixture.DensityBreakdown(nullptr, record);
        ASSERT_EQ(contributions.size(), 2u);
        EXPECT_DOUBLE_EQ(contributions[0], 0.5);
        EXPECT_DOUBLE_EQ(contributions[1], 5.25);
    }
}

TEST(ScatteringMeasureConversion, UsesIncomingAndOutgoingCmMomenta) {
    namespace J = siren::injection::phase_space_jacobian;
    constexpr double s = 25.0;
    constexpr double m_beam = 0.5;
    constexpr double m_target = 1.0;
    constexpr double m_outgoing = 2.0;
    constexpr double m_recoil = 0.75;

    double p_in_sq = siren::injection::Kallen(
        s, m_beam * m_beam, m_target * m_target) / (4.0 * s);
    double p_out_sq = siren::injection::Kallen(
        s, m_outgoing * m_outgoing, m_recoil * m_recoil) / (4.0 * s);
    double expected = 2.0 * std::sqrt(p_in_sq * p_out_sq);

    EXPECT_DOUBLE_EQ(J::SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target, m_outgoing, m_recoil), expected);
    EXPECT_NE(J::SolidAngleRestToMandelstamQ2AbsJacobian(
        s, m_beam, m_target), expected);
}

TEST(ScatteringMeasureConversion, MultiChannelUsesInelasticJacobian) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            3.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            5.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    record.primary_mass = 0.5;
    record.target_mass = 1.0;
    record.primary_momentum = {12.0, 0.0, 0.0, 0.0};
    record.secondary_masses = {2.0, 0.75};
    double s = record.primary_mass * record.primary_mass
             + record.target_mass * record.target_mass
             + 2.0 * record.target_mass * record.primary_momentum[0];
    double jacobian = siren::injection::phase_space_jacobian::
        SolidAngleRestToMandelstamQ2AbsJacobian(
            s, record.primary_mass, record.target_mass,
            record.secondary_masses[0], record.secondary_masses[1]);

    double expected = 0.5 * 3.0 + 0.5 * 5.0 * jacobian / (2.0 * M_PI);
    EXPECT_NEAR(mixture.Density(nullptr, record), expected, 1e-14);
}


TEST(ScatteringMeasureConversion, MixedFixedMassYAndQ2HasNoInverseYInflation) {
    constexpr double target_mass = 0.020;
    constexpr double incident_energy = 0.300;
    constexpr double y_density = 3.0;
    constexpr double q2_density = 5.0;
    double jacobian = 2.0 * target_mass * incident_energy;

    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            y_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::FixedMassY()),
        std::make_shared<ConstantChannel>(
            q2_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    record.target_mass = target_mass;
    record.primary_momentum = {incident_energy, 0.0, 0.0, incident_energy};
    double expected = 0.5 * y_density / jacobian + 0.5 * q2_density;

    record.interaction_parameters["bjorken_y"] = 0.2;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), expected);
    record.interaction_parameters["bjorken_y"] = 0.8;
    EXPECT_DOUBLE_EQ(mixture.Density(nullptr, record), expected);
}

TEST(ScatteringMeasureConversion, ExplicitAzimuthWinsAndLiftsMarginal) {
    constexpr double marginal_density = 6.0;
    constexpr double joint_density = 2.0;

    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            marginal_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            marginal_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2()),
        std::make_shared<ConstantChannel>(
            joint_density, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::MandelstamQ2Phi())};
    mixture.weights = {0.25, 0.25, 0.5};

    InteractionRecord record;
    EXPECT_EQ(mixture.CommonMeasure(),
              PhaseSpaceMeasure::MandelstamQ2Phi());
    EXPECT_NEAR(
        mixture.Density(nullptr, record),
        0.5 * marginal_density / (2.0 * M_PI) + 0.5 * joint_density,
        1e-14);

    siren::injection::PhaseSpaceConvention common = mixture.CommonConvention();
    EXPECT_EQ(common.topology, PhaseSpaceTopology::Scatter2to2);
    EXPECT_EQ(common.measure, PhaseSpaceMeasure::MandelstamQ2Phi());
    EXPECT_DOUBLE_EQ(mixture.DensityIn(nullptr, record, common),
                     mixture.Density(nullptr, record));

    siren::injection::PhaseSpaceConvention wrong_topology = common;
    wrong_topology.topology = PhaseSpaceTopology::Decay2Body;
    EXPECT_THROW(
        mixture.DensityIn(nullptr, record, wrong_topology),
        siren::utilities::MeasureCompatibilityError);
}

TEST(ScatteringMeasureConversion, ExplicitAzimuthCannotBeMarginalizedPointwise) {
    InteractionRecord record;
    EXPECT_THROW(
        siren::injection::ConvertDensity(
            1.0,
            PhaseSpaceMeasure::FixedMassYPhi(),
            PhaseSpaceMeasure::FixedMassY(),
            PhaseSpaceTopology::Scatter2to2,
            record),
        siren::utilities::MeasureCompatibilityError);
    EXPECT_THROW(
        siren::injection::ConvertDensity(
            0.0,
            PhaseSpaceMeasure::FixedMassYPhi(),
            PhaseSpaceMeasure::FixedMassY(),
            PhaseSpaceTopology::Scatter2to2,
            record),
        siren::utilities::MeasureCompatibilityError);
}

TEST(ScatteringMeasureConversion, BjorkenXYConvertsOnlyToQ2Y) {
    InteractionRecord record;
    record.target_mass = 0.938;
    record.primary_momentum = {5.0, 0.0, 0.0, 5.0};
    record.interaction_parameters["bjorken_y"] = 0.4;

    double converted = siren::injection::ConvertDensity(
        3.0,
        PhaseSpaceMeasure::BjorkenXY(),
        PhaseSpaceMeasure::MandelstamQ2Y(),
        PhaseSpaceTopology::Scatter2to2,
        record);
    EXPECT_GT(converted, 0.0);
    EXPECT_THROW(
        siren::injection::ConvertDensity(
            3.0,
            PhaseSpaceMeasure::BjorkenXY(),
            PhaseSpaceMeasure::MandelstamQ2(),
            PhaseSpaceTopology::Scatter2to2,
            record),
        siren::utilities::MeasureCompatibilityError);
}

TEST(ScatteringMeasureConversion, RejectsDecayStyleLabBoost) {
    MultiChannelPhaseSpace mixture;
    mixture.channels = {
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleRest()),
        std::make_shared<ConstantChannel>(
            1.0, PhaseSpaceTopology::Scatter2to2,
            PhaseSpaceMeasure::SolidAngleLab())};
    mixture.weights = {0.5, 0.5};

    InteractionRecord record;
    EXPECT_THROW(mixture.Density(nullptr, record), std::runtime_error);
}


TEST(KinematicInjectionFailure, IsotropicTwoBodyRejectsSubThresholdDecay) {
    siren::injection::Isotropic2BodyChannel channel(0);
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_mass = 0.9;
    record.primary_momentum = {0.9, 0.0, 0.0, 0.0};
    record.secondary_masses = {0.5, 0.5};
    auto momenta_before = record.secondary_momenta;
    auto random = std::make_shared<siren::utilities::SIREN_random>(264575);

    EXPECT_DOUBLE_EQ(channel.Density(nullptr, record), 0.0);
    EXPECT_THROW(channel.Sample(random, nullptr, record),
                 siren::utilities::InjectionFailure);
    EXPECT_EQ(record.secondary_momenta, momenta_before);
}

TEST(KinematicInjectionFailure, IsotropicTwoBodyAcceptsExactThreshold) {
    siren::injection::Isotropic2BodyChannel channel(1);
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_mass = 1.0;
    record.primary_momentum = {1.0, 0.0, 0.0, 0.0};
    record.secondary_masses = {0.4, 0.6};
    auto random = std::make_shared<siren::utilities::SIREN_random>(271828);

    EXPECT_NO_THROW(channel.Sample(random, nullptr, record));
    EXPECT_NEAR(channel.Density(nullptr, record), 1.0 / (4.0 * M_PI), 1e-15);
    for (size_t i = 0; i < 2; ++i) {
        EXPECT_NEAR(record.secondary_momenta[i][0],
                    record.secondary_masses[i], 1e-15);
        EXPECT_DOUBLE_EQ(record.secondary_momenta[i][1], 0.0);
        EXPECT_DOUBLE_EQ(record.secondary_momenta[i][2], 0.0);
        EXPECT_DOUBLE_EQ(record.secondary_momenta[i][3], 0.0);
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


TEST(SharedInteractionRecordUtils, ReadWriteAndValidationUseOneLayout) {
    InteractionRecord record = TwoBodyDecayRecord();
    record.primary_momentum = {5.0, 1.0, 2.0, 3.0};
    record.interaction_vertex = {4.0, 5.0, 6.0};

    ASSERT_TRUE(siren::injection::detail::HasSecondaryStorage(record, 2));
    auto primary = siren::injection::detail::ReadPrimary(record);
    EXPECT_DOUBLE_EQ(primary.e, 5.0);
    EXPECT_DOUBLE_EQ(primary.p.GetX(), 1.0);
    EXPECT_DOUBLE_EQ(primary.p.GetY(), 2.0);
    EXPECT_DOUBLE_EQ(primary.p.GetZ(), 3.0);
    auto vertex = siren::injection::detail::ReadVertex(record);
    EXPECT_DOUBLE_EQ(vertex.GetX(), 4.0);
    EXPECT_DOUBLE_EQ(vertex.GetY(), 5.0);
    EXPECT_DOUBLE_EQ(vertex.GetZ(), 6.0);

    siren::injection::detail::WriteSecondary(record, 1, {
        7.0, siren::math::Vector3D(8.0, 9.0, 10.0)});
    auto secondary = siren::injection::detail::ReadSecondary(record, 1);
    EXPECT_DOUBLE_EQ(secondary.e, 7.0);
    EXPECT_DOUBLE_EQ(secondary.p.GetX(), 8.0);
    EXPECT_DOUBLE_EQ(secondary.p.GetY(), 9.0);
    EXPECT_DOUBLE_EQ(secondary.p.GetZ(), 10.0);

    record.secondary_momenta.pop_back();
    EXPECT_FALSE(siren::injection::detail::HasSecondaryStorage(record, 2));
    EXPECT_THROW(
        siren::injection::detail::RequireSecondaryStorage(
            record, 2, "test"),
        std::runtime_error);
}


TEST(TwoBodyLabAngle, RejectsNaNAndClampsRoundoffAtAngularBoundary) {
    constexpr double parent_mass = 2.0;
    constexpr double daughter_mass = 0.5;
    constexpr double other_mass = 0.5;
    constexpr double beta = 0.5;
    const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    const double p_rest = siren::injection::TwoBodyRestMomentum(
        parent_mass, daughter_mass, other_mass);
    const double E_rest = siren::injection::TwoBodyRestEnergy(
        parent_mass, daughter_mass, other_mass);

    auto invalid = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, E_rest, daughter_mass,
        std::numeric_limits<double>::quiet_NaN());
    for (auto const & solution : invalid) {
        EXPECT_FALSE(solution.valid);
        EXPECT_TRUE(std::isfinite(solution.cos_theta_rest));
        EXPECT_TRUE(std::isfinite(solution.p_lab));
        EXPECT_TRUE(std::isfinite(solution.jacobian));
    }

    auto rounded = siren::injection::SolveLabAngle(
        beta, gamma, p_rest, E_rest, daughter_mass,
        std::nextafter(1.0, 2.0));
    EXPECT_TRUE(rounded[0].valid || rounded[1].valid);
    for (auto const & solution : rounded) {
        if (!solution.valid) continue;
        EXPECT_TRUE(std::isfinite(solution.cos_theta_rest));
        EXPECT_TRUE(std::isfinite(solution.p_lab));
        EXPECT_TRUE(std::isfinite(solution.jacobian));
    }
}

TEST(SharedKinematics, StableBreakupMomentumBacksInjectionWrapper) {
    constexpr double mass_a = 0.001;
    constexpr double mass_b = 300.0;
    const double parent_mass = mass_a + mass_b + 1e-10;
    double canonical = siren::math::TwoBodyRestMomentum(
        parent_mass, mass_a, mass_b);

    EXPECT_TRUE(std::isfinite(canonical));
    EXPECT_GT(canonical, 0.0);
    EXPECT_DOUBLE_EQ(
        siren::injection::TwoBodyRestMomentum(
            parent_mass, mass_a, mass_b),
        canonical);
    EXPECT_DOUBLE_EQ(
        siren::injection::Kallen(7.0, 2.0, 1.0),
        siren::math::Kallen(7.0, 2.0, 1.0));
    EXPECT_DOUBLE_EQ(
        siren::math::TwoBodyRestMomentum(
            mass_a + mass_b, mass_a, mass_b),
        0.0);
}

TEST(SharedLorentzBoost, RkAdapterMatchesAnalyticBoost) {
    constexpr double parent_mass = 2.0;
    constexpr double gamma = 1.5;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    const double parent_energy = gamma * parent_mass;
    const double parent_pz = gamma * beta * parent_mass;
    constexpr double daughter_mass = 0.4;
    constexpr double rest_px = 0.3;
    constexpr double rest_py = -0.2;
    constexpr double rest_pz = 0.7;
    const double rest_energy = std::sqrt(
        daughter_mass * daughter_mass
        + rest_px * rest_px + rest_py * rest_py + rest_pz * rest_pz);

    auto lab = siren::injection::detail::BoostRestFrameToLab(
        parent_energy, 0.0, 0.0, parent_pz,
        rest_energy, rest_px, rest_py, rest_pz);
    EXPECT_NEAR(lab[0], gamma * (rest_energy + beta * rest_pz), 1e-14);
    EXPECT_NEAR(lab[1], rest_px, 1e-14);
    EXPECT_NEAR(lab[2], rest_py, 1e-14);
    EXPECT_NEAR(lab[3], gamma * (rest_pz + beta * rest_energy), 1e-14);
    EXPECT_NEAR(
        lab[0] * lab[0]
        - lab[1] * lab[1] - lab[2] * lab[2] - lab[3] * lab[3],
        daughter_mass * daughter_mass, 1e-13);
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

TEST(CoreMixture, InvalidComponentCannotBeMaskedByPositivePeer) {
    for (double bad : {-1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        MultiChannelPhaseSpace mixture({std::make_shared<ConstantChannel>(bad),
                                       std::make_shared<ConstantChannel>(10.0)});
        EXPECT_THROW(mixture.Density(nullptr, TwoBodyDecayRecord()),
                     siren::utilities::WeightCalculationError);
    }
    EXPECT_THROW(MultiChannelPhaseSpace({nullptr}), siren::utilities::ConfigurationError);
    EXPECT_THROW(MultiChannelPhaseSpace({std::make_shared<ConstantChannel>(),
        std::make_shared<ConstantChannel>()}, {1e308, 1e308}), siren::utilities::ConfigurationError);
}

TEST(CoreMeasure, RestCosThetaLiftsUniformAzimuthButCannotMarginalize) {
    auto integrated = PhaseSpaceMeasure::CosThetaRest();
    auto joint = PhaseSpaceMeasure::SolidAngleRest();
    auto topology = PhaseSpaceTopology::Decay2Body;
    auto record = TwoBodyDecayRecord();
    EXPECT_NEAR(siren::injection::ConvertDensity(0.5, integrated, joint, topology, record),
                1.0 / (4.0 * M_PI), 1e-16);
    EXPECT_TRUE(siren::injection::PhaseSpaceDensityConvertible(topology, integrated, joint));
    EXPECT_FALSE(siren::injection::PhaseSpaceDensityConvertible(topology, joint, integrated));
    EXPECT_THROW(siren::injection::ConvertDensity(0.5, joint, integrated, topology, record),
                 siren::utilities::MeasureCompatibilityError);
}

} // namespace
