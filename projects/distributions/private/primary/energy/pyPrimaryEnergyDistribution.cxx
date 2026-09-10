#include "SIREN/distributions/primary/energy/pyPrimaryEnergyDistribution.h"

#include <memory>
#include <set>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector

#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

double pyPrimaryEnergyDistribution::SampleEnergy(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        PrimaryEnergyDistribution,
        double,
        SampleEnergy,
        "SampleEnergy",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyPrimaryEnergyDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_REF(
        self,
        PrimaryEnergyDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyPrimaryEnergyDistribution::SetNormalization(double norm) {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        void,
        SetNormalization,
        "SetNormalization",
        norm
    )
}

double pyPrimaryEnergyDistribution::GetNormalization() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        double,
        GetNormalization,
        "GetNormalization"
    )
}

bool pyPrimaryEnergyDistribution::IsNormalizationSet() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        bool,
        IsNormalizationSet,
        "IsNormalizationSet"
    )
}

double pyPrimaryEnergyDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        PrimaryEnergyDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyPrimaryEnergyDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::set<DistributionVariable> pyPrimaryEnergyDistribution::SetVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        std::set<DistributionVariable>,
        SetVariables,
        "SetVariables"
    )
}

std::set<DistributionVariable> pyPrimaryEnergyDistribution::RequiredVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        std::set<DistributionVariable>,
        RequiredVariables,
        "RequiredVariables"
    )
}

std::string pyPrimaryEnergyDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        PrimaryEnergyDistribution
    )
}

std::shared_ptr<PrimaryInjectionDistribution> pyPrimaryEnergyDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        PrimaryEnergyDistribution,
        PrimaryInjectionDistribution
    )
}

bool pyPrimaryEnergyDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDistribution,
        bool,
        AreEquivalent,
        "AreEquivalent",
        detector_model,
        interactions,
        distribution,
        second_detector_model,
        second_interactions
    )
}

bool pyPrimaryEnergyDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        PrimaryEnergyDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyPrimaryEnergyDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        PrimaryEnergyDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
