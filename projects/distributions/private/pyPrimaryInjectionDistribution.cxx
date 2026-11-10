#include "SIREN/distributions/pyPrimaryInjectionDistribution.h"

#include <memory>
#include <set>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector

#include "SIREN/distributions/Distributions.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

void pyPrimaryInjectionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        PrimaryInjectionDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

double pyPrimaryInjectionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        PrimaryInjectionDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyPrimaryInjectionDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryInjectionDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::set<DistributionVariable> pyPrimaryInjectionDistribution::SetVariables() const {
    SELF_OVERRIDE_PURE(
        self,
        PrimaryInjectionDistribution,
        std::set<DistributionVariable>,
        SetVariables,
        "SetVariables"
    )
}

std::set<DistributionVariable> pyPrimaryInjectionDistribution::RequiredVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryInjectionDistribution,
        std::set<DistributionVariable>,
        RequiredVariables,
        "RequiredVariables"
    )
}

std::string pyPrimaryInjectionDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        PrimaryInjectionDistribution
    )
}

std::shared_ptr<PrimaryInjectionDistribution> pyPrimaryInjectionDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        PrimaryInjectionDistribution,
        PrimaryInjectionDistribution
    )
}

bool pyPrimaryInjectionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        PrimaryInjectionDistribution,
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

bool pyPrimaryInjectionDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        PrimaryInjectionDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyPrimaryInjectionDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        PrimaryInjectionDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
