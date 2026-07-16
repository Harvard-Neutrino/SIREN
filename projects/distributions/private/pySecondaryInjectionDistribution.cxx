#include "SIREN/distributions/pySecondaryInjectionDistribution.h"

#include <memory>                                 // for shared_ptr
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

void pySecondaryInjectionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        SecondaryInjectionDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

double pySecondaryInjectionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        SecondaryInjectionDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pySecondaryInjectionDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        SecondaryInjectionDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::string pySecondaryInjectionDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        SecondaryInjectionDistribution
    )
}

std::shared_ptr<SecondaryInjectionDistribution> pySecondaryInjectionDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        SecondaryInjectionDistribution,
        SecondaryInjectionDistribution
    )
}

bool pySecondaryInjectionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        SecondaryInjectionDistribution,
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

bool pySecondaryInjectionDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        SecondaryInjectionDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pySecondaryInjectionDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        SecondaryInjectionDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
