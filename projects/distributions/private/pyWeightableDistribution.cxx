#include "SIREN/distributions/pyWeightableDistribution.h"

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector

#include "SIREN/distributions/Distributions.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros

namespace siren {
namespace distributions {

double pyWeightableDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        WeightableDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyWeightableDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        WeightableDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::string pyWeightableDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        WeightableDistribution
    )
}

bool pyWeightableDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        WeightableDistribution,
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

bool pyWeightableDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        WeightableDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyWeightableDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        WeightableDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
