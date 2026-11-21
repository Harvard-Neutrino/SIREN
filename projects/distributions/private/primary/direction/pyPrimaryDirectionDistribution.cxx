#include "SIREN/distributions/primary/direction/pyPrimaryDirectionDistribution.h"

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector

#include "SIREN/distributions/primary/direction/PrimaryDirectionDistribution.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/math/Vector3D.h"                    // for Vector3D
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

siren::math::Vector3D pyPrimaryDirectionDistribution::SampleDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        PrimaryDirectionDistribution,
        siren::math::Vector3D,
        SampleDirection,
        "SampleDirection",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyPrimaryDirectionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_REF(
        self,
        PrimaryDirectionDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

double pyPrimaryDirectionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        PrimaryDirectionDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyPrimaryDirectionDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryDirectionDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::string pyPrimaryDirectionDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        PrimaryDirectionDistribution
    )
}

std::shared_ptr<PrimaryInjectionDistribution> pyPrimaryDirectionDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        PrimaryDirectionDistribution,
        PrimaryInjectionDistribution
    )
}

bool pyPrimaryDirectionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        PrimaryDirectionDistribution,
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

bool pyPrimaryDirectionDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        PrimaryDirectionDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyPrimaryDirectionDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        PrimaryDirectionDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
