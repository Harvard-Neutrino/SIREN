#include "SIREN/distributions/primary/area/pyPrimaryAreaDistribution.h"

#include <memory>
#include <set>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector

#include "SIREN/distributions/primary/area/PrimaryAreaDistribution.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/math/Vector3D.h"                    // for Vector3D
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

siren::math::Vector3D pyPrimaryAreaDistribution::SamplePointOfClosestApproach(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        PrimaryAreaDistribution,
        siren::math::Vector3D,
        SamplePointOfClosestApproach,
        "SamplePointOfClosestApproach",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyPrimaryAreaDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_REF(
        self,
        PrimaryAreaDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

double pyPrimaryAreaDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        PrimaryAreaDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyPrimaryAreaDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryAreaDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::set<DistributionVariable> pyPrimaryAreaDistribution::SetVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryAreaDistribution,
        std::set<DistributionVariable>,
        SetVariables,
        "SetVariables"
    )
}

std::set<DistributionVariable> pyPrimaryAreaDistribution::RequiredVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryAreaDistribution,
        std::set<DistributionVariable>,
        RequiredVariables,
        "RequiredVariables"
    )
}

std::string pyPrimaryAreaDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        PrimaryAreaDistribution
    )
}

std::shared_ptr<PrimaryInjectionDistribution> pyPrimaryAreaDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        PrimaryAreaDistribution,
        PrimaryInjectionDistribution
    )
}

bool pyPrimaryAreaDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        PrimaryAreaDistribution,
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

bool pyPrimaryAreaDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        PrimaryAreaDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyPrimaryAreaDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        PrimaryAreaDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
