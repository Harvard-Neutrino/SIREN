#include "SIREN/distributions/primary/energy_direction/pyPrimaryEnergyDirectionDistribution.h"

#include <memory>
#include <set>                                 // for shared_ptr
#include <string>                                 // for string
#include <utility>                                // for pair
#include <vector>                                 // for vector

#include "SIREN/distributions/primary/energy_direction/PrimaryEnergyDirectionDistribution.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/math/Vector3D.h"                    // for Vector3D
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

// Alias with no un-parenthesized comma so it can be passed as a single
// preprocessor macro argument to the SELF_OVERRIDE_PURE_REF invocation below.
using PrimaryEnergyDirectionSample = std::pair<double, siren::math::Vector3D>;

std::pair<double, siren::math::Vector3D> pyPrimaryEnergyDirectionDistribution::SampleEnergyAndDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        PrimaryEnergyDirectionDistribution,
        PrimaryEnergyDirectionSample,
        SampleEnergyAndDirection,
        "SampleEnergyAndDirection",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyPrimaryEnergyDirectionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_REF(
        self,
        PrimaryEnergyDirectionDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyPrimaryEnergyDirectionDistribution::SetNormalization(double norm) {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
        void,
        SetNormalization,
        "SetNormalization",
        norm
    )
}

double pyPrimaryEnergyDirectionDistribution::GetNormalization() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
        double,
        GetNormalization,
        "GetNormalization"
    )
}

bool pyPrimaryEnergyDirectionDistribution::IsNormalizationSet() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
        bool,
        IsNormalizationSet,
        "IsNormalizationSet"
    )
}

double pyPrimaryEnergyDirectionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        PrimaryEnergyDirectionDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyPrimaryEnergyDirectionDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::set<DistributionVariable> pyPrimaryEnergyDirectionDistribution::SetVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
        std::set<DistributionVariable>,
        SetVariables,
        "SetVariables"
    )
}

std::set<DistributionVariable> pyPrimaryEnergyDirectionDistribution::RequiredVariables() const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
        std::set<DistributionVariable>,
        RequiredVariables,
        "RequiredVariables"
    )
}

std::string pyPrimaryEnergyDirectionDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        PrimaryEnergyDirectionDistribution
    )
}

std::shared_ptr<PrimaryInjectionDistribution> pyPrimaryEnergyDirectionDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        PrimaryEnergyDirectionDistribution,
        PrimaryInjectionDistribution
    )
}

bool pyPrimaryEnergyDirectionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        PrimaryEnergyDirectionDistribution,
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

bool pyPrimaryEnergyDirectionDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        PrimaryEnergyDirectionDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyPrimaryEnergyDirectionDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        PrimaryEnergyDirectionDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
