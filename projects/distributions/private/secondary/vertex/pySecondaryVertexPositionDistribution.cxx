#include "SIREN/distributions/secondary/vertex/pySecondaryVertexPositionDistribution.h"

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <tuple>                                  // for tuple
#include <vector>                                 // for vector

#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/math/Vector3D.h"                    // for Vector3D
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

// Alias with no un-parenthesized comma so it can be passed as a single
// preprocessor macro argument to the SELF_OVERRIDE_PURE invocation below.
using SecondaryVertexPositionBounds = std::tuple<siren::math::Vector3D, siren::math::Vector3D>;

void pySecondaryVertexPositionDistribution::SampleVertex(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        SecondaryVertexPositionDistribution,
        void,
        SampleVertex,
        "SampleVertex",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pySecondaryVertexPositionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const {
    SELF_OVERRIDE_REF(
        self,
        SecondaryVertexPositionDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

double pySecondaryVertexPositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        SecondaryVertexPositionDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pySecondaryVertexPositionDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        SecondaryVertexPositionDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::string pySecondaryVertexPositionDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        SecondaryVertexPositionDistribution
    )
}

std::shared_ptr<SecondaryInjectionDistribution> pySecondaryVertexPositionDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        SecondaryVertexPositionDistribution,
        SecondaryInjectionDistribution
    )
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> pySecondaryVertexPositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        SecondaryVertexPositionDistribution,
        SecondaryVertexPositionBounds,
        InjectionBounds,
        "InjectionBounds",
        detector_model,
        interactions,
        interaction
    )
}

bool pySecondaryVertexPositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        SecondaryVertexPositionDistribution,
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

bool pySecondaryVertexPositionDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        SecondaryVertexPositionDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pySecondaryVertexPositionDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        SecondaryVertexPositionDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
