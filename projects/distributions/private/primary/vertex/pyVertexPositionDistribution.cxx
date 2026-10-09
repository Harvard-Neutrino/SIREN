#include "SIREN/distributions/primary/vertex/pyVertexPositionDistribution.h"

#include <memory>
#include <set>                                 // for shared_ptr
#include <string>                                 // for string
#include <tuple>                                  // for tuple
#include <vector>                                 // for vector

#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"

#include "SIREN/dataclasses/InteractionRecord.h"    // for InteractionRecord
#include "SIREN/detector/DetectorModel.h"           // for DetectorModel
#include "SIREN/interactions/InteractionCollection.h" // for InteractionCollection
#include "SIREN/math/Vector3D.h"                    // for Vector3D
#include "SIREN/utilities/Pybind11Trampoline.h"     // for SELF_OVERRIDE macros
#include "SIREN/utilities/Random.h"                 // for SIREN_random

namespace siren {
namespace distributions {

// Alias with no un-parenthesized comma so it can be passed as a single
// preprocessor macro argument to the SELF_OVERRIDE invocations below.
using VertexPositionBounds = std::tuple<siren::math::Vector3D, siren::math::Vector3D>;

std::tuple<siren::math::Vector3D, siren::math::Vector3D> pyVertexPositionDistribution::SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_PURE_REF(
        self,
        VertexPositionDistribution,
        VertexPositionBounds,
        SamplePosition,
        "SamplePosition",
        rand,
        detector_model,
        interactions,
        record
    )
}

void pyVertexPositionDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    SELF_OVERRIDE_REF(
        self,
        VertexPositionDistribution,
        void,
        Sample,
        "Sample",
        rand,
        detector_model,
        interactions,
        record
    )
}

double pyVertexPositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        VertexPositionDistribution,
        double,
        GenerationProbability,
        "GenerationProbability",
        detector_model,
        interactions,
        record
    )
}

std::vector<std::string> pyVertexPositionDistribution::DensityVariables() const {
    SELF_OVERRIDE(
        self,
        VertexPositionDistribution,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

std::set<DistributionVariable> pyVertexPositionDistribution::SetVariables() const {
    SELF_OVERRIDE(
        self,
        VertexPositionDistribution,
        std::set<DistributionVariable>,
        SetVariables,
        "SetVariables"
    )
}

std::set<DistributionVariable> pyVertexPositionDistribution::RequiredVariables() const {
    SELF_OVERRIDE(
        self,
        VertexPositionDistribution,
        std::set<DistributionVariable>,
        RequiredVariables,
        "RequiredVariables"
    )
}

std::string pyVertexPositionDistribution::Name() const {
    SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(
        self,
        VertexPositionDistribution
    )
}

std::shared_ptr<PrimaryInjectionDistribution> pyVertexPositionDistribution::clone() const {
    SELF_OVERRIDE_CLONE_REQUIRED(
        self,
        VertexPositionDistribution,
        PrimaryInjectionDistribution
    )
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> pyVertexPositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        VertexPositionDistribution,
        VertexPositionBounds,
        InjectionBounds,
        "InjectionBounds",
        detector_model,
        interactions,
        interaction
    )
}

bool pyVertexPositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    SELF_OVERRIDE(
        self,
        VertexPositionDistribution,
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

bool pyVertexPositionDistribution::equal(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(
        self,
        VertexPositionDistribution,
        WeightableDistribution,
        distribution
    )
}

bool pyVertexPositionDistribution::less(WeightableDistribution const & distribution) const {
    SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(
        self,
        VertexPositionDistribution,
        WeightableDistribution,
        distribution
    )
}

} // namespace distributions
} // namespace siren
