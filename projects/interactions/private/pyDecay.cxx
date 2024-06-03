#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/pyDecay.h"

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/dataclasses/InteractionRecord.h" // for InteractionRecord
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline
#include "SIREN/utilities/Random.h" // for SIREN_random

namespace siren {
namespace interactions {

bool pyDecay::equal(Decay const & other) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        bool,
        equal,
        "equal",
        other
    )
}

double pyDecay::TotalDecayLength(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
        self,
        Decay,
        double,
        TotalDecayLength,
        "TotalDecayLength",
        interaction
    )
}

double pyDecay::TotalDecayLengthForFinalState(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
        self,
        Decay,
        double,
        TotalDecayLengthForFinalState,
        "TotalDecayLengthForFinalState",
        interaction
    )
}

double pyDecay::TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        double,
        TotalDecayWidth,
        "TotalDecayWidth",
        interaction
    )
}

double pyDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        double,
        TotalDecayWidthForFinalState,
        "TotalDecayWidthForFinalState",
        interaction
    )
}

double pyDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        double,
        TotalDecayWidth,
        "TotalDecayWidth",
        primary
    )
}

double pyDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        double,
        DifferentialDecayWidth,
        "DifferentialDecayWidth",
        interaction
    )
}

void pyDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        void,
        SampleFinalState,
        "SampleFinalState",
        record,
        random
    )
}

std::vector<siren::dataclasses::InteractionSignature> pyDecay::GetPossibleSignatures() const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        std::vector<siren::dataclasses::InteractionSignature>,
        GetPossibleSignatures,
        "GetPossibleSignatures"
    )
}

std::vector<siren::dataclasses::InteractionSignature> pyDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary_type) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        std::vector<siren::dataclasses::InteractionSignature>,
        GetPossibleSignaturesFromParents,
        "GetPossibleSignaturesFromParents",
        primary_type
    )
}

std::vector<std::string> pyDecay::DensityVariables() const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        std::vector<std::string>,
        DensityVariables,
        "DensityVariables"
    )
}

double pyDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        Decay,
        double,
        FinalStateProbability,
        "FinalStateProbability",
        record
    )
}

} // namespace interactions
} // namespace siren

