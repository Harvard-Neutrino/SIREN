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

#include "SIREN/interactions/pyCrossSection.h"
#include "SIREN/interactions/CrossSection.h"

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/dataclasses/InteractionRecord.h" // for InteractionRecord
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline
#include "SIREN/utilities/Random.h" // for SIREN_random

namespace siren {
namespace interactions {

bool pyCrossSection::equal(CrossSection const & other) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        bool,
        equal,
        "equal",
        other
    )
}

double pyCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        double,
        TotalCrossSection,
        "TotalCrossSection",
        interaction
    )
}

double pyCrossSection::TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE(
        self,
        CrossSection,
        double,
        TotalCrossSectionAllFinalStates,
        "TotalCrossSectionAllFinalStates",
        record
    )
}

double pyCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        double,
        DifferentialCrossSection,
        "DifferentialCrossSection",
        interaction
    )
}

double pyCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        double,
        InteractionThreshold,
        "InteractionThreshold",
        interaction
    )
}

void pyCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        void,
        SampleFinalState,
        "SampleFinalState",
        record,
        random
    )
}

std::vector<siren::dataclasses::ParticleType> pyCrossSection::GetPossibleTargets() const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        std::vector<siren::dataclasses::ParticleType>,
        GetPossibleTargets,
        "GetPossibleTargets"
    )
}

std::vector<siren::dataclasses::ParticleType> pyCrossSection::GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        std::vector<siren::dataclasses::ParticleType>,
        GetPossibleTargetsFromPrimary,
        "GetPossibleTargetsFromPrimary",
        primary_type
    )
}

std::vector<siren::dataclasses::ParticleType> pyCrossSection::GetPossiblePrimaries() const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        std::vector<siren::dataclasses::ParticleType>,
        GetPossiblePrimaries,
        "GetPossiblePrimaries"
    )
}

std::vector<siren::dataclasses::InteractionSignature> pyCrossSection::GetPossibleSignatures() const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        std::vector<siren::dataclasses::InteractionSignature>,
        GetPossibleSignatures,
        "GetPossibleSignatures"
    )
}

std::vector<siren::dataclasses::InteractionSignature> pyCrossSection::GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        std::vector<siren::dataclasses::InteractionSignature>,
        GetPossibleSignaturesFromParents,
        "GetPossibleSignaturesFromParents",
        primary_type,
        target_type
    )
}

double pyCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE_PURE(
        self,
        CrossSection,
        double,
        FinalStateProbability,
        "FinalStateProbability",
        record
    )
}

std::vector<std::string> pyCrossSection::DensityVariables() const {
    PYBIND11_OVERRIDE_PURE(
        std::vector<std::string>,
        CrossSection,
        DensityVariables
    );
}

} // namespace interactions
} // namespace siren

