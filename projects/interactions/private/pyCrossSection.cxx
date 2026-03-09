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
        const CrossSection * ref;
        if(self) {
            ref = self.cast<CrossSection *>();
        } else {
            ref = this;
        }
        do {
            do {
                pybind11::gil_scoped_acquire gil;
                pybind11::function override
                    = pybind11::get_override(static_cast<const CrossSection *>(ref), "SampleFinalState");
                if (override) {
                auto o = override.operator()<pybind11::return_value_policy::reference>(record, random);
                    if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
                        static pybind11::detail::override_caster_t<void> caster;
                        return pybind11::detail::cast_ref<void>(std::move(o), caster);
                    }
                    return pybind11::detail::cast_safe<void>(std::move(o));
                }
            } while (false);
            pybind11::pybind11_fail(
                "Tried to call pure virtual function \"CrossSection::SampleFinalState\"");
        } while (false);
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

