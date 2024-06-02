#include <memory>
#include <string>
#include <vector>                                       // for vector
#include <cstdint>                                      // for uint32_t
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/interactions/pyDarkNewsCrossSection.h"  // for DarkNewsCrossSection
#include "SIREN/interactions/DarkNewsCrossSection.h"  // for DarkNewsCrossSection

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

pyDarkNewsCrossSection::pyDarkNewsCrossSection(DarkNewsCrossSection && parent) : DarkNewsCrossSection(std::move(parent)) {
    self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsCrossSection)))));
}
pyDarkNewsCrossSection::pyDarkNewsCrossSection(DarkNewsCrossSection const & parent) : DarkNewsCrossSection(parent) {
    self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsCrossSection)))));
}

double pyDarkNewsCrossSection::TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE(
            self,
            CrossSection,
            double,
            TotalCrossSectionAllFinalStates,
            "TotalCrossSectionAllFinalStates",
            std::cref(record)
            )
}

double pyDarkNewsCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            std::cref(interaction)
            )
}

double pyDarkNewsCrossSection::TotalCrossSection(siren::dataclasses::ParticleType primary, double energy, siren::dataclasses::ParticleType target) const {
    SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            primary,
            energy,
            target
            )
}

double pyDarkNewsCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            std::cref(interaction)
            )
}

double pyDarkNewsCrossSection::DifferentialCrossSection(siren::dataclasses::ParticleType primary, siren::dataclasses::ParticleType target, double energy, double Q2) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            primary,
            target,
            energy,
            Q2
            )
}

double pyDarkNewsCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            std::cref(interaction)
            )
}

double pyDarkNewsCrossSection::Q2Min(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Min,
            "Q2Min",
            std::cref(interaction)
            )
}

double pyDarkNewsCrossSection::Q2Max(dataclasses::InteractionRecord const & interaction) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Max,
            "Q2Max",
            std::cref(interaction)
            )
}

double pyDarkNewsCrossSection::TargetMass(dataclasses::ParticleType const & target_type) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TargetMass,
            "TargetMass",
            std::cref(target_type)
            )
}

std::vector<double> pyDarkNewsCrossSection::SecondaryMasses(std::vector<dataclasses::ParticleType> const & secondary_types) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            std::vector<double>,
            SecondaryMasses,
            "SecondaryMasses",
            std::cref(secondary_types)
            )
}

std::vector<double> pyDarkNewsCrossSection::SecondaryHelicities(dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            std::vector<double>,
            SecondaryHelicities,
            "SecondaryHelicities",
            std::cref(record)
            )
}

void pyDarkNewsCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            std::ref(record),
            random
            )
}

std::vector<siren::dataclasses::ParticleType> pyDarkNewsCrossSection::GetPossibleTargets() const {
    SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
            )
}

std::vector<siren::dataclasses::ParticleType> pyDarkNewsCrossSection::GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const {
    SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
            )
}

std::vector<siren::dataclasses::ParticleType> pyDarkNewsCrossSection::GetPossiblePrimaries() const {
    SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
            )
}

std::vector<siren::dataclasses::InteractionSignature> pyDarkNewsCrossSection::GetPossibleSignatures() const {
    SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
            )
}

std::vector<siren::dataclasses::InteractionSignature> pyDarkNewsCrossSection::GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const {
    SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
            )
}

double pyDarkNewsCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            std::cref(record)
            )
}

} // namespace interactions
} // namespace siren

