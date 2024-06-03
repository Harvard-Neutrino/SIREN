#pragma once
#ifndef SIREN_pyCrossSection_H
#define SIREN_pyCrossSection_H

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

#include "SIREN/interactions/CrossSection.h"

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline
#include "SIREN/utilities/Random.h" // for SIREN_random

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Trampoline class for CrossSection
class pyCrossSection : public CrossSection {
public:
    using CrossSection::CrossSection;
    pyCrossSection(CrossSection && parent) : CrossSection(std::move(parent)) {}

    bool equal(CrossSection const & other) const override;
    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override;
    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;
    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    std::vector<std::string> DensityVariables() const override;

    Pybind11TrampolineCerealMethods(CrossSection, pyCrossSection);
}; // class pyCrossSection

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::pyCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::pyCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::pyCrossSection);

#endif // SIREN_pyCrossSection_H

