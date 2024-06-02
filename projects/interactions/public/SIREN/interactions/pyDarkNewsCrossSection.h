#pragma once
#ifndef SIREN_pyDarkNewsCrossSection_H
#define SIREN_pyDarkNewsCrossSection_H

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

#include "SIREN/interactions/DarkNewsCrossSection.h"  // for DarkNewsCrossSection

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Trampoline class for DarkNewsCrossSection
class pyDarkNewsCrossSection : public DarkNewsCrossSection {
public:
    using DarkNewsCrossSection::DarkNewsCrossSection;
    pyDarkNewsCrossSection(DarkNewsCrossSection && parent);
    pyDarkNewsCrossSection(DarkNewsCrossSection const & parent);
    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override;
    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override;
    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy, siren::dataclasses::ParticleType target) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override;
    double DifferentialCrossSection(siren::dataclasses::ParticleType primary, siren::dataclasses::ParticleType target, double energy, double Q2) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override;
    double Q2Min(dataclasses::InteractionRecord const & interaction) const override;
    double Q2Max(dataclasses::InteractionRecord const & interaction) const override;
    double TargetMass(dataclasses::ParticleType const & target_type) const override;
    std::vector<double> SecondaryMasses(std::vector<dataclasses::ParticleType> const & secondary_types) const override;
    std::vector<double> SecondaryHelicities(dataclasses::InteractionRecord const & record) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;
    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    Pybind11TrampolineCerealMethods(DarkNewsCrossSection, pyDarkNewsCrossSection);
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::pyDarkNewsCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::pyDarkNewsCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::DarkNewsCrossSection, siren::interactions::pyDarkNewsCrossSection);

//CEREAL_FORCE_DYNAMIC_INIT(pyDarkNewsCrossSection);

#endif // SIREN_pyDarkNewsCrossSection_H

