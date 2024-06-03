#pragma once
#ifndef SIREN_pyDecay_H
#define SIREN_pyDecay_H

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

// Trampoline class for Decay
class pyDecay : public Decay {
public:
    using Decay::Decay;
    pyDecay(Decay && parent) : Decay(std::move(parent)) {}

    bool equal(Decay const & other) const override;
    double TotalDecayLength(dataclasses::InteractionRecord const & interaction) const override;
    double TotalDecayLengthForFinalState(dataclasses::InteractionRecord const & interaction) const override;
    double TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const override;
    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override;
    double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override;
    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<std::string> DensityVariables() const override;
    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    Pybind11TrampolineCerealMethods(Decay, pyDecay);
}; // class pyDecay

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::pyDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::pyDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::pyDecay);

#endif // SIREN_Decay_H

