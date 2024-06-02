#pragma once
#ifndef SIREN_pyDarkNewsDecay_H
#define SIREN_pyDarkNewsDecay_H

#include <set>                                    // for set
#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t
#include <stdexcept>                              // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include "SIREN/interactions/DarkNewsDecay.h"

#include "SIREN/interactions/Decay.h"   // for Decay
#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Trampoline class for DarkNewsDecay
class pyDarkNewsDecay : public DarkNewsDecay {
public:
    using DarkNewsDecay::DarkNewsDecay;
    pyDarkNewsDecay(DarkNewsDecay && parent);
    pyDarkNewsDecay(DarkNewsDecay const & parent);

    double TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const override;
    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override;
    double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override;
    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override;
    void SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<std::string> DensityVariables() const override;
    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    Pybind11TrampolineCerealMethods(DarkNewsDecay, pyDarkNewsDecay);

}; // class pyDarkNewsDecay

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::pyDarkNewsDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::pyDarkNewsDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::DarkNewsDecay, siren::interactions::pyDarkNewsDecay);

#endif // SIREN_pyDarkNewsDecay_H

