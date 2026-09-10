#pragma once
#ifndef SIREN_pyPrimaryInjectionDistribution_H
#define SIREN_pyPrimaryInjectionDistribution_H

#include <memory>
#include <set>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/Distributions.h"

#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class PrimaryDistributionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

// Trampoline class for PrimaryInjectionDistribution
class pyPrimaryInjectionDistribution : public PrimaryInjectionDistribution {
public:
    pyPrimaryInjectionDistribution() = default;

    void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const override;
    double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    std::vector<std::string> DensityVariables() const override;
    std::set<DistributionVariable> SetVariables() const override;
    std::set<DistributionVariable> RequiredVariables() const override;
    std::string Name() const override;
    std::shared_ptr<PrimaryInjectionDistribution> clone() const override;
    bool AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const override;
    bool equal(WeightableDistribution const & distribution) const override;
    bool less(WeightableDistribution const & distribution) const override;

    Pybind11TrampolineCerealMethods(PrimaryInjectionDistribution, pyPrimaryInjectionDistribution);
}; // class pyPrimaryInjectionDistribution

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::pyPrimaryInjectionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::pyPrimaryInjectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::PrimaryInjectionDistribution, siren::distributions::pyPrimaryInjectionDistribution);

#endif // SIREN_pyPrimaryInjectionDistribution_H
