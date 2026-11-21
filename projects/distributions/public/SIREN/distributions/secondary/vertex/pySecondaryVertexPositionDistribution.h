#pragma once
#ifndef SIREN_pySecondaryVertexPositionDistribution_H
#define SIREN_pySecondaryVertexPositionDistribution_H

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <tuple>                                  // for tuple
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"

#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class SecondaryDistributionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace math { class Vector3D; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

// Trampoline class for SecondaryVertexPositionDistribution
class pySecondaryVertexPositionDistribution : public SecondaryVertexPositionDistribution {
public:
    pySecondaryVertexPositionDistribution() = default;

    void SampleVertex(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const override;
    void Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::SecondaryDistributionRecord & record) const override;
    double GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const override;
    std::vector<std::string> DensityVariables() const override;
    std::string Name() const override;
    std::shared_ptr<SecondaryInjectionDistribution> clone() const override;
    std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const override;
    bool AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const override;
    bool equal(WeightableDistribution const & distribution) const override;
    bool less(WeightableDistribution const & distribution) const override;

    Pybind11TrampolineCerealMethods(SecondaryVertexPositionDistribution, pySecondaryVertexPositionDistribution);
}; // class pySecondaryVertexPositionDistribution

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::pySecondaryVertexPositionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::pySecondaryVertexPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::distributions::SecondaryVertexPositionDistribution, siren::distributions::pySecondaryVertexPositionDistribution);

#endif // SIREN_pySecondaryVertexPositionDistribution_H
