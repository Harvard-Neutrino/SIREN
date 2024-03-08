#pragma once
#ifndef LI_DummyCrossSection_H
#define LI_DummyCrossSection_H

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

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace SI { namespace dataclasses { struct InteractionSignature; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace interactions {

class DummyCrossSection : public CrossSection {
friend cereal::access;
private:
public:
    DummyCrossSection();

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(SI::dataclasses::ParticleType primary, double energy, SI::dataclasses::ParticleType target) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<SI::utilities::LI_random> random) const override;

    std::vector<SI::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<SI::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(SI::dataclasses::ParticleType primary_type) const override;
    std::vector<SI::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(SI::dataclasses::ParticleType primary_type, SI::dataclasses::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DummyCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DummyCrossSection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace SI

CEREAL_CLASS_VERSION(SI::interactions::DummyCrossSection, 0);
CEREAL_REGISTER_TYPE(SI::interactions::DummyCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::interactions::CrossSection, SI::interactions::DummyCrossSection);

#endif // LI_DummyCrossSection_H
