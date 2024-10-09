#pragma once
#ifndef SIREN_DMesonELoss_H
#define SIREN_DMesonELoss_H

#include <set>                                                // for set
#include <map>                                                // for map
#include <memory>
#include <vector>                                             // for vector
#include <cstdint>                                            // for uint32_t
#include <utility>                                            // for pair
#include <algorithm>
#include <stdexcept>                                          // for runtime...

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include "SIREN/interactions/CrossSection.h"        // for CrossSe...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class DMesonELoss : public CrossSection {
friend cereal::access;
private:
    std::set<siren::dataclasses::Particle::ParticleType> primary_types_ = {siren::dataclasses::Particle::ParticleType::D0, siren::dataclasses::Particle::ParticleType::DPlus};
    std::set<siren::dataclasses::Particle::ParticleType> target_types_ = {siren::dataclasses::Particle::ParticleType::PPlus};

public:
    DMesonELoss();

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(siren::dataclasses::Particle::ParticleType primary, double energy) const;
    // double TotalCrossSection(siren::dataclasses::Particle::ParticleType primary, double energy, siren::dataclasses::Particle::ParticleType target) const override;
    
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random> random) const override;

    std::vector<siren::dataclasses::Particle::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::Particle::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::Particle::ParticleType primary_type, siren::dataclasses::Particle::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DMesonELoss only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(cereal::virtual_base_class<CrossSection>(this));
            InitializeSignatures();
        } else {
            throw std::runtime_error("DMesonELoss only supports version <= 0!");
        }
    }
private:
    void InitializeSignatures();
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::DMesonELoss, 0);
CEREAL_REGISTER_TYPE(siren::interactions::DMesonELoss);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::DMesonELoss);

#endif // SIREN_DMesonELoss_H
