#pragma once
#ifndef LI_DarkNewsCrossSection_H
#define LI_DarkNewsCrossSection_H

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

#include "LeptonInjector/crosssections/CrossSection.h"  // for CrossSection
#include "LeptonInjector/dataclasses/Particle.h"        // for Particle

namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace dataclasses { struct InteractionSignature; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace crosssections {

class DarkNewsCrossSection : public CrossSection {
friend cereal::access;
private:
public:
    DarkNewsCrossSection();

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const override = 0; // Assumes Python-side implementation
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary, LI::dataclasses::Particle::ParticleType target, double energy, double Q2) const; 
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    double Q2Min(dataclasses::InteractionRecord const &) const;
    double Q2Max(dataclasses::InteractionRecord const &) const;
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random> random) const override;

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const override;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const override = 0; // Assumes Python-side implementation
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const override = 0; // Assumes Python-side implementation

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DarkNewsCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DarkNewsCrossSection only supports version <= 0!");
        }
    }
};

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::crosssections::DarkNewsCrossSection, 0);
CEREAL_REGISTER_TYPE(LI::crosssections::DarkNewsCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::crosssections::CrossSection, LI::crosssections::DarkNewsCrossSection);

#endif // LI_DarkNewsCrossSection_H