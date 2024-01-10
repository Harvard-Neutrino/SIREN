#pragma once
#ifndef LI_ElasticScattering_H
#define LI_ElasticScattering_H

#include <set>                                          // for set
#include <memory>
#include <string>
#include <vector>                                       // for vector
#include <cstdint>                                      // for uint32_t
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/set.hpp>
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

// For details, see appendix A of 1906.00111v4
class ElasticScattering : public CrossSection {
friend cereal::access;
protected:
private:
    const double CLR = 0.2334; // at one loop
    const std::set<LI::dataclasses::Particle::ParticleType> primary_types = {LI::dataclasses::Particle::ParticleType::NuE, LI::dataclasses::Particle::ParticleType::NuMu};
public:
    ElasticScattering() {};
    ElasticScattering(std::set<LI::dataclasses::Particle::ParticleType> const & primary_types) : primary_types(primary_types) {};
    virtual bool equal(CrossSection const & other) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, double y) const;
    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const override;

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const override;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const override;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("ElasticScattering only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<ElasticScattering> & construct, std::uint32_t version) {
        if(version == 0) {
            std::set<LI::dataclasses::Particle::ParticleType> _primary_types;
            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            construct(_primary_types);
            archive(cereal::virtual_base_class<CrossSection>(construct.ptr()));
        } else {
            throw std::runtime_error("ElasticScattering only supports version <= 0!");
        }
    }
};

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::crosssections::ElasticScattering, 0);
CEREAL_REGISTER_TYPE(LI::crosssections::ElasticScattering);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::crosssections::CrossSection, LI::crosssections::ElasticScattering);

#endif // LI_ElasticScattering_H
