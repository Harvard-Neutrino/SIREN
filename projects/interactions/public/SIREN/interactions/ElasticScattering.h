#pragma once
#ifndef SIREN_ElasticScattering_H
#define SIREN_ElasticScattering_H

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

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// For details, see appendix A of 1906.00111v4
class ElasticScattering : public CrossSection {
friend cereal::access;
protected:
private:
    const double CLR = 0.2334; // at one loop
    const std::set<siren::dataclasses::ParticleType> primary_types = {siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu};
public:
    ElasticScattering() {};
    ElasticScattering(std::set<siren::dataclasses::ParticleType> const & primary_types) : primary_types(primary_types) {};
    virtual bool equal(CrossSection const & other) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(siren::dataclasses::ParticleType primary_type, double primary_energy, double y) const;
    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy, siren::dataclasses::ParticleType target) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;

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
            std::set<siren::dataclasses::ParticleType> _primary_types;
            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            construct(_primary_types);
            archive(cereal::virtual_base_class<CrossSection>(construct.ptr()));
        } else {
            throw std::runtime_error("ElasticScattering only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::ElasticScattering, 0);
CEREAL_REGISTER_TYPE(siren::interactions::ElasticScattering);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::ElasticScattering);

#endif // SIREN_ElasticScattering_H
