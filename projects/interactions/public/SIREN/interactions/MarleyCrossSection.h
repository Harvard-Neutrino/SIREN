#pragma once
#ifndef SIREN_MarleyCrossSection_H
#define SIREN_MarleyCrossSection_H

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

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/utilities/Random.h" // for SIREN_random
#include "SIREN/interactions/Interaction.h" // for Interaction
#include "SIREN/interactions/CrossSection.h"

#include "marley/Generator.hh"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class MarleyCrossSection : public CrossSection {
friend cereal::access;
private:
    //void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;
    std::string marley_config_;
    marley::Generator marley_generator_;
    void InitializeMarley(const std::string& marley_config);
    bool has_nu_cc;
    bool has_nubar_cc;
    bool has_nc;
    bool has_elastic;

public:
    MarleyCrossSection(std::string marley_config, std::string search_path);
    virtual ~MarleyCrossSection() {};
    virtual bool equal(CrossSection const & other) const override;
    virtual double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    virtual double TotalCrossSectionAllFinalStates(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    virtual double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;

    virtual std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    virtual std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    virtual std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("MarleyConfig", marley_config_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("MarleyCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<MarleyCrossSection> & construct, std::uint32_t const version) {
        if(version == 0) {
            std::string marley_config;
            std::string search_path;
            archive(::cereal::make_nvp("MarleyConfig", marley_config));
            construct(marley_config, search_path);
            archive(cereal::virtual_base_class<CrossSection>(construct.ptr()));
        } else {
            throw std::runtime_error("MarleyCrossSection only supports version <= 0!");
        }
    }
}; // class MarleyCrossSection

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::MarleyCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::MarleyCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::MarleyCrossSection);

#endif // SIREN_MarleyCrossSection_H
