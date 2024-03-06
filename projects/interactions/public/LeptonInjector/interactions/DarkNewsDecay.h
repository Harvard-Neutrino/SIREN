#pragma once
#ifndef LI_DarkNewsDecay_H
#define LI_DarkNewsDecay_H

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

#include "LeptonInjector/interactions/Decay.h"   // for Decay
#include "LeptonInjector/dataclasses/Particle.h"  // for Particle

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace LI { namespace dataclasses { struct InteractionSignature; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace interactions {

class DarkNewsDecay : public Decay {
friend cereal::access;
public:
    
    DarkNewsDecay();

    virtual pybind11::object get_representation();
    
    virtual bool equal(Decay const & other) const override;
    
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual double TotalDecayWidth(LI::dataclasses::ParticleType primary) const override;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual void SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<LI::utilities::LI_random>) const;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<LI::utilities::LI_random>) const override;
    
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override = 0; // Requires python-side implementation
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::ParticleType primary) const override = 0; // Requires python-side implementation
    
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
public:
    virtual std::vector<std::string> DensityVariables() const override = 0; // Requires python-side implementation
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("DarkNewsDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<DarkNewsDecay> & construct, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(construct.ptr())));
        } else {
            throw std::runtime_error("DarkNewsDecay only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace LI

CEREAL_CLASS_VERSION(LI::interactions::DarkNewsDecay, 0);
CEREAL_REGISTER_TYPE(LI::interactions::DarkNewsDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::interactions::Decay, LI::interactions::DarkNewsDecay);

#endif // LI_DarkNewsDecay_H
