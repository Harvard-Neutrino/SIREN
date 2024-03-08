#pragma once
#ifndef LI_Decay_H
#define LI_Decay_H

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

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace SI { namespace dataclasses { struct InteractionSignature; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace interactions {

class Decay{
friend cereal::access;
private:
public:
    Decay();
    virtual ~Decay() {};
    bool operator==(Decay const & other) const;
    virtual bool equal(Decay const & other) const = 0;
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayWidth(SI::dataclasses::ParticleType primary) const = 0;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayLength(SI::dataclasses::InteractionRecord const & record) const;
    virtual double TotalDecayLengthForFinalState(SI::dataclasses::InteractionRecord const & record) const;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<SI::utilities::LI_random>) const = 0;
    virtual std::vector<SI::dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;
    virtual std::vector<SI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(SI::dataclasses::ParticleType primary) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};

}; // class Decay

} // namespace interactions
} // namespace SI

CEREAL_CLASS_VERSION(SI::interactions::Decay, 0);

#endif // LI_Decay_H
