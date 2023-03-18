#pragma once
#ifndef LI_Decay_H
#define LI_Decay_H

#include <memory>
#include <string>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

namespace LI {
namespace utilities {
class LI_random;
}
}

namespace LI {
namespace crosssections {

class Decay{
  friend cereal::access;
  private:
  public: 
    Decay();
    virtual ~Decay() {};
    bool operator==(Decay const & other) const;
    virtual bool equal(Decay const & other) const = 0;
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const = 0;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayLength(LI::dataclasses::InteractionRecord const & record) const;
    virtual double TotalDecayLengthForFinalState(LI::dataclasses::InteractionRecord const & record) const;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const = 0;
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};

}; // class Decay

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::crosssections::Decay, 0);

#endif // LI_Decay_H
