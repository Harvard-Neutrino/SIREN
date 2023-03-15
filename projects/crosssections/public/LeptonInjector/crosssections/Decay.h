#pragma once
#ifndef LI_Decay_H
#define LI_Decay_H

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

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
    virtual double TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const = 0;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const = 0;
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;

}; // class Decay

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::crosssections::Decay, 0);

#endif // LI_Decay_H
