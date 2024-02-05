#include "LeptonInjector/dataclasses/Particle.h"

#include <iosfwd>
#include <math.h>
#include <string>
#include <cstdint>
#include <stdint.h>
#include <assert.h>
#include <stdexcept>

#include <cereal/cereal.hpp>

#include "LeptonInjector/utilities/Constants.h"

namespace LI {
namespace dataclasses {

Particle::Particle(ParticleID id, ParticleType type, double mass, std::array<double, 4> momentum, std::array<double, 3> position, double length, double helicity) : id(id), type(type), mass(mass), momentum(momentum), position(position), length(length), helicity(helicity) {}

Particle::Particle(ParticleType type, double mass, std::array<double, 4> momentum, std::array<double, 3> position, double length, double helicity) : type(type), mass(mass), momentum(momentum), position(position), length(length), helicity(helicity) {}

ParticleID & Particle::GenerateID() {
    id = ParticleID::GenerateID();
    return id;
}

// Helper functions for dealing with particle types
bool isNeutrino(Particle::ParticleType p) {
    return (
           p==Particle::ParticleType::NuE      || p==Particle::ParticleType::NuEBar ||
           p==Particle::ParticleType::NuMu     || p==Particle::ParticleType::NuMuBar ||
           p==Particle::ParticleType::NuTau    || p==Particle::ParticleType::NuTauBar
           );
}

// returns true if a particle is a Lepton. False if not
bool isLepton(Particle::ParticleType p){
    return(p==Particle::ParticleType::EMinus   || p==Particle::ParticleType::EPlus ||
           p==Particle::ParticleType::MuMinus  || p==Particle::ParticleType::MuPlus ||
           p==Particle::ParticleType::TauMinus || p==Particle::ParticleType::TauPlus ||
           p==Particle::ParticleType::NuE      || p==Particle::ParticleType::NuEBar ||
           p==Particle::ParticleType::NuMu     || p==Particle::ParticleType::NuMuBar ||
           p==Particle::ParticleType::NuTau    || p==Particle::ParticleType::NuTauBar);
}

// returns true if the particle is either
//        a charged lepton
//   (OR) a "hadrons" particle
bool isCharged(Particle::ParticleType p){
    if( !(isLepton(p) || p==Particle::ParticleType::Hadrons) ){
        throw std::runtime_error("You should only be using Leptons or Hadrons!");
    }

    // keeps this within scope. Shouldn't be getting some other kind of charged particle
    return(p==Particle::ParticleType::EMinus   || p==Particle::ParticleType::EPlus ||
           p==Particle::ParticleType::MuMinus  || p==Particle::ParticleType::MuPlus ||
           p==Particle::ParticleType::TauMinus || p==Particle::ParticleType::TauPlus ||
           p==Particle::ParticleType::Hadrons);
}

} // namespace utilities
} // namespace LI
