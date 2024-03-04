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

std::ostream& operator<<(std::ostream& os, LI::dataclasses::Particle const& p) {
    os << "Particle (" << &p << ")\n";

    std::stringstream ss;
    ss << p.id;
    std::string id_str = ss.str();
    std::string from = "\n";
    std::string to = "\n    ";
    size_t start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }

    os << "ID: " << id_str << "\n";
    os << "Type: " << p.type << "\n";
    os << "Mass: " << p.mass << "\n";
    os << "Momentum: " << p.momentum.at(0) << " " << p.momentum.at(1) << " " << p.momentum.at(2) << " " << p.momentum.at(3) << "\n";
    os << "Position: " << p.position.at(0) << " " << p.position.at(1) << " " << p.position.at(2) << "\n";
    os << "Length: " << p.length << "\n";
    os << "Helicity: " << p.helicity;

    return os;
}

namespace LI {
namespace dataclasses {

Particle::Particle(ParticleID id, ParticleType type, double mass, std::array<double, 4> momentum, std::array<double, 3> position, double length, double helicity) : id(id), type(type), mass(mass), momentum(momentum), position(position), length(length), helicity(helicity) {}

Particle::Particle(ParticleType type, double mass, std::array<double, 4> momentum, std::array<double, 3> position, double length, double helicity) : type(type), mass(mass), momentum(momentum), position(position), length(length), helicity(helicity) {}

ParticleID & Particle::GenerateID() {
    id = ParticleID::GenerateID();
    return id;
}

// Helper functions for dealing with particle types
bool isNeutrino(ParticleType p) {
    return (
           p==ParticleType::NuE      || p==ParticleType::NuEBar ||
           p==ParticleType::NuMu     || p==ParticleType::NuMuBar ||
           p==ParticleType::NuTau    || p==ParticleType::NuTauBar
           );
}

// returns true if a particle is a Lepton. False if not
bool isLepton(ParticleType p){
    return(p==ParticleType::EMinus   || p==ParticleType::EPlus ||
           p==ParticleType::MuMinus  || p==ParticleType::MuPlus ||
           p==ParticleType::TauMinus || p==ParticleType::TauPlus ||
           p==ParticleType::NuE      || p==ParticleType::NuEBar ||
           p==ParticleType::NuMu     || p==ParticleType::NuMuBar ||
           p==ParticleType::NuTau    || p==ParticleType::NuTauBar);
}

// returns true if the particle is either
//        a charged lepton
//   (OR) a "hadrons" particle
bool isCharged(ParticleType p){
    if( !(isLepton(p) || p==ParticleType::Hadrons) ){
        throw std::runtime_error("You should only be using Leptons or Hadrons!");
    }

    // keeps this within scope. Shouldn't be getting some other kind of charged particle
    return(p==ParticleType::EMinus   || p==ParticleType::EPlus ||
           p==ParticleType::MuMinus  || p==ParticleType::MuPlus ||
           p==ParticleType::TauMinus || p==ParticleType::TauPlus ||
           p==ParticleType::Hadrons);
}

} // namespace utilities
} // namespace LI
