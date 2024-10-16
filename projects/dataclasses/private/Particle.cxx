#include "SIREN/dataclasses/Particle.h"

#include <iosfwd>
#include <math.h>
#include <string>
#include <cstdint>
#include <sstream>
#include <stdint.h>
#include <assert.h>
#include <stdexcept>

#include <cereal/cereal.hpp>

#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/StringManipulation.h"

std::ostream& operator<<(std::ostream& os, siren::dataclasses::Particle const& p) {
    os << to_repr(p);
    return os;
}

std::string to_str(siren::dataclasses::Particle const& p) {
    using siren::utilities::tab;
    std::stringstream ss;
    ss << "[ Particle (" << &p << "):\n";
    ss << tab << "ID: " << to_repr(p.id) << '\n';
    ss << tab << "Type: " << p.type << '\n';
    ss << tab << "Mass: " << p.mass << '\n';
    ss << tab << "Momentum: " << p.momentum.at(0) << ' ' << p.momentum.at(1) << ' ' << p.momentum.at(2) << ' ' << p.momentum.at(3) << '\n';
    ss << tab << "Position: " << p.position.at(0) << ' ' << p.position.at(1) << ' ' << p.position.at(2) << '\n';
    ss << tab << "Length: " << p.length << '\n';
    ss << tab << "Helicity: " << p.helicity << '\n';
    ss << ']';

    return ss.str();
}

std::string to_repr(siren::dataclasses::Particle const& p) {
    std::stringstream ss;
    ss << "Particle(";
    ss << "id=" << to_repr(p.id) << ", ";
    ss << "type=" << p.type << ", ";
    ss << "mass=" << p.mass << ", ";
    ss << "momentum=(" << p.momentum.at(0) << ", " << p.momentum.at(1) << ", " << p.momentum.at(2) << ", " << p.momentum.at(3) << "), ";
    ss << "position=(" << p.position.at(0) << ", " << p.position.at(1) << ", " << p.position.at(2) << "), ";
    ss << "length=" << p.length << ", ";
    ss << "helicity=" << p.helicity;
    ss << ')';

    return ss.str();
}

namespace siren {
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
} // namespace siren
