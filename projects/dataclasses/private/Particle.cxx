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

Particle::Particle(ParticleID id, ParticleType type, double mass, std::array<double, 4> momentum, std::array<double, 3> position, double helicity) : id(id), type(type), mass(mass), momentum(momentum), position(position), helicity(helicity) {}

Particle::Particle(ParticleType type, double mass, std::array<double, 4> momentum, std::array<double, 3> position, double helicity) : type(type), mass(mass), momentum(momentum), position(position), helicity(helicity) {}

ParticleID & Particle::GenerateID() {
    particle_id = ParticleID::GenerateID();
    return particle_id;
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


// returns string of particle's name
std::string particleName(Particle::ParticleType p){
    return(Particle(p).GetTypeString());
}


// gets the mass of a particle for a given type
double particleMass(Particle::ParticleType type){
    Particle p(type);
    if(!p.HasMass()){
        return(0);
    }
    return(p.GetMass());
}

// Uses a particle's type (mass) and total energy to calculate kinetic energy
double kineticEnergy(Particle::ParticleType type, double totalEnergy){
    double mass=particleMass(type);
    if(totalEnergy<mass){
        return(0.);
    }
    return(sqrt(totalEnergy*totalEnergy-mass*mass));
}

// uses the particle type and kinetic energy to calculate the speed of the particle
// relies on the constants!
double particleSpeed(Particle::ParticleType type, double kineticEnergy){
    Particle p=Particle(type);
    if(!p.HasMass()){
        return(LI::utilities::Constants::c);
    }
    double mass=p.GetMass();
    if(kineticEnergy<0){
        return(0.);
    }

    // IF mass>0 THEN return mass/(stuff) ... ELSE return 0
    double r=(mass>0 ? mass/(kineticEnergy+mass) : 0.);
    return(LI::utilities::Constants::c*sqrt(1-r*r));
}

Particle::ParticleShape decideShape(Particle::ParticleType t){
    switch(t){
        case Particle::ParticleType::MuMinus:  case Particle::ParticleType::MuPlus:
        case Particle::ParticleType::TauMinus: case Particle::ParticleType::TauPlus:
        case Particle::ParticleType::NuE:      case Particle::ParticleType::NuEBar:
        case Particle::ParticleType::NuMu:     case Particle::ParticleType::NuMuBar:
        case Particle::ParticleType::NuTau:    case Particle::ParticleType::NuTauBar:
            return(Particle::ParticleShape::MCTrack);
        case Particle::ParticleType::EMinus: case Particle::ParticleType::EPlus:
        case Particle::ParticleType::Hadrons:
            return(Particle::ParticleShape::Cascade);
        case Particle::ParticleType::unknown:
            return(Particle::ParticleShape::unknown);
        default:
            throw "BadShape"; // this replaces the previous fatal log
//				log_fatal_stream("Unable to decide shape for unexpected particle type: " << particleName(t));
    }
}

// This function returns the primary particle type given the final state particles
// returns a particle type object
Particle::ParticleType deduceInitialType(Particle::ParticleType pType1, Particle::ParticleType pType2){
    //only accept certain particle types in general
    if(!isLepton(pType1) && pType1!=Particle::ParticleType::Hadrons)
        throw std::runtime_error("BadParticle");
    if(!isLepton(pType2) && pType2!=Particle::ParticleType::Hadrons)
        throw std::runtime_error("BadParticle");

    bool c1=isCharged(pType1);
    bool c2=isCharged(pType2);
    bool l1=isLepton(pType1);
    bool l2=isLepton(pType2);

    //at least one particle should be charged
    if(!c1 && !c2)
        throw std::runtime_error("Final state should have at least one charged particle");
    //first particle is charged, second is not
    if(c1 && !c2){
        //valid cases are charged lepton + matching antineutrino for GR
        if(l1){
            //!c2 => pType2 is a neutrino
            if(!((pType1==Particle::ParticleType::EMinus   && pType2==Particle::ParticleType::NuEBar) ||
                 (pType1==Particle::ParticleType::EPlus    && pType2==Particle::ParticleType::NuE) ||
                 (pType1==Particle::ParticleType::MuMinus  && pType2==Particle::ParticleType::NuMuBar) ||
                 (pType1==Particle::ParticleType::MuPlus   && pType2==Particle::ParticleType::NuMu) ||
                 (pType1==Particle::ParticleType::TauMinus && pType2==Particle::ParticleType::NuTauBar) ||
                 (pType1==Particle::ParticleType::TauPlus  && pType2==Particle::ParticleType::NuTau)))
                 throw std::runtime_error("Final states with a charged lepton must have an anti-matching neutrino.");
            return(Particle::ParticleType::NuEBar);
        }
        throw std::runtime_error("BadFinal");
    }

    //first particle is neutral, second is charged
    if(!c1 && c2){
        if(l1 && pType2==Particle::ParticleType::Hadrons){
            //particle 1 is a neutral lepton, so it must be a neutrino
            return(pType1); //the incoming neutrino type is the same as the outgoing
        }
        throw std::runtime_error("BadFinal");
    }

    //have two charged particles
    if(c1 && c2){
        //no two charged lepton states
        if(l1 && l2)
            throw std::runtime_error("BadFinal");
        //lepton should be given first
        if(!l1 && l2)
            throw std::runtime_error("BadFinal");
        if(l1 && !l2){ //valid: charged lepton + Hadrons for CC
            switch(pType1){
                case Particle::ParticleType::EMinus: return(Particle::ParticleType::NuE);
                case Particle::ParticleType::EPlus: return(Particle::ParticleType::NuEBar);
                case Particle::ParticleType::MuMinus: return(Particle::ParticleType::NuMu);
                case Particle::ParticleType::MuPlus: return(Particle::ParticleType::NuMuBar);
                case Particle::ParticleType::TauMinus: return(Particle::ParticleType::NuTau);
                case Particle::ParticleType::TauPlus: return(Particle::ParticleType::NuTauBar);
                default: assert(false && "This point should be unreachable");
            }
        }
        if(!l1 && !l2){ //valid: two hadrons (for GR)
            return(Particle::ParticleType::NuEBar);
        }
    }
    throw std::runtime_error("You must be a wizard: this point should be unreachable");
}
} // namespace utilities
} // namespace LI
