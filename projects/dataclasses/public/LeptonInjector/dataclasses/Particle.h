#pragma once
#ifndef LI_Particle_H
#define LI_Particle_H

// Used to define the Particle class
// Partiles have a type, energy, position, and direction

#include <set>
#include <string>
#include <utility>
#include <stdint.h>

#include "LeptonInjector/dataclasses/ParticleType.h"

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include <cereal/types/set.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>

namespace LI {
namespace dataclasses {

// simple data structure for particles
class Particle {
public:
    typedef LI::dataclasses::ParticleType ParticleType;

    Particle() = default;

    ParticleID id;
    ParticleType type = ParticleType::unknown;
    double mass = 0;
    std::array<double, 4> momentum = {0, 0, 0, 0};
    std::array<double, 3> position = {0, 0, 0};
    double helicity = 0;

};

// prototype some of the particle helper functions

bool isLepton(Particle::ParticleType p);
bool isCharged(Particle::ParticleType p);
bool isNeutrino(Particle::ParticleType p);
std::string particleName( Particle::ParticleType p);
double particleMass( Particle::ParticleType type);
double kineticEnergy( Particle::ParticleType type, double totalEnergy);
double particleSpeed( Particle::ParticleType type, double kineticEnergy);
Particle::ParticleShape decideShape(Particle::ParticleType t);
Particle::ParticleType deduceInitialType( Particle::ParticleType pType1, Particle::ParticleType pType2);

} // namespace utilities
} // namespace LI

#endif // LI_Particle_H
