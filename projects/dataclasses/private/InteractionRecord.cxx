#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include <cmath>
#include <tuple>    // for tie, operator==, tuple
#include <ostream>  // for operator<<, basic_ostream, char_traits, endl, ost...

namespace LI {
namespace dataclasses {

PrimaryDistributionRecord::PrimaryDistributionRecord(InteractionRecord const & record) :
    record(record),
    signature(record.signature),
    id((record.primary_id) ? (record.primary_id) : (ParticleID::GenerateID())),
    type(record.signature.primary_type)
{}

InteractionSignature const & PrimaryDistributionRecord::GetSignature() const {
    return signature;
}

Particle PrimaryDistributionRecord::GetParticle() {
    Particle p;
    p.id = id;
    p.type = signature.primary_type;
    try {
        p.mass = GetMass();
    } catch(...) {
        p.mass = 0;
    }
    try {
        p.momentum = GetFourMomentum();
    } catch(...) {
        p.momentum = {0, 0, 0, 0};
    }
    try {
        p.position = GetInitialPosition();
    } catch(...) {
        p.position = {0, 0, 0};
    }
    try {
        p.length = GetLength();
    } catch(...) {
        p.length = 0;
    }
    try {
        p.helicity = GetHelicity();
    } catch(...) {
        p.helicity = 0;
    }
    return p;
}

ParticleID const & PrimaryDistributionRecord::GetID() const {
    return id;
}

ParticleType const & PrimaryDistributionRecord::GetType() const {
    return type;
}

double PrimaryDistributionRecord::GetMass() const {
    if(mass_set) {
        return mass;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateMass();
        return non_const_this.mass;
    }
}

double const & PrimaryDistributionRecord::GetMass() {
    if(not mass_set) {
        UpdateMass();
    }
    return mass;
}

double PrimaryDistributionRecord::GetEnergy() const {
    if(energy_set) {
        return energy;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateEnergy();
        return non_const_this.energy;
    }
}

double const & PrimaryDistributionRecord::GetEnergy() {
    if(not energy_set) {
        UpdateEnergy();
    }
    return energy;
}

double PrimaryDistributionRecord::GetKineticEnergy() const {
    if(kinetic_energy_set) {
        return kinetic_energy;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateKineticEnergy();
        return non_const_this.kinetic_energy;
    }
}

double const & PrimaryDistributionRecord::GetKineticEnergy() {
    if(not kinetic_energy_set) {
        UpdateKineticEnergy();
    }
    return kinetic_energy;
}


std::array<double, 3> PrimaryDistributionRecord::GetDirection() const {
    if(direction_set) {
        return direction;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateDirection();
        return non_const_this.direction;
    }
}

std::array<double, 3> const & PrimaryDistributionRecord::GetDirection() {
    if(not direction_set) {
        UpdateDirection();
    }
    return direction;
}

std::array<double, 3> PrimaryDistributionRecord::GetThreeMomentum() const {
    if(momentum_set) {
        return momentum;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateMomentum();
        return non_const_this.momentum;
    }
}

std::array<double, 3> const & PrimaryDistributionRecord::GetThreeMomentum() {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return momentum;
}

std::array<double, 4> PrimaryDistributionRecord::GetFourMomentum() const {
    if(momentum_set) {
        return {energy, momentum.at(0), momentum.at(1), momentum.at(2)};
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateMomentum();
        return {non_const_this.energy, non_const_this.momentum.at(0), non_const_this.momentum.at(1), non_const_this.momentum.at(2)};
    }
}

std::array<double, 4> PrimaryDistributionRecord::GetFourMomentum() {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return {momentum.at(0), momentum.at(1), momentum.at(2), GetEnergy()};
}

double PrimaryDistributionRecord::GetLength() const {
    if(length_set) {
        return length;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateLength();
        return non_const_this.length;
    }
}

double const & PrimaryDistributionRecord::GetLength() {
    if(not length_set) {
        UpdateLength();
    }
    return length;
}

std::array<double, 3> PrimaryDistributionRecord::GetInitialPosition() const {
    if(initial_position_set) {
        return initial_position;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateInitialPosition();
        return non_const_this.initial_position;
    }
}

std::array<double, 3> const & PrimaryDistributionRecord::GetInitialPosition() {
    if(not initial_position_set) {
        UpdateInitialPosition();
    }
    return initial_position;
}

std::array<double, 3> PrimaryDistributionRecord::GetInteractionVertex() const {
    if(interaction_vertex_set) {
        return interaction_vertex;
    } else {
        PrimaryDistributionRecord non_const_this = *this;
        non_const_this.UpdateInteractionVertex();
        return non_const_this.interaction_vertex;
    }
}

std::array<double, 3> const & PrimaryDistributionRecord::GetInteractionVertex() {
    if(not interaction_vertex_set) {
        UpdateInteractionVertex();
    }
    return interaction_vertex;
}

double const & PrimaryDistributionRecord::GetHelicity() const {
    return helicity;
}

void PrimaryDistributionRecord::SetParticle(Particle const & particle) {
    if(particle.id != id) {
        throw std::runtime_error("Cannot set particle with different ID!");
    }
    if(particle.type != type) {
        throw std::runtime_error("Cannot set particle with different type!");
    }

    mass_set = true;
    mass = particle.mass;

    momentum_set = true;
    momentum = {particle.momentum.at(1), particle.momentum.at(2), particle.momentum.at(3)};

    energy_set = true;
    energy = particle.momentum.at(0);

    initial_position_set = true;
    initial_position = particle.position;

    length_set = true;
    length = particle.length;

    helicity_set = true;
    helicity = particle.helicity;
}

void PrimaryDistributionRecord::SetMass(double mass) {
    mass_set = true;
    this->mass = mass;
}

void PrimaryDistributionRecord::SetEnergy(double energy) {
    energy_set = true;
    this->energy = energy;
}

void PrimaryDistributionRecord::SetKineticEnergy(double kinetic_energy) {
    kinetic_energy_set = true;
    this->kinetic_energy = kinetic_energy;
}

void PrimaryDistributionRecord::SetDirection(std::array<double, 3> direction) {
    direction_set = true;
    this->direction = direction;
}

void PrimaryDistributionRecord::SetThreeMomentum(std::array<double, 3> momentum) {
    momentum_set = true;
    this->momentum = momentum;
}

void PrimaryDistributionRecord::SetFourMomentum(std::array<double, 4> momentum) {
    momentum_set = true;
    this->momentum = {momentum.at(1), momentum.at(2), momentum.at(3)};
    energy_set = true;
    this->energy = momentum.at(0);
}

void PrimaryDistributionRecord::SetLength(double length) {
    length_set = true;
    this->length = length;
}

void PrimaryDistributionRecord::SetInitialPosition(std::array<double, 3> initial_position) {
    initial_position_set = true;
    this->initial_position = initial_position;
}

void PrimaryDistributionRecord::SetInteractionVertex(std::array<double, 3> interaction_vertex) {
    interaction_vertex_set = true;
    this->interaction_vertex = interaction_vertex;
}

void PrimaryDistributionRecord::SetHelicity(double helicity) {
    helicity_set = true;
    this->helicity = helicity;
}

void PrimaryDistributionRecord::UpdateMass() {
    if(mass_set)
        return;
    if(energy_set and momentum_set) {
        mass = std::sqrt(energy*energy - momentum.at(0)*momentum.at(0) - momentum.at(1)*momentum.at(1) - momentum.at(2)*momentum.at(2));
    } else if(energy_set and kinetic_energy_set) {
        mass = std::sqrt(energy*energy - kinetic_energy*kinetic_energy);
    } else {
        throw std::runtime_error("Cannot calculate mass without energy and momentum or energy and kinetic energy!");
    }
}

void PrimaryDistributionRecord::UpdateEnergy() {
    if(energy_set)
        return;
    if(mass_set and momentum_set) {
        energy = std::sqrt(mass*mass + momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
    } else if(mass_set and kinetic_energy_set) {
        energy = std::sqrt(mass*mass + kinetic_energy*kinetic_energy);
    } else {
        throw std::runtime_error("Cannot calculate energy without mass and momentum or mass and kinetic energy!");
    }
}

void PrimaryDistributionRecord::UpdateKineticEnergy() {
    if(kinetic_energy_set)
        return;
    if(mass_set and energy_set) {
        kinetic_energy = std::sqrt(energy*energy - mass*mass);
    } else if(momentum_set) {
        kinetic_energy = std::sqrt(momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
    } else {
        throw std::runtime_error("Cannot calculate kinetic energy without mass and energy or momentum!");
    }
}

void PrimaryDistributionRecord::UpdateDirection() {
    if(direction_set)
        return;
    if(momentum_set) {
        double magnitude = std::sqrt(momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
        direction = {momentum.at(0)/magnitude, momentum.at(1)/magnitude, momentum.at(2)/magnitude};
    } else if(initial_position_set and interaction_vertex_set) {
        direction = {interaction_vertex.at(0) - initial_position.at(0), interaction_vertex.at(1) - initial_position.at(1), interaction_vertex.at(2) - initial_position.at(2)};
        double magnitude = std::sqrt(direction.at(0)*direction.at(0) + direction.at(1)*direction.at(1) + direction.at(2)*direction.at(2));
        direction = {direction.at(0)/magnitude, direction.at(1)/magnitude, direction.at(2)/magnitude};
    } else {
        throw std::runtime_error("Cannot calculate direction without momentum or initial position and interaction vertex!");
    }
}

void PrimaryDistributionRecord::UpdateMomentum() {
    if(momentum_set)
        return;
    if(energy_set and mass_set and direction_set) {
        double magnitude = std::sqrt(energy*energy - mass*mass);
        momentum = {magnitude*direction.at(0), magnitude*direction.at(1), magnitude*direction.at(2)};
    } else if(kinetic_energy_set and direction_set) {
        double magnitude = kinetic_energy;
        momentum = {magnitude*direction.at(0), magnitude*direction.at(1), magnitude*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate momentum without energy and mass and direction or kinetic energy and direction!");
    }
}

void PrimaryDistributionRecord::UpdateLength() {
    if(length_set)
        return;
    if(initial_position_set and interaction_vertex_set) {
        length = std::sqrt(
            (interaction_vertex.at(0) - initial_position.at(0))*(interaction_vertex.at(0) - initial_position.at(0)) +
            (interaction_vertex.at(1) - initial_position.at(1))*(interaction_vertex.at(1) - initial_position.at(1)) +
            (interaction_vertex.at(2) - initial_position.at(2))*(interaction_vertex.at(2) - initial_position.at(2))
        );
    } else {
        throw std::runtime_error("Cannot calculate length without initial position and interaction vertex!");
    }
}

void PrimaryDistributionRecord::UpdateInitialPosition() {
    if(initial_position_set)
        return;
    if(interaction_vertex_set and direction_set and length_set) {
        initial_position = {interaction_vertex.at(0) - length*direction.at(0), interaction_vertex.at(1) - length*direction.at(1), interaction_vertex.at(2) - length*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate initial position without interaction vertex and direction and length!");
    }
}

void PrimaryDistributionRecord::UpdateInteractionVertex() {
    if(interaction_vertex_set)
        return;
    if(initial_position_set and direction_set and length_set) {
        interaction_vertex = {initial_position.at(0) + length*direction.at(0), initial_position.at(1) + length*direction.at(1), initial_position.at(2) + length*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate interaction vertex without initial position and direction and length!");
    }
}

void PrimaryDistributionRecord::Finalize(InteractionRecord & record) {
    record.signature.primary_type = signature.primary_type;
    record.primary_id = GetID();
    record.primary_initial_position = GetInitialPosition();
    record.primary_mass = GetMass();
    record.primary_momentum = GetFourMomentum();
    record.primary_helicity = GetHelicity();
}

/////////////////////////////////////////

SecondaryParticleRecord::SecondaryParticleRecord(InteractionRecord const & record, size_t secondary_index) :
    id((record.secondary_ids.size() > secondary_index and record.secondary_ids.at(secondary_index)) ? (record.secondary_ids.at(secondary_index)) : (ParticleID::GenerateID())),
    type(record.signature.secondary_types.at(secondary_index)),
    initial_position(record.interaction_vertex)
{}

Particle SecondaryParticleRecord::GetParticle() {
    Particle p;
    p.id = id;
    p.type = type;
    try {
        p.mass = GetMass();
    } catch(...) {
        p.mass = 0;
    }
    try {
        p.momentum = GetFourMomentum();
    } catch(...) {
        p.momentum = {0, 0, 0, 0};
    }
    try {
        p.position = GetInitialPosition();
    } catch(...) {
        p.position = {0, 0, 0};
    }
    try {
        p.helicity = GetHelicity();
    } catch(...) {
        p.helicity = 0;
    }
    return p;
}

void SecondaryParticleRecord::SetParticle(Particle const & particle) {
    if(particle.id != id) {
        throw std::runtime_error("Cannot set particle with different ID!");
    }
    if(particle.type != type) {
        throw std::runtime_error("Cannot set particle with different type!");
    }

    mass_set = true;
    mass = particle.mass;

    momentum_set = true;
    momentum = {particle.momentum.at(1), particle.momentum.at(2), particle.momentum.at(3)};

    energy_set = true;
    energy = particle.momentum.at(0);

    helicity_set = true;
    helicity = particle.helicity;
}

ParticleID const & SecondaryParticleRecord::GetID() const {
    return id;
}

ParticleType const & SecondaryParticleRecord::GetType() const {
    return type;
}

double SecondaryParticleRecord::GetMass() const {
    if(mass_set) {
        return mass;
    } else {
        SecondaryParticleRecord non_const_this = *this;
        non_const_this.UpdateMass();
        return non_const_this.mass;
    }
}

double const & SecondaryParticleRecord::GetMass() {
    if(not mass_set) {
        UpdateMass();
    }
    return mass;
}

double SecondaryParticleRecord::GetEnergy() const {
    if(energy_set) {
        return energy;
    } else {
        SecondaryParticleRecord non_const_this = *this;
        non_const_this.UpdateEnergy();
        return non_const_this.energy;
    }
}

double const & SecondaryParticleRecord::GetEnergy() {
    if(not energy_set) {
        UpdateEnergy();
    }
    return energy;
}

double SecondaryParticleRecord::GetKineticEnergy() const {
    if(kinetic_energy_set) {
        return kinetic_energy;
    } else {
        SecondaryParticleRecord non_const_this = *this;
        non_const_this.UpdateKineticEnergy();
        return non_const_this.kinetic_energy;
    }
}

double const & SecondaryParticleRecord::GetKineticEnergy() {
    if(not kinetic_energy_set) {
        UpdateKineticEnergy();
    }
    return kinetic_energy;
}

std::array<double, 3> SecondaryParticleRecord::GetDirection() const {
    if(direction_set) {
        return direction;
    } else {
        SecondaryParticleRecord non_const_this = *this;
        non_const_this.UpdateDirection();
        return non_const_this.direction;
    }
}

std::array<double, 3> const & SecondaryParticleRecord::GetDirection() {
    if(not direction_set) {
        UpdateDirection();
    }
    return direction;
}

std::array<double, 3> SecondaryParticleRecord::GetThreeMomentum() const {
    if(momentum_set) {
        return momentum;
    } else {
        SecondaryParticleRecord non_const_this = *this;
        non_const_this.UpdateMomentum();
        return non_const_this.momentum;
    }
}

std::array<double, 3> const & SecondaryParticleRecord::GetThreeMomentum() {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return momentum;
}

std::array<double, 4> SecondaryParticleRecord::GetFourMomentum() const {
    if(momentum_set) {
        return {energy, momentum.at(0), momentum.at(1), momentum.at(2)};
    } else {
        SecondaryParticleRecord non_const_this = *this;
        non_const_this.UpdateMomentum();
        return {non_const_this.energy, non_const_this.momentum.at(0), non_const_this.momentum.at(1), non_const_this.momentum.at(2)};
    }
}

std::array<double, 4> SecondaryParticleRecord::GetFourMomentum() {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return {momentum.at(0), momentum.at(1), momentum.at(2), GetEnergy()};
}

std::array<double, 3> const & SecondaryParticleRecord::GetInitialPosition() const {
    return initial_position;
}

double const & SecondaryParticleRecord::GetHelicity() const {
    return helicity;
}

void SecondaryParticleRecord::SetMass(double mass) {
    mass_set = true;
    this->mass = mass;
}

void SecondaryParticleRecord::SetEnergy(double energy) {
    energy_set = true;
    this->energy = energy;
}

void SecondaryParticleRecord::SetKineticEnergy(double kinetic_energy) {
    kinetic_energy_set = true;
    this->kinetic_energy = kinetic_energy;
}

void SecondaryParticleRecord::SetDirection(std::array<double, 3> direction) {
    direction_set = true;
    this->direction = direction;
}

void SecondaryParticleRecord::SetThreeMomentum(std::array<double, 3> momentum) {
    momentum_set = true;
    this->momentum = momentum;
}

void SecondaryParticleRecord::SetFourMomentum(std::array<double, 4> momentum) {
    momentum_set = true;
    this->momentum = {momentum.at(1), momentum.at(2), momentum.at(3)};
    energy_set = true;
    this->energy = momentum.at(0);
}

void SecondaryParticleRecord::SetHelicity(double helicity) {
    helicity_set = true;
    this->helicity = helicity;
}

void SecondaryParticleRecord::UpdateMass() {
    if(mass_set)
        return;
    if(energy_set and momentum_set) {
        mass = std::sqrt(energy*energy - momentum.at(0)*momentum.at(0) - momentum.at(1)*momentum.at(1) - momentum.at(2)*momentum.at(2));
    } else if(energy_set and kinetic_energy_set) {
        mass = std::sqrt(energy*energy - kinetic_energy*kinetic_energy);
    } else {
        throw std::runtime_error("Cannot calculate mass without energy and momentum or energy and kinetic energy!");
    }
}

void SecondaryParticleRecord::UpdateEnergy() {
    if(energy_set)
        return;
    if(mass_set and momentum_set) {
        energy = std::sqrt(mass*mass + momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
    } else if(mass_set and kinetic_energy_set) {
        energy = std::sqrt(mass*mass + kinetic_energy*kinetic_energy);
    } else {
        throw std::runtime_error("Cannot calculate energy without mass and momentum or mass and kinetic energy!");
    }
}

void SecondaryParticleRecord::UpdateKineticEnergy() {
    if(kinetic_energy_set)
        return;
    if(mass_set and energy_set) {
        kinetic_energy = std::sqrt(energy*energy - mass*mass);
    } else if(momentum_set) {
        kinetic_energy = std::sqrt(momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
    } else {
        throw std::runtime_error("Cannot calculate kinetic energy without mass and energy or momentum!");
    }
}

void SecondaryParticleRecord::UpdateDirection() {
    if(direction_set)
        return;
    if(momentum_set) {
        double magnitude = std::sqrt(momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
        direction = {momentum.at(0)/magnitude, momentum.at(1)/magnitude, momentum.at(2)/magnitude};
    } else {
        throw std::runtime_error("Cannot calculate direction without momentum or initial position and interaction vertex!");
    }
}

void SecondaryParticleRecord::UpdateMomentum() {
    if(momentum_set)
        return;
    if(energy_set and mass_set and direction_set) {
        double magnitude = std::sqrt(energy*energy - mass*mass);
        momentum = {magnitude*direction.at(0), magnitude*direction.at(1), magnitude*direction.at(2)};
    } else if(kinetic_energy_set and direction_set) {
        double magnitude = kinetic_energy;
        momentum = {magnitude*direction.at(0), magnitude*direction.at(1), magnitude*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate momentum without energy and mass and direction or kinetic energy and direction!");
    }
}

void SecondaryParticleRecord::Finalize(InteractionRecord & record) {
    record.primary_id = id;
    record.signature.primary_type = type;
    record.primary_initial_position = initial_position;
    record.primary_mass = GetMass();
    record.primary_momentum = GetFourMomentum();
    record.primary_helicity = GetHelicity();
}

/////////////////////////////////////////

CrossSectionDistributionRecord::CrossSectionDistributionRecord(InteractionRecord const & record) :
    record(record),
    signature(record.signature),
    primary_id(record.primary_id),
    primary_type(record.signature.primary_type),
    primary_initial_position(record.primary_initial_position),
    primary_mass(record.primary_mass),
    primary_momentum(record.primary_momentum),
    primary_helicity(record.primary_helicity),
    interaction_vertex(record.interaction_vertex),
    target_id((record.target_id) ? (record.target_id) : (ParticleID::GenerateID())),
    target_type(record.signature.target_type),
    target_mass(record.target_mass),
    target_helicity(record.target_helicity) {

    secondary_particles.reserve(record.signature.secondary_types.size());
    for(size_t i = 0; i < record.signature.secondary_types.size(); ++i) {
        secondary_particles.emplace_back(record, i);
    }
}

InteractionSignature const & CrossSectionDistributionRecord::GetSignature() const {
    return signature;
}

ParticleID const & CrossSectionDistributionRecord::GetPrimaryID() const {
    return primary_id;
}

ParticleType const & CrossSectionDistributionRecord::GetPrimaryType() const {
    return primary_type;
}

std::array<double, 3> const & CrossSectionDistributionRecord::GetPrimaryInitialPosition() const {
    return primary_initial_position;
}

double const & CrossSectionDistributionRecord::GetPrimaryMass() const {
    return primary_mass;
}

std::array<double, 4> const & CrossSectionDistributionRecord::GetPrimaryMomentum() const {
    return primary_momentum;
}

double const & CrossSectionDistributionRecord::GetPrimaryHelicity() const {
    return primary_helicity;
}

Particle CrossSectionDistributionRecord::GetPrimaryParticle() const {
    Particle p;
    p.id = primary_id;
    p.type = primary_type;
    p.mass = primary_mass;
    p.momentum = primary_momentum;
    p.position = primary_initial_position;
    p.length = std::sqrt((interaction_vertex.at(0) - primary_initial_position.at(0))*(interaction_vertex.at(0) - primary_initial_position.at(0)) +
                         (interaction_vertex.at(1) - primary_initial_position.at(1))*(interaction_vertex.at(1) - primary_initial_position.at(1)) +
                         (interaction_vertex.at(2) - primary_initial_position.at(2))*(interaction_vertex.at(2) - primary_initial_position.at(2)));
    p.helicity = primary_helicity;
    return p;
}

std::array<double, 3> const & CrossSectionDistributionRecord::GetInteractionVertex() const {
    return interaction_vertex;
}

ParticleID const & CrossSectionDistributionRecord::GetTargetID() const {
    return target_id;
}

ParticleType const & CrossSectionDistributionRecord::GetTargetType() const {
    return target_type;
}

double const & CrossSectionDistributionRecord::GetTargetMass() const {
    return target_mass;
}

double const & CrossSectionDistributionRecord::GetTargetHelicity() const {
    return target_helicity;
}

double & CrossSectionDistributionRecord::GetTargetMass() {
    return target_mass;
}

double & CrossSectionDistributionRecord::GetTargetHelicity() {
    return target_helicity;
}

std::map<std::string, double> & CrossSectionDistributionRecord::GetInteractionParameters() {
    return interaction_parameters;
}

void CrossSectionDistributionRecord::SetTargetMass(double mass) {
    target_mass = mass;
}

void CrossSectionDistributionRecord::SetTargetHelicity(double helicity) {
    target_helicity = helicity;
}

void CrossSectionDistributionRecord::SetInteractionParameters(std::map<std::string, double> const & parameters) {
    interaction_parameters = parameters;
}

void CrossSectionDistributionRecord::SetInteractionParameter(std::string const & name, double value) {
    interaction_parameters[name] = value;
}

Particle CrossSectionDistributionRecord::GetTargetParticle() const {
    Particle p;
    p.id = target_id;
    p.type = target_type;
    p.mass = target_mass;
    p.helicity = target_helicity;
    return p;
}

void CrossSectionDistributionRecord::SetTargetParticle(Particle const & particle) {
    if(particle.id != target_id) {
        throw std::runtime_error("Cannot set particle with different ID!");
    }
    if(particle.type != target_type) {
        throw std::runtime_error("Cannot set particle with different type!");
    }

    target_mass = particle.mass;
    target_helicity = particle.helicity;
}

SecondaryParticleRecord & CrossSectionDistributionRecord::GetSecondaryParticleRecord(size_t index) {
    return secondary_particles.at(index);
}

Particle CrossSectionDistributionRecord::GetSecondaryParticle(size_t index) {
    return secondary_particles.at(index).GetParticle();
}

std::vector<Particle> CrossSectionDistributionRecord::GetSecondaryParticles() {
    std::vector<Particle> particles;
    for(SecondaryParticleRecord & secondary: secondary_particles) {
        particles.push_back(secondary.GetParticle());
    }
    return particles;
}

void CrossSectionDistributionRecord::SetSecondaryParticle(size_t index, Particle const & particle) {
    secondary_particles.at(index).SetParticle(particle);
}

void CrossSectionDistributionRecord::SetSecondaryParticles(std::vector<Particle> const & particles) {
    if(particles.size() != secondary_particles.size()) {
        throw std::runtime_error("Cannot set particles with different size!");
    }
    for(size_t i = 0; i < particles.size(); ++i) {
        secondary_particles.at(i).SetParticle(particles.at(i));
    }
}

void CrossSectionDistributionRecord::Finalize(InteractionRecord & record) {
    record.target_id = target_id;
    record.target_mass = target_mass;
    record.target_helicity = target_helicity;

    record.interaction_parameters = interaction_parameters;

    record.secondary_ids.clear();
    record.secondary_masses.clear();
    record.secondary_momenta.clear();
    record.secondary_helicities.clear();

    for(SecondaryParticleRecord const & secondary: secondary_particles) {
        record.secondary_ids.push_back(secondary.GetID());
        record.secondary_masses.push_back(secondary.GetMass());
        record.secondary_momenta.push_back(secondary.GetFourMomentum());
        record.secondary_helicities.push_back(secondary.GetHelicity());
    }
}

/////////////////////////////////////////

SecondaryDistributionRecord::SecondaryDistributionRecord(InteractionRecord const & record, size_t secondary_index) :
    parent_record(record),
    secondary_index(secondary_index),
    id((record.secondary_ids.at(secondary_index)) ? (record.secondary_ids.at(secondary_index)) : (ParticleID::GenerateID())),
    type(record.signature.secondary_types.at(secondary_index)),
    mass(record.secondary_masses.at(secondary_index)),
    direction(),
    momentum(record.secondary_momenta.at(secondary_index)),
    helicity(record.secondary_helicities.at(secondary_index)),
    initial_position(record.primary_initial_position) {
    double magnitude = std::sqrt(momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2) + momentum.at(3)*momentum.at(3));
    direction = {momentum.at(1)/magnitude, momentum.at(2)/magnitude, momentum.at(3)/magnitude};
}

double & SecondaryDistributionRecord::GetLength() {
    return length;
}

void SecondaryDistributionRecord::SetLength(double const & length) {
    this->length = length;
}

void SecondaryDistributionRecord::Finalize(InteractionRecord & record) {
    record.primary_id = id;
    record.primary_initial_position = initial_position;
    record.primary_mass = mass;
    record.primary_momentum = momentum;
    record.primary_helicity = helicity;
    record.interaction_vertex = initial_position;
    record.interaction_vertex.at(0) += length*direction.at(0);
    record.interaction_vertex.at(1) += length*direction.at(1);
    record.interaction_vertex.at(2) += length*direction.at(2);
}

/////////////////////////////////////////

bool InteractionRecord::operator==(InteractionRecord const & other) const {
    return std::tie(
        signature,
        primary_id,
        primary_initial_position,
        primary_mass,
        primary_momentum,
        primary_helicity,
        target_id,
        target_mass,
        target_helicity,
        interaction_vertex,
        secondary_ids,
        secondary_masses,
        secondary_momenta,
        secondary_helicities,
        interaction_parameters)
        ==
        std::tie(
        other.signature,
        other.primary_id,
        other.primary_initial_position,
        other.primary_mass,
        other.primary_momentum,
        other.primary_helicity,
        other.target_id,
        other.target_mass,
        other.target_helicity,
        other.interaction_vertex,
        other.secondary_ids,
        other.secondary_masses,
        other.secondary_momenta,
        other.secondary_helicities,
        other.interaction_parameters);
}

bool InteractionRecord::operator<(InteractionRecord const & other) const {
    return std::tie(
        signature,
        primary_id,
        primary_initial_position,
        primary_mass,
        primary_momentum,
        primary_helicity,
        target_id,
        target_mass,
        target_helicity,
        interaction_vertex,
        secondary_ids,
        secondary_masses,
        secondary_momenta,
        secondary_helicities,
        interaction_parameters)
        <
        std::tie(
        other.signature,
        other.primary_id,
        other.primary_initial_position,
        other.primary_mass,
        other.primary_momentum,
        other.primary_helicity,
        other.target_id,
        other.target_mass,
        other.target_helicity,
        other.interaction_vertex,
        other.secondary_ids,
        other.secondary_masses,
        other.secondary_momenta,
        other.secondary_helicities,
        other.interaction_parameters);
}

} // namespace dataclasses
} // namespace LI

std::ostream& operator<<(std::ostream& os, LI::dataclasses::InteractionRecord const& record) {
    std::stringstream ss;
    ss << "InteractionRecord (" << &record << ") ";
    os << ss.str() << '\n';
    os << "Signature(" << &record.signature << "): " << record.signature.primary_type << " + " << record.signature.target_type << " ->";
    for(auto secondary: record.signature.secondary_types) {
        os << " " << secondary;
    }
    os << "\n";

    os << "PrimaryID: " << record.primary_id << "\n";
    os << "PrimaryInitialPosition: " << record.primary_initial_position.at(0) << " " << record.primary_initial_position.at(1) << " " << record.primary_initial_position.at(2) << "\n";
    os << "InteractionVertex: " << record.interaction_vertex.at(0) << " " << record.interaction_vertex.at(1) << " " << record.interaction_vertex.at(2) << "\n";
    os << "PrimaryMass: " << record.primary_mass << "\n";
    os << "PrimaryMomentum: " << record.primary_momentum.at(0) << " " << record.primary_momentum.at(1) << " " << record.primary_momentum.at(2) << " " << record.primary_momentum.at(3) << "\n";
    os << "TargetID: " << record.target_id << "\n";
    os << "TargetMass: " << record.target_mass << "\n";
    os << "SecondaryIDs:\n";
    for(auto const & secondary: record.secondary_ids) {
        os << "\t" << secondary << "\n";
    }
    os << "SecondaryMomenta:\n";
    for(auto const & secondary: record.secondary_momenta) {
        os << "\t" << secondary.at(0) << " " << secondary.at(1) << " " << secondary.at(2) << " " << secondary.at(3) << "\n";
    }
    os << "SecondaryMasses:\n";
    for(auto const & secondary: record.secondary_masses) {
        os << "\t" << secondary << "\n";
    }
    os << "InteractionParameters:\n";
    for(std::pair<std::string const, double> const & param : record.interaction_parameters) {
        os << "\t\"" << param.first << "\": " << param.second << "\n";
    }
    os << std::endl;

    return os;
}

