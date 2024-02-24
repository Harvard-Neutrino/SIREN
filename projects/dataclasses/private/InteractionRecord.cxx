#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include <cmath>
#include <tuple>    // for tie, operator==, tuple
#include <cassert>
#include <ostream>  // for operator<<, basic_ostream, char_traits, endl, ost...

std::ostream& operator<<(std::ostream& os, LI::dataclasses::InteractionRecord const& record);
std::ostream& operator<<(std::ostream& os, LI::dataclasses::PrimaryDistributionRecord const& record);
std::ostream& operator<<(std::ostream& os, LI::dataclasses::SecondaryParticleRecord const& record);
std::ostream& operator<<(std::ostream& os, LI::dataclasses::CrossSectionDistributionRecord const& record);
std::ostream& operator<<(std::ostream& os, LI::dataclasses::SecondaryDistributionRecord const& record);

namespace LI {
namespace dataclasses {

PrimaryDistributionRecord::PrimaryDistributionRecord(ParticleType type) :
    id(ParticleID::GenerateID()),
    type(type)
{}

Particle PrimaryDistributionRecord::GetParticle() const {
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

double const & PrimaryDistributionRecord::GetMass() const {
    if(not mass_set) {
        UpdateMass();
    }
    return mass;
}

double const & PrimaryDistributionRecord::GetEnergy() const {
    if(not energy_set) {
        UpdateEnergy();
    }
    return energy;
}

double const & PrimaryDistributionRecord::GetKineticEnergy() const {
    if(not kinetic_energy_set) {
        UpdateKineticEnergy();
    }
    return kinetic_energy;
}


std::array<double, 3> const & PrimaryDistributionRecord::GetDirection() const {
    if(not direction_set) {
        UpdateDirection();
    }
    return direction;
}

std::array<double, 3> const & PrimaryDistributionRecord::GetThreeMomentum() const {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return momentum;
}

std::array<double, 4> PrimaryDistributionRecord::GetFourMomentum() const {
    if(not (momentum_set and energy_set)) {
        UpdateMomentum();
        UpdateEnergy();
    }
    return {energy, momentum.at(0), momentum.at(1), momentum.at(2)};
}

double const & PrimaryDistributionRecord::GetLength() const {
    if(not length_set) {
        UpdateLength();
    }
    return length;
}

std::array<double, 3> const & PrimaryDistributionRecord::GetInitialPosition() const {
    if(not initial_position_set) {
        UpdateInitialPosition();
    }
    return initial_position;
}

std::array<double, 3> const & PrimaryDistributionRecord::GetInteractionVertex() const {
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

void PrimaryDistributionRecord::UpdateMass() const {
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

void PrimaryDistributionRecord::UpdateEnergy() const {
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

void PrimaryDistributionRecord::UpdateKineticEnergy() const {
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

void PrimaryDistributionRecord::UpdateDirection() const {
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

void PrimaryDistributionRecord::UpdateMomentum() const {
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

void PrimaryDistributionRecord::UpdateLength() const {
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

void PrimaryDistributionRecord::UpdateInitialPosition() const {
    if(initial_position_set)
        return;
    if(interaction_vertex_set and direction_set and length_set) {
        initial_position = {interaction_vertex.at(0) - length*direction.at(0), interaction_vertex.at(1) - length*direction.at(1), interaction_vertex.at(2) - length*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate initial position without interaction vertex and direction and length!");
    }
}

void PrimaryDistributionRecord::UpdateInteractionVertex() const {
    if(interaction_vertex_set)
        return;
    if(initial_position_set and direction_set and length_set) {
        interaction_vertex = {initial_position.at(0) + length*direction.at(0), initial_position.at(1) + length*direction.at(1), initial_position.at(2) + length*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate interaction vertex without initial position and direction and length!");
    }
}

void PrimaryDistributionRecord::FinalizeAvailable(InteractionRecord & record) const {
    record.signature.primary_type = type;
    record.primary_id = GetID();
    try {
        record.primary_initial_position = GetInitialPosition();
    } catch(std::runtime_error e) {}
    try {
        record.interaction_vertex = GetInteractionVertex();
    } catch(std::runtime_error e) {}
    try {
        record.primary_mass = GetMass();
    } catch(std::runtime_error e) {}
    try {
        record.primary_momentum = GetFourMomentum();
    } catch(std::runtime_error e) {}
    try {
        record.primary_helicity = GetHelicity();
    } catch(std::runtime_error e) {}
}

void PrimaryDistributionRecord::Finalize(InteractionRecord & record) const {
    record.signature.primary_type = type;
    record.primary_id = GetID();
    record.interaction_vertex = GetInteractionVertex();
    record.primary_initial_position = GetInitialPosition();
    record.primary_mass = GetMass();
    record.primary_momentum = GetFourMomentum();
    record.primary_helicity = GetHelicity();
}

/////////////////////////////////////////

SecondaryParticleRecord::SecondaryParticleRecord(InteractionRecord const & record, size_t secondary_index) :
    secondary_index(secondary_index),
    id((record.secondary_ids.size() > secondary_index and record.secondary_ids.at(secondary_index)) ? (record.secondary_ids.at(secondary_index)) : (ParticleID::GenerateID())),
    type(record.signature.secondary_types.at(secondary_index)),
    initial_position(record.interaction_vertex)
{}

Particle SecondaryParticleRecord::GetParticle() const {
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

double const & SecondaryParticleRecord::GetMass() const {
    if(not mass_set) {
        UpdateMass();
    }
    return mass;
}

double const & SecondaryParticleRecord::GetEnergy() const {
    if(not energy_set) {
        UpdateEnergy();
    }
    return energy;
}

double const & SecondaryParticleRecord::GetKineticEnergy() const {
    if(not kinetic_energy_set) {
        UpdateKineticEnergy();
    }
    return kinetic_energy;
}

std::array<double, 3> const & SecondaryParticleRecord::GetDirection() const {
    if(not direction_set) {
        UpdateDirection();
    }
    return direction;
}

std::array<double, 3> const & SecondaryParticleRecord::GetThreeMomentum() const {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return momentum;
}

std::array<double, 4> SecondaryParticleRecord::GetFourMomentum() const {
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

void SecondaryParticleRecord::UpdateMass() const {
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

void SecondaryParticleRecord::UpdateEnergy() const {
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

void SecondaryParticleRecord::UpdateKineticEnergy() const {
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

void SecondaryParticleRecord::UpdateDirection() const {
    if(direction_set)
        return;
    if(momentum_set) {
        double magnitude = std::sqrt(momentum.at(0)*momentum.at(0) + momentum.at(1)*momentum.at(1) + momentum.at(2)*momentum.at(2));
        direction = {momentum.at(0)/magnitude, momentum.at(1)/magnitude, momentum.at(2)/magnitude};
    } else {
        throw std::runtime_error("Cannot calculate direction without momentum or initial position and interaction vertex!");
    }
}

void SecondaryParticleRecord::UpdateMomentum() const {
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

void SecondaryParticleRecord::Finalize(InteractionRecord & record) const {
    assert(record.signature.secondary_types.at(secondary_index) == type);
    record.secondary_ids.at(secondary_index) = GetID();
    record.secondary_masses.at(secondary_index) = GetMass();
    record.secondary_momenta.at(secondary_index) = GetFourMomentum();
    record.secondary_helicities.at(secondary_index) = GetHelicity();
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

std::vector<SecondaryParticleRecord> & CrossSectionDistributionRecord::GetSecondaryParticleRecords() {
    return secondary_particles;
}

std::vector<SecondaryParticleRecord> const & CrossSectionDistributionRecord::GetSecondaryParticleRecords() const {
    return secondary_particles;
}

SecondaryParticleRecord & CrossSectionDistributionRecord::GetSecondaryParticleRecord(size_t index) {
    return secondary_particles.at(index);
}

SecondaryParticleRecord const & CrossSectionDistributionRecord::GetSecondaryParticleRecord(size_t index) const {
    return secondary_particles.at(index);
}

void CrossSectionDistributionRecord::Finalize(InteractionRecord & record) const {
    record.target_id = target_id;
    record.target_mass = target_mass;
    record.target_helicity = target_helicity;

    record.interaction_parameters = interaction_parameters;

    record.secondary_ids.resize(secondary_particles.size());
    record.secondary_masses.resize(secondary_particles.size());
    record.secondary_momenta.resize(secondary_particles.size());
    record.secondary_helicities.resize(secondary_particles.size());

    for(SecondaryParticleRecord const & secondary : secondary_particles) {
        secondary.Finalize(record);
    }
}

/////////////////////////////////////////

InteractionRecord SecondaryDistributionRecord::CreateSecondaryRecord(InteractionRecord const & parent_record, size_t secondary_index) {
    InteractionRecord record;

    record.primary_id = ((parent_record.secondary_ids.at(secondary_index)) ? (parent_record.secondary_ids.at(secondary_index)) : (ParticleID::GenerateID())),
    record.signature.primary_type = parent_record.signature.secondary_types.at(secondary_index);
    record.primary_mass = parent_record.secondary_masses.at(secondary_index);
    record.primary_momentum = parent_record.secondary_momenta.at(secondary_index);
    record.primary_helicity = parent_record.secondary_helicities.at(secondary_index);
    record.primary_initial_position = parent_record.interaction_vertex;
    return record;
}

SecondaryDistributionRecord::SecondaryDistributionRecord(InteractionRecord & record) :
    record([](InteractionRecord & record) -> InteractionRecord & {
        record.primary_id = ((record.primary_id) ? (record.primary_id) : (ParticleID::GenerateID()));
        return record;
    }(record)),
    id(record.primary_id),
    type(record.signature.primary_type),
    mass(record.primary_mass),
    direction([](InteractionRecord const & record) -> std::array<double, 3> {
        std::array<double, 3> direction;
        if(record.primary_momentum.at(0) != 0) {
            double magnitude = std::sqrt(record.primary_momentum.at(1)*record.primary_momentum.at(1) + record.primary_momentum.at(2)*record.primary_momentum.at(2) + record.primary_momentum.at(3)*record.primary_momentum.at(3));
            direction = {record.primary_momentum.at(1)/magnitude, record.primary_momentum.at(2)/magnitude, record.primary_momentum.at(3)/magnitude};
        } else {
            direction = {0, 0, 0};
        }
        return direction;
    }(record)),
    momentum(record.primary_momentum),
    helicity(record.primary_helicity),
    initial_position(record.primary_initial_position) {}

SecondaryDistributionRecord::SecondaryDistributionRecord(InteractionRecord const & parent_record, size_t secondary_index) :
    secondary_index(secondary_index),
    record(SecondaryDistributionRecord::CreateSecondaryRecord(parent_record, secondary_index)),
    id(record.primary_id),
    type(record.signature.primary_type),
    mass(record.primary_mass),
    direction([](InteractionRecord const & record) -> std::array<double, 3> {
        std::array<double, 3> direction;
        if(record.primary_momentum.at(0) != 0) {
            double magnitude = std::sqrt(record.primary_momentum.at(1)*record.primary_momentum.at(1) + record.primary_momentum.at(2)*record.primary_momentum.at(2) + record.primary_momentum.at(3)*record.primary_momentum.at(3));
            direction = {record.primary_momentum.at(1)/magnitude, record.primary_momentum.at(2)/magnitude, record.primary_momentum.at(3)/magnitude};
        } else {
            direction = {0, 0, 0};
        }
        return direction;
    }(record)),
    momentum(record.primary_momentum),
    helicity(record.primary_helicity),
    initial_position(record.primary_initial_position) {}

double const & SecondaryDistributionRecord::GetLength() const {
    if(not length_set) {
        throw std::runtime_error("Length not set!");
    }
    return length;
}

void SecondaryDistributionRecord::SetLength(double const & length) {
    length_set = true;
    this->length = length;
}

void SecondaryDistributionRecord::Finalize(InteractionRecord & record) const {
    record.signature.primary_type = type;
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

std::ostream & operator<<(std::ostream & os, LI::dataclasses::PrimaryDistributionRecord const & record) {
    std::stringstream ss;
    ss << "PrimaryDistributionRecord (" << &record << ") ";
    os << ss.str() << '\n';

    ss.str(std::string());
    std::string id_str;
    ss << record.GetID();
    id_str = ss.str();
    std::string from = "\n";
    std::string to = "\n    ";
    size_t start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    os << "ID: " << id_str << "\n";

    os << "Type: " << record.GetType() << "\n";

    if(record.mass_set) {
        os << "Mass: " << record.GetMass() << "\n";
    } else {
        os << "Mass: " << "None" << "\n";
    }

    if(record.energy_set) {
        os << "Energy: " << record.GetEnergy() << "\n";
    } else {
        os << "Energy: " << "None" << "\n";
    }

    if(record.kinetic_energy_set) {
        os << "KineticEnergy: " << record.GetKineticEnergy() << "\n";
    } else {
        os << "KineticEnergy: " << "None" << "\n";
    }

    if(record.direction_set) {
        os << "Direction: " << record.GetDirection().at(0) << " " << record.GetDirection().at(1) << " " << record.GetDirection().at(2) << "\n";
    } else {
        os << "Direction: " << "None" << "\n";
    }

    if(record.momentum_set) {
        os << "Momentum: " << record.GetThreeMomentum().at(0) << " " << record.GetThreeMomentum().at(1) << " " << record.GetThreeMomentum().at(2) << "\n";
    } else {
        os << "Momentum: " << "None" << "\n";
    }

    if(record.length_set) {
        os << "Length: " << record.GetLength() << "\n";
    } else {
        os << "Length: " << "None" << "\n";
    }

    if(record.initial_position_set) {
        os << "InitialPosition: " << record.GetInitialPosition().at(0) << " " << record.GetInitialPosition().at(1) << " " << record.GetInitialPosition().at(2) << "\n";
    } else {
        os << "InitialPosition: " << "None" << "\n";
    }

    if(record.interaction_vertex_set) {
        os << "InteractionVertex: " << record.GetInteractionVertex().at(0) << " " << record.GetInteractionVertex().at(1) << " " << record.GetInteractionVertex().at(2) << "\n";
    } else {
        os << "InteractionVertex: " << "None" << "\n";
    }

    if(record.helicity_set) {
        os << "Helicity: " << record.GetHelicity() << "\n";
    } else {
        os << "Helicity: " << "None" << "\n";
    }

    return os;
}

std::ostream & operator<<(std::ostream & os, LI::dataclasses::CrossSectionDistributionRecord const & record) {
    std::stringstream ss;
    ss << "CrossSectionDistributionRecord (" << &record << ") ";
    os << ss.str() << '\n';

    ss.str(std::string());
    std::string id_str;
    ss << record.GetPrimaryID();
    id_str = ss.str();
    std::string from = "\n";
    std::string to = "\n    ";
    size_t start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    os << "PrimaryID: " << id_str << "\n";

    os << "PrimaryType: " << record.primary_type << "\n";

    os << "PrimaryInitialPosition: " << record.primary_initial_position.at(0) << " " << record.primary_initial_position.at(1) << " " << record.primary_initial_position.at(2) << "\n";

    os << "PrimaryMass: " << record.primary_mass << "\n";

    os << "PrimaryMomentum: " << record.primary_momentum.at(0) << " " << record.primary_momentum.at(1) << " " << record.primary_momentum.at(2) << " " << record.primary_momentum.at(3) << "\n";

    os << "PrimaryHelicity: " << record.primary_helicity << "\n";

    os << "InteractionVertex: " << record.interaction_vertex.at(0) << " " << record.interaction_vertex.at(1) << " " << record.interaction_vertex.at(2) << "\n";

    ss.str(std::string());
    ss << record.GetTargetID();
    id_str = ss.str();
    start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    os << "TargetID: " << id_str << "\n";

    os << "TargetType: " << record.target_type << "\n";

    os << "TargetMass: " << record.target_mass << "\n";

    os << "TargetHelicity: " << record.target_helicity << "\n";

    if(record.interaction_parameters.size() > 0) {
        os << "InteractionParameters:\n";
        for(auto const & parameter: record.interaction_parameters) {
            os << "\t" << parameter.first << ": " << parameter.second << "\n";
        }
    } else {
        os << "InteractionParameters: " << "None" << "\n";
    }

    os << "SecondaryParticles:\n";
    std::string secondary_str;
    for(size_t i = 0; i < record.signature.secondary_types.size(); ++i) {
        ss.str(std::string());
        ss << record.GetSecondaryParticleRecord(i);
        secondary_str = ss.str();
        start_pos = 0;
        while((start_pos = secondary_str.find(from, start_pos)) != std::string::npos) {
            secondary_str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
        }
        os << secondary_str << "\n";
    }

    return os;
}

std::ostream & operator<<(std::ostream & os, LI::dataclasses::SecondaryParticleRecord const & record) {
    std::stringstream ss;
    ss << "SecondaryParticleRecord (" << &record << ") ";
    os << ss.str() << '\n';

    ss.str(std::string());
    std::string id_str;
    ss << record.GetID();
    id_str = ss.str();
    std::string from = "\n";
    std::string to = "\n    ";
    size_t start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    os << "ID: " << id_str << "\n";

    os << "Type: " << record.GetType() << "\n";

    if(record.mass_set) {
        os << "Mass: " << record.mass << "\n";
    } else {
        os << "Mass: " << "None" << "\n";
    }

    if(record.energy_set) {
        os << "Energy: " << record.energy << "\n";
    } else {
        os << "Energy: " << "None" << "\n";
    }

    if(record.kinetic_energy_set) {
        os << "KineticEnergy: " << record.kinetic_energy << "\n";
    } else {
        os << "KineticEnergy: " << "None" << "\n";
    }

    if(record.direction_set) {
        os << "Direction: " << record.direction.at(0) << " " << record.direction.at(1) << " " << record.direction.at(2) << "\n";
    } else {
        os << "Direction: " << "None" << "\n";
    }

    if(record.momentum_set) {
        os << "Momentum: " << record.momentum.at(0) << " " << record.momentum.at(1) << " " << record.momentum.at(2) << "\n";
    } else {
        os << "Momentum: " << "None" << "\n";
    }

    os << "InitialPosition: " << record.initial_position.at(0) << " " << record.initial_position.at(1) << " " << record.initial_position.at(2) << "\n";

    if(record.helicity_set) {
        os << "Helicity: " << record.helicity << "\n";
    } else {
        os << "Helicity: " << "None" << "\n";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, LI::dataclasses::SecondaryDistributionRecord const& record) {
    std::stringstream ss;
    ss << "SecondaryDistributionRecord (" << &record << ") ";
    os << ss.str() << '\n';

    ss.str(std::string());
    std::string id_str;
    ss << record.id;
    id_str = ss.str();
    std::string from = "\n";
    std::string to = "\n    ";
    size_t start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    os << "ID: " << id_str << "\n";

    os << "Type: " << record.type << "\n";

    os << "Mass: " << record.mass << "\n";

    os << "Direction: " << record.direction.at(0) << " " << record.direction.at(1) << " " << record.direction.at(2) << "\n";

    os << "Momentum: " << record.momentum.at(0) << " " << record.momentum.at(1) << " " << record.momentum.at(2) << " " << record.momentum.at(3) << "\n";

    os << "Helicity: " << record.helicity << "\n";

    os << "InitialPosition: " << record.initial_position.at(0) << " " << record.initial_position.at(1) << " " << record.initial_position.at(2) << "\n";

    if(record.length_set) {
        os << "Length: " << record.GetLength() << "\n";
    } else {
        os << "Length: " << "None" << "\n";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, LI::dataclasses::InteractionRecord const& record) {
    std::stringstream ss;
    ss << "InteractionRecord (" << &record << ") ";
    os << ss.str() << '\n';
    os << "Signature(" << &record.signature << "): " << record.signature.primary_type << " + " << record.signature.target_type << " ->";
    for(auto secondary: record.signature.secondary_types) {
        os << " " << secondary;
    }
    os << "\n";

    ss.str(std::string());
    std::string id_str;
    ss << record.primary_id;
    id_str = ss.str();
    std::string from = "\n";
    std::string to = "\n    ";
    size_t start_pos = 0;
    while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
        id_str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    ss << "PrimaryID: " << id_str << "\n";
    os << "PrimaryInitialPosition: " << record.primary_initial_position.at(0) << " " << record.primary_initial_position.at(1) << " " << record.primary_initial_position.at(2) << "\n";
    os << "InteractionVertex: " << record.interaction_vertex.at(0) << " " << record.interaction_vertex.at(1) << " " << record.interaction_vertex.at(2) << "\n";
    os << "PrimaryMass: " << record.primary_mass << "\n";
    os << "PrimaryMomentum: " << record.primary_momentum.at(0) << " " << record.primary_momentum.at(1) << " " << record.primary_momentum.at(2) << " " << record.primary_momentum.at(3) << "\n";
    os << "TargetID: " << record.target_id << "\n";
    os << "TargetMass: " << record.target_mass << "\n";
    os << "SecondaryIDs:\n";
    for(auto const & secondary: record.secondary_ids) {
        ss.str(std::string());
        id_str.clear();
        ss << secondary;
        id_str = ss.str();
        start_pos = 0;
        while((start_pos = id_str.find(from, start_pos)) != std::string::npos) {
            id_str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
        }
        os << "\t" << id_str << "\n";
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

