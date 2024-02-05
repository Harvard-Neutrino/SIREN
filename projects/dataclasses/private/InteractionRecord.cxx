#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include <tuple>    // for tie, operator==, tuple
#include <ostream>  // for operator<<, basic_ostream, char_traits, endl, ost...

namespace LI {
namespace dataclasses {

PrimaryRecord::PrimaryRecord(ParticleType const & type) : type(type) {}

PrimaryRecord & PrimaryRecord::operator=(PrimaryRecord const & record) {
    type = record.type;
    mass_set = record.mass_set;
    energy_set = record.energy_set;
    kinetic_energy_set = record.kinetic_energy_set;
    direction_set = record.direction_set;
    momentum_set = record.momentum_set;
    length_set = record.length_set;
    initial_position_set = record.initial_position_set;
    interaction_vertex_set = record.interaction_vertex_set;
    mass = record.mass;
    energy = record.energy;
    kinetic_energy = record.kinetic_energy;
    direction = record.direction;
    momentum = record.momentum;
    length = record.length;
    initial_position = record.initial_position;
    interaction_vertex = record.interaction_vertex;
    return *this;
}

ParticleType const & PrimaryRecord::GetType() const {
    return type;
}

double PrimaryRecord::GetMass() const {
    if(mass_set) {
        return mass;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateMass();
        return non_const_this.mass;
    }
}

double const & PrimaryRecord::GetMass() {
    if(not mass_set) {
        UpdateMass();
    }
    return mass;
}

double PrimaryRecord::GetEnergy() const {
    if(energy_set) {
        return energy;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateEnergy();
        return non_const_this.energy;
    }
}

double const & PrimaryRecord::GetEnergy() {
    if(not energy_set) {
        UpdateEnergy();
    }
    return energy;
}

double PrimaryRecord::GetKineticEnergy() const {
    if(kinetic_energy_set) {
        return kinetic_energy;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateKineticEnergy();
        return non_const_this.kinetic_energy;
    }
}

double const & PrimaryRecord::GetKineticEnergy() {
    if(not kinetic_energy_set) {
        UpdateKineticEnergy();
    }
    return kinetic_energy;
}


std::array<double, 3> PrimaryRecord::GetDirection() const {
    if(direction_set) {
        return direction;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateDirection();
        return non_const_this.direction;
    }
}

std::array<double, 3> const & PrimaryRecord::GetDirection() {
    if(not direction_set) {
        UpdateDirection();
    }
    return direction;
}

std::array<double, 3> PrimaryRecord::GetThreeMomentum() const {
    if(momentum_set) {
        return momentum;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateMomentum();
        return non_const_this.momentum;
    }
}

std::array<double, 3> const & PrimaryRecord::GetThreeMomentum() {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return momentum;
}

std::array<double, 4> PrimaryRecord::GetFourMomentum() const {
    if(momentum_set) {
        return {energy, momentum.at(0), momentum.at(1), momentum.at(2)};
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateMomentum();
        return {non_const_this.energy, non_const_this.momentum.at(0), non_const_this.momentum.at(1), non_const_this.momentum.at(2)};
    }
}

std::array<double, 4> PrimaryRecord::GetFourMomentum() {
    if(not momentum_set) {
        UpdateMomentum();
    }
    return {momentum.at(0), momentum.at(1), momentum.at(2), GetEnergy()};
}

double PrimaryRecord::GetLength() const {
    if(length_set) {
        return length;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateLength();
        return non_const_this.length;
    }
}

double const & PrimaryRecord::GetLength() {
    if(not length_set) {
        UpdateLength();
    }
    return length;
}

std::array<double, 3> PrimaryRecord::GetIntialPosition() const {
    if(initial_position_set) {
        return initial_position;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateIntialPosition();
        return non_const_this.initial_position;
    }
}

std::array<double, 3> const & PrimaryRecord::GetIntialPosition() {
    if(not initial_position_set) {
        UpdateIntialPosition();
    }
    return initial_position;
}

std::array<double, 3> PrimaryRecord::GetInteractionVertex() const {
    if(interaction_vertex_set) {
        return interaction_vertex;
    } else {
        PrimaryRecord non_const_this = *this;
        non_const_this.UpdateInteractionVertex();
        return non_const_this.interaction_vertex;
    }
}

std::array<double, 3> const & PrimaryRecord::GetInteractionVertex() {
    if(not interaction_vertex_set) {
        UpdateInteractionVertex();
    }
    return interaction_vertex;
}

// void PrimaryRecord::SetType(ParticleType type) {
//     this->type = type;
// }

void PrimaryRecord::SetMass(double mass) {
    mass_set = true;
    this->mass = mass;
}

void PrimaryRecord::SetEnergy(double energy) {
    energy_set = true;
    this->energy = energy;
}

void PrimaryRecord::SetKineticEnergy(double kinetic_energy) {
    kinetic_energy_set = true;
    this->kinetic_energy = kinetic_energy;
}

void PrimaryRecord::SetDirection(std::array<double, 3> direction) {
    direction_set = true;
    this->direction = direction;
}

void PrimaryRecord::SetThreeMomentum(std::array<double, 3> momentum) {
    momentum_set = true;
    this->momentum = momentum;
}

void PrimaryRecord::SetFourMomentum(std::array<double, 4> momentum) {
    momentum_set = true;
    this->momentum = {momentum.at(1), momentum.at(2), momentum.at(3)};
    energy_set = true;
    this->energy = momentum.at(0);
}

void PrimaryRecord::SetLength(double length) {
    length_set = true;
    this->length = length;
}

void PrimaryRecord::SetIntialPosition(std::array<double, 3> initial_position) {
    initial_position_set = true;
    this->initial_position = initial_position;
}

void PrimaryRecord::SetInteractionVertex(std::array<double, 3> interaction_vertex) {
    interaction_vertex_set = true;
    this->interaction_vertex = interaction_vertex;
}

void PrimaryRecord::UpdateMass() {
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

void PrimaryRecord::UpdateEnergy() {
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

void PrimaryRecord::UpdateKineticEnergy() {
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

void PrimaryRecord::UpdateDirection() {
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

void PrimaryRecord::UpdateMomentum() {
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

void PrimaryRecord::UpdateLength() {
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

void PrimaryRecord::UpdateIntialPosition() {
    if(initial_position_set)
        return;
    if(interaction_vertex_set and direction_set and length_set) {
        initial_position = {interaction_vertex.at(0) - length*direction.at(0), interaction_vertex.at(1) - length*direction.at(1), interaction_vertex.at(2) - length*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate initial position without interaction vertex and direction and length!");
    }
}

void PrimaryRecord::UpdateInteractionVertex() {
    if(interaction_vertex_set)
        return;
    if(initial_position_set and direction_set and length_set) {
        interaction_vertex = {initial_position.at(0) + length*direction.at(0), initial_position.at(1) + length*direction.at(1), initial_position.at(2) + length*direction.at(2)};
    } else {
        throw std::runtime_error("Cannot calculate interaction vertex without initial position and direction and length!");
    }
}

/////////////////////////////////////////

SecondaryRecord::SecondaryRecord(ParticleType const & type, std::array<double, 3> initial_position, double mass, std::array<double, 4> momentum) : PrimaryRecord(type) {
    SetIntialPosition(initial_position);
    SetMass(mass);
    SetFourMomentum(momentum);
}

SecondaryRecord & SecondaryRecord::operator=(SecondaryRecord const & record) {
    PrimaryRecord::operator=(record);
    return *this;
}

ParticleType const & SecondaryRecord::GetType() const {
    return PrimaryRecord::GetType();
}

double const & SecondaryRecord::GetMass() {
    return PrimaryRecord::GetMass();
}

double const & SecondaryRecord::GetEnergy() {
    return PrimaryRecord::GetEnergy();
}

double const & SecondaryRecord::GetKineticEnergy() {
    return PrimaryRecord::GetKineticEnergy();
}

std::array<double, 3> const & SecondaryRecord::GetDirection() {
    return PrimaryRecord::GetDirection();
}

std::array<double, 3> const & SecondaryRecord::GetThreeMomentum() {
    return PrimaryRecord::GetThreeMomentum();
}

std::array<double, 4> SecondaryRecord::GetFourMomentum() {
    return PrimaryRecord::GetFourMomentum();
}

double const & SecondaryRecord::GetLength() {
    return PrimaryRecord::GetLength();
}

std::array<double, 3> const & SecondaryRecord::GetIntialPosition() {
    return PrimaryRecord::GetIntialPosition();
}

std::array<double, 3> const & SecondaryRecord::GetInteractionVertex() {
    return PrimaryRecord::GetInteractionVertex();
}

void SecondaryRecord::SetLength(double length) {
    PrimaryRecord::SetLength(length);
}

/////////////////////////////////////////

InteractionRecord::InteractionRecord(InteractionSignature const & signature) : signature_set(true), signature(signature) {}

PrimaryRecord InteractionRecord::GetPrimaryRecord() const {
    PrimaryRecord record(signature.primary_type);
    record.SetMass(primary_mass);
    record.SetEnergy(primary_momentum.at(0));
    record.SetThreeMomentum({primary_momentum.at(1), primary_momentum.at(2), primary_momentum.at(3)});
    record.SetIntialPosition(primary_initial_position);
    record.SetInteractionVertex(interaction_vertex);
    return record;
}

void InteractionRecord::SetPrimaryRecord(PrimaryRecord const & record) {
    if(signature_set) {
        if(record.GetType() != signature.primary_type) {
            throw std::runtime_error("Primary particle type does not match signature!");
        }
    } else {
        signature.primary_type = record.GetType();
    }
    primary_mass = record.GetMass();
    primary_momentum = record.GetFourMomentum();
    primary_initial_position = record.GetIntialPosition();
    interaction_vertex = record.GetInteractionVertex();
}

InteractionSignature const & InteractionRecord::GetSignature() const {
    return signature;
}

Particle InteractionRecord::GetPrimary() const {
    double primary_length = std::sqrt(
        (interaction_vertex.at(0) - primary_initial_position.at(0))*(interaction_vertex.at(0) - primary_initial_position.at(0)) +
        (interaction_vertex.at(1) - primary_initial_position.at(1))*(interaction_vertex.at(1) - primary_initial_position.at(1)) +
        (interaction_vertex.at(2) - primary_initial_position.at(2))*(interaction_vertex.at(2) - primary_initial_position.at(2))
    );
    Particle primary(primary_id, signature.primary_type, primary_mass, primary_momentum, primary_initial_position, primary_length, primary_helicity);
    return primary;
}

ParticleID const & InteractionRecord::GetPrimaryID() const {
    return primary_id;
}

ParticleType const & InteractionRecord::GetPrimaryType() const {
    return signature.primary_type;
}

std::array<double, 3> const & InteractionRecord::GetPrimaryInitialPosition() const {
    return primary_initial_position;
}

double const & InteractionRecord::GetPrimaryMass() const {
    return primary_mass;
}

std::array<double, 4> const & InteractionRecord::GetPrimaryMomentum() const {
    return primary_momentum;
}

double const & InteractionRecord::GetPrimaryHelicity() const {
    return primary_helicity;
}

Particle InteractionRecord::GetTarget() const {
    Particle target(target_id, signature.target_type, target_mass, {0, 0, 0}, {0, 0, 0}, 0, target_helicity);
    return target;
}

ParticleID const & InteractionRecord::GetTargetID() const {
    return target_id;
}

ParticleType const & InteractionRecord::GetTargetType() const {
    return signature.target_type;
}

double const & InteractionRecord::GetTargetMass() const {
    return target_mass;
}

double const & InteractionRecord::GetTargetHelicity() const {
    return target_helicity;
}

std::array<double, 3> const & InteractionRecord::GetInteractionVertex() const {
    return interaction_vertex;
}

std::vector<Particle> InteractionRecord::GetSecondaries() const {
    std::vector<Particle> secondaries;
    for(size_t i=0; i<secondary_ids.size(); ++i) {
        secondaries.emplace_back(secondary_ids.at(i), signature.secondary_types.at(i), secondary_masses.at(i), secondary_momenta.at(i), std::array<double, 3>{0, 0, 0}, 0, secondary_helicities.at(i));
    }
    return secondaries;
}

std::vector<ParticleID> const & InteractionRecord::GetSecondaryIDs() const {
    return secondary_ids;
}

std::vector<Particle::ParticleType> const & InteractionRecord::GetSecondaryTypes() const {
    return signature.secondary_types;
}

std::vector<double> const & InteractionRecord::GetSecondaryMasses() const {
    return secondary_masses;
}

std::vector<std::array<double, 4>> const & InteractionRecord::GetSecondaryMomenta() const {
    return secondary_momenta;
}

std::vector<double> const & InteractionRecord::GetSecondaryHelicities() const {
    return secondary_helicities;
}

std::map<std::string, double> const & InteractionRecord::GetInteractionParameters() const {
    return interaction_parameters;
}

Particle InteractionRecord::GetSecondary(size_t const & index) const {
    if(index >= secondary_ids.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    Particle secondary(secondary_ids.at(index), signature.secondary_types.at(index), secondary_masses.at(index), secondary_momenta.at(index), std::array<double, 3>{0, 0, 0}, 0, secondary_helicities.at(index));
    return secondary;
}

ParticleID const & InteractionRecord::GetSecondaryID(size_t const & index) const {
    if(index >= secondary_ids.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    return secondary_ids.at(index);
}

ParticleType const & InteractionRecord::GetSecondaryType(size_t const & index) const {
    if(index >= signature.secondary_types.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    return signature.secondary_types.at(index);
}

double const & InteractionRecord::GetSecondaryMass(size_t const & index) const {
    if(index >= secondary_masses.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    return secondary_masses.at(index);
}

std::array<double, 4> const & InteractionRecord::GetSecondaryMomentum(size_t const & index) const {
    if(index >= secondary_momenta.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    return secondary_momenta.at(index);
}

double const & InteractionRecord::GetSecondaryHelicity(size_t const & index) const {
    if(index >= secondary_helicities.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    return secondary_helicities.at(index);
}

double const & InteractionRecord::GetInteractionParameter(std::string const & key) const {
    return interaction_parameters.at(key);
}

void InteractionRecord::SetSignature(InteractionSignature const & signature) {
    signature_set = true;
    this->signature = signature;
}

ParticleID InteractionRecord::SetPrimary(Particle & primary) {
    if(signature_set) {
        if(primary.type != ParticleType::unknown and primary.type != signature.primary_type) {
            throw std::runtime_error("Primary particle type does not match signature!");
        }
        reference_signature.primary_type = signature.primary_type;
    } else {
        signature.primary_type = primary.type;
    }

    if(primary.id) {
        primary_id = primary.id;
    } else if(not primary_id) {
        primary_id = primary.GenerateID();
    }

    std::array<double, 3> direction = {primary.momentum.at(1), primary.momentum.at(2), primary.momentum.at(3)};
    double magnitude = std::sqrt(direction.at(0)*direction.at(0) + direction.at(1)*direction.at(1) + direction.at(2)*direction.at(2));
    direction = {direction.at(0)/magnitude, direction.at(1)/magnitude, direction.at(2)/magnitude};

    primary_initial_position = primary.position;
    interaction_vertex = primary.position;
    interaction_vertex.at(0) += primary.length * direction.at(0);
    interaction_vertex.at(1) += primary.length * direction.at(1);
    interaction_vertex.at(2) += primary.length * direction.at(2);
    primary_mass = primary.mass;
    primary_momentum = primary.momentum;
    primary_helicity = primary.helicity;
    return primary_id;
}

ParticleID InteractionRecord::SetPrimary(Particle const & primary) {
    Particle primary_copy = primary;
    return SetPrimary(primary_copy);
}

ParticleID InteractionRecord::SetPrimaryID(ParticleID const & primary_id) {
    this->primary_id = primary_id;
    return primary_id;
}

void InteractionRecord::SetPrimaryType(Particle::ParticleType const & primary_type) {
    if(signature_set) {
        if(primary_type != ParticleType::unknown and primary_type != signature.primary_type) {
            throw std::runtime_error("Primary particle type does not match signature!");
        }
        reference_signature.primary_type = signature.primary_type;
    } else {
        signature.primary_type = primary_type;
    }
}

void InteractionRecord::SetPrimaryInitialPosition(std::array<double, 3> const & primary_initial_position) {
    this->primary_initial_position = primary_initial_position;
}

void InteractionRecord::SetPrimaryMass(double const & primary_mass) {
    this->primary_mass = primary_mass;
}

void InteractionRecord::SetPrimaryMomentum(std::array<double, 4> const & primary_momentum) {
    this->primary_momentum = primary_momentum;
}

void InteractionRecord::SetPrimaryHelicity(double const & primary_helicity) {
    this->primary_helicity = primary_helicity;
}

ParticleID InteractionRecord::SetTarget(Particle & target) {
    if(signature_set) {
        if(target.type != ParticleType::unknown and target.type != signature.target_type) {
            throw std::runtime_error("Target particle type does not match signature!");
        }
        reference_signature.target_type = signature.target_type;
    } else {
        signature.target_type = target.type;
    }

    if(target.id) {
        target_id = target.id;
    } else if(not target_id) {
        target_id = target.GenerateID();
    }

    target_mass = target.mass;
    target_helicity = target.helicity;
    return target_id;
}

ParticleID InteractionRecord::SetTarget(Particle const & target) {
    Particle target_copy = target;
    return SetTarget(target_copy);
}

ParticleID InteractionRecord::SetTargetID(ParticleID const & target_id) {
    this->target_id = target_id;
    return target_id;
}

void InteractionRecord::SetTargetType(Particle::ParticleType const & target_type) {
    if(signature_set) {
        if(target_type != ParticleType::unknown and target_type != signature.target_type) {
            throw std::runtime_error("Target particle type does not match signature!");
        }
        reference_signature.target_type = signature.target_type;
    } else {
        signature.target_type = target_type;
    }
}

void InteractionRecord::SetTargetMass(double const & target_mass) {
    this->target_mass = target_mass;
}

void InteractionRecord::SetTargetHelicity(double const & target_helicity) {
    this->target_helicity = target_helicity;
}

void InteractionRecord::SetInteractionVertex(std::array<double, 3> const & interaction_vertex) {
    this->interaction_vertex = interaction_vertex;
}

std::vector<ParticleID> InteractionRecord::SetSecondaries(std::vector<Particle> & secondaries) {
    secondary_ids.clear();
    secondary_masses.clear();
    secondary_momenta.clear();
    secondary_helicities.clear();
    for(Particle & secondary: secondaries) {
        AddSecondary(secondary);
    }
    return secondary_ids;
}

std::vector<ParticleID> InteractionRecord::SetSecondaries(std::vector<Particle> const & secondaries) {
    secondary_ids.clear();
    secondary_masses.clear();
    secondary_momenta.clear();
    secondary_helicities.clear();
    for(Particle const & secondary: secondaries) {
        AddSecondary(secondary);
    }
    return secondary_ids;
}

std::vector<ParticleID> InteractionRecord::SetSecondaryIDs(std::vector<ParticleID> const & secondary_ids) {
    this->secondary_ids = secondary_ids;
    return this->secondary_ids;
}

void InteractionRecord::SetSecondaryTypes(std::vector<Particle::ParticleType> const & secondary_types) {
    if(signature_set) {
        if(signature.secondary_types.size() != reference_signature.secondary_types.size()) {
            throw std::runtime_error("Secondary types do not match signature!");
        }
        for(size_t i=0; i<secondary_types.size(); ++i) {
            if(secondary_types.at(i) != ParticleType::unknown and secondary_types.at(i) != signature.secondary_types.at(i)) {
                throw std::runtime_error("Secondary types do not match signature!");
            }
            reference_signature.secondary_types.at(i) = signature.secondary_types.at(i);
        }
    } else {
        signature.secondary_types = secondary_types;
    }
}

void InteractionRecord::SetSecondaryMasses(std::vector<double> const & secondary_masses) {
    this->secondary_masses = secondary_masses;
}

void InteractionRecord::SetSecondaryMomenta(std::vector<std::array<double, 4>> const & secondary_momenta) {
    this->secondary_momenta = secondary_momenta;
}

void InteractionRecord::SetSecondaryHelicities(std::vector<double> const & secondary_helicities) {
    this->secondary_helicities = secondary_helicities;
}

void InteractionRecord::SetInteractionParameters(std::map<std::string, double> const & interaction_parameters) {
    this->interaction_parameters = interaction_parameters;
}

ParticleID InteractionRecord::SetSecondary(size_t const & index, Particle & secondary) {
    SetSecondaryType(index, secondary.type);
    SetSecondaryID(index, secondary.id);
    SetSecondaryMass(index, secondary.mass);
    SetSecondaryMomentum(index, secondary.momentum);
    SetSecondaryHelicity(index, secondary.helicity);
    return secondary_ids.at(index);
}

ParticleID InteractionRecord::SetSecondary(size_t const & index, Particle const & secondary) {
    Particle secondary_copy = secondary;
    return SetSecondary(index, secondary_copy);
}

ParticleID InteractionRecord::SetSecondaryID(size_t const & index, ParticleID const & secondary_id) {
    if(secondary_ids.size() < index) {
        throw std::runtime_error("Secondary index out of range!");
    } else if(secondary_ids.size() == index) {
        secondary_ids.resize(index+1);
    }
    secondary_ids.at(index) = secondary_id;
    return secondary_ids.at(index);
}

void InteractionRecord::SetSecondaryType(size_t const & index, Particle::ParticleType const & secondary_type) {
    if(signature_set) {
        if(secondary_type != ParticleType::unknown and secondary_type != signature.secondary_types.at(index)) {
            throw std::runtime_error("Secondary particle type does not match signature!");
        }
        if(reference_signature.secondary_types.size() == index) {
            reference_signature.secondary_types.push_back(signature.secondary_types.at(index));
        } else if(reference_signature.secondary_types.size() > index) {
            reference_signature.secondary_types.at(index) = signature.secondary_types.at(index);
        } else {
            throw std::runtime_error("Secondary index out of range!");
        }
    } else {
        signature.secondary_types.at(index) = secondary_type;
    }
}

void InteractionRecord::SetSecondaryMass(size_t const & index, double const & secondary_mass) {
    if(secondary_masses.size() < index) {
        throw std::runtime_error("Secondary index out of range!");
    } else if(secondary_masses.size() == index) {
        secondary_masses.resize(index+1);
    }
    secondary_masses.at(index) = secondary_mass;
}

void InteractionRecord::SetSecondaryMomentum(size_t const & index, std::array<double, 4> const & secondary_momentum) {
    if(secondary_momenta.size() < index) {
        throw std::runtime_error("Secondary index out of range!");
    } else if(secondary_momenta.size() == index) {
        secondary_momenta.resize(index+1);
    }
    secondary_momenta.at(index) = secondary_momentum;
}

void InteractionRecord::SetSecondaryHelicity(size_t const & index, double const & secondary_helicity) {
    if(this->secondary_helicities.size() < index) {
        throw std::runtime_error("Secondary index out of range!");
    } else if(this->secondary_helicities.size() == index) {
        this->secondary_helicities.resize(index+1);
    }
    this->secondary_helicities.at(index) = secondary_helicity;
}

void InteractionRecord::SetInteractionParameter(std::string const & key, double const & value) {
    interaction_parameters.at(key) = value;
}

ParticleID InteractionRecord::AddSecondary() {
    ParticleType secondary_type = ParticleType::unknown;
    if(signature_set) {
        secondary_type = signature.secondary_types.at(secondary_ids.size());
    }
    Particle secondary(secondary_type, 0, {0, 0, 0, 0}, {0, 0, 0}, 0, 0);
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(ParticleID const & secondary_id) {
    if(signature_set) {
        if(signature.secondary_types.size() < reference_signature.secondary_types.size() + 1) {
            throw std::runtime_error("Too many secondary particles!");
        }
        reference_signature.secondary_types.push_back(signature.secondary_types.at(secondary_ids.size()));
    } else {
        signature.secondary_types.push_back(ParticleType::unknown);
    }

    secondary_ids.push_back(secondary_id);
    secondary_masses.push_back(0);
    secondary_momenta.push_back({0, 0, 0, 0});
    secondary_helicities.push_back(0);
    return secondary_ids.back();
}

ParticleID InteractionRecord::AddSecondary(Particle & secondary) {
    if(signature_set) {
        if(signature.secondary_types.size() < reference_signature.secondary_types.size() + 1) {
            throw std::runtime_error("Too many secondary particles!");
        }
        if(secondary.type != ParticleType::unknown and secondary.type != signature.secondary_types.at(secondary_ids.size())) {
            throw std::runtime_error("Secondary particle type does not match signature!");
        }
        reference_signature.secondary_types.push_back(signature.secondary_types.at(secondary_ids.size()));
    } else {
        signature.secondary_types.push_back(secondary.type);
    }

    if(secondary.id) {
        secondary_ids.push_back(secondary.id);
    } else {
        secondary_ids.push_back(secondary.GenerateID());
    }

    secondary_masses.push_back(secondary.mass);
    secondary_momenta.push_back(secondary.momentum);
    secondary_helicities.push_back(secondary.helicity);
    return secondary_ids.back();
}

ParticleID InteractionRecord::AddSecondary(Particle const & secondary) {
    Particle secondary_copy = secondary;
    return AddSecondary(secondary_copy);
}

ParticleID InteractionRecord::AddSecondary(Particle::ParticleType const & secondary_type, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(secondary_type, secondary_mass, secondary_momentum, {0, 0, 0}, 0, secondary_helicity);
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(ParticleID const & secondary_id, Particle::ParticleType const & secondary_type, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(secondary_type, secondary_mass, secondary_momentum, {0, 0, 0}, 0, secondary_helicity);
    secondary.id = secondary_id;
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(Particle::ParticleType::unknown, secondary_mass, secondary_momentum, {0, 0, 0}, 0, secondary_helicity);
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(ParticleID const & secondary_id, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(Particle::ParticleType::unknown, secondary_mass, secondary_momentum, {0, 0, 0}, 0, secondary_helicity);
    secondary.id = secondary_id;
    return AddSecondary(secondary);
}

bool InteractionRecord::CheckSignature() const {
    if(signature.primary_type == ParticleType::unknown) {
        return false;
    }
    if(signature.target_type == ParticleType::unknown) {
        return false;
    }
    for(ParticleType secondary_type: signature.secondary_types) {
        if(secondary_type == ParticleType::unknown) {
            return false;
        }
    }
    if(signature_set) {
        if(signature.primary_type != reference_signature.primary_type) {
            return false;
        }
        if(signature.target_type != reference_signature.target_type) {
            return false;
        }
        if(signature.secondary_types.size() != reference_signature.secondary_types.size()) {
            return false;
        }
        for(size_t i=0; i<signature.secondary_types.size(); ++i) {
            if(signature.secondary_types.at(i) != reference_signature.secondary_types.at(i)) {
                return false;
            }
        }
    }
    return true;
}

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

