#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include <tuple>    // for tie, operator==, tuple
#include <ostream>  // for operator<<, basic_ostream, char_traits, endl, ost...

namespace LI {
namespace dataclasses {

InteractionRecord::InteractionRecord(InteractionSignature const & signature) : signature_set(true), signature(signature) {}

InteractionSignature const & InteractionRecord::GetSignature() const {
    return signature;
}

Particle InteractionRecord::GetPrimary() const {
    Particle primary(primary_id, signature.primary_type, primary_mass, primary_momentum, primary_initial_position, primary_helicity);
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
    Particle target(target_id, signature.target_type, target_mass, target_momentum, {0, 0, 0}, target_helicity);
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

std::array<double, 4> const & InteractionRecord::GetTargetMomentum() const {
    return target_momentum;
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
        secondaries.emplace_back(secondary_ids.at(i), signature.secondary_types.at(i), secondary_masses.at(i), secondary_momenta.at(i), std::array<double, 3>{0, 0, 0}, secondary_helicity.at(i));
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

std::vector<double> const & InteractionRecord::GetSecondaryHelicity() const {
    return secondary_helicity;
}

std::map<std::string, double> const & InteractionRecord::GetInteractionParameters() const {
    return interaction_parameters;
}

Particle InteractionRecord::GetSecondary(size_t const & index) const {
    if(index >= secondary_ids.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    Particle secondary(secondary_ids.at(index), signature.secondary_types.at(index), secondary_masses.at(index), secondary_momenta.at(index), std::array<double, 3>{0, 0, 0}, secondary_helicity.at(index));
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
    if(index >= secondary_helicity.size()) {
        throw std::runtime_error("Secondary index out of range!");
    }
    return secondary_helicity.at(index);
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

    primary_initial_position = primary.position;
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
    target_momentum = target.momentum;
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

void InteractionRecord::SetTargetMomentum(std::array<double, 4> const & target_momentum) {
    this->target_momentum = target_momentum;
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
    secondary_helicity.clear();
    for(Particle & secondary: secondaries) {
        AddSecondary(secondary);
    }
    return secondary_ids;
}

std::vector<ParticleID> InteractionRecord::SetSecondaries(std::vector<Particle> const & secondaries) {
    secondary_ids.clear();
    secondary_masses.clear();
    secondary_momenta.clear();
    secondary_helicity.clear();
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

void InteractionRecord::SetSecondaryHelicity(std::vector<double> const & secondary_helicity) {
    this->secondary_helicity = secondary_helicity;
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
    if(this->secondary_helicity.size() < index) {
        throw std::runtime_error("Secondary index out of range!");
    } else if(this->secondary_helicity.size() == index) {
        this->secondary_helicity.resize(index+1);
    }
    this->secondary_helicity.at(index) = secondary_helicity;
}

void InteractionRecord::SetInteractionParameter(std::string const & key, double const & value) {
    interaction_parameters.at(key) = value;
}

ParticleID InteractionRecord::AddSecondary() {
    ParticleType secondary_type = ParticleType::unknown;
    if(signature_set) {
        secondary_type = signature.secondary_types.at(secondary_ids.size());
    }
    Particle secondary(secondary_type, 0, {0, 0, 0, 0}, {0, 0, 0}, 0);
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
    secondary_helicity.push_back(0);
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
    secondary_helicity.push_back(secondary.helicity);
    return secondary_ids.back();
}

ParticleID InteractionRecord::AddSecondary(Particle const & secondary) {
    Particle secondary_copy = secondary;
    return AddSecondary(secondary_copy);
}

ParticleID InteractionRecord::AddSecondary(Particle::ParticleType const & secondary_type, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(secondary_type, secondary_mass, secondary_momentum, {0, 0, 0}, secondary_helicity);
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(ParticleID const & secondary_id, Particle::ParticleType const & secondary_type, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(secondary_type, secondary_mass, secondary_momentum, {0, 0, 0}, secondary_helicity);
    secondary.id = secondary_id;
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(Particle::ParticleType::unknown, secondary_mass, secondary_momentum, {0, 0, 0}, secondary_helicity);
    return AddSecondary(secondary);
}

ParticleID InteractionRecord::AddSecondary(ParticleID const & secondary_id, double const & secondary_mass, std::array<double, 4> const & secondary_momentum, double const & secondary_helicity) {
    Particle secondary(Particle::ParticleType::unknown, secondary_mass, secondary_momentum, {0, 0, 0}, secondary_helicity);
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
        target_momentum,
        target_helicity,
        interaction_vertex,
        secondary_ids,
        secondary_masses,
        secondary_momenta,
        secondary_helicity,
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
        other.target_momentum,
        other.target_helicity,
        other.interaction_vertex,
        other.secondary_ids,
        other.secondary_masses,
        other.secondary_momenta,
        other.secondary_helicity,
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
        target_momentum,
        target_helicity,
        interaction_vertex,
        secondary_ids,
        secondary_masses,
        secondary_momenta,
        secondary_helicity,
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
        other.target_momentum,
        other.target_helicity,
        other.interaction_vertex,
        other.secondary_ids,
        other.secondary_masses,
        other.secondary_momenta,
        other.secondary_helicity,
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
    os << "TargetMomentum: " << record.target_momentum.at(0) << " " << record.target_momentum.at(1) << " " << record.target_momentum.at(2) << " " << record.target_momentum.at(3) << "\n";
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

