#include "LeptonInjector/dataclasses/InteractionRecord.h"

namespace LI {
namespace dataclasses {

bool InteractionRecord::operator==(InteractionRecord const & other) const {
    return std::tie(
        signature,
        primary_mass,
        primary_momentum,
        primary_helicity,
        target_mass,
        target_momentum,
        target_helicity,
        interaction_vertex,
        secondary_masses,
        secondary_momenta,
        secondary_helicity,
        interaction_parameters)
        ==
        std::tie(
        other.signature,
        other.primary_mass,
        other.primary_momentum,
        other.primary_helicity,
        other.target_mass,
        other.target_momentum,
        other.target_helicity,
        other.interaction_vertex,
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

    os << "InteractionVertex: " << record.interaction_vertex[0] << " " << record.interaction_vertex[1] << " " << record.interaction_vertex[2] << "\n";
    os << "PrimaryMass: " << record.primary_mass << "\n";
    os << "PrimaryMomentum: " << record.primary_momentum[0] << " " << record.primary_momentum[1] << " " << record.primary_momentum[2] << " " << record.primary_momentum[3] << "\n";
    os << "TargetMass: " << record.target_mass << "\n";
    os << "TargetMomentum: " << record.target_momentum[0] << " " << record.target_momentum[1] << " " << record.target_momentum[2] << " " << record.target_momentum[3] << "\n";
    os << "SecondaryMomenta:\n";
    for(auto const & secondary: record.secondary_momenta) {
        os << "\t" << secondary[0] << " " << secondary[1] << " " << secondary[2] << " " << secondary[3] << "\n";
    }
    os << "SecondaryMasses:\n";
    for(auto const & secondary: record.secondary_masses) {
        os << "\t" << secondary << "\n";
    }
    os << "InteractionParameters:";
    for(auto param: record.interaction_parameters) {
        os << " " << param;
    }
    os << std::endl;

    return os;
}

