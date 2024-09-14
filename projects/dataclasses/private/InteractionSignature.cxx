#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/StringManipulation.h"

#include <tuple>    // for tie, operator<, operator==, tuple
#include <ostream>  // for operator<<, char_traits, basic_ostream, endl, ost...

namespace siren {
namespace dataclasses {

bool InteractionSignature::operator==(InteractionSignature const & other) const {
    return
        std::tie(primary_type, target_type, secondary_types)
        ==
        std::tie(other.primary_type, other.target_type, other.secondary_types);
}

bool InteractionSignature::operator<(InteractionSignature const & other) const {
    return
        std::tie(primary_type, target_type, secondary_types)
        <
        std::tie(other.primary_type, other.target_type, other.secondary_types);
}

} // namespace dataclasses
} // namespace siren

std::ostream& operator<<(std::ostream& os, siren::dataclasses::InteractionSignature const & signature) {
    os << to_repr(signature);
    return os;
}

std::string to_str(siren::dataclasses::InteractionSignature const & signature) {
    using siren::utilities::tab;
    std::stringstream ss;
    ss << "[ InteractionSignature (" << &signature << ") \n";
    ss << tab << "PrimaryType: " << signature.primary_type << '\n';
    ss << tab << "TargetType: " << signature.target_type << '\n';
    ss << tab << "SecondaryTypes:";
    for(auto secondary: signature.secondary_types) {
        ss << ' ' << secondary;
    }
    ss << "\n]";

    return ss.str();
}

std::string to_repr(siren::dataclasses::InteractionSignature const & signature) {
    using siren::dataclasses::ParticleType;
    std::stringstream ss;
    ss << "InteractionSignature( ";
    ss << signature.primary_type << " ";
    if(signature.primary_type == ParticleType::unknown or signature.target_type != ParticleType::unknown) {
        ss << signature.target_type << " ";
    }
    ss << "-> ";
    for(auto const & secondary : signature.secondary_types) {
        ss << secondary << " ";
    }
    ss << ")";
    return ss.str();
}
