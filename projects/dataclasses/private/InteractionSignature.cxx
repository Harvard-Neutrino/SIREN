#include "SIREN/dataclasses/InteractionSignature.h"

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

std::ostream& operator<<(std::ostream& os, siren::dataclasses::InteractionSignature const& signature) {
    std::stringstream ss;
    ss << "InteractionSignature (" << &signature << ") ";
    os << ss.str() << '\n';


    os << "PrimaryType: " << signature.primary_type << "\n";
    os << "TargetType: " << signature.target_type << "\n";
    os << "SecondaryTypes:";
    for(auto secondary: signature.secondary_types) {
        os << " " << secondary;
    }
    os << std::endl;

    return os;
}

