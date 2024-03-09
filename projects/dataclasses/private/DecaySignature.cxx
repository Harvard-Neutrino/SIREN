#include "SIREN/dataclasses/DecaySignature.h"

#include <tuple>

namespace siren {
namespace dataclasses {

bool DecaySignature::operator==(DecaySignature const & other) const {
    return
        std::tie(primary_type, secondary_types)
        ==
        std::tie(other.primary_type, other.secondary_types);
}

} // namespace dataclasses
} // namespace siren
