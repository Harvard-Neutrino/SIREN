#include <tuple>

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/DecaySignature.h"

namespace LI {
namespace dataclasses {

bool DecaySignature::operator==(DecaySignature const & other) const {
    return
        std::tie(primary_type, secondary_types)
        ==
        std::tie(other.primary_type, other.secondary_types);
}

} // namespace dataclasses
} // namespace LI
