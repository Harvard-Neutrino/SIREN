#include "SIREN/dataclasses/ParticleType.h"

#include <map>
#include <string>
#include <iostream>
#include <inttypes.h>

std::ostream & operator<<(std::ostream & os, SI::dataclasses::ParticleType && p) {
    SI::dataclasses::ParticleType pp = p;
    return operator<<(os, pp);
}

std::ostream & operator<<(std::ostream & os, SI::dataclasses::ParticleType const & p) {
    if (SI::dataclasses::ParticleTypeNames.find(p) != SI::dataclasses::ParticleTypeNames.end())
        os << SI::dataclasses::ParticleTypeNames.at(p);
    else
        os << static_cast<int32_t>(p);
    return os;
}
