#include "SIREN/dataclasses/ParticleType.h"

#include <map>
#include <string>
#include <iostream>
#include <inttypes.h>

std::ostream & operator<<(std::ostream & os, siren::dataclasses::ParticleType && p) {
    siren::dataclasses::ParticleType pp = p;
    return operator<<(os, pp);
}

std::ostream & operator<<(std::ostream & os, siren::dataclasses::ParticleType const & p) {
    if (siren::dataclasses::ParticleTypeNames.find(p) != siren::dataclasses::ParticleTypeNames.end())
        os << siren::dataclasses::ParticleTypeNames.at(p);
    else
        os << static_cast<int32_t>(p);
    return os;
}
