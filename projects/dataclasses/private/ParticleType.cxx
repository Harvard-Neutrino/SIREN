#include "LeptonInjector/dataclasses/ParticleType.h"

#include <map>
#include <string>
#include <iostream>
#include <inttypes.h>

std::ostream & operator<<(std::ostream & os, LI::dataclasses::ParticleType && p) {
    LI::dataclasses::ParticleType pp = p;
    return operator<<(os, pp);
}

std::ostream & operator<<(std::ostream & os, LI::dataclasses::ParticleType const & p) {
    if (LI::dataclasses::ParticleTypeNames.find(p) != LI::dataclasses::ParticleTypeNames.end())
        os << LI::dataclasses::ParticleTypeNames.at(p);
    else
        os << static_cast<int32_t>(p);
    return os;
}
