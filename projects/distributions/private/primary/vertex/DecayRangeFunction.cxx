#include "SIREN/distributions/primary/vertex/DecayRangeFunction.h"

#include <cmath>
#include <tuple>              // for tie, operator<, operator==, tuple
#include <math.h>             // for sqrt
#include <algorithm>          // for min

namespace siren { namespace dataclasses { enum class ParticleType : int32_t; } }

namespace siren {
namespace distributions {

//---------------
// class DecayRangeFunction
//---------------

double DecayRangeFunction::DecayLength(double particle_mass, double decay_width, double energy) {
    double beta = std::sqrt(energy*energy - particle_mass*particle_mass) / energy;
    double gamma = energy / particle_mass;
    double time_in_rest_frame = 1.0 / decay_width; // inverse GeV
    double time_in_lab_frame = time_in_rest_frame * gamma; // inverse GeV
    constexpr double iGeV_in_m = 1.973269804593025e-16; // meters per inverse GeV
    double length = time_in_lab_frame * beta * iGeV_in_m; // meters = ((inverse GeV * dimensionless) * (meters per inverse GeV))
    return length; // meters
}

double DecayRangeFunction::DecayLength(siren::dataclasses::ParticleType const & primary_type, double energy) const {
    return DecayRangeFunction::DecayLength(particle_mass, decay_width, energy);
}

double DecayRangeFunction::Range(siren::dataclasses::ParticleType const & primary_type, double energy) const {
    return std::min(DecayLength(primary_type, energy) * multiplier, max_distance);
}

double DecayRangeFunction::operator()(siren::dataclasses::ParticleType const & primary_type, double energy) const {
    return Range(primary_type, energy);
}

double DecayRangeFunction::Multiplier() const {
    return multiplier;
}

double DecayRangeFunction::ParticleMass() const {
    return particle_mass;
}

double DecayRangeFunction::DecayWidth() const {
    return decay_width;
}

double DecayRangeFunction::MaxDistance() const {
    return max_distance;
}

DecayRangeFunction::DecayRangeFunction(double particle_mass, double decay_width, double multiplier, double max_distance) : particle_mass(particle_mass), decay_width(decay_width), multiplier(multiplier), max_distance(max_distance) {}

bool DecayRangeFunction::equal(RangeFunction const & other) const {
    const DecayRangeFunction* x = dynamic_cast<const DecayRangeFunction*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(particle_mass, decay_width, multiplier, max_distance)
            ==
            std::tie(x->particle_mass, x->decay_width, x->multiplier, x->max_distance);
}

bool DecayRangeFunction::less(RangeFunction const & other) const {
    const DecayRangeFunction* x = dynamic_cast<const DecayRangeFunction*>(&other);

    return
        std::tie(particle_mass, decay_width, multiplier, max_distance)
        <
        std::tie(x->particle_mass, x->decay_width, x->multiplier, x->max_distance);
}

} // namespace distributions
} // namespace sirenREN
