#include "LeptonInjector/crosssections/Decay.h"

namespace LI {
namespace crosssections {

Decay::Decay() {}

bool Decay::operator==(Decay const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

double Decay::TotalDecayLength(dataclasses::InteractionRecord const & record) {
  return TotalDecayLength(record.signature.primary_type,record.primary_momentum[0]);
}

double Decay::TotalDecayLength(LI::dataclasses::Particle::ParticleType primary, double energy) {
    double tau = 1./TotalDecayWidth(primary,energy); // in inverse GeV
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    return p1.beta() * p1.gamma() * tau * LI::utilities::Constants::hbarc * LI::utilities::constants::cm;
}

} // namespace crosssections
} // namespace LI

