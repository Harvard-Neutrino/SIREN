#include "LeptonInjector/crosssections/Decay.h"

#include <rk/rk.hh>

#include "LeptonInjector/utilities/Constants.h"

namespace LI {
namespace crosssections {

Decay::Decay() {}

bool Decay::operator==(Decay const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

double Decay::TotalDecayLength(dataclasses::InteractionRecord const & interaction) const {
    double tau = 1./TotalDecayWidth(interaction); // in inverse GeV
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    return p1.beta() * p1.gamma() * tau * LI::utilities::Constants::hbarc;
}

double Decay::TotalDecayLengthForFinalState(dataclasses::InteractionRecord const & interaction) const {
    double tau = 1./TotalDecayWidthForFinalState(interaction); // in inverse GeV
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    return p1.beta() * p1.gamma() * tau * LI::utilities::Constants::hbarc;
}

} // namespace crosssections
} // namespace LI

