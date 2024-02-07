#include "LeptonInjector/interactions/Decay.h"

#include <array>                                           // for array

#include <rk/rk.hh>                                        // for P4
#include <rk/geom3.hh>                                     // for Vector3

#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/utilities/Constants.h"

namespace LI {
namespace interactions {

Decay::Decay() {}

bool Decay::operator==(Decay const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

double Decay::TotalDecayLength(dataclasses::InteractionRecord const & interaction) const {
    double tau = 1./TotalDecayWidth(interaction); // in inverse GeV
    std::array<double, 4> const & p4 = interaction.primary_momentum;
    double const & mass = interaction.primary_mass;
    rk::P4 p1(geom3::Vector3(p4[1], p4[2], p4[3]), mass);
    return p1.beta() * p1.gamma() * tau * LI::utilities::Constants::hbarc;
}

double Decay::TotalDecayLengthForFinalState(dataclasses::InteractionRecord const & interaction) const {
    double tau = 1./TotalDecayWidthForFinalState(interaction); // in inverse GeV
    std::array<double, 4> const & p4 = interaction.primary_momentum;
    double const & mass = interaction.primary_mass;
    rk::P4 p1(geom3::Vector3(p4[1], p4[2], p4[3]), mass);
    return p1.beta() * p1.gamma() * tau * LI::utilities::Constants::hbarc;
}

} // namespace interactions
} // namespace LI

