#include "SIREN/interactions/Decay.h"

#include <array>                                           // for array

#include <rk/rk.hh>                                        // for P4
#include <rk/geom3.hh>                                     // for Vector3

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleMasses.h"
#include "SIREN/utilities/Constants.h"

namespace siren {
namespace interactions {

Decay::Decay() {}

bool Decay::operator==(Decay const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

double Decay::TotalDecayLengthAllFinalStates(dataclasses::InteractionRecord const & interaction) const {
    double tau = 1./TotalDecayWidthAllFinalStates(interaction); // in inverse GeV
    std::array<double, 4> const & p4 = interaction.primary_momentum;
    double const & mass = interaction.primary_mass;
    rk::P4 p1(geom3::Vector3(p4[1], p4[2], p4[3]), mass);
    return p1.beta() * p1.gamma() * tau * siren::utilities::Constants::hbarc;
}

double Decay::TotalDecayLength(dataclasses::InteractionRecord const & interaction) const {
    double tau = 1./TotalDecayWidth(interaction); // in inverse GeV
    std::array<double, 4> const & p4 = interaction.primary_momentum;
    double const & mass = interaction.primary_mass;
    rk::P4 p1(geom3::Vector3(p4[1], p4[2], p4[3]), mass);
    return p1.beta() * p1.gamma() * tau * siren::utilities::Constants::hbarc;
}

std::vector<double> Decay::SecondaryMasses(std::vector<siren::dataclasses::ParticleType> const & secondary_types) const {
    std::vector<double> masses;
    masses.reserve(secondary_types.size());
    for(auto const & type : secondary_types) {
        masses.push_back(siren::dataclasses::GetParticleMass(type));
    }
    return masses;
}

std::vector<double> Decay::SecondaryHelicities(dataclasses::InteractionRecord const & record) const {
    return std::vector<double>(record.signature.secondary_types.size(), 0.0);
}

} // namespace interactions
} // namespace siren
