#include "SIREN/interactions/Decay.h"

#include <array>                                           // for array

#include <rk/rk.hh>                                        // for P4
#include <rk/geom3.hh>                                     // for Vector3

#include "SIREN/dataclasses/InteractionRecord.h"
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

double Decay::TotalDecayWidthAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {
    return TotalDecayWidthAllFinalStates(record.signature.primary_type);
}

double Decay::TotalDecayWidthAllFinalStates(siren::dataclasses::ParticleType const & primary) const {
    std::vector<siren::dataclasses::InteractionSignature> signatures = this->GetPossibleSignaturesFromParent(primary);
    siren::dataclasses::InteractionRecord fake_record;
    fake_record.signature.primary_type = primary;
    double total_decay_width = 0;
    for(auto signature : signatures) {
        fake_record.signature = signature;
        total_decay_width += this->TotalDecayWidth(fake_record);
    }
    return total_decay_width;
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

} // namespace interactions
} // namespace siren

