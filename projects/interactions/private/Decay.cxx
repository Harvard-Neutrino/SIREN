#include "SIREN/interactions/Decay.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>

#include <rk/rk.hh>                                        // for P4
#include <rk/geom3.hh>                                     // for Vector3

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleMasses.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"
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

siren::dataclasses::PhaseSpaceConvention Decay::Convention() const {
    using C = siren::dataclasses::PhaseSpaceConvention;
    auto variables = DensityVariables();

    auto lower_contains = [](std::vector<std::string> const & vars,
                             std::string const & needle) {
        for (auto const & v : vars) {
            std::string lv = v;
            std::transform(lv.begin(), lv.end(), lv.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (lv.find(needle) != std::string::npos) return true;
        }
        return false;
    };

    auto signatures = GetPossibleSignatures();
    if (signatures.empty()) {
        std::cerr << "Warning: Decay subclass does not override Convention(); "
                  << "auto-detected Custom from DensityVariables(). "
                  << "Override Convention() to silence this warning."
                  << std::endl;
        return C::Custom;
    }
    size_t n = signatures.front().secondary_types.size();
    for (auto const & sig : signatures) {
        if (sig.secondary_types.size() != n) {
            std::cerr << "Warning: Decay subclass does not override Convention(); "
                      << "auto-detected Custom (mixed arities) from "
                      << "DensityVariables(). Override Convention() to silence "
                      << "this warning." << std::endl;
            return C::Custom;
        }
    }

    C result = C::Custom;
    if (n == 2) {
        if (lower_contains(variables, "energy") ||
            lower_contains(variables, "e_") ||
            lower_contains(variables, "lab") ||
            lower_contains(variables, "cone") ||
            lower_contains(variables, "biased")) {
            result = C::Custom;
        } else {
            result = C::RestFrameSolidAngle;
        }
    } else if (n == 3) {
        if (lower_contains(variables, "dalitz") ||
            lower_contains(variables, "s12") ||
            lower_contains(variables, "s13") ||
            lower_contains(variables, "s_12")) {
            result = C::Dalitz;
        } else if (lower_contains(variables, "theta") ||
                   lower_contains(variables, "cos")) {
            result = C::HelicityAngles;
        }
    }

    std::cerr << "Warning: Decay subclass does not override Convention(); "
              << "auto-detected "
              << siren::dataclasses::PhaseSpaceConventionName(result)
              << " from DensityVariables(). Override Convention() to silence "
              << "this warning." << std::endl;
    return result;
}

} // namespace interactions
} // namespace siren
