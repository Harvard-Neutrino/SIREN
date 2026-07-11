#include "SIREN/interactions/Decay.h"

#include <algorithm>
#include <array>
#include <atomic>
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

double Decay::SampleDecayTime(dataclasses::CrossSectionDistributionRecord const & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // Identity default: keep the flight-time value already on the record.
    return record.GetInteractionTime();
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

siren::dataclasses::PhaseSpaceTopology Decay::Topology() const {
    using T = siren::dataclasses::PhaseSpaceTopology;
    auto signatures = GetPossibleSignatures();
    if (signatures.empty()) return T::Unspecified;
    size_t n = signatures.front().secondary_types.size();
    for (auto const & sig : signatures) {
        if (sig.secondary_types.size() != n) return T::Unspecified;
    }
    if (n == 2) return T::Decay2Body;
    if (n == 3) return T::Decay3Body;
    if (n > 3)  return T::DecayNBody;
    return T::Unspecified;
}

siren::dataclasses::PhaseSpaceMeasure Decay::Measure() const {
    using M = siren::dataclasses::PhaseSpaceMeasure;
    auto variables = DensityVariables();

    // Warn at most once per process that this base fallback is being used.
    static std::atomic<bool> warned{false};
    auto warn_once = [](std::string const & detail) {
        if (!warned.exchange(true)) {
            std::cerr << "Warning: Decay subclass does not override Measure(); "
                      << detail << " Override Measure() to silence this warning."
                      << std::endl;
        }
    };

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
        warn_once("auto-detected Unspecified (no signatures defined).");
        return M::Unspecified();
    }
    size_t n = signatures.front().secondary_types.size();
    for (auto const & sig : signatures) {
        if (sig.secondary_types.size() != n) {
            warn_once("auto-detected Unspecified (mixed secondary counts).");
            return M::Unspecified();
        }
    }

    M result = M::Unspecified();
    if (n == 2) {
        if (lower_contains(variables, "energy") ||
            lower_contains(variables, "e_") ||
            lower_contains(variables, "lab") ||
            lower_contains(variables, "cone") ||
            lower_contains(variables, "biased")) {
            result = M::Unspecified();
        } else {
            result = M::SolidAngleRest();
        }
    } else if (n == 3) {
        if (lower_contains(variables, "dalitz") ||
            lower_contains(variables, "s12") ||
            lower_contains(variables, "s13") ||
            lower_contains(variables, "s_12")) {
            result = M::DalitzPair();
        } else if (lower_contains(variables, "theta") ||
                   lower_contains(variables, "cos")) {
            result = M::HelicityAngles();
        }
    }

    warn_once("auto-detected " +
              siren::dataclasses::PhaseSpaceMeasureName(result) +
              " from DensityVariables().");
    return result;
}

siren::dataclasses::PhaseSpaceTopology Decay::TopologyForSignature(
    siren::dataclasses::InteractionSignature const & signature) const
{
    using T = siren::dataclasses::PhaseSpaceTopology;
    size_t n = signature.secondary_types.size();
    if (n == 2) return T::Decay2Body;
    if (n == 3) return T::Decay3Body;
    if (n > 3)  return T::DecayNBody;
    return T::Unspecified;
}

siren::dataclasses::PhaseSpaceMeasure Decay::MeasureForSignature(
    siren::dataclasses::InteractionSignature const & signature) const
{
    using M = siren::dataclasses::PhaseSpaceMeasure;
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

    size_t n = signature.secondary_types.size();
    if (n == 2) {
        if (lower_contains(variables, "energy") ||
            lower_contains(variables, "e_") ||
            lower_contains(variables, "lab") ||
            lower_contains(variables, "cone") ||
            lower_contains(variables, "biased")) {
            return M::Unspecified();
        }
        return M::SolidAngleRest();
    }
    if (n == 3) {
        if (lower_contains(variables, "dalitz") ||
            lower_contains(variables, "s12") ||
            lower_contains(variables, "s13") ||
            lower_contains(variables, "s_12")) {
            return M::DalitzPair();
        }
        if (lower_contains(variables, "theta") ||
            lower_contains(variables, "cos")) {
            return M::HelicityAngles();
        }
    }
    return M::Unspecified();
}

} // namespace interactions
} // namespace siren
