#pragma once
#ifndef SIREN_CharmCrossSectionHelpers_H
#define SIREN_CharmCrossSectionHelpers_H

// Shared, header-only helpers for the charm-DIS cross sections
// (QuarkDISFromSpline and PythiaDISCrossSection). These two classes are
// intentionally independent (no common base); this file holds only the byte-
// identical pieces, moved verbatim, so the values and behavior are unchanged.

#include <set>
#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <cctype>
#include <stdexcept>

#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/utilities/Constants.h"

namespace siren {
namespace interactions {
namespace charm_xsec {

// D0:D+/-:Ds = 0.60:0.23:0.15, renormalized to sum to 1.0. The Lambda_c channel
// (~0.02) is not modeled, so its fraction is redistributed into the implemented
// D species (each /0.98); otherwise summing TotalCrossSection over the three
// registered D signatures recovers only 0.98*sigma_inclusive and under-counts
// charm production. Values are load-bearing for interaction-depth normalization.
inline double FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) {
    if (secondary == siren::dataclasses::Particle::ParticleType::D0 || secondary == siren::dataclasses::Particle::ParticleType::D0Bar) {
        return 0.6 / 0.98;
    } else if (secondary == siren::dataclasses::Particle::ParticleType::DPlus || secondary == siren::dataclasses::Particle::ParticleType::DMinus) {
        return 0.23 / 0.98;
    } else if (secondary == siren::dataclasses::Particle::ParticleType::DsPlus || secondary == siren::dataclasses::Particle::ParticleType::DsMinus) {
        return 0.15 / 0.98;
    }
    return 0;
}

// Map a neutrino primary to its charged-lepton product (CC signature building).
// Throws on a non-neutrino input. The interaction-type policy (which classes
// accept CC/NC) stays in the caller; only this flavor map is shared.
inline siren::dataclasses::ParticleType ChargedLeptonProduct(siren::dataclasses::ParticleType primary_type) {
    if(primary_type == siren::dataclasses::ParticleType::NuE) {
        return siren::dataclasses::ParticleType::EMinus;
    } else if(primary_type == siren::dataclasses::ParticleType::NuEBar) {
        return siren::dataclasses::ParticleType::EPlus;
    } else if(primary_type == siren::dataclasses::ParticleType::NuMu) {
        return siren::dataclasses::ParticleType::MuMinus;
    } else if(primary_type == siren::dataclasses::ParticleType::NuMuBar) {
        return siren::dataclasses::ParticleType::MuPlus;
    } else if(primary_type == siren::dataclasses::ParticleType::NuTau) {
        return siren::dataclasses::ParticleType::TauMinus;
    } else if(primary_type == siren::dataclasses::ParticleType::NuTauBar) {
        return siren::dataclasses::ParticleType::TauPlus;
    } else {
        throw std::runtime_error("InitializeSignatures: Unknown parent neutrino type!");
    }
}

// Lepton mass by |PDG| (neutrinos -> 0); throws on a non-lepton input.
inline double GetLeptonMass(siren::dataclasses::ParticleType lepton_type) {
    int32_t lepton_number = std::abs(static_cast<int32_t>(lepton_type));
    switch(lepton_number) {
        case 11: return siren::utilities::Constants::electronMass;
        case 13: return siren::utilities::Constants::muonMass;
        case 15: return siren::utilities::Constants::tauMass;
        case 12: case 14: case 16: return 0;
        default: throw std::runtime_error("Unknown lepton type!");
    }
}

// Cross-section unit multiplier: "cm" -> 1, "m" -> 10000 (cm^2 internal).
inline double UnitForString(std::string units) {
    std::transform(units.begin(), units.end(), units.begin(),
        [](unsigned char c){ return std::tolower(c); });
    if(units == "cm") {
        return 1.0;
    } else if(units == "m") {
        return 10000.0;
    } else {
        throw std::runtime_error("Cross section units not supported!");
    }
}

inline std::vector<siren::dataclasses::ParticleType> ToVector(std::set<siren::dataclasses::ParticleType> const & types) {
    return std::vector<siren::dataclasses::ParticleType>(types.begin(), types.end());
}

inline std::vector<dataclasses::InteractionSignature> SignaturesForParents(
    std::map<std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType>, std::vector<dataclasses::InteractionSignature>> const & by_parent,
    siren::dataclasses::ParticleType primary_type,
    siren::dataclasses::ParticleType target_type) {
    std::pair<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType> key(primary_type, target_type);
    auto it = by_parent.find(key);
    if(it != by_parent.end()) {
        return it->second;
    }
    return std::vector<dataclasses::InteractionSignature>();
}

} // namespace charm_xsec
} // namespace interactions
} // namespace siren

#endif // SIREN_CharmCrossSectionHelpers_H
