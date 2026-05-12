#pragma once
#ifndef SIREN_ParticleMasses_H
#define SIREN_ParticleMasses_H

#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/utilities/Constants.h"
#include <map>

namespace siren {
namespace dataclasses {

static const std::map<ParticleType, double> ParticleMasses = {
    {ParticleType::d,          siren::utilities::Constants::downMass},
    {ParticleType::dBar,       siren::utilities::Constants::downMass},
    {ParticleType::u,          siren::utilities::Constants::upMass},
    {ParticleType::uBar,       siren::utilities::Constants::upMass},
    {ParticleType::s,          siren::utilities::Constants::strangeMass},
    {ParticleType::sBar,       siren::utilities::Constants::strangeMass},
    {ParticleType::c,          siren::utilities::Constants::charmMass},
    {ParticleType::cBar,       siren::utilities::Constants::charmMass},
    {ParticleType::b,          siren::utilities::Constants::bottomMass},
    {ParticleType::bBar,       siren::utilities::Constants::bottomMass},
    {ParticleType::t,          siren::utilities::Constants::topMass},
    {ParticleType::tBar,       siren::utilities::Constants::topMass},
    {ParticleType::PPlus,      siren::utilities::Constants::protonMass},
    {ParticleType::PMinus,     siren::utilities::Constants::protonMass},
    {ParticleType::Neutron,    siren::utilities::Constants::neutronMass},
    {ParticleType::Gamma,      0.0},
    {ParticleType::EMinus,     siren::utilities::Constants::electronMass},
    {ParticleType::EPlus,      siren::utilities::Constants::electronMass},
    {ParticleType::MuMinus,    siren::utilities::Constants::muonMass},
    {ParticleType::MuPlus,     siren::utilities::Constants::muonMass},
    {ParticleType::TauMinus,   siren::utilities::Constants::tauMass},
    {ParticleType::TauPlus,    siren::utilities::Constants::tauMass},
    {ParticleType::Pi0,        siren::utilities::Constants::Pi0Mass},
    {ParticleType::PiPlus,     siren::utilities::Constants::PiPlusMass},
    {ParticleType::PiMinus,    siren::utilities::Constants::PiMinusMass},
    {ParticleType::K0_Long,    siren::utilities::Constants::K0Mass},
    {ParticleType::KPlus,      siren::utilities::Constants::KPlusMass},
    {ParticleType::KMinus,     siren::utilities::Constants::KMinusMass},
    {ParticleType::KPrime0,    siren::utilities::Constants::KPrime0Mass},
    {ParticleType::KPrimePlus, siren::utilities::Constants::KPrimePlusMass},
    {ParticleType::KPrimeMinus,siren::utilities::Constants::KPrimeMinusMass},
    {ParticleType::D0,         siren::utilities::Constants::D0Mass},
    {ParticleType::D0Bar,      siren::utilities::Constants::D0Mass},
    {ParticleType::DPlus,      siren::utilities::Constants::DPlusMass},
    {ParticleType::DMinus,     siren::utilities::Constants::DMinusMass},
    {ParticleType::DsPlus,     siren::utilities::Constants::DsPlusMass},
    {ParticleType::DsMinus,    siren::utilities::Constants::DsMinusMass},
    {ParticleType::Eta,        siren::utilities::Constants::EtaMass},
    {ParticleType::EtaPrime,   siren::utilities::Constants::EtaPrimeMass},
    {ParticleType::Rho0,       siren::utilities::Constants::Rho0Mass},
    {ParticleType::RhoPlus,    siren::utilities::Constants::RhoPlusMass},
    {ParticleType::RhoMinus,   siren::utilities::Constants::RhoMinusMass},
    {ParticleType::Omega,      siren::utilities::Constants::OmegaMass},
    {ParticleType::Phi,        siren::utilities::Constants::PhiMass},
    {ParticleType::BPlus,      siren::utilities::Constants::BPlusMass},
    {ParticleType::BMinus,     siren::utilities::Constants::BMinusMass},
    {ParticleType::WMinus,     siren::utilities::Constants::wMass},
    {ParticleType::WPlus,      siren::utilities::Constants::wMass},
    {ParticleType::Z0,         siren::utilities::Constants::zMass},
    {ParticleType::NuE,        0.0}, // Neutrinos are massless in this context
    {ParticleType::NuEBar,     0.0}, // Neutrinos are massless in this context
    {ParticleType::NuMu,       0.0}, // Neutrinos are massless in this context
    {ParticleType::NuMuBar,    0.0}, // Neutrinos are massless in this context
    {ParticleType::NuTau,      0.0}, // Neutrinos are massless in this context
    {ParticleType::NuTauBar,   0.0}, // Neutrinos are massless in this context
    {ParticleType::NuLight,    0.0}, // Neutrinos are massless in this context
    {ParticleType::NuLightBar, 0.0}, // Neutrinos are massless in this context
    // ...add more as needed
};

inline double GetParticleMass(ParticleType type) {
    auto it = ParticleMasses.find(type);
    if (it != ParticleMasses.end())
        return it->second;
    if(ParticleTypeNames.find(type) != ParticleTypeNames.end())
        throw std::invalid_argument("Particle type " + ParticleTypeNames.at(type) + " not found in mass map.");
    throw std::invalid_argument("Particle type with value " + std::to_string(static_cast<int>(type)) + " not found in mass map.");
}

} // namespace dataclasses
} // namespace siren

#endif // SIREN_ParticleMasses_H
