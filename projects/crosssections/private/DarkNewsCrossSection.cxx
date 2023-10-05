#include "LeptonInjector/crosssections/DarkNewsCrossSection.h"

#include <array>                                              // for array
#include <cmath>                                              // for sqrt, M_PI
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <stddef.h>                                           // for size_t

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "LeptonInjector/crosssections/CrossSection.h"        // for CrossSe...
#include "LeptonInjector/dataclasses/InteractionRecord.h"     // for Interac...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle
#include "LeptonInjector/utilities/Random.h"                  // for LI_random
#include "LeptonInjector/utilities/Errors.h"                  // for PythonImplementationError

namespace LI {
namespace crosssections {

DarkNewsCrossSection::DarkNewsCrossSection() {}

bool DarkNewsCrossSection::equal(CrossSection const & other) const {
    const DarkNewsCrossSection* x = dynamic_cast<const DarkNewsCrossSection*>(&other);

    if(!x)
        return false;
    else
        return true;
}

double DarkNewsCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    LI::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;
    LI::dataclasses::Particle::ParticleType target_type = interaction.signature.target_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double primary_energy;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        rk::P4 p1_lab = boost_start_to_lab * p1;
        primary_energy = p1_lab.e();
    }
    return TotalCrossSection(primary_type, primary_energy, target_type);
}

double DarkNewsCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    rk::P4 p3(geom3::Vector3(interaction.secondary_momenta[0][1], interaction.secondary_momenta[0][2], interaction.secondary_momenta[0][3]), interaction.secondary_masses[0]);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * target_mass_ + 2 * target_mass_ * primary_energy;
    double s = std::pow(rk::invMass(p1, p2), 2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p1_lab = p1;
        p2_lab = p2;
        primary_energy = p1_lab.e();
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p3_lab = boost_start_to_lab * p3;
        primary_energy = p1_lab.e();
    }
    double Q2 = -1*(std::pow(p1_lab.m(),2) + std::pow(p3_lab.m(),2) - 2.0*p1_lab.dot(p3_lab));

    return DifferentialCrossSection(interaction.signature.primary_type interaction.signature.target_type, primary_energy, Q2);
}

double DarkNewsCrossSection::DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary, LI::dataclasses::Particle::ParticleType target, double energy, double Q2
) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::DifferentialCrossSection should be implemented in Python!"));
    return 0;
}


double DarkNewsCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::InteractionThreshold should be implemented in Python!"));
    return 0;
}

double DarkNewsCrossSection::Q2Min(dataclasses::InteractionRecord const & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::Q2Min should be implemented in Python!"));
    return 0;
}

double DarkNewsCrossSection::Q2Max(dataclasses::InteractionRecord const & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::Q2Max should be implemented in Python!"));
    return 0;
}


void DarkNewsCrossSection::SampleFinalState(dataclasses::InteractionRecord & interaction) const {
    
    // Uses Metropolis-Hastings Algorithm
    // Assumes we have the differential xsec v.s. Q^2

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p1_lab = p1;
        p2_lab = p2;
        primary_energy = p1_lab.e();
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p2_lab = boost_start_to_lab * p2;
        primary_energy = p1_lab.e();
    }
    double Q2 = -1*(std::pow(p1_lab.m(),2) + std::pow(p2_lab.m(),2) - 2.0*p1_lab.dot(p2_lab));
    double minQ2 = Q2Min(interaction);
    double maxQ2 = Q2Max(interaction);
    double log_minQ2 = log10(Q2Min(interaction));
    double log_maxQ2 = log10(Q2Max(interaction));

    bool accept;

    // kin_vars and its twin are 2-vectors containing [nu-energy, Q2]
    std::array<double,2> kin_vars, test_kin_vars;

    // values of cross_section from the table
    double cross_section, test_cross_section;

    // No matter what, we're evaluating at this specific energy.
    kin_vars[0] = test_kin_vars[0] = primary_energy;

    // sample an intial point
    kin_vars[1] = std::pow(10,random->Uniform(log_minQ2,log_maxQ2));
    

    
}

} // namespace crosssections
} // namespace LI