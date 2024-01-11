#include "LeptonInjector/interactions/DarkNewsCrossSection.h"

#include <array>                                              // for array
#include <cmath>                                              // for sqrt, M_PI
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <stddef.h>                                           // for size_t

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "LeptonInjector/interactions/CrossSection.h"        // for CrossSe...
#include "LeptonInjector/dataclasses/InteractionRecord.h"     // for Interac...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle
#include "LeptonInjector/utilities/Random.h"                  // for LI_random
#include "LeptonInjector/utilities/Errors.h"                  // for PythonImplementationError

namespace LI {
namespace interactions {

DarkNewsCrossSection::DarkNewsCrossSection() {}

pybind11::object DarkNewsCrossSection::get_self() {
    return pybind11::cast<pybind11::none>(Py_None);
}

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

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p3_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p1_lab = p1;
        p3_lab = p2;
        primary_energy = p1_lab.e();
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p3_lab = boost_start_to_lab * p3;
        primary_energy = p1_lab.e();
    }
    double Q2 = -1*(std::pow(p1_lab.m(),2) + std::pow(p3_lab.m(),2) - 2.0*p1_lab.dot(p3_lab));

    return DifferentialCrossSection(interaction.signature.primary_type, interaction.signature.target_type, primary_energy, Q2);
}

double DarkNewsCrossSection::TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target
) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow TotalCrossSection to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::TotalCrossSection should be implemented in Python!"));
    return 0;
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

void DarkNewsCrossSection::SetUpscatteringMasses(dataclasses::InteractionRecord & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::SetUpscatteringMasses should be implemented in Python!"));
    return;
}

void DarkNewsCrossSection::SetUpscatteringHelicities(dataclasses::InteractionRecord & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsCrossSection::SetUpscatteringHelicities should be implemented in Python!"));
    return;
}


void DarkNewsCrossSection::SampleFinalState(dataclasses::InteractionRecord & interaction, std::shared_ptr<LI::utilities::LI_random> random) const {
    // Set our upscattering masses and helicities using values from DarkNews
    SetUpscatteringMasses(interaction);
    SetUpscatteringHelicities(interaction);
    interaction.primary_mass = 0;
    interaction.target_mass = m_target;
    interaction.secondary_masses.push_back(m_ups);
    interaction.secondary_masses.push_back(m_target);
    interaction.secondary_helicity.push_back(h_ups);
    interaction.secondary_helicity.push_back(h_target);
    // Uses Metropolis-Hastings Algorithm
    // Assumes we have the differential xsec v.s. Q^2

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = std::pow(rk::invMass(p1, p2), 2);

    // define masses that we will use
    double m1 = interaction.primary_mass;
    double m2 = interaction.target_mass;
    double m3 = interaction.secondary_masses[0];
    double m4 = interaction.secondary_masses[1];

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
    double minQ2 = Q2Min(interaction);
    double maxQ2 = Q2Max(interaction);
    double log_minQ2 = log10(minQ2);
    double log_maxQ2 = log10(maxQ2);

    bool accept;

    // kin_vars and its twin are 2-vectors containing [nu-energy, Q2]
    std::array<double,2> kin_vars, test_kin_vars;

    // values of cross_section from the table
    double cross_section, test_cross_section;

    // No matter what, we're evaluating at this specific energy.
    kin_vars[0] = test_kin_vars[0] = primary_energy;

    // sample an intial point
    kin_vars[1] = std::pow(10,random->Uniform(log_minQ2,log_maxQ2));
    
    test_cross_section = DifferentialCrossSection(interaction.signature.primary_type, interaction.signature.target_type, primary_energy, kin_vars[1]);
    cross_section = test_cross_section;

    // this is the magic part. Metropolis Hastings Algorithm.
    // MCMC method!
    const size_t burnin = 40; // converges to the correct distribution over multiple samplings.
                              // big number means more accurate, but slower
    for(size_t j = 0; j <= burnin; j++) {
        // repeat the sampling from above to get a new valid point
        test_kin_vars[1] = std::pow(10,random->Uniform(log_minQ2,log_maxQ2));
        test_cross_section = DifferentialCrossSection(interaction.signature.primary_type, interaction.signature.target_type, primary_energy, kin_vars[1]);
        double odds = (test_cross_section / cross_section);
        accept = (cross_section == 0 || (odds > 1.) || random->Uniform(0, 1) < odds);

        if(accept) {
            kin_vars = test_kin_vars;
            cross_section = test_cross_section;
        }
    }
    double final_Q2 = kin_vars[1];

    // // Working in center of mass frame, assuming 2 -> 2 scattering
    // double E1CM = (s + pow(interaction.primary_mass,2) - pow(interaction.target_mass,2)) / (2*sqrt(s));
    // double E3CM = (s + pow(interaction.secondary_masses[0],2) - pow(interaction.secondary_masses[1],2)) / (2*sqrt(s));
    // double P1CM = sqrt(E1CM*E1CM - pow(interaction.primary_mass,2));
    // double P3CM = sqrt(E3CM*E3CM - pow(interaction.secondary_masses[0],2));
                 
    // double CosThetaCM = (final_Q2 
    //                      + pow(interaction.primary_mass,2) 
    //                      + pow(interaction.secondary_masses[0],2)
    //                      - 2*E1CM*E3CM)
    //                     / (-2*P1CM*P3CM);

    // Leverage the fact that the target is at rest
    double E4Lab = (final_Q2 + m2*m2 + m4*m4) / (2*m2);
    double E3Lab = primary_energy + p2_lab.e() - E4Lab;
    double P1Lab = sqrt(primary_energy*primary_energy - m1*m1);
    double P3Lab = sqrt(E3Lab*E3Lab - m3*m3);
    double CosThetaLab = (primary_energy*E3Lab - 0.5*(final_Q2 + m1*m1 + m3*m3))/(P1Lab*P3Lab);
    double PhiLab = random->Uniform(0, 2.0 * M_PI);

    // Apply scattering angle in the x direction, then rotate to
    // incoming particle direction
    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);
    geom3::Rotation3 rand_rot(p1_lab_dir, PhiLab);

    rk::P4 p3_lab(E3Lab,P3Lab*geom3::Vector3(CosThetaLab,sqrt(1-CosThetaLab*CosThetaLab),0));
    p3_lab.rotate(x_to_p1_lab_rot);
    p3_lab.rotate(rand_rot);
    rk::P4 p4_lab = p1_lab + p2_lab - p3_lab;

    // Rotate back to whatever frame the traget was in originally.
    // I believe we have identified the lab frame as the target 
    // rest ferame in this function
    rk::P4 p3;
    rk::P4 p4;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p3 = p3_lab;
        p4 = p4_lab;
    } else {
        rk::Boost boost_lab_to_start = p2.labBoost();
        p3 = boost_lab_to_start * p3_lab;
        p4 = boost_lab_to_start * p4_lab;
    }

    // TODO: helicity update for secondary particles
    interaction.secondary_momenta.resize(2);

    interaction.secondary_momenta[0][0] = p3.e(); // p3_energy
    interaction.secondary_momenta[0][1] = p3.px(); // p3_x
    interaction.secondary_momenta[0][2] = p3.py(); // p3_y
    interaction.secondary_momenta[0][3] = p3.pz(); // p3_z

    interaction.secondary_momenta[1][0] = p4.e(); // p4_energy
    interaction.secondary_momenta[1][1] = p4.px(); // p4_x
    interaction.secondary_momenta[1][2] = p4.py(); // p4_y
    interaction.secondary_momenta[1][3] = p4.pz(); // p4_z

}

double DarkNewsCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    double dxs = DifferentialCrossSection(record);
    double txs = TotalCrossSection(record);
    if(dxs == 0) {
        return 0.0;
    } else if (txs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<std::string> DarkNewsCrossSection::DensityVariables() const {
    return std::vector<std::string>{"Q2"};
}

} // namespace interactions
} // namespace LI