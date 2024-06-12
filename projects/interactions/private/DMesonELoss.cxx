#include "SIREN/interactions/DMesonELoss.h"

#include <map>                                             // for map, opera...
#include <set>                                             // for set, opera...
#include <array>                                           // for array
#include <cmath>                                           // for pow, log10
#include <tuple>                                           // for tie, opera...
#include <memory>                                          // for allocator
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <assert.h>                                        // for assert
#include <stddef.h>                                        // for size_t

#include <rk/rk.hh>                                        // for P4, Boost
#include <rk/geom3.hh>                                     // for Vector3

#include <photospline/splinetable.h>                       // for splinetable
//#include <photospline/cinter/splinetable.h>

#include "SIREN/interactions/CrossSection.h"     // for CrossSection
#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/dataclasses/Particle.h"           // for Particle
#include "SIREN/utilities/Random.h"               // for SIREN_random
#include "SIREN/utilities/Constants.h"            // for electronMass
#include "SIREN/utilities/Integration.h"          // for rombergInt...


namespace siren {
namespace interactions {

DMesonELoss::DMesonELoss() {
}

    
bool DMesonELoss::equal(CrossSection const & other) const {
    const DMesonELoss* x = dynamic_cast<const DMesonELoss*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(
            primary_types_,
            target_types_)
            ==
            std::tie(
            x->primary_types_,
            x->target_types_);
}



std::vector<siren::dataclasses::Particle::ParticleType> DMesonELoss::GetPossiblePrimaries() const {
    return std::vector<siren::dataclasses::Particle::ParticleType>(primary_types_.begin(), primary_types_.end());
}

// getting target should be the same for all the primary types
std::vector<siren::dataclasses::Particle::ParticleType> DMesonELoss::GetPossibleTargetsFromPrimary(siren::dataclasses::Particle::ParticleType primary_type) const {
    return std::vector<siren::dataclasses::Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<siren::dataclasses::Particle::ParticleType> DMesonELoss::GetPossibleTargets() const {
        return std::vector<siren::dataclasses::Particle::ParticleType>(target_types_.begin(), target_types_.end());
}

std::vector<dataclasses::InteractionSignature> DMesonELoss::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types_) {
        // hardcode the target type here, this should be fine
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParents(primary, siren::dataclasses::Particle::ParticleType::PPlus);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end()); 
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> DMesonELoss::GetPossibleSignaturesFromParents(siren::dataclasses::Particle::ParticleType primary_type, siren::dataclasses::Particle::ParticleType target_type) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary_type;
    signature.target_type = target_type;

    // first we deal with semileptonic decays where there are 3 final state particles
    signature.secondary_types.resize(1);
    signature.secondary_types[0] = primary_type; // same particle comes out
    signatures.push_back(signature);
    return signatures;
}

// i am here

double DMesonELoss::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    siren::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    // rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    // double primary_energy;

    // // make sure of the reference frame before assigning energy
    // if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
    //     primary_energy = interaction.primary_momentum[0];
    // } else {
    //     rk::Boost boost_start_to_lab = p2.restBoost();
    //     rk::P4 p1_lab = boost_start_to_lab * p1;
    //     primary_energy = p1_lab.e();
    // }
    double primary_energy = interaction.primary_momentum[0];


    return TotalCrossSection(primary_type, primary_energy);
}

double DMesonELoss::TotalCrossSection(siren::dataclasses::Particle::ParticleType primary_type, double primary_energy) const {
    if(not primary_types_.count(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    double log_energy = log10(primary_energy);
    double mb_to_cm2 = 1e-27;

    // current implementation uses only > 1PeV data
    double xsec = exp(1.891 + 0.205 * log_energy) - 2.157 + 1.264 * log_energy;

    return xsec * mb_to_cm2;
}

// double DMesonELoss::TotalCrossSection(siren::dataclasses::Particle::ParticleType primary_type, double primary_energy, siren::dataclasses::Particle::ParticleType target) const {
// 		return DMesonELoss::TotalCrossSection(primary_type,primary_energy);
// }


double DMesonELoss::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    // rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double primary_energy;
    rk::P4 p1_lab;
    // rk::P4 p2_lab;
    // if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
    primary_energy = interaction.primary_momentum[0];
    p1_lab = p1;
    // p2_lab = p2;
    // } else {
    //     rk::Boost boost_start_to_lab = p2.restBoost();
    //     p1_lab = boost_start_to_lab * p1;
    //     p2_lab = boost_start_to_lab * p2;
    //     primary_energy = p1_lab.e();
    //     std::cout << "D Meson Diff Xsec: not in lab frame???" << std::endl;
    // }

    double final_energy = interaction.secondary_momenta[0][0];
    double z = 1 - final_energy / primary_energy;
    
    // now normalize the gaussian
    double total_xsec = TotalCrossSection(interaction.signature.primary_type, primary_energy);
    double z0 = 0.56;
    double sigma = 0.2;
    std::function<double(double)> integrand = [&] (double z) -> double {
            return exp(-(pow(z - z0, 2))/(2 * pow(sigma, 2)));
        };
    double unnormalized = siren::utilities::rombergIntegrate(integrand, 0.001, 0.999);
    double normalization = total_xsec / unnormalized;

    double diff_xsec = normalization * exp(-(pow(z - z0, 2))/(2 * pow(sigma, 2)));

    return diff_xsec;
}


double DMesonELoss::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    // Consider implementing DIS thershold at some point
    return 0;
}

void DMesonELoss::SampleFinalState(dataclasses::CrossSectionDistributionRecord& interaction, std::shared_ptr<siren::utilities::SIREN_random> random) const {

    // std::cout << "In D Meson E Loss Sample Final State" << std::endl;

    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    // rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * tinteraction.secondary_momentarget_mass_ + 2 * target_mass_ * primary_energy;
    // double s = std::pow(rk::invMass(p1, p2), 2);

    double primary_energy;
    double Dmass = interaction.primary_mass;
    rk::P4 p1_lab;
    // rk::P4 p2_lab;
    // if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
    p1_lab = p1;
    // p2_lab = p2;
    primary_energy = p1_lab.e();
    // } else {
    //     // this is currently not implemented
    //     // Rest frame of p2 will be our "lab" frame
    //     rk::Boost boost_start_to_lab = p2.restBoost();
    //     p1_lab = boost_start_to_lab * p1;
    //     p2_lab = boost_start_to_lab * p2;
    //     primary_energy = p1_lab.e();
    //     // std::cout << "D Meson Energy Loss: not in lab frame???" << std::endl;
    // }
    // following line is wrong but i dont want to change it now fuck it.
    // std::cout << " " << interaction.primary_momentum[0] << " " << interaction.primary_momentum[1] << " " << interaction.primary_momentum[2] << " " << interaction.primary_momentum[3];
    // std::cout << primary_energy << " " << pow(primary_energy, 2) - pow(Dmass, 2) << " " <<
    //         sqrt(pow(interaction.primary_momentum[1], 2) +pow(interaction.primary_momentum[2], 2) +pow(interaction.primary_momentum[3], 2)) << std::endl;
    // sample an inelasticity from gaussian using Box-Muller Transform
    double sigma = 0.2;
    double z0 = 0.56; // for mesons only, for baryons it's 0.59, but not implemented yet
    double u1, u2;
    double final_energy;
    bool accept;


    do {
        do
        {
            u1 = random->Uniform(0, 1);
        }
        while (u1 == 0);
        u2 = random->Uniform(0, 1);
        double z = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2) + z0;
        // std::cout << z<< std::endl;

        // now modify the energy of the charm hadron and the corresponding momentum
        final_energy = primary_energy * (1-z);
        if (pow(final_energy, 2) - pow(Dmass, 2) >= 0) {
            accept = true;
        } else {
            accept = false;
        }
    } while (!accept);
    // this might be an infinite loop???????
    // need to check if the cross section length is good enough, how to make some relevant plots?
    // std::cout << final_energy << std::endl;
    double p3f = sqrt(pow(final_energy, 2) - pow(Dmass, 2));
    double p3i = std::sqrt(std::pow(p1.px(), 2) + std::pow(p1.py(), 2) + std::pow(p1.pz(), 2));
    double p_ratio = p3f / p3i;

    // std::cout << " " << p3f << " " << p3i << " " << p_ratio << std::endl;
    rk::P4 pf(p_ratio * geom3::Vector3(p1.px(), p1.py(), p1.pz()), Dmass);

    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = interaction.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & dmeson = secondaries[0];


    dmeson.SetFourMomentum({pf.e(), pf.px(), pf.py(), pf.pz()});
    dmeson.SetMass(pf.m());
    dmeson.SetHelicity(interaction.primary_helicity);

    // interaction.secondary_momenta.resize(1);
    // interaction.secondary_masses.resize(1);
    // interaction.secondary_helicity.resize(1);

    // interaction.secondary_momenta[0][0] = pf.e(); // p3_energy
    // interaction.secondary_momenta[0][1] = pf.px(); // p3_x
    // interaction.secondary_momenta[0][2] = pf.py(); // p3_y
    // interaction.secondary_momenta[0][3] = pf.pz(); // p3_z 
    // interaction.secondary_masses[0] = pf.m();

    // interaction.secondary_helicity[0] = interaction.primary_helicity;
}

double DMesonELoss::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<std::string> DMesonELoss::DensityVariables() const {
    return std::vector<std::string>{"Bjorken y"};
}

void DMesonELoss::InitializeSignatures() {
    return;
}

} // namespace interactions
} // namespace siren
