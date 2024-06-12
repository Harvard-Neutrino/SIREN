#include "SIREN/interactions/CharmHadronization.h"

#include <array>                                              // for array
#include <cmath>                                              // for sqrt, M_PI
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <stddef.h>                                           // for size_t

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "SIREN/interactions/Hadronization.h"        // for Hadronization
#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/utilities/Random.h"                  // for SIREN_random
#include "SIREN/utilities/Errors.h"                  // for PythonImplementationError
#include "SIREN/utilities/Constants.h"            // for electronMass



namespace siren {
namespace interactions {

CharmHadronization::CharmHadronization() {}

// pybind11::object CharmHadronization::get_self() {
//     return pybind11::cast<pybind11::none>(Py_None);
// }

bool CharmHadronization::equal(Hadronization const & other) const {
    const CharmHadronization* x = dynamic_cast<const CharmHadronization*>(&other);

    if(!x)
        return false;
    else
        return primary_types == x->primary_types;
}

double CharmHadronization::getHadronMass(siren::dataclasses::ParticleType hadron_type) {
    switch(hadron_type){
			case siren::dataclasses::ParticleType::D0:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::D0Bar:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::DPlus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::DMinus:
				return( siren::utilities::Constants::DPlusMass);	
			case siren::dataclasses::ParticleType::Charm:
				return( siren::utilities::Constants::CharmMass);
			case siren::dataclasses::ParticleType::CharmBar:
				return( siren::utilities::Constants::CharmMass);	
            default:
                return(0.0);
        }
}


std::vector<dataclasses::InteractionSignature> CharmHadronization::GetPossibleSignatures() const
{
    std::vector<dataclasses::InteractionSignature> signatures;
    for(auto primary : primary_types) {
      std::vector<dataclasses::InteractionSignature> new_signatures = GetPossibleSignaturesFromParent(primary);
      signatures.insert(signatures.end(),new_signatures.begin(),new_signatures.end()); 
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> CharmHadronization::GetPossibleSignaturesFromParent(siren::dataclasses::Particle::ParticleType primary) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.primary_type = primary;
    signature.target_type = siren::dataclasses::Particle::ParticleType::Decay;
    
    signature.secondary_types.resize(2);
    signature.secondary_types[0] = siren::dataclasses::Particle::ParticleType::Hadrons;
    if(primary==siren::dataclasses::Particle::ParticleType::Charm) {
        signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::D0;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::DPlus;
        signatures.push_back(signature);
      }
    else if(primary==siren::dataclasses::Particle::ParticleType::CharmBar) {
        signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::D0Bar;
        signatures.push_back(signature);
        signature.secondary_types[1] = siren::dataclasses::Particle::ParticleType::DMinus;
        signatures.push_back(signature);
    }
    return signatures;
}

void CharmHadronization::SampleFinalState(dataclasses::CrossSectionDistributionRecord & interaction, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // Uses Perterson distribution with epsilon = 0.2
    // charm quark momentum
    rk::P4 pc(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    double p3c = std::sqrt(std::pow(interaction.primary_momentum[1], 2) + std::pow(interaction.primary_momentum[2], 2) + std::pow(interaction.primary_momentum[3], 2));
    double Ec = pc.e(); //energy of primary charm
    // std::cout << "hadronization sample final state" << std::endl;
    // double peterson_distribution(double z) {
    //     return 0.8 / x * std::pow((1 - 1 / x - 0.2 / (1 - x)), 2)
    // }

    double z = 0.6; // replace by actual sampling: inverse CDF
    double ECH = z * Ec;

    // is it ok to compute everything in lab frame?
    double mCH = getHadronMass(interaction.signature.secondary_types[1]); // obtain charmed hadron mass
    double p3CH = std::sqrt(std::pow(ECH, 2) - std::pow(mCH, 2)); //obtain charmed hadron 3-momentum
    double rCH = p3CH/p3c; // ratio of momentum carried away by the charmed hadron, assume collinearity
    rk::P4 p4CH(geom3::Vector3(rCH * interaction.primary_momentum[1], rCH * interaction.primary_momentum[2], rCH * interaction.primary_momentum[3]), mCH);

    double EX = (1 - z) * Ec; // energy of the hadronic shower
    double p3X = EX; // assume no hadronic mass
    double rX = p3X/p3c; // assume collinear
    rk::P4 p4X(geom3::Vector3(rX * interaction.primary_momentum[1], rX * interaction.primary_momentum[2], rX * interaction.primary_momentum[3]), 0);


    // update interaction parameters: to be added here later 

    // new implementation of updateing outgoing particles
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondaries = interaction.GetSecondaryParticleRecords();
    siren::dataclasses::SecondaryParticleRecord & hadronic_vertex = secondaries[0];
    siren::dataclasses::SecondaryParticleRecord & d_meson = secondaries[1]; // these indices are hard-coded, should be automated in a future time

    hadronic_vertex.SetFourMomentum({p4X.e(), p4X.px(), p4X.py(), p4X.pz()});
    hadronic_vertex.SetMass(p4X.m());
    hadronic_vertex.SetHelicity(interaction.primary_helicity);

    d_meson.SetFourMomentum({p4CH.e(), p4CH.px(), p4CH.py(), p4CH.pz()});
    d_meson.SetMass(p4CH.m());
    d_meson.SetHelicity(interaction.primary_helicity);

    // interaction.secondary_momenta.resize(2);
    // interaction.secondary_masses.resize(2);
    // interaction.secondary_helicities.resize(2);

    // // the hadronic shower
    // interaction.secondary_momenta[0][0] = p4X.e();
    // interaction.secondary_momenta[0][1] = p4X.px();
    // interaction.secondary_momenta[0][2] = p4X.py();
    // interaction.secondary_momenta[0][3] = p4X.pz();
    // interaction.secondary_masses[0] = 0;
    // interaction.secondary_helicities[0] = interaction.primary_helicity; // true?

    // // the charmed hadron
    // interaction.secondary_momenta[1][0] = p4CH.e();
    // interaction.secondary_momenta[1][1] = p4CH.px();
    // interaction.secondary_momenta[1][2] = p4CH.py();
    // interaction.secondary_momenta[1][3] = p4CH.pz();
    // interaction.secondary_masses[1] = p4CH.m();
    // interaction.secondary_helicities[1] = interaction.primary_helicity; // true?

}

double CharmHadronization::FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) const {
    if (secondary == siren::dataclasses::Particle::ParticleType::D0 || secondary == siren::dataclasses::Particle::ParticleType::D0Bar) {
        return 0.6;
    } else if (secondary == siren::dataclasses::Particle::ParticleType::DPlus || secondary == siren::dataclasses::Particle::ParticleType::DMinus) {
        return 0.23;
    } // D_s and Lambda^+ not yet implemented
    return 0;
}

} // namespace interactions
} // namespace siren