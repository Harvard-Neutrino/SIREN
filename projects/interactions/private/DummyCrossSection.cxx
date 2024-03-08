#include "SIREN/interactions/DummyCrossSection.h"

#include <array>                                              // for array
#include <cmath>                                              // for sqrt, M_PI
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <stddef.h>                                           // for size_t

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "SIREN/interactions/CrossSection.h"        // for CrossSe...
#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/utilities/Random.h"                  // for LI_random

namespace SI {
namespace interactions {

DummyCrossSection::DummyCrossSection() {}


bool DummyCrossSection::equal(CrossSection const & other) const {
    const DummyCrossSection* x = dynamic_cast<const DummyCrossSection*>(&other);

    if(!x)
        return false;
    else
        return true;
}

double DummyCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    SI::dataclasses::ParticleType primary_type = interaction.signature.primary_type;
    SI::dataclasses::ParticleType target_type = interaction.signature.target_type;
    double primary_energy = interaction.primary_momentum[0];
    return TotalCrossSection(primary_type, primary_energy, target_type);
}

double DummyCrossSection::TotalCrossSection(SI::dataclasses::ParticleType primary_type, double primary_energy, SI::dataclasses::ParticleType target_type) const {
    double interaction_length_m = 1e6;
    double interaction_length_cm = interaction_length_m * 100;
    double density = 1.0; // g/cm^3
    double mol = 6.0221415e23; // particles/g
    double per_cm2 = density * interaction_length_cm * mol; // particles / cm^2
    double cm2 = 1.0 / per_cm2; // cm2/particle
    return cm2 * primary_energy / 1e5; // Interaction length was computed for 100TeV
}


double DummyCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    SI::dataclasses::ParticleType primary_type = interaction.signature.primary_type;
    SI::dataclasses::ParticleType target_type = interaction.signature.target_type;
    double primary_energy = interaction.primary_momentum[0];
    return TotalCrossSection(primary_type, primary_energy, target_type);
}

double DummyCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & interaction) const {
    return 0;
}

void DummyCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<SI::utilities::LI_random> random) const {
    rk::P4 p1(geom3::Vector3(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]), record.primary_mass);
    rk::P4 p2(geom3::Vector3(0, 0, 0), record.target_mass);

    // we assume that:
    // the target is stationary so its energy is just its mass
    // the incoming neutrino is massless, so its kinetic energy is its total energy
    // double s = target_mass_ * trecord.secondary_momentarget_mass_ + 2 * target_mass_ * primary_energy;
    // double s = std::pow(rk::invMass(p1, p2), 2);

    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    p1_lab = p1;
    p2_lab = p2;
    primary_energy = p1_lab.e();

    double final_x = random->Uniform(0, 1);
    double final_y = random->Uniform(0, 1);

    record.interaction_parameters.clear();
    record.interaction_parameters["energy"] = primary_energy;
    record.interaction_parameters["bjorken_x"] = final_x;
    record.interaction_parameters["bjorken_y"] = final_y;

    double m1 = record.primary_mass;
    double m3 = 0;
    double E1_lab = p1_lab.e();
    double E2_lab = p2_lab.e();
    double Q2 = 2 * E1_lab * E2_lab * final_x * final_y;
    double p1x_lab = std::sqrt(p1_lab.px() * p1_lab.px() + p1_lab.py() * p1_lab.py() + p1_lab.pz() * p1_lab.pz());
    double pqx_lab = (m1*m1 + m3*m3 + 2 * p1x_lab * p1x_lab + Q2 + 2 * E1_lab * E1_lab * (final_y - 1)) / (2.0 * p1x_lab);
    double momq_lab = std::sqrt(m1*m1 + p1x_lab*p1x_lab + Q2 + E1_lab * E1_lab * (final_y * final_y - 1));
    double pqy_lab = std::sqrt(momq_lab*momq_lab - pqx_lab *pqx_lab);
    double Eq_lab = E1_lab * final_y;

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_lab_dir, phi);

    rk::P4 pq_lab(Eq_lab, geom3::Vector3(pqx_lab, pqy_lab, 0));
    pq_lab.rotate(x_to_p1_lab_rot);
    pq_lab.rotate(rand_rot);

    rk::P4 p3_lab((p1_lab - pq_lab).momentum(), m3);
    rk::P4 p4_lab = p2_lab + pq_lab;

    rk::P4 p3;
    rk::P4 p4;
    p3 = p3_lab;
    p4 = p4_lab;

    size_t lepton_index = 0;
    size_t other_index = 1;

    std::vector<SI::dataclasses::SecondaryParticleRecord> & secondaries = record.GetSecondaryParticleRecords();
    SI::dataclasses::SecondaryParticleRecord & lepton = secondaries[lepton_index];
    SI::dataclasses::SecondaryParticleRecord & other = secondaries[other_index];

    lepton.SetFourMomentum({p3.e(), p3.px(), p3.py(), p3.pz()});
    lepton.SetMass(p3.m());
    lepton.SetHelicity(record.primary_helicity);

    other.SetFourMomentum({p4.e(), p4.px(), p4.py(), p4.pz()});
    other.SetMass(p4.m());
    other.SetHelicity(record.target_helicity);
}

double DummyCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<SI::dataclasses::ParticleType> DummyCrossSection::GetPossibleTargets() const {
    return {
        SI::dataclasses::ParticleType::Nucleon
    };
}

std::vector<SI::dataclasses::ParticleType> DummyCrossSection::GetPossibleTargetsFromPrimary(SI::dataclasses::ParticleType primary_type) const {
    return {
        SI::dataclasses::ParticleType::Nucleon
    };
}

std::vector<SI::dataclasses::ParticleType> DummyCrossSection::GetPossiblePrimaries() const {
    return {
        SI::dataclasses::ParticleType::NuE,
        SI::dataclasses::ParticleType::NuEBar,
        SI::dataclasses::ParticleType::NuMu,
        SI::dataclasses::ParticleType::NuMuBar,
        SI::dataclasses::ParticleType::NuTau,
        SI::dataclasses::ParticleType::NuTauBar
    };
}

std::vector<dataclasses::InteractionSignature> DummyCrossSection::GetPossibleSignatures() const {
    std::vector<SI::dataclasses::InteractionSignature> sigs;
    SI::dataclasses::InteractionSignature signature;
    std::vector<SI::dataclasses::ParticleType> targets = GetPossibleTargets();
    std::vector<SI::dataclasses::ParticleType> primaries = GetPossiblePrimaries();
    for(auto target : targets) {
        signature.target_type = target;
        for(auto primary : primaries) {
            signature.primary_type = primary;
            signature.secondary_types = {primary, target};
            sigs.push_back(signature);
        }
    }
    return sigs;
}

std::vector<dataclasses::InteractionSignature> DummyCrossSection::GetPossibleSignaturesFromParents(SI::dataclasses::ParticleType primary_type, SI::dataclasses::ParticleType target_type) const {
    std::vector<dataclasses::InteractionSignature> sigs = DummyCrossSection::GetPossibleSignatures();
    std::vector<dataclasses::InteractionSignature> these_sigs;
    for(auto sig : sigs) {
        if(sig.primary_type == primary_type)
            these_sigs.push_back(sig);
    }
    return these_sigs;
}

std::vector<std::string> DummyCrossSection::DensityVariables() const {
    return std::vector<std::string>{"Bjorken x", "Bjorken y"};
}

} // namespace interactions
} // namespace SI

