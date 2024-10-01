#include "SIREN/interactions/MarleyCrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/Particle.h"
#include "marley/Generator.hh"  // MARLEY generator
#include "marley/JSON.hh"  // For handling MARLEY JSON configuration


namespace siren {
namespace interactions {

//MarleyCrossSection::MarleyCrossSection() {}
MarleyCrossSection::MarleyCrossSection(const std::string& marley_config) {
    InitializeMarley(marley_config);
}

void MarleyCrossSection::InitializeMarley(const std::string& marley_config) {
    marley::JSON marley_json = marley::JSON::load(marley_config);
    marley::JSONConfig config(marley_json);
    marley_generator_ = config.create_generator();
    marley_config_ = marley_config;

    std::vector<std::unique_ptr<marley::Reaction>> & reactions = marley_generator_.get_reactions();

    has_nu_cc = false;
    has_nubar_cc = false;
    has_nc = false;
    has_elastic = false;

    for(std::unique_ptr<Reaction> & reaction : reactions) {
        marley::Reaction::ProcessType process = reactions->process_type();
        switch (process) {
            case marley::Reaction::ProcessType::NeutrinoCC:
                has_nu_cc = true;
                break;
            case marley::Reaction::ProcessType::AntiNeutrinoCC:
                has_nubar_cc = true;
                break;
            case marley::Reaction::ProcessType::NC:
                has_nc = true;
                break;
            case marley::Reaction::ProcessType::NuElectronElastic:
                has_elastic = true;
                break;
            default:
                std::cout << "Unknown process type" << std::endl;
                break;
        }
    }
}



double MarleyCrossSection::TotalCrossSection(siren::dataclasses::InteractionRecord const & record) const {

    double px = record.primary_momentum[1];  // x component of the primary particle momentum
    double py = record.primary_momentum[2];  // y component of the primary particle momentum
    double pz = record.primary_momentum[3];  // z component of the primary particle momentum

    int pdg_a = record.signature.primary_type;  // PDG code of the primary particle
    double KEa = std::sqrt(px*px + py*py + pz*pz);  // Kinetic energy of the primary particle
    int pdg_atom = record.signature.target_type;  // PDG code of the target atom

    double KEa_MeV = KEa * 1e3;  // Convert kinetic energy to MeV

    double total_xs_in_MeV2 = marley_generator_.total_xs(pdg_a, KEa, pdg_atom); //in MeV^(-2)

    double conversion_factor = 3.893793721718594e-22;  // MeV^(-2) to cm^2
    double total_xs_in_cm2 = total_xs_in_MeV2 * conversion_factor;  // in cm^2
    return total_xs_in_cm2;
}


double MarleyCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & record) const {
    return 1.0 ;
}

double MarleyCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & record) const {
    return 0.0;
}

void MarleyCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord &csdr, std::shared_ptr<siren::utilities::SIREN_random> rand) const {
    std::vector<SecondaryParticleRecord> & secondary_particles = csdr.GetSecondaryParticleRecords();
    size_t lepton_index;

    if(isLepton(secondary_particles[0].type)) {
        lepton_index = 0;
    } else {
        lepton_index = 1;
    }
    SecondaryParticleRecord & lepton = secondary_particles[lepton_index];
    std::array<double,3>lepton_direction = {csdr.primary_momentum[1], csdr.primary_momentum[2], csdr.primary_momentum[3]};
    double norm = std::sqrt(lepton_direction[0]*lepton_direction[0] + lepton_direction[1]*lepton_direction[1] + lepton_direction[2]*lepton_direction[2]);
    lepton_direction[0] /= norm;
    lepton_direction[1] /= norm;
    lepton_direction[2] /= norm;
    lepton.SetDirection(lepton_direction);
    lepton.SetEnergy(csdr.primary_momentum[0]);
    double mass = 0.0;
    if(not isNeutrino(csdr.signature.primary_type)) {
        int32_t abs_pdg_code = std::abs(static_cast<int32_t>(csdr.signature.primary_type));
        switch(abs_pdg_code) {
            case 11:
                mass = siren::utilities::Constants::electronMass;
                break;
            case 13:
                mass = siren::utilities::Constants::muonMass;
                break;
            case 15:
                mass = siren::utilities::Constants::tauMass;
                break;
            default:
                throw std::runtime_error("Unknown lepton type");
        }
    }
    lepton.SetMass(mass);
    lepton.SetHelicity(1.0);

    double molar_mass = siren::detector::MaterialModel::GetMolarMass(csdr.signature.target_type);
    double target_mass = molar_mass * siren::utilities::Constants::GeV_per_amu;

    SecondaryParticleRecord & nucleus = secondary_particles[1 - lepton_index];
    std::array<double,3> nucleus_direction = {-lepton_direction[0], -lepton_direction[1], -lepton_direction[2]};
    nucleus.SetDirection(nucleus_direction);
    double nucleus_energy = csdr.primary_momentum[0] - lepton.GetEnergy();
    nucleus.SetEnergy(nucleus_energy);
    nucleus.SetMass(target_mass);
    nucleus.SetHelicity(1.0);
}



bool MarleyCrossSection::equal(MarleyCrossSection const & other) const {
    return marley_config_ == other.marley_config_;
}

std::vector<siren::dataclasses::ParticleType> MarleyCrossSection::GetPossibleTargets() const {
    return {siren::dataclasses::ParticleType::Ar40Nucleus};
}

std::vector<siren::dataclasses::ParticleType> MarleyCrossSection::GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const {
    if(isNeutrino(primary_type)) {
        return {siren::dataclasses::ParticleType::Ar40Nucleus};
    } else {
        return {};
    }
}

std::vector<siren::dataclasses::ParticleType> MarleyCrossSection::GetPossiblePrimaries() const {
    if(has_nu_cc and (not has_nubar_cc) and (not has_nc) and (not has_elastic)) {
        return {siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuTau};
    } else if(has_nubar_cc and (not has_nu_cc) and (not has_nc) and (not has_elastic) ) {
        return {siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTauBar};
    } else if(has_nu_cc or has_nubar_cc or has_nc or has_elastic) {
        return {siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::NuEBar, siren::dataclasses::ParticleType::NuMu, siren::dataclasses::ParticleType::NuMuBar, siren::dataclasses::ParticleType::NuTau, siren::dataclasses::ParticleType::NuTauBar};
    } else {
        throw std::runtime_error("No valid primary particles found");
        return {};
    }
}

std::vector<dataclasses::InteractionSignature> MarleyCrossSection::GetPossibleSignatures() const {
    std::vector<dataclasses::InteractionSignature> signatures;
    //Possible signatures are all combinations of primary and target particles

    if(has_nu_cc) {
        //CC: Nu_e + Ar40 -> e- + K40
        dataclasses::InteractionSignature cc_interaction_electron;
        cc_interaction_electron.primary_type = siren::dataclasses::ParticleType::NuE;
        cc_interaction_electron.target_type = siren::dataclasses::ParticleType::Ar40Nucleus;
        cc_interaction_electron.secondary_types = {siren::dataclasses::ParticleType::EMinus, siren::dataclasses::ParticleType::K40Nucleus};
        for(size_t i=0; i<3; ++i) {
            signatures.push_back(cc_interaction_electron);
            cc_interaction_electron.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cc_interaction_electron.primary_type) + 2);
            cc_interaction_electron.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cc_interaction_electron.secondary_types[0]) + 2),
                cc_interaction_electron.secondary_types[1]
            };
        }
    }

    if(has_nubar_cc) {
        //CC: Nu_e_bar + Ar40 -> e+ + Cl39
        dataclasses::InteractionSignature cc_interaction_positron;
        cc_interaction_positron.primary_type = siren::dataclasses::ParticleType::NuEBar;
        cc_interaction_positron.target_type = siren::dataclasses::ParticleType::Ar40Nucleus;
        cc_interaction_positron.secondary_types = {siren::dataclasses::ParticleType::EPlus, static_cast<int>(siren::dataclasses::ParticleType::1000170390)};
        for(size_t i=0; i<3; ++i) {
            signatures.push_back(cc_interaction_positron);
            cc_interaction_positron.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cc_interaction_positron.primary_type) - 2);
            cc_interaction_positron.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cc_interaction_positron.secondary_types[0]) - 2),
                cc_interaction_positron.secondary_types[1]
            };
        }
    }


    if(has_elastic) {
        //NC: Nu_e + e -> Nu_e + e
        dataclasses::InteractionSignature es_interaction;
        es_interaction.primary_type = siren::dataclasses::ParticleType::NuE;
        es_interaction.target_type = siren::dataclasses::ParticleType::EMinus;
        es_interaction.secondary_types = {siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::EMinus};
        for (size_t i=0; i<3; ++i) {
            signatures.push_back(es_interaction);
            es_interaction.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(es_interaction.primary_type) * -1);
            es_interaction.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(es_interaction.secondary_types[0]) * -1),
                es_interaction.secondary_types[1]
            };
            signatures.push_back(es_interaction);

            es_interaction.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(es_interaction.primary_type) * -1);
            es_interaction.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(es_interaction.secondary_types[0]) * -1),
                es_interaction.secondary_types[1]
            };

            es_interaction.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(es_interaction.primary_type) + 2);
            es_interaction.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(es_interaction.secondary_types[0]) + 2),
                es_interaction.secondary_types[1]
            };
        }
    }

    if(has_nc) {
        //CEvNS: Nu_e + Ar40 -> Nu_e + Ar40
        dataclasses::InteractionSignature cevns_interaction;
        cevns_interaction.primary_type = siren::dataclasses::ParticleType::NuE;
        cevns_interaction.target_type = siren::dataclasses::ParticleType::Ar40Nucleus;
        cevns_interaction.secondary_types = {siren::dataclasses::ParticleType::NuE, siren::dataclasses::ParticleType::Ar40Nucleus};
        for (size_t i=0; i<3; ++i) {
            signatures.push_back(cevns_interaction);
            cevns_interaction.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cevns_interaction.primary_type) * -1);
            cevns_interaction.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cevns_interaction.secondary_types[0]) * -1),
                cevns_interaction.secondary_types[1]
            };
            signatures.push_back(cevns_interaction);

            cevns_interaction.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cevns_interaction.primary_type) * -1);
            cevns_interaction.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cevns_interaction.secondary_types[0]) * -1),
                cevns_interaction.secondary_types[1]
            };

            cevns_interaction.primary_type = static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cevns_interaction.primary_type) + 2);
            cevns_interaction.secondary_types = {
                static_cast<siren::dataclasses::ParticleType>(static_cast<int32_t>(cevns_interaction.secondary_types[0]) + 2),
                cevns_interaction.secondary_types[1]
            };
        }

    }

    return signatures;
}

std::vector<dataclasses::InteractionSignature> MarleyCrossSection::GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const {
    std::vector<dataclasses::InteractionSignature> signatures;
    for(dataclasses::InteractionSignature const & signature : GetPossibleSignatures()) {
        if(signature.primary_type == primary_type and signature.target_type == target_type) {
            signatures.push_back(signature);
        }
    }
    return signatures;
}

double MarleyCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    return 1.0;
}

std::vector<std::string> MarleyCrossSection::DensityVariables() const {
    return {};
}


} // namespace interactions
} // namespace siren

