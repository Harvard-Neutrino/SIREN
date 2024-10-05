#include <filesystem>  // for path, create_directories

#include "SIREN/interactions/MarleyCrossSection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/detector/MaterialModel.h"

#include "marley/Generator.hh"  // MARLEY generator
#include "marley/JSON.hh"  // For handling MARLEY JSON configuration
#include "marley/JSONConfig.hh"

namespace siren {
namespace interactions {

MarleyCrossSection::MarleyCrossSection(std::string marley_react_file, std::string marley_nuclide_index_file, std::string marley_nuclide_file, std::string marley_masses_file, std::string marley_gs_parity_file) {
    std::filesystem::path tmp_dir_path {std::filesystem::temp_directory_path() /= std::tmpnam(nullptr)};
    std::filesystem::create_directories(tmp_dir_path);

    marley_react_fname_ = marley_react_file;
    marley_nuclide_index_fname_ = marley_nuclide_index_file;
    marley_nuclide_fname_ = marley_nuclide_file;
    marley_masses_fname_ = marley_masses_file;
    marley_gs_parity_fname_ = marley_gs_parity_file;

    std::vector<std::tuple<std::string, std::vector<char> *> > files = {
        {marley_react_file, &marley_react_data_},
        {marley_nuclide_index_file, &marley_nuclide_index_data_},
        {marley_nuclide_file, &marley_nuclide_data_},
        {marley_masses_file, &marley_masses_data_},
        {marley_gs_parity_file, &marley_gs_parity_data_}
    };

    for(auto const & file : files) {
        std::string filename = std::get<0>(file);
        std::vector<char> & data = *std::get<1>(file);

        std::string basename = std::filesystem::path(filename).filename();
        std::string dest_file = (tmp_dir_path / basename).string();

        std::ifstream ifs(filename, std::ios::binary);
        std::stringstream ss;
        ss << ifs.rdbuf();
        ifs.close();
        data.resize(ss.str().size());
        std::string string_data = ss.str();
        std::copy(string_data.begin(), string_data.end(), data.begin());

        std::ofstream ofs(dest_file, std::ios::binary);
        ofs << ss.str();
    }

    std::string search_path = tmp_dir_path.string();
    std::string react_file = (tmp_dir_path / std::filesystem::path(marley_react_file).filename()).string();

    setenv("MARLEY", "", 0);
    setenv("MARLEY_SEARCH_PATH", search_path.c_str(), 0);
    marley::FileManager::Instance();
    siren::interactions::marley_::FileManager_::set_search_path(search_path);
    InitializeMarley(react_file);
}

MarleyCrossSection::MarleyCrossSection(std::array<std::vector<char>, 5> const & data, std::array<std::string, 5> const & fnames) {
    marley_react_data_ = data[0];
    marley_nuclide_index_data_ = data[1];
    marley_nuclide_data_ = data[2];
    marley_masses_data_ = data[3];
    marley_gs_parity_data_ = data[4];
    marley_react_fname_ = fnames[0];
    marley_nuclide_index_fname_ = fnames[1];
    marley_nuclide_fname_ = fnames[2];
    marley_masses_fname_ = fnames[3];
    marley_gs_parity_fname_ = fnames[4];

    std::filesystem::path tmp_dir_path {std::filesystem::temp_directory_path() /= std::tmpnam(nullptr)};
    std::filesystem::create_directories(tmp_dir_path);

    for(size_t i=0; i<5; ++i) {
        std::vector<char> const & d = data[i];
        std::string const & filename = fnames[i];
        std::string dest_file = (tmp_dir_path / filename).string();
        std::string string_data(d.begin(), d.end());
        std::ofstream ofs(dest_file, std::ios::binary);
        ofs << string_data;
    }

    std::string search_path = tmp_dir_path.string();
    std::string react_file = (tmp_dir_path / marley_react_fname_).string();

    setenv("MARLEY", "", 0);
    setenv("MARLEY_SEARCH_PATH", search_path.c_str(), 0);
    marley::FileManager::Instance();
    siren::interactions::marley_::FileManager_::set_search_path(search_path);
    InitializeMarley(react_file);
}

void MarleyCrossSection::InitializeMarley(std::string const & marley_react_file) {
    structure_database_ = std::make_unique<marley::StructureDatabase>();
    reactions_ = marley::Reaction::load_from_file(marley_react_file, *structure_database_);

    std::vector<std::unique_ptr<marley::Reaction>> const & reactions = reactions_;

    has_nu_cc = false;
    has_nubar_cc = false;
    has_nc = false;
    has_elastic = false;

    for(std::unique_ptr<marley::Reaction> const & reaction : reactions) {
        marley::Reaction::ProcessType process = reaction->process_type();
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

    double energy = record.primary_momentum[0];  // Energy of the primary particle
    double mass = record.primary_mass;  // Mass of the primary particle
    int pdg_a = static_cast<int32_t>(record.signature.primary_type); // PDG code of the primary particle
    double KEa = std::sqrt(energy*energy - mass*mass);  // Kinetic energy of the primary particle
    int pdg_atom = static_cast<int32_t>(record.signature.target_type);  // PDG code of the target atom
    double KEa_MeV = KEa * 1e3;  // Convert kinetic energy to MeV

    std::vector<std::unique_ptr<marley::Reaction>> const & reactions = reactions_;
    std::vector<marley::Reaction const *> the_reactions;

    marley::Reaction::ProcessType desired_process_type;

    size_t lepton_index;

    if(isLepton(record.signature.secondary_types[0])) {
        lepton_index = 0;
    } else {
        lepton_index = 1;
    }

    siren::dataclasses::ParticleType lepton_type = record.signature.secondary_types[lepton_index];

    if(isNeutrino(lepton_type)) {
        // NC or elastic
        siren::dataclasses::ParticleType target_type = record.signature.target_type;
        if(target_type == siren::dataclasses::ParticleType::EMinus) {
            desired_process_type = marley::Reaction::ProcessType::NuElectronElastic;
        } else {
            desired_process_type = marley::Reaction::ProcessType::NC;
        }
    } else {
        // CC
        siren::dataclasses::ParticleType primary_type = record.signature.primary_type;
        int32_t primary_pdg_code = static_cast<int32_t>(primary_type);
        if(primary_pdg_code > 0) {
            desired_process_type = marley::Reaction::ProcessType::NeutrinoCC;
        } else {
            desired_process_type = marley::Reaction::ProcessType::AntiNeutrinoCC;
        }
    }

    for(std::unique_ptr<marley::Reaction> const & reaction : reactions) {
        marley::Reaction::ProcessType process = reaction->process_type();
        // Skip reactions which involve a different target atom
        if(pdg_atom != reaction->atomic_target().pdg())
            continue;
        // Skip reactions which involve a different projectile
        if(pdg_a != reaction->pdg_a())
            continue;
        // Skip reactions which involve a different final state (determined by the process type)
        if(process != desired_process_type)
            continue;
        the_reactions.push_back(reaction.get());
    }

    double total_xs_in_MeV2 = 0.0;

    for(marley::Reaction const * reaction : the_reactions) {
        total_xs_in_MeV2 += reaction->total_xs(pdg_a, KEa_MeV); //in MeV^(-2)
    }

    double conversion_factor = 3.893793721718594e-22;  // MeV^(-2) to cm^2
    double total_xs_in_cm2 = total_xs_in_MeV2 * conversion_factor;  // in cm^2
    return total_xs_in_cm2;
}

double MarleyCrossSection::TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {

    double energy = record.primary_momentum[0];  // Energy of the primary particle
    double mass = record.primary_mass;  // Mass of the primary particle
    int pdg_a = static_cast<int32_t>(record.signature.primary_type); // PDG code of the primary particle
    double KEa = std::sqrt(energy*energy - mass*mass);  // Kinetic energy of the primary particle
    int pdg_atom = static_cast<int32_t>(record.signature.target_type);  // PDG code of the target atom
    double KEa_MeV = KEa * 1e3;  // Convert kinetic energy to MeV

    std::vector<std::unique_ptr<marley::Reaction>> const & reactions = reactions_;
    std::vector<marley::Reaction const *> the_reactions;

    for(std::unique_ptr<marley::Reaction> const & reaction : reactions) {
        // Skip reactions which involve a different target atom
        if(pdg_atom != reaction->atomic_target().pdg())
            continue;
        // Skip reactions which involve a different projectile
        if(pdg_a != reaction->pdg_a())
            continue;
        the_reactions.push_back(reaction.get());
    }
    double total_xs_in_MeV2 = 0.0;


    for(marley::Reaction const * reaction : the_reactions) {
        total_xs_in_MeV2 += reaction->total_xs(pdg_a, KEa_MeV); //in MeV^(-2)
    }

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
    std::vector<siren::dataclasses::SecondaryParticleRecord> & secondary_particles = csdr.GetSecondaryParticleRecords();
    size_t lepton_index;

    if(isLepton(secondary_particles[0].type)) {
        lepton_index = 0;
    } else {
        lepton_index = 1;
    }
    siren::dataclasses::SecondaryParticleRecord & lepton = secondary_particles[lepton_index];
    std::array<double,3>lepton_direction = {csdr.primary_momentum[1], csdr.primary_momentum[2], csdr.primary_momentum[3]};
    double norm = std::sqrt(lepton_direction[0]*lepton_direction[0] + lepton_direction[1]*lepton_direction[1] + lepton_direction[2]*lepton_direction[2]);
    lepton_direction[0] /= norm;
    lepton_direction[1] /= norm;
    lepton_direction[2] /= norm;
    lepton.SetDirection(lepton_direction);
    lepton.SetEnergy(csdr.primary_momentum[0]);
    double mass = 0.0;
    if(not isNeutrino(lepton.type)) {
        int32_t abs_pdg_code = std::abs(static_cast<int32_t>(lepton.type));
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

    siren::dataclasses::SecondaryParticleRecord & nucleus = secondary_particles[1 - lepton_index];
    std::array<double,3> nucleus_direction = {-lepton_direction[0], -lepton_direction[1], -lepton_direction[2]};
    nucleus.SetDirection(nucleus_direction);
    double nucleus_energy = csdr.primary_momentum[0] - lepton.GetEnergy();
    nucleus.SetEnergy(nucleus_energy);
    nucleus.SetMass(target_mass);
    nucleus.SetHelicity(1.0);
}



bool MarleyCrossSection::equal(CrossSection const & other) const {
    const MarleyCrossSection* x = dynamic_cast<const MarleyCrossSection*>(&other);

    if(!x)
        return false;
    else
        return std::tie(marley_react_data_, marley_nuclide_index_data_, marley_nuclide_data_, marley_masses_data_, marley_gs_parity_data_, marley_react_fname_, marley_nuclide_index_fname_, marley_nuclide_fname_, marley_masses_fname_, marley_gs_parity_fname_) == std::tie(x->marley_react_data_, x->marley_nuclide_index_data_, x->marley_nuclide_data_, x->marley_masses_data_, x->marley_gs_parity_data_, x->marley_react_fname_, x->marley_nuclide_index_fname_, x->marley_nuclide_fname_, x->marley_masses_fname_, x->marley_gs_parity_fname_);
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
        cc_interaction_positron.secondary_types = {siren::dataclasses::ParticleType::EPlus, static_cast<siren::dataclasses::ParticleType>(1000170390)};
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

