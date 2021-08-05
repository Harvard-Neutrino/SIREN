#include "LeptonInjector/CrossSection.h"

#include <array>

namespace LeptonInjector {

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types) primary_types_(primary_types), target_types(target_types), minimum_Q2_(minimum_Q2), target_mass_(target_mass), interaction_type_(interaction) {
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types) primary_types_(primary_types), target_types(target_types) {
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types) primary_types_(primary_types.begin(), primary_types.end()), target_types(target_types.begin(), target_types.end()), minimum_Q2_(minimum_Q2), target_mass_(target_mass), interaction_type_(interaction) {
    LoadFromFile(differential_filename, total_filename);
    InitializeSignatures();
}

DISFromSpline::DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types) primary_types_(primary_types.begin(), primary_types.end()), target_types(target_types.begin(), target_types.end()) {
    LoadFromFile(differential_filename, total_filename);
    ReadParamsFromSplineTable();
    InitializeSignatures();
}

void DISFromSpline::LoadFromFile(std::string dd_crossSectionFile, std::string total_crossSectionFile) {

    crossSection = photospline::splinetable<>(dd_crossSectionFile.c_str());

    if(differential_cross_section_.get_ndim()!=3 && differential_cross_section_.get_ndim()!=2)
        throw("cross section spline has " + std::to_string(differential_cross_section_.get_ndim())
                + " dimensions, should have either 3 (log10(E), log10(x), log10(y)) or 2 (log10(E), log10(y))");

    total_cross_section_ = photospline::splinetable<>(total_crossSectionFile.c_str());

    if(total_cross_section_.get_ndim() != 1)
        throw("Total cross section spline has " + std::to_string(total_cross_section_.get_ndim())
                + " dimensions, should have 1, log10(E)");
}

void DISFromSpline::ReadParamsFromSplineTable() {
    // returns true if successfully read target mass
    bool mass_good = differential_cross_section_.read_key("TARGETMASS", target_mass_);
    // returns true if successfully read interaction type
    bool int_good = differential_cross_section_.read_key("INTERACTION", interaction_type_);
    // returns true if successfully read minimum Q2
    bool q2_good = differential_cross_section_.read_key("Q2MIN", minimum_Q2_);

    if(!int_good) {
        // assume DIS to preserve compatability with previous versions
        interaction = 1;
    }

    if(!q2_good) {
        // assume 1 GeV^2
        Q2Min = 1;
    }

    if(!mass_good) {
        if(int_good) {
            if(interaction_type_ == 1 or interaction_type_ == 2) {
                target_mass_ = (LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::PPlus)+
                        LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::Neutron))/2;
            } else if(interaction_type_ == 3) {
                target_mass_ = LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::EMinus);
            } else {
                throw("Logic error. Interaction type is not 1, 2, or 3!");
            }

        } else {
            if(differential_cross_section_.get_ndim() == 3) {
                target_mass_ = (LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::PPlus)+
                        LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::Neutron))/2;
            } else if(crossSection.get_ndim() == 2) {
                target_mass_ = LeptonInjector::particleMass(LeptonInjector::Particle::ParticleType::EMinus);
            } else {
                throw("Logic error. Spline dimensionality is not 2, or 3!");
            }
        }
    }
}

DISFromSpline::InitializeSignatures() {
    signatures_.clear();
    for(auto primary_type : primary_types_) {
        InteractionSignature signature;
        signature.primary_type = primary_type;

        if(not isNeutrino(primary_type)) {
            throw std::runtime_error("This DIS implementation only supports neutrinos as primaries!");
        }

        Particle::ParticleType charged_lepton_product = Particle::ParticleType::unknown;
        Particle::ParticleType neutral_lepton_product = primary_type;

        if(primary_type == Particle::ParticleType::NuE) {
            charged_lepton_product = Particle::ParticleType::EMinus;
        } else if(primary_type == Particle::ParticleType::NuEBar) {
            charged_lepton_product = Particle::ParticleType::EPlus;
        } else if(primary_type == Particle::ParticleType::NuMu) {
            charged_lepton_product = Particle::ParticleType::MuMinus;
        } else if(primary_type == Particle::ParticleType::NuMuBar) {
            charged_lepton_product = Particle::ParticleType::MuPlus;
        } else if(primary_type == Particle::ParticleType::NuTau) {
            charged_lepton_product = Particle::ParticleType::TauMinus;
        } else if(primary_type == Particle::ParticleType::NuTauBar) {
            charged_lepton_product = Particle::ParticleType::TauPlus;
        } else {
            throw std::runtime_error("InitializeSignatures: Unkown parent neutrino type!");
        }

        if(interaction_type_ == 1) {
            secondary_types.append(charged_lepton_product);
        } else if(interaction_type == 2) {
            secondary_types.append(neutral_lepton_product);
        } else if(interaction_type == 3) {
            secondary_types.append(Particle::ParticleType::Hadrons);
        } else {
            throw std::runtime_error("InitializeSignatures: Unkown interaction type!");
        }
        for(auto target_type : target_types_) {
            std::pair<Particle::ParticleType, Particle::ParticleType> key(primary_type, target_type);
            signature.target_type = target_type;
            signatures.emplace(key, signature);
        }
    }
}

DISFromSpline::TotalCrossSection(InteractionRecord const & interaction) const {
    LeptonInjector::Particle::ParticleType primary_type = interaction.primary_type;
    double primary_energy;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
    } else {
        throw std::runtime_error("Lorentz boost not implemented!");
    }
    return TotalCrossSection(primary_type, primary_energy);
}

DISFromSpline::TotalCrossSection(LeptonInjector::Particle::ParticleType primary_type, double primary_energy) const {
    if(not primary_types.contains(primary_type)) {
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    double target_mass = target_mass_;
    double log_energy = log10(primary_energy);

    if(log_energy < totalCrossSection.lower_extent(0)
            or log_energy > totalCrossSection.upper_extent(0)) {
        throw("Interaction energy ("+ std::to_string(energy) +
                ") out of cross section table range: ["
                + std::to_string(pow(10.,totalCrossSection.lower_extent(0))) + " GeV,"
                + std::to_string(pow(10.,totalCrossSection.upper_extent(0))) + " GeV]");
    }

    int center;
    total_cross_section_.searchcenters(&log_energy, &center);
    double log_xs = total_cross_section_.ndsplineeval(&log_energy, &center, 0);

    return std::pow(10.0, log_xs);
}


} // namespace LeptonInjector

