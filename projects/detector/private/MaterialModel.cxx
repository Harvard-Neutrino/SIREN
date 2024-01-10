#include "LeptonInjector/detector/MaterialModel.h"

#include <map>                                    // for map, operator!=
#include <cmath>                                  // for pow, exp, log, sqrt
#include <string>                                 // for basic_string, opera...
#include <vector>                                 // for vector, operator==
#include <fstream>                                // for basic_istream, getline
#include <sstream>
#include <stdio.h>                                // for snprintf, sscanf
#include <stdlib.h>                               // for abs

#include "LeptonInjector/dataclasses/Particle.h"  // for Particle
#include "LeptonInjector/utilities/Constants.h"   // for avogadro, GeV_per_amu

using namespace LI::detector ;

MaterialModel::Component::Component(LI::dataclasses::Particle::ParticleType type)
    : type(type)
{
    if(type == LI::dataclasses::Particle::ParticleType::PPlus) {
        neutron_count = 0;
        proton_count = 1;
        nucleon_count = 1;
        is_atom = false;
        molar_mass = atomic_masses.at({strange_count, neutron_count, proton_count, nucleon_count});
    } else if (type == LI::dataclasses::Particle::ParticleType::Neutron) {
        neutron_count = 1;
        proton_count = 0;
        nucleon_count = 1;
        is_atom = false;
        molar_mass = atomic_masses.at({strange_count, neutron_count, proton_count, nucleon_count});
    } else if (type == LI::dataclasses::Particle::ParticleType::Nucleon) {
        neutron_count = 0;
        proton_count = 0;
        nucleon_count = 1;
        is_atom = false;
        molar_mass = (atomic_masses.at({0, 1, 0, 1}) + atomic_masses.at({0, 0, 1, 1})) / 2.0;
    } else if (type == LI::dataclasses::Particle::ParticleType::EMinus) {
        neutron_count = 0;
        proton_count = 0;
        nucleon_count = 0;
        is_atom = false;
        molar_mass = electron_molar_mass;
    } else {
        try {
            GetNucleonContent(static_cast<int>(type), strange_count, neutron_count, proton_count, nucleon_count);
            is_atom = true;
            std::tuple<int, int, int, int> key(strange_count, neutron_count, proton_count, nucleon_count);
            auto it = atomic_masses.find(key);
            if(it != atomic_masses.end())
                molar_mass = it->second;
            else {
                double binding_energy = GetEmpericalNuclearBindingEnergy(strange_count, neutron_count, proton_count, nucleon_count);
                molar_mass = strange_count * LI::utilities::Constants::lambda0Mass / LI::utilities::Constants::GeV_per_amu
                    + neutron_count * atomic_masses.at({0, 1, 0, 1})
                    + proton_count * atomic_masses.at({0, 0, 1, 1})
                    - binding_energy / LI::utilities::Constants::GeV_per_amu;
            }
        } catch (std::runtime_error const & e) {
            strange_count = 0;
            neutron_count = 0;
            proton_count = 0;
            nucleon_count = 0;
        }
    }
}

bool MaterialModel::Component::operator==(MaterialModel::Component const & component) const {
    return type == component.type;
}

bool MaterialModel::MaterialComponent::operator==(MaterialModel::MaterialComponent const & other) const {
    return
        std::tie(
                component,
                mass_density_over_total_mass_density,
                particle_density_over_total_mass_density)
        ==
        std::tie(
                other.component,
                other.mass_density_over_total_mass_density,
                other.particle_density_over_total_mass_density);
}

bool MaterialModel::operator==(MaterialModel const & other) const {
    return material_components_ == other.material_components_;
}

MaterialModel::MaterialModel() {}

MaterialModel::MaterialModel(std::string const & file) {
    AddModelFile(file);
}

MaterialModel::MaterialModel(std::string const & path, std::string const & file) : path_(path) {
    AddModelFile(file);
}

MaterialModel::MaterialModel(std::vector<std::string> const & files) {
    AddModelFiles(files);
}

MaterialModel::MaterialModel(std::string const & path, std::vector<std::string> const & files) : path_(path) {
    AddModelFiles(files);
}

void MaterialModel::SetPath(std::string const & path) {
    path_ = path;
}

void MaterialModel::AddMaterial(std::string const & material_name, std::map<int, double> const & component_mass_fractions) {
    int material_id;
    bool new_material = material_ids_.find(material_name) == material_ids_.end();
    if(new_material) {
        material_id = material_names_.size();
        material_components_.resize(material_components_.size() + 1);
        material_ids_.insert({material_name, material_id});
        material_names_.push_back(material_name);
    } else {
        material_id = material_ids_[material_name];
    }

    double neutrons_per_amu = 0;
    double nucleons_per_amu = 0;
    double protons_per_amu = 0;

    std::vector<MaterialComponent> material_components;

    for(auto const & mass_frac_it : component_mass_fractions) {
        int component_id = mass_frac_it.first;
        double component_mass_fraction = mass_frac_it.second;
        Component component(static_cast<LI::dataclasses::Particle::ParticleType>(component_id));
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = component_mass_fraction;
        material_component.particle_density_over_total_mass_density = LI::utilities::Constants::avogadro * component_mass_fraction / material_component.component.molar_mass;
        material_components.push_back(material_component);
        neutrons_per_amu += component_mass_fraction / material_component.component.molar_mass * material_component.component.neutron_count;
        nucleons_per_amu += component_mass_fraction / material_component.component.molar_mass * material_component.component.nucleon_count;
        protons_per_amu += component_mass_fraction / material_component.component.molar_mass * material_component.component.proton_count;
    }

    if(neutrons_per_amu > 0) {
        Component component(LI::dataclasses::Particle::ParticleType::Neutron);
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = neutrons_per_amu * component.molar_mass;
        material_component.particle_density_over_total_mass_density = neutrons_per_amu * LI::utilities::Constants::avogadro;
        material_components.push_back(material_component);
    }

    if(nucleons_per_amu > 0) {
        Component component(LI::dataclasses::Particle::ParticleType::Nucleon);
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = nucleons_per_amu * component.molar_mass;
        material_component.particle_density_over_total_mass_density = nucleons_per_amu * LI::utilities::Constants::avogadro;
        material_components.push_back(material_component);
    }

    if(protons_per_amu > 0) {
        Component component(LI::dataclasses::Particle::ParticleType::PPlus);
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = protons_per_amu * component.molar_mass;
        material_component.particle_density_over_total_mass_density = protons_per_amu * LI::utilities::Constants::avogadro;
        material_components.push_back(material_component);
        component = Component(LI::dataclasses::Particle::ParticleType::EMinus);
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = protons_per_amu * component.molar_mass;
        material_component.particle_density_over_total_mass_density = protons_per_amu * LI::utilities::Constants::avogadro;
        material_components.push_back(material_component);
    }

    // Store material component list for material
    material_components_[material_id] = material_components;

    // Compute radiation length
    double rad_length = ComputeMaterialRadiationLength(material_id);
    if(new_material) {
        material_radiation_length_.push_back(rad_length);
        // Store material components by material id and component particle type
        for(auto const & material_component : material_components) {
            material_components_by_id_.insert({std::make_pair(material_id, material_component.component.type), material_component});
        }
    } else {
        material_radiation_length_[material_id] = rad_length;
        // Store material components by material id and component particle type
        for(auto const & material_component : material_components) {
            material_components_by_id_[std::make_pair(material_id, material_component.component.type)] = material_component;
        }
    }
}

void MaterialModel::AddModelFiles(std::vector<std::string> const & matratios) {
    for(auto matratio : matratios)
        AddModelFile(matratio);
}

namespace {
    bool fexists(const std::string filename)
    {
        std::ifstream ifile(filename.c_str());
        return (bool)ifile;
    }
}

void MaterialModel::AddModelFile(std::string matratio) {
    std::string fname;

    if(matratio.empty())
        throw(std::runtime_error("Received empty matratio filename!"));

    if(fexists(matratio)) {
        fname = matratio;
    }
    else if(fexists(matratio + ".dat")) {
        fname = matratio + ".dat";
    }
    else if(fexists(path_ + "/materials/" + matratio)) {
        fname = path_ + "/materials/" + matratio;
    }
    else if(fexists(path_ + "/materials/" + matratio + ".dat")) {
        fname = path_ + "/materials/" + matratio + ".dat";
    }
    else if(fexists(path_ + "/" + matratio)) {
        fname = path_ + "/" + matratio;
    }
    else if(fexists(path_ + "/" + matratio + ".dat")) {
        fname = path_ + "/" + matratio + ".dat";
    }
    else {
        throw(std::runtime_error("Cannot open matratio file!"));
    }

    // check earthmodel file
    std::ifstream in(fname.c_str(), std::ifstream::in);

    if (in.fail())
        throw(std::runtime_error("Failed to open " + fname + ". Set correct path."));

    // read the file
    std::string buffer;
    std::string material_name;
    int component_pdg_code, num_components;
    double mass_fraction;
    int nread = 0;

    while(not (in.eof() or in.fail())) {

        std::getline(in, buffer);
        nread = buffer.size();

        if(nread == 0 || buffer[0] == ' ' || buffer[0] == '\t' || buffer[0] == '#') {
            // new line, start from white space, or comment line.
            continue;
        } else {
            // material density data
            std::stringstream ss(buffer);

            // read material name and number of components
            ss >> material_name >> num_components;

            std::map<int, double> component_mass_fractions;
            for(int i=0; i<num_components; ++i) {
                std::getline(in, buffer);
                nread = buffer.size();
                if(nread == 0 || buffer[0] == ' ' || buffer[0] == '\t' || buffer[0] == '#') {
                    // new line, start from white space, or comment line.
                    --i;
                    continue;
                } else {
                    std::stringstream ss2(buffer);
                    ss2 >> component_pdg_code >> mass_fraction;
                    component_mass_fractions[component_pdg_code] = mass_fraction;
                }
            }
            AddMaterial(material_name, component_mass_fractions);
        }

    } // end of the while loop
    model_files_.push_back(matratio);

    in.close();
}

std::vector<LI::dataclasses::Particle::ParticleType> MaterialModel::GetMaterialTargets(int material_id) const {
    std::vector<LI::dataclasses::Particle::ParticleType> targets;
    std::vector<MaterialComponent> const & components = material_components_[material_id];
    targets.reserve(components.size());
    for(auto const & comp : components) {
        targets.push_back(comp.component.type);
    }
    return targets;
}

double MaterialModel::GetMaterialRadiationLength(int id) const {
    return material_radiation_length_.at(id);
}


std::string MaterialModel::GetMaterialName(int material_id) const {
    return material_names_.at(material_id);
}

int MaterialModel::GetMaterialId(std::string const & material_name) const {
    return material_ids_.at(material_name);
}

bool MaterialModel::HasMaterial(std::string const & name) const {
    return material_ids_.count(name) > 0;
}

bool MaterialModel::HasMaterial(int id) const {
    if(id < 0) {
        return false;
    }
    return material_names_.size() > (unsigned int)(id);
}

double MaterialModel::GetTargetMassFraction(int material_id, LI::dataclasses::Particle::ParticleType particle_type) const {
    std::pair<int, LI::dataclasses::Particle::ParticleType> key(material_id, particle_type);
    if(material_components_by_id_.find(key) != material_components_by_id_.end())
        return material_components_by_id_.at(key).mass_density_over_total_mass_density;
    else
        return 0.0;
}

double MaterialModel::GetTargetParticleFraction(int material_id, LI::dataclasses::Particle::ParticleType particle_type) const {
    std::pair<int, LI::dataclasses::Particle::ParticleType> key(material_id, particle_type);
    if(material_components_by_id_.find(key) != material_components_by_id_.end())
        return material_components_by_id_.at(key).particle_density_over_total_mass_density;
    else
        return 0.0;
}

std::vector<double> MaterialModel::GetTargetMassFraction(int material_id, std::vector<LI::dataclasses::Particle::ParticleType> const & particle_types) const {
    std::vector<double> fractions;
    fractions.reserve(particle_types.size());
    for(auto const & particle_type : particle_types) {
        std::pair<int, LI::dataclasses::Particle::ParticleType> key(material_id, particle_type);
        if(material_components_by_id_.find(key) != material_components_by_id_.end())
            fractions.push_back(material_components_by_id_.at(key).mass_density_over_total_mass_density);
        else
            fractions.push_back(0.0);
    }
    return fractions;
}

std::vector<double> MaterialModel::GetTargetParticleFraction(int material_id, std::vector<LI::dataclasses::Particle::ParticleType> const & particle_types) const {
    std::vector<double> fractions;
    fractions.reserve(particle_types.size());
    for(auto const & particle_type : particle_types) {
        std::pair<int, LI::dataclasses::Particle::ParticleType> key(material_id, particle_type);
        if(material_components_by_id_.find(key) != material_components_by_id_.end())
            fractions.push_back(material_components_by_id_.at(key).particle_density_over_total_mass_density);
        else
            fractions.push_back(0.0);
    }
    return fractions;
}

double MaterialModel::ComputeMaterialRadiationLength(int material_id) const {
    // This function calculates the radiation length of a given material in g/cm^2
    // Takes screening effects into account
    // Averages over constituent materials in a composite
    // See Page 21 of Particle Detectors by Grupen and Shwartz
    double X0inv = 0;
    for (auto const & component : material_components_[material_id]) {
        if(not component.component.is_atom)
            continue;
        int A = component.component.nucleon_count;
        int Z = component.component.proton_count;
        double f = component.mass_density_over_total_mass_density;
        double X0i = 716.4 * A / ( Z*(Z + 1) * std::log(287./std::sqrt(Z))); // g/cm^2, Grupen eq 1.59
        X0inv += f / X0i;
    }
    return 1.0 / X0inv;
}

std::vector<double> MaterialModel::GetTargetRadiationFraction(int material_id, std::vector<LI::dataclasses::Particle::ParticleType> const & particle_types) const {
    double X0inv = 0;
    std::vector<double> fractions;
    fractions.reserve(particle_types.size());
    for(auto const & particle_type : particle_types) {
        std::pair<int, LI::dataclasses::Particle::ParticleType> key(material_id, particle_type);
        if(material_components_by_id_.find(key) != material_components_by_id_.end()) {
            fractions.push_back(0.0);
            continue;
        }
        MaterialComponent const & component = material_components_by_id_.at(key);
        if(not component.component.is_atom) {
            fractions.push_back(0.0);
            continue;
        }
        int A = component.component.nucleon_count;
        int Z = component.component.proton_count;
        double f = component.mass_density_over_total_mass_density;
        double X0i = 716.4 * A / ( Z*(Z + 1) * std::log(287./std::sqrt(Z))); // g/cm^2, Grupen eq 1.59
        double frac = f / X0i;
        fractions.push_back(frac);
        X0inv += frac;
    }
    for(unsigned int i=0; i<fractions.size(); ++i) {
        fractions[i] /= X0inv;
    }
    return fractions;
}

int MaterialModel::GetNucleonContent(int code, int & strange_count, int & neutron_count, int & proton_count, int & nucleon_count) {
    int prefix = 0;
    int excitation = 0;

    char buf[CHAR_BUF_SIZE];
    snprintf(buf, CHAR_BUF_SIZE, "%d", code);
    int nread = sscanf(buf, "%2d%1d%3d%3d%1d", &prefix, &strange_count, &proton_count, &nucleon_count, &excitation);
    if (nread != 5) {
        throw std::runtime_error("Failed to convert nuclear pdg to 10LZZZAAAI "
                "prefix "+std::to_string(prefix)+", L "+std::to_string(strange_count)+", Z "+std::to_string(proton_count)+", A "+std::to_string(nucleon_count)+", I "+std::to_string(excitation));
    }
    neutron_count = nucleon_count - proton_count - strange_count;
    return 0;
}

std::vector<LI::dataclasses::Particle::ParticleType> MaterialModel::GetMaterialConstituents(int material_id) const {
    std::vector<LI::dataclasses::Particle::ParticleType> particles;;
    particles.reserve(material_components_[material_id].size());
    for (auto const & component : material_components_[material_id]) {
        particles.push_back(component.component.type);
    }
    return particles;
}

double MaterialModel::GetMolarMass(LI::dataclasses::Particle::ParticleType particle) {
    return Component(particle).molar_mass;
}

int MaterialModel::GetStrangeCount(LI::dataclasses::Particle::ParticleType particle) {
    return Component(particle).strange_count;
}

int MaterialModel::GetNucleonCount(LI::dataclasses::Particle::ParticleType particle) {
    return Component(particle).nucleon_count;
}

int MaterialModel::GetNeutronCount(LI::dataclasses::Particle::ParticleType particle) {
    return Component(particle).neutron_count;
}

int MaterialModel::GetProtonCount(LI::dataclasses::Particle::ParticleType particle) {
    return Component(particle).proton_count;
}

int MaterialModel::GetEmpericalNuclearBindingEnergy(int strange_count, int neutron_count, int proton_count, int nucleon_count) {
    // Generalized mass formula and parameters from https://arxiv.org/abs/nucl-th/0504085
    // Nucleus mass formula parameters comes from least squares fit to experimental data
    // Hypernucleus correction parameters in the paper result from a two parameter fit of c1 and c2 while keeping c0 fixed
    // PDG particle codes assume the strange contribution to the hypernucleus comes from lambdas
    // Here we assume the strange contribution comes from lambda0
    constexpr const double a_nu = 15.777; // MeV
    constexpr const double a_s = 18.34; // MeV
    constexpr const double a_c = 0.71; // MeV
    constexpr const double a_sym = 23.21; // MeV
    constexpr const double k = 17;
    constexpr const double c = 30;
    const double lambda0_mass = LI::utilities::Constants::lambda0Mass * 1e3; // GeV --> MeV

    constexpr const double c0 = 0.0335;
    constexpr const double c1 = 26.7;
    constexpr const double c2 = 48.7;
    constexpr const double S = -1; // Strangeness of lambda


    double N = neutron_count;
    double Z = proton_count;
    double A = nucleon_count;
    double L = strange_count;

    double delta = 12 * std::pow(A, -0.5);
    if(proton_count % 2 == 0 and neutron_count % 2 == 0) {
    } else if (proton_count % 2 == 1 and neutron_count % 2 == 1) {
        delta = -delta;
    } else {
        delta = 0;
    }

    double delta_new = (1.0 - std::exp(-A/c))*delta;

    double binding_energy_in_MeV =
        a_nu * A
        - a_s * std::pow(A, 2.0/3.0)
        - a_c * Z * (Z - 1) / std::pow(A, 1.0/3.0)
        - a_sym * std::pow(N - Z, 2) / ((1 + std::exp(-A/k)) * A)
        + delta_new
        + L * (c0 * lambda0_mass - c1 - c2 * std::abs(S) / std::pow(A, 2.0/3.0));
    return binding_energy_in_MeV * 1e-3; // GeV
}
