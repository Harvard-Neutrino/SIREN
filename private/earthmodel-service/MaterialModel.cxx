#include <map>
#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "LeptonInjector/Particle.h"
#include "LeptonInjector/Constants.h"

#include "earthmodel-service/MaterialModel.h"

using namespace earthmodel;

MaterialModel::Component::Component(LeptonInjector::Particle::ParticleType type)
    : type(type)
{
    if(type == LeptonInjector::Particle::ParticleType::PPlus) {
        neutron_count = 0;
        proton_count = 1;
        nucleon_count = 1;
        is_atom = false;
        molar_mass = atomic_masses.at({neutron_count, proton_count, nucleon_count});
    } else if (type == LeptonInjector::Particle::ParticleType::Neutron) {
        neutron_count = 1;
        proton_count = 0;
        nucleon_count = 1;
        is_atom = false;
        molar_mass = atomic_masses.at({neutron_count, proton_count, nucleon_count});
    } else if (type == LeptonInjector::Particle::ParticleType::Nucleon) {
        neutron_count = 0;
        proton_count = 0;
        nucleon_count = 1;
        is_atom = false;
        molar_mass = (atomic_masses.at({1, 0, 1}) + atomic_masses.at({0, 1, 1})) / 2.0;
    } else if (type == LeptonInjector::Particle::ParticleType::EMinus) {
        neutron_count = 0;
        proton_count = 0;
        nucleon_count = 0;
        is_atom = false;
        molar_mass = electron_molar_mass;
    } else {
        GetNucleonContent(static_cast<int>(type), neutron_count, proton_count, nucleon_count);
        is_atom = true;
        molar_mass = atomic_masses.at({neutron_count, proton_count, nucleon_count});
    }
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

    double neutrons_per_gram = 0;
    double nucleons_per_gram = 0;
    double protons_per_gram = 0;

    std::vector<MaterialComponent> material_components;

    for(auto const & mass_frac_it : component_mass_fractions) {
        int component_id = mass_frac_it.first;
        double component_mass_fraction = mass_frac_it.second;
        Component component(static_cast<LeptonInjector::Particle::ParticleType>(component_id));
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = component_mass_fraction;
        material_component.particle_density_over_total_mass_density = LeptonInjector::Constants::avogadro * component_mass_fraction / material_component.component.molar_mass;
        material_components.push_back(material_component);
        neutrons_per_gram += material_component.particle_density_over_total_mass_density * material_component.component.neutron_count;
        nucleons_per_gram += material_component.particle_density_over_total_mass_density * material_component.component.nucleon_count;
        protons_per_gram += material_component.particle_density_over_total_mass_density * material_component.component.proton_count;
    }

    if(neutrons_per_gram > 0) {
        Component component(LeptonInjector::Particle::ParticleType::Neutron);
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = neutrons_per_gram * component.molar_mass;
        material_component.particle_density_over_total_mass_density = neutrons_per_gram;
        material_components.push_back(material_component);
    }

    if(nucleons_per_gram > 0) {
        Component component(LeptonInjector::Particle::ParticleType::Nucleon);
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = nucleons_per_gram * component.molar_mass;
        material_component.particle_density_over_total_mass_density = nucleons_per_gram;
        material_components.push_back(material_component);
    }

    if(protons_per_gram > 0) {
        Component component(LeptonInjector::Particle::ParticleType::PPlus);
        MaterialComponent material_component;
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = protons_per_gram * component.molar_mass;
        material_component.particle_density_over_total_mass_density = protons_per_gram;
        material_components.push_back(material_component);
        component = Component(LeptonInjector::Particle::ParticleType::EMinus);
        material_component.component = component;
        material_component.mass_density_over_total_mass_density = protons_per_gram * component.molar_mass;
        material_component.particle_density_over_total_mass_density = protons_per_gram;
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
    bool fexists(const char *filename)
    {
        std::ifstream ifile(filename);
        return (bool)ifile;
    }
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

std::vector<LeptonInjector::Particle::ParticleType> MaterialModel::GetMaterialTargets(int material_id) const {
    std::vector<LeptonInjector::Particle::ParticleType> targets;
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
    return material_names_.size() > id;
}

double MaterialModel::GetTargetMassFraction(int material_id, LeptonInjector::Particle::ParticleType particle_type) const {
    std::pair<int, LeptonInjector::Particle::ParticleType> key(material_id, particle_type);
    if(material_components_by_id_.find(key) != material_components_by_id_.end())
        return material_components_by_id_.at(key).mass_density_over_total_mass_density;
    else
        return 0.0;
}

double MaterialModel::GetTargetParticleFraction(int material_id, LeptonInjector::Particle::ParticleType particle_type) const {
    std::pair<int, LeptonInjector::Particle::ParticleType> key(material_id, particle_type);
    if(material_components_by_id_.find(key) != material_components_by_id_.end())
        return material_components_by_id_.at(key).particle_density_over_total_mass_density;
    else
        return 0.0;
}

std::vector<double> MaterialModel::GetTargetMassFraction(int material_id, std::vector<LeptonInjector::Particle::ParticleType> const & particle_types) const {
    std::vector<double> fractions;
    fractions.reserve(particle_types.size());
    for(auto const & particle_type : particle_types) {
        std::pair<int, LeptonInjector::Particle::ParticleType> key(material_id, particle_type);
        if(material_components_by_id_.find(key) != material_components_by_id_.end())
            fractions.push_back(material_components_by_id_.at(key).mass_density_over_total_mass_density);
        else
            fractions.push_back(0.0);
    }
    return fractions;
}

std::vector<double> MaterialModel::GetTargetParticleFraction(int material_id, std::vector<LeptonInjector::Particle::ParticleType> const & particle_types) const {
    std::vector<double> fractions;
    fractions.reserve(particle_types.size());
    for(auto const & particle_type : particle_types) {
        std::pair<int, LeptonInjector::Particle::ParticleType> key(material_id, particle_type);
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

std::vector<double> MaterialModel::GetTargetRadiationFraction(int material_id, std::vector<LeptonInjector::Particle::ParticleType> const & particle_types) const {
    double X0inv = 0;
    std::vector<double> fractions;
    fractions.reserve(particle_types.size());
    for(auto const & particle_type : particle_types) {
        std::pair<int, LeptonInjector::Particle::ParticleType> key(material_id, particle_type);
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

void MaterialModel::GetNucleonContent(int code, int & neutron_count, int & proton_count, int & nucleon_count) {
    int prefix = 0;
    int suffix = 0;

    char buf[CHAR_BUF_SIZE];
    sprintf(buf, "%d", code);
    int nread = sscanf(buf, "%3d%3d%3d%1d", &prefix, &proton_count, &nucleon_count, &suffix);
    if (nread != 4) {
        throw std::runtime_error("Failed to convert nuclear pdg to A and Z "
                "prefix "+std::to_string(prefix)+", A "+std::to_string(nucleon_count)+", Z "+std::to_string(proton_count)+", suffix "+std::to_string(suffix));
    }
    neutron_count = nucleon_count - proton_count;
}

std::vector<LeptonInjector::Particle::ParticleType> MaterialModel::GetMaterialConstituents(int material_id) const {
    std::vector<LeptonInjector::Particle::ParticleType> particles;;
    particles.reserve(material_components_[material_id].size());
    for (auto const & component : material_components_[material_id]) {
        particles.push_back(component.component.type);
    }
    return particles;
}

double MaterialModel::GetMolarMass(LeptonInjector::Particle::ParticleType particle) {
    return Component(particle).molar_mass;
}

int MaterialModel::GetNucleonCount(LeptonInjector::Particle::ParticleType particle) {
    return Component(particle).nucleon_count;
}

int MaterialModel::GetNeutronCount(LeptonInjector::Particle::ParticleType particle) {
    return Component(particle).neutron_count;
}

int MaterialModel::GetProtonCount(LeptonInjector::Particle::ParticleType particle) {
    return Component(particle).proton_count;
}

