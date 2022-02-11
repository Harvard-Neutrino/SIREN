#include <map>
#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "earthmodel-service/MaterialModel.h"

using namespace earthmodel;

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

void MaterialModel::AddMaterial(std::string const & name, std::map<int, double> matratios) {
    double pne_ratio = ComputePNERatio(matratios);
    if(material_ids_.find(name) == material_ids_.end()) {
        int id = material_names_.size();

        // Fill basic information
        material_ids_.insert({name, id});
        material_names_.push_back(name);

        // Fill proton/neutron/nucleon number maps
        int np, nn;
        std::map<int, int> num_protons_, num_neutrons_, num_nucleons_;
        for (auto & k: matratios) {
            int pdgcode = k.first;
            GetNucleonContent(pdgcode, np, nn);
            num_protons_[pdgcode] = np;
            num_neutrons_[pdgcode] = nn;
            num_nucleons_[pdgcode] = np+nn;
        }
        material_num_protons_.insert({id, num_protons_});
        material_num_neutrons_.insert({id, num_neutrons_});
        material_num_nucleons_.insert({id, num_nucleons_});

        // Fill mass fraction, molar mass, and atomic fraction maps
        material_mass_frac_.insert({id, matratios});
        std::map<int,double> molar_masses = GetMolarMasses(num_protons_);
        material_molar_mass_.insert({id, molar_masses});
        double nfrac_denom = 0;
        for (auto& k : matratios) {nfrac_denom += matratios[k.first]/molar_masses[k.first];}
        std::map<int,double> atom_fracs;
        for (auto& k : matratios) {
            atom_fracs[k.first] = (matratios[k.first]/molar_masses[k.first])/nfrac_denom;
        }
        material_atom_frac_.insert({id,atom_fracs});


        // Fill particle type map
        std::set<LeptonInjector::Particle::ParticleType> ptypes;
        for (auto& k : matratios) {ptypes.insert(static_cast<LeptonInjector::Particle::ParticleType>(k.first));}
        material_constituents_.insert({id,ptypes});

        // Fill proton:electron ratio map
        pne_ratios_.insert({id, pne_ratio});
        
        // Compute radiation length 
        material_rad_length_.insert({id, ComputeRadLength(id)});
    }
    else {
        int id = material_ids_[name];
        material_mass_frac_[id] = matratios;
        pne_ratios_[id] = pne_ratio;
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
        throw("Received empty matratio filename!");

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
        throw("Cannot open matratio file!");
    }

    // check earthmodel file
    std::ifstream in(fname.c_str(), std::ifstream::in);

    if (in.fail())
        throw("Failed to open " + fname + ". Set correct path.");

    // read the file
    std::string buf;
    std::string medtype;
    int matpdg, nmats;
    double weight;
    int nread = 0;

    while(not (in.eof() or in.fail())) {

        std::getline(in, buf);
        nread = buf.size();

        if(nread == 0 || buf[0] == ' ' || buf[0] == '\t' || buf[0] == '#') {
            // new line, start from white space, or comment line.
            continue;
        } else {
            // material density data
            std::stringstream ss(buf);

            ss >> medtype >> nmats;

            std::map<int, double> matratio;
            for(int i=0; i<nmats; ++i) {
                std::getline(in, buf);
                nread = buf.size();
                if(nread == 0 || buf[0] == ' ' || buf[0] == '\t' || buf[0] == '#') {
                    // new line, start from white space, or comment line.
                    --i;
                    continue;
                } else {
                    std::stringstream ss2(buf);
                    ss2 >> matpdg >> weight;
                    matratio[matpdg] = weight;
                }
            }
            AddMaterial(medtype, matratio);
        }

    } // end of the while loop
    model_files_.push_back(matratio);

    in.close();
}

double MaterialModel::ComputePNERatio(std::map<int, double> const & mats) const {
    // calculate P, N, E ratio
    double tot_np = 0;
    double tot_nn = 0;
    int np, nn;
    for(auto const & it : mats) {
        int pdg = it.first;
        GetNucleonContent(pdg, np, nn);
        tot_np += np*it.second;
        tot_nn += nn*it.second;
    }

    double tot_A = tot_np + tot_nn;
    if(tot_A==0)
        tot_A=1; //avoid division by zero

    double nw_proton = tot_np / tot_A;
    double nw_electron = nw_proton;
    //double nw_neutron = tot_nn / tot_z;

    return nw_electron;
}

double MaterialModel::ComputeRadLength(int id) {
    // This function calculates the radiation length of a given material in g/cm^2
    // Takes screening effects into account
    // Averages over constituent materials in a composite
    // See Page 21 of Particle Detectors by Grupen and Shwartz
    
    double X0inv = 0;
    int i;
    double X0i, Z, A, f;
    std::map<int, double>::iterator it;
    for (it = material_molar_mass_[id].begin(); it != material_molar_mass_[id].end(); it++) {
        i = it->first;
        A = it->second;
        Z = (material_num_protons_[id])[i];  
        f = (material_mass_frac_[id])[i];  
        X0i = 716.4 * A / ( Z*(Z + 1) * std::log(287./std::sqrt(Z))); // g/cm^2, Grupen eq 1.59
        X0inv += f/X0i;
    }
    return 1/X0inv;
}

std::map<int, double> MaterialModel::GetMolarMasses(std::map<int, int> const & pnums) const {
    std::map<int, double> molar_masses;
    for(auto const & it : pnums) {
        molar_masses[it.first] = molar_mass_table.at(it.second);
    }
    return molar_masses;
}

double MaterialModel::GetPNERatio(int id) const {
    return pne_ratios_.at(id);
}

std::string MaterialModel::GetMaterialName(int id) const {
    return material_names_.at(id);
}

int MaterialModel::GetMaterialId(std::string const & name) const {
    return material_ids_.at(name);
}

bool MaterialModel::HasMaterial(std::string const & name) const {
    return material_ids_.count(name) > 0;
}

bool MaterialModel::HasMaterial(int id) const {
    return material_names_.size() > id;
}

std::set<LeptonInjector::Particle::ParticleType> MaterialModel::GetMaterialConstituents(int id) const {
    return material_constituents_.at(id);
}

std::map<int, double> MaterialModel::GetMaterialMassFracs(int id) const {
    return material_mass_frac_.at(id);
}

std::map<int, double> MaterialModel::GetMaterialAtomFracs(int id) const {
    return material_atom_frac_.at(id);
}

std::map<int, int> MaterialModel::GetMaterialNumNucleons(int id) const {
    return material_num_nucleons_.at(id);
}

std::map<int, int> MaterialModel::GetMaterialNumProtons(int id) const {
    return material_num_protons_.at(id);
}

std::map<int, int> MaterialModel::GetMaterialNumNeutrons(int id) const {
    return material_num_neutrons_.at(id);
}

double MaterialModel::GetMaterialRadLength(int id) const {
    return material_rad_length_.at(id);
}

double MaterialModel::GetTargetListMassFrac(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {

    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_mass_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second;
    }
    return sum;
}

double MaterialModel::GetTargetListAtomFrac(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second;
    }
    return sum;
}

double MaterialModel::GetTargetListNucleonFrac(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double num = 0, dem = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        dem += it.second*material_num_nucleons_.at(id).at(it.first);
        if(t.count(pdg) > 0)
            num += it.second;
    }
    return num/dem;
}

double MaterialModel::GetTargetListProtonFrac(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double num = 0, dem = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        dem += it.second*material_num_protons_.at(id).at(it.first);
        if(t.count(pdg) > 0)
            num += it.second;
    }
    return num/dem;
}

double MaterialModel::GetTargetListNeutronFrac(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double num = 0, dem = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        dem += it.second*material_num_neutrons_.at(id).at(it.first);
        if(t.count(pdg) > 0)
            num += it.second;
    }
    return num/dem;
}

double MaterialModel::GetTargetListAtomsToMass(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_mass_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second/material_molar_mass_.at(id).at(it.first);
    }
    return NA*sum;
}

double MaterialModel::GetTargetListNucleonsToMass(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_mass_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second/material_molar_mass_.at(id).at(it.first) * material_num_nucleons_.at(id).at(it.first);
    }
    return NA*sum;
}

double MaterialModel::GetTargetListProtonsToMass(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_mass_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second/material_molar_mass_.at(id).at(it.first) * material_num_protons_.at(id).at(it.first);
    }
    return NA*sum;
}

double MaterialModel::GetTargetListNeutronsToMass(int id, std::set<LeptonInjector::Particle::ParticleType> const & targets) const {
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_mass_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second/material_molar_mass_.at(id).at(it.first) * material_num_neutrons_.at(id).at(it.first);
    }
    return NA*sum;
}

void MaterialModel::GetNucleonContent(int code, int & np, int & nn) {
    np = 0;
    int z = 0;
    int prefix = 0;
    int suffix = 0;

    char buf[CHAR_BUF_SIZE];
    sprintf(buf, "%d", code);
    int nread = sscanf(buf, "%3d%3d%3d%1d", &prefix, &np, &z, &suffix);
    if (nread != 4) {
        throw std::runtime_error("Failed to convert nuclear pdg to A and Z "
                "prefix "+std::to_string(prefix)+", A "+std::to_string(np)+", Z "+std::to_string(z)+", suffix "+std::to_string(suffix));
    }
    nn = z - np;
}

const std::map<int, double> MaterialModel::molar_mass_table {

                    {1, 1.007},
                    {2, 4.002},
                    {3, 6.941},
                    {4, 9.012},
                    {5, 10.811},
                    {6, 12.011},
                    {7, 14.007},
                    {8, 15.999},
                    {9, 18.998},
                    {10, 20.18},
                    {11, 22.99},
                    {12, 24.305},
                    {13, 26.982},
                    {14, 28.086},
                    {15, 30.974},
                    {16, 32.065},
                    {17, 35.453},
                    {18, 39.948},
                    {19, 39.098},
                    {20, 40.078},
                    {21, 44.956},
                    {22, 47.867},
                    {23, 50.942},
                    {24, 51.996},
                    {25, 54.938},
                    {26, 55.845},
                    {27, 58.933},
                    {28, 58.693},
                    {29, 63.546},
                    {30, 65.38},
                    {31, 69.723},
                    {32, 72.64},
                    {33, 74.922},
                    {34, 78.96},
                    {35, 79.904},
                    {36, 83.798},
                    {37, 85.468},
                    {38, 87.62},
                    {39, 88.906},
                    {40, 91.224},
                    {41, 92.906},
                    {42, 95.96},
                    {43, 98.0},
                    {44, 101.07},
                    {45, 102.906},
                    {46, 106.42},
                    {47, 107.868},
                    {48, 112.411},
                    {49, 114.818},
                    {50, 118.71},
                    {51, 121.76},
                    {52, 127.6},
                    {53, 126.904},
                    {54, 131.293},
                    {55, 132.905},
                    {56, 137.327},
                    {57, 138.905},
                    {58, 140.116},
                    {59, 140.908},
                    {60, 144.242},
                    {61, 145.0},
                    {62, 150.36},
                    {63, 151.964},
                    {64, 157.25},
                    {65, 158.925},
                    {66, 162.5},
                    {67, 164.93},
                    {68, 167.259},
                    {69, 168.934},
                    {70, 173.054},
                    {71, 174.967},
                    {72, 178.49},
                    {73, 180.948},
                    {74, 183.84},
                    {75, 186.207},
                    {76, 190.23},
                    {77, 192.217},
                    {78, 195.084},
                    {79, 196.967},
                    {80, 200.59},
                    {81, 204.383},
                    {82, 207.2},
                    {83, 208.98},
                    {84, 210.0},
                    {85, 210.0},
                    {86, 222.0},
                    {87, 223.0},
                    {88, 226.0},
                    {89, 227.0},
                    {90, 232.038},
                    {91, 231.036},
                    {92, 238.029},
                    {93, 237.0},
                    {94, 244.0},
                    {95, 243.0},
                    {96, 247.0},
                    {97, 247.0},
                    {98, 251.0},
                    {99, 252.0},
                    {100, 257.0},
                    {101, 258.0},
                    {102, 259.0},
                    {103, 262.0},
                    {104, 261.0},
                    {105, 262.0},
                    {106, 266.0},
                    {107, 264.0},
                    {108, 267.0},
                    {109, 268.0},
                    {110, 271.0},
                    {111, 272.0},
                    {112, 285.0},
                    {113, 284.0},
                    {114, 289.0},
                    {115, 288.0},
                    {116, 292.0},
                    {117, 295.0},
                    {118, 294.0} };
