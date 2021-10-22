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
        int np,nn;
        std::map<int, int> num_protons_,num_neutrons_,num_nucleons_;
        for (auto& k: matratios)
        {
					GetAZ(k.first,np,nn);
					num_protons_[k.first] = np;
					num_neutrons_[k.first] = nn;
					num_nucleons_[k.first] = np+nn;
				}
				material_num_protons_.insert({id,num_protons_});
				material_num_neutrons_.insert({id,num_neutrons_});
				material_num_nucleons_.insert({id,num_nucleons_});
        
        // Fill mass fraction, molar mass, and atomic fraction maps
        material_mass_frac_.insert({id, matratios});
        std::map<int,double> molar_masses = GetMolarMasses(num_protons_);
        material_molar_mass_.insert({id, molar_masses});
        double nfrac_denom = 0;
        for (auto& k : matratios) {nfrac_denom += matratios[k.first]/molar_masses[k.first];}
        std::map<int,double> atom_fracs;
        for (auto& k : matratios) {atom_fracs[k.first] = (matratios[k.first]/molar_masses[k.first])/nfrac_denom;}
        material_atom_frac_.insert({id,atom_fracs});


        // Fill particle type map
        std::vector<LeptonInjector::Particle::ParticleType> ptypes;
        for (auto& k : matratios) {ptypes.push_back(static_cast<LeptonInjector::Particle::ParticleType>(k.first));}
        material_constituents_.insert({id,ptypes});
        
        // Fill proton:electron ratio map
        pne_ratios_.insert({id, pne_ratio});
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
        GetAZ(pdg, np, nn);
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

std::map<int, double> MaterialModel::GetMolarMasses(std::map<int, int> const & pnums) const {
		
		std::ifstream ifs("AtomicData.csv");
		std::string line;
		int nproton; double molmass;
		std::map<int, double> molar_masses;

		while(std::getline(ifs, line))
		{
				std::stringstream linestream(line);
				linestream >> nproton >> molmass;
				auto result = std::find_if(pnums.begin(),pnums.end(),[nproton](const auto& mo) {return mo.second==nproton; });
				if(result != pnums.end()) molar_masses[result->first] = nproton;
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

std::vector<LeptonInjector::Particle::ParticleType> MaterialModel::GetMaterialConstituents(int id) const {
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
    return material_num_protons.at(id);
}

std::map<int, int> MaterialModel::GetMaterialNumNeutrons(int id) const {
    return material_num_neutrons.at(id);
}

double MaterialModel::GetTargetListMassFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const {
    
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_mass_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second;
    }
    return sum;
}

double MaterialModel::GetTargetListAtomFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const {
    double sum = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        if(t.count(pdg) > 0)
            sum += it.second;
    }
    return sum;
}

double MaterialModel::GetTargetListNucleonFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const {
    double num = 0, dem = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        dem += it.second*material_num_nucleons_.at(id)[it.first];
        if(t.count(pdg) > 0)
            num += it.second;
    }
    return num/dem;
}

double MaterialModel::GetTargetListProtonFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const {
    double num = 0, dem = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        dem += it.second*material_num_protons_.at(id)[it.first];
        if(t.count(pdg) > 0)
            num += it.second;
    }
    return num/dem;
}

double MaterialModel::GetTargetListNeutronFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const {
    double num = 0, dem = 0;
    std::set<LeptonInjector::Particle::ParticleType> t(targets.begin(), targets.end());
    for(auto const & it : material_atom_frac_.at(id)) {
        LeptonInjector::Particle::ParticleType pdg = (LeptonInjector::Particle::ParticleType)it.first;
        dem += it.second*material_num_neutrons_.at(id)[it.first];
        if(t.count(pdg) > 0)
            num += it.second;
    }
    return num/dem;
}

void MaterialModel::GetAZ(int code, int & np, int & nn) {
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

