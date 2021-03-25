#include <map>
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
        material_ids_.insert({name, id});
        material_names_.push_back(name);
        material_maps_.insert({id, matratios});
        pne_ratios_.insert({id, pne_ratio});
    }
    else {
        int id = material_ids_[name];
        material_maps_[id] = matratios;
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

double MaterialModel::ComputePNERatio(std::map<int, double> const & mats) {
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

double MaterialModel::GetPNERatio(int id) {
    return pne_ratios_[id];
}

std::string MaterialModel::GetMaterialName(int id) {
    return material_names_[id];
}

int MaterialModel::GetMaterialId(std::string const & name) {
    return material_ids_[name];
}

bool MaterialModel::HasMaterial(std::string const & name) {
    return material_ids_.count(name) > 0;
}

bool MaterialModel::HasMaterial(int id) {
    return material_names_.size() > id;
}

std::map<int, double> MaterialModel::GetMaterialMap(int id) {
    return material_maps_[id];
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

