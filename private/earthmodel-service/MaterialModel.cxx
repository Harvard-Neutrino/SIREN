#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include <earthmodel-service/MaterialModel.h>

using namespace earthmodel;

MaterialModel::MaterialModel() {}

MaterialModel::MaterialModel(std::string const & path) : path_(path) {}

MaterialModel::MaterialModel(std::string const & path, std::string const & matratio) : path_(path) {
    AddModelFile(matratio);
}

MaterialModel::MaterialModel(std::string const & path, std::vector<std::string> const & matratios) : path_(path) {
    AddModelFiles(matratios);
}

void MaterialModel::SetPath(std::string const & path) {
    path_ = path;
}

void MaterialModel::AddMaterial(std::string const & name, int matpdg, std::map<int, double> matratios) {
    AddMaterial(name, ComputePNERatio(matpdg), matratios);
}

void MaterialModel::AddMaterial(std::string const & name, double pne_ratio, std::map<int, double> matratios) {
    int id = material_names_.size();
    material_ids_.insert({name, id});
    material_names_.push_back(name);
    material_maps_.insert({id, matratios});
    pne_ratios_.insert({id, pne_ratio});
}

void MaterialModel::AddModelFiles(std::vector<std::string> const & matratios) {
    for(auto matratio : matratios)
        AddModelFile(matratio);
}

void MaterialModel::AddModelFile(std::string matratio) {
    if (matratio.find(".dat") == std::string::npos)
        matratio += ".dat";

    // check earthmodel file
    std::string fname = (matratio.find('/')==std::string::npos ? path_ + "materials/" + matratio : matratio);
    std::ifstream in(fname.c_str(), std::ifstream::in);

    if (in.fail())
        throw("failed to open " + fname + ". Set correct path.");

    // read the file
    const int bufsize = CHAR_BUF_SIZE;
    char buf[bufsize];
    std::string medtype;
    int matpdg, nmats;
    double weight;
    int nread = 0;

    while(!in.eof()) {

        in.getline(buf, bufsize);
        nread = in.gcount();

        if(nread == -1) {
            throw("getline failed");

        } else if(nread == 1 || buf[0] == ' ' || buf[0] == '#') {
            // new line, start from white space, or comment line.
            continue;
        } else {
            // material density data
            std::stringstream ss(buf);

            ss >> medtype >> nmats;

            std::map<int, double> matratio;
            for(int i=0; i<nmats; ++i) {
                in.getline(buf, bufsize);
                nread = in.gcount();
                if(nread == -1) {
                    throw("getline failed");
                } else if(nread == 1 || buf[0] == ' ' || buf[0] == '#') {
                    // new line, start from white space, or comment line.
                    --i;
                    continue;
                } else {
                    std::stringstream ss2(buf);
                    ss2 >> matpdg >> weight;
                    matratio[matpdg] = weight;
                }
            }
            AddMaterial(medtype, matpdg, matratio);
        }

    } // end of the while loop
    model_files_.push_back(matratio);

    in.close();
}

double MaterialModel::ComputePNERatio(int id) {
    // calculate P, N, E ratio
    std::map<int, double> &mats = material_maps_[id];
    double tot_np = 0;
    double tot_nn = 0;
    int np, nn;
    for(auto const & it : mats) {
        int pdg = it.first;
        GetAZ(pdg, np, nn);
        tot_np += np*it.second;
        tot_nn += nn*it.second;
    }

    int tot_z = tot_np + tot_nn;
    if(tot_z==0)
        tot_z=1; //avoid division by zero

    double nw_proton = tot_np / tot_z;
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

