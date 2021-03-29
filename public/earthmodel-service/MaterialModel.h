#ifndef LI_MaterialModel_H
#define LI_MaterialModel_H

#include <map>
#include <string>
#include <vector>

namespace earthmodel {

class MaterialModel {
    static constexpr const int CHAR_BUF_SIZE = 8196;
private:
    std::string path_;
    std::vector<std::string> model_files_;

    std::vector<std::string> material_names_;
    std::map<std::string, int> material_ids_;
    std::map<int, std::map<int, double> > material_maps_;
    std::map<int, double> pne_ratios_;
public:
    MaterialModel();
    MaterialModel(std::string const & file);
    MaterialModel(std::string const & path, std::string const & file);
    MaterialModel(std::vector<std::string> const & files);
    MaterialModel(std::string const & path, std::vector<std::string> const & files);

    void SetPath(std::string const & path);
    void AddMaterial(std::string const & name, std::map<int, double> matratios);
    //void AddMaterial(std::string const & name, double pne_ratio);
    void AddModelFiles(std::vector<std::string> const & matratios);
    void AddModelFile(std::string matratio);

    double GetPNERatio(int id) const;
    std::string GetMaterialName(int id) const;
    int GetMaterialId(std::string const & name) const;
    bool HasMaterial(std::string const & name) const;
    bool HasMaterial(int) const;
    std::map<int, double> GetMaterialMap(int id) const;
private:
    double ComputePNERatio(std::map<int, double> const & matratios) const;
public:
    static void GetAZ(int code, int & np, int & nn);
};

} // namespace earthmodel

# endif // LI_MaterialModel_H

