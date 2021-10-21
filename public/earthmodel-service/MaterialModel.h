#ifndef LI_MaterialModel_H
#define LI_MaterialModel_H

#include <map>
#include <string>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include "serialization/array.h"

#include "LeptonInjector/Particle.h"

namespace earthmodel {

class MaterialModel {
    static constexpr const int CHAR_BUF_SIZE = 8196;
private:
    std::string path_;
    std::vector<std::string> model_files_;

    std::vector<std::string> material_names_;
    std::map<std::string, int> material_ids_;
    std::map<int, std::map<int, double> > material_maps_;
    std::map<int, std::vector<LeptonInjector::Particle::ParticleType> > material_constituents_;
    std::map<int, double> pne_ratios_;
public:
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Path", path_));
            archive(cereal::make_nvp("MaterialNames", material_names_));
            archive(cereal::make_nvp("MaterialIDs", material_ids_));
            archive(cereal::make_nvp("MaterialMaps", material_maps_));
            archive(cereal::make_nvp("MaterialConstituents", material_constituents_));
            archive(cereal::make_nvp("PNERatios", pne_ratios_));
        } else {
            throw std::runtime_error("MaterialModel only supports version <= 0!");
        }
    }
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
    double GetTargetComposition(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    std::string GetMaterialName(int id) const;
    int GetMaterialId(std::string const & name) const;
    bool HasMaterial(std::string const & name) const;
    bool HasMaterial(int) const;
    std::map<int, double> GetMaterialMap(int id) const;
    std::vector<LeptonInjector::Particle::ParticleType> GetMaterialConstituents(int id) const;
private:
    double ComputePNERatio(std::map<int, double> const & matratios) const;
public:
    static void GetAZ(int code, int & np, int & nn);
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::MaterialModel, 0);

# endif // LI_MaterialModel_H

