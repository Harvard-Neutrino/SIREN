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
#include "LeptonInjector/Constants.h"

namespace earthmodel {

class MaterialModel {
    static constexpr const int CHAR_BUF_SIZE = 8196;
    static constexpr const double NA = 6.02e23;
    static const std::map<int, double> molar_mass_table;
private:
    std::string path_;
    std::vector<std::string> model_files_;

    std::vector<std::string> material_names_;
    std::map<std::string, int> material_ids_;
    std::map<int, std::map<int, double> > material_mass_frac_;
    std::map<int, std::map<int, double> > material_atom_frac_;
    std::map<int, std::map<int, double> > material_molar_mass_;
    std::map<int, std::map<int, int> > material_num_protons_;
    std::map<int, std::map<int, int> > material_num_neutrons_;
    std::map<int, std::map<int, int> > material_num_nucleons_;
    std::map<int, double > material_rad_length_;
    std::map<int, std::vector<LeptonInjector::Particle::ParticleType> > material_constituents_;
    std::map<int, double> pne_ratios_;


public:
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Path", path_));
            archive(cereal::make_nvp("MaterialNames", material_names_));
            archive(cereal::make_nvp("MaterialIDs", material_ids_));
            archive(cereal::make_nvp("MaterialMaps", material_mass_frac_));
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
    std::string GetMaterialName(int id) const;
    int GetMaterialId(std::string const & name) const;
    bool HasMaterial(std::string const & name) const;
    bool HasMaterial(int) const;
    std::vector<LeptonInjector::Particle::ParticleType> GetMaterialConstituents(int id) const;
    std::map<int, double> GetMaterialMassFracs(int id) const;
    std::map<int, double> GetMaterialAtomFracs(int id) const;
    std::map<int, int> GetMaterialNumNucleons(int id) const;
    std::map<int, int> GetMaterialNumProtons(int id) const;
    std::map<int, int> GetMaterialNumNeutrons(int id) const;
    double GetMaterialRadLength(int id) const;
    
    double GetTargetListMassFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListAtomFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListNucleonFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListProtonFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListNeutronFrac(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    
    double GetTargetListAtomsToMass(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListNucleonsToMass(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListProtonsToMass(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
    double GetTargetListNeutronsToMass(int id, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;
private:
    double ComputePNERatio(std::map<int, double> const & matratios) const;
    double ComputeRadLength(int id);
    std::map<int, double> GetMolarMasses(std::map<int, int> const & pnums) const;
public:
    static void GetNucleonContent(int code, int & np, int & nn);
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::MaterialModel, 0);

# endif // LI_MaterialModel_H

