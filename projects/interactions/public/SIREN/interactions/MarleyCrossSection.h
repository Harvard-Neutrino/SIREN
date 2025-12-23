#pragma once
#ifndef SIREN_MarleyCrossSection_H
#define SIREN_MarleyCrossSection_H

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/utilities/Random.h" // for SIREN_random
#include "SIREN/interactions/Interaction.h" // for Interaction
#include "SIREN/interactions/CrossSection.h"

#include "marley/Generator.hh"
#include "marley/FileManager.hh"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }
namespace siren { namespace interactions { class MarleyCrossSection; } }

namespace siren {
namespace interactions {

namespace marley_ {

class FileManager_ : public ::marley::FileManager {
friend class siren::interactions::MarleyCrossSection;
    static void set_search_path(std::string const & search_path) {
        ::marley::FileManager::default_search_path_ = search_path;
    }
    static void set_marley_path(std::string const & marley_path) {
        std::string default_search_path;
        default_search_path = marley_path + "/data";
        default_search_path += ':' + marley_path + "/data/react";
        default_search_path += ':' + marley_path + "/data/structure";
        set_search_path(default_search_path);
    }
};

} // namespace marley

class MarleyCrossSection : public CrossSection {
//friend cereal::access;
private:
    std::vector<char> marley_react_data_;
    std::vector<char> marley_nuclide_index_data_;
    std::vector<std::vector<char>> marley_nuclide_data_;
    std::vector<char> marley_masses_data_;
    std::vector<char> marley_gs_parity_data_;
    std::string marley_react_fname_;
    std::string marley_nuclide_index_fname_;
    std::vector<std::string> marley_nuclide_fnames_;
    std::string marley_masses_fname_;
    std::string marley_gs_parity_fname_;
    std::vector<std::unique_ptr<marley::Reaction>> reactions_;
    std::unique_ptr<marley::StructureDatabase> structure_database_;
    void InitializeMarley(const std::string& marley_config);
    bool has_nu_cc;
    bool has_nubar_cc;
    bool has_nc;
    bool has_elastic;

public:
    MarleyCrossSection(std::string marley_react_file, std::string marley_nuclide_index_file, std::vector<std::string> marley_nuclide_files, std::string marley_masses_file, std::string marley_gs_parity_file);
    MarleyCrossSection(std::array<std::vector<char>, 4> const & data, std::vector<std::vector<char>> const & nuclide_data, std::array<std::string, 4> const & fnames, std::vector<std::string> const & nuclide_fnames);
    virtual ~MarleyCrossSection() {};
    virtual bool equal(CrossSection const & other) const override;
    virtual double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    virtual double TotalCrossSectionAllFinalStates(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    virtual double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;

    virtual std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    virtual std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    virtual std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("MarleyReactData", marley_react_data_));
            archive(::cereal::make_nvp("MarleyNuclideIndexData", marley_nuclide_index_data_));
            archive(::cereal::make_nvp("MarleyNuclideData", marley_nuclide_data_));
            archive(::cereal::make_nvp("MarleyMassesData", marley_masses_data_));
            archive(::cereal::make_nvp("MarleyGSParityData", marley_gs_parity_data_));
            archive(::cereal::make_nvp("MarleyReactFname", marley_react_fname_));
            archive(::cereal::make_nvp("MarleyNuclideIndexFname", marley_nuclide_index_fname_));
            archive(::cereal::make_nvp("MarleyNuclideFnames", marley_nuclide_fnames_));
            archive(::cereal::make_nvp("MarleyMassesFname", marley_masses_fname_));
            archive(::cereal::make_nvp("MarleyGSParityFname", marley_gs_parity_fname_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("MarleyCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<MarleyCrossSection> & construct, std::uint32_t const version) {
        if(version == 0) {
            std::array<std::vector<char>, 4> data;
            std::array<std::string, 4> fnames;
            std::vector<std::string> nuclide_fnames;
            std::vector<std::vector<char>> nuclide_data;
            archive(::cereal::make_nvp("MarleyReactData", data[0]));
            archive(::cereal::make_nvp("MarleyNuclideIndexData", data[1]));
            archive(::cereal::make_nvp("MarleyNuclideData", nuclide_data));
            archive(::cereal::make_nvp("MarleyMassesData", data[2]));
            archive(::cereal::make_nvp("MarleyGSParityData", data[3]));
            archive(::cereal::make_nvp("MarleyReactFname", fnames[0]));
            archive(::cereal::make_nvp("MarleyNuclideIndexFname", fnames[1]));
            archive(::cereal::make_nvp("MarleyNuclideFnames", nuclide_fnames));
            archive(::cereal::make_nvp("MarleyMassesFname", fnames[2]));
            archive(::cereal::make_nvp("MarleyGSParityFname", fnames[3]));
            construct(data, nuclide_data, fnames, nuclide_fnames);
            archive(cereal::virtual_base_class<CrossSection>(construct.ptr()));
        } else {
            throw std::runtime_error("MarleyCrossSection only supports version <= 0!");
        }
    }
}; // class MarleyCrossSection

} // namespace interactions
} // namespace siren


CEREAL_CLASS_VERSION(siren::interactions::MarleyCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::MarleyCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::MarleyCrossSection);

#endif // SIREN_MarleyCrossSection_H
