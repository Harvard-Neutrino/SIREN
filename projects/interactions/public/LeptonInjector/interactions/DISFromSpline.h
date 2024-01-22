#pragma once
#ifndef LI_DISFromSpline_H
#define LI_DISFromSpline_H

#include <set>                                                // for set
#include <map>                                                // for map
#include <memory>
#include <vector>                                             // for vector
#include <cstdint>                                            // for uint32_t
#include <utility>                                            // for pair
#include <algorithm>
#include <stdexcept>                                          // for runtime...

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include "LeptonInjector/interactions/CrossSection.h"        // for CrossSe...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle

namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace interactions {

class DISFromSpline : public CrossSection {
friend cereal::access;
private:
    photospline::splinetable<> differential_cross_section_;
    photospline::splinetable<> total_cross_section_;

    std::vector<dataclasses::InteractionSignature> signatures_;
    std::set<LI::dataclasses::Particle::ParticleType> primary_types_;
    std::set<LI::dataclasses::Particle::ParticleType> target_types_;
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<LI::dataclasses::Particle::ParticleType>> targets_by_primary_types_;
    std::map<std::pair<LI::dataclasses::Particle::ParticleType, LI::dataclasses::Particle::ParticleType>, std::vector<dataclasses::InteractionSignature>> signatures_by_parent_types_;

    int interaction_type_;
    double target_mass_;
    double minimum_Q2_;
    
    double unit;

public:
    DISFromSpline();
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::set<LI::dataclasses::Particle::ParticleType> primary_types, std::set<LI::dataclasses::Particle::ParticleType> target_types, std::string units = "cm");
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::vector<LI::dataclasses::Particle::ParticleType> primary_types, std::vector<LI::dataclasses::Particle::ParticleType> target_types, std::string units = "cm");
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LI::dataclasses::Particle::ParticleType> primary_types, std::set<LI::dataclasses::Particle::ParticleType> target_types, std::string units = "cm");
    DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LI::dataclasses::Particle::ParticleType> primary_types, std::set<LI::dataclasses::Particle::ParticleType> target_types, std::string units = "cm");
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LI::dataclasses::Particle::ParticleType> primary_types, std::vector<LI::dataclasses::Particle::ParticleType> target_types, std::string units = "cm");
    DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LI::dataclasses::Particle::ParticleType> primary_types, std::vector<LI::dataclasses::Particle::ParticleType> target_types, std::string units = "cm");
    
    void SetUnits(std::string units);

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy) const;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass, double Q2=std::numeric_limits<double>::quiet_NaN()) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random> random) const override;

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const override;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const override;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    void LoadFromFile(std::string differential_filename, std::string total_filename);
    void LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data);

    double GetMinimumQ2() const {return minimum_Q2_;};
    double GetTargetMass() const {return target_mass_;};
    int GetInteractionType() const {return interaction_type_;};

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            splinetable_buffer buf;
            buf.size = 0;
            auto result_obj = differential_cross_section_.write_fits_mem();
            buf.data = result_obj.first;
            buf.size = result_obj.second;

            std::vector<char> diff_blob;
            diff_blob.resize(buf.size);
            std::copy((char*)buf.data, (char*)buf.data + buf.size, &diff_blob[0]);

            archive(::cereal::make_nvp("DifferentialCrossSectionSpline", diff_blob));

            buf.size = 0;
            result_obj = total_cross_section_.write_fits_mem();
            buf.data = result_obj.first;
            buf.size = result_obj.second;

            std::vector<char> total_blob;
            total_blob.resize(buf.size);
            std::copy((char*)buf.data, (char*)buf.data + buf.size, &total_blob[0]);

            archive(::cereal::make_nvp("TotalCrossSectionSpline", total_blob));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(::cereal::make_nvp("TargetTypes", target_types_));
            archive(::cereal::make_nvp("InteractionType", interaction_type_));
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DISFromSpline only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            std::vector<char> differential_data;
            std::vector<char> total_data;
            archive(::cereal::make_nvp("DifferentialCrossSectionSpline", differential_data));
            archive(::cereal::make_nvp("TotalCrossSectionSpline", total_data));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(::cereal::make_nvp("TargetTypes", target_types_));
            archive(::cereal::make_nvp("InteractionType", interaction_type_));
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(cereal::virtual_base_class<CrossSection>(this));
            LoadFromMemory(differential_data, total_data);
            InitializeSignatures();
        } else {
            throw std::runtime_error("DISFromSpline only supports version <= 0!");
        }
    }
private:
    void ReadParamsFromSplineTable();
    void InitializeSignatures();
};

} // namespace interactions
} // namespace LI

CEREAL_CLASS_VERSION(LI::interactions::DISFromSpline, 0);
CEREAL_REGISTER_TYPE(LI::interactions::DISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::interactions::CrossSection, LI::interactions::DISFromSpline);

#endif // LI_DISFromSpline_H
