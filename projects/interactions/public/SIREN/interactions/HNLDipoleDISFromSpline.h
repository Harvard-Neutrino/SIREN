#pragma once
#ifndef SIREN_HNLDipoleDISFromSpline_H
#define SIREN_HNLDipoleDISFromSpline_H

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

#include "SIREN/interactions/CrossSection.h"        // for CrossSe...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle

namespace siren { namespace dataclasses { struct InteractionRecord; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class HNLDipoleDISFromSpline : public CrossSection {
friend cereal::access;
private:
    photospline::splinetable<> differential_cross_section_;
    photospline::splinetable<> total_cross_section_;

    double hnl_mass_;
    std::vector<double> dipole_coupling_;  // d_e, d_mu, d_tau
    double target_mass_;
    double minimum_Q2_;

    std::vector<dataclasses::InteractionSignature> signatures_;
    std::set<siren::dataclasses::Particle::ParticleType> primary_types_;
    std::set<siren::dataclasses::Particle::ParticleType> target_types_;
    std::map<siren::dataclasses::Particle::ParticleType, std::vector<siren::dataclasses::Particle::ParticleType>> targets_by_primary_types_;
    std::map<std::pair<siren::dataclasses::Particle::ParticleType, siren::dataclasses::Particle::ParticleType>, std::vector<dataclasses::InteractionSignature>> signatures_by_parent_types_;


    double unit;

public:
    HNLDipoleDISFromSpline();
    HNLDipoleDISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, double hnl_mass, std::vector<double> dipole_coupling, double target_mass, double minumum_Q2, std::set<siren::dataclasses::ParticleType> primary_types, std::set<siren::dataclasses::ParticleType> target_types, std::string units = "cm");
    HNLDipoleDISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, double hnl_mass, std::vector<double> dipole_coupling, double target_mass, double minumum_Q2, std::vector<siren::dataclasses::ParticleType> primary_types, std::vector<siren::dataclasses::ParticleType> target_types, std::string units = "cm");
    HNLDipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> dipole_coupling, double target_mass, double minumum_Q2, std::set<siren::dataclasses::ParticleType> primary_types, std::set<siren::dataclasses::ParticleType> target_types, std::string units = "cm");
    HNLDipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> dipole_coupling, std::set<siren::dataclasses::ParticleType> primary_types, std::set<siren::dataclasses::ParticleType> target_types, std::string units = "cm");
    HNLDipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> dipole_coupling, double target_mass, double minumum_Q2, std::vector<siren::dataclasses::ParticleType> primary_types, std::vector<siren::dataclasses::ParticleType> target_types, std::string units = "cm");
    HNLDipoleDISFromSpline(std::string differential_filename, std::string total_filename, double hnl_mass, std::vector<double> dipole_coupling, std::vector<siren::dataclasses::ParticleType> primary_types, std::vector<siren::dataclasses::ParticleType> target_types, std::string units = "cm");
    
    void SetUnits(std::string units);

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(siren::dataclasses::Particle::ParticleType primary_type, double energy) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(siren::dataclasses::Particle::ParticleType primary_type, double energy, double x, double y, double Q2=std::numeric_limits<double>::quiet_NaN()) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random> random) const override;

    std::vector<siren::dataclasses::Particle::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::Particle::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::Particle::ParticleType primary_type, siren::dataclasses::Particle::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    void LoadFromFile(std::string differential_filename, std::string total_filename);
    void LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data);

    double GetMinimumQ2() const {return minimum_Q2_;};
    double GetTargetMass() const {return target_mass_;};
    double GetHNLMass() const {return hnl_mass_;};

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
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("HNLMass", hnl_mass_));
            archive(::cereal::make_nvp("DipoleCoupling", dipole_coupling_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("HNLDipoleDISFromSpline only supports version <= 0!");
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
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("HNLMass", hnl_mass_));
            archive(::cereal::make_nvp("DipoleCoupling", dipole_coupling_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(cereal::virtual_base_class<CrossSection>(this));
            LoadFromMemory(differential_data, total_data);
            InitializeSignatures();
        } else {
            throw std::runtime_error("HNLDipoleDISFromSpline only supports version <= 0!");
        }
    }
private:
    void ReadParamsFromSplineTable();
    void InitializeSignatures();
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::HNLDipoleDISFromSpline, 0);
CEREAL_REGISTER_TYPE(siren::interactions::HNLDipoleDISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::HNLDipoleDISFromSpline);

#endif // SIREN_HNLDipoleDISFromSpline_H