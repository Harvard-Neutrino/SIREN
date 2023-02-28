#pragma once
#ifndef LI_CrossSection_H
#define LI_CrossSection_H

#include <map>
#include <set>
#include <array>
#include <memory>
#include <string>
#include <stdexcept>

#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

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

#include "LeptonInjector/serialization/array.h"
#include "LeptonInjector/utilities/Interpolator.h"

#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

namespace LI {
namespace utilities {
class LI_random;
}
}

namespace LI {
namespace crosssections {

class CrossSection {
friend cereal::access;
private:
public:
    CrossSection();
    bool operator==(CrossSection const & other) const;
    virtual bool equal(CrossSection const & other) const = 0;
    virtual double TotalCrossSection(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const = 0;
    virtual double DifferentialCrossSection(dataclasses::InteractionRecord const &) const = 0;
    virtual double InteractionThreshold(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const = 0;

    virtual std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const = 0;
    virtual std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;

    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};
};

class CrossSectionCollection {
private:
    LI::dataclasses::Particle::ParticleType primary_type;
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> cross_sections_by_target;
    std::set<LI::dataclasses::Particle::ParticleType> target_types;
    static const std::vector<std::shared_ptr<CrossSection>> empty;
    void InitializeTargetTypes();
public:
    CrossSectionCollection();
    CrossSectionCollection(LI::dataclasses::Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections);
    bool operator==(CrossSectionCollection const & other) const;
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSections() const {return cross_sections;};
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSectionsForTarget(LI::dataclasses::Particle::ParticleType p) const;
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & GetCrossSectionsByTarget() const {
        return cross_sections_by_target;
    };
    std::set<LI::dataclasses::Particle::ParticleType> const & TargetTypes() const {
        return target_types;
    };
    virtual bool MatchesPrimary(dataclasses::InteractionRecord const & record) const;
public:
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("CrossSections", cross_sections));
        } else {
            throw std::runtime_error("CrossSectionCollection only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("CrossSections", cross_sections));
        } else {
            throw std::runtime_error("CrossSectionCollection only supports version <= 0!");
        }
    }
};

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

public:
    DISFromSpline();
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::set<LI::dataclasses::Particle::ParticleType> primary_types, std::set<LI::dataclasses::Particle::ParticleType> target_types);
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::vector<LI::dataclasses::Particle::ParticleType> primary_types, std::vector<LI::dataclasses::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LI::dataclasses::Particle::ParticleType> primary_types, std::set<LI::dataclasses::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LI::dataclasses::Particle::ParticleType> primary_types, std::set<LI::dataclasses::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LI::dataclasses::Particle::ParticleType> primary_types, std::vector<LI::dataclasses::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LI::dataclasses::Particle::ParticleType> primary_types, std::vector<LI::dataclasses::Particle::ParticleType> target_types);

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy) const;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const;
    double DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const;
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random> random) const;

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const;

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
            LoadFromMemory(differential_data, total_data);
            InitializeSignatures();
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DISFromSpline only supports version <= 0!");
        }
    }
private:
    void ReadParamsFromSplineTable();
    void InitializeSignatures();
};

// For details, see appendix A of 1906.00111v4
class ElasticScattering : public CrossSection {
friend cereal::access;
protected:
private:
		const double CLR = 0.2334; // at one loop
    const std::set<LI::dataclasses::Particle::ParticleType> primary_types = {LI::dataclasses::Particle::ParticleType::NuE, LI::dataclasses::Particle::ParticleType::NuMu};
public:
		ElasticScattering() {};
		virtual bool equal(CrossSection const & other) const override;
		double DifferentialCrossSection(dataclasses::InteractionRecord const &) const;
    double DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, double y) const;
    double TotalCrossSection(dataclasses::InteractionRecord const &) const;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const;
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const;

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const;
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("ElasticScattering only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            std::set<LI::dataclasses::Particle::ParticleType> prim;
            archive(::cereal::make_nvp("PrimaryTypes", prim));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("ElasticScattering only supports version <= 0!");
        }
    }
};

class DipoleFromTable : public CrossSection {
friend cereal::access;
protected:
DipoleFromTable() {};
public:
    enum HelicityChannel {Conserving, Flipping};
private:
    bool z_samp = true;
    bool in_invGeV = true;
    bool inelastic = true;
    std::map<LI::dataclasses::Particle::ParticleType, LI::utilities::Interpolator2D<double>> differential;
    std::map<LI::dataclasses::Particle::ParticleType, LI::utilities::Interpolator1D<double>> total;
    const std::set<LI::dataclasses::Particle::ParticleType> primary_types = {LI::dataclasses::Particle::ParticleType::NuE, LI::dataclasses::Particle::ParticleType::NuMu, LI::dataclasses::Particle::ParticleType::NuTau, LI::dataclasses::Particle::ParticleType::NuEBar, LI::dataclasses::Particle::ParticleType::NuMuBar, LI::dataclasses::Particle::ParticleType::NuTauBar};
    double hnl_mass;
    double dipole_coupling;
    HelicityChannel channel;
public:
    virtual bool equal(CrossSection const & other) const override;
    double GetHNLMass() const {return hnl_mass;};
    static double DipoleyMin(double Enu, double mHNL, double target_mass);
    static double DipoleyMax(double Enu, double mHNL, double target_mass);
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV) : z_samp(z_samp), in_invGeV(in_invGeV), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV, bool inelastic) : z_samp(z_samp), in_invGeV(in_invGeV), inelastic(inelastic), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, std::set<LI::dataclasses::Particle::ParticleType> const & primary_types) : primary_types(primary_types), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV, std::set<LI::dataclasses::Particle::ParticleType> const & primary_types) : z_samp(z_samp), in_invGeV(in_invGeV), primary_types(primary_types), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV, bool inelastic, std::set<LI::dataclasses::Particle::ParticleType> const & primary_types) : z_samp(z_samp), in_invGeV(in_invGeV), inelastic(inelastic), primary_types(primary_types), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    double TotalCrossSection(dataclasses::InteractionRecord const &) const;
    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const;
    double DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type, double target_mass, double y) const;
    double DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type, double target_mass, double y, double thresh) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const;
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const;

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const;
    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const;

    void AddDifferentialCrossSectionFile(std::string filename, LI::dataclasses::Particle::ParticleType target);
    void AddTotalCrossSectionFile(std::string filename, LI::dataclasses::Particle::ParticleType target);
    void AddDifferentialCrossSection(LI::dataclasses::Particle::ParticleType target, LI::utilities::Interpolator2D<double>);
    void AddTotalCrossSection(LI::dataclasses::Particle::ParticleType target, LI::utilities::Interpolator1D<double>);
public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("DifferentialCrossSection", differential));
            archive(::cereal::make_nvp("TotalCrossSection", total));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("HelicityChannel", static_cast<int>(channel)));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DifferentialCrossSection", differential));
            archive(::cereal::make_nvp("TotalCrossSection", total));
            std::set<LI::dataclasses::Particle::ParticleType> prim;
            archive(::cereal::make_nvp("PrimaryTypes", prim));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("HelicityChannel", channel));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
};

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::crosssections::CrossSection, 0);

CEREAL_CLASS_VERSION(LI::crosssections::DISFromSpline, 0);
CEREAL_REGISTER_TYPE(LI::crosssections::DISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::crosssections::CrossSection, LI::crosssections::DISFromSpline);

CEREAL_CLASS_VERSION(LI::crosssections::DipoleFromTable, 0);
CEREAL_REGISTER_TYPE(LI::crosssections::DipoleFromTable);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::crosssections::CrossSection, LI::crosssections::DipoleFromTable);

CEREAL_CLASS_VERSION(LI::crosssections::CrossSectionCollection, 0);

#endif // LI_CrossSection_H

