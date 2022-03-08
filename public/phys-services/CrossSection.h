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
#include "serialization/array.h"

#include "phys-services/Interpolator.h"

#include "LeptonInjector/Particle.h"

namespace LeptonInjector {
// #include "LeptonInjector/Random.h"
class LI_random;
}

namespace LeptonInjector {

struct InteractionSignature {
    LeptonInjector::Particle::ParticleType primary_type;
    LeptonInjector::Particle::ParticleType target_type;
    std::vector<LeptonInjector::Particle::ParticleType> secondary_types;
    bool operator==(InteractionSignature const & other) const;
    bool operator<(InteractionSignature const & other) const;
    friend std::ostream& operator<<(std::ostream& os, InteractionSignature const& signature);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("TargetType", target_type));
            archive(cereal::make_nvp("SecondaryTypes", secondary_types));
        } else {
            throw std::runtime_error("InteractionSignature only supports version <= 0!");
        }
    }
};

struct InteractionRecord {
    InteractionSignature signature;
    double primary_mass = 0;
    std::array<double, 4> primary_momentum = {0, 0, 0, 0};
    double primary_helicity = 0;
    double target_mass = 0;
    std::array<double, 4> target_momentum = {0, 0, 0, 0};
    double target_helicity = 0;
    std::array<double, 3> interaction_vertex = {0, 0, 0};
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::vector<double> interaction_parameters;
    bool operator==(InteractionRecord const & other) const;
    friend std::ostream& operator<<(std::ostream& os, InteractionRecord const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InteractionSignature", signature));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("TargetMass", target_mass));
            archive(::cereal::make_nvp("TargetMomentum", target_momentum));
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("InteractionParameters", interaction_parameters));
        } else {
            throw std::runtime_error("InteractionRecord only supports version <= 0!");
        }
    };
};

struct DecaySignature {
    LeptonInjector::Particle::ParticleType primary_type;
    std::vector<LeptonInjector::Particle::ParticleType> secondary_types;
    bool operator==(DecaySignature const & other) const;
    friend std::ostream& operator<<(std::ostream& os, DecaySignature const& signature);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("SecondaryTypes", secondary_types));
        } else {
            throw std::runtime_error("DecaySignature only supports version <= 0!");
        }
    }
};

struct DecayRecord {
    DecaySignature signature;
    double primary_mass;
    std::array<double, 4> primary_momentum;
    double primary_helicity;
    std::array<double, 3> decay_vertex = {0, 0, 0};
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::vector<double> decay_parameters;
    bool operator==(DecayRecord const & other) const;
    friend std::ostream& operator<<(std::ostream& os, DecayRecord const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DecaySignature", signature));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("DecayVertex", decay_vertex));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("DecayParameters", decay_parameters));
        } else {
            throw std::runtime_error("DecayRecord only supports version <= 0!");
        }
    };
};

class CrossSection {
friend cereal::access;
private:
public:
    CrossSection();
    bool operator==(CrossSection const & other) const;
    virtual bool equal(CrossSection const & other) const = 0;
    virtual double TotalCrossSection(InteractionRecord const &) const = 0;
    virtual double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy, Particle::ParticleType target) const = 0;
    virtual double DifferentialCrossSection(InteractionRecord const &) const = 0;
    virtual double InteractionThreshold(InteractionRecord const &) const = 0;
    virtual void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random>) const = 0;

    virtual std::vector<Particle::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<Particle::ParticleType> GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const = 0;
    virtual std::vector<Particle::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<InteractionSignature> GetPossibleSignatures() const = 0;

    virtual std::vector<InteractionSignature> GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const = 0;
    virtual double FinalStateProbability(InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};
};

class CrossSectionCollection {
private:
    Particle::ParticleType primary_type;
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> cross_sections_by_target;
    std::set<Particle::ParticleType> target_types;
    static const std::vector<std::shared_ptr<CrossSection>> empty;
    void InitializeTargetTypes();
public:
    CrossSectionCollection();
    CrossSectionCollection(Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections);
    bool operator==(CrossSectionCollection const & other) const;
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSections() const {return cross_sections;};
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSectionsForTarget(Particle::ParticleType p) const;
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & GetCrossSectionsByTarget() const {
        return cross_sections_by_target;
    };
    std::set<Particle::ParticleType> const & TargetTypes() const {
        return target_types;
    };
    virtual bool MatchesPrimary(InteractionRecord const & record) const;
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

    std::vector<InteractionSignature> signatures_;
    std::set<LeptonInjector::Particle::ParticleType> primary_types_;
    std::set<LeptonInjector::Particle::ParticleType> target_types_;
    std::map<LeptonInjector::Particle::ParticleType, std::vector<LeptonInjector::Particle::ParticleType>> targets_by_primary_types_;
    std::map<std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>, std::vector<InteractionSignature>> signatures_by_parent_types_;

    int interaction_type_;
    double target_mass_;
    double minimum_Q2_;

public:
    DISFromSpline();
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy, Particle::ParticleType target) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    double DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass) const;
    double InteractionThreshold(InteractionRecord const &) const;
    void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random> random) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;
    std::vector<InteractionSignature> GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(InteractionRecord const & record) const;

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


class DipoleFromTable : public CrossSection {
friend cereal::access;
protected:
DipoleFromTable() {};
public:
    enum HelicityChannel {Conserving, Flipping};
private:
    bool z_samp = true;
    bool in_invGeV = true;
    std::map<Particle::ParticleType, Interpolator2D<double>> differential;
    std::map<Particle::ParticleType, Interpolator1D<double>> total;
    const std::set<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau, Particle::ParticleType::NuEBar, Particle::ParticleType::NuMuBar, Particle::ParticleType::NuTauBar};
    double hnl_mass;
    double dipole_coupling;
    HelicityChannel channel;
public:
    virtual bool equal(CrossSection const & other) const override;
    double GetHNLMass() const {return hnl_mass;};
    static double DipoleyMin(double Enu, double mHNL, double target_mass);
    static double DipoleyMax(double Enu, double mHNL, double target_mass);
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel), z_samp(z_samp), in_invGeV(in_invGeV) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, std::set<Particle::ParticleType> const & primary_types) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel), primary_types(primary_types) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV, std::set<Particle::ParticleType> const & primary_types) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel), z_samp(z_samp), in_invGeV(in_invGeV), primary_types(primary_types) {};
    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy, Particle::ParticleType target) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    double DifferentialCrossSection(Particle::ParticleType primary_type, double primary_energy, Particle::ParticleType target_type, double target_mass, double y) const;
    double DifferentialCrossSection(Particle::ParticleType primary_type, double primary_energy, Particle::ParticleType target_type, double target_mass, double y, double thresh) const;
    double InteractionThreshold(InteractionRecord const &) const;
    void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random>) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;
    std::vector<InteractionSignature> GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(InteractionRecord const & record) const;

    void AddDifferentialCrossSectionFile(std::string filename, Particle::ParticleType target);
    void AddTotalCrossSectionFile(std::string filename, Particle::ParticleType target);
    void AddDifferentialCrossSection(Particle::ParticleType target, Interpolator2D<double>);
    void AddTotalCrossSection(Particle::ParticleType target, Interpolator1D<double>);
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
            std::set<LeptonInjector::Particle::ParticleType> prim;
            archive(::cereal::make_nvp("PrimaryTypes", prim));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("HelicityChannel", channel));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
};

} // namespace LeptonInjector

CEREAL_CLASS_VERSION(LeptonInjector::CrossSection, 0);

CEREAL_CLASS_VERSION(LeptonInjector::DISFromSpline, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::DISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::CrossSection, LeptonInjector::DISFromSpline);

CEREAL_CLASS_VERSION(LeptonInjector::DipoleFromTable, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::DipoleFromTable);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::CrossSection, LeptonInjector::DipoleFromTable);

CEREAL_CLASS_VERSION(LeptonInjector::InteractionSignature, 0);
CEREAL_CLASS_VERSION(LeptonInjector::InteractionRecord, 0);
CEREAL_CLASS_VERSION(LeptonInjector::CrossSectionCollection, 0);


#endif // LI_CrossSection_H

