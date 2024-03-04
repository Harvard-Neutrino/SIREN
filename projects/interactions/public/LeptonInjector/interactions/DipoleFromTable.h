#pragma once
#ifndef LI_DipoleFromTable_H
#define LI_DipoleFromTable_H

#include <map>                                          // for map
#include <set>                                          // for set
#include <memory>                                       // for shared_ptr
#include <string>                                       // for string
#include <vector>                                       // for vector
#include <cstdint>                                      // for uint32_t
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/interactions/CrossSection.h"  // for CrossSection
#include "LeptonInjector/dataclasses/Particle.h"        // for Particle
#include "LeptonInjector/utilities/Interpolator.h"      // for Interpolator1D

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace LI { namespace dataclasses { struct InteractionSignature; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace interactions {

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
    std::map<LI::dataclasses::ParticleType, LI::utilities::Interpolator2D<double>> differential;
    std::map<LI::dataclasses::ParticleType, LI::utilities::Interpolator1D<double>> total;
    const std::set<LI::dataclasses::ParticleType> primary_types = {LI::dataclasses::ParticleType::NuE, LI::dataclasses::ParticleType::NuMu, LI::dataclasses::ParticleType::NuTau, LI::dataclasses::ParticleType::NuEBar, LI::dataclasses::ParticleType::NuMuBar, LI::dataclasses::ParticleType::NuTauBar};
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
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, std::set<LI::dataclasses::ParticleType> const & primary_types) : primary_types(primary_types), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV, std::set<LI::dataclasses::ParticleType> const & primary_types) : z_samp(z_samp), in_invGeV(in_invGeV), primary_types(primary_types), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    DipoleFromTable(double hnl_mass, double dipole_coupling, HelicityChannel channel, bool z_samp, bool in_invGeV, bool inelastic, std::set<LI::dataclasses::ParticleType> const & primary_types) : z_samp(z_samp), in_invGeV(in_invGeV), inelastic(inelastic), primary_types(primary_types), hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), channel(channel) {};
    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(LI::dataclasses::ParticleType primary, double energy, LI::dataclasses::ParticleType target) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double DifferentialCrossSection(LI::dataclasses::ParticleType primary_type, double primary_energy, LI::dataclasses::ParticleType target_type, double target_mass, double y) const;
    double DifferentialCrossSection(LI::dataclasses::ParticleType primary_type, double primary_energy, LI::dataclasses::ParticleType target_type, double target_mass, double y, double thresh) const;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<LI::utilities::LI_random>) const override;

    std::vector<LI::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<LI::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::ParticleType primary_type) const override;
    std::vector<LI::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::ParticleType primary_type, LI::dataclasses::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    void AddDifferentialCrossSectionFile(std::string filename, LI::dataclasses::ParticleType target);
    void AddTotalCrossSectionFile(std::string filename, LI::dataclasses::ParticleType target);
    void AddDifferentialCrossSection(LI::dataclasses::ParticleType target, LI::utilities::Interpolator2D<double>);
    void AddTotalCrossSection(LI::dataclasses::ParticleType target, LI::utilities::Interpolator1D<double>);
public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("ZSampling", z_samp));
            archive(::cereal::make_nvp("InvGeVUnits", in_invGeV));
            archive(::cereal::make_nvp("Inelastic", inelastic));
            archive(::cereal::make_nvp("DifferentialCrossSection", differential));
            archive(::cereal::make_nvp("TotalCrossSection", total));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("DipoleCoupling", dipole_coupling));
            archive(::cereal::make_nvp("HelicityChannel", static_cast<int>(channel)));
            archive(::cereal::make_nvp("CrossSection", cereal::virtual_base_class<CrossSection>(this)));
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<DipoleFromTable> & construct, std::uint32_t version) {
        if(version == 0) {
            bool _z_samp = true;
            bool _in_invGeV = true;
            bool _inelastic = true;
            std::map<LI::dataclasses::ParticleType, LI::utilities::Interpolator2D<double>> _differential;
            std::map<LI::dataclasses::ParticleType, LI::utilities::Interpolator1D<double>> _total;
            std::set<LI::dataclasses::ParticleType> _primary_types;
            double _hnl_mass;
            double _dipole_coupling;
            HelicityChannel _channel;

            archive(::cereal::make_nvp("ZSampling", _z_samp));
            archive(::cereal::make_nvp("InvGeVUnits", _in_invGeV));
            archive(::cereal::make_nvp("Inelastic", _inelastic));
            archive(::cereal::make_nvp("DifferentialCrossSection", _differential));
            archive(::cereal::make_nvp("TotalCrossSection", _total));
            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            archive(::cereal::make_nvp("HNLMass", _hnl_mass));
            archive(::cereal::make_nvp("DipoleCoupling", _dipole_coupling));
            archive(::cereal::make_nvp("HelicityChannel", _channel));
            construct(_hnl_mass, _dipole_coupling, _channel, _z_samp, _in_invGeV, _inelastic, _primary_types);
            archive(::cereal::make_nvp("CrossSection", cereal::virtual_base_class<CrossSection>(construct.ptr())));
            construct.ptr()->differential = _differential;
            construct.ptr()->total = _total;
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace LI

CEREAL_CLASS_VERSION(LI::interactions::DipoleFromTable, 0);
CEREAL_REGISTER_TYPE(LI::interactions::DipoleFromTable);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::interactions::CrossSection, LI::interactions::DipoleFromTable);

#endif // LI_DipoleFromTable_H
