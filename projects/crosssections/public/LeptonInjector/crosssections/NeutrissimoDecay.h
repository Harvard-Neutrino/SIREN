#pragma once
#ifndef LI_NeutrissimoDecay_H
#define LI_NeutrissimoDecay_H

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>


#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include "LeptonInjector/crosssections/Decay.h"

namespace LI {
namespace crosssections {

class NeutrissimoDecay : public Decay {
friend cereal::access;
protected:
NeutrissimoDecay() {};
public:
    enum ChiralNature {Dirac, Majorana};
private:
    double hnl_mass;
    std::vector<double> dipole_coupling; // d_e, d_mu, d_tau
    ChiralNature nature;
    const std::set<LI::dataclasses::Particle::ParticleType> primary_types = {LI::dataclasses::Particle::ParticleType::NuF4, LI::dataclasses::Particle::ParticleType::NuF4Bar};
public:
    NeutrissimoDecay(double hnl_mass, std::vector<double> dipole_coupling, ChiralNature nature) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), nature(nature) {};
    NeutrissimoDecay(double hnl_mass, std::vector<double> dipole_coupling, ChiralNature nature, std::set<LI::dataclasses::Particle::ParticleType> const & primary_types) : hnl_mass(hnl_mass), dipole_coupling(dipole_coupling), nature(nature), primary_types(primary_types) {};
    virtual bool equal(Decay const & other) const override;
    double GetHNLMass() const {return hnl_mass;};
    // if only one dipole coupling provided, assume it is d_mu
    NeutrissimoDecay(double hnl_mass, double dipole_coupling, ChiralNature nature) : hnl_mass(hnl_mass), dipole_coupling(std::vector<double>{0,dipole_coupling,0}), nature(nature) {};
    NeutrissimoDecay(double hnl_mass, double dipole_coupling, ChiralNature nature, std::set<LI::dataclasses::Particle::ParticleType> const & primary_types) : hnl_mass(hnl_mass), dipole_coupling(std::vector<double>{0,dipole_coupling,0}), nature(nature), primary_types(primary_types) {};
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual double TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const override;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const override;
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    virtual std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("DipoleCoupling", dipole_coupling));
            archive(::cereal::make_nvp("ChiralNature", static_cast<int>(nature)));
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("NeutrissimoDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<NeutrissimoDecay> & construct, std::uint32_t version) {
        if(version == 0) {
            std::set<LI::dataclasses::Particle::ParticleType> _primary_types;
            double _hnl_mass;
            double _dipole_coupling;
            ChiralNature _nature;

            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            archive(::cereal::make_nvp("HNLMass", _hnl_mass));
            archive(::cereal::make_nvp("DipoleCoupling", _dipole_coupling));
            archive(::cereal::make_nvp("ChiralNature", _nature));
            construct(_hnl_mass, _dipole_coupling, _nature, _primary_types);
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(construct.ptr())));
        } else {
            throw std::runtime_error("NeutrissimoDecay only supports version <= 0!");
        }
    }

};

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::crosssections::NeutrissimoDecay, 0);
CEREAL_REGISTER_TYPE(LI::crosssections::NeutrissimoDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::crosssections::Decay, LI::crosssections::NeutrissimoDecay);

#endif // LI_NeutrissimoDecay_H
