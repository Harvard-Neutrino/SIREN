#pragma once
#ifndef SIREN_HNLDecay_H
#define SIREN_HNLDecay_H

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <math>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>


#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"

#include "SIREN/interactions/Decay.h"

namespace siren {
namespace interactions {

class HNLDecay : public Decay {
friend cereal::access;
protected:
HNLDecay() {};
public:
    enum ChiralNature {Dirac, Majorana};
private:
    double hnl_mass;
    std::vector<double> mixing; // Ue4, Um4, Ut4
    ChiralNature nature;
    const std::set<siren::dataclasses::ParticleType> primary_types = {siren::dataclasses::ParticleType::NuF4, siren::dataclasses::ParticleType::NuF4Bar};
public:
    HNLDecay(double hnl_mass, std::vector<double> mixing, ChiralNature nature) : hnl_mass(hnl_mass), mixing(mixing), nature(nature) {};
    HNLDecay(double hnl_mass, std::vector<double> mixing, ChiralNature nature, std::set<siren::dataclasses::ParticleType> const & primary_types) : hnl_mass(hnl_mass), mixing(mixing), nature(nature), primary_types(primary_types) {};
    virtual bool equal(Decay const & other) const override;
    double GetHNLMass() const {return hnl_mass;};
    // if only one dipole coupling provided, assume it is U4t
    HNLDecay(double hnl_mass, double mixing, ChiralNature nature) : hnl_mass(hnl_mass), mixing(std::vector<double>{0,0,mixing}), nature(nature) {};
    HNLDecay(double hnl_mass, double mixing, ChiralNature nature, std::set<siren::dataclasses::ParticleType> const & primary_types) : hnl_mass(hnl_mass), mixing(std::vector<double>{0,0,mixing}), nature(nature), primary_types(primary_types) {};
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("DipoleCoupling", mixing));
            archive(::cereal::make_nvp("ChiralNature", static_cast<int>(nature)));
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("HNLDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<HNLDecay> & construct, std::uint32_t version) {
        if(version == 0) {
            std::set<siren::dataclasses::ParticleType> _primary_types;
            double _hnl_mass;
            double _mixing;
            ChiralNature _nature;

            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            archive(::cereal::make_nvp("HNLMass", _hnl_mass));
            archive(::cereal::make_nvp("DipoleCoupling", _mixing));
            archive(::cereal::make_nvp("ChiralNature", _nature));
            construct(_hnl_mass, _mixing, _nature, _primary_types);
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(construct.ptr())));
        } else {
            throw std::runtime_error("HNLDecay only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::HNLDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::HNLDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::HNLDecay);

#endif // SIREN_HNLDecay_H
