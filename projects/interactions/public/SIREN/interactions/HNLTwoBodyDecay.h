#pragma once
#ifndef SIREN_HNLTwoBodyDecay_H
#define SIREN_HNLTwoBodyDecay_H

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>

#include <gsl/gsl_integration.h>

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

#include <CRunDec3.1/CRunDec.h>

namespace siren {
namespace interactions {

class HNLTwoBodyDecay : public Decay {
friend cereal::access;
protected:
HNLTwoBodyDecay() {crundec = new CRunDec();};
public:
    enum ChiralNature {Dirac, Majorana};
private:
    double hnl_mass;
    std::vector<double> mixing; // Ue4, Um4, Ut4
    ChiralNature nature;
    const std::set<siren::dataclasses::ParticleType> primary_types = {siren::dataclasses::ParticleType::N4, siren::dataclasses::ParticleType::N4Bar};
    CRunDec * crundec;
    double _GammaHadronsCC = 0;
    double _GammaHadronsNC = 0;
    std::vector<siren::dataclasses::ParticleType> MinusChargedMesons = {siren::dataclasses::ParticleType::PiMinus,
                                                                        siren::dataclasses::ParticleType::KMinus,
                                                                        siren::dataclasses::ParticleType::RhoMinus,
                                                                        siren::dataclasses::ParticleType::KPrimeMinus,
                                                                        siren::dataclasses::ParticleType::Hadrons,
                                                                        siren::dataclasses::ParticleType::DMinus,
                                                                        siren::dataclasses::ParticleType::DsMinus};
    std::vector<siren::dataclasses::ParticleType> PlusChargedMesons = {siren::dataclasses::ParticleType::PiPlus,
                                                                        siren::dataclasses::ParticleType::KPlus,
                                                                        siren::dataclasses::ParticleType::RhoPlus,
                                                                        siren::dataclasses::ParticleType::KPrimePlus,
                                                                        siren::dataclasses::ParticleType::DPlus,
                                                                        siren::dataclasses::ParticleType::DsPlus};
    std::vector<siren::dataclasses::ParticleType> NeutralMesons = {siren::dataclasses::ParticleType::Pi0,
                                                                   siren::dataclasses::ParticleType::Eta,
                                                                   siren::dataclasses::ParticleType::Rho0,
                                                                   siren::dataclasses::ParticleType::Omega,
                                                                   siren::dataclasses::ParticleType::EtaPrime,
                                                                   siren::dataclasses::ParticleType::KPrime0,
                                                                   siren::dataclasses::ParticleType::Phi};


public:
    HNLTwoBodyDecay(double hnl_mass, std::vector<double> mixing, ChiralNature nature) : hnl_mass(hnl_mass), mixing(mixing), nature(nature) {crundec = new CRunDec();};
    HNLTwoBodyDecay(double hnl_mass, std::vector<double> mixing, ChiralNature nature, std::set<siren::dataclasses::ParticleType> const & primary_types) : hnl_mass(hnl_mass), mixing(mixing), nature(nature), primary_types(primary_types) {crundec = new CRunDec();};
    HNLTwoBodyDecay(double hnl_mass, double mixing, ChiralNature nature) : hnl_mass(hnl_mass), mixing(std::vector<double>{0,0,mixing}), nature(nature) {crundec = new CRunDec();};
    HNLTwoBodyDecay(double hnl_mass, double mixing, ChiralNature nature, std::set<siren::dataclasses::ParticleType> const & primary_types) : hnl_mass(hnl_mass), mixing(std::vector<double>{0,0,mixing}), nature(nature), primary_types(primary_types) {crundec = new CRunDec();};
    virtual bool equal(Decay const & other) const override;
    double GetHNLMass() const {return hnl_mass;};
    // if only one coupling provided, assume it is U4t
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    double DeltaQCD(int nloops=5) const;
    double GetAlpha(dataclasses::ParticleType const & secondary) const;
    double GetMass(dataclasses::ParticleType const & secondary) const;
    void SetGammaHadrons(double Gamma, std::string mode);
    double CCMesonDecayWidth(dataclasses::InteractionRecord const &) const;
    double NCMesonDecayWidth(dataclasses::InteractionRecord const &) const;
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
            throw std::runtime_error("HNLTwoBodyDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<HNLTwoBodyDecay> & construct, std::uint32_t version) {
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
            throw std::runtime_error("HNLTwoBodyDecay only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::HNLTwoBodyDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::HNLTwoBodyDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::HNLTwoBodyDecay);

#endif // SIREN_HNLTwoBodyDecay_H
