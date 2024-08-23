#pragma once
#ifndef SIREN_ElectroweakDecay_H
#define SIREN_ElectroweakDecay_H

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>

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

class ElectroweakDecay : public Decay {
friend cereal::access;
protected:
ElectroweakDecay() {};
private:
    const std::set<siren::dataclasses::ParticleType> primary_types = {siren::dataclasses::ParticleType::WPlus,
                                                                      siren::dataclasses::ParticleType::WMinus,
                                                                      siren::dataclasses::ParticleType::Z0};
    std::vector<siren::dataclasses::ParticleType> Nus =              {siren::dataclasses::ParticleType::NuE,
                                                                      siren::dataclasses::ParticleType::NuMu,
                                                                      siren::dataclasses::ParticleType::NuTau};
    std::vector<siren::dataclasses::ParticleType> AntiNus =          {siren::dataclasses::ParticleType::NuEBar,
                                                                      siren::dataclasses::ParticleType::NuMuBar,
                                                                      siren::dataclasses::ParticleType::NuTauBar};
    std::vector<siren::dataclasses::ParticleType> Leptons =          {siren::dataclasses::ParticleType::EMinus,
                                                                      siren::dataclasses::ParticleType::MuMinus,
                                                                      siren::dataclasses::ParticleType::TauMinus};
    std::vector<siren::dataclasses::ParticleType> AntiLeptons =      {siren::dataclasses::ParticleType::EPlus,
                                                                      siren::dataclasses::ParticleType::MuPlus,
                                                                      siren::dataclasses::ParticleType::TauPlus};
    std::vector<siren::dataclasses::ParticleType> UpQuarks =         {siren::dataclasses::ParticleType::u,
                                                                      siren::dataclasses::ParticleType::c,
                                                                      siren::dataclasses::ParticleType::t};
    std::vector<siren::dataclasses::ParticleType> UpAntiQuarks =     {siren::dataclasses::ParticleType::uBar,
                                                                      siren::dataclasses::ParticleType::cBar,
                                                                      siren::dataclasses::ParticleType::tBar};
    std::vector<siren::dataclasses::ParticleType> DownQuarks =       {siren::dataclasses::ParticleType::d,
                                                                      siren::dataclasses::ParticleType::s,
                                                                      siren::dataclasses::ParticleType::b};
    std::vector<siren::dataclasses::ParticleType> DownAntiQuarks =   {siren::dataclasses::ParticleType::dBar,
                                                                      siren::dataclasses::ParticleType::sBar,
                                                                      siren::dataclasses::ParticleType::bBar};


public:
    ElectroweakDecay(std::set<siren::dataclasses::ParticleType> const & primary_types) :  primary_types(primary_types) {};
    virtual bool equal(Decay const & other) const override;
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
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("ElectroweakDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<ElectroweakDecay> & construct, std::uint32_t version) {
        if(version == 0) {
            std::set<siren::dataclasses::ParticleType> _primary_types;
            construct(_primary_types);
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(construct.ptr())));
        } else {
            throw std::runtime_error("ElectroweakDecay only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::ElectroweakDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::ElectroweakDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::ElectroweakDecay);

#endif // SIREN_ElectroweakDecay_H