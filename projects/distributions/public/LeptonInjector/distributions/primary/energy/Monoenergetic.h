#pragma once
#ifndef LI_Monoenergetic_H
#define LI_Monoenergetic_H

#include <memory>
#include <string>
#include <cstdint>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace distributions {

class Monoenergetic : virtual public PrimaryEnergyDistribution {
friend cereal::access;
protected:
    Monoenergetic() {};
private:
    double gen_energy;
public:
    Monoenergetic(double gen_energy);
    double pdf(double energy) const;
    double SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const override;
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("GenEnergy", gen_energy));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("Monoenergetic only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<Monoenergetic> & construct, std::uint32_t const version) {
        if(version == 0) {
            double energy;
            archive(::cereal::make_nvp("GenEnergy", energy));
            construct(energy);
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("Monoenergetic only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::Monoenergetic, 0);
CEREAL_REGISTER_TYPE(LI::distributions::Monoenergetic);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryEnergyDistribution, LI::distributions::Monoenergetic);

#endif // LI_Monoenergetic_H
