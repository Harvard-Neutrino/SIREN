#pragma once
#ifndef LI_PrimaryInjector_H
#define LI_PrimaryInjector_H

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include <cereal/access.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/utilities/Particle.h"
#include "LeptonInjector/utilities/Interpolator.h"

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace detector {
class EarthModel;
} // namespace detector

namespace crosssections {
struct InteractionRecord;
struct InteractionSignature;
class CrossSectionCollection;
} // namespace crosssections
} // namespace LeptonInjector

namespace LI {
namespace distributions {

class PrimaryInjector : virtual public InjectionDistribution {
friend cereal::access;
protected:
    PrimaryInjector() {};
private:
    LI::utilities::Particle::ParticleType primary_type;
    double primary_mass;
public:
    PrimaryInjector(LI::utilities::Particle::ParticleType primary_type, double primary_mass = 0);
    LI::utilities::Particle::ParticleType PrimaryType() const;
    double PrimaryMass() const;
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryType", primary_type));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryInjector only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PrimaryInjector> & construct, std::uint32_t const version) {
        if(version == 0) {
            LI::utilities::Particle::ParticleType t;
            double m;
            archive(::cereal::make_nvp("PrimaryType", t));
            archive(::cereal::make_nvp("PrimaryMass", m));
            construct(t, m);
            archive(cereal::virtual_base_class<InjectionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PrimaryInjector only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PrimaryInjector, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryInjector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::PrimaryInjector);

#endif // LI_PrimaryInjector_H
