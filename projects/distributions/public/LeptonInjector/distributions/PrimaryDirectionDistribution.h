#pragma once
#ifndef LI_PrimaryDirectionDistribution_H
#define LI_PrimaryDirectionDistribution_H

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

class PrimaryDirectionDistribution : virtual public InjectionDistribution {
friend cereal::access;
protected:
    PrimaryDirectionDistribution() {};
private:
    virtual LI::math::Vector3D SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class IsotropicDirection : virtual public PrimaryDirectionDistribution {
friend cereal::access;
private:
    LI::math::Vector3D SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("IsotropicDirection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("IsotropicDirection only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class FixedDirection : virtual public PrimaryDirectionDistribution {
friend cereal::access;
protected:
    FixedDirection() {};
private:
    LI::math::Vector3D dir;
public:
    FixedDirection(LI::math::Vector3D dir) : dir(dir) {};
private:
    LI::math::Vector3D SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Direction", dir));
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("FixedDirection only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<FixedDirection> & construct, std::uint32_t const version) {
        if(version == 0) {
            LI::math::Vector3D d;
            archive(::cereal::make_nvp("Direction", d));
            construct(d);
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("FixedDirection only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class Cone : virtual public PrimaryDirectionDistribution {
friend cereal::access;
protected:
    Cone() {};
private:
    LI::math::Vector3D dir;
    LI::math::Quaternion rotation;
    double opening_angle;
public:
    Cone(LI::math::Vector3D dir, double opening_angle);
private:
    LI::math::Vector3D SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Direction", dir));
            archive(::cereal::make_nvp("OpeningAngle", opening_angle));
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("Cone only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<Cone> & construct, std::uint32_t const version) {
        if(version == 0) {
            LI::math::Vector3D d;
            double angle;
            archive(::cereal::make_nvp("Direction", d));
            archive(::cereal::make_nvp("OpeningAngle", angle));
            construct(d, angle);
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("Cone only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::PrimaryDirectionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::PrimaryDirectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::InjectionDistribution, LI::distributions::PrimaryDirectionDistribution);

CEREAL_CLASS_VERSION(LI::distributions::IsotropicDirection, 0);
CEREAL_REGISTER_TYPE(LI::distributions::IsotropicDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryDirectionDistribution, LI::distributions::IsotropicDirection);

CEREAL_CLASS_VERSION(LI::distributions::FixedDirection, 0);
CEREAL_REGISTER_TYPE(LI::distributions::FixedDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryDirectionDistribution, LI::distributions::FixedDirection);

CEREAL_CLASS_VERSION(LI::distributions::Cone, 0);
CEREAL_REGISTER_TYPE(LI::distributions::Cone);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::PrimaryDirectionDistribution, LI::distributions::Cone);

#endif // LI_PrimaryDirectionDistribution_H
