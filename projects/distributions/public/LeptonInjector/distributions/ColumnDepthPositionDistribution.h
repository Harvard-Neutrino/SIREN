#pragma once
#ifndef LI_ColumnDepthPositionDistribution_H
#define LI_ColumnDepthPositionDistribution_H

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
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/VertexPositionDistribution.h"

namespace LI {
namespace utilities {
class LI_random;
} // namespace utilities

namespace crosssections {
struct InteractionRecord;
struct InteractionSignature;
class CrossSectionCollection;
} // namespace crosssections
} // namespace LeptonInjector

namespace LI {
namespace distributions {

class DepthFunction {
friend cereal::access;
public:
    DepthFunction();
    virtual double operator()(LI::crosssections::InteractionSignature const & signature, double energy) const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            throw std::runtime_error("DepthFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("DepthFunction only supports version <= 0!");
        }
    }
    bool operator==(DepthFunction const & distribution) const;
    bool operator<(DepthFunction const & distribution) const;
protected:
    virtual bool equal(DepthFunction const & distribution) const = 0;
    virtual bool less(DepthFunction const & distribution) const = 0;
};

class LeptonDepthFunction : virtual public DepthFunction {
friend cereal::access;
private:
    double mu_alpha = 1.76666667e-3;
    double mu_beta = 2.0916666667e-6;
    double tau_alpha = 1.473684210526e1;
    double tau_beta = 2.6315789473684212e-7;
    double scale = 1.0;
    double max_depth = 3e7;
    std::set<LI::utilities::Particle::ParticleType> tau_primaries = {LI::utilities::Particle::ParticleType::NuTau, LI::utilities::Particle::ParticleType::NuTauBar};
public:
    LeptonDepthFunction();
    void SetMuParams(double mu_alpha, double mu_beta);
    void SetTauParams(double tau_alpha, double tau_beta);
    void SetScale(double scale);
    void SetMaxDepth(double max_depth);
    double GetMuAlpha() const;
    double GetMuBeta() const;
    double GetTauAlpha() const;
    double GetTauBeta() const;
    double GetScale() const;
    double GetMaxDepth() const;
    virtual double operator()(LI::crosssections::InteractionSignature const & signature, double energy) const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            archive(::cereal::make_nvp("MuAlpha", mu_alpha));
            archive(::cereal::make_nvp("MuBeta", mu_beta));
            archive(::cereal::make_nvp("TauAlpha", tau_alpha));
            archive(::cereal::make_nvp("TauBeta", tau_beta));
            archive(::cereal::make_nvp("Scale", scale));
            archive(::cereal::make_nvp("MaxDepth", max_depth));
            archive(::cereal::make_nvp("TauPrimaries", tau_primaries));
            throw std::runtime_error("LeptonDepthFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            archive(::cereal::make_nvp("MuAlpha", mu_alpha));
            archive(::cereal::make_nvp("MuBeta", mu_beta));
            archive(::cereal::make_nvp("TauAlpha", tau_alpha));
            archive(::cereal::make_nvp("TauBeta", tau_beta));
            archive(::cereal::make_nvp("Scale", scale));
            archive(::cereal::make_nvp("MaxDepth", max_depth));
            archive(::cereal::make_nvp("TauPrimaries", tau_primaries));
            throw std::runtime_error("LeptonDepthFunction only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(DepthFunction const & distribution) const override;
    virtual bool less(DepthFunction const & distribution) const override;
};

class ColumnDepthPositionDistribution : virtual public VertexPositionDistribution {
friend cereal::access;
protected:
    ColumnDepthPositionDistribution() {};
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DepthFunction> depth_function;
    std::set<LI::utilities::Particle::ParticleType> target_types;

    LI::math::Vector3D SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const;

    LI::math::Vector3D SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const override;
    ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::set<LI::utilities::Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("DepthFunction", depth_function));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("ColumnDepthPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<ColumnDepthPositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::set<LI::utilities::Particle::ParticleType> t;
            std::shared_ptr<DepthFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("DepthFunction", f));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, f, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("ColumnDepthPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace LI

CEREAL_CLASS_VERSION(LI::distributions::DepthFunction, 0);
CEREAL_CLASS_VERSION(LI::distributions::LeptonDepthFunction, 0);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::DepthFunction, LI::distributions::LeptonDepthFunction);

CEREAL_CLASS_VERSION(LI::distributions::ColumnDepthPositionDistribution, 0);
CEREAL_REGISTER_TYPE(LI::distributions::ColumnDepthPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::distributions::VertexPositionDistribution, LI::distributions::ColumnDepthPositionDistribution);

#endif // LI_ColumnDepthPositionDistribution_H
