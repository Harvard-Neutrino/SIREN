#pragma once
#ifndef LI_Distributions_H
#define LI_Distributions_H

#include <cereal/access.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include "serialization/array.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"
#include "LeptonInjector/Constants.h"
#include "LeptonInjector/DataWriter.h"
#include "LeptonInjector/EventProps.h"
#include "LeptonInjector/Coordinates.h"
#include "LeptonInjector/BasicInjectionConfiguration.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"

namespace LeptonInjector {

class WeightableDistribution {
friend cereal::access;
public:
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const;
    virtual std::string Name() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            throw std::runtime_error("WeightableDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("WeightableDistribution only supports version <= 0!");
        }
    }
    bool operator==(WeightableDistribution const & distribution) const;
    bool operator<(WeightableDistribution const & distribution) const;
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class InjectionDistribution : public WeightableDistribution {
friend cereal::access;
private:
public:
    virtual void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("InjectionDistribution only supports version <= 0!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<WeightableDistribution>(this));
        } else {
            throw std::runtime_error("InjectionDistribution only supports version <= 0!");
        }
    }
};

class PrimaryInjector : public InjectionDistribution {
friend cereal::access;
protected:
    PrimaryInjector() {};
private:
    LeptonInjector::Particle::ParticleType primary_type;
    double primary_mass;
public:
    PrimaryInjector(LeptonInjector::Particle::ParticleType primary_type, double primary_mass = 0);
    LeptonInjector::Particle::ParticleType PrimaryType() const;
    double PrimaryMass() const;
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
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
            LeptonInjector::Particle::ParticleType t;
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

class TargetMomentumDistribution : public InjectionDistribution {
friend cereal::access;
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("TargetMomentumDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("TargetMomentumDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class TargetAtRest : public TargetMomentumDistribution {
friend cereal::access;
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const;
public:
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    std::string Name() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<TargetMomentumDistribution>(this));
        } else {
            throw std::runtime_error("TargetAtRest only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<TargetMomentumDistribution>(this));
        } else {
            throw std::runtime_error("TargetAtRest only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class PrimaryEnergyDistribution : public InjectionDistribution {
friend cereal::access;
private:
    virtual double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const override;
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class PowerLaw : public PrimaryEnergyDistribution {
friend cereal::access;
protected:
    PowerLaw() {};
private:
    double powerLawIndex;
    double energyMin;
    double energyMax;
public:
    PowerLaw(double powerLawIndex, double energyMin, double energyMax);
    double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PowerLawIndex", powerLawIndex));
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("PowerLaw only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<PowerLaw> & construct, std::uint32_t const version) {
        if(version == 0) {
            double gamma, min, max;
            archive(::cereal::make_nvp("PowerLawIndex", gamma));
            archive(::cereal::make_nvp("EnergyMin", min));
            archive(::cereal::make_nvp("EnergyMax", max));
            construct(gamma, min, max);
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("PowerLaw only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class ModifiedMoyalPlusExponentialEnergyDistribution : public PrimaryEnergyDistribution {
friend cereal::access;
protected:
    ModifiedMoyalPlusExponentialEnergyDistribution() {};
private:
    double energyMin;
    double energyMax;
    double mu;
    double sigma;
    double A;
    double l;
    double B;
    double integral;
    const size_t burnin = 40;
    double unnormed_pdf(double energy) const ;
    double pdf(double energy) const;
public:
    double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    ModifiedMoyalPlusExponentialEnergyDistribution(double energyMin, double energyMax, double mu, double sigma, double A, double l, double B);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(::cereal::make_nvp("ParameterMu", mu));
            archive(::cereal::make_nvp("ParameterSigma", sigma));
            archive(::cereal::make_nvp("ParameterA", A));
            archive(::cereal::make_nvp("ParameterL", l));
            archive(::cereal::make_nvp("ParameteB", B));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("ModifiedMoyalPlusExponentialEnergyDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<ModifiedMoyalPlusExponentialEnergyDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double min, max, mu, s, a, l, b;
            archive(::cereal::make_nvp("EnergyMin", min));
            archive(::cereal::make_nvp("EnergyMax", max));
            archive(::cereal::make_nvp("ParameterMu", mu));
            archive(::cereal::make_nvp("ParameterSigma", s));
            archive(::cereal::make_nvp("ParameterA", a));
            archive(::cereal::make_nvp("ParameterL", l));
            archive(::cereal::make_nvp("ParameteB", b));
            construct(min, max, mu, s, a, l, b);
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("ModifiedMoyalPlusExponentialEnergyDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class PrimaryDirectionDistribution : public InjectionDistribution {
friend cereal::access;
protected:
    PrimaryDirectionDistribution() {};
private:
    virtual earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
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

class IsotropicDirection : public PrimaryDirectionDistribution {
friend cereal::access;
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
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

class FixedDirection : public PrimaryDirectionDistribution {
friend cereal::access;
protected:
    FixedDirection() {};
private:
    earthmodel::Vector3D dir;
public:
    FixedDirection(earthmodel::Vector3D dir) : dir(dir) {};
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
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
            earthmodel::Vector3D d;
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

class Cone : public PrimaryDirectionDistribution {
friend cereal::access;
protected:
    Cone() {};
private:
    earthmodel::Vector3D dir;
    earthmodel::Quaternion rotation;
    double opening_angle;
public:
    Cone(earthmodel::Vector3D dir, double opening_angle);
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
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
            earthmodel::Vector3D d;
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

class VertexPositionDistribution : public InjectionDistribution {
friend cereal::access;
private:
    virtual earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const;
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    virtual std::pair<earthmodel::Vector3D, earthmodel::Vector3D> InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & interaction) const = 0;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("VertexPositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("VertexPositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override = 0;
    virtual bool less(WeightableDistribution const & distribution) const = 0;
};

class CylinderVolumePositionDistribution : public VertexPositionDistribution {
friend cereal::access;
protected:
    CylinderVolumePositionDistribution() {};
private:
    earthmodel::Cylinder cylinder;
    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    CylinderVolumePositionDistribution(earthmodel::Cylinder);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    virtual std::pair<earthmodel::Vector3D, earthmodel::Vector3D> InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Cylinder", cylinder));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("CylinderVolumePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<LeptonInjector::CylinderVolumePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            earthmodel::Cylinder c;
            archive(::cereal::make_nvp("Cylinder", c));
            construct(c);
            archive(cereal::virtual_base_class<LeptonInjector::VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("CylinderVolumePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class DepthFunction {
friend cereal::access;
public:
    DepthFunction();
    virtual double operator()(InteractionSignature const & signature, double energy) const;
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

class RangeFunction {
friend cereal::access;
public:
    RangeFunction();
    virtual double operator()(InteractionSignature const & signature, double energy) const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
        } else {
            throw std::runtime_error("RangeFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("RangeFunction only supports version <= 0!");
        }
    }
    bool operator==(RangeFunction const & distribution) const;
    bool operator<(RangeFunction const & distribution) const;
protected:
    virtual bool equal(RangeFunction const & distribution) const = 0;
    virtual bool less(RangeFunction const & distribution) const = 0;
};

class DecayRangeFunction : public RangeFunction {
friend cereal::access;
protected:
    DecayRangeFunction() {};
private:
    double particle_mass; // GeV
    double decay_width; // GeV
    double multiplier;
public:
    DecayRangeFunction(double particle_mass, double decay_width, double multiplier);
    static double DecayLength(double mass, double width, double energy);
    double operator()(InteractionSignature const & signature, double energy) const override;
    double DecayLength(InteractionSignature const & signature, double energy) const;
    double Range(InteractionSignature const & signature, double energy) const;
    double Multiplier() const;
    double ParticleMass() const;
    double DecayWidth() const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("ParticleMass", particle_mass));
            archive(::cereal::make_nvp("DecayWidth", decay_width));
            archive(cereal::virtual_base_class<RangeFunction>(this));
        } else {
            throw std::runtime_error("DecayRangeFunction only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<DecayRangeFunction> & construct, std::uint32_t const version) {
        if(version == 0) {
            double mass;
            double width;
            double multiplier;
            archive(::cereal::make_nvp("ParticleMass", mass));
            archive(::cereal::make_nvp("DecayWidth", width));
            archive(::cereal::make_nvp("Multiplier", multiplier));
            construct(mass, width, multiplier);
            archive(cereal::virtual_base_class<RangeFunction>(construct.ptr()));
        } else {
            throw std::runtime_error("DecayRangeFunction only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(RangeFunction const & distribution) const override;
    virtual bool less(RangeFunction const & distribution) const override;
};

class ColumnDepthPositionDistribution : public VertexPositionDistribution {
friend cereal::access;
protected:
    ColumnDepthPositionDistribution() {};
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DepthFunction> depth_function;
    std::set<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const;

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::set<Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    virtual std::pair<earthmodel::Vector3D, earthmodel::Vector3D> InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & interaction) const override;
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
            std::set<Particle::ParticleType> t;
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

class RangePositionDistribution : public VertexPositionDistribution {
friend cereal::access;
private:
    double radius;
    double endcap_length;
    std::shared_ptr<RangeFunction> range_function;
    std::set<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const;

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    RangePositionDistribution();
    RangePositionDistribution(const RangePositionDistribution &) = default;
    RangePositionDistribution(double radius, double endcap_length, std::shared_ptr<RangeFunction> range_function, std::set<Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::pair<earthmodel::Vector3D, earthmodel::Vector3D> InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("RangeFunction", range_function));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("RangePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<RangePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::set<Particle::ParticleType> t;
            std::shared_ptr<RangeFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("RangeFunction", f));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, f, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("RangePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class DecayRangePositionDistribution : public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DecayRangeFunction> range_function;
    std::set<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const;

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
public:
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
    DecayRangePositionDistribution();
    DecayRangePositionDistribution(const DecayRangePositionDistribution &) = default;
    DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function, std::set<Particle::ParticleType> target_types);
    std::string Name() const override;
    virtual std::pair<earthmodel::Vector3D, earthmodel::Vector3D> InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & interaction) const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Radius", radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("DecayRangeFunction", range_function));
            archive(::cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::virtual_base_class<VertexPositionDistribution>(this));
        } else {
            throw std::runtime_error("DecayRangePositionDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<DecayRangePositionDistribution> & construct, std::uint32_t const version) {
        if(version == 0) {
            double r;
            double l;
            std::set<Particle::ParticleType> t;
            std::shared_ptr<DecayRangeFunction> f;
            archive(::cereal::make_nvp("Radius", r));
            archive(::cereal::make_nvp("EndcapLength", l));
            archive(::cereal::make_nvp("DecayRangeFunction", f));
            archive(::cereal::make_nvp("TargetTypes", t));
            construct(r, l, f, t);
            archive(cereal::virtual_base_class<VertexPositionDistribution>(construct.ptr()));
        } else {
            throw std::runtime_error("DecayRangePositionDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

class PrimaryNeutrinoHelicityDistribution : public InjectionDistribution {
friend cereal::access;
public:
    virtual void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override;
    virtual double GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override;
	PrimaryNeutrinoHelicityDistribution();
	PrimaryNeutrinoHelicityDistribution(const PrimaryNeutrinoHelicityDistribution &) = default;
    virtual std::vector<std::string> DensityVariables() const;
    std::string Name() const override;
    virtual std::shared_ptr<InjectionDistribution> clone() const;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryNeutrinoHelicityDistribution only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryNeutrinoHelicityDistribution only supports version <= 0!");
        }
    }
protected:
    virtual bool equal(WeightableDistribution const & distribution) const override;
    virtual bool less(WeightableDistribution const & distribution) const override;
};

} // namespace LeptonInjector

CEREAL_CLASS_VERSION(LeptonInjector::WeightableDistribution, 0);

CEREAL_CLASS_VERSION(LeptonInjector::InjectionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::InjectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::WeightableDistribution, LeptonInjector::InjectionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::PrimaryInjector, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::PrimaryInjector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::InjectionDistribution, LeptonInjector::PrimaryInjector);

CEREAL_CLASS_VERSION(LeptonInjector::TargetMomentumDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::TargetMomentumDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::InjectionDistribution, LeptonInjector::TargetMomentumDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::TargetAtRest, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::TargetAtRest);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::TargetMomentumDistribution, LeptonInjector::TargetAtRest);

CEREAL_CLASS_VERSION(LeptonInjector::PrimaryEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::PrimaryEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::InjectionDistribution, LeptonInjector::PrimaryEnergyDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::PowerLaw, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::PowerLaw);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::PrimaryEnergyDistribution, LeptonInjector::PowerLaw);

CEREAL_CLASS_VERSION(LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::PrimaryEnergyDistribution, LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::PrimaryDirectionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::PrimaryDirectionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::InjectionDistribution, LeptonInjector::PrimaryDirectionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::IsotropicDirection, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::IsotropicDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::PrimaryDirectionDistribution, LeptonInjector::IsotropicDirection);

CEREAL_CLASS_VERSION(LeptonInjector::FixedDirection, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::FixedDirection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::PrimaryDirectionDistribution, LeptonInjector::FixedDirection);

CEREAL_CLASS_VERSION(LeptonInjector::Cone, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::Cone);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::PrimaryDirectionDistribution, LeptonInjector::Cone);

CEREAL_CLASS_VERSION(LeptonInjector::VertexPositionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::VertexPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::InjectionDistribution, LeptonInjector::VertexPositionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::CylinderVolumePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::CylinderVolumePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::VertexPositionDistribution, LeptonInjector::CylinderVolumePositionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::DepthFunction, 0);

CEREAL_CLASS_VERSION(LeptonInjector::RangeFunction, 0);

CEREAL_CLASS_VERSION(LeptonInjector::DecayRangeFunction, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::DecayRangeFunction);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::RangeFunction, LeptonInjector::DecayRangeFunction);

CEREAL_CLASS_VERSION(LeptonInjector::ColumnDepthPositionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::ColumnDepthPositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::VertexPositionDistribution, LeptonInjector::ColumnDepthPositionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::RangePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::RangePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::VertexPositionDistribution, LeptonInjector::RangePositionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::DecayRangePositionDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::DecayRangePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::VertexPositionDistribution, LeptonInjector::DecayRangePositionDistribution);

CEREAL_CLASS_VERSION(LeptonInjector::PrimaryNeutrinoHelicityDistribution, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::PrimaryNeutrinoHelicityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::InjectionDistribution, LeptonInjector::PrimaryNeutrinoHelicityDistribution);

#endif // LI_Distributions_H


