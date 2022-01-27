#ifndef LI_LeptonInjector_H
#define LI_LeptonInjector_H

#include <queue>
#include <memory> // adds shared pointer
#include <iostream>

#include <photospline/bspline.h>
#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
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

#include "stga3/STGA3.h"
#include "stga3/Typedefs.h"
#include "stga3/Utilities.h"

namespace LeptonInjector {

class InjectionDistribution {
private:
public:
    virtual void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {};
    virtual std::vector<std::string> DensityVariables() const {return {};};
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("InjectionDistribution only supports version <= 0!");
        }
    }
};

class TargetMomentumDistribution : public InjectionDistribution {
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        record.target_momentum = SampleMomentum(rand, earth_model, cross_sections, record);
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"TargetMomentum"};};
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("TargetMomentumDistribution only supports version <= 0!");
        }
    }
};

class TargetAtRest : public TargetMomentumDistribution {
private:
public:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
        return std::array<double, 4>{record.target_mass, 0, 0, 0};
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>();};
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<TargetMomentumDistribution>(new TargetAtRest(*this));
    };
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<TargetMomentumDistribution>(this));
        } else {
            throw std::runtime_error("TargetAtRest only supports version <= 0!");
        }
    }
};

class PrimaryEnergyDistribution : public InjectionDistribution {
private:
    virtual double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        record.primary_momentum[0] = SampleEnergy(rand, earth_model, cross_sections, record);
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"PrimaryEnergy"};};
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryEnergyDistribution only supports version <= 0!");
        }
    }
};

class PowerLaw : public PrimaryEnergyDistribution {
public:
    double powerLawIndex;
    double energyMin;
    double energyMax;
    double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
        if(energyMin == energyMax)
            return energyMin; //return the only allowed energy

        if(powerLawIndex == 1.0) //sample uniformly in log space
            return pow(10.0, rand->Uniform(log10(energyMin), log10(energyMax)));
        else {
            double u = rand->Uniform();
            double energyP = (1 - u) * pow(energyMin, 1 - powerLawIndex) + u * pow(energyMax, 1 - powerLawIndex);
            return pow(energyP, 1 / (1 - powerLawIndex));
        }
    };
    std::string Name() const override {
        return "PowerLaw";
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new PowerLaw(*this));
    };
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PowerLawIndex", powerLawIndex));
            archive(::cereal::make_nvp("EnergyMin", energyMin));
            archive(::cereal::make_nvp("EnergyMax", energyMax));
            archive(cereal::virtual_base_class<PrimaryEnergyDistribution>(this));
        } else {
            throw std::runtime_error("PowerLaw only supports version <= 0!");
        }
    }
};

class PrimaryDirectionDistribution : public InjectionDistribution {
private:
    virtual earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        earthmodel::Vector3D dir = SampleDirection(rand, earth_model, cross_sections, record);
        double energy = record.primary_momentum[0];
        double mass = record.primary_mass;
        double momentum = std::sqrt(energy*energy - mass*mass);
        record.primary_momentum[1] = momentum * dir.GetX();
        record.primary_momentum[2] = momentum * dir.GetY();
        record.primary_momentum[3] = momentum * dir.GetZ();
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"PrimaryDirection"};};
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<InjectionDistribution>(this));
        } else {
            throw std::runtime_error("PrimaryDirectionDistribution only supports version <= 0!");
        }
    }
};

class IsotropicDirection : public PrimaryDirectionDistribution {
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
        double nx = rand->Uniform(0, 1);
        double ny = rand->Uniform(0, 1);
        double nz = rand->Uniform(0, 1);
        earthmodel::Vector3D res(nx, ny, nz);
        res.normalize();
        return res;
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new IsotropicDirection(*this));
    };
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<PrimaryDirectionDistribution>(this));
        } else {
            throw std::runtime_error("IsotropicDirection only supports version <= 0!");
        }
    }
};

class FixedDirection : public PrimaryDirectionDistribution {
friend cereal::access;
private:
    earthmodel::Vector3D dir;
public:
    FixedDirection(earthmodel::Vector3D dir) : dir(dir) {};
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
        return dir;
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>();};
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new FixedDirection(*this));
    };
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
};

class Cone : public PrimaryDirectionDistribution {
friend cereal::access;
private:
    earthmodel::Vector3D dir;
    earthmodel::Quaternion rotation;
    double opening_angle;
public:
    Cone(earthmodel::Vector3D dir, double opening_angle) : dir(dir), opening_angle(opening_angle) {
        this->dir.normalize();
        if(this->dir == earthmodel::Vector3D(0,0,1)) {
            rotation = earthmodel::Quaternion(0,0,0,1);
        } else {
            earthmodel::Vector3D r = cross_product(earthmodel::Vector3D(0, 0, 1), dir);
            r.normalize();
            rotation = earthmodel::Quaternion(r);
            rotation.SetW(1.0 + dir.GetZ());
        }
    };
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
        double theta = rand->Uniform(0, opening_angle);
        double phi = rand->Uniform(0, 2.0 * M_PI);
        earthmodel::Quaternion q;
        q.SetEulerAnglesZXZr(phi, theta, 0.0);
        return rotation.rotate(q.rotate(earthmodel::Vector3D(0,0,1), false), false);
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new Cone(*this));
    };
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
};

class VertexPositionDistribution : public InjectionDistribution {
private:
    virtual earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
        earthmodel::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, record);
        record.interaction_vertex[0] = pos.GetX();
        record.interaction_vertex[1] = pos.GetY();
        record.interaction_vertex[2] = pos.GetZ();
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"VertexPosition"};};
    virtual std::string Name() const = 0;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
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
};

class CylinderVolumePositionDistribution : public VertexPositionDistribution {
friend cereal::access;
private:
    earthmodel::Cylinder cylinder;
    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        double t = rand->Uniform(0, 2 * M_PI);
        const double outer_radius = cylinder.GetRadius();
        const double inner_radius = cylinder.GetInnerRadius();
        const double height = cylinder.GetZ();
        double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
        double z = rand->Uniform(-height/2.0, height/2.0);
        earthmodel::Vector3D pos(r * cos(t), r * sin(t), z);
        return cylinder.LocalToGlobalPosition(pos);
    };
public:
    CylinderVolumePositionDistribution(earthmodel::Cylinder) : cylinder(cylinder) {};
    std::string Name() const override {
        return "CylinderVolumePositionDistribution";
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new CylinderVolumePositionDistribution(*this));
    };
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
};

class DepthFunction {
public:
    DepthFunction() {};
    virtual double operator()(InteractionSignature const & signature, double energy) const {
        return 0.0;
    };
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
};

class RangeFunction {
public:
    RangeFunction() {};
    virtual double operator()(InteractionSignature const & signature, double energy) const {
        return 0.0;
    };
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
};

class DecayRangeFunction : public RangeFunction {
friend cereal::access;
private:
    double particle_mass; // GeV
    double decay_width; // GeV
    double multiplier;
public:
    DecayRangeFunction(double particle_mass, double decay_width, double multiplier) : particle_mass(particle_mass), decay_width(decay_width), multiplier(multiplier) {};
    double operator()(InteractionSignature const & signature, double energy) const override {
        stga3::FourVector<double> lab_momentum{energy, 0.0, 0.0, sqrt(energy*energy - particle_mass*particle_mass)}; // GeV
        stga3::Beta<double> beta = stga3::beta_to_rest_frame_of(lab_momentum); // dimensionless
        double decay_time = 1.0 / decay_width; // inverse GeV
        stga3::FourVector<double> time_in_rest_frame{decay_time, 0,0,0}; // inverse GeV
        stga3::FourVector<double> time_in_lab_frame = stga3::apply_boost(-beta, time_in_rest_frame); // inverse GeV
        constexpr double iGeV_in_m = 1.973269804593025e-16; // meters per inverse GeV
        double length = (time_in_lab_frame.e0() * beta.norm()) * iGeV_in_m; // meters = ((inverse GeV | dimensionless) * (meters per inverse GeV))
        return length * multiplier; // meters
    };
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
};


class ColumnDepthPositionDistribution : public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    std::shared_ptr<DepthFunction> depth_function;
    std::vector<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
        double t = rand->Uniform(0, 2 * M_PI);
        double r = radius * std::sqrt(rand->Uniform());
        earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
        earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
        return q.rotate(pos, false);
    };

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
        dir.normalize();
        earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

        double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

        earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
        earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

        earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
        path.ExtendFromStartByColumnDepth(lepton_depth, target_types);
        path.ClipToOuterBounds();

        double totalColumnDepth = path.GetColumnDepthInBounds(target_types);

        double traversedColumnDepth = totalColumnDepth * rand->Uniform();
        double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth);
        earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

        return vertex;
    };
public:
    ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::vector<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), depth_function(depth_function), target_types(target_types) {};
    std::string Name() const override {
        return "ColumnDepthPositionDistribution";
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new ColumnDepthPositionDistribution(*this));
    };
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
            std::vector<Particle::ParticleType> t;
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
};

class RangePositionDistribution : public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    std::shared_ptr<RangeFunction> range_function;
    std::vector<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
        double t = rand->Uniform(0, 2 * M_PI);
        double r = radius * std::sqrt(rand->Uniform());
        earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
        earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
        return q.rotate(pos, false);
    };

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
        dir.normalize();
        earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

        double lepton_range = (*range_function)(record.signature, record.primary_momentum[0]);
        record.decay_length = lepton_range;

        earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
        earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

        
        earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
        path.ExtendFromStartByDistance(lepton_range);
        path.ClipToOuterBounds();




        double totalColumnDepth = path.GetColumnDepthInBounds(target_types);

        double traversedColumnDepth = totalColumnDepth * rand->Uniform();
        double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth, target_types);
        earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

        return vertex;
    };
public:
    RangePositionDistribution(){};
    RangePositionDistribution(const RangePositionDistribution &) = default;
    RangePositionDistribution(double radius, double endcap_length, std::shared_ptr<RangeFunction> range_function, std::vector<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {};
    std::string Name() const override {
        return "RangePositionDistribution";
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new RangePositionDistribution(*this));
    };
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
            std::vector<Particle::ParticleType> t;
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
};

class InjectorBase {
protected:
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<LI_random> random;
    Particle::ParticleType primary_type;
    CrossSectionCollection cross_sections;
    std::shared_ptr<earthmodel::EarthModel> earth_model;
    std::vector<std::shared_ptr<InjectionDistribution>> distributions;
public:
    InjectorBase(unsigned int events_to_inject, Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::shared_ptr<earthmodel::EarthModel> earth_model, std::vector<std::shared_ptr<InjectionDistribution>> distributions, std::shared_ptr<LI_random> random) : events_to_inject(events_to_inject), primary_type(primary_type), cross_sections(primary_type, cross_sections), earth_model(earth_model), distributions(distributions), random(random) {};
    InjectorBase(unsigned int events_to_inject, Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<LI_random> random) : events_to_inject(events_to_inject), primary_type(primary_type), cross_sections(primary_type, cross_sections), earth_model(earth_model), random(random) {};
    InjectorBase(unsigned int events_to_inject, CrossSectionCollection cross_sections) : events_to_inject(events_to_inject), cross_sections(cross_sections) {};
    virtual InteractionRecord NewRecord() const {
        InteractionRecord record;
        record.signature.primary_type = primary_type;
        record.primary_mass = Particle(primary_type).GetMass();
        return record;
    };
    void SetRandom(std::shared_ptr<LI_random> random) {
        this->random = random;
    };
    virtual void SampleCrossSection(InteractionRecord & record) const {
        std::vector<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();
        std::vector<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
        std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

        earthmodel::Vector3D intVertex(record.interaction_vertex[0],record.interaction_vertex[1],record.interaction_vertex[2]);

        earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(intVertex, earthmodel::Vector3D());

        double total_prob = 0.0;
        std::vector<Particle::ParticleType> single;
        std::vector<double> probs;
        std::vector<Particle::ParticleType> matching_targets;
        std::vector<InteractionSignature> matching_signatures;
        std::vector<std::shared_ptr<CrossSection>> matching_cross_sections;
        for(auto const target : possible_targets) {
            if(available_targets.find(target) != available_targets.end()) {
                // Get target density
                single.push_back(target);
                double target_density = earth_model->GetDensity(intVertex, single);
                single.clear();
                // Loop over cross sections that have this target
                std::vector<std::shared_ptr<CrossSection>> const & target_cross_sections = cross_sections.GetCrossSectionsForTarget(target);
                for(auto const & cross_section : target_cross_sections) {
                    // Loop over cross section signatures with the same target
                    std::vector<InteractionSignature> signatures = cross_section->GetPossibleSignatures();
                    for(auto const & signature : signatures) {
                        record.signature = signature;
                        // Add total cross section times density to the total prob
                        total_prob += target_density * cross_section->TotalCrossSection(record);
                        // Add total prob to probs
                        probs.push_back(total_prob);
                        // Add target and cross section pointer to the lists
                        matching_targets.push_back(target);
                        matching_cross_sections.push_back(cross_section);
                        matching_signatures.push_back(signature);
                    }
                }
            }
        }
        // Throw a random number
        double r = random->Uniform(0, total_prob);
        // Choose the target and cross section
        unsigned int index = 0;
        for(; (index < probs.size()-1) and (r > probs[index]); ++index) {}
        record.signature.target_type = matching_targets[index];
        record.signature = matching_signatures[index];
        record.target_mass = earth_model->GetTargetMass(record.signature.target_type);
        record.target_momentum = {record.target_mass,0,0,0};
        matching_cross_sections[index]->SampleFinalState(record, random);
    }
    virtual void SampleSecondaryDecay(InteractionRecord & record) const {
        // This function takes an interaction record containing an HNL and simulates the decay to a photon
        // Currently assumes Majorana HNL
        // Final state photon added to secondary particle vectors in InteractionRecord

        // Find the HNL in the secondar particle vector and save its momentum/cartesian direction
        unsigned int lepton_index = (record.signature.secondary_types[0] == Particle::ParticleType::NuF4 or record.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar) ? 0 : 1;
        std::array<double, 4> hnl_momentum = record.secondary_momenta[lepton_index];
        stga3::FourVector<double> pHNL_lab{hnl_momentum[0], hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]};
        double hnl_mass = std::sqrt(pHNL_lab | pHNL_lab);
        earthmodel::Vector3D hnl_dir(hnl_momentum[1],hnl_momentum[2],hnl_momentum[3]);
        hnl_dir.normalize();
        
        // Calculate the decay location of the HNL
        double decay_loc = -1 * record.decay_length * std::log(random->Uniform(0,1));
        record.decay_vertex = {record.interaction_vertex[0] + decay_loc*hnl_dir.GetX(),
                               record.interaction_vertex[1] + decay_loc*hnl_dir.GetY(),
                               record.interaction_vertex[2] + decay_loc*hnl_dir.GetZ()};

        // Sample decay angles
        // Majorana Case: Isotropic Decay
        double costh = random->Uniform(-1,1);
        double theta = std::acos(costh);
        double phi = random->Uniform(0,2*Constants::pi);
        stga3::FourVector<double> pGamma_HNLrest{hnl_mass/2.0,
                                                 hnl_mass/2.0*std::cos(phi)*std::sin(theta),
                                                 hnl_mass/2.0*std::sin(phi)*std::sin(theta),
                                                 hnl_mass/2.0*costh};

        // Boost gamma to lab frame
        stga3::Beta<double> beta_to_hnl_rest = stga3::beta_to_rest_frame_of(pHNL_lab);
        stga3::Boost<double> boost_to_lab = stga3::boost_from_beta(-beta_to_hnl_rest);
        stga3::FourVector<double> pGamma_lab = stga3::apply_boost(boost_to_lab,pGamma_HNLrest);

        std::array<double,4> gamma_momentum;
        gamma_momentum[0] = pGamma_lab.e0();
        gamma_momentum[1] = pGamma_lab.e1();
        gamma_momentum[2] = pGamma_lab.e2();
        gamma_momentum[3] = pGamma_lab.e3();
        record.secondary_momenta.push_back(gamma_momentum);
        record.signature.secondary_types.push_back(Particle::ParticleType::Gamma);
    }
    virtual void SamplePairProduction(InteractionRecord & record) {
        // Nick TODO: finish implementing this function which samples the photon pair produciton location
        earthmodel::Vector3D decay_vtx(record.decay_vertex[0],
                                       record.decay_vertex[1],
                                       record.decay_vertex[2]);
        unsigned int gamma_index = record.secondary_momenta.size() - 1; 
        earthmodel::Vector3D decay_dir(record.secondary_momenta[gamma_index][1],
                                       record.secondary_momenta[gamma_index][2],
                                       record.secondary_momenta[gamma_index][3]);
        decay_dir.normalize();
        earthmodel::Path path(earth_model, decay_vtx, decay_dir, 0);
        path.ComputeIntersections();
        std::vector<double> P;
        double N = 0;
        for
        earthmodel::Geometry::IntersectionList ilist = path.GetIntersections();
        std::cout << "HNL DECAY VERTEX: " << decay_vtx.GetX() << " " << decay_vtx.GetY() << " " << decay_vtx.GetZ() << std::endl;
        std::cout << "HNL DECAY DIRECTION: " << decay_dir.GetX() << " " << decay_dir.GetY() << " " << decay_dir.GetZ() << std::endl;
        std::cout << "INTERSECTIONS:\n";
        for(auto& i : ilist.intersections){
            std::cout << i.hierarchy << " " << i.distance << " " << i.position.GetX() << " " << i.position.GetY() << " " << i.position.GetZ() << std::endl;
        }
        std::cout << "\n";
    
    }
    virtual InteractionRecord GenerateEvent() {
        InteractionRecord record = this->NewRecord();
        for(auto & distribution : distributions) {
            distribution->Sample(random, earth_model, cross_sections, record);
        }
        SampleCrossSection(record);
        injected_events += 1;
        return record;
    };
    virtual std::string Name() const {return("InjectorBase");}
    operator bool() const {
        return injected_events < events_to_inject;
    };
    /*
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<LI_random> random;
    Particle::ParticleType primary_type;
    CrossSectionCollection cross_sections;
    std::shared_ptr<earthmodel::EarthModel> earth_model;
    std::vector<std::shared_ptr<InjectionDistribution>> distributions;
    */
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            archive(::cereal::make_nvp("PrimaryType", primary_type));
            archive(::cereal::make_nvp("CrossSectrions", cross_sections));
            archive(::cereal::make_nvp("EarthModel", earth_model));
            archive(::cereal::make_nvp("InjectionDistributions", distributions));
        } else {
            throw std::runtime_error("InjectorBase only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<InjectorBase> & construct, std::uint32_t const version) {
        if(version == 0) {
            unsigned int events_to_inject;
            unsigned int injected_events;
            Particle::ParticleType primary_type;
            CrossSectionCollection cross_sections;
            std::shared_ptr<earthmodel::EarthModel> earth_model;
            std::vector<std::shared_ptr<InjectionDistribution>> distributions;
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            archive(::cereal::make_nvp("PrimaryType", primary_type));
            archive(::cereal::make_nvp("CrossSectrions", cross_sections));
            archive(::cereal::make_nvp("EarthModel", earth_model));
            archive(::cereal::make_nvp("InjectionDistributions", distributions));
            construct(events_to_inject, cross_sections);
            construct.ptr()->injected_events = injected_events;
            construct.ptr()->primary_type = primary_type;
            construct.ptr()->earth_model = earth_model;
            construct.ptr()->distributions = distributions;
        } else {
            throw std::runtime_error("InjectorBase only supports version <= 0!");
        }
    }
};

class RangedLeptonInjector : public InjectorBase {
    private:
        std::shared_ptr<PrimaryEnergyDistribution> energy_distribution;
        std::shared_ptr<PrimaryDirectionDistribution> direction_distribution;
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution;
        std::shared_ptr<RangeFunction> range_func;
        double disk_radius;
        double endcap_length;
        std::shared_ptr<RangePositionDistribution> position_distribution;
    public:
        RangedLeptonInjector(unsigned int events_to_inject, Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<LI_random> random, std::shared_ptr<PrimaryEnergyDistribution> edist, std::shared_ptr<PrimaryDirectionDistribution> ddist, std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution, std::shared_ptr<RangeFunction> range_func, double disk_radius, double endcap_length) : energy_distribution(edist), direction_distribution(ddist), target_momentum_distribution(target_momentum_distribution), disk_radius(disk_radius), endcap_length(endcap_length), InjectorBase(events_to_inject, primary_type, cross_sections, earth_model, random) {
            std::vector<Particle::ParticleType> target_types = this->cross_sections.TargetTypes();
            position_distribution = std::make_shared<RangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
        };
        virtual InteractionRecord GenerateEvent() override;
        std::string Name() const override {return("RangedInjector");}
};

class VolumeLeptonInjector : public InjectorBase {
    private:
        std::shared_ptr<PrimaryEnergyDistribution> energy_distribution;
        std::shared_ptr<PrimaryDirectionDistribution> direction_distribution;
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution;
        std::shared_ptr<CylinderVolumePositionDistribution> position_distribution;
    public:
        VolumeLeptonInjector(unsigned int events_to_inject, Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<LI_random> random, std::shared_ptr<PrimaryEnergyDistribution> edist, std::shared_ptr<PrimaryDirectionDistribution> ddist, std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution, earthmodel::Cylinder cylinder) : energy_distribution(edist), direction_distribution(ddist), target_momentum_distribution(target_momentum_distribution), position_distribution(std::make_shared<CylinderVolumePositionDistribution>(cylinder)), InjectorBase(events_to_inject, primary_type, cross_sections, earth_model, random) {};
        virtual InteractionRecord GenerateEvent() override;
        std::string Name() const override {return("VolumeInjector");}
};

//----


///Construct a new direction with the given relative angles with respect to
///an existing direction.
///\param base the existing base direction
///\param zenith the angle of the new direction with respect to the base
///\param azimuth the rotation of the new direction about the base
//std::pair<double,double> rotateRelative(std::pair<double,double> base, double zenith, double azimuth);

} //namespace LeptonInjector

CEREAL_CLASS_VERSION(LeptonInjector::InjectionDistribution, 0);

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

CEREAL_CLASS_VERSION(LeptonInjector::InjectorBase, 0);

#endif // LI_LeptonInjector_H

