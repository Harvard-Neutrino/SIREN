#ifndef LI_LeptonInjector_H
#define LI_LeptonInjector_H

#include <queue>
#include <memory> // adds shared pointer
#include <iostream>

#include <photospline/bspline.h>
#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"
#include "LeptonInjector/Constants.h"
#include "LeptonInjector/DataWriter.h"
#include "LeptonInjector/EventProps.h"
#include "LeptonInjector/Coordinates.h"
#include "LeptonInjector/BasicInjectionConfiguration.h"

#include "phys-services/CrossSection.h"
#include "earthmodel-service/EarthModel.h"
#include "earthmodel-service/Vector3D.h"

namespace LeptonInjector {

class InjectionDistribution {
private:
public:
    virtual void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord & record) const;
    virtual std::vector<std::string> DensityVariables() const;
};

class PrimaryEnergyDistribution : public InjectionDistribution {
private:
    virtual double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord & record) const override {
        record.primary_momentum[0] = SampleEnergy(rand, earth_model, record);
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"PrimaryEnergy"};};
    virtual std::string Name() const;
};

class PowerLaw : public PrimaryEnergyDistribution {
public:
    double powerLawIndex;
    double energyMin;
    double energyMax;
    double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const override {
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
};

class PrimaryDirectionDistribution : public InjectionDistribution {
private:
    virtual earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord & record) const override {
        earthmodel::Vector3D dir = SampleDirection(rand, earth_model, record);
        double energy = record.primary_momentum[0];
        double mass = record.primary_mass;
        double momentum = std::sqrt(energy*energy - mass*mass);
        record.primary_momentum[1] = momentum * dir.GetX();
        record.primary_momentum[2] = momentum * dir.GetY();
        record.primary_momentum[3] = momentum * dir.GetZ();
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"PrimaryDirection"};};
};

class IsotropicDirection : public PrimaryDirectionDistribution {
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const override {
        double nx = rand->Uniform(0, 1);
        double ny = rand->Uniform(0, 1);
        double nz = rand->Uniform(0, 1);
        earthmodel::Vector3D res(nx, ny, nz);
        res.normalize();
        return res;
    };
};

class FixedDirection : public PrimaryDirectionDistribution {
private:
    earthmodel::Vector3D dir;
public:
    FixedDirection(earthmodel::Vector3D dir) : dir(dir) {};
private:
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const override {
        return dir;
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>();};
};

class Cone : public PrimaryDirectionDistribution {
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
    earthmodel::Vector3D SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const override {
        double theta = rand->Uniform(0, opening_angle);
        double phi = rand->Uniform(0, 2.0 * M_PI);
        earthmodel::Quaternion q;
        q.SetEulerAnglesZXZr(phi, theta, 0.0);
        return rotation.rotate(q.rotate(earthmodel::Vector3D(0,0,1), false), false);
    };
};

class VertexPositionDistribution : public InjectionDistribution {
private:
    virtual earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord & record) const {
        earthmodel::Vector3D pos = SamplePosition(rand, earth_model, record);
        record.interaction_vertex[0] = pos.GetX();
        record.interaction_vertex[1] = pos.GetY();
        record.interaction_vertex[2] = pos.GetZ();
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"VertexPosition"};};
    virtual std::string Name() const;
};

class VolumePositionDistribution : public VertexPositionDistribution {
private:
    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const override {

    };
public:
    VolumePositionDistribution();
    std::string Name() const override {
        return "VolumePositionDistribution";
    };
};

class RangePositionDistribution : public VertexPositionDistribution {
private:
    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, InteractionRecord const & record) const override {

    };
public:
    std::string Name() const override {
        return "RangePositionDistribution";
    };
};

class InjectorBase {
protected:
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<LI_random> random;
    Particle::ParticleType primary_type;
    CrossSectionCollection cross_sections;
    std::shared_ptr<earthmodel::EarthModel> earth_model;
    std::vector<InjectionDistribution> distributions;
public:
    InjectorBase(Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections) : primary_type(primary_type), cross_sections(primary_type, cross_sections) {};
    InjectorBase(CrossSectionCollection cross_sections) : cross_sections(cross_sections) {};
    InteractionRecord NewRecord() const {
        InteractionRecord record;
        record.signature.primary_type = primary_type;
        record.primary_mass = Particle(primary_type).GetMass();
    };
    virtual InteractionRecord GenerateEvent(std::shared_ptr<LI_random>) {
        InteractionRecord record = this->NewRecord();
        for(auto & distribution : distributions) {
            distribution.Sample(random, earth_model, record);
        }
        std::vector<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();
        std::vector<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
        std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());
        double total_prob = 0.0;
        std::vector<double> probs;
        std::vector<Particle::ParticleTypes> matching_targets;
        std::vector<std::shared_ptr<CrossSection>> matching_cross_sections;
        for(auto const & target : targets) {
            if(available_targets.find(target) != available_targets.end()) {
                // Get target density
                // Loop over cross sections that have this target
                // Add total cross section times density to the total prob
                // Add total prob to probs
                // Add target and cross section pointer to the lists
            }
        }
        // Throw a random number
        // Choose the target and cross section
        return record;
    };
};

class RangedLeptonInjector : public InjectorBase {
    public:
        RangedLeptonInjector();
        RangedLeptonInjector(BasicInjectionConfiguration config, std::shared_ptr<earthmodel::EarthModel> earth, std::shared_ptr<LI_random> random_);
        bool Generate() override;
        std::string Name() const override {return("RangedInjector");}
        bool isRanged() const override {return(true);}

        // the earthmodel will just be a null poitner at instantiation
        std::shared_ptr<earthmodel::EarthModel> earthModel;

};

class VolumeLeptonInjector : public InjectorBase {
    public:
        VolumeLeptonInjector();
        VolumeLeptonInjector(BasicInjectionConfiguration config, std::shared_ptr<earthmodel::EarthModel> earth, std::shared_ptr<LI_random> random_);
        bool Generate() override;
        std::string Name() const override {return("VolumeInjector");}
        bool isRanged() const override{return(false);}

        // the earthmodel will just be a null poitner at instantiation
        std::shared_ptr<earthmodel::EarthModel> earthModel;
};

//----


///Construct a new direction with the given relative angles with respect to
///an existing direction.
///\param base the existing base direction
///\param zenith the angle of the new direction with respect to the base
///\param azimuth the rotation of the new direction about the base
//std::pair<double,double> rotateRelative(std::pair<double,double> base, double zenith, double azimuth);

} //namespace LeptonInjector

#endif // LI_LeptonInjector_H

