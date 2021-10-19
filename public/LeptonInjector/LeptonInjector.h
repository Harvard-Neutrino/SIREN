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
    virtual void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const;
    virtual std::vector<std::string> DensityVariables() const;
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::make_shared<InjectionDistribution>(InjectionDistribution(*this));
    };
};

class TargetMomentumDistribution : public InjectionDistribution {
private:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        record.target_momentum = SampleMomentum(rand, earth_model, cross_sections, record);
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"TargetMomentum"};};
};

class TargetAtRest : TargetMomentumDistribution {
private:
public:
    virtual std::array<double, 4> SampleMomentum(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
        return std::array<double, 4>{record.target_mass, 0, 0, 0};
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>();};
};

class PrimaryEnergyDistribution : public InjectionDistribution {
private:
    virtual double SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const override {
        record.primary_momentum[0] = SampleEnergy(rand, earth_model, cross_sections, record);
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"PrimaryEnergy"};};
    virtual std::string Name() const;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
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
};

class FixedDirection : public PrimaryDirectionDistribution {
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
};

class VertexPositionDistribution : public InjectionDistribution {
private:
    virtual earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const = 0;
public:
    void Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
        earthmodel::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, record);
        record.interaction_vertex[0] = pos.GetX();
        record.interaction_vertex[1] = pos.GetY();
        record.interaction_vertex[2] = pos.GetZ();
    };
    virtual std::vector<std::string> DensityVariables() const {return std::vector<std::string>{"VertexPosition"};};
    virtual std::string Name() const;
    virtual std::shared_ptr<InjectionDistribution> clone() const = 0;
};

class CylinderVolumePositionDistribution : public VertexPositionDistribution {
private:
    earthmodel::Cylinder cylinder;
    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
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
};

class RangeFunction {
public:
    RangeFunction();
    virtual double operator()(InteractionSignature const & signature, double energy) const {
        return 0.0;
    };
};

class DecayRangeFunction : public RangeFunction {
private:
    double particle_mass; // GeV
    double decay_width; // GeV
public:
    DecayRangeFunction(double particle_mass, double decay_width) : particle_mass(particle_mass), decay_width(decay_width) {};
    double operator()(InteractionSignature const & signature, double energy) const override {
        stga3::FourVector<double> lab_momentum{energy, energy*energy - particle_mass*particle_mass, 0.0, 0.0}; // GeV
        stga3::Beta<double> beta = stga3::beta_to_rest_frame_of(lab_momentum); // dimensionless
        double decay_time = 1.0 / decay_width; // inverse GeV
        stga3::FourVector<double> time_in_rest_frame{decay_time, 0,0,0}; // inverse GeV
        stga3::FourVector<double> time_in_lab_frame = stga3::apply_boost(-beta, time_in_rest_frame); // inverse GeV
        constexpr double iGeV_in_m = 1.973269804593025e-16; // meters per GeV
        double length = (time_in_lab_frame.e0() * beta.norm()) * iGeV_in_m; // meters = ((inverse GeV | dimensionless) * (meters per GeV))
        return length; // meters
    };
};


class ColumnDepthPositionDistribution : public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    std::function<double(InteractionSignature const &, double)> depth_function;
    std::vector<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
        double t = rand->Uniform(0, 2 * M_PI);
        double r = radius * std::sqrt(rand->Uniform());
        earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
        earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
        return q.rotate(pos, false);
    };

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
        earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
        dir.normalize();
        earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

        double lepton_depth = depth_function(record.signature, record.primary_momentum[0]);

        earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
        earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

        earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
        // TODO method should take vector of targets
        // path.ExtendFromStartByColumnDepth(lepton_depth, target_types);
        path.ClipToOuterBounds();

        // TODO method should take vector of targets
        // double totalColumnDepth = path.GetColumnDepthInBounds(target_types);
        double totalColumnDepth = path.GetColumnDepthInBounds(false);

        double traversedColumnDepth = totalColumnDepth * rand->Uniform();
        double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth);
        earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

        return vertex;
    };
public:
    ColumnDepthPositionDistribution(double radius, double endcap_length, std::function<double(InteractionSignature const &, double)> depth_function, std::vector<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), depth_function(depth_function), target_types(target_types) {};
    std::string Name() const override {
        return "ColumnDepthPositionDistribution";
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new ColumnDepthPositionDistribution(*this));
    };
};

class RangePositionDistribution : public VertexPositionDistribution {
private:
    double radius;
    double endcap_length;
    RangeFunction range_function;
    std::vector<Particle::ParticleType> target_types;

    earthmodel::Vector3D SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
        double t = rand->Uniform(0, 2 * M_PI);
        double r = radius * std::sqrt(rand->Uniform());
        earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
        earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
        return q.rotate(pos, false);
    };

    earthmodel::Vector3D SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const override {
        earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
        dir.normalize();
        earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

        double lepton_range = range_function(record.signature, record.primary_momentum[0]);

        earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
        earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

        earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
        path.ExtendFromStartByDistance(lepton_range);
        path.ClipToOuterBounds();

        // TODO method should take vector of targets
        // double totalColumnDepth = path.GetColumnDepthInBounds(target_types);
        double totalColumnDepth = path.GetColumnDepthInBounds(false);

        double traversedColumnDepth = totalColumnDepth * rand->Uniform();
        double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth);
        //double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth, target_types);
        earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

        return vertex;
    };
public:
    RangePositionDistribution(){};
    RangePositionDistribution(const RangePositionDistribution &) = default;
    RangePositionDistribution(double radius, double endcap_length, RangeFunction range_function, std::vector<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {};
    std::string Name() const override {
        return "RangePositionDistribution";
    };
    virtual std::shared_ptr<InjectionDistribution> clone() const {
        return std::shared_ptr<InjectionDistribution>(new RangePositionDistribution(*this));
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
    virtual void SampleCrossSection(InteractionRecord & record) const {
        std::vector<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();
        std::vector<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
        std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());
        double total_prob = 0.0;
        std::vector<double> probs;
        std::vector<Particle::ParticleType> matching_targets;
        std::vector<InteractionSignature> matching_signatures;
        std::vector<std::shared_ptr<CrossSection>> matching_cross_sections;
        for(auto target : possible_targets) {
            if(available_targets.find(target) != available_targets.end()) {
                // // Get target density
                // double target_density = earth_model->GetTargetDensity(record.interaction_vertex);
                // // Loop over cross sections that have this target
                // std::vector<std::shared_ptr<CrossSection>> const & target_cross_sections = cross_sections.GetCrossSectionsForTarget(target);
                // for(auto const & cross_section : target_cross_sections) {
                //     // Loop over cross section signatures with the same target
                //     std::vector<InteractionSignature> signatures = cross_section->
                //     // Add total cross section times density to the total prob
                //     total_prob += target_density * cross_section->TotalCrossSection(record);
                //     // Add total prob to probs
                //     probs.push_back(total_prob);
                //     // Add target and cross section pointer to the lists
                //     matching_targets.push_back(target);
                //     matching_cross_sections.push_back(cross_section);
                // }
                //
            }
        }
        // Throw a random number
        double r = random->Uniform(0, total_prob);
        // Choose the target and cross section
        unsigned int index = 0;
        for(; (index < probs.size()-1) and (r > probs[index]); ++index) {}
        record.signature.target_type = matching_targets[index];
        matching_cross_sections[index]->SampleFinalState(record, random);
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
};

class RangedLeptonInjector : public InjectorBase {
    private:
        std::shared_ptr<PrimaryEnergyDistribution> energy_distribution;
        std::shared_ptr<PrimaryDirectionDistribution> direction_distribution;
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution;
        RangeFunction range_func;
        double disk_radius;
        double endcap_length;
        std::shared_ptr<RangePositionDistribution> position_distribution;
    public:
        RangedLeptonInjector(unsigned int events_to_inject, Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<LI_random> random, std::shared_ptr<PrimaryEnergyDistribution> edist, std::shared_ptr<PrimaryDirectionDistribution> ddist, std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution, RangeFunction range_func, double disk_radius, double endcap_length) : energy_distribution(edist), direction_distribution(ddist), target_momentum_distribution(target_momentum_distribution), disk_radius(disk_radius), endcap_length(endcap_length), InjectorBase(events_to_inject, primary_type, cross_sections, earth_model, random) {
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

#endif // LI_LeptonInjector_H

