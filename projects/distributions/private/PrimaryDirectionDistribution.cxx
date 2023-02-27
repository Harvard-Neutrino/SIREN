#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/detector/EarthModelCalculator.h"

#include "LeptonInjector/crosssections/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Particle.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/PrimaryDirectionDistribution.h"

#include "LeptonInjector/utilities/Errors.h"

namespace LI {
namespace distributions {

//---------------
// class PrimaryDirectionDistribution : InjectionDistribution
//---------------
void PrimaryDirectionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const {
    LI::math::Vector3D dir = SampleDirection(rand, earth_model, cross_sections, record);
    double energy = record.primary_momentum[0];
    double mass = record.primary_mass;
    double momentum = std::sqrt(energy*energy - mass*mass);
    record.primary_momentum[1] = momentum * dir.GetX();
    record.primary_momentum[2] = momentum * dir.GetY();
    record.primary_momentum[3] = momentum * dir.GetZ();
}

std::vector<std::string> PrimaryDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryDirection"};
}

//---------------
// class IsotropicDirection : PrimaryDirectionDistribution
//---------------
LI::math::Vector3D IsotropicDirection::SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    double nx = rand->Uniform(0, 1);
    double ny = rand->Uniform(0, 1);
    double nz = rand->Uniform(0, 1);
    LI::math::Vector3D res(nx, ny, nz);
    res.normalize();
    return res;
}

double IsotropicDirection::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    return 1.0 / (4.0 * M_PI);
}

std::shared_ptr<InjectionDistribution> IsotropicDirection::clone() const {
    return std::shared_ptr<InjectionDistribution>(new IsotropicDirection(*this));
}

std::string IsotropicDirection::Name() const {
    return "IsotropicDirection";
}

bool IsotropicDirection::equal(WeightableDistribution const & other) const {
    const IsotropicDirection* x = dynamic_cast<const IsotropicDirection*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool IsotropicDirection::less(WeightableDistribution const & other) const {
    return false;
}

//---------------
// class FixedDirection : PrimaryDirectionDistribution
//---------------
LI::math::Vector3D FixedDirection::SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    return dir;
}

double FixedDirection::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    LI::math::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    if(abs(1.0 - LI::math::scalar_product(dir, event_dir)) < 1e-9)
        return 1.0;
    else
        return 0.0;
}

std::vector<std::string> FixedDirection::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<InjectionDistribution> FixedDirection::clone() const {
    return std::shared_ptr<InjectionDistribution>(new FixedDirection(*this));
}

std::string FixedDirection::Name() const {
    return "FixedDirection";
}

bool FixedDirection::equal(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - LI::math::scalar_product(dir, x->dir)) < 1e-9);
}

bool FixedDirection::less(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);
    if(abs(1.0 - LI::math::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        double X = dir.GetX();
        double Y = dir.GetY();
        double Z = dir.GetZ();
        double other_X = dir.GetX();
        double other_Y = dir.GetY();
        double other_Z = dir.GetZ();
        return
            std::tie(X, Y, Z)
            <
            std::tie(other_X, other_Y, other_Z);
    }
}

//---------------
// class Cone : PrimaryDirectionDistribution
//---------------
Cone::Cone(LI::math::Vector3D dir, double opening_angle) : dir(dir), opening_angle(opening_angle) {
    this->dir.normalize();
    if(this->dir == LI::math::Vector3D(0,0,1)) {
        rotation = LI::math::Quaternion(0,0,0,1);
    } else {
        LI::math::Vector3D r = cross_product(LI::math::Vector3D(0, 0, 1), dir);
        r.normalize();
        rotation = LI::math::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
    }
}

LI::math::Vector3D Cone::SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const{
    double theta = cos(rand->Uniform(acos(opening_angle), 1));
    double phi = rand->Uniform(0, 2.0 * M_PI);
    LI::math::Quaternion q;
    q.SetEulerAnglesZXZr(phi, theta, 0.0);
    return rotation.rotate(q.rotate(LI::math::Vector3D(0,0,1), false), false);
}

double Cone::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    LI::math::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    double theta = acos(LI::math::scalar_product(dir, event_dir));
    if(theta < opening_angle)
        return 1.0 / (2.0 * M_PI * (1.0 - cos(opening_angle)));
    else
        return 0.0;
}

std::shared_ptr<InjectionDistribution> Cone::clone() const {
    return std::shared_ptr<InjectionDistribution>(new Cone(*this));
}

std::string Cone::Name() const {
    return "Cone";
}

bool Cone::equal(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - LI::math::scalar_product(dir, x->dir)) < 1e-9
            and opening_angle == x->opening_angle);
}

bool Cone::less(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);
    if(abs(1.0 - LI::math::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        return opening_angle < x->opening_angle;
    }
}

} // namespace distributions
} // namespace LeptonInjector
