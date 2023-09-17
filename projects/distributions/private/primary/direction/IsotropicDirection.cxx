#include "LeptonInjector/distributions/primary/direction/IsotropicDirection.h"

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/direction/PrimaryDirectionDistribution.h"

namespace LI {
namespace distributions {

//---------------
// class IsotropicDirection : PrimaryDirectionDistribution
//---------------
LI::math::Vector3D IsotropicDirection::SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    double nz = rand->Uniform(-1, 1);
    double nr = sqrt(1.0 - nz*nz);
    double phi = rand->Uniform(-M_PI, M_PI);
    double nx = nr * cos(phi);
    double ny = nr * sin(phi);
    LI::math::Vector3D res(nx, ny, nz);
    res.normalize();
    return res;
}

double IsotropicDirection::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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

} // namespace distributions
} // namespace LI
