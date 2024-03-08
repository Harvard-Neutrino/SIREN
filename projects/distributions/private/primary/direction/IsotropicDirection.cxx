#include "SIREN/distributions/primary/direction/IsotropicDirection.h"

#include <math.h>                                        // for M_PI, cos, sin
#include <string>                                        // for basic_string

#include "SIREN/distributions/Distributions.h"  // for InjectionDis...
#include "SIREN/math/Vector3D.h"                // for Vector3D
#include "SIREN/utilities/Random.h"             // for LI_random

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }

namespace SI {
namespace distributions {

//---------------
// class IsotropicDirection : PrimaryDirectionDistribution
//---------------
SI::math::Vector3D IsotropicDirection::SampleDirection(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const {
    double nz = rand->Uniform(-1, 1);
    double nr = sqrt(1.0 - nz*nz);
    double phi = rand->Uniform(-M_PI, M_PI);
    double nx = nr * cos(phi);
    double ny = nr * sin(phi);
    SI::math::Vector3D res(nx, ny, nz);
    res.normalize();
    return res;
}

double IsotropicDirection::GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const {
    return 1.0 / (4.0 * M_PI);
}

std::shared_ptr<PrimaryInjectionDistribution> IsotropicDirection::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new IsotropicDirection(*this));
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
} // namespace SI
