#include "SIREN/distributions/primary/direction/Cone.h"

#include <array>                                           // for array
#include <math.h>                                          // for acos, cos
#include <string>                                          // for basic_string
#include <stdlib.h>                                        // for abs

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/math/Quaternion.h"                // for Quaternion
#include "SIREN/math/Vector3D.h"                  // for Vector3D
#include "SIREN/utilities/Random.h"               // for LI_random

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace detector { class DetectorModel; } }

namespace SI {
namespace distributions {

//---------------
// class Cone : PrimaryDirectionDistribution
//---------------
Cone::Cone(SI::math::Vector3D dir, double opening_angle) : dir(dir), opening_angle(opening_angle) {
    this->dir.normalize();
    if(this->dir == SI::math::Vector3D(0,0,1)) {
        rotation = SI::math::Quaternion(0,0,0,1);
    } else if(this->dir == SI::math::Vector3D(0,0,-1)) {
        rotation = SI::math::Quaternion(0,1,0,0);
    } else {
        SI::math::Vector3D r = cross_product(SI::math::Vector3D(0, 0, 1), dir);
        rotation = SI::math::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
        rotation.normalize();
    }
}

SI::math::Vector3D Cone::SampleDirection(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const{
    double theta = acos(rand->Uniform(cos(opening_angle), 1));
    double phi = rand->Uniform(0, 2.0 * M_PI);
    SI::math::Quaternion q;
    q.SetEulerAnglesZXZr(phi, theta, 0.0);
    return rotation.rotate(q.rotate(SI::math::Vector3D(0,0,1), false), false);
}

double Cone::GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const {
    SI::math::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    double c = SI::math::scalar_product(dir, event_dir);
    double theta;
    if(c > 1)
        theta = 0;
    else
        theta = acos(c);
    if(theta < opening_angle)
        return 1.0 / (2.0 * M_PI * (1.0 - cos(opening_angle)));
    else
        return 0.0;
}

std::shared_ptr<PrimaryInjectionDistribution> Cone::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new Cone(*this));
}

std::string Cone::Name() const {
    return "Cone";
}

bool Cone::equal(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - SI::math::scalar_product(dir, x->dir)) < 1e-9
            and opening_angle == x->opening_angle);
}

bool Cone::less(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);
    if(abs(1.0 - SI::math::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        return opening_angle < x->opening_angle;
    }
}

} // namespace distributions
} // namespace SI

