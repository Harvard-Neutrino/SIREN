#include "SIREN/distributions/primary/direction/Cone.h"

#include <array>                                           // for array
#include <math.h>                                          // for acos, cos
#include <string>                                          // for basic_string
#include <stdlib.h>                                        // for abs

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/math/Quaternion.h"                // for Quaternion
#include "SIREN/math/Vector3D.h"                  // for Vector3D
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace distributions {

//---------------
// class Cone : PrimaryDirectionDistribution
//---------------
Cone::Cone(siren::math::Vector3D dir, double opening_angle) : dir(dir), opening_angle(opening_angle) {
    this->dir.normalize();
    if(this->dir == siren::math::Vector3D(0,0,1)) {
        rotation = siren::math::Quaternion(0,0,0,1);
    } else if(this->dir == siren::math::Vector3D(0,0,-1)) {
        rotation = siren::math::Quaternion(0,1,0,0);
    } else {
        siren::math::Vector3D r = cross_product(siren::math::Vector3D(0, 0, 1), dir);
        rotation = siren::math::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
        rotation.normalize();
    }
}

siren::math::Vector3D Cone::SampleDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const{
    double theta = acos(rand->Uniform(cos(opening_angle), 1));
    double phi = rand->Uniform(0, 2.0 * M_PI);
    siren::math::Quaternion q;
    q.SetEulerAnglesZXZr(phi, theta, 0.0);
    return rotation.rotate(q.rotate(siren::math::Vector3D(0,0,1), false), false);
}

double Cone::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    double c = siren::math::scalar_product(dir, event_dir);
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
        return (abs(1.0 - siren::math::scalar_product(dir, x->dir)) < 1e-9
            and opening_angle == x->opening_angle);
}

bool Cone::less(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);
    if(abs(1.0 - siren::math::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        return opening_angle < x->opening_angle;
    }
}

} // namespace distributions
} // namespace siren

