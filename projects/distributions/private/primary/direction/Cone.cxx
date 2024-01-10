#include "LeptonInjector/distributions/primary/direction/Cone.h"

#include <array>                                           // for array
#include <math.h>                                          // for acos, cos
#include <string>                                          // for basic_string
#include <stdlib.h>                                        // for abs

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/distributions/Distributions.h"    // for InjectionD...
#include "LeptonInjector/math/Quaternion.h"                // for Quaternion
#include "LeptonInjector/math/Vector3D.h"                  // for Vector3D
#include "LeptonInjector/utilities/Random.h"               // for LI_random

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace detector { class DetectorModel; } }

namespace LI {
namespace distributions {

//---------------
// class Cone : PrimaryDirectionDistribution
//---------------
Cone::Cone(LI::math::Vector3D dir, double opening_angle) : dir(dir), opening_angle(opening_angle) {
    this->dir.normalize();
    if(this->dir == LI::math::Vector3D(0,0,1)) {
        rotation = LI::math::Quaternion(0,0,0,1);
    } else if(this->dir == LI::math::Vector3D(0,0,-1)) {
        rotation = LI::math::Quaternion(0,1,0,0);
    } else {
        LI::math::Vector3D r = cross_product(LI::math::Vector3D(0, 0, 1), dir);
        rotation = LI::math::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
        rotation.normalize();
    }
}

LI::math::Vector3D Cone::SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const{
    double theta = acos(rand->Uniform(cos(opening_angle), 1));
    double phi = rand->Uniform(0, 2.0 * M_PI);
    LI::math::Quaternion q;
    q.SetEulerAnglesZXZr(phi, theta, 0.0);
    return rotation.rotate(q.rotate(LI::math::Vector3D(0,0,1), false), false);
}

double Cone::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    double c = LI::math::scalar_product(dir, event_dir);
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
} // namespace LI

