#include "SIREN/distributions/primary/direction/FixedDirection.h"

#include <array>                                           // for array
#include <tuple>                                           // for tie, opera...
#include <string>                                          // for basic_string
#include <stdlib.h>                                        // for abs

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/math/Vector3D.h"                  // for Vector3D

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

//---------------
// class FixedDirection : PrimaryDirectionDistribution
//---------------
siren::math::Vector3D FixedDirection::SampleDirection(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    return dir;
}

double FixedDirection::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    if(abs(1.0 - siren::math::scalar_product(dir, event_dir)) < 1e-9)
        return 1.0;
    else
        return 0.0;
}

std::vector<std::string> FixedDirection::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<PrimaryInjectionDistribution> FixedDirection::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new FixedDirection(*this));
}

std::string FixedDirection::Name() const {
    return "FixedDirection";
}

bool FixedDirection::equal(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - siren::math::scalar_product(dir, x->dir)) < 1e-9);
}

bool FixedDirection::less(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);
    if(abs(1.0 - siren::math::scalar_product(dir, x->dir)) < 1e-9) {
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

} // namespace distributions
} // namespace siren

