#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/detector/EarthModelCalculator.h"

#include "LeptonInjector/crosssections/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "LeptonInjector/distributions/primary/direction/FixedDirection.h"

namespace LI {
namespace distributions {

//---------------
// class FixedDirection : PrimaryDirectionDistribution
//---------------
LI::math::Vector3D FixedDirection::SampleDirection(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    return dir;
}

double FixedDirection::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
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

} // namespace distributions
} // namespace LI

