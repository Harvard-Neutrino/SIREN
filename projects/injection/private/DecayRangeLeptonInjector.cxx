#include <map>
#include <set>
#include <string>
#include <memory>
#include <vector>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/injection/InjectorBase.h"
#include "LeptonInjector/injection/DecayRangeLeptonInjector.h"

namespace LI {
namespace injection {

//---------------
// class DecayRangeLeptonInjector : InjectorBase
//---------------
DecayRangeLeptonInjector::DecayRangeLeptonInjector() {}

DecayRangeLeptonInjector::DecayRangeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::distributions::PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> cross_sections,
        std::shared_ptr<LI::detector::EarthModel> earth_model,
        std::shared_ptr<LI::utilities::LI_random> random,
        std::shared_ptr<LI::distributions::PrimaryEnergyDistribution> edist,
        std::shared_ptr<LI::distributions::PrimaryDirectionDistribution> ddist,
        std::shared_ptr<LI::distributions::TargetMomentumDistribution> target_momentum_distribution,
        std::shared_ptr<LI::distributions::DecayRangeFunction> range_func,
        double disk_radius,
        double endcap_length,
        std::shared_ptr<LI::distributions::PrimaryNeutrinoHelicityDistribution> helicity_distribution) :
    InjectorBase(events_to_inject, primary_injector, cross_sections, earth_model, random),
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    range_func(range_func),
    helicity_distribution(helicity_distribution),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    std::set<LI::dataclasses::Particle::ParticleType> target_types = this->cross_sections->TargetTypes();
    position_distribution = std::make_shared<LI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
    distributions = {target_momentum_distribution, energy_distribution, helicity_distribution, direction_distribution, position_distribution};
}

std::string DecayRangeLeptonInjector::Name() const {
    return("DecayRangeInjector");
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> DecayRangeLeptonInjector::InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

} // namespace injection
} // namespace LI
