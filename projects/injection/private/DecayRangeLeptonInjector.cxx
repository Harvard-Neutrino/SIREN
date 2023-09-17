#include "LeptonInjector/injection/DecayRangeLeptonInjector.h"

#include <map>
#include <set>
#include <string>
#include <memory>
#include <vector>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/injection/InjectorBase.h"

namespace LI {
namespace injection {

//---------------
// class DecayRangeLeptonInjector : InjectorBase
//---------------
DecayRangeLeptonInjector::DecayRangeLeptonInjector() {}

DecayRangeLeptonInjector::DecayRangeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::EarthModel> earth_model,
        std::shared_ptr<injection::InjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes,
        std::shared_ptr<LI::utilities::LI_random> random,
        std::shared_ptr<LI::distributions::DecayRangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    InjectorBase(events_to_inject, earth_model, random),
    range_func(range_func),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    cross_sections = primary_process->cross_sections;
    std::set<LI::dataclasses::Particle::ParticleType> target_types = cross_sections->TargetTypes();
    position_distribution = std::make_shared<LI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
    primary_process->injection_distributions.push_back(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      target_types = sec_process->cross_sections->TargetTypes();
      sec_process->injection_distributions.push_back(std::make_shared<LI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types));
      */
    }
}

std::string DecayRangeLeptonInjector::Name() const {
    return("DecayRangeInjector");
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> DecayRangeLeptonInjector::InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

} // namespace injection
} // namespace LI
