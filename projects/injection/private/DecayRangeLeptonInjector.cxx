#include "LeptonInjector/injection/DecayRangeLeptonInjector.h"

#include <set>
#include <string>
#include <vector>

#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangePositionDistribution.h"
#include "LeptonInjector/injection/Injector.h"
#include "LeptonInjector/injection/Process.h"
#include "LeptonInjector/math/Vector3D.h"

namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }

namespace LI {
namespace injection {

//---------------
// class DecayRangeLeptonInjector : Injector
//---------------
DecayRangeLeptonInjector::DecayRangeLeptonInjector() {}

DecayRangeLeptonInjector::DecayRangeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::DetectorModel> earth_model,
        std::shared_ptr<injection::InjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes,
        std::shared_ptr<LI::utilities::LI_random> random,
        std::shared_ptr<LI::distributions::DecayRangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    Injector(events_to_inject, earth_model, random),
    range_func(range_func),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    cross_sections = primary_process->GetInteractions();
    position_distribution = std::make_shared<LI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func);
    primary_process->AddInjectionDistribution(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      target_types = sec_process->GetInteractions()->TargetTypes();
      sec_process->GetInjectionDistributions().push_back(std::make_shared<LI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types));
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
