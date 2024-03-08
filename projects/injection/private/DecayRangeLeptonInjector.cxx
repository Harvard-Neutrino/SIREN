#include "SIREN/injection/DecayRangeSIREN.h"

#include <set>
#include <string>
#include <vector>

#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/distributions/primary/vertex/DecayRangeFunction.h"
#include "SIREN/distributions/primary/vertex/DecayRangePositionDistribution.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }

namespace SI {
namespace injection {

//---------------
// class DecayRangeSIREN : Injector
//---------------
DecayRangeSIREN::DecayRangeSIREN() {}

DecayRangeSIREN::DecayRangeSIREN(
        unsigned int events_to_inject,
        std::shared_ptr<SI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<SI::utilities::LI_random> random,
        std::shared_ptr<SI::distributions::DecayRangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    Injector(events_to_inject, detector_model, random),
    range_func(range_func),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    interactions = primary_process->GetInteractions();
    position_distribution = std::make_shared<SI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func);
    primary_process->AddPrimaryInjectionDistribution(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      target_types = sec_process->GetInteractions()->TargetTypes();
      sec_process->GetPrimaryInjectionDistributions().push_back(std::make_shared<SI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types));
      */
    }
}

std::string DecayRangeSIREN::Name() const {
    return("DecayRangeInjector");
}

std::tuple<SI::math::Vector3D, SI::math::Vector3D> DecayRangeSIREN::PrimaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(detector_model, interactions, interaction);
}

} // namespace injection
} // namespace SI
