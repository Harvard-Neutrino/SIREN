#include "SIREN/injection/RangedSIREN.h"

#include <set>
#include <string>
#include <vector>

#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/distributions/primary/vertex/RangeFunction.h"
#include "SIREN/distributions/primary/vertex/RangePositionDistribution.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }

namespace SI {
namespace injection {

//---------------
// class RangedSIREN : Injector
//---------------
RangedSIREN::RangedSIREN() {}

RangedSIREN::RangedSIREN(
        unsigned int events_to_inject,
        std::shared_ptr<SI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<SI::utilities::LI_random> random,
        std::shared_ptr<SI::distributions::RangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    Injector(events_to_inject, detector_model, random),
    range_func(range_func),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    interactions = primary_process->GetInteractions();
    std::set<SI::dataclasses::ParticleType> target_types = interactions->TargetTypes();
    position_distribution = std::make_shared<SI::distributions::RangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
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

std::string RangedSIREN::Name() const {
    return("RangedInjector");
}

std::tuple<SI::math::Vector3D, SI::math::Vector3D> RangedSIREN::PrimaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(detector_model, interactions, interaction);
}

} // namespace injection
} // namespace SI
