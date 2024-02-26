#include "LeptonInjector/injection/RangedLeptonInjector.h"

#include <set>
#include <string>
#include <vector>

#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/distributions/primary/vertex/RangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/RangePositionDistribution.h"
#include "LeptonInjector/injection/Injector.h"
#include "LeptonInjector/injection/Process.h"
#include "LeptonInjector/math/Vector3D.h"

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }

namespace LI {
namespace injection {

//---------------
// class RangedLeptonInjector : Injector
//---------------
RangedLeptonInjector::RangedLeptonInjector() {}

RangedLeptonInjector::RangedLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<LI::utilities::LI_random> random,
        std::shared_ptr<LI::distributions::RangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    Injector(events_to_inject, detector_model, random),
    range_func(range_func),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    interactions = primary_process->GetInteractions();
    std::set<LI::dataclasses::Particle::ParticleType> target_types = interactions->TargetTypes();
    position_distribution = std::make_shared<LI::distributions::RangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
    primary_process->AddPrimaryInjectionDistribution(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      target_types = sec_process->GetInteractions()->TargetTypes();
      sec_process->GetPrimaryInjectionDistributions().push_back(std::make_shared<LI::distributions::DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types));
      */
    }
}

std::string RangedLeptonInjector::Name() const {
    return("RangedInjector");
}

std::tuple<LI::math::Vector3D, LI::math::Vector3D> RangedLeptonInjector::PrimaryInjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(detector_model, interactions, interaction);
}

} // namespace injection
} // namespace LI
