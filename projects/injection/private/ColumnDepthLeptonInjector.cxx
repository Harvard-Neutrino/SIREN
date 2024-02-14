#include "LeptonInjector/injection/ColumnDepthLeptonInjector.h"

#include <set>
#include <tuple>
#include <string>
#include <vector>

#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/distributions/primary/vertex/ColumnDepthPositionDistribution.h"
#include "LeptonInjector/injection/Injector.h"
#include "LeptonInjector/injection/Process.h"
#include "LeptonInjector/math/Vector3D.h"

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class DepthFunction; } }

namespace LI {
namespace injection {

//---------------
// class ColumnDepthLeptonInjector : Injector
//---------------
ColumnDepthLeptonInjector::ColumnDepthLeptonInjector() {}

ColumnDepthLeptonInjector::ColumnDepthLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<LI::utilities::LI_random> random,
        std::shared_ptr<LI::distributions::DepthFunction> depth_func,
        double disk_radius,
        double endcap_length) :
    Injector(events_to_inject, detector_model, random),
    depth_func(depth_func),
    disk_radius(disk_radius),
    endcap_length(endcap_length)
{
    interactions = primary_process->GetInteractions();
    std::set<LI::dataclasses::Particle::ParticleType> target_types = interactions->TargetTypes();
    position_distribution = std::make_shared<LI::distributions::ColumnDepthPositionDistribution>(disk_radius, endcap_length, depth_func, target_types);
    primary_process->AddPrimaryInjectionDistribution(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      target_types = sec_process->GetInteractions()->TargetTypes();
      sec_process->GetPrimaryInjectionDistributions().push_back(std::make_shared<LI::distributions::ColumnDepthPositionDistribution>(disk_radius, endcap_length, depth_func, target_types));
      */
    }
}

std::string ColumnDepthLeptonInjector::Name() const {
    return("ColumnDepthInjector");
}

std::tuple<LI::math::Vector3D, LI::math::Vector3D> ColumnDepthLeptonInjector::PrimaryInjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(detector_model, interactions, interaction);
}

} // namespace injection
} // namespace LI
