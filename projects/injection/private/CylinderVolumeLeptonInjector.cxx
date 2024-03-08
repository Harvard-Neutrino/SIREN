#include "SIREN/injection/CylinderVolumeSIREN.h"

#include <string>
#include <vector>

#include "SIREN/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/math/Vector3D.h"

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }

namespace SI {
namespace injection {

//---------------
// class CylinderVolumeSIREN : Injector
//---------------
CylinderVolumeSIREN::CylinderVolumeSIREN() {}

CylinderVolumeSIREN::CylinderVolumeSIREN(
        unsigned int events_to_inject,
        std::shared_ptr<SI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<SI::utilities::LI_random> random,
        SI::geometry::Cylinder cylinder) :
    Injector(events_to_inject, detector_model, random),
    position_distribution(std::make_shared<SI::distributions::CylinderVolumePositionDistribution>(cylinder)) {
    interactions = primary_process->GetInteractions();
    primary_process->AddPrimaryInjectionDistribution(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      sec_process->GetPrimaryInjectionDistributions().push_back(position_distribution);
      */
    }
}

std::string CylinderVolumeSIREN::Name() const {
    return("VolumeInjector");
}

std::tuple<SI::math::Vector3D, SI::math::Vector3D> CylinderVolumeSIREN::PrimaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(detector_model, interactions, interaction);
}

} // namespace injection
} // namespace SI
