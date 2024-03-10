#include "SIREN/injection/CylinderVolumeSIREN.h"

#include <string>
#include <vector>

#include "SIREN/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace injection {

//---------------
// class CylinderVolumeSIREN : Injector
//---------------
CylinderVolumeSIREN::CylinderVolumeSIREN() {}

CylinderVolumeSIREN::CylinderVolumeSIREN(
        unsigned int events_to_inject,
        std::shared_ptr<siren::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<siren::utilities::SIREN_random> random,
        siren::geometry::Cylinder cylinder) :
    Injector(events_to_inject, detector_model, random),
    position_distribution(std::make_shared<siren::distributions::CylinderVolumePositionDistribution>(cylinder)) {
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

std::tuple<siren::math::Vector3D, siren::math::Vector3D> CylinderVolumeSIREN::PrimaryInjectionBounds(siren::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(detector_model, interactions, interaction);
}

} // namespace injection
} // namespace siren
