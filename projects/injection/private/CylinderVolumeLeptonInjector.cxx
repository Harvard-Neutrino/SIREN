#include <map>
#include <set>
#include <string>
#include <memory>
#include <vector>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/injection/InjectorBase.h"
#include "LeptonInjector/injection/CylinderVolumeLeptonInjector.h"

namespace LI {
namespace injection {

//---------------
// class CylinderVolumeLeptonInjector : InjectorBase
//---------------
CylinderVolumeLeptonInjector::CylinderVolumeLeptonInjector() {}

CylinderVolumeLeptonInjector::CylinderVolumeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::distributions::PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> cross_sections,
        std::shared_ptr<LI::detector::EarthModel> earth_model,
        std::shared_ptr<LI::utilities::LI_random> random,
        std::shared_ptr<LI::distributions::PrimaryEnergyDistribution> edist,
        std::shared_ptr<LI::distributions::PrimaryDirectionDistribution> ddist,
        std::shared_ptr<LI::distributions::TargetMomentumDistribution> target_momentum_distribution,
        LI::geometry::Cylinder cylinder,
        std::shared_ptr<LI::distributions::PrimaryNeutrinoHelicityDistribution> helicity_distribution) :
    InjectorBase(events_to_inject, primary_injector, cross_sections, earth_model, random),
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    position_distribution(std::make_shared<LI::distributions::CylinderVolumePositionDistribution>(cylinder)),
    helicity_distribution(helicity_distribution) {
    distributions = {target_momentum_distribution, energy_distribution, helicity_distribution, direction_distribution, position_distribution};
}

std::string CylinderVolumeLeptonInjector::Name() const {
    return("VolumeInjector");
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> CylinderVolumeLeptonInjector::InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

} // namespace injection
} // namespace LI
