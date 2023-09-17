#include "LeptonInjector/distributions/primary/direction/PrimaryDirectionDistribution.h"

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace distributions {

//---------------
// class PrimaryDirectionDistribution : InjectionDistribution
//---------------
void PrimaryDirectionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D dir = SampleDirection(rand, earth_model, cross_sections, record);
    double energy = record.primary_momentum[0];
    double mass = record.primary_mass;
    double momentum = std::sqrt(energy*energy - mass*mass);
    record.primary_momentum[1] = momentum * dir.GetX();
    record.primary_momentum[2] = momentum * dir.GetY();
    record.primary_momentum[3] = momentum * dir.GetZ();
}

std::vector<std::string> PrimaryDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryDirection"};
}

} // namespace distributions
} // namespace LI

