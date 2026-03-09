#include "SIREN/distributions/primary/area/PrimaryAreaDistribution.h"

#include <array>                                                  // for array
#include <string>                                                 // for bas...

#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...
#include "SIREN/math/Vector3D.h"                         // for Vec...

namespace siren {
namespace distributions {

//---------------
// class PrimaryAreaDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryAreaDistribution::Sample(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    siren::math::Vector3D point_of_closest_approach = SamplePointOfClosestApproach(rand, detector_model, interactions, record);
    record.SetPointOfClosestApproach((std::array<double, 3>)point_of_closest_approach);
}

std::vector<std::string> PrimaryAreaDistribution::DensityVariables() const {
    return std::vector<std::string>{"PointOfClosestApproach"};
}

bool PrimaryAreaDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution) and detector_model->operator==(*second_detector_model) and interactions->operator==(*second_interactions);
}

} // namespace distributions
} // namespace sirenREN
