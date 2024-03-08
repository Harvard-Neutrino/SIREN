#include "SIREN/distributions/primary/vertex/OrientedCylinderPositionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt, cos

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/math/Quaternion.h"                // for rotation_b...
#include "SIREN/math/Vector3D.h"                  // for Vector3D
#include "SIREN/utilities/Random.h"               // for LI_random

namespace SI { namespace distributions { class WeightableDistribution; } }

namespace SI {
namespace distributions {

//---------------
// class OrientedCylinderPositionDistribution : VertexPositionDistribution
//---------------
//
SI::math::Vector3D OrientedCylinderPositionDistribution::SampleFromDisk(std::shared_ptr<SI::utilities::LI_random> rand, SI::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    SI::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    SI::math::Quaternion q = rotation_between(SI::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

std::tuple<SI::math::Vector3D, SI::math::Vector3D> OrientedCylinderPositionDistribution::SamplePosition(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const {
    SI::math::Vector3D dir(record.GetDirection());
    dir.normalize();
    SI::math::Vector3D pca = SampleFromDisk(rand, dir);

    /*
    std::tuple<SI::math::Vector3D, SI::math::Vector3D> GetBounds(detector_model, interactions, pca);

    SI::math::Vector3D p0;
    SI::math::Vector3D p1;

    SI::detector::Path path(detector_model, detector_model->GetEarthCoordPosFromDetCoordPos(endcap_0), detector_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    */
    return {SI::math::Vector3D(), SI::math::Vector3D()};
}

double OrientedCylinderPositionDistribution::GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const {
    return 0.0;
}

std::tuple<SI::math::Vector3D, SI::math::Vector3D> OrientedCylinderPositionDistribution::InjectionBounds(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & interaction) const {
    return std::make_tuple(SI::math::Vector3D(), SI::math::Vector3D());
}

bool OrientedCylinderPositionDistribution::AreEquivalent(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<SI::detector::DetectorModel const> second_detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> second_interactions) const {
    return false;
}

} // namespace distributions
} // namespace SIREN
