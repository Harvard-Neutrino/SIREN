#include "LeptonInjector/distributions/primary/vertex/OrientedCylinderPositionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt, cos

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/math/Quaternion.h"                // for rotation_b...
#include "LeptonInjector/math/Vector3D.h"                  // for Vector3D
#include "LeptonInjector/utilities/Random.h"               // for LI_random

namespace LI { namespace distributions { class WeightableDistribution; } }

namespace LI {
namespace distributions {

//---------------
// class OrientedCylinderPositionDistribution : VertexPositionDistribution
//---------------
//
LI::math::Vector3D OrientedCylinderPositionDistribution::SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    LI::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    LI::math::Quaternion q = rotation_between(LI::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

LI::math::Vector3D OrientedCylinderPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D pca = SampleFromDisk(rand, dir);

    /*
    std::pair<LI::math::Vector3D, LI::math::Vector3D> GetBounds(earth_model, interactions, pca);

    LI::math::Vector3D p0;
    LI::math::Vector3D p1;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    */
    return LI::math::Vector3D();
}

double OrientedCylinderPositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    return 0.0;
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> OrientedCylinderPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & interaction) const {
    return std::make_pair(LI::math::Vector3D(), LI::math::Vector3D());
}

bool OrientedCylinderPositionDistribution::AreEquivalent(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::DetectorModel const> second_earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> second_interactions) const {
    return false;
}

} // namespace distributions
} // namespace LeptonInjector
