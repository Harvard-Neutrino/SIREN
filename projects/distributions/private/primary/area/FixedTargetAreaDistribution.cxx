#include "SIREN/distributions/primary/area/FixedTargetAreaDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt, cos
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <stdlib.h>                                        // for abs

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/detector/DetectorModel.h"            // for DetectorModel
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/geometry/Geometry.h"              // for Geometry
#include "SIREN/math/Vector3D.h"                  // for Vector3D
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }

namespace siren {
namespace distributions {

using detector::DetectorPosition;
using detector::DetectorDirection;

//---------------
// class FixedTargetAreaDistribution : public VertexPositionDistribution
//---------------
siren::math::Vector3D FixedTargetAreaDistribution::SamplePointOfClosestApproach(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    // sample uniformly within the cylinder
    double t = rand->Uniform(0, 2 * M_PI);
    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
    double z = rand->Uniform(-height/2.0, height/2.0);
    siren::math::Vector3D cylinder_pos(r * cos(t), r * sin(t), z);
    siren::math::Vector3D init_pos = cylinder.LocalToGlobalPosition(cylinder_pos);

    siren::math::Vector3D dir = record.GetDirection();
    
    siren::math::Vector3D point_of_closest_approach = init_pos - dir * (dir *init_pos);

    return point_of_closest_approach;
}

double FixedTargetAreaDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {

    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex);

    std::vector<siren::geometry::Geometry::Intersection> intersections = cylinder.Intersections(vertex, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);

    if(intersections.size() == 0) {
        return 0.0;
    }

    double distance = 0.0;
    for(size_t i=0; i*2<intersections.size(); ++i) {
        distance += intersections[2*i+1].distance - intersections[2*i].distance;
    }
    distance = std::abs(distance);

    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double volume = M_PI * ((outer_radius * outer_radius) - (inner_radius * inner_radius)) * height;

    return distance / volume;
}

FixedTargetAreaDistribution::FixedTargetAreaDistribution(siren::geometry::Cylinder cylinder) : cylinder(cylinder) {}

std::string FixedTargetAreaDistribution::Name() const {
    return "FixedTargetAreaDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> FixedTargetAreaDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new FixedTargetAreaDistribution(*this));
}

bool FixedTargetAreaDistribution::equal(WeightableDistribution const & other) const {
    const FixedTargetAreaDistribution* x = dynamic_cast<const FixedTargetAreaDistribution*>(&other);

    if(!x)
        return false;
    else
        return (cylinder == x->cylinder);
}

bool FixedTargetAreaDistribution::less(WeightableDistribution const & other) const {
    const FixedTargetAreaDistribution* x = dynamic_cast<const FixedTargetAreaDistribution*>(&other);
    return cylinder < x->cylinder;
}

bool FixedTargetAreaDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace sirenREN
