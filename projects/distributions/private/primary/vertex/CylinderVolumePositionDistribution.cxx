#include "SIREN/distributions/primary/vertex/CylinderVolumePositionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt, cos
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <stdlib.h>                                        // for abs

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/detector/DetectorModel.h"            // for DetectorModel
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/geometry/Geometry.h"              // for Geometry
#include "SIREN/math/Vector3D.h"                  // for Vector3D
#include "SIREN/utilities/Random.h"               // for LI_random

namespace siren { namespace interactions { class InteractionCollection; } }

namespace siren {
namespace distributions {

//---------------
// class CylinderVolumePositionDistribution : public VertexPositionDistribution
//---------------
std::tuple<siren::math::Vector3D, siren::math::Vector3D> CylinderVolumePositionDistribution::SamplePosition(std::shared_ptr<siren::utilities::LI_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    double t = rand->Uniform(0, 2 * M_PI);
    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
    double z = rand->Uniform(-height/2.0, height/2.0);
    siren::math::Vector3D pos(r * cos(t), r * sin(t), z);
    siren::math::Vector3D final_pos = cylinder.LocalToGlobalPosition(pos);

    siren::math::Vector3D dir = record.GetDirection();
    std::vector<siren::geometry::Geometry::Intersection> intersections = cylinder.Intersections(final_pos, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);

    siren::math::Vector3D init_pos;

    if(intersections.size() == 0) {
        init_pos = final_pos;
    } else if(intersections.size() >= 2) {
        init_pos = intersections.front().position;
    } else {
        throw std::runtime_error("Only found one cylinder intersection!");
    }

    return {init_pos, final_pos};
}

double CylinderVolumePositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D pos(cylinder.GlobalToLocalPosition(record.interaction_vertex));
    double z = pos.GetZ();
    double r = sqrt(pos.GetX() * pos.GetX() + pos.GetY() * pos.GetY());
    if(abs(z) >= 0.5 * cylinder.GetZ()
            or r <= cylinder.GetInnerRadius()
            or r >= cylinder.GetRadius()) {
        return 0.0;
    } else {
        return 1.0 / (M_PI*(cylinder.GetRadius() * cylinder.GetRadius() - cylinder.GetInnerRadius() * cylinder.GetInnerRadius()) * cylinder.GetZ());
    }
}


CylinderVolumePositionDistribution::CylinderVolumePositionDistribution(siren::geometry::Cylinder cylinder) : cylinder(cylinder) {}

std::string CylinderVolumePositionDistribution::Name() const {
    return "CylinderVolumePositionDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> CylinderVolumePositionDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new CylinderVolumePositionDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> CylinderVolumePositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const {
    siren::math::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D pos(interaction.interaction_vertex);
    std::vector<siren::geometry::Geometry::Intersection> intersections = cylinder.Intersections(pos, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);
    if(intersections.size() == 0) {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));
    } else if(intersections.size() >= 2) {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(intersections.front().position, intersections.back().position);
    } else {
        throw std::runtime_error("Only found one cylinder intersection!");
    }
}

bool CylinderVolumePositionDistribution::equal(WeightableDistribution const & other) const {
    const CylinderVolumePositionDistribution* x = dynamic_cast<const CylinderVolumePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (cylinder == x->cylinder);
}

bool CylinderVolumePositionDistribution::less(WeightableDistribution const & other) const {
    const CylinderVolumePositionDistribution* x = dynamic_cast<const CylinderVolumePositionDistribution*>(&other);
    return cylinder < x->cylinder;
}

bool CylinderVolumePositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace sirenREN
