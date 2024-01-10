#include "LeptonInjector/distributions/primary/vertex/CylinderVolumePositionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt, cos
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <stdlib.h>                                        // for abs

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/detector/EarthModel.h"            // for EarthModel
#include "LeptonInjector/distributions/Distributions.h"    // for InjectionD...
#include "LeptonInjector/geometry/Geometry.h"              // for Geometry
#include "LeptonInjector/math/Vector3D.h"                  // for Vector3D
#include "LeptonInjector/utilities/Random.h"               // for LI_random

namespace LI { namespace crosssections { class InteractionCollection; } }

namespace LI {
namespace distributions {

//---------------
// class CylinderVolumePositionDistribution : public VertexPositionDistribution
//---------------
LI::math::Vector3D CylinderVolumePositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    double t = rand->Uniform(0, 2 * M_PI);
    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
    double z = rand->Uniform(-height/2.0, height/2.0);
    LI::math::Vector3D pos(r * cos(t), r * sin(t), z);
    return cylinder.LocalToGlobalPosition(pos);
}

double CylinderVolumePositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D pos(record.interaction_vertex);
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


CylinderVolumePositionDistribution::CylinderVolumePositionDistribution(LI::geometry::Cylinder cylinder) : cylinder(cylinder) {}

std::string CylinderVolumePositionDistribution::Name() const {
    return "CylinderVolumePositionDistribution";
}

std::shared_ptr<InjectionDistribution> CylinderVolumePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new CylinderVolumePositionDistribution(*this));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> CylinderVolumePositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & interaction) const {
    LI::math::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D pos(interaction.interaction_vertex);
    std::vector<LI::geometry::Geometry::Intersection> intersections = cylinder.Intersections(pos, dir);
    LI::detector::EarthModel::SortIntersections(intersections);
    if(intersections.size() == 0) {
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));
    } else if(intersections.size() >= 2) {
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(intersections.front().position, intersections.back().position);
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

bool CylinderVolumePositionDistribution::AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::crosssections::InteractionCollection const> second_cross_sections) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace LeptonInjector
