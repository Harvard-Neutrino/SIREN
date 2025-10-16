#include "SIREN/distributions/primary/vertex/SphereVolumePositionDistribution.h"

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
#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren { namespace interactions { class InteractionCollection; } }

namespace siren {
namespace distributions {

//---------------
// class SphereVolumePositionDistribution : public VertexPositionDistribution
//---------------
std::tuple<siren::math::Vector3D, siren::math::Vector3D> SphereVolumePositionDistribution::SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {

    // Sample uniformly within a spherical shell
    const double outer_radius = sphere.GetRadius();
    const double inner_radius = sphere.GetInnerRadius();

    // Sample radius with proper volume weighting for uniform distribution
    double outer_vol = outer_radius * outer_radius * outer_radius;
    double inner_vol = inner_radius * inner_radius * inner_radius;
    double r_cubed = rand->Uniform(inner_vol, outer_vol);
    double r = std::cbrt(r_cubed);

    // Sample angles uniformly on sphere surface
    double phi = rand->Uniform(0, 2 * M_PI);  // azimuthal angle
    double cos_theta = rand->Uniform(-1, 1);  // uniform in cos(theta)
    double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

    // Convert to Cartesian coordinates
    siren::math::Vector3D pos(
        r * sin_theta * std::cos(phi),
        r * sin_theta * std::sin(phi),
        r * cos_theta
    );

    siren::math::Vector3D final_pos = sphere.LocalToGlobalPosition(pos);

    siren::math::Vector3D dir = record.GetDirection();
    std::vector<siren::geometry::Geometry::Intersection> intersections = sphere.Intersections(final_pos, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);

    siren::math::Vector3D init_pos;

    if(intersections.size() == 0) {
        init_pos = final_pos;
    } else if(intersections.size() >= 2) {
        init_pos = intersections.front().position;
    } else {
        throw std::runtime_error("Only found one sphere intersection!");
    }

    return {init_pos, final_pos};
}

double SphereVolumePositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D pos(sphere.GlobalToLocalPosition(record.interaction_vertex));
    double r = pos.magnitude();  // Distance from sphere center

    if(r <= sphere.GetInnerRadius() || r >= sphere.GetRadius()) {
        return 0.0;
    } else {
        // Volume of spherical shell: (4/3)π(R³ - r_inner³)
        double outer_volume = (4.0/3.0) * M_PI * sphere.GetRadius() * sphere.GetRadius() * sphere.GetRadius();
        double inner_volume = (4.0/3.0) * M_PI * sphere.GetInnerRadius() * sphere.GetInnerRadius() * sphere.GetInnerRadius();
        return 1.0 / (outer_volume - inner_volume);
    }
}

SphereVolumePositionDistribution::SphereVolumePositionDistribution(siren::geometry::Sphere sphere) : sphere(sphere) {}

std::string SphereVolumePositionDistribution::Name() const {
    return "SphereVolumePositionDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> SphereVolumePositionDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new SphereVolumePositionDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> SphereVolumePositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const {
    siren::math::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D pos(interaction.interaction_vertex);
    std::vector<siren::geometry::Geometry::Intersection> intersections = sphere.Intersections(pos, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);
    if(intersections.size() == 0) {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));
    } else if(intersections.size() >= 2) {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(intersections.front().position, intersections.back().position);
    } else {
        throw std::runtime_error("Only found one sphere intersection!");
    }
}

bool SphereVolumePositionDistribution::equal(WeightableDistribution const & other) const {
    const SphereVolumePositionDistribution* x = dynamic_cast<const SphereVolumePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (sphere == x->sphere);
}

bool SphereVolumePositionDistribution::less(WeightableDistribution const & other) const {
    const SphereVolumePositionDistribution* x = dynamic_cast<const SphereVolumePositionDistribution*>(&other);
    return sphere < x->sphere;
}

bool SphereVolumePositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace siren
