#include "SIREN/distributions/primary/vertex/FixedTargetPositionDistribution.h"

#include <array>                                           // for array
#include <cmath>                                           // for sqrt, cos
#include <string>                                          // for basic_string
#include <vector>                                          // for vector
#include <stdlib.h>                                        // for abs

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/detector/DetectorModel.h"            // for DetectorModel
#include "SIREN/detector/Path.h"                         // for Path
#include "SIREN/distributions/Distributions.h"    // for InjectionD...
#include "SIREN/geometry/Geometry.h"              // for Geometry
#include "SIREN/math/Vector3D.h"                  // for Vector3D
#include "SIREN/utilities/Random.h"               // for SIREN_random
#include "SIREN/utilities/Errors.h"                      // for Sec...

namespace siren { namespace interactions { class InteractionCollection; } }

namespace siren {
namespace distributions {

using detector::DetectorPosition;
using detector::DetectorDirection;

namespace {
double log_one_minus_exp_of_negative(double x) {
    if(x < 1e-1) {
        return std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0;
    } else if(x > 3) {
        double ex = std::exp(-x);
        double ex2 = ex * ex;
        double ex3 = ex2 * ex;
        double ex4 = ex3 * ex;
        double ex5 = ex4 * ex;
        double ex6 = ex5 * ex;
        return -(ex + ex2 / 2.0 + ex3 / 3.0 + ex4 / 4.0 + ex5 / 5.0 + ex6 / 6.0);
    } else {
        return std::log(1.0 - std::exp(-x));
    }
}
}

//---------------
// class FixedTargetPositionDistribution : public VertexPositionDistribution
//---------------
std::tuple<siren::math::Vector3D, siren::math::Vector3D> FixedTargetPositionDistribution::SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    // sample from the cylinder of the target to the get the neutrino production point
    double t = rand->Uniform(0, 2 * M_PI);
    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
    double z = rand->Uniform(-height/2.0, height/2.0);
    siren::math::Vector3D cylinder_pos(r * cos(t), r * sin(t), z);
    siren::math::Vector3D init_pos = cylinder.LocalToGlobalPosition(cylinder_pos);

    siren::math::Vector3D dir = record.GetDirection();

    std::vector<siren::geometry::Geometry::Intersection> intersections = cylinder.Intersections(init_pos, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);

    // now sample from fiducial volume to get interaction vertex
    siren::math::Vector3D endcap_0 = intersections.back().position;
    siren::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), max_length);
    path.ClipToOuterBounds();

    // Check if fiducial volume is provided
    if(fiducial_volume) {
        std::vector<siren::geometry::Geometry::Intersection> fid_intersections = fiducial_volume->Intersections(endcap_0, dir);
        // If the path intersects the fiducial volume, restrict position to that volume
        if(!fid_intersections.empty()) {
            // make sure the first intersection happens before the maximum generation length
            // and the last intersection happens in front of the generation point
            bool update_path = (fid_intersections.front().distance < max_length
                    && fid_intersections.back().distance > 0);
            if(update_path) {
                siren::math::Vector3D first_point((fid_intersections.front().distance > 0) ? fid_intersections.front().position : endcap_0);
                siren::math::Vector3D last_point((fid_intersections.back().distance < max_length) ? fid_intersections.back().position : endcap_1);
                path.SetPoints(DetectorPosition(first_point), DetectorPosition(last_point));
            }
        }
    }

    std::vector<siren::dataclasses::ParticleType> targets(interactions->TargetTypes().begin(), interactions->TargetTypes().end());

    std::vector<double> total_cross_sections(targets.size(), 0.0);
    siren::dataclasses::InteractionRecord fake_record;
    record.FinalizeAvailable(fake_record);
    double total_decay_length = interactions->TotalDecayLength(fake_record);
    for(unsigned int i=0; i<targets.size(); ++i) {
        siren::dataclasses::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = detector_model->GetTargetMass(target);
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSectionAllFinalStates(fake_record);
        }
    }

    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);
    if(total_interaction_depth == 0) {
        throw(siren::utilities::InjectionFailure("No available interactions along path!"));
    }

    double traversed_interaction_depth;
    if(total_interaction_depth < 1e-6) {
        traversed_interaction_depth = rand->Uniform() * total_interaction_depth;
    } else {
        double exp_m_total_interaction_depth = exp(-total_interaction_depth);

        double y = rand->Uniform();
        traversed_interaction_depth = -log(y * exp_m_total_interaction_depth + (1.0 - y));
    }

    double dist = path.GetDistanceFromStartAlongPath(traversed_interaction_depth, targets, total_cross_sections, total_decay_length);
    siren::math::Vector3D vertex = path.GetFirstPoint() + dist * path.GetDirection();

    return {init_pos, vertex};
}

double FixedTargetPositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {

    //pdf = chord length * volume * secondary bounded vertex distriution pdf

    // first copying code from the secondary bounded vertex distribution generation propability function
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex);

    std::vector<siren::geometry::Geometry::Intersection> intersections = cylinder.Intersections(vertex, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);

    // now sample from fiducial volume to get interaction vertex
    siren::math::Vector3D endcap_0 = intersections.back().position;
    siren::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), max_length);
    path.ClipToOuterBounds();

    // Check if fiducial volume is provided
    if(fiducial_volume) {
        std::vector<siren::geometry::Geometry::Intersection> fid_intersections = fiducial_volume->Intersections(endcap_0, dir);
        // If the path intersects the fiducial volume, restrict position to that volume
        if(!fid_intersections.empty()) {
            // make sure the first intersection happens before the maximum generation length
            // and the last intersection happens in front of the generation point
            bool update_path = (fid_intersections.front().distance < max_length
                    && fid_intersections.back().distance > 0);
            if(update_path) {
                DetectorPosition first_point((fid_intersections.front().distance > 0) ? fid_intersections.front().position : endcap_0);
                DetectorPosition last_point((fid_intersections.back().distance < max_length) ? fid_intersections.back().position : endcap_1);
                path.SetPoints(first_point, last_point);
            }
        }
    }

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return 0.0;

    std::set<siren::dataclasses::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<siren::dataclasses::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(record);
    siren::dataclasses::InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        siren::dataclasses::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = detector_model->GetTargetMass(target);
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSectionAllFinalStates(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);

    path.SetPointsWithRay(path.GetFirstPoint(), path.GetDirection(), path.GetDistanceFromStartInBounds(DetectorPosition(vertex)));

    double traversed_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);

    double interaction_density = detector_model->GetInteractionDensity(path.GetIntersections(), DetectorPosition(vertex), targets, total_cross_sections, total_decay_length);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    // now time to calculate the chord length through the target cylinder* volume
    double length = 0.0;
    if(intersections.size() == 0) {
    } else if(intersections.size() >= 2) {
        siren::math::Vector3D delta_intersections = intersections.back().position - intersections.front().position;
        length = delta_intersections.magnitude();
    } else {
        throw std::runtime_error("Only found one cylinder intersection!");
    }

    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetRadius();
    const double height = cylinder.GetZ();
    double target_volume = M_PI * ((outer_radius * outer_radius) - (inner_radius * inner_radius)) * height;

    return prob_density * length / target_volume;
}


FixedTargetPositionDistribution::FixedTargetPositionDistribution(siren::geometry::Cylinder cylinder, std::shared_ptr<siren::geometry::Geometry> fiducial_volume) : cylinder(cylinder), fiducial_volume(fiducial_volume) {}
FixedTargetPositionDistribution::FixedTargetPositionDistribution(siren::geometry::Cylinder cylinder, std::shared_ptr<siren::geometry::Geometry> fiducial_volume, double max_length) : cylinder(cylinder), fiducial_volume(fiducial_volume), max_length(max_length) {}

FixedTargetPositionDistribution::FixedTargetPositionDistribution(siren::geometry::Cylinder cylinder, double max_length) : cylinder(cylinder), max_length(max_length) {}

std::string FixedTargetPositionDistribution::Name() const {
    return "FixedTargetPositionDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> FixedTargetPositionDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new FixedTargetPositionDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> FixedTargetPositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & interaction) const {
    siren::math::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D pos(interaction.interaction_vertex);
    std::vector<siren::geometry::Geometry::Intersection> intersections = cylinder.Intersections(pos, dir);
    siren::detector::DetectorModel::SortIntersections(intersections);

    if(intersections.size() == 0) {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));
    } else if (intersections.size() < 2) {
        throw std::runtime_error("Only found one cylinder intersection!");
    }

    // now sample from fiducial volume to get interaction vertex
    siren::math::Vector3D endcap_0 = intersections.back().position;
    siren::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), max_length);
    path.ClipToOuterBounds();

    // Check if fiducial volume is provided
    if(fiducial_volume) {
        std::vector<siren::geometry::Geometry::Intersection> fid_intersections = fiducial_volume->Intersections(endcap_0, dir);
        // If the path intersects the fiducial volume, restrict position to that volume
        if(!fid_intersections.empty()) {
            // make sure the first intersection happens before the maximum generation length
            // and the last intersection happens in front of the generation point
            bool update_path = (fid_intersections.front().distance < max_length
                    && fid_intersections.back().distance > 0);
            if(update_path) {
                DetectorPosition first_point((fid_intersections.front().distance > 0) ? fid_intersections.front().position : endcap_0);
                DetectorPosition last_point((fid_intersections.back().distance < max_length) ? fid_intersections.back().position : endcap_1);
                path.SetPoints(first_point, last_point);
            }
        }
    }

    if(not path.IsWithinBounds(DetectorPosition(pos))){
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));
    } else {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
    }
}

bool FixedTargetPositionDistribution::equal(WeightableDistribution const & other) const {
    const FixedTargetPositionDistribution* x = dynamic_cast<const FixedTargetPositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (cylinder == x->cylinder);
}

bool FixedTargetPositionDistribution::less(WeightableDistribution const & other) const {
    const FixedTargetPositionDistribution* x = dynamic_cast<const FixedTargetPositionDistribution*>(&other);
    return cylinder < x->cylinder;
}

bool FixedTargetPositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace sirenREN
