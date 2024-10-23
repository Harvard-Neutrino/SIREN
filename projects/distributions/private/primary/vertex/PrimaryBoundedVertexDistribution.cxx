#include "SIREN/distributions/primary/vertex/PrimaryBoundedVertexDistribution.h"

#include <set>                                                    // for set
#include <array>                                                  // for array
#include <cmath>                                                  // for exp
#include <tuple>                                                  // for tie
#include <string>                                                 // for bas...
#include <vector>                                                 // for vector

#include "SIREN/interactions/CrossSection.h"            // for Cro...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/dataclasses/InteractionSignature.h"      // for Int...
#include "SIREN/dataclasses/Particle.h"                  // for Par...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...
#include "SIREN/detector/Path.h"                         // for Path
#include "SIREN/detector/Coordinates.h"
#include "SIREN/distributions/Distributions.h"           // for Inj...
#include "SIREN/geometry/Geometry.h"                     // for Geo...
#include "SIREN/math/Vector3D.h"                         // for Vec...
#include "SIREN/utilities/Errors.h"                      // for Sec...
#include "SIREN/utilities/Random.h"                      // for SIREN_...

namespace siren {
namespace distributions {

using detector::DetectorPosition;
using detector::DetectorDirection;
using detector::GeometryPosition;
using detector::GeometryDirection;

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
// class PrimaryBoundedVertexDistribution : public VertexPositionDistribution
//---------------


std::tuple<siren::math::Vector3D, siren::math::Vector3D> PrimaryBoundedVertexDistribution::SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    siren::math::Vector3D pos = record.GetInitialPosition();
    siren::math::Vector3D dir = record.GetDirection();

    siren::math::Vector3D endcap_0 = pos;
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

    siren::dataclasses::InteractionRecord fake_record;
    record.FinalizeAvailable(fake_record);
    std::vector<double> total_cross_sections(targets.size(), 0.0);
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

    return {pos, vertex};
}

double PrimaryBoundedVertexDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex);

    siren::math::Vector3D endcap_0 = record.primary_initial_position;
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

    return prob_density;
}

PrimaryBoundedVertexDistribution::PrimaryBoundedVertexDistribution() {}

PrimaryBoundedVertexDistribution::PrimaryBoundedVertexDistribution(double max_length) : max_length(max_length) {}

PrimaryBoundedVertexDistribution::PrimaryBoundedVertexDistribution(std::shared_ptr<siren::geometry::Geometry> fiducial_volume) : fiducial_volume(fiducial_volume) {}

PrimaryBoundedVertexDistribution::PrimaryBoundedVertexDistribution(std::shared_ptr<siren::geometry::Geometry> fiducial_volume, double max_length) : fiducial_volume(fiducial_volume), max_length(max_length) {}

std::string PrimaryBoundedVertexDistribution::Name() const {
    return "PrimaryBoundedVertexDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> PrimaryBoundedVertexDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new PrimaryBoundedVertexDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> PrimaryBoundedVertexDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex);

    siren::math::Vector3D endcap_0 = record.primary_initial_position;
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
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));
    return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool PrimaryBoundedVertexDistribution::equal(WeightableDistribution const & other) const {
    const PrimaryBoundedVertexDistribution* x = dynamic_cast<const PrimaryBoundedVertexDistribution*>(&other);

    if(!x)
        return false;
    else
        return (max_length == x->max_length);
}

bool PrimaryBoundedVertexDistribution::less(WeightableDistribution const & other) const {
    const PrimaryBoundedVertexDistribution* x = dynamic_cast<const PrimaryBoundedVertexDistribution*>(&other);
    return
        std::tie(max_length)
        <
        std::tie(x->max_length);
}

} // namespace distributions
} // namespace sirenREN
