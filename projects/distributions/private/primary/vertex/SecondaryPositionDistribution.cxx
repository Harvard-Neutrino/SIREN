#include "LeptonInjector/distributions/primary/vertex/SecondaryPositionDistribution.h"

#include <set>                                                    // for set
#include <array>                                                  // for array
#include <cmath>                                                  // for exp
#include <tuple>                                                  // for tie
#include <string>                                                 // for bas...
#include <vector>                                                 // for vector

#include "LeptonInjector/interactions/CrossSection.h"            // for Cro...
#include "LeptonInjector/interactions/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/dataclasses/InteractionSignature.h"      // for Int...
#include "LeptonInjector/dataclasses/Particle.h"                  // for Par...
#include "LeptonInjector/detector/DetectorModel.h"                   // for Ear...
#include "LeptonInjector/detector/Path.h"                         // for Path
#include "LeptonInjector/distributions/Distributions.h"           // for Inj...
#include "LeptonInjector/geometry/Geometry.h"                     // for Geo...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...
#include "LeptonInjector/utilities/Errors.h"                      // for Sec...
#include "LeptonInjector/utilities/Random.h"                      // for LI_...

namespace LI {
namespace distributions {

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
// class SecondaryPositionDistribution : public VertexPositionDistribution
//---------------

LI::math::Vector3D SecondaryPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const {
  throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::SamplePosition without a datum to access the parent"));
  return LI::math::Vector3D(0,0,0);
}

void SecondaryPositionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionTreeDatum & datum) const {
    LI::math::Vector3D pos = SamplePosition(rand, earth_model, interactions, datum);
    datum.record.interaction_vertex[0] = pos.GetX();
    datum.record.interaction_vertex[1] = pos.GetY();
    datum.record.interaction_vertex[2] = pos.GetZ();
}

void SecondaryPositionDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const {
    throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::Sample without a datum to access the parent"));
}

LI::math::Vector3D SecondaryPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionTreeDatum & datum) const {
    LI::math::Vector3D dir(datum.record.primary_momentum[1], datum.record.primary_momentum[2], datum.record.primary_momentum[3]);
    dir.normalize();



    LI::math::Vector3D endcap_0 = LI::math::Vector3D(datum.parent->record.interaction_vertex);
    LI::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), max_length);
    path.ClipToOuterBounds();

    // Check if fiducial volume is provided
    if(fiducial) {
      std::vector<LI::geometry::Geometry::Intersection> fid_intersections = fiducial->Intersections(earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0),
                                                                                                    earth_model->GetEarthCoordDirFromDetCoordDir(dir));
      // If the path intersects the fiducial volume, restrict position to that volume
      if(!fid_intersections.empty()) {
        // make sure the first intersection happens before the maximum generation length
        // and the last intersection happens in front of the generation point
        bool update_path = (fid_intersections.front().distance < max_length
                         && fid_intersections.back().distance > 0);
        if(update_path) {
          LI::math::Vector3D first_point = (fid_intersections.front().distance > 0) ? fid_intersections.front().position : endcap_0;
          LI::math::Vector3D last_point = (fid_intersections.back().distance < max_length) ? fid_intersections.back().position : endcap_1;
          path.SetPoints(first_point,last_point);
        }
      }
    }

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(datum.record);
    LI::dataclasses::InteractionRecord fake_record = datum.record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);
    if(total_interaction_depth == 0) {
        throw(LI::utilities::InjectionFailure("No available interactions along path!"));
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
    LI::math::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double SecondaryPositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
  throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::GenerationProbability without a datum to access the parent"));
  return 0;
}

double SecondaryPositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionTreeDatum const & datum) const {
    LI::math::Vector3D dir(datum.record.primary_momentum[1], datum.record.primary_momentum[2], datum.record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(datum.record.interaction_vertex);

    LI::math::Vector3D endcap_0 = LI::math::Vector3D(datum.parent->record.interaction_vertex);
    LI::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), max_length);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)))
        return 0.0;

    // Check if fiducial volume is provided
    if(fiducial) {
      std::vector<LI::geometry::Geometry::Intersection> fid_intersections = fiducial->Intersections(earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0),
                                                                                                    earth_model->GetEarthCoordDirFromDetCoordDir(dir));
      // If the path intersects the fiducial volume, restrict position to that volume
      if(!fid_intersections.empty()) {
        // make sure the first intersection happens before the maximum generation length
        // and the last intersection happens in front of the generation point
        bool update_path = (fid_intersections.front().distance < max_length
                         && fid_intersections.back().distance > 0);
        if(update_path) {
          LI::math::Vector3D first_point = (fid_intersections.front().distance > 0) ? fid_intersections.front().position : endcap_0;
          LI::math::Vector3D last_point = (fid_intersections.back().distance < max_length) ? fid_intersections.back().position : endcap_1;
          path.SetPoints(first_point,last_point);
        }
      }
    }

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(datum.record);
    LI::dataclasses::InteractionRecord fake_record = datum.record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);

    path.SetPointsWithRay(path.GetFirstPoint(), path.GetDirection(), path.GetDistanceFromStartInBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)));

    double traversed_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);

    double interaction_density = earth_model->GetInteractionDensity(path.GetIntersections(), earth_model->GetEarthCoordPosFromDetCoordPos(vertex), targets, total_cross_sections, total_decay_length);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    return prob_density;
}

SecondaryPositionDistribution::SecondaryPositionDistribution() {}

SecondaryPositionDistribution::SecondaryPositionDistribution(double max_length) : max_length(max_length) {}

SecondaryPositionDistribution::SecondaryPositionDistribution(double max_length, std::shared_ptr<LI::geometry::Geometry> fiducial) :
  max_length(max_length),
  fiducial(fiducial) {}

SecondaryPositionDistribution::SecondaryPositionDistribution(std::shared_ptr<const LI::geometry::Geometry> fiducial) : fiducial(fiducial) {}

std::string SecondaryPositionDistribution::Name() const {
    return "SecondaryPositionDistribution";
}

std::shared_ptr<InjectionDistribution> SecondaryPositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new SecondaryPositionDistribution(*this));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> SecondaryPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
  throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::InjectionBounds without a datum to access the parent"));
  return std::make_pair(LI::math::Vector3D(0,0,0),LI::math::Vector3D(0,0,0));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> SecondaryPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionTreeDatum const & datum) const {
    LI::math::Vector3D dir(datum.record.primary_momentum[1], datum.record.primary_momentum[2], datum.record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(datum.record.interaction_vertex);

    LI::math::Vector3D endcap_0 = LI::math::Vector3D(datum.parent->record.interaction_vertex);
    LI::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), max_length);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));
    return std::pair<LI::math::Vector3D, LI::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool SecondaryPositionDistribution::equal(WeightableDistribution const & other) const {
    const SecondaryPositionDistribution* x = dynamic_cast<const SecondaryPositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (max_length == x->max_length);
}

bool SecondaryPositionDistribution::less(WeightableDistribution const & other) const {
    const SecondaryPositionDistribution* x = dynamic_cast<const SecondaryPositionDistribution*>(&other);
    return
        std::tie(max_length)
        <
        std::tie(x->max_length);
}

} // namespace distributions
} // namespace LeptonInjector
