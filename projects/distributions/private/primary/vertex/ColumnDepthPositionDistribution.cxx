#include "SIREN/distributions/primary/vertex/ColumnDepthPositionDistribution.h"

#include <array>
#include <cmath>
#include <tuple>
#include <string>
#include <vector>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/Path.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/primary/vertex/DepthFunction.h"
#include "SIREN/math/Quaternion.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

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
// class ColumnDepthPositionDistribution : VertexPositionDistribution
//---------------
siren::math::Vector3D ColumnDepthPositionDistribution::SampleFromDisk(std::shared_ptr<siren::utilities::SIREN_random> rand, siren::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    siren::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    siren::math::Quaternion q = rotation_between(siren::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> ColumnDepthPositionDistribution::SamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    siren::math::Vector3D dir(record.GetDirection());
    dir.normalize();
    siren::math::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_depth = (*depth_function)(record.type, record.GetEnergy());//note: return is in cgs units!!!

    siren::math::Vector3D endcap_0 = pca - endcap_length * dir;
    siren::math::Vector3D endcap_1 = pca + endcap_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    std::set<siren::dataclasses::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<siren::dataclasses::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    siren::dataclasses::InteractionRecord fake_record;
    fake_record.signature.primary_type = record.type;
    fake_record.primary_mass = record.GetMass();
    fake_record.primary_momentum[0] = record.GetEnergy();
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
    siren::math::Vector3D init_pos = path.GetFirstPoint();
    siren::math::Vector3D vertex = path.GetFirstPoint() + dist * path.GetDirection();

    return {init_pos, vertex};
}

// public getter function for the private SamplePosition function (for debugging)
std::tuple<siren::math::Vector3D, siren::math::Vector3D> ColumnDepthPositionDistribution::GetSamplePosition(std::shared_ptr<siren::utilities::SIREN_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) {

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> samplepos = ColumnDepthPositionDistribution::SamplePosition(rand, detector_model, interactions, record);

    return samplepos;
}

double ColumnDepthPositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex); // m
    siren::math::Vector3D pca = vertex - dir * siren::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_depth = (*depth_function)(record.signature.primary_type, record.primary_momentum[0]);

    siren::math::Vector3D endcap_0 = pca - (endcap_length * dir);
    siren::math::Vector3D endcap_1 = pca + (endcap_length * dir);

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

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
    
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

ColumnDepthPositionDistribution::ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function) : radius(radius), endcap_length(endcap_length), depth_function(depth_function) {}

std::string ColumnDepthPositionDistribution::Name() const {
    return "ColumnDepthPositionDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> ColumnDepthPositionDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new ColumnDepthPositionDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> ColumnDepthPositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex); // m
    siren::math::Vector3D pca = vertex - dir * siren::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));

    double lepton_depth = (*depth_function)(record.signature.primary_type, record.primary_momentum[0]);

    siren::math::Vector3D endcap_0 = pca - endcap_length * dir;
    siren::math::Vector3D endcap_1 = pca + endcap_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool ColumnDepthPositionDistribution::equal(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (depth_function and x->depth_function and *depth_function == *x->depth_function)
                    or (!depth_function and !x->depth_function)
                )
        );
}

bool ColumnDepthPositionDistribution::less(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);
    bool depth_less =
        (!depth_function and x->depth_function) // this->NULL and other->(not NULL)
        or (depth_function and x->depth_function // both not NULL
                and *depth_function < *x->depth_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f)
        <
        std::tie(radius, x->endcap_length, depth_less);
}

} // namespace distributions
} // namespace siren
