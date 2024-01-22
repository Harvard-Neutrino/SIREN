#include "LeptonInjector/distributions/primary/vertex/RangePositionDistribution.h"

#include <array>
#include <cmath>
#include <tuple>
#include <string>
#include <vector>

#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/detector/DetectorModel.h"
#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/detector/Coordinates.h"
#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/RangeFunction.h"
#include "LeptonInjector/math/Quaternion.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/utilities/Errors.h"
#include "LeptonInjector/utilities/Random.h"

namespace LI {
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
// class RangePositionDistribution : public VertexPositionDistribution
//---------------
LI::math::Vector3D RangePositionDistribution::SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    LI::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    LI::math::Quaternion q = rotation_between(LI::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

LI::math::Vector3D RangePositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(record);
    LI::dataclasses::InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = detector_model->GetTargetMass(target);
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
    LI::math::Vector3D vertex = path.GetFirstPoint() + dist * path.GetDirection();

    return vertex;
}

double RangePositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex); // m
    LI::math::Vector3D pca = vertex - dir * LI::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return 0.0;

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = interactions->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = interactions->TotalDecayLength(record);
    LI::dataclasses::InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = detector_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
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

RangePositionDistribution::RangePositionDistribution() {}

RangePositionDistribution::RangePositionDistribution(double radius, double endcap_length, std::shared_ptr<RangeFunction> range_function, std::set<LI::dataclasses::Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {}

std::string RangePositionDistribution::Name() const {
    return "RangePositionDistribution";
}

std::shared_ptr<InjectionDistribution> RangePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new RangePositionDistribution(*this));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> RangePositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex); // m
    LI::math::Vector3D pca = vertex - dir * LI::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));
    return std::pair<LI::math::Vector3D, LI::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool RangePositionDistribution::equal(WeightableDistribution const & other) const {
    const RangePositionDistribution* x = dynamic_cast<const RangePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (range_function and x->range_function and *range_function == *x->range_function)
                    or (!range_function and !x->range_function)
                )
            and target_types == x->target_types);
}

bool RangePositionDistribution::less(WeightableDistribution const & other) const {
    const RangePositionDistribution* x = dynamic_cast<const RangePositionDistribution*>(&other);
    bool range_less =
        (!range_function and x->range_function) // this->NULL and other->(not NULL)
        or (range_function and x->range_function // both not NULL
                and *range_function < *x->range_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, range_less, x->target_types);
}

} // namespace distributions
} // namespace LeptonInjector
