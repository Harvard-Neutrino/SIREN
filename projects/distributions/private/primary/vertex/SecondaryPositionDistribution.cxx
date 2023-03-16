#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/SecondaryPositionDistribution.h"

#include "LeptonInjector/utilities/Errors.h"

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

LI::math::Vector3D SecondaryPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
  throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::SamplePosition without a datum to access the parent"));
  return LI::math::Vector3D(0,0,0);
}

void SecondaryPositionDistribution:: Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum & datum) {
    LI::math::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, datum);
    datum.record.interaction_vertex[0] = pos.GetX();
    datum.record.interaction_vertex[1] = pos.GetY();
    datum.record.interaction_vertex[2] = pos.GetZ();
}

LI::math::Vector3D SecondaryPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum & datum) const {
    LI::math::Vector3D dir(datum.record.primary_momentum[1], datum.record.primary_momentum[2], datum.record.primary_momentum[3]);
    dir.normalize();


    LI::math::Vector3D endcap_0 = LI::math::Vector3D(datum.parent->record.interaction_vertex);
    LI::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), max_length);
    path.ClipToOuterBounds();

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = cross_sections->TotalDecayLength(datum.record);
    LI::dataclasses::InteractionRecord fake_record = datum.record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
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

double SecondaryPositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
  throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::GenerationProbability without a datum to access the parent"));
  return 0;
}

double SecondaryPositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum const & datum) const {
    LI::math::Vector3D dir(datum.record.primary_momentum[1], datum.record.primary_momentum[2], datum.record.primary_momentum[3]);
    dir.normalize();
    
    LI::math::Vector3D endcap_0 = LI::math::Vector3D(datum.parent->record.interaction_vertex);
    LI::math::Vector3D endcap_1 = endcap_0 + max_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), max_length);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)))
        return 0.0;

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LI::dataclasses::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    double total_decay_length = cross_sections->TotalDecayLength(datum.record);
    LI::dataclasses::InteractionRecord fake_record = datum.record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::dataclasses::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
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
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

SecondaryPositionDistribution::SecondaryPositionDistribution() {}

SecondaryPositionDistribution::SecondaryPositionDistribution(double max_length) : max_length(max_length) {}

std::string SecondaryPositionDistribution::Name() const {
    return "SecondaryPositionDistribution";
}

std::shared_ptr<InjectionDistribution> SecondaryPositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new SecondaryPositionDistribution(*this));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> SecondaryPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
  throw(LI::utilities::SecondaryProcessFailure("Cannot call SecondaryPositionDistribution::InjectionBounds without a datum to access the parent"));
  return std::make_pair(LI::math::Vector3D(0,0,0),LI::math::Vector3D(0,0,0));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> SecondaryPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionTreeDatum const & datum) const {
    LI::math::Vector3D dir(datum.record.primary_momentum[1], datum.record.primary_momentum[2], datum.record.primary_momentum[3]);
    dir.normalize();

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
        return (max_length == x->max_length)
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
