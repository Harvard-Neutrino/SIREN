#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/detector/EarthModelCalculator.h"

#include "LeptonInjector/crosssections/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Particle.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/ColumnDepthPositionDistribution.h"

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
// class DepthFunction
//---------------
DepthFunction::DepthFunction() {}

double DepthFunction::operator()(LI::crosssections::InteractionSignature const & signature, double energy) const {
    return 0.0;
}

bool DepthFunction::operator==(DepthFunction const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool DepthFunction::operator<(DepthFunction const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

//---------------
// class LeptonDepthFunction : DepthFunction
//---------------
LeptonDepthFunction::LeptonDepthFunction() {}

void LeptonDepthFunction::SetMuParams(double mu_alpha, double mu_beta) {
    this->mu_alpha = mu_alpha;
    this->mu_beta = mu_beta;
}

void LeptonDepthFunction::SetTauParams(double tau_alpha, double tau_beta) {
    this->tau_alpha = tau_alpha;
    this->tau_beta = tau_beta;

}

void LeptonDepthFunction::SetScale(double scale) {
    this->scale = scale;
}

void LeptonDepthFunction::SetMaxDepth(double max_depth) {
    this->max_depth = max_depth;
}

double LeptonDepthFunction::GetMuAlpha() const {
    return mu_alpha;
}

double LeptonDepthFunction::GetMuBeta() const {
    return mu_beta;
}

double LeptonDepthFunction::GetTauAlpha() const {
    return tau_alpha;
}

double LeptonDepthFunction::GetTauBeta() const {
    return tau_beta;
}

double LeptonDepthFunction::GetScale() const {
    return scale;
}

double LeptonDepthFunction::GetMaxDepth() const {
    return max_depth;
}

double LeptonDepthFunction::operator()(LI::crosssections::InteractionSignature const & signature, double energy) const {
    double range = log(1.0 + energy * mu_beta / mu_alpha) / mu_beta;
    if(tau_primaries.count(signature.primary_type) > 0)
        range += log(1.0 + energy * tau_beta / tau_alpha) / tau_beta;
    return std::min(range, max_depth);
}

bool LeptonDepthFunction::equal(DepthFunction const & other) const {
    const LeptonDepthFunction* x = dynamic_cast<const LeptonDepthFunction*>(&other);

    if(not x)
        return false;

    return
        std::tie(mu_alpha, mu_beta, tau_alpha, tau_beta, scale, max_depth, tau_primaries)
        ==
        std::tie(x->mu_alpha, x->mu_beta, x->tau_alpha, x->tau_beta, x->scale, x->max_depth, x->tau_primaries);
}

bool LeptonDepthFunction::less(DepthFunction const & other) const {
    const LeptonDepthFunction* x = dynamic_cast<const LeptonDepthFunction*>(&other);

    if(not x)
        return false;

    return
        std::tie(mu_alpha, mu_beta, tau_alpha, tau_beta, scale, max_depth, tau_primaries)
        <
        std::tie(x->mu_alpha, x->mu_beta, x->tau_alpha, x->tau_beta, x->scale, x->max_depth, x->tau_primaries);
}

//---------------
// class ColumnDepthPositionDistribution : VertexPositionDistribution
//---------------
LI::math::Vector3D ColumnDepthPositionDistribution::SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    LI::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    LI::math::Quaternion q = rotation_between(LI::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

LI::math::Vector3D ColumnDepthPositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    std::set<LI::utilities::Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LI::utilities::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    LI::crosssections::InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::utilities::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);
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

    double dist = path.GetDistanceFromStartAlongPath(traversed_interaction_depth, targets, total_cross_sections);
    LI::math::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double ColumnDepthPositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex); // m
    LI::math::Vector3D pca = vertex - dir * LI::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    std::set<LI::utilities::Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LI::utilities::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    LI::crosssections::InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LI::utilities::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    path.SetPointsWithRay(path.GetFirstPoint(), path.GetDirection(), path.GetDistanceFromStartInBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)));

    double traversed_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    double interaction_density = earth_model->GetInteractionDensity(path.GetIntersections(), earth_model->GetEarthCoordPosFromDetCoordPos(vertex), targets, total_cross_sections);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

ColumnDepthPositionDistribution::ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::set<LI::utilities::Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), depth_function(depth_function), target_types(target_types) {}

std::string ColumnDepthPositionDistribution::Name() const {
    return "ColumnDepthPositionDistribution";
}

std::shared_ptr<InjectionDistribution> ColumnDepthPositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ColumnDepthPositionDistribution(*this));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> ColumnDepthPositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::crosssections::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex); // m
    LI::math::Vector3D pca = vertex - dir * LI::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    return std::pair<LI::math::Vector3D, LI::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
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
            and target_types == x->target_types);
}

bool ColumnDepthPositionDistribution::less(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);
    bool depth_less =
        (!depth_function and x->depth_function) // this->NULL and other->(not NULL)
        or (depth_function and x->depth_function // both not NULL
                and *depth_function < *x->depth_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, depth_less, x->target_types);
}

} // namespace distributions
} // namespace LeptonInjector
