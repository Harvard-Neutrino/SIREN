#include "LeptonInjector/distributions/primary/vertex/DecayRangePositionDistribution.h"

#include <array>
#include <cmath>
#include <tuple>
#include <string>

#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/RangeFunction.h"
#include "LeptonInjector/math/Quaternion.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/utilities/Random.h"

namespace LI {
namespace distributions {

//---------------
// class DecayRangePositionDistribution : public VertexPositionDistribution
//---------------
LI::math::Vector3D DecayRangePositionDistribution::SampleFromDisk(std::shared_ptr<LI::utilities::LI_random> rand, LI::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    LI::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    LI::math::Quaternion q = rotation_between(LI::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

LI::math::Vector3D DecayRangePositionDistribution::SamplePosition(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D pca = SampleFromDisk(rand, dir);

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    double y = rand->Uniform();
    double total_distance = path.GetDistance();
    double dist = -decay_length * log(y * (exp(-total_distance/decay_length) - 1) + 1);

    LI::math::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double DecayRangePositionDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex); // m
    LI::math::Vector3D pca = vertex - dir * LI::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    LI::math::Vector3D earth_vertex = earth_model->GetEarthCoordPosFromDetCoordPos(vertex);

    if(not path.IsWithinBounds(earth_vertex))
        return 0.0;

    double total_distance = path.GetDistance();
    double dist = LI::math::scalar_product(path.GetDirection(), earth_vertex - path.GetFirstPoint());

    double prob_density = exp(-dist / decay_length) / (decay_length * (1.0 - exp(-total_distance / decay_length))); // m^-1
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3
    return prob_density;
}

DecayRangePositionDistribution::DecayRangePositionDistribution() {}

DecayRangePositionDistribution::DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function) : radius(radius), endcap_length(endcap_length), range_function(range_function) {}

std::string DecayRangePositionDistribution::Name() const {
    return "DecayRangePositionDistribution";
}

std::shared_ptr<InjectionDistribution> DecayRangePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new DecayRangePositionDistribution(*this));
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> DecayRangePositionDistribution::InjectionBounds(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    LI::math::Vector3D vertex(record.interaction_vertex); // m
    LI::math::Vector3D pca = vertex - dir * LI::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    LI::math::Vector3D endcap_0 = pca - endcap_length * dir;
    LI::math::Vector3D endcap_1 = pca + endcap_length * dir;

    LI::detector::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    LI::math::Vector3D earth_vertex = earth_model->GetEarthCoordPosFromDetCoordPos(vertex);

    if(not path.IsWithinBounds(earth_vertex))
        return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));

    return std::pair<LI::math::Vector3D, LI::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool DecayRangePositionDistribution::equal(WeightableDistribution const & other) const {
    const DecayRangePositionDistribution* x = dynamic_cast<const DecayRangePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (range_function and x->range_function and *range_function == *x->range_function)
                    or (!range_function and !x->range_function)
                )
            );
}

bool DecayRangePositionDistribution::less(WeightableDistribution const & other) const {
    const DecayRangePositionDistribution* x = dynamic_cast<const DecayRangePositionDistribution*>(&other);
    bool range_less =
        (!range_function and x->range_function) // this->NULL and other->(not NULL)
        or (range_function and x->range_function // both not NULL
                and *range_function < *x->range_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f)
        <
        std::tie(radius, x->endcap_length, range_less);
}

bool DecayRangePositionDistribution::AreEquivalent(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<LI::detector::EarthModel const> second_earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> second_cross_sections) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace LeptonInjector
