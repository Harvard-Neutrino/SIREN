#include "SIREN/distributions/primary/vertex/DecayRangePositionDistribution.h"

#include <array>
#include <cmath>
#include <tuple>
#include <string>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/Path.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/primary/vertex/DecayRangeFunction.h"
#include "SIREN/distributions/primary/vertex/RangeFunction.h"
#include "SIREN/math/Quaternion.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

namespace siren {
namespace distributions {

using detector::DetectorPosition;
using detector::DetectorDirection;

//---------------
// class DecayRangePositionDistribution : public VertexPositionDistribution
//---------------
siren::math::Vector3D DecayRangePositionDistribution::SampleFromDisk(std::shared_ptr<siren::utilities::LI_random> rand, siren::math::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    siren::math::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    siren::math::Quaternion q = rotation_between(siren::math::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> DecayRangePositionDistribution::SamplePosition(std::shared_ptr<siren::utilities::LI_random> rand, std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::PrimaryDistributionRecord & record) const {
    siren::math::Vector3D dir(record.GetDirection());
    dir.normalize();
    siren::math::Vector3D pca = SampleFromDisk(rand, dir);

    double decay_length = range_function->DecayLength(record.type, record.GetEnergy());

    siren::math::Vector3D endcap_0 = pca - endcap_length * dir;
    siren::math::Vector3D endcap_1 = pca + endcap_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    double y = rand->Uniform();
    double total_distance = path.GetDistance();
    double dist = -decay_length * log(y * (exp(-total_distance/decay_length) - 1) + 1);

    siren::math::Vector3D init_pos = path.GetFirstPoint();
    siren::math::Vector3D vertex = path.GetFirstPoint() + dist * path.GetDirection();

    return {init_pos, vertex};
}

double DecayRangePositionDistribution::GenerationProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex); // m
    siren::math::Vector3D pca = vertex - dir * siren::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double decay_length = range_function->DecayLength(record.signature.primary_type, record.primary_momentum[0]);

    siren::math::Vector3D endcap_0 = pca - endcap_length * dir;
    siren::math::Vector3D endcap_1 = pca + endcap_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return 0.0;

    double total_distance = path.GetDistance();
    double dist = siren::math::scalar_product(path.GetDirection(), vertex - path.GetFirstPoint());

    double prob_density = exp(-dist / decay_length) / (decay_length * (1.0 - exp(-total_distance / decay_length))); // m^-1
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3
    return prob_density;
}

DecayRangePositionDistribution::DecayRangePositionDistribution() {}

DecayRangePositionDistribution::DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function) : radius(radius), endcap_length(endcap_length), range_function(range_function) {}

std::string DecayRangePositionDistribution::Name() const {
    return "DecayRangePositionDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> DecayRangePositionDistribution::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new DecayRangePositionDistribution(*this));
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D> DecayRangePositionDistribution::InjectionBounds(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) const {
    siren::math::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    siren::math::Vector3D vertex(record.interaction_vertex); // m
    siren::math::Vector3D pca = vertex - dir * siren::math::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));

    double decay_length = range_function->DecayLength(record.signature.primary_type, record.primary_momentum[0]);

    siren::math::Vector3D endcap_0 = pca - endcap_length * dir;
    siren::math::Vector3D endcap_1 = pca + endcap_length * dir;

    siren::detector::Path path(detector_model, DetectorPosition(endcap_0), DetectorDirection(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(DetectorPosition(vertex)))
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));

    return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
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

bool DecayRangePositionDistribution::AreEquivalent(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<siren::detector::DetectorModel const> second_detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> second_interactions) const {
    return this->operator==(*distribution);
}

} // namespace distributions
} // namespace sirenREN
