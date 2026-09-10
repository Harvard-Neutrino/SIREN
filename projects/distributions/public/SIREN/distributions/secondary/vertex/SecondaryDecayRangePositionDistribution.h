#pragma once
#ifndef SIREN_SecondaryDecayRangePositionDistribution_H
#define SIREN_SecondaryDecayRangePositionDistribution_H

#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>

#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace geometry { class Geometry; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace distributions {

/**
 * Bias a secondary interaction vertex toward positions from which a proxy
 * daughter is likely to interact or decay inside a fiducial volume.
 *
 * The optimized density is proportional to the physical first-interaction
 * density of the current particle times the conditional probability that the
 * proxy daughter's first process occurs in the fiducial volume. The daughter
 * is approximated as collinear with the current particle and assigned a fixed
 * fraction of its energy. Both legs use detector interaction depth, so layered
 * density profiles, cross sections, and decays are handled consistently.
 * For forward fiducial intervals [a_i, b_i], the continuous target is
 *
 *   q(s) = kappa_1(s) exp(-tau_1(0,s)) P_2(s),
 *   P_2(s) = sum_i exp(-(tau_2(max(s,a_i))-tau_2(s)))
 *                    (1-exp(-(tau_2(b_i)-tau_2(max(s,a_i))))).
 *
 * The parent interaction depth is used as the sampling coordinate. Detector
 * paths provide the forward and inverse depth maps, while the remaining
 * one-dimensional conditional normalization is evaluated adaptively.
 */
class SecondaryDecayRangePositionDistribution
    : virtual public SecondaryVertexPositionDistribution {
    friend cereal::access;

private:
    struct DensityContext;

    std::shared_ptr<siren::geometry::Geometry> fiducial_volume_ = nullptr;
    std::shared_ptr<siren::interactions::InteractionCollection> daughter_interactions_ = nullptr;
    double daughter_mass_ = 0.0;
    double daughter_energy_fraction_ = 1.0;
    double max_length_ = std::numeric_limits<double>::infinity();

    std::shared_ptr<DensityContext> BuildDensityContext(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::InteractionRecord const & record) const;
    double LogDaughterSuccessProbability(
        DensityContext const & context, double distance) const;

public:
    SecondaryDecayRangePositionDistribution();
    SecondaryDecayRangePositionDistribution(
        std::shared_ptr<siren::geometry::Geometry> fiducial_volume,
        std::shared_ptr<siren::interactions::InteractionCollection> daughter_interactions,
        double daughter_mass,
        double daughter_energy_fraction = 1.0,
        double max_length = std::numeric_limits<double>::infinity());
    SecondaryDecayRangePositionDistribution(
        SecondaryDecayRangePositionDistribution const &) = default;

    void SampleVertex(
        std::shared_ptr<siren::utilities::SIREN_random> rand,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::SecondaryDistributionRecord & record) const override;

    double GenerationProbability(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::InteractionRecord const & record) const override;

    std::string Name() const override;
    std::shared_ptr<SecondaryInjectionDistribution> clone() const override;

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
        siren::dataclasses::InteractionRecord const & record) const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version != 0) {
            throw std::runtime_error(
                "SecondaryDecayRangePositionDistribution only supports version <= 0!");
        }
        archive(::cereal::make_nvp("FiducialVolume", fiducial_volume_));
        archive(::cereal::make_nvp("DaughterInteractions", daughter_interactions_));
        archive(::cereal::make_nvp("DaughterMass", daughter_mass_));
        archive(::cereal::make_nvp("DaughterEnergyFraction", daughter_energy_fraction_));
        archive(::cereal::make_nvp("MaxLength", max_length_));
        archive(cereal::virtual_base_class<SecondaryVertexPositionDistribution>(this));
    }

    template<typename Archive>
    static void load_and_construct(
        Archive & archive,
        cereal::construct<SecondaryDecayRangePositionDistribution> & construct,
        std::uint32_t const version) {
        if(version != 0) {
            throw std::runtime_error(
                "SecondaryDecayRangePositionDistribution only supports version <= 0!");
        }
        std::shared_ptr<siren::geometry::Geometry> fiducial_volume;
        std::shared_ptr<siren::interactions::InteractionCollection> daughter_interactions;
        double daughter_mass;
        double daughter_energy_fraction;
        double max_length;
        archive(::cereal::make_nvp("FiducialVolume", fiducial_volume));
        archive(::cereal::make_nvp("DaughterInteractions", daughter_interactions));
        archive(::cereal::make_nvp("DaughterMass", daughter_mass));
        archive(::cereal::make_nvp("DaughterEnergyFraction", daughter_energy_fraction));
        archive(::cereal::make_nvp("MaxLength", max_length));
        construct(fiducial_volume, daughter_interactions, daughter_mass,
                  daughter_energy_fraction, max_length);
        archive(cereal::virtual_base_class<SecondaryVertexPositionDistribution>(construct.ptr()));
    }

protected:
    bool equal(WeightableDistribution const & distribution) const override;
    bool less(WeightableDistribution const & distribution) const override;
};

} // namespace distributions
} // namespace siren

CEREAL_CLASS_VERSION(siren::distributions::SecondaryDecayRangePositionDistribution, 0);
CEREAL_REGISTER_TYPE(siren::distributions::SecondaryDecayRangePositionDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::distributions::SecondaryVertexPositionDistribution,
    siren::distributions::SecondaryDecayRangePositionDistribution);

#endif // SIREN_SecondaryDecayRangePositionDistribution_H
