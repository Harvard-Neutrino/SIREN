#pragma once
#ifndef SIREN_DetectorDirectedAngularSectorChannel_H
#define SIREN_DetectorDirectedAngularSectorChannel_H

#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

namespace siren {
namespace injection {

// Biased 2-body decay channel that tiles the cone subtended by a target
// geometry into angular sectors and samples one daughter into a single
// (u, phi) sector.
//
//   u    in [0,1] -- fraction of the per-event bounding half-angle
//                    (u=0 is the to-target axis, u=1 the cone edge);
//   phi  in [0,2pi] -- azimuth around the to-target axis.
//
// Within the active regime (the bounding cone fits inside the daughter
// kinematic cone, BoundInKin, or the parent is at rest with the vertex
// outside the target) the sector is a uniform-in-solid-angle proposal
// with lab-angular density 1/Omega_bin, boosted to the lab via the
// standard 2-body Jacobian.  Outside that regime (collimated /
// unreachable / partial-overlap) the sector falls back to isotropic
// rest-frame sampling, exactly like the volume-directed channel where
// directing is inert.
//
// A complete set of sectors (a partition of [0,1] x [0,2pi]) tiles the
// reachable cone disjointly, giving the Kleiss-Pittau optimizer a
// non-degenerate per-sector weight.  Sample and Density share the bin
// test, Omega_bin, axis, and Jacobian, so Sample == Density (Contract C1).
class DetectorDirectedAngularSectorChannel : public PhaseSpaceChannel {
public:
    DetectorDirectedAngularSectorChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        double u_lo,
        double u_hi,
        double phi_lo,
        double phi_hi,
        int daughter_index = 0
    );

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override { return "DetectorDirectedAngularSector"; }
    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Decay2Body;
    }
    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::SolidAngleRest();
    }

    // True when this sector genuinely directs (the active regime), false when it
    // falls back to isotropic 1/4pi.  Lets the optimizer separate genuine
    // directing from the shared fallback.
    bool DirectingActive(
        siren::dataclasses::InteractionRecord const & record) const override;

private:
    friend class cereal::access;

    DetectorDirectedAngularSectorChannel() = default;

    std::shared_ptr<siren::geometry::Geometry const> target_;
    double u_lo_;
    double u_hi_;
    double phi_lo_;
    double phi_hi_;
    int daughter_index_;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Target", target_));
            archive(::cereal::make_nvp("ULower", u_lo_));
            archive(::cereal::make_nvp("UUpper", u_hi_));
            archive(::cereal::make_nvp("PhiLower", phi_lo_));
            archive(::cereal::make_nvp("PhiUpper", phi_hi_));
            archive(::cereal::make_nvp("DaughterIndex", daughter_index_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "DetectorDirectedAngularSectorChannel only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Target", target_));
            archive(::cereal::make_nvp("ULower", u_lo_));
            archive(::cereal::make_nvp("UUpper", u_hi_));
            archive(::cereal::make_nvp("PhiLower", phi_lo_));
            archive(::cereal::make_nvp("PhiUpper", phi_hi_));
            archive(::cereal::make_nvp("DaughterIndex", daughter_index_));
            archive(::cereal::virtual_base_class<PhaseSpaceChannel>(this));
        } else {
            throw std::runtime_error(
                "DetectorDirectedAngularSectorChannel only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(
    siren::injection::DetectorDirectedAngularSectorChannel, 0);
CEREAL_REGISTER_TYPE(
    siren::injection::DetectorDirectedAngularSectorChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::DetectorDirectedAngularSectorChannel);

#endif // SIREN_DetectorDirectedAngularSectorChannel_H
