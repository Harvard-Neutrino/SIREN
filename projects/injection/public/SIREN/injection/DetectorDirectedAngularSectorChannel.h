#pragma once
#ifndef SIREN_DetectorDirectedAngularSectorChannel_H
#define SIREN_DetectorDirectedAngularSectorChannel_H

#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>

namespace siren { namespace geometry { class Geometry; } }

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
        siren::dataclasses::InteractionRecord const & record) const;

private:
    std::shared_ptr<siren::geometry::Geometry const> target_;
    double u_lo_;
    double u_hi_;
    double phi_lo_;
    double phi_hi_;
    int daughter_index_;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirectedAngularSectorChannel_H
