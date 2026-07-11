#pragma once
#ifndef SIREN_Isotropic2BodyChannel_H
#define SIREN_Isotropic2BodyChannel_H

#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>

namespace siren {
namespace injection {

// Physical (unbiased) 2-body decay channel.
//
// Samples cos(theta*) uniformly in the parent rest frame and phi*
// uniformly in [0, 2*pi].  Boosts the result to the lab frame.
//
// Density: 1/(4*pi) per unit rest-frame solid angle.
// Kinematically closed decays have zero density, and Sample reports a
// retryable InjectionFailure without modifying the secondary momenta.
//
// Supports an angular distribution via the wrapped Decay's
// DifferentialDecayWidth, but defaults to isotropic.
//
// The daughter_index parameter selects which secondary particle
// from the InteractionRecord to treat as daughter A (the other
// is daughter B).
class Isotropic2BodyChannel : public PhaseSpaceChannel {
public:
    // daughter_index: which secondary (0 or 1) in the record is
    // the "daughter A" whose direction we track.
    explicit Isotropic2BodyChannel(int daughter_index = 0);

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord & record
    ) const override;

    double Density(
        std::shared_ptr<siren::detector::DetectorModel const> detector_model,
        siren::dataclasses::InteractionRecord const & record
    ) const override;

    std::string Name() const override { return "Isotropic2Body"; }
    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Decay2Body;
    }
    PhaseSpaceMeasure Measure() const override {
        return PhaseSpaceMeasure::SolidAngleRest();
    }

private:
    int daughter_index_;
};

} // namespace injection
} // namespace siren

#endif // SIREN_Isotropic2BodyChannel_H
