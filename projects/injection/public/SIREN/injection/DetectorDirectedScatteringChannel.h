#pragma once
#ifndef SIREN_DetectorDirectedScatteringChannel_H
#define SIREN_DetectorDirectedScatteringChannel_H

#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

// Detector-directed 2 -> 2 scattering channel for target-at-rest models.
//
// A lab direction for one outgoing particle fixes the remaining
// kinematics. Density() reports the proposal density in the same
// scalar variable used by common SIREN scattering models: Q2, energy
// loss y = 1 - E_out/E_in, or recoil y = (E_out - m_out)/E_in.
class DetectorDirectedScatteringChannel : public PhaseSpaceChannel {
public:
    enum class Variable { Q2, BjorkenY, RecoilY };

    DetectorDirectedScatteringChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int directed_index = 0,
        Variable variable = Variable::Q2,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume
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

    std::string Name() const override { return "DetectorDirectedScattering"; }
    PhaseSpaceTopology Topology() const override {
        return PhaseSpaceTopology::Scatter2to2;
    }
    PhaseSpaceMeasure Measure() const override;

    void SetVolume(double volume);

private:
    std::shared_ptr<siren::geometry::Geometry const> target_;
    int directed_index_;
    Variable variable_;
    DetectorDirected2BodyChannel::Mode mode_;
    double target_volume_;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirectedScatteringChannel_H
