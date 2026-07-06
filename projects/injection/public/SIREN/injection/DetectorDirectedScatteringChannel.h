#pragma once
#ifndef SIREN_DetectorDirectedScatteringChannel_H
#define SIREN_DetectorDirectedScatteringChannel_H

#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <memory>
#include <string>
#include <vector>

namespace siren { namespace geometry { class Geometry; } }

namespace siren {
namespace injection {

// Detector-directed 2 -> 2 scattering channel for target-at-rest models.
//
// A lab direction for one outgoing particle fixes the remaining
// kinematics. Density() reports the proposal density in the same
// scalar variable used by common SIREN scattering models: Q2, energy
// loss y = 1 - E_out/E_in, or recoil y = (E_out - m_out)/E_in.
//
// Q2 sampling mode (Variable::Q2 only):
//
//   Geometry (default): the outgoing direction is sampled toward the
//     target geometry and Q2 is read off the resulting polar angle.
//     Density() reports |d cos_theta / d var| via a finite-difference
//     derivative.  Good when geometric pointing dominates acceptance.
//
//   Propagator / Tabulated: Q2 is sampled directly from a 1-D mapping
//     (t-channel propagator 1/(Q2 + m_med^2)^2, or a tabulated CDF of
//     the exact dsigma/dQ2) with a uniform azimuth about the beam axis.
//     Because the 2->2 final state is azimuthally symmetric and the Q2
//     measure marginalizes phi, Density() then returns the mapping's
//     analytic density in Q2 directly -- no finite-difference Jacobian.
//     Use this to make a biased channel follow the physical Q2 shape
//     (Sample == Density), removing the dominant upscatter weight
//     variance and the numerical-derivative noise.
class DetectorDirectedScatteringChannel : public PhaseSpaceChannel {
public:
    enum class Variable { Q2, BjorkenY, RecoilY };
    enum class Q2Mode { Geometry, Propagator, Tabulated };

    DetectorDirectedScatteringChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int directed_index = 0,
        Variable variable = Variable::Q2,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume,
        Q2Mode q2_mode = Q2Mode::Geometry,
        double mediator_mass = 0.0,
        std::vector<double> q2_cdf_nodes = {},
        std::vector<double> q2_cdf_values = {}
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

    Q2Mode q2_mode_;
    double mediator_mass_;                  // for Q2Mode::Propagator
    std::vector<double> q2_cdf_nodes_;      // for Q2Mode::Tabulated
    std::vector<double> q2_cdf_values_;

    // Sample / evaluate Q2 directly from the configured 1-D mapping
    // (Propagator or Tabulated), bypassing geometric direction sampling.
    void SampleFromQ2Mapping(
        std::shared_ptr<siren::utilities::SIREN_random> random,
        siren::dataclasses::InteractionRecord & record) const;

    double DensityFromQ2Mapping(
        siren::dataclasses::InteractionRecord const & record) const;
};

} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirectedScatteringChannel_H
