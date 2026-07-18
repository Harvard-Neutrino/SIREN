#pragma once
#ifndef SIREN_DetectorDirectedScatteringChannel_H
#define SIREN_DetectorDirectedScatteringChannel_H

#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/PhaseSpaceChannel.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

namespace siren {
namespace injection {

// Detector-directed 2 -> 2 scattering channel for target-at-rest models.
//
// A lab direction for one outgoing particle fixes the remaining
// kinematics. Density() reports the full proposal density in Q2 and azimuth,
// or in y and azimuth, where y is energy loss 1 - E_out/E_in or recoil
// (E_out - m_out)/E_in.
//
// Q2 sampling mode (Variable::Q2 only):
//
//   Geometry (default): the outgoing direction is sampled toward the
//     intersection of the target's bounding cone and the accessible
//     kinematic cone. If those cones are disjoint, or the target cone adds no
//     useful restriction, sampling falls back deterministically to uniform CM
//     solid angle. Density() applies the matching analytic CM-to-Q2/y
//     Jacobian in every regime.
//
//   Propagator / Tabulated: Q2 is sampled directly from a 1-D mapping
//     (t-channel propagator 1/(Q2 + m_med^2)^2, or a tabulated CDF of
//     the exact dsigma/dQ2) with a uniform azimuth about the beam axis.
//     Density() returns the mapping's analytic density in Q2 times the uniform
//     1/(2*pi) azimuthal conditional -- no finite-difference Jacobian.
//     Use this to make a biased channel follow the physical Q2 shape
//     (Sample == Density), removing the dominant upscatter weight
//     variance and the numerical-derivative noise.
class DetectorDirectedScatteringChannel : public PhaseSpaceChannel {
public:
    enum class Variable { Q2, BjorkenY, RecoilY };
    enum class Q2Mode { Geometry, Propagator, Tabulated };

    // A positive `volume` supplies the exact target volume for composite or
    // very thin geometries in Volume mode.
    DetectorDirectedScatteringChannel(
        std::shared_ptr<siren::geometry::Geometry const> target,
        int directed_index = 0,
        Variable variable = Variable::Q2,
        DetectorDirected2BodyChannel::Mode mode = DetectorDirected2BodyChannel::Mode::Volume,
        Q2Mode q2_mode = Q2Mode::Geometry,
        double mediator_mass = 0.0,
        std::vector<double> q2_cdf_nodes = {},
        std::vector<double> q2_cdf_values = {},
        double volume = -1.0
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

    bool DirectingActive(
        siren::dataclasses::InteractionRecord const & record) const override;

    void SetVolume(double volume);

private:
    friend cereal::access;

    DetectorDirectedScatteringChannel() = default;

    std::shared_ptr<siren::geometry::Geometry const> target_;
    int directed_index_;
    Variable variable_;
    DetectorDirected2BodyChannel::Mode mode_;
    double target_volume_;

    Q2Mode q2_mode_;
    double mediator_mass_;                  // for Q2Mode::Propagator
    std::shared_ptr<TabulatedMappingTable const> q2_cdf_table_;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (version != 0) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel only supports version <= 0!");
        }

        archive(::cereal::make_nvp("Target", target_));
        archive(::cereal::make_nvp("DirectedIndex", directed_index_));
        archive(::cereal::make_nvp("Variable", static_cast<int>(variable_)));
        archive(::cereal::make_nvp("Mode", static_cast<int>(mode_)));
        archive(::cereal::make_nvp("TargetVolume", target_volume_));
        archive(::cereal::make_nvp("Q2Mode", static_cast<int>(q2_mode_)));
        archive(::cereal::make_nvp("MediatorMass", mediator_mass_));
        archive(::cereal::make_nvp("Q2CdfTable", q2_cdf_table_));
        archive(cereal::virtual_base_class<PhaseSpaceChannel>(this));
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if (version != 0) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel only supports version <= 0!");
        }

        int variable_int;
        int mode_int;
        int q2_mode_int;
        archive(::cereal::make_nvp("Target", target_));
        archive(::cereal::make_nvp("DirectedIndex", directed_index_));
        archive(::cereal::make_nvp("Variable", variable_int));
        archive(::cereal::make_nvp("Mode", mode_int));
        archive(::cereal::make_nvp("TargetVolume", target_volume_));
        archive(::cereal::make_nvp("Q2Mode", q2_mode_int));
        archive(::cereal::make_nvp("MediatorMass", mediator_mass_));
        archive(::cereal::make_nvp("Q2CdfTable", q2_cdf_table_));
        archive(cereal::virtual_base_class<PhaseSpaceChannel>(this));

        if (variable_int != static_cast<int>(Variable::Q2) &&
            variable_int != static_cast<int>(Variable::BjorkenY) &&
            variable_int != static_cast<int>(Variable::RecoilY)) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: invalid Variable value " +
                std::to_string(variable_int) + " in archive");
        }
        if (mode_int != static_cast<int>(DetectorDirected2BodyChannel::Mode::Cone) &&
            mode_int != static_cast<int>(DetectorDirected2BodyChannel::Mode::Volume)) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: invalid Mode value " +
                std::to_string(mode_int) + " in archive");
        }
        if (q2_mode_int != static_cast<int>(Q2Mode::Geometry) &&
            q2_mode_int != static_cast<int>(Q2Mode::Propagator) &&
            q2_mode_int != static_cast<int>(Q2Mode::Tabulated)) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: invalid Q2Mode value " +
                std::to_string(q2_mode_int) + " in archive");
        }

        variable_ = static_cast<Variable>(variable_int);
        mode_ = static_cast<DetectorDirected2BodyChannel::Mode>(mode_int);
        q2_mode_ = static_cast<Q2Mode>(q2_mode_int);
        // Tabulated mode owns exactly one immutable table.
        if ((q2_mode_ == Q2Mode::Tabulated) !=
            static_cast<bool>(q2_cdf_table_)) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: invalid tabulated table state");
        }
        if (q2_mode_ != Q2Mode::Geometry && variable_ != Variable::Q2) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: Q2 mapping requires Variable::Q2");
        }
        if (q2_mode_ == Q2Mode::Propagator && mediator_mass_ <= 0.0) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: invalid propagator mediator mass");
        }
    }

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

CEREAL_CLASS_VERSION(siren::injection::DetectorDirectedScatteringChannel, 0);
CEREAL_REGISTER_TYPE(siren::injection::DetectorDirectedScatteringChannel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
    siren::injection::PhaseSpaceChannel,
    siren::injection::DetectorDirectedScatteringChannel);

#endif // SIREN_DetectorDirectedScatteringChannel_H
