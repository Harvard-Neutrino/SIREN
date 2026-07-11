#include "SIREN/injection/DetectorDirectedScatteringChannel.h"

#include "DetectorDirectedChannelUtils.h"
#include "InteractionRecordUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace siren {
namespace injection {

namespace {

using detail::FourVector;
using detail::ReadPrimary;
using detail::ReadSecondary;
using detail::WriteSecondary;

double MomentumMagnitude(double e, double m) {
    double p2 = e * e - m * m;
    if (p2 <= 0.0) return 0.0;
    return std::sqrt(p2);
}

double ComputeQ2(
    FourVector const & p1,
    FourVector const & p3,
    double m1,
    double m3)
{
    return 2.0 * (p1.e * p3.e - p1.p * p3.p) - m1 * m1 - m3 * m3;
}

double CosThetaFromQ2(
    double Q2,
    double E1,
    double p1_mag,
    double m1,
    double m2,
    double m3,
    double m4)
{
    if (m2 <= 0.0 || p1_mag <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double E4 = (Q2 + m2 * m2 + m4 * m4) / (2.0 * m2);
    double E3 = E1 + m2 - E4;
    double p3_mag = MomentumMagnitude(E3, m3);
    if (p3_mag <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return (E1 * E3 - 0.5 * (Q2 + m1 * m1 + m3 * m3)) / (p1_mag * p3_mag);
}

// Kinematic Q^2 limits for the 2->2 process 1 + 2(rest) -> 3 + 4, from
// the CM scattering angle endpoints (cos_cm = +1 forward -> Q2min,
// cos_cm = -1 backward -> Q2max).  Q^2 = -(p1 - p3)^2.
bool Q2Range(double E1, double m1, double m2, double m3, double m4,
             double & q2min, double & q2max) {
    double s = m1 * m1 + m2 * m2 + 2.0 * m2 * E1;
    if (s <= 0.0) return false;
    double root_s = std::sqrt(s);
    if (root_s <= m3 + m4) return false;
    double E1cm = (s + m1 * m1 - m2 * m2) / (2.0 * root_s);
    double E3cm = (s + m3 * m3 - m4 * m4) / (2.0 * root_s);
    double p1cm = MomentumMagnitude(E1cm, m1);
    double p3cm = MomentumMagnitude(E3cm, m3);
    if (p1cm <= 0.0 || p3cm <= 0.0) return false;
    // Q^2 = -(m1^2 + m3^2 - 2 E1cm E3cm + 2 p1cm p3cm cos_cm)
    double base = m1 * m1 + m3 * m3 - 2.0 * E1cm * E3cm;
    q2min = -(base + 2.0 * p1cm * p3cm);   // cos_cm = +1
    q2max = -(base - 2.0 * p1cm * p3cm);   // cos_cm = -1
    // Do not clamp the forward endpoint to zero. For unequal initial/final
    // masses, physical downscattering can have a small timelike-transfer
    // region where Q^2 = -(p1 - p3)^2 is negative.
    return q2max > q2min;
}

// The propagator-shaped proposal 1/(Q^2 + m_med^2)^2 is not normalizable
// across the t-channel pole at Q^2 = -m_med^2, which lies inside the
// kinematic window for mass-changing scattering with a light mediator.
// Sample and Density apply the same restriction, so the sliver below the
// pole carries zero proposal density from this channel (a mixture's
// physical channel covers it).  Returns false when no part of the window
// lies above the pole.
bool RestrictToAbovePropagatorPole(
    double mediator_mass_squared, double & q2min, double q2max) {
    double span_above = q2max + mediator_mass_squared;
    if (!(span_above > 0.0)) return false;
    // The relative margin keeps Q^2 reconstructed from the record's
    // momenta resolvable against the pole offset; the importance weight
    // is flat in the offset because the physical amplitude shares the
    // propagator-squared shape.
    double q2_floor = -mediator_mass_squared + 1e-9 * span_above;
    if (q2min < q2_floor) q2min = q2_floor;
    return q2max > q2min;
}

// Construct a unit direction at polar angle theta (cos_theta) from the
// beam axis, with the given azimuth phi about that axis.
siren::math::Vector3D DirectionAtPolarAngle(
    siren::math::Vector3D const & beam, double cos_theta, double phi) {
    siren::math::Vector3D z = beam.normalized();
    siren::math::Vector3D x, y;
    detail::SectorPerpFrame(z, x, y);
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    return (z * cos_theta
            + x * (sin_theta * std::cos(phi))
            + y * (sin_theta * std::sin(phi))).normalized();
}

} // anonymous namespace

DetectorDirectedScatteringChannel::DetectorDirectedScatteringChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int directed_index,
    Variable variable,
    DetectorDirected2BodyChannel::Mode mode,
    Q2Mode q2_mode,
    double mediator_mass,
    std::vector<double> q2_cdf_nodes,
    std::vector<double> q2_cdf_values,
    double volume)
    : target_(std::move(target))
    , directed_index_(directed_index)
    , variable_(variable)
    , mode_(mode)
    , q2_mode_(q2_mode)
    , mediator_mass_(mediator_mass)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirectedScatteringChannel requires a target geometry");
    }
    target_volume_ = ResolveDetectorDirectedVolume(
        *target_, mode_ == DetectorDirected2BodyChannel::Mode::Volume, volume);

    if (q2_mode_ != Q2Mode::Geometry) {
        if (variable_ != Variable::Q2) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: Propagator/Tabulated Q2 mode "
                "requires Variable::Q2");
        }
        if (q2_mode_ == Q2Mode::Propagator && mediator_mass_ <= 0.0) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: Propagator Q2 mode requires "
                "mediator_mass > 0");
        }
        if (q2_mode_ == Q2Mode::Tabulated &&
            (q2_cdf_nodes.size() < 2 || q2_cdf_nodes.size() != q2_cdf_values.size())) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: Tabulated Q2 mode requires "
                "matching q2_cdf_nodes / q2_cdf_values arrays of length >= 2");
        }
        if (q2_mode_ == Q2Mode::Tabulated) {
            q2_cdf_table_ = std::make_shared<TabulatedMappingTable>(
                std::move(q2_cdf_nodes), std::move(q2_cdf_values));
        }
    }
}

void DetectorDirectedScatteringChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

PhaseSpaceMeasure DetectorDirectedScatteringChannel::Measure() const {
    if (variable_ == Variable::Q2) {
        return PhaseSpaceMeasure::MandelstamQ2Phi();
    }
    // BjorkenY and RecoilY are one-dimensional fixed-mass y variables. They
    // differ by an affine map with unit absolute slope and both satisfy
    // |dQ2/dy| = 2 M_target E_in. They are not the DIS dx dy measure.
    return PhaseSpaceMeasure::FixedMassYPhi();
}

bool DetectorDirectedScatteringChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const
{
    if (q2_mode_ != Q2Mode::Geometry) return true;
    if (!detail::HasSecondaryStorage(record, 2) ||
        directed_index_ < 0 || directed_index_ > 1) {
        return false;
    }

    int other_index = 1 - directed_index_;
    FourVector p1 = ReadPrimary(record);
    double m2 = record.target_mass;
    double m3 = record.secondary_masses[directed_index_];
    double m4 = record.secondary_masses[other_index];
    FourVector total{p1.e + m2, p1.p};
    double s = total.e * total.e - total.p * total.p;
    if (s <= 0.0 || std::sqrt(s) <= m3 + m4) return false;

    siren::math::Vector3D vertex = detail::ReadVertex(record);
    auto geometry = detail::ClassifyDirectedRegime(
        total.e, total.p.GetX(), total.p.GetY(), total.p.GetZ(),
        std::sqrt(s), m3, m4, vertex, *target_);
    return detail::IsDirectedStepActive(
        geometry.regime, geometry.inside_geometry);
}

void DetectorDirectedScatteringChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    detail::RequireSecondaryStorage(
        record, 2, "DetectorDirectedScatteringChannel");
    if (directed_index_ < 0 || directed_index_ > 1) {
        throw std::runtime_error("directed_index must be 0 or 1");
    }

    if (q2_mode_ != Q2Mode::Geometry) {
        SampleFromQ2Mapping(random, record);
        return;
    }

    int other_index = 1 - directed_index_;
    FourVector p1 = ReadPrimary(record);
    double m1 = record.primary_mass;
    double m2 = record.target_mass;
    double m3 = record.secondary_masses[directed_index_];
    double m4 = record.secondary_masses[other_index];

    FourVector p2{m2, siren::math::Vector3D(0, 0, 0)};
    FourVector total{p1.e + p2.e, p1.p + p2.p};
    double s = total.e * total.e - total.p * total.p;
    if (s <= 0.0 || std::sqrt(s) <= m3 + m4) {
        throw siren::utilities::InjectionFailure(
            "DetectorDirectedScatteringChannel: no physical 2->2 phase space");
    }

    siren::math::Vector3D vertex = detail::ReadVertex(record);

    // A target-at-rest 2->2 scatter is a two-body decay of the total
    // four-momentum. The shared regime machinery intersects the detector and
    // kinematic cones and uses a deterministic uniform-CM fallback when
    // pointing is impossible or adds no support restriction.
    auto step = detail::SampleDirectedStep(
        total.e, total.p.GetX(), total.p.GetY(), total.p.GetZ(),
        std::sqrt(s), m3, m4,
        vertex, *target_, target_volume_, mode_, random,
        detail::DirectedBranchSelection::Uniform);

    FourVector directed{step.E_lab, step.lab_dir * step.p_lab};
    FourVector other{total.e - directed.e, total.p - directed.p};

    WriteSecondary(record, directed_index_, directed);
    WriteSecondary(record, other_index, other);

    double Q2 = ComputeQ2(p1, directed, m1, m3);
    double energy_loss_y = 1.0 - directed.e / p1.e;
    double recoil_y = (directed.e - m3) / p1.e;
    record.interaction_parameters["energy"] = p1.e;
    record.interaction_parameters["Q2"] = Q2;
    record.interaction_parameters["bjorken_y"] = energy_loss_y;
    record.interaction_parameters["recoil_y"] = recoil_y;
}

double DetectorDirectedScatteringChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 2) ||
        directed_index_ < 0 || directed_index_ > 1) {
        return 0.0;
    }

    if (q2_mode_ != Q2Mode::Geometry) {
        return DensityFromQ2Mapping(record);
    }

    int other_index = 1 - directed_index_;
    FourVector p1 = ReadPrimary(record);
    FourVector directed = ReadSecondary(record, directed_index_);
    double p1_mag = p1.p.magnitude();
    double directed_mag = directed.p.magnitude();
    if (p1.e <= 0.0 || p1_mag <= 0.0 || directed_mag <= 0.0) return 0.0;

    double m1 = record.primary_mass;
    double m2 = record.target_mass;
    double m3 = record.secondary_masses[directed_index_];
    double m4 = record.secondary_masses[other_index];

    FourVector p2{m2, siren::math::Vector3D(0, 0, 0)};
    FourVector total{p1.e + p2.e, p1.p + p2.p};
    double s = total.e * total.e - total.p * total.p;
    if (s <= 0.0 || std::sqrt(s) <= m3 + m4) return 0.0;

    siren::math::Vector3D vertex = detail::ReadVertex(record);

    double rest_density = detail::DensityDirectedStep(
        total.e, total.p.GetX(), total.p.GetY(), total.p.GetZ(),
        std::sqrt(s), m3, m4,
        directed.e,
        directed.p.GetX(), directed.p.GetY(), directed.p.GetZ(),
        vertex, *target_, target_volume_, mode_,
        detail::DirectedBranchSelection::Uniform);
    if (rest_density <= 0.0) return 0.0;

    namespace J = siren::injection::phase_space_jacobian;
    double q2_density = J::SolidAngleRestDensityToMandelstamQ2PhiDensity(
        rest_density, s, m1, m2, m3, m4);
    if (variable_ == Variable::Q2) return q2_density;
    return J::MandelstamQ2DensityToFixedMassYDensity(
        q2_density, m2, p1.e);
}

void DetectorDirectedScatteringChannel::SampleFromQ2Mapping(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::dataclasses::InteractionRecord & record) const
{
    int other_index = 1 - directed_index_;
    FourVector p1 = ReadPrimary(record);
    double m1 = record.primary_mass;
    double m2 = record.target_mass;
    double m3 = record.secondary_masses[directed_index_];
    double m4 = record.secondary_masses[other_index];
    double p1_mag = p1.p.magnitude();

    double q2min = 0.0, q2max = 0.0;
    if (p1_mag <= 0.0 ||
        !Q2Range(p1.e, m1, m2, m3, m4, q2min, q2max)) {
        throw siren::utilities::InjectionFailure(
            "DetectorDirectedScatteringChannel: no valid Q2 range for mapping mode");
    }

    double Q2;
    if (q2_mode_ == Q2Mode::Propagator) {
        double m_med2 = mediator_mass_ * mediator_mass_;
        if (!RestrictToAbovePropagatorPole(m_med2, q2min, q2max)) {
            throw siren::utilities::InjectionFailure(
                "Propagator Q2 window lies entirely at or below the "
                "t-channel pole for this event");
        }
        Q2 = PropagatorMapping(m_med2, q2min, q2max)
                 .Forward(random->Uniform(0.0, 1.0));
    } else {
        TabulatedMapping map(q2_cdf_table_, q2min, q2max);
        if (!map.HasSupport()) {
            throw siren::utilities::InjectionFailure(
                "Tabulated Q2 mapping has no support in the event's kinematic range");
        }
        Q2 = map.Forward(random->Uniform(0.0, 1.0));
    }

    double cos_theta = CosThetaFromQ2(Q2, p1.e, p1_mag, m1, m2, m3, m4);
    if (!std::isfinite(cos_theta)) {
        // Numerical edge: clamp into the physical range.
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
    }
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    double phi = random->Uniform(0.0, detail::kTwoPi);
    siren::math::Vector3D dir = DirectionAtPolarAngle(p1.p, cos_theta, phi);

    // With the polar angle (hence Q^2) fixed, the outgoing energy follows
    // from 2->2 kinematics: E3 = E1 + m2 - E4, E4 = (Q2 + m2^2 + m4^2)/(2 m2).
    double E4 = (Q2 + m2 * m2 + m4 * m4) / (2.0 * m2);
    double E3 = p1.e + m2 - E4;
    double p3_mag = MomentumMagnitude(E3, m3);
    if (p3_mag <= 0.0) {
        throw siren::utilities::InjectionFailure(
            "DetectorDirectedScatteringChannel: unphysical outgoing momentum in mapping mode");
    }

    FourVector p3{E3, dir * p3_mag};
    FourVector p2{m2, siren::math::Vector3D(0, 0, 0)};
    FourVector p4{p1.e + p2.e - p3.e, p1.p + p2.p - p3.p};

    WriteSecondary(record, directed_index_, p3);
    WriteSecondary(record, other_index, p4);

    double Q2_actual = ComputeQ2(p1, p3, m1, m3);
    record.interaction_parameters["energy"] = p1.e;
    record.interaction_parameters["Q2"] = Q2_actual;
    record.interaction_parameters["bjorken_y"] = 1.0 - p3.e / p1.e;
    record.interaction_parameters["recoil_y"] = (p3.e - m3) / p1.e;
}

double DetectorDirectedScatteringChannel::DensityFromQ2Mapping(
    siren::dataclasses::InteractionRecord const & record) const
{
    int other_index = 1 - directed_index_;
    FourVector p1 = ReadPrimary(record);
    FourVector directed = ReadSecondary(record, directed_index_);
    if (p1.e <= 0.0 || p1.p.magnitude() <= 0.0 || directed.p.magnitude() <= 0.0) {
        return 0.0;
    }

    double m1 = record.primary_mass;
    double m2 = record.target_mass;
    double m3 = record.secondary_masses[directed_index_];
    double m4 = record.secondary_masses[other_index];

    double q2min = 0.0, q2max = 0.0;
    if (!Q2Range(p1.e, m1, m2, m3, m4, q2min, q2max)) return 0.0;

    double Q2 = ComputeQ2(p1, directed, m1, m3);
    if (Q2 < q2min || Q2 > q2max) return 0.0;

    // The channel reports the full dQ2 dphi proposal density. The mapping is
    // normalized in Q2 and SampleFromQ2Mapping draws a uniform azimuth, so the
    // joint density carries the explicit 1/(2*pi) conditional.
    if (q2_mode_ == Q2Mode::Propagator) {
        double m_med2 = mediator_mass_ * mediator_mass_;
        if (!RestrictToAbovePropagatorPole(m_med2, q2min, q2max))
            return 0.0;
        return PropagatorMapping(m_med2, q2min, q2max)
                   .Density(Q2) / detail::kTwoPi;
    }
    return TabulatedMapping(q2_cdf_table_, q2min, q2max)
               .Density(Q2) / detail::kTwoPi;
}

} // namespace injection
} // namespace siren
