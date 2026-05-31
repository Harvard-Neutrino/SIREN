#include "SIREN/injection/DetectorDirectedScatteringChannel.h"

#include "DetectorDirectedChannelUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace siren {
namespace injection {

namespace {

struct FourVector {
    double e;
    siren::math::Vector3D p;
};

FourVector ReadPrimary(siren::dataclasses::InteractionRecord const & record) {
    return FourVector{
        record.primary_momentum[0],
        siren::math::Vector3D(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3])
    };
}

FourVector ReadSecondary(siren::dataclasses::InteractionRecord const & record, int index) {
    return FourVector{
        record.secondary_momenta[index][0],
        siren::math::Vector3D(
            record.secondary_momenta[index][1],
            record.secondary_momenta[index][2],
            record.secondary_momenta[index][3])
    };
}

void WriteSecondary(
    siren::dataclasses::InteractionRecord & record,
    int index,
    FourVector const & p4)
{
    record.secondary_momenta[index] = {
        p4.e,
        p4.p.GetX(),
        p4.p.GetY(),
        p4.p.GetZ()
    };
}

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

bool SolveDirectedScattering(
    FourVector const & p1,
    double m1,
    double m2,
    double m3,
    double m4,
    siren::math::Vector3D const & direction,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    FourVector & p3,
    FourVector & p4)
{
    if (m2 <= 0.0) return false;

    FourVector p2{m2, siren::math::Vector3D(0, 0, 0)};
    FourVector total{p1.e + p2.e, p1.p + p2.p};
    double total_p2 = total.p * total.p;
    double s = total.e * total.e - total_p2;
    double C = total.p * direction;
    double A = 0.5 * (s + m3 * m3 - m4 * m4);
    double a = total.e * total.e - C * C;
    double b = -2.0 * A * C;
    double c = total.e * total.e * m3 * m3 - A * A;
    double disc = b * b - 4.0 * a * c;
    if (a <= 0.0 || disc < 0.0) return false;

    double sqrt_disc = std::sqrt(std::max(0.0, disc));
    std::array<double, 2> roots = {
        (-b + sqrt_disc) / (2.0 * a),
        (-b - sqrt_disc) / (2.0 * a)
    };

    std::array<double, 2> valid_roots = {0.0, 0.0};
    int n_valid = 0;
    for (double p : roots) {
        if (p <= 0.0) continue;
        double E3 = std::sqrt(p * p + m3 * m3);
        double lhs = total.e * E3 - C * p;
        if (std::abs(lhs - A) > 1e-7 * std::max(1.0, std::abs(A))) continue;
        FourVector candidate3{E3, direction * p};
        FourVector candidate4{total.e - candidate3.e, total.p - candidate3.p};
        double m4_sq = candidate4.e * candidate4.e - candidate4.p * candidate4.p;
        if (std::abs(m4_sq - m4 * m4) > 1e-6 * std::max(1.0, m4 * m4)) continue;
        valid_roots[n_valid++] = p;
    }
    if (n_valid == 0) return false;

    int chosen = 0;
    if (n_valid == 2 && random->Uniform(0, 1) >= 0.5) chosen = 1;

    double p = valid_roots[chosen];
    p3 = FourVector{std::sqrt(p * p + m3 * m3), direction * p};
    p4 = FourVector{total.e - p3.e, total.p - p3.p};
    return true;
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

double Q2FromEnergyLossY(
    double y,
    double E1,
    double m2,
    double m4)
{
    double E3 = E1 * (1.0 - y);
    double E4 = E1 + m2 - E3;
    return 2.0 * m2 * E4 - m2 * m2 - m4 * m4;
}

double CosThetaFromEnergyLossY(
    double y,
    double E1,
    double p1_mag,
    double m1,
    double m2,
    double m3,
    double m4)
{
    double E3 = E1 * (1.0 - y);
    double p3_mag = MomentumMagnitude(E3, m3);
    if (p3_mag <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    double Q2 = Q2FromEnergyLossY(y, E1, m2, m4);
    return (E1 * E3 - 0.5 * (Q2 + m1 * m1 + m3 * m3)) / (p1_mag * p3_mag);
}

double CosThetaFromRecoilY(
    double y,
    double E1,
    double p1_mag,
    double m1,
    double m2,
    double m3,
    double m4)
{
    double E3 = m3 + E1 * y;
    double p3_mag = MomentumMagnitude(E3, m3);
    if (p3_mag <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    double E4 = E1 + m2 - E3;
    double Q2 = 2.0 * m2 * E4 - m2 * m2 - m4 * m4;
    return (E1 * E3 - 0.5 * (Q2 + m1 * m1 + m3 * m3)) / (p1_mag * p3_mag);
}

template <typename F>
double NumericDerivative(F f, double x) {
    double h = std::max(1e-8, 1e-5 * std::max(1.0, std::abs(x)));
    double fp = f(x + h);
    double fm = f(x - h);
    if (std::isfinite(fp) && std::isfinite(fm)) {
        return (fp - fm) / (2.0 * h);
    }
    double f0 = f(x);
    if (std::isfinite(fp) && std::isfinite(f0)) {
        return (fp - f0) / h;
    }
    if (std::isfinite(fm) && std::isfinite(f0)) {
        return (f0 - fm) / h;
    }
    return 0.0;
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
    if (q2min < 0.0) q2min = 0.0;
    return q2max > q2min;
}

// Construct a unit direction at polar angle theta (cos_theta) from the
// beam axis, with the given azimuth phi about that axis.
siren::math::Vector3D DirectionAtPolarAngle(
    siren::math::Vector3D const & beam, double cos_theta, double phi) {
    siren::math::Vector3D z = beam.normalized();
    siren::math::Vector3D ref =
        std::abs(z.GetZ()) < 0.9 ? siren::math::Vector3D(0, 0, 1)
                                 : siren::math::Vector3D(1, 0, 0);
    siren::math::Vector3D x = cross_product(z, ref).normalized();
    siren::math::Vector3D y = cross_product(z, x);
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
    std::vector<double> q2_cdf_values)
    : target_(std::move(target))
    , directed_index_(directed_index)
    , variable_(variable)
    , mode_(mode)
    , q2_mode_(q2_mode)
    , mediator_mass_(mediator_mass)
    , q2_cdf_nodes_(std::move(q2_cdf_nodes))
    , q2_cdf_values_(std::move(q2_cdf_values))
{
    if (!target_) {
        throw std::runtime_error("DetectorDirectedScatteringChannel requires a target geometry");
    }
    double aabb_volume = detail::AABBVolume(*target_);
    target_volume_ = detail::ComputeGeometryVolume(*target_, aabb_volume);

    if (mode_ == DetectorDirected2BodyChannel::Mode::Volume) {
        if (aabb_volume <= 0.0) {
            throw std::runtime_error("Target bounding box has zero volume");
        }
        if (target_volume_ / aabb_volume < 1e-4) {
            throw std::runtime_error("Target volume is too small relative to its bounding box for sampling to be viable");
        }
    }

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
            (q2_cdf_nodes_.size() < 2 || q2_cdf_nodes_.size() != q2_cdf_values_.size())) {
            throw std::runtime_error(
                "DetectorDirectedScatteringChannel: Tabulated Q2 mode requires "
                "matching q2_cdf_nodes / q2_cdf_values arrays of length >= 2");
        }
    }
}

void DetectorDirectedScatteringChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

PhaseSpaceMeasure DetectorDirectedScatteringChannel::Measure() const {
    if (variable_ == Variable::Q2) {
        return PhaseSpaceMeasure::MandelstamQ2();
    }
    return PhaseSpaceMeasure::BjorkenXY();
}

void DetectorDirectedScatteringChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 2 ||
        record.secondary_masses.size() != 2 ||
        record.secondary_momenta.size() != 2) {
        throw std::runtime_error("DetectorDirectedScatteringChannel requires exactly 2 secondaries");
    }
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

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);
    siren::math::Vector3D direction = detail::SampleDirectedDirection(
        *target_, mode_, random, vertex);

    FourVector directed;
    FourVector other;
    if (!SolveDirectedScattering(p1, m1, m2, m3, m4, direction, random,
                                 directed, other)) {
        throw std::runtime_error("DetectorDirectedScatteringChannel could not solve 2->2 kinematics");
    }

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
    if (record.signature.secondary_types.size() != 2 ||
        record.secondary_masses.size() != 2 ||
        record.secondary_momenta.size() != 2 ||
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

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);
    siren::math::Vector3D direction = directed.p / directed_mag;
    double g_angular = detail::AngularDensity(*target_, target_volume_, mode_, vertex, direction);
    if (g_angular <= 0.0) return 0.0;

    double variable = 0.0;
    double derivative = 0.0;
    if (variable_ == Variable::Q2) {
        variable = ComputeQ2(p1, directed, m1, m3);
        derivative = NumericDerivative(
            [&](double x) {
                return CosThetaFromQ2(x, p1.e, p1_mag, m1, m2, m3, m4);
            },
            variable);
    } else if (variable_ == Variable::BjorkenY) {
        variable = 1.0 - directed.e / p1.e;
        derivative = NumericDerivative(
            [&](double x) {
                return CosThetaFromEnergyLossY(x, p1.e, p1_mag, m1, m2, m3, m4);
            },
            variable);
    } else {
        variable = (directed.e - m3) / p1.e;
        derivative = NumericDerivative(
            [&](double x) {
                return CosThetaFromRecoilY(x, p1.e, p1_mag, m1, m2, m3, m4);
            },
            variable);
    }

    (void) variable;
    if (!std::isfinite(derivative)) return 0.0;

    // Existing SIREN scattering probabilities omit the uniform azimuthal
    // factor. Multiply by 2*pi so this density is comparable to dP/dQ2
    // or dP/dy rather than dP/(dvar dphi).
    return detail::kTwoPi * g_angular * std::abs(derivative);
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
        throw std::runtime_error(
            "DetectorDirectedScatteringChannel: no valid Q2 range for mapping mode");
    }

    double Q2;
    if (q2_mode_ == Q2Mode::Propagator) {
        Q2 = PropagatorMapping(mediator_mass_ * mediator_mass_, q2min, q2max)
                 .Forward(random->Uniform(0.0, 1.0));
    } else {
        Q2 = TabulatedMapping(q2_cdf_nodes_, q2_cdf_values_, q2min, q2max)
                 .Forward(random->Uniform(0.0, 1.0));
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
        throw std::runtime_error(
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

    // The measure is the 1-D Q^2 marginal: the azimuth is uniform and
    // already integrated out, so the proposal density in Q^2 is exactly
    // the mapping density -- no Jacobian, no finite-difference noise.
    if (q2_mode_ == Q2Mode::Propagator) {
        return PropagatorMapping(mediator_mass_ * mediator_mass_, q2min, q2max)
                   .Density(Q2);
    }
    return TabulatedMapping(q2_cdf_nodes_, q2_cdf_values_, q2min, q2max)
               .Density(Q2);
}

} // namespace injection
} // namespace siren
