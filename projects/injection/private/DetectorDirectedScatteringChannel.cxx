#include "SIREN/injection/DetectorDirectedScatteringChannel.h"

#include "DetectorDirectedChannelUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
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

} // anonymous namespace

DetectorDirectedScatteringChannel::DetectorDirectedScatteringChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int directed_index,
    Variable variable,
    DetectorDirected2BodyChannel::Mode mode)
    : target_(std::move(target))
    , directed_index_(directed_index)
    , variable_(variable)
    , mode_(mode)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirectedScatteringChannel requires a target geometry");
    }
    double aabb_volume = detail::AABBVolume(*target_);
    target_volume_ = detail::ComputeGeometryVolume(*target_, aabb_volume);
}

void DetectorDirectedScatteringChannel::SetVolume(double volume) {
    target_volume_ = volume;
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

} // namespace injection
} // namespace siren
