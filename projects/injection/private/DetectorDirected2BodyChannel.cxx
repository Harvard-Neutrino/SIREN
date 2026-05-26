#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace siren {
namespace injection {

static const double TWO_PI = 2.0 * M_PI;
static const double FOUR_PI = 4.0 * M_PI;

namespace {

double ComputeVolume(siren::geometry::Geometry const & geo, double aabb_volume) {
    // Analytic volume for known geometry types
    if (auto const * cyl = dynamic_cast<siren::geometry::Cylinder const *>(&geo)) {
        double R = cyl->GetRadius();
        double r = cyl->GetInnerRadius();
        double H = cyl->GetZ();
        return M_PI * (R * R - r * r) * H;
    }
    if (auto const * sph = dynamic_cast<siren::geometry::Sphere const *>(&geo)) {
        double R = sph->GetRadius();
        double r = sph->GetInnerRadius();
        return (4.0 / 3.0) * M_PI * (R * R * R - r * r * r);
    }
    if (auto const * box = dynamic_cast<siren::geometry::Box const *>(&geo)) {
        return box->GetX() * box->GetY() * box->GetZ();
    }

    // Monte Carlo fallback: estimate V from AABB rejection rate
    auto aabb = geo.GetWorldBoundingBox();
    int N = 10000;
    int n_inside = 0;
    // Deterministic scan (not random — no RNG available here)
    for (int i = 0; i < N; ++i) {
        double fi = (i + 0.5) / N;
        // 3D Halton-like quasi-random sequence
        int ix = i;
        double fx = 0, fy = 0, fz = 0;
        double bx = 1.0/2, by = 1.0/3, bz = 1.0/5;
        for (int j = 0; j < 15; ++j) {
            fx += (ix % 2) * bx; bx /= 2; ix /= 2;
        }
        ix = i;
        for (int j = 0; j < 10; ++j) {
            fy += (ix % 3) * by; by /= 3; ix /= 3;
        }
        ix = i;
        for (int j = 0; j < 7; ++j) {
            fz += (ix % 5) * bz; bz /= 5; ix /= 5;
        }
        double x = aabb.min_corner.GetX() + fx * (aabb.max_corner.GetX() - aabb.min_corner.GetX());
        double y = aabb.min_corner.GetY() + fy * (aabb.max_corner.GetY() - aabb.min_corner.GetY());
        double z = aabb.min_corner.GetZ() + fz * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ());
        if (geo.IsInside(siren::math::Vector3D(x, y, z))) {
            ++n_inside;
        }
    }
    if (n_inside == 0) return aabb_volume;
    return aabb_volume * static_cast<double>(n_inside) / N;
}

} // anonymous namespace

DetectorDirected2BodyChannel::DetectorDirected2BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int daughter_index,
    Mode mode)
    : target_(std::move(target))
    , daughter_index_(daughter_index)
    , mode_(mode)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected2BodyChannel requires a target geometry");
    }
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    aabb_volume_ = extent.GetX() * extent.GetY() * extent.GetZ();
    target_volume_ = ComputeVolume(*target_, aabb_volume_);

    if (mode_ == Mode::Volume) {
        if (aabb_volume_ <= 0.0) {
            throw std::runtime_error("Target bounding box has zero volume");
        }
        if (target_volume_ / aabb_volume_ < 1e-4) {
            throw std::runtime_error("Target volume is too small relative to its bounding box for sampling to be viable");
        }
    }
}

void DetectorDirected2BodyChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

// ================================================================ //
//  Cone mode helpers                                                 //
// ================================================================ //

double DetectorDirected2BodyChannel::ConeSolidAngle(
    siren::math::Vector3D const & position) const
{
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D offset = center - position;
    double dist = offset.magnitude();

    if (dist < bounding_radius) return FOUR_PI;

    double sin_theta = bounding_radius / dist;
    double cos_theta = std::sqrt(1.0 - sin_theta * sin_theta);
    return TWO_PI * (1.0 - cos_theta);
}

std::pair<siren::math::Vector3D, double>
DetectorDirected2BodyChannel::SampleConeDirection(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position) const
{
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D to_center = center - position;
    double dist = to_center.magnitude();

    if (dist < bounding_radius) {
        double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
        double phi = TWO_PI * random->Uniform(0, 1);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        return {
            siren::math::Vector3D(
                sin_theta * std::cos(phi),
                sin_theta * std::sin(phi),
                cos_theta),
            FOUR_PI
        };
    }

    double sin_cone = bounding_radius / dist;
    double cos_cone = std::sqrt(1.0 - sin_cone * sin_cone);
    double solid_angle = TWO_PI * (1.0 - cos_cone);

    double cos_theta = cos_cone + (1.0 - cos_cone) * random->Uniform(0, 1);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = TWO_PI * random->Uniform(0, 1);

    siren::math::Vector3D axis = to_center;
    axis.normalize();

    siren::math::Vector3D perp1, perp2;
    if (std::abs(axis.GetX()) < 0.9)
        perp1 = siren::math::Vector3D(1, 0, 0);
    else
        perp1 = siren::math::Vector3D(0, 1, 0);
    perp2 = siren::math::vector_product(axis, perp1);
    perp2.normalize();
    perp1 = siren::math::vector_product(perp2, axis);
    perp1.normalize();

    siren::math::Vector3D direction = axis * cos_theta
        + perp1 * (sin_theta * std::cos(phi))
        + perp2 * (sin_theta * std::sin(phi));
    direction.normalize();

    return {direction, solid_angle};
}

bool DetectorDirected2BodyChannel::DirectionHitsTarget(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const
{
    auto intersections = target_->Intersections(position, direction);
    for (auto const & isec : intersections) {
        if (isec.distance > 0 && isec.entering) return true;
    }
    return false;
}

// ================================================================ //
//  Volume mode helpers                                               //
// ================================================================ //

siren::math::Vector3D DetectorDirected2BodyChannel::SampleVolumePoint(
    std::shared_ptr<siren::utilities::SIREN_random> random) const
{
    auto aabb = target_->GetWorldBoundingBox();

    for (int attempt = 0; attempt < 10000; ++attempt) {
        double x = aabb.min_corner.GetX()
            + random->Uniform(0, 1) * (aabb.max_corner.GetX() - aabb.min_corner.GetX());
        double y = aabb.min_corner.GetY()
            + random->Uniform(0, 1) * (aabb.max_corner.GetY() - aabb.min_corner.GetY());
        double z = aabb.min_corner.GetZ()
            + random->Uniform(0, 1) * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ());
        siren::math::Vector3D point(x, y, z);
        if (target_->IsInside(point)) return point;
    }
    throw std::runtime_error("Failed to sample a point inside target volume after 10000 attempts");
}

double DetectorDirected2BodyChannel::SolidAngleDensity(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const
{
    auto intersections = target_->Intersections(position, direction);

    std::vector<siren::geometry::Geometry::Intersection> pos_intersections;
    for (auto const & isec : intersections) {
        if (isec.distance > 0) {
            pos_intersections.push_back(isec);
        }
    }
    std::sort(pos_intersections.begin(), pos_intersections.end(),
              [](auto const & a, auto const & b) { return a.distance < b.distance; });

    // g(omega) = (1/V) * integral of r^2 dr over intersection segments
    //          = (1/V) * sum_segments (r_exit^3 - r_enter^3) / 3
    double r2_integral = 0.0;
    if (!pos_intersections.empty()) {
        size_t start_idx = 0;
        if (!pos_intersections[0].entering) {
            // Origin is inside the volume.
            // First segment is from 0 to the first exit.
            double r_exit = pos_intersections[0].distance;
            r2_integral += (r_exit * r_exit * r_exit) / 3.0;
            start_idx = 1;
        }

        for (size_t i = start_idx; i < pos_intersections.size(); ++i) {
            if (pos_intersections[i].entering) {
                double r_enter = pos_intersections[i].distance;
                double r_exit = r_enter;
                for (size_t j = i + 1; j < pos_intersections.size(); ++j) {
                    if (!pos_intersections[j].entering) {
                        r_exit = pos_intersections[j].distance;
                        break;
                    }
                }
                if (r_exit > r_enter) {
                    r2_integral += (r_exit * r_exit * r_exit
                                  - r_enter * r_enter * r_enter) / 3.0;
                }
            }
        }
    }

    if (r2_integral <= 0) return 0.0;
    return r2_integral / target_volume_;
}

// ================================================================ //
//  Sample and Density                                                //
// ================================================================ //

void DetectorDirected2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 2) {
        throw std::runtime_error(
            "DetectorDirected2BodyChannel requires exactly 2 secondaries");
    }

    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_A_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];
    double p_parent = std::sqrt(
        px_parent * px_parent + py_parent * py_parent + pz_parent * pz_parent);

    if (p_parent < 1e-15) {
        Isotropic2BodyChannel fallback(daughter_index_);
        fallback.Sample(random, nullptr, record);
        return;
    }

    double beta = p_parent / E_parent;
    double gamma = E_parent / M_parent;

    siren::math::Vector3D parent_dir(
        px_parent / p_parent, py_parent / p_parent, pz_parent / p_parent);

    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    // Compute the kinematic cone: the maximum lab angle the daughter
    // can reach.  Beyond this angle, no rest-frame solution exists.
    double cos_critical = CriticalCosTheta(
        beta, gamma, p_rest, E_A_rest, m_A);

    // Sample a lab direction toward the target, rejecting directions
    // outside the kinematic cone.  For long-baseline geometries the
    // target is inside the cone and rejection is rare or never needed.
    siren::math::Vector3D lab_dir;
    int n_valid = 0;
    std::array<TwoBodyLabSolution, 2> solutions;

    for (int attempt = 0; attempt < 1000; ++attempt) {
        if (mode_ == Mode::Cone) {
            auto [dir, sa] = SampleConeDirection(random, decay_pos);
            lab_dir = dir;
        } else {
            siren::math::Vector3D target_point;
            siren::math::Vector3D diff;
            double diff_mag = 0.0;
            for (int inner = 0; inner < 100; ++inner) {
                target_point = SampleVolumePoint(random);
                diff = target_point - decay_pos;
                diff_mag = diff.magnitude();
                if (diff_mag > 1e-6) break;
            }
            if (diff_mag <= 1e-6) continue;
            lab_dir = diff;
            lab_dir.normalize();
        }

        double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);
        if (cos_theta_lab < cos_critical) continue;

        solutions = SolveLabAngle(
            beta, gamma, p_rest, E_A_rest, m_A, cos_theta_lab);
        n_valid = (solutions[0].valid ? 1 : 0)
                + (solutions[1].valid ? 1 : 0);
        if (n_valid > 0) break;
    }

    if (n_valid == 0) {
        Isotropic2BodyChannel fallback(daughter_index_);
        fallback.Sample(random, nullptr, record);
        return;
    }

    int chosen = 0;
    if (n_valid == 2) {
        double w0 = solutions[0].jacobian;
        double w1 = solutions[1].jacobian;
        double r = random->Uniform(0, 1) * (w0 + w1);
        chosen = (r < w0) ? 0 : 1;
    } else {
        chosen = solutions[0].valid ? 0 : 1;
    }

    double p_lab = solutions[chosen].p_lab;
    double E_A_lab = std::sqrt(p_lab * p_lab + m_A * m_A);

    record.secondary_momenta[daughter_index_] = {
        E_A_lab,
        p_lab * lab_dir.GetX(),
        p_lab * lab_dir.GetY(),
        p_lab * lab_dir.GetZ()
    };

    record.secondary_momenta[1 - daughter_index_] = {
        E_parent - E_A_lab,
        px_parent - p_lab * lab_dir.GetX(),
        py_parent - p_lab * lab_dir.GetY(),
        pz_parent - p_lab * lab_dir.GetZ()
    };
}

double DetectorDirected2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    // The density is reported in LabFrameSolidAngle convention.
    //
    // Sample() picks a lab direction from the target geometry,
    // rejection-sampling against the kinematic cone.  The density
    // is the proposal density (g_angular) scaled by the branch
    // selection probability (J_match / J_total) and divided by the
    // acceptance probability (fraction of the proposal that falls
    // within the kinematic cone).
    //
    // For Cone mode:  g_lab = (1/Omega_effective) * J_match / J_total
    //   where Omega_effective = min(Omega_cone, Omega_kinematic) when
    //   the target and parent directions are nearly aligned.
    // For Volume mode: g_lab = SolidAngleDensity / acceptance * J_match / J_total

    if (record.signature.secondary_types.size() != 2) return 0.0;

    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_A_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];
    double p_parent = std::sqrt(
        px_parent * px_parent + py_parent * py_parent + pz_parent * pz_parent);

    if (p_parent < 1e-15) return 1.0 / FOUR_PI;

    double beta = p_parent / E_parent;
    double gamma = E_parent / M_parent;

    double px_A = record.secondary_momenta[daughter_index_][1];
    double py_A = record.secondary_momenta[daughter_index_][2];
    double pz_A = record.secondary_momenta[daughter_index_][3];
    double p_A = std::sqrt(px_A * px_A + py_A * py_A + pz_A * pz_A);
    if (p_A < 1e-15) return 0.0;

    siren::math::Vector3D lab_dir(px_A / p_A, py_A / p_A, pz_A / p_A);
    siren::math::Vector3D parent_dir(
        px_parent / p_parent, py_parent / p_parent, pz_parent / p_parent);
    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    // Check kinematic cone
    double cos_critical = CriticalCosTheta(
        beta, gamma, p_rest, E_A_rest, m_A);
    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);
    if (cos_theta_lab < cos_critical) return 0.0;

    double g_angular;
    if (mode_ == Mode::Cone) {
        // Bounding cone check
        auto aabb = target_->GetWorldBoundingBox();
        siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
        siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
        double bounding_radius = 0.5 * extent.magnitude();
        siren::math::Vector3D to_center = center - decay_pos;
        double dist = to_center.magnitude();
        if (dist > bounding_radius) {
            to_center.normalize();
            double cos_to_center = siren::math::scalar_product(lab_dir, to_center);
            double cos_cone = std::sqrt(1.0 - (bounding_radius / dist) * (bounding_radius / dist));
            if (cos_to_center < cos_cone) return 0.0;
        }
        // The effective solid angle is the bounding cone clipped by
        // the kinematic cone.  For the typical case where target and
        // parent directions are nearly aligned, this is approximately
        // min(Omega_cone, Omega_kinematic).  Use the kinematic cone
        // solid angle as a cap on the bounding cone.
        double omega_bounding = ConeSolidAngle(decay_pos);
        double omega_kinematic = TWO_PI * (1.0 - cos_critical);
        double omega_effective = std::min(omega_bounding, omega_kinematic);
        g_angular = 1.0 / omega_effective;
    } else {
        g_angular = SolidAngleDensity(decay_pos, lab_dir);
        if (g_angular <= 0) return 0.0;
    }

    auto solutions = SolveLabAngle(
        beta, gamma, p_rest, E_A_rest, m_A, cos_theta_lab);

    double E_A_lab = record.secondary_momenta[daughter_index_][0];
    double p_par_lab = p_A * cos_theta_lab;
    double p_par_rest = gamma * (p_par_lab - beta * E_A_lab);
    double cos_theta_rest_actual = p_par_rest / p_rest;

    double J_total = 0.0;
    for (auto const & sol : solutions) {
        if (sol.valid) J_total += sol.jacobian;
    }
    if (J_total <= 0.0) return 0.0;

    double best_jacobian = 0.0;
    double best_distance = 1e30;
    for (auto const & sol : solutions) {
        if (!sol.valid) continue;
        double d = std::abs(sol.cos_theta_rest - cos_theta_rest_actual);
        if (d < best_distance) {
            best_distance = d;
            best_jacobian = sol.jacobian;
        }
    }

    if (best_distance > 0.01) return 0.0;
    return g_angular * best_jacobian / J_total;
}

} // namespace injection
} // namespace siren
