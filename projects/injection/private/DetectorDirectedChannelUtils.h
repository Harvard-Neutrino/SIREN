#pragma once
#ifndef SIREN_DetectorDirectedChannelUtils_H
#define SIREN_DetectorDirectedChannelUtils_H

#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <array>
#include <cmath>
#include <memory>
#include <utility>
#include <algorithm>

namespace siren {
namespace injection {
namespace detail {

static const double kTwoPi = 2.0 * M_PI;
static const double kFourPi = 4.0 * M_PI;

inline double AABBVolume(siren::geometry::Geometry const & geo) {
    auto aabb = geo.GetWorldBoundingBox();
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    return extent.GetX() * extent.GetY() * extent.GetZ();
}

inline double ComputeGeometryVolume(siren::geometry::Geometry const & geo, double aabb_volume) {
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

    auto aabb = geo.GetWorldBoundingBox();
    int N = 10000;
    int n_inside = 0;
    for (int i = 0; i < N; ++i) {
        int ix = i;
        double fx = 0.0, fy = 0.0, fz = 0.0;
        double bx = 1.0 / 2.0, by = 1.0 / 3.0, bz = 1.0 / 5.0;
        for (int j = 0; j < 15; ++j) {
            fx += (ix % 2) * bx;
            bx /= 2.0;
            ix /= 2;
        }
        ix = i;
        for (int j = 0; j < 10; ++j) {
            fy += (ix % 3) * by;
            by /= 3.0;
            ix /= 3;
        }
        ix = i;
        for (int j = 0; j < 7; ++j) {
            fz += (ix % 5) * bz;
            bz /= 5.0;
            ix /= 5;
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

inline double ConeSolidAngle(
    siren::geometry::Geometry const & target,
    siren::math::Vector3D const & position)
{
    auto aabb = target.GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D offset = center - position;
    double dist = offset.magnitude();
    if (dist < bounding_radius) return kFourPi;

    double sin_theta = bounding_radius / dist;
    double cos_theta = std::sqrt(std::max(0.0, 1.0 - sin_theta * sin_theta));
    return kTwoPi * (1.0 - cos_theta);
}

inline std::pair<siren::math::Vector3D, double> SampleConeDirection(
    siren::geometry::Geometry const & target,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position)
{
    auto aabb = target.GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double bounding_radius = 0.5 * extent.magnitude();

    siren::math::Vector3D to_center = center - position;
    double dist = to_center.magnitude();

    if (dist < bounding_radius) {
        double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
        double phi = kTwoPi * random->Uniform(0, 1);
        double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        return {
            siren::math::Vector3D(
                sin_theta * std::cos(phi),
                sin_theta * std::sin(phi),
                cos_theta),
            kFourPi
        };
    }

    double sin_cone = bounding_radius / dist;
    double cos_cone = std::sqrt(std::max(0.0, 1.0 - sin_cone * sin_cone));
    double solid_angle = kTwoPi * (1.0 - cos_cone);

    double cos_theta = cos_cone + (1.0 - cos_cone) * random->Uniform(0, 1);
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double phi = kTwoPi * random->Uniform(0, 1);

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

inline bool DirectionHitsTarget(
    siren::geometry::Geometry const & target,
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction)
{
    auto intersections = target.Intersections(position, direction);
    for (auto const & isec : intersections) {
        if (isec.distance > 0 && isec.entering) return true;
    }
    return false;
}

inline siren::math::Vector3D SampleVolumePoint(
    siren::geometry::Geometry const & target,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    auto aabb = target.GetWorldBoundingBox();
    for (int attempt = 0; attempt < 10000; ++attempt) {
        double x = aabb.min_corner.GetX()
            + random->Uniform(0, 1) * (aabb.max_corner.GetX() - aabb.min_corner.GetX());
        double y = aabb.min_corner.GetY()
            + random->Uniform(0, 1) * (aabb.max_corner.GetY() - aabb.min_corner.GetY());
        double z = aabb.min_corner.GetZ()
            + random->Uniform(0, 1) * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ());
        siren::math::Vector3D point(x, y, z);
        if (target.IsInside(point)) return point;
    }
    throw std::runtime_error("Failed to sample a point inside target volume after 10000 attempts");
}

inline double SolidAngleDensity(
    siren::geometry::Geometry const & target,
    double target_volume,
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction)
{
    auto intersections = target.Intersections(position, direction);
    std::vector<siren::geometry::Geometry::Intersection> pos_intersections;
    for (auto const & isec : intersections) {
        if (isec.distance > 0) {
            pos_intersections.push_back(isec);
        }
    }
    std::sort(pos_intersections.begin(), pos_intersections.end(),
              [](auto const & a, auto const & b) { return a.distance < b.distance; });

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

    if (r2_integral <= 0.0 || target_volume <= 0.0) return 0.0;
    return r2_integral / target_volume;
}

inline siren::math::Vector3D SampleDirectedDirection(
    siren::geometry::Geometry const & target,
    DetectorDirected2BodyChannel::Mode mode,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position)
{
    if (mode == DetectorDirected2BodyChannel::Mode::Cone) {
        return SampleConeDirection(target, random, position).first;
    }
    siren::math::Vector3D target_point;
    siren::math::Vector3D direction;
    double diff_mag = 0.0;
    for (int attempt = 0; attempt < 100; ++attempt) {
        target_point = SampleVolumePoint(target, random);
        direction = target_point - position;
        diff_mag = direction.magnitude();
        if (diff_mag > 1e-6) {
            break;
        }
    }
    if (diff_mag <= 1e-6) {
        throw std::runtime_error("Failed to sample a target point sufficiently separated from decay vertex after 100 attempts");
    }
    direction.normalize();
    return direction;
}

inline double AngularDensity(
    siren::geometry::Geometry const & target,
    double target_volume,
    DetectorDirected2BodyChannel::Mode mode,
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction)
{
    if (mode == DetectorDirected2BodyChannel::Mode::Cone) {
        if (!DirectionHitsTarget(target, position, direction)) return 0.0;
        return 1.0 / ConeSolidAngle(target, position);
    }
    return SolidAngleDensity(target, target_volume, position, direction);
}

// ================================================================ //
//  Regime classification for directed 2-body sub-steps              //
// ================================================================ //

enum class DirectedRegime {
    Rest,           // parent at rest (lab = rest frame)
    Disjoint,       // cones don't overlap
    KinInBound,     // kinematic cone fully inside bounding cone
    BoundInKin,     // bounding cone fully inside kinematic cone
    Overlap,        // partial overlap
};

struct DirectedGeometry {
    DirectedRegime regime;
    double theta_bound;
    double theta_kin;
    double axis_sep;
    double omega_bound;
    double omega_kin;
    double omega_eff;
    siren::math::Vector3D to_center;
    bool inside_geometry;
};

// Solid angle of the intersection of two circular cones with
// half-angles theta1, theta2 and axis separation angle alpha.
// From Mazonka (2012), arXiv:1205.1396, Eqs 45, 48-50.
inline double ConeIntersectionSolidAngle(
    double theta1, double theta2, double alpha)
{
    double omega1 = kTwoPi * (1.0 - std::cos(theta1));
    double omega2 = kTwoPi * (1.0 - std::cos(theta2));

    if (alpha <= 0.0) return std::min(omega1, omega2);
    if (alpha >= M_PI) return std::max(omega1 + omega2 - kFourPi, 0.0);
    if (theta1 <= 0.0 || theta2 <= 0.0) return 0.0;

    auto clamp01 = [](double x) -> double {
        if (x > 1.0) return 1.0;
        if (x < -1.0) return -1.0;
        return x;
    };

    auto segment = [&](double theta, double theta_other) -> double {
        double ct = std::cos(theta);
        double st = std::sin(theta);
        double co = std::cos(theta_other);
        double ca = std::cos(alpha);
        double sa = std::sin(alpha);

        double ty = co - ca * ct;
        double tx = sa * ct;

        double cos_phi = clamp01((ty * ct) / (tx * st));
        double cos_beta = clamp01(ty / (st * std::sqrt(tx * tx + ty * ty)));
        double phi = std::acos(cos_phi);
        double beta = std::acos(cos_beta);
        return 2.0 * (beta - phi * ct);
    };

    return segment(theta1, theta2) + segment(theta2, theta1);
}

// Classify the geometric relationship between the bounding cone
// (around the target) and the kinematic cone (around the parent
// direction) for a single 2-body sub-step.
inline DirectedGeometry ClassifyDirectedRegime(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    siren::math::Vector3D const & decay_pos,
    siren::geometry::Geometry const & target)
{
    DirectedGeometry geo;
    geo.inside_geometry = false;

    double p_rest = TwoBodyRestMomentum(parent_mass, daughter_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, daughter_mass, other_mass);
    double p_parent = std::sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);

    // Parent at rest
    if (p_parent < 1e-15) {
        geo.regime = DirectedRegime::Rest;
        geo.inside_geometry = target.IsInside(decay_pos);
        if (geo.inside_geometry) {
            geo.theta_bound = M_PI;
            geo.omega_bound = kFourPi;
        } else {
            auto aabb = target.GetWorldBoundingBox();
            siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
            siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
            double br = 0.5 * extent.magnitude();
            geo.to_center = center - decay_pos;
            double dist = geo.to_center.magnitude();
            geo.to_center.normalize();
            geo.theta_bound = std::asin(std::min(br / dist, 1.0));
            geo.omega_bound = kTwoPi * (1.0 - std::cos(geo.theta_bound));
        }
        geo.theta_kin = M_PI;
        geo.omega_kin = kFourPi;
        geo.axis_sep = 0;
        geo.omega_eff = geo.omega_bound;
        return geo;
    }

    double beta = p_parent / parent_E;
    double gamma = parent_E / parent_mass;
    siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);

    // Kinematic cone
    double cos_critical = CriticalCosTheta(beta, gamma, p_rest, E_rest, daughter_mass);
    geo.theta_kin = std::acos(cos_critical);
    geo.omega_kin = kTwoPi * (1.0 - cos_critical);

    // Bounding cone
    auto aabb = target.GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double br = 0.5 * extent.magnitude();
    geo.to_center = center - decay_pos;
    double dist = geo.to_center.magnitude();

    if (dist <= br) {
        geo.inside_geometry = true;
        geo.theta_bound = M_PI;
        geo.omega_bound = kFourPi;
        geo.axis_sep = 0;
    } else {
        geo.to_center.normalize();
        geo.theta_bound = std::asin(br / dist);
        geo.omega_bound = kTwoPi * (1.0 - std::cos(geo.theta_bound));
        double cos_sep = siren::math::scalar_product(geo.to_center, parent_dir);
        geo.axis_sep = std::acos(std::clamp(cos_sep, -1.0, 1.0));
    }

    // Intersection: if the kinematic cone is the full sphere
    // (all lab angles accessible), the intersection is the bounding cone.
    if (cos_critical <= -1.0 + 1e-10) {
        geo.omega_eff = geo.omega_bound;
    } else {
        geo.omega_eff = ConeIntersectionSolidAngle(
            geo.theta_bound, geo.theta_kin, geo.axis_sep);
    }

    // Classify
    double tol = 1e-12;
    if (geo.omega_eff <= tol) {
        geo.regime = DirectedRegime::Disjoint;
    } else if (geo.omega_eff >= geo.omega_kin - tol) {
        geo.regime = DirectedRegime::KinInBound;
    } else if (geo.omega_eff >= geo.omega_bound - tol) {
        geo.regime = DirectedRegime::BoundInKin;
    } else {
        geo.regime = DirectedRegime::Overlap;
    }

    return geo;
}

// ================================================================ //
//  Per-step sample and density for directed 2-body sub-problems    //
// ================================================================ //

struct DirectedStepResult {
    siren::math::Vector3D lab_dir;
    double p_lab;
    double E_lab;
    double rest_density;
};

// Sample one directed 2-body sub-step and return the lab direction,
// momentum, and rest-frame density.  Handles all five regimes.
inline DirectedStepResult SampleDirectedStep(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    siren::math::Vector3D const & decay_pos,
    siren::geometry::Geometry const & target,
    double target_volume,
    DetectorDirected2BodyChannel::Mode mode,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    DirectedStepResult result;
    double p_rest = TwoBodyRestMomentum(parent_mass, daughter_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, daughter_mass, other_mass);
    double p_parent = std::sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);

    DirectedGeometry geo = ClassifyDirectedRegime(
        parent_E, parent_px, parent_py, parent_pz,
        parent_mass, daughter_mass, other_mass,
        decay_pos, target);

    // ---- Isotropic regimes: sample uniformly in rest frame ----
    if (geo.regime == DirectedRegime::Disjoint ||
        geo.regime == DirectedRegime::KinInBound ||
        (geo.regime == DirectedRegime::Rest && geo.inside_geometry)) {

        double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
        double phi = kTwoPi * random->Uniform(0, 1);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

        siren::math::Vector3D rest_dir(
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta);

        if (p_parent < 1e-15) {
            result.lab_dir = rest_dir;
            result.p_lab = p_rest;
        } else {
            double beta = p_parent / parent_E;
            double gamma = parent_E / parent_mass;
            siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);

            // Build perpendicular frame for parent direction
            siren::math::Vector3D perp1, perp2;
            if (std::abs(parent_dir.GetX()) < 0.9)
                perp1 = siren::math::Vector3D(1, 0, 0);
            else
                perp1 = siren::math::Vector3D(0, 1, 0);
            perp2 = siren::math::vector_product(parent_dir, perp1);
            perp2.normalize();
            perp1 = siren::math::vector_product(perp2, parent_dir);
            perp1.normalize();

            // Rest-frame momentum in parent's frame
            double p_rest_par = p_rest * cos_theta;
            double p_rest_perp1 = p_rest * sin_theta * std::cos(phi);
            double p_rest_perp2 = p_rest * sin_theta * std::sin(phi);

            // Boost to lab
            double E_lab_val = gamma * (E_rest + beta * p_rest_par);
            double p_lab_par = gamma * (p_rest_par + beta * E_rest);

            siren::math::Vector3D p_lab_vec =
                parent_dir * p_lab_par + perp1 * p_rest_perp1 + perp2 * p_rest_perp2;
            result.p_lab = p_lab_vec.magnitude();
            result.lab_dir = (result.p_lab > 1e-15) ? p_lab_vec / result.p_lab
                           : parent_dir;
        }

        result.E_lab = std::sqrt(result.p_lab * result.p_lab + daughter_mass * daughter_mass);
        result.rest_density = 1.0 / kFourPi;
        return result;
    }

    // ---- Directed regimes ----
    double beta = p_parent / parent_E;
    double gamma = parent_E / parent_mass;
    siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);
    double cos_critical = CriticalCosTheta(beta, gamma, p_rest, E_rest, daughter_mass);

    siren::math::Vector3D lab_dir;
    std::array<TwoBodyLabSolution, 2> solutions;
    int n_valid = 0;

    if (geo.regime == DirectedRegime::Rest) {
        // Rest, outside geometry: sample from bounding cone
        auto [dir, sa] = SampleConeDirection(target, random, decay_pos);
        lab_dir = dir;
        n_valid = 1;
        solutions[0].valid = true;
        solutions[0].p_lab = p_rest;
        solutions[0].jacobian = 1.0;
        solutions[0].cos_theta_rest = siren::math::scalar_product(lab_dir, geo.to_center);
    } else if (geo.regime == DirectedRegime::BoundInKin) {
        if (mode == DetectorDirected2BodyChannel::Mode::Volume) {
            siren::math::Vector3D diff;
            double diff_mag = 0.0;
            for (int attempt = 0; attempt < 100; ++attempt) {
                siren::math::Vector3D pt = SampleVolumePoint(target, random);
                diff = pt - decay_pos;
                diff_mag = diff.magnitude();
                if (diff_mag > 1e-6) break;
            }
            if (diff_mag > 1e-6) {
                lab_dir = diff / diff_mag;
            } else {
                lab_dir = SampleConeDirection(target, random, decay_pos).first;
            }
        } else {
            lab_dir = SampleConeDirection(target, random, decay_pos).first;
        }
        double cos_lab = siren::math::scalar_product(lab_dir, parent_dir);
        solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_lab);
        n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);
    } else {
        // Overlap: rejection from smaller cone
        bool bound_smaller = geo.omega_bound <= geo.omega_kin;
        for (int attempt = 0; attempt < 10000; ++attempt) {
            if (bound_smaller) {
                lab_dir = SampleConeDirection(target, random, decay_pos).first;
                double cos_lab = siren::math::scalar_product(lab_dir, parent_dir);
                if (cos_lab < cos_critical) continue;
            } else {
                double cos_theta = cos_critical + (1.0 - cos_critical) * random->Uniform(0, 1);
                double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                double phi_val = kTwoPi * random->Uniform(0, 1);
                siren::math::Vector3D perp1, perp2;
                if (std::abs(parent_dir.GetX()) < 0.9)
                    perp1 = siren::math::Vector3D(1, 0, 0);
                else
                    perp1 = siren::math::Vector3D(0, 1, 0);
                perp2 = siren::math::vector_product(parent_dir, perp1);
                perp2.normalize();
                perp1 = siren::math::vector_product(perp2, parent_dir);
                perp1.normalize();
                lab_dir = parent_dir * cos_theta
                    + perp1 * (sin_theta * std::cos(phi_val))
                    + perp2 * (sin_theta * std::sin(phi_val));
                lab_dir.normalize();
                double cos_to = siren::math::scalar_product(lab_dir, geo.to_center);
                if (cos_to < std::cos(geo.theta_bound)) continue;
            }
            double cos_lab = siren::math::scalar_product(lab_dir, parent_dir);
            solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_lab);
            n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);
            if (n_valid > 0) break;
        }
    }

    if (n_valid == 0) {
        // Defensive fallback (should not happen in classified regimes)
        result.lab_dir = parent_dir;
        result.p_lab = p_rest;
        result.E_lab = E_rest;
        result.rest_density = 1.0 / kFourPi;
        return result;
    }

    // Choose branch
    int chosen = 0;
    if (n_valid == 2) {
        double w0 = solutions[0].jacobian;
        double w1 = solutions[1].jacobian;
        chosen = (random->Uniform(0, 1) * (w0 + w1) < w0) ? 0 : 1;
    } else {
        chosen = solutions[0].valid ? 0 : 1;
    }

    result.lab_dir = lab_dir;
    result.p_lab = solutions[chosen].p_lab;
    result.E_lab = std::sqrt(result.p_lab * result.p_lab + daughter_mass * daughter_mass);

    // Compute rest-frame density
    double J_chosen = solutions[chosen].jacobian;
    double J_total = 0.0;
    for (auto const & sol : solutions) {
        if (sol.valid) J_total += sol.jacobian;
    }

    double g_angular;
    if (geo.regime == DirectedRegime::Rest) {
        g_angular = 1.0 / geo.omega_bound;
        // At rest J=1, J_total=1, so rest_density = g_angular
        result.rest_density = g_angular;
    } else {
        if (geo.regime == DirectedRegime::BoundInKin &&
            mode == DetectorDirected2BodyChannel::Mode::Volume) {
            g_angular = SolidAngleDensity(target, target_volume, decay_pos, lab_dir);
            if (g_angular <= 0.0) g_angular = 1.0 / geo.omega_bound;
        } else if (geo.regime == DirectedRegime::BoundInKin) {
            g_angular = 1.0 / geo.omega_bound;
        } else {
            g_angular = (geo.omega_eff > 0.0) ? 1.0 / geo.omega_eff : 1.0 / kFourPi;
        }
        result.rest_density = g_angular * J_chosen * J_chosen / J_total;
    }

    return result;
}

// Compute the rest-frame density for one directed 2-body sub-step
// at the phase-space point described by the daughter's observed
// lab-frame 4-momentum.
inline double DensityDirectedStep(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    double daughter_E,
    double daughter_px, double daughter_py, double daughter_pz,
    siren::math::Vector3D const & decay_pos,
    siren::geometry::Geometry const & target,
    double target_volume,
    DetectorDirected2BodyChannel::Mode mode)
{
    double p_rest = TwoBodyRestMomentum(parent_mass, daughter_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, daughter_mass, other_mass);
    double p_parent = std::sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);

    DirectedGeometry geo = ClassifyDirectedRegime(
        parent_E, parent_px, parent_py, parent_pz,
        parent_mass, daughter_mass, other_mass,
        decay_pos, target);

    // ---- Isotropic regimes ----
    if (geo.regime == DirectedRegime::Disjoint ||
        geo.regime == DirectedRegime::KinInBound ||
        (geo.regime == DirectedRegime::Rest && geo.inside_geometry)) {
        return 1.0 / kFourPi;
    }

    double p_A = std::sqrt(daughter_px*daughter_px + daughter_py*daughter_py + daughter_pz*daughter_pz);
    if (p_A < 1e-15) return 0.0;
    siren::math::Vector3D lab_dir(daughter_px/p_A, daughter_py/p_A, daughter_pz/p_A);

    // Rest-frame parent, outside geometry
    if (geo.regime == DirectedRegime::Rest) {
        double cos_to = siren::math::scalar_product(lab_dir, geo.to_center);
        if (cos_to < std::cos(geo.theta_bound)) return 0.0;
        return 1.0 / geo.omega_bound;
    }

    // ---- Boosted directed regimes ----
    siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);
    double beta = p_parent / parent_E;
    double gamma = parent_E / parent_mass;

    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);
    double cos_critical = CriticalCosTheta(beta, gamma, p_rest, E_rest, daughter_mass);
    if (cos_theta_lab < cos_critical) return 0.0;

    // In BoundInKin, Sample() always produces directions inside the
    // bounding cone, so skip the boundary check (avoids floating-point
    // edge rejections).  In Overlap, the direction must be inside both cones.
    if (geo.regime != DirectedRegime::BoundInKin) {
        double cos_to = siren::math::scalar_product(lab_dir, geo.to_center);
        if (cos_to < std::cos(geo.theta_bound)) return 0.0;
    }

    auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_theta_lab);

    double p_par_lab = p_A * cos_theta_lab;
    double p_par_rest = gamma * (p_par_lab - beta * daughter_E);
    double cos_rest_actual = p_par_rest / p_rest;

    double J_total = 0.0;
    for (auto const & sol : solutions) {
        if (sol.valid) J_total += sol.jacobian;
    }
    if (J_total <= 0.0) return 0.0;

    double J_chosen = 0.0;
    double best_dist = 1e30;
    for (auto const & sol : solutions) {
        if (!sol.valid) continue;
        double d = std::abs(sol.cos_theta_rest - cos_rest_actual);
        if (d < best_dist) {
            best_dist = d;
            J_chosen = sol.jacobian;
        }
    }
    if (best_dist > 0.1) return 0.0;

    double g_angular;
    if (geo.regime == DirectedRegime::BoundInKin &&
        mode == DetectorDirected2BodyChannel::Mode::Volume) {
        g_angular = SolidAngleDensity(target, target_volume, decay_pos, lab_dir);
        if (g_angular <= 0.0) g_angular = 1.0 / geo.omega_bound;
    } else if (geo.regime == DirectedRegime::BoundInKin) {
        g_angular = 1.0 / geo.omega_bound;
    } else {
        if (geo.omega_eff <= 0.0) return 0.0;
        g_angular = 1.0 / geo.omega_eff;
    }

    return g_angular * J_chosen * J_chosen / J_total;
}

} // namespace detail
} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirectedChannelUtils_H
