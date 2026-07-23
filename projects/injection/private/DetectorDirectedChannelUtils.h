#pragma once
#ifndef SIREN_DetectorDirectedChannelUtils_H
#define SIREN_DetectorDirectedChannelUtils_H

#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "LorentzBoostUtils.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <algorithm>

namespace siren {
namespace injection {
namespace detail {

static const double kTwoPi = 2.0 * M_PI;
static const double kFourPi = 4.0 * M_PI;

// Deterministic orthonormal frame perpendicular to a unit axis. All directed
// samplers use the same convention so their azimuth references cannot drift.
inline void SectorPerpFrame(
    siren::math::Vector3D const & axis,
    siren::math::Vector3D & perp1,
    siren::math::Vector3D & perp2)
{
    perp1 = std::abs(axis.GetX()) < 0.9
        ? siren::math::Vector3D(1, 0, 0)
        : siren::math::Vector3D(0, 1, 0);
    perp2 = siren::math::vector_product(axis, perp1);
    perp2.normalize();
    perp1 = siren::math::vector_product(perp2, axis);
    perp1.normalize();
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
    SectorPerpFrame(axis, perp1, perp2);

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

// Sample a direction uniformly within a spherical cap centered on
// `axis` with half-angle whose cosine is `cos_half`.
inline siren::math::Vector3D SampleCapDirection(
    siren::math::Vector3D const & axis,
    double cos_half,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    double cos_theta = cos_half + (1.0 - cos_half) * random->Uniform(0, 1);
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double phi = kTwoPi * random->Uniform(0, 1);

    siren::math::Vector3D perp1, perp2;
    SectorPerpFrame(axis, perp1, perp2);

    siren::math::Vector3D dir = axis * cos_theta
        + perp1 * (sin_theta * std::cos(phi))
        + perp2 * (sin_theta * std::sin(phi));
    dir.normalize();
    return dir;
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
inline DirectedGeometry ComputeDirectedRegime(
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
            if (dist > 1e-15) {
                geo.to_center.normalize();
            } else {
                // The axis is arbitrary for full-sphere support.
                geo.to_center = siren::math::Vector3D(0, 0, 1);
            }
            if (dist < br) {
                // SampleConeDirection uses the full sphere here because no
                // forward cone can enclose an AABB whose bounding sphere
                // contains the vertex. Classification and density must use
                // the same support even when the vertex is outside the exact
                // target geometry.
                geo.theta_bound = M_PI;
                geo.omega_bound = kFourPi;
            } else {
                geo.theta_bound = std::asin(std::min(br / dist, 1.0));
                geo.omega_bound = kTwoPi * (1.0 - std::cos(geo.theta_bound));
            }
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
    if (!std::isfinite(cos_critical)) {
        throw std::runtime_error(
            "Directed two-body classification produced a non-finite critical angle");
    }
    cos_critical = std::clamp(cos_critical, -1.0, 1.0);
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

    // Classify.
    // Near-zero overlaps (< 0.1% of the smaller cone) are treated as
    // Disjoint: the directed contribution is negligible and the
    // rejection sampler cannot reliably find the overlap region.
    double tol = 1e-12;
    double min_omega = std::min(geo.omega_kin, geo.omega_bound);
    if (geo.omega_eff <= tol ||
        (min_omega > tol && geo.omega_eff < 1e-3 * min_omega)) {
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

struct DirectedRegimeCacheState {
    bool valid = false;
    std::size_t epoch = 0;
    siren::geometry::Geometry const * target = nullptr;
    std::array<double, 10> inputs{};
    DirectedGeometry geometry{};
    std::size_t misses = 0;
};

inline DirectedRegimeCacheState & DirectedRegimeCache() {
    static thread_local DirectedRegimeCacheState cache;
    return cache;
}

// Cache the most recent classification per thread. Multi-channel density and
// optimizer-credit loops evaluate sibling sector channels consecutively with
// bit-identical inputs, so one classification serves the entire sector family.
inline DirectedGeometry ClassifyDirectedRegime(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    siren::math::Vector3D const & decay_pos,
    siren::geometry::Geometry const & target)
{
    std::array<double, 10> inputs{
        parent_E, parent_px, parent_py, parent_pz,
        parent_mass, daughter_mass, other_mass,
        decay_pos.GetX(), decay_pos.GetY(), decay_pos.GetZ()};
    auto & cache = DirectedRegimeCache();
    bool cache_active = phase_space_detail::EvaluationCacheActive();
    std::size_t epoch = phase_space_detail::EvaluationEpoch();
    if (cache_active && cache.valid && cache.epoch == epoch &&
        cache.target == &target && cache.inputs == inputs) {
        return cache.geometry;
    }
    DirectedGeometry geometry = ComputeDirectedRegime(
        parent_E, parent_px, parent_py, parent_pz,
        parent_mass, daughter_mass, other_mass, decay_pos, target);
    ++cache.misses;
    if (cache_active) {
        cache.geometry = geometry;
        cache.target = &target;
        cache.inputs = inputs;
        cache.epoch = epoch;
        cache.valid = true;
    } else {
        cache.valid = false;
    }
    return geometry;
}

inline void ResetDirectedRegimeCacheForTesting() {
    DirectedRegimeCache() = DirectedRegimeCacheState{};
}

inline std::size_t DirectedRegimeCacheMissesForTesting() {
    return DirectedRegimeCache().misses;
}

// Regimes where the directed 2-body step produces a DIRECTIONAL proposal,
// vs the isotropic 1/4pi fallback (Disjoint / KinInBound / Rest-inside).
// Used to attribute a directed channel's variance into directing vs fallback.
inline bool IsDirectedStepActive(DirectedRegime regime, bool inside_geometry) {
    if (regime == DirectedRegime::BoundInKin) return true;
    if (regime == DirectedRegime::Overlap) return true;
    if (regime == DirectedRegime::Rest && !inside_geometry) return true;
    return false;
}

// ================================================================ //
//  Per-step sample and density for directed 2-body sub-problems    //
// ================================================================ //

// A two-body split P -> A + B only has phase space for finite masses with
// parent_mass >= daughter_mass + other_mass (the acceptance rule shared with
// Isotropic2BodyChannel). The directed sub-steps must reject a sub-threshold
// parent before TwoBodyRestMomentum turns the mass deficit into NaN
// kinematics: samplers fail the attempt, densities report zero support.
inline bool HasTwoBodyStepPhaseSpace(
    double parent_mass, double daughter_mass, double other_mass)
{
    return std::isfinite(parent_mass) &&
           std::isfinite(daughter_mass) &&
           std::isfinite(other_mass) &&
           parent_mass > 0.0 &&
           daughter_mass >= 0.0 &&
           other_mass >= 0.0 &&
           parent_mass >= daughter_mass + other_mass;
}

struct DirectedStepResult {
    siren::math::Vector3D lab_dir;
    double p_lab;
    double E_lab;
    double rest_density;
};

enum class DirectedBranchSelection {
    JacobianWeighted,
    Uniform
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
    std::shared_ptr<siren::utilities::SIREN_random> random,
    DirectedBranchSelection branch_selection =
        DirectedBranchSelection::JacobianWeighted)
{
    if (!HasTwoBodyStepPhaseSpace(parent_mass, daughter_mass, other_mass)) {
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::KinematicallyForbidden,
            "Directed 2-body step has no kinematically allowed phase space");
    }

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
            siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);

            siren::math::Vector3D perp1, perp2;
            SectorPerpFrame(parent_dir, perp1, perp2);
            siren::math::Vector3D p_rest_vec =
                parent_dir * (p_rest * cos_theta)
                + perp1 * (p_rest * sin_theta * std::cos(phi))
                + perp2 * (p_rest * sin_theta * std::sin(phi));
            auto lab = BoostRestFrameToLab(
                parent_E, parent_px, parent_py, parent_pz,
                E_rest, p_rest_vec.GetX(), p_rest_vec.GetY(), p_rest_vec.GetZ());
            siren::math::Vector3D p_lab_vec(lab[1], lab[2], lab[3]);
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
    // The Rest regime populates only the physical solution at index zero.
    // Value-initialize both entries so the unused branch is deterministically
    // invalid before the common Jacobian accumulation below.
    std::array<TwoBodyLabSolution, 2> solutions{};
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
        for (int retry = 0; retry < 100; ++retry) {
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
            if (n_valid > 0) {
                if (mode == DetectorDirected2BodyChannel::Mode::Volume &&
                    SolidAngleDensity(target, target_volume, decay_pos, lab_dir) <= 0.0) {
                    n_valid = 0;
                    continue;
                }
                break;
            }
        }
    } else {
        // Overlap: the lens is contained in each cone, so uniform draws from
        // the smaller cone are a proven superset of the lens; the two
        // half-space tests below make accepted directions uniform on the lens,
        // matching DensityDirectedStep's 1/omega_eff over the same support.
        double cos_bound = std::cos(geo.theta_bound);
        bool bound_smaller = geo.omega_bound <= geo.omega_kin;

        for (int attempt = 0; attempt < 10000; ++attempt) {
            if (bound_smaller) {
                lab_dir = SampleConeDirection(target, random, decay_pos).first;
            } else {
                lab_dir = SampleCapDirection(
                    parent_dir, cos_critical, random);
            }
            double cos_lab = siren::math::scalar_product(lab_dir, parent_dir);
            if (cos_lab < cos_critical) continue;
            double cos_to = siren::math::scalar_product(lab_dir, geo.to_center);
            if (cos_to < cos_bound) continue;

            solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_lab);
            n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);
            if (n_valid > 0) break;
        }
    }

    if (n_valid == 0) {
        // Exhaustion must fail loudly: any junk direction returned here would
        // claim 1/4pi while DensityDirectedStep reports 1/omega_eff over the
        // lens, breaking Sample/Density closure.
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::KinematicallyForbidden,
            "Directed 2-body rejection sampler exhausted its attempts without a "
            "valid lab direction in the cone-intersection lens");
    }

    // Choose branch
    int chosen = 0;
    if (n_valid == 2) {
        if (branch_selection == DirectedBranchSelection::Uniform) {
            chosen = random->Uniform(0, 1) < 0.5 ? 0 : 1;
        } else {
            double w0 = solutions[0].jacobian;
            double w1 = solutions[1].jacobian;
            chosen = (random->Uniform(0, 1) * (w0 + w1) < w0) ? 0 : 1;
        }
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
        } else if (geo.regime == DirectedRegime::BoundInKin) {
            g_angular = 1.0 / geo.omega_bound;
        } else {
            g_angular = (geo.omega_eff > 0.0) ? 1.0 / geo.omega_eff : 1.0 / kFourPi;
        }
        double branch_probability =
            branch_selection == DirectedBranchSelection::Uniform
            ? 1.0 / n_valid
            : J_chosen / J_total;
        result.rest_density = g_angular * J_chosen * branch_probability;
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
    DetectorDirected2BodyChannel::Mode mode,
    DirectedBranchSelection branch_selection =
        DirectedBranchSelection::JacobianWeighted)
{
    if (!HasTwoBodyStepPhaseSpace(parent_mass, daughter_mass, other_mass)) {
        return 0.0;
    }

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

    {
        double cos_to = siren::math::scalar_product(lab_dir, geo.to_center);
        if (cos_to < std::cos(geo.theta_bound)) return 0.0;
    }

    auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_theta_lab);

    double p_par_lab = p_A * cos_theta_lab;
    double p_par_rest = gamma * (p_par_lab - beta * daughter_E);
    double cos_rest_actual = p_par_rest / p_rest;

    double J_total = 0.0;
    int n_valid = 0;
    for (auto const & sol : solutions) {
        if (sol.valid) {
            J_total += sol.jacobian;
            ++n_valid;
        }
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
        if (g_angular <= 0.0) return 0.0;
    } else if (geo.regime == DirectedRegime::BoundInKin) {
        g_angular = 1.0 / geo.omega_bound;
    } else {
        if (geo.omega_eff <= 0.0) return 0.0;
        g_angular = 1.0 / geo.omega_eff;
    }

    double branch_probability =
        branch_selection == DirectedBranchSelection::Uniform
        ? 1.0 / n_valid
        : J_chosen / J_total;
    return g_angular * J_chosen * branch_probability;
}

// ================================================================ //
//  Angular-sector directed 2-body sub-step                          //
// ================================================================ //
//
// Tiles the cone subtended by the target (around the to-center axis)
// into (u, phi) bins, where u in [0,1] is the fraction of the bounding
// half-angle (u=0 the axis, u=1 the cone edge) and phi is the azimuth
// around the axis.  In the active regime (BoundInKin, or Rest with the
// vertex outside the target) each bin is a uniform-in-solid-angle
// proposal with lab-angular density 1/Omega_bin; elsewhere
// (KinInBound / Disjoint / Overlap / Rest-inside) the sector is
// isotropic, mirroring the volume-directed channel where directing is
// inert.  Sample and Density share the bin test, Omega_bin, axis, perp
// frame, and boost Jacobian, so Sample == Density (Contract C1).

struct AngularSectorBin {
    double u_lo;
    double u_hi;
    double phi_lo;
    double phi_hi;
};

inline double WrapSectorPhi(double phi) {
    phi = std::fmod(phi, kTwoPi);
    if (phi < 0.0) phi += kTwoPi;
    return phi;
}

inline bool PhiInAngularSector(double phi, AngularSectorBin const & bin) {
    double width = bin.phi_hi - bin.phi_lo;
    if (width >= kTwoPi) return true;
    if (!(width > 0.0)) return false;
    double offset = WrapSectorPhi(phi - bin.phi_lo);
    return offset <= width + 1e-12;
}

inline siren::math::Vector3D NormalizedSectorAxis(
    siren::math::Vector3D axis)
{
    double magnitude = axis.magnitude();
    if (magnitude > 1e-15) return axis / magnitude;
    // At the AABB center a full-sphere Rest proposal has no preferred axis.
    // Choose one deterministically so sector coordinates still produce unit
    // directions instead of zero three-momenta.
    return siren::math::Vector3D(0.0, 0.0, 1.0);
}

// Regimes where the sector tiling presents a bin proposal (vs isotropic).
inline bool SectorActive(DirectedRegime regime, bool inside_geometry) {
    if (regime == DirectedRegime::BoundInKin) return true;
    if (regime == DirectedRegime::Rest && !inside_geometry) return true;
    return false;
}

// Isotropic rest-frame sample with proper boost to lab; rest_density 1/4pi.
// Used only when the sector proposal is deterministically inactive.
inline DirectedStepResult IsotropicBoostStep(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    DirectedStepResult result;
    double p_rest = TwoBodyRestMomentum(parent_mass, daughter_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, daughter_mass, other_mass);
    double p_parent = std::sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);

    double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
    double phi = kTwoPi * random->Uniform(0, 1);
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

    if (p_parent < 1e-15) {
        result.lab_dir = siren::math::Vector3D(
            sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
        result.p_lab = p_rest;
    } else {
        siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);
        siren::math::Vector3D perp1, perp2;
        SectorPerpFrame(parent_dir, perp1, perp2);
        siren::math::Vector3D p_rest_vec =
            parent_dir * (p_rest * cos_theta)
            + perp1 * (p_rest * sin_theta * std::cos(phi))
            + perp2 * (p_rest * sin_theta * std::sin(phi));
        auto lab = BoostRestFrameToLab(
            parent_E, parent_px, parent_py, parent_pz,
            E_rest, p_rest_vec.GetX(), p_rest_vec.GetY(), p_rest_vec.GetZ());
        siren::math::Vector3D p_lab_vec(lab[1], lab[2], lab[3]);
        result.p_lab = p_lab_vec.magnitude();
        result.lab_dir = (result.p_lab > 1e-15) ? p_lab_vec / result.p_lab : parent_dir;
    }
    result.E_lab = std::sqrt(result.p_lab * result.p_lab + daughter_mass * daughter_mass);
    result.rest_density = 1.0 / kFourPi;
    return result;
}

inline DirectedStepResult SampleAngularSectorStep(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    siren::math::Vector3D const & decay_pos,
    siren::geometry::Geometry const & target,
    AngularSectorBin const & bin,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    if (!HasTwoBodyStepPhaseSpace(parent_mass, daughter_mass, other_mass)) {
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::KinematicallyForbidden,
            "Directed angular-sector step has no kinematically allowed "
            "phase space");
    }

    double p_rest = TwoBodyRestMomentum(parent_mass, daughter_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, daughter_mass, other_mass);
    double p_parent = std::sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);

    DirectedGeometry geo = ClassifyDirectedRegime(
        parent_E, parent_px, parent_py, parent_pz,
        parent_mass, daughter_mass, other_mass, decay_pos, target);

    if (!SectorActive(geo.regime, geo.inside_geometry)) {
        return IsotropicBoostStep(parent_E, parent_px, parent_py, parent_pz,
                                  parent_mass, daughter_mass, other_mass, random);
    }

    // Active: sample within the (u, phi) bin of the bounding cone.
    siren::math::Vector3D axis = NormalizedSectorAxis(geo.to_center);
    double cos_bound = std::cos(geo.theta_bound);
    double c_hi = 1.0 - bin.u_lo * (1.0 - cos_bound);
    double c_lo = 1.0 - bin.u_hi * (1.0 - cos_bound);
    double omega_bin = (c_hi - c_lo) * (bin.phi_hi - bin.phi_lo);
    double g_angular = (omega_bin > 1e-300) ? 1.0 / omega_bin : 1.0 / kFourPi;

    double ct = c_lo + (c_hi - c_lo) * random->Uniform(0, 1);
    double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
    double ph = bin.phi_lo + (bin.phi_hi - bin.phi_lo) * random->Uniform(0, 1);
    siren::math::Vector3D perp1, perp2;
    SectorPerpFrame(axis, perp1, perp2);
    siren::math::Vector3D lab_dir =
        axis * ct + perp1 * (st * std::cos(ph)) + perp2 * (st * std::sin(ph));
    lab_dir.normalize();

    DirectedStepResult result;
    if (geo.regime == DirectedRegime::Rest) {
        result.lab_dir = lab_dir;
        result.p_lab = p_rest;
        result.E_lab = std::sqrt(p_rest * p_rest + daughter_mass * daughter_mass);
        result.rest_density = g_angular;
        return result;
    }

    // BoundInKin: boost to lab along this direction.
    double beta = p_parent / parent_E;
    double gamma = parent_E / parent_mass;
    siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);
    double cos_lab = siren::math::scalar_product(lab_dir, parent_dir);
    auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_lab);
    int n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);
    if (n_valid == 0) {
        // An isotropic fallback would add support that DensityAngularSectorStep
        // does not model. Reject this injection attempt instead of returning a
        // sample with an incorrect proposal density.
        throw siren::utilities::InjectionFailure(
            "Directed angular-sector sampler found no valid kinematic branch");
    }
    int chosen = 0;
    if (n_valid == 2) {
        double w0 = solutions[0].jacobian;
        double w1 = solutions[1].jacobian;
        chosen = (random->Uniform(0, 1) * (w0 + w1) < w0) ? 0 : 1;
    } else {
        chosen = solutions[0].valid ? 0 : 1;
    }
    double J_chosen = solutions[chosen].jacobian;
    double J_total = 0.0;
    for (auto const & sol : solutions) {
        if (sol.valid) J_total += sol.jacobian;
    }
    result.lab_dir = lab_dir;
    result.p_lab = solutions[chosen].p_lab;
    result.E_lab = std::sqrt(result.p_lab * result.p_lab + daughter_mass * daughter_mass);
    result.rest_density = (J_total > 0.0)
        ? g_angular * J_chosen * J_chosen / J_total
        : g_angular;
    return result;
}

inline double DensityAngularSectorStep(
    double parent_E,
    double parent_px, double parent_py, double parent_pz,
    double parent_mass,
    double daughter_mass,
    double other_mass,
    double daughter_E,
    double daughter_px, double daughter_py, double daughter_pz,
    siren::math::Vector3D const & decay_pos,
    siren::geometry::Geometry const & target,
    AngularSectorBin const & bin)
{
    if (!HasTwoBodyStepPhaseSpace(parent_mass, daughter_mass, other_mass)) {
        return 0.0;
    }

    double p_rest = TwoBodyRestMomentum(parent_mass, daughter_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, daughter_mass, other_mass);
    double p_parent = std::sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);

    DirectedGeometry geo = ClassifyDirectedRegime(
        parent_E, parent_px, parent_py, parent_pz,
        parent_mass, daughter_mass, other_mass, decay_pos, target);

    if (!SectorActive(geo.regime, geo.inside_geometry)) {
        return 1.0 / kFourPi;
    }

    double p_A = std::sqrt(daughter_px*daughter_px + daughter_py*daughter_py + daughter_pz*daughter_pz);
    if (p_A < 1e-15) return 0.0;
    siren::math::Vector3D lab_dir(daughter_px/p_A, daughter_py/p_A, daughter_pz/p_A);

    siren::math::Vector3D axis = NormalizedSectorAxis(geo.to_center);
    double cos_to = siren::math::scalar_product(lab_dir, axis);
    double cos_bound = std::cos(geo.theta_bound);
    double denom = 1.0 - cos_bound;
    if (denom < 1e-300) return 0.0;
    if (cos_to < cos_bound) return 0.0;             // outside the bounding cone
    double u = (1.0 - cos_to) / denom;
    if (u < bin.u_lo || u > bin.u_hi) return 0.0;   // outside this polar bin

    siren::math::Vector3D perp1, perp2;
    SectorPerpFrame(axis, perp1, perp2);
    siren::math::Vector3D perp = lab_dir - axis * cos_to;
    double phi = std::atan2(siren::math::scalar_product(perp, perp2),
                            siren::math::scalar_product(perp, perp1));
    if (!PhiInAngularSector(phi, bin)) return 0.0;  // outside azimuth bin

    double c_hi = 1.0 - bin.u_lo * denom;
    double c_lo = 1.0 - bin.u_hi * denom;
    double omega_bin = (c_hi - c_lo) * (bin.phi_hi - bin.phi_lo);
    double g_angular = (omega_bin > 1e-300) ? 1.0 / omega_bin : 1.0 / kFourPi;

    if (geo.regime == DirectedRegime::Rest) {
        return g_angular;
    }

    // BoundInKin: boost Jacobian (mirrors DensityDirectedStep).
    siren::math::Vector3D parent_dir(parent_px/p_parent, parent_py/p_parent, parent_pz/p_parent);
    double beta = p_parent / parent_E;
    double gamma = parent_E / parent_mass;
    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);
    double cos_critical = CriticalCosTheta(beta, gamma, p_rest, E_rest, daughter_mass);
    if (cos_theta_lab < cos_critical) return 0.0;

    auto solutions = SolveLabAngle(beta, gamma, p_rest, E_rest, daughter_mass, cos_theta_lab);
    double p_par_lab = p_A * cos_theta_lab;
    double p_par_rest = gamma * (p_par_lab - beta * daughter_E);
    double cos_rest_actual = (p_rest > 0.0) ? p_par_rest / p_rest : 0.0;

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

    return g_angular * J_chosen * J_chosen / J_total;
}

} // namespace detail
} // namespace injection
} // namespace siren

#endif // SIREN_DetectorDirectedChannelUtils_H
