#include "SIREN/geometry/Polycone.h"

#include <cmath>
#include <tuple>
#include <limits>
#include <string>
#include <vector>
#include <ostream>
#include <cassert>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"
#include "GeometryMacros.h"
#include "PhiUtils.h"
#include "RayIntervals.h"

namespace siren {
namespace geometry {

void Polycone::validate() {
    if(z_planes_.size() < 2) {
        throw std::runtime_error("Polycone requires at least 2 z-planes!");
    }
    if(z_planes_.size() != rmin_.size() || z_planes_.size() != rmax_.size()) {
        throw std::runtime_error("Polycone z_planes, rmin, and rmax vectors must have the same size!");
    }
    // Accept ascending or descending z-order (including duplicate z-values).
    // If descending, reverse all three vectors so the geometry math sees
    // ascending order. Non-monotonic sequences are rejected.
    // This must happen before the rmin/rmax check since the reversal
    // reorders the radii arrays.
    bool ascending = true, descending = true;
    for(size_t i = 1; i < z_planes_.size(); ++i) {
        if(z_planes_[i] < z_planes_[i - 1]) ascending = false;
        if(z_planes_[i] > z_planes_[i - 1]) descending = false;
    }
    if(!ascending && !descending) {
        throw std::runtime_error("Polycone z-planes must be monotonic (ascending or descending)!");
    }
    if(descending && !ascending) {
        std::reverse(z_planes_.begin(), z_planes_.end());
        std::reverse(rmin_.begin(), rmin_.end());
        std::reverse(rmax_.begin(), rmax_.end());
    }
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        if(rmin_[i] < 0 || rmax_[i] < 0) {
            throw std::runtime_error("Polycone radii must be non-negative!");
        }
        if(rmin_[i] > rmax_[i]) {
            throw std::runtime_error("Polycone inner radius must not exceed outer radius at each z-plane!");
        }
    }
}

Polycone::Polycone() : Geometry("Polycone"), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Polycone::Polycone(std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax) : Geometry("Polycone"), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { validate(); RecomputeWorldAABB(); }
Polycone::Polycone(Placement const & p) : Geometry("Polycone", p), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Polycone::Polycone(Placement const & p, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax) : Geometry("Polycone", p), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { validate(); RecomputeWorldAABB(); }
Polycone::Polycone(std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax, double start_phi, double delta_phi) : Geometry("Polycone"), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(start_phi), delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    RecomputeWorldAABB();
}
Polycone::Polycone(Placement const & p, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax, double start_phi, double delta_phi) : Geometry("Polycone", p), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(start_phi), delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    RecomputeWorldAABB();
}
Polycone::Polycone(const Polycone& o) : Geometry(o), z_planes_(o.z_planes_), rmin_(o.rmin_), rmax_(o.rmax_), start_phi_(o.start_phi_), delta_phi_(o.delta_phi_), has_phi_cut_(o.has_phi_cut_) { RecomputeWorldAABB(); }

// Swap includes derived has_phi_cut_; equal/less use only primary fields
SIREN_GEOMETRY_SWAP(Polycone, z_planes_, rmin_, rmax_, start_phi_, delta_phi_, has_phi_cut_)
SIREN_GEOMETRY_ASSIGN(Polycone)
SIREN_GEOMETRY_EQUAL(Polycone, z_planes_, rmin_, rmax_, start_phi_, delta_phi_)
SIREN_GEOMETRY_LESS(Polycone, z_planes_, rmin_, rmax_, start_phi_, delta_phi_)

void Polycone::print(std::ostream& os) const {
    os << "Polycone(" << z_planes_.size() << " z-planes";
    if(has_phi_cut_) os << ", " << start_phi_ << ", " << delta_phi_;
    os << ")\n";
}

// ------------------------------------------------------------------------- //
// ComputeIntersections
//
// Interval (slab) method. The polycone is the union over sections of
//     (z-slab) INTERSECT (outer cone) MINUS (inner cone),
// one conical frustum section per pair of adjacent z-planes. Each
// section's in-solid ray-parameter set is a short list of disjoint
// intervals from exact quadratic sign conditions, and the joint plane
// ray parameters are computed once per plane index, so adjacent
// sections share their boundary parameter bit-exactly. The union over
// sections therefore merges seamlessly where material continues across
// a joint and splits exactly where it does not (step joints, end caps,
// grazes), and emitting each connected piece as an (enter, exit) pair
// makes the hit list parity-consistent by construction. End caps and
// internal step annuli need no dedicated surface tests: they are
// exactly the piece boundaries that land on joint plane parameters.
//
// Testing each surface independently with tolerance windows (the
// previous approach) let a crossing near a joint or corner circle pass
// one window and fail another, emitting unpaired hits that broke every
// consumer relying on enter/exit alternation
// (DetectorModel::SectorLoop in particular).
//
// All coordinates are in local frame. Distance can be negative
// (full-line intersections, same as Cone.cxx).
// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Polycone::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {

    if(z_planes_.size() < 2) {
        return {};
    }

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    size_t n = z_planes_.size();

    // Quick rejects against the z-range and the bounding cylinder keep
    // the common miss path allocation-free.
    if(dz == 0 && (pz < z_planes_.front() || pz >= z_planes_.back())) {
        return {};
    }
    double C_xy = dx * dx + dy * dy;
    double cross_z = px * dy - py * dx;
    double max_r = *std::max_element(rmax_.begin(), rmax_.end());
    if(C_xy != 0) {
        // det = C*r^2 - (px*dy - py*dx)^2 avoids catastrophic cancellation
        // at large distances.
        if(C_xy * max_r * max_r - cross_z * cross_z <= 0) {
            return {}; // miss, or tangent line (zero measure)
        }
    } else if(px * px + py * py > max_r * max_r) {
        return {};
    }

    // Origin shift: quadratics are evaluated about the closest approach to
    // the coordinate origin. This keeps their coefficients small for
    // far-field rays, avoiding catastrophic cancellation in the
    // discriminant B^2 - 4AC.
    double t_shift = -(px * dx + py * dy + pz * dz);
    double qx = px + t_shift * dx;
    double qy = py + t_shift * dy;
    double qz = pz + t_shift * dz;

    // Ray parameter of the crossing with the z-plane of index k. Adjacent
    // sections call this for the same k, giving bit-identical shared
    // boundary parameters: that is what lets the piece union below merge
    // exactly at continuous joints.
    double inv_dz = (dz != 0) ? 1.0 / dz : 0.0;
    auto plane_t = [&](size_t k) { return (z_planes_[k] - pz) * inv_dz; };

    using ray_intervals::IntervalSet;
    double const inf = std::numeric_limits<double>::infinity();

    std::vector<std::pair<double, double>> pieces;
    pieces.reserve(2 * (n - 1));

    for(size_t seg = 0; seg + 1 < n; ++seg) {
        double z_lo = z_planes_[seg];
        double z_hi = z_planes_[seg + 1];

        // Zero-height sections encode step joints; they carry no volume.
        if(z_hi <= z_lo) {
            continue;
        }

        // Ray-parameter interval inside this section's z-slab.
        IntervalSet slab;
        if(dz != 0) {
            double ta = plane_t(seg);
            double tb = plane_t(seg + 1);
            if(ta > tb) std::swap(ta, tb);
            slab.Add(ta, tb);
        } else {
            // Half-open [z_lo, z_hi): a horizontal ray at an internal
            // plane belongs to exactly one section.
            if(pz < z_lo || pz >= z_hi) {
                continue;
            }
            slab.Add(-inf, inf);
        }

        double dz_sec = z_hi - z_lo;

        // Outer cone r_outer(z) = a + b*z: the solid needs
        // x^2 + y^2 <= (a + b*z)^2. r_outer is non-negative throughout the
        // slab, so the mirror nappe of the squared form lies outside it.
        double b_outer = (rmax_[seg + 1] - rmax_[seg]) / dz_sec;
        double a_outer = rmax_[seg] - b_outer * z_lo;
        double r_q_outer = a_outer + b_outer * qz;
        IntervalSet sec = ray_intervals::QuadraticLEQ(
            C_xy - b_outer * b_outer * dz * dz,
            qx * dx + qy * dy - b_outer * dz * r_q_outer,
            qx * qx + qy * qy - r_q_outer * r_q_outer);
        sec.Shift(t_shift);
        sec = ray_intervals::Intersect(slab, sec);

        // Inner cone (hollow sections): subtract the bore.
        if(sec.n > 0 && (rmin_[seg] > 0 || rmin_[seg + 1] > 0)) {
            double b_inner = (rmin_[seg + 1] - rmin_[seg]) / dz_sec;
            double a_inner = rmin_[seg] - b_inner * z_lo;
            double r_q_inner = a_inner + b_inner * qz;
            IntervalSet bore = ray_intervals::QuadraticLEQ(
                C_xy - b_inner * b_inner * dz * dz,
                qx * dx + qy * dy - b_inner * dz * r_q_inner,
                qx * qx + qy * qy - r_q_inner * r_q_inner);
            bore.Shift(t_shift);
            sec = ray_intervals::Subtract(sec, bore);
        }

        for(int i = 0; i < sec.n; ++i) {
            pieces.emplace_back(sec.lo[i], sec.hi[i]);
        }
    }

    if(pieces.empty()) {
        return {};
    }

    // Sections were visited in ascending z; for dz < 0 the ray traverses
    // them in descending parameter order.
    std::sort(pieces.begin(), pieces.end());

    // Snap piece boundaries within GEOMETRY_PRECISION ahead of the ray
    // origin to distance zero (the on-border convention), drop pieces the
    // snap collapses, and merge pieces that touch or overlap: continuous
    // joints share their plane parameter bit-exactly and coalesce here,
    // so only real material boundaries survive as (enter, exit) pairs.
    std::vector<std::pair<double, double>> solid;
    solid.reserve(pieces.size());
    for(auto const & piece : pieces) {
        double lo = piece.first;
        double hi = piece.second;
        if(lo > 0 && lo < GEOMETRY_PRECISION) lo = 0;
        if(hi > 0 && hi < GEOMETRY_PRECISION) hi = 0;
        if(!(lo < hi)) continue;
        if(!solid.empty() && lo <= solid.back().second) {
            if(hi > solid.back().second) solid.back().second = hi;
            continue;
        }
        solid.emplace_back(lo, hi);
    }
    if(solid.empty()) {
        return {};
    }

    auto make_hit = [&](double t, bool entering) {
        Intersection isect;
        isect.distance = t;
        isect.hierarchy = 0;
        isect.entering = entering;
        isect.position = siren::math::Vector3D(px + t * dx, py + t * dy, pz + t * dz);
        return isect;
    };

    if(!has_phi_cut_) {
        // No phi cut: the interval endpoints are the final result
        std::vector<Intersection> result;
        result.reserve(2 * solid.size());
        for(auto const & piece : solid) {
            result.push_back(make_hit(piece.first, true));
            result.push_back(make_hit(piece.second, false));
        }
        return result;
    }

    // Phi cut: merge surface hits with infinite wedge hits and run CSG walk.
    // Method: intersect the ray with the full-rotation solid AND an infinite
    // wedge (two half-planes from the z-axis). A sorted walk over both hit
    // lists produces the CSG intersection. This pattern is duplicated across
    // Polycone, GenericPolycone, Sphere, Torus, Cylinder, Cone, CutTube.
    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    std::vector<TaggedHit> all_hits;
    all_hits.reserve(2 * solid.size() + 2);
    for(auto const & piece : solid) {
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? piece.first : piece.second;
            all_hits.push_back({t, siren::math::Vector3D(px + t * dx, py + t * dy, pz + t * dz), k == 0, 0});
        }
    }

    // Compute infinite wedge intersections (two half-planes from z-axis)
    bool z_axis_hit_emitted = false;
    for(int face = 0; face < 2; ++face) {
        double alpha = start_phi_ + face * delta_phi_;
        double ca = std::cos(alpha), sa = std::sin(alpha);
        double nx, ny;
        if(face == 0) { nx = sa; ny = -ca; }
        else { nx = -sa; ny = ca; }
        double n_dot_d = nx*dx + ny*dy;
        if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue;
        double n_dot_p = nx*px + ny*py;
        double t = -n_dot_p / n_dot_d;
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

        double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;

        bool entering = (n_dot_d < 0);
        if(hx*hx + hy*hy < GEOMETRY_PRECISION * 1e3 * GEOMETRY_PRECISION * 1e3) {
            if(z_axis_hit_emitted) continue;
            z_axis_hit_emitted = true;
            entering = phi_utils::ZAxisWedgeEntering(dx, dy, start_phi_, delta_phi_);
        }
        all_hits.push_back({t, siren::math::Vector3D(hx, hy, hz), entering, 1});
    }

    if(all_hits.empty()) return {};

    std::sort(all_hits.begin(), all_hits.end(), [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    bool in_surface = false;
    bool in_wedge = phi_utils::InitialPhiState(px, py, dx, dy, start_phi_, delta_phi_);

    bool was_inside = in_surface && in_wedge;

    std::vector<Intersection> result;
    for(size_t i = 0; i < all_hits.size(); ++i) {
        if(all_hits[i].source == 0) {
            in_surface = all_hits[i].entering;
        } else {
            in_wedge = all_hits[i].entering;
        }
        bool now_inside = in_surface && in_wedge;
        if(now_inside != was_inside) {
            Intersection isect;
            isect.distance = all_hits[i].distance;
            isect.hierarchy = 0;
            isect.entering = now_inside;
            isect.position = all_hits[i].position;
            result.push_back(isect);
        }
        was_inside = now_inside;
    }
    return result;
}

// ------------------------------------------------------------------------- //
AABB Polycone::GetBoundingBox() const {
    if(z_planes_.empty()) {
        return AABB(math::Vector3D(0, 0, 0), math::Vector3D(0, 0, 0));
    }
    double max_r = *std::max_element(rmax_.begin(), rmax_.end());
    double z_lo = z_planes_.front();
    double z_hi = z_planes_.back();

    if(!has_phi_cut_) {
        return AABB(
            math::Vector3D(-max_r, -max_r, z_lo),
            math::Vector3D( max_r,  max_r, z_hi)
        );
    }

    // Phi sector: compute bounding box from the two edge rays and any
    // cardinal directions (0, pi/2, pi, 3pi/2) that fall within the sector.
    double sp = phi_utils::NormalizePhi(start_phi_);
    double ep = sp + delta_phi_;

    double x_min = 0, x_max = 0, y_min = 0, y_max = 0;

    // Check the two edge directions
    double angles[2] = {sp, sp + delta_phi_};
    for(int i = 0; i < 2; ++i) {
        double a = angles[i];
        double cx = std::cos(a) * max_r;
        double cy = std::sin(a) * max_r;
        if(cx < x_min) x_min = cx;
        if(cx > x_max) x_max = cx;
        if(cy < y_min) y_min = cy;
        if(cy > y_max) y_max = cy;
    }

    // Check cardinal directions if they fall in the sector
    double cardinals[4] = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    for(int i = 0; i < 4; ++i) {
        double c = cardinals[i];
        // Check if cardinal is within [sp, ep] (with wraparound)
        double cn = c;
        if(cn < sp) cn += phi_utils::TWO_PI;
        if(cn <= ep + 1e-9) {
            double cx = std::cos(cardinals[i]) * max_r;
            double cy = std::sin(cardinals[i]) * max_r;
            if(cx < x_min) x_min = cx;
            if(cx > x_max) x_max = cx;
            if(cy < y_min) y_min = cy;
            if(cy > y_max) y_max = cy;
        }
    }

    return AABB(
        math::Vector3D(x_min, y_min, z_lo),
        math::Vector3D(x_max, y_max, z_hi)
    );
}

} // namespace geometry
} // namespace siren
