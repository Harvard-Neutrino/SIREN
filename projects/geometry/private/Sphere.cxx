#include "SIREN/geometry/Sphere.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <math.h>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"
#include "GeometryMacros.h"
#include "PhiUtils.h"
#include "RayIntervals.h"

namespace siren {
namespace geometry {

namespace {
using phi_utils::InitialPhiState;
using phi_utils::ZAxisWedgeEntering;
} // anonymous namespace

Sphere::Sphere() : Geometry("Sphere"), radius_(0), inner_radius_(0), start_phi_(0), delta_phi_(2.0 * M_PI), start_theta_(0), delta_theta_(M_PI), has_phi_cut_(false), has_theta_cut_(false) { RecomputeWorldAABB(); }
Sphere::Sphere(double radius, double inner_radius) : Geometry("Sphere"), radius_(radius), inner_radius_(inner_radius), start_phi_(0), delta_phi_(2.0 * M_PI), start_theta_(0), delta_theta_(M_PI), has_phi_cut_(false), has_theta_cut_(false) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    RecomputeWorldAABB();
}
Sphere::Sphere(Placement const & p) : Geometry("Sphere", p), radius_(0), inner_radius_(0), start_phi_(0), delta_phi_(2.0 * M_PI), start_theta_(0), delta_theta_(M_PI), has_phi_cut_(false), has_theta_cut_(false) { RecomputeWorldAABB(); }
Sphere::Sphere(Placement const & p, double radius, double inner_radius) : Geometry("Sphere", p), radius_(radius), inner_radius_(inner_radius), start_phi_(0), delta_phi_(2.0 * M_PI), start_theta_(0), delta_theta_(M_PI), has_phi_cut_(false), has_theta_cut_(false) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    RecomputeWorldAABB();
}
Sphere::Sphere(double radius, double inner_radius, double start_phi, double delta_phi, double start_theta, double delta_theta)
    : Geometry("Sphere"), radius_(radius), inner_radius_(inner_radius), start_phi_(start_phi), delta_phi_(delta_phi), start_theta_(start_theta), delta_theta_(delta_theta) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    if(start_theta_ < -1e-9 || start_theta_ > M_PI + 1e-9) throw std::invalid_argument("Sphere start_theta must be in [0, pi]!");
    if(delta_theta_ <= 0 || start_theta_ + delta_theta_ > M_PI + 1e-9) throw std::invalid_argument("Sphere start_theta + delta_theta must be in (0, pi]!");
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
    RecomputeWorldAABB();
}
Sphere::Sphere(Placement const & p, double radius, double inner_radius, double start_phi, double delta_phi, double start_theta, double delta_theta)
    : Geometry("Sphere", p), radius_(radius), inner_radius_(inner_radius), start_phi_(start_phi), delta_phi_(delta_phi), start_theta_(start_theta), delta_theta_(delta_theta) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    if(start_theta_ < -1e-9 || start_theta_ > M_PI + 1e-9) throw std::invalid_argument("Sphere start_theta must be in [0, pi]!");
    if(delta_theta_ <= 0 || start_theta_ + delta_theta_ > M_PI + 1e-9) throw std::invalid_argument("Sphere start_theta + delta_theta must be in (0, pi]!");
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
    RecomputeWorldAABB();
}
Sphere::Sphere(const Sphere& o) : Geometry(o), radius_(o.radius_), inner_radius_(o.inner_radius_), start_phi_(o.start_phi_), delta_phi_(o.delta_phi_), start_theta_(o.start_theta_), delta_theta_(o.delta_theta_), has_phi_cut_(o.has_phi_cut_), has_theta_cut_(o.has_theta_cut_) { RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(Sphere, radius_, inner_radius_, start_phi_, delta_phi_, start_theta_, delta_theta_, has_phi_cut_, has_theta_cut_)
SIREN_GEOMETRY_ASSIGN(Sphere)
SIREN_GEOMETRY_EQUAL(Sphere, radius_, inner_radius_, start_phi_, delta_phi_, start_theta_, delta_theta_)
SIREN_GEOMETRY_LESS(Sphere, radius_, inner_radius_, start_phi_, delta_phi_, start_theta_, delta_theta_)

// ------------------------------------------------------------------------- //
void Sphere::print(std::ostream& os) const {
    os << "Radius: " << radius_ << "\tInner radius: " << inner_radius_;
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    if(has_theta_cut_) os << "\tStartTheta: " << start_theta_ << "\tDeltaTheta: " << delta_theta_;
    os << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Sphere::ComputeIntersections(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const {

    double px = position.GetX(), py = position.GetY(), pz = position.GetZ();
    double dx = direction.GetX(), dy = direction.GetY(), dz = direction.GetZ();
    double pd = px*dx + py*dy + pz*dz;
    double dd = dx*dx + dy*dy + dz*dz;
    if(dd == 0) return {};

    // Interval (slab) method. The solid is
    //     (ball MINUS inner ball) INTERSECT (theta band) INTERSECT (phi wedge).
    // The first two factors are solids of revolution about z, so their
    // in-solid ray-parameter sets are short lists of disjoint intervals
    // computed with exact set algebra; emitting each connected piece of
    // their intersection as an (enter, exit) pair makes the surface hit
    // list parity-consistent by construction. Only the phi wedge remains
    // a CSG state walk, consuming the interval endpoints as surface hits
    // (same pattern as Cylinder and Polycone). The previous per-root
    // theta-cone tests filtered each root through tolerance windows
    // (nappe sign, apex proximity) with sign-derived entering flags,
    // which could emit inconsistent band sequences at the cone/sphere
    // seam circles and silently corrupt the state walk.

    using ray_intervals::IntervalSet;

    // ---- Radial shell: {|x| <= radius} minus the inner ball ----
    // det = dd*r^2 - |p x d|^2 (squared perpendicular distance from the
    // origin to the ray, times dd). Computing |p x d|^2 directly avoids
    // catastrophic cancellation when the ray origin is far from the sphere.
    double cx = py * dz - pz * dy;
    double cy = pz * dx - px * dz;
    double cz = px * dy - py * dx;
    double perp2 = cx*cx + cy*cy + cz*cz;
    double inv_dd = 1.0 / dd;
    double det_o = dd * radius_ * radius_ - perp2;
    if(det_o <= 0) return {}; // miss, or tangent line (zero measure)
    double sq_o = std::sqrt(det_o);
    IntervalSet solid;
    solid.Add((-pd - sq_o) * inv_dd, (-pd + sq_o) * inv_dd);
    if(inner_radius_ > 0) {
        double det_i = dd * inner_radius_ * inner_radius_ - perp2;
        if(det_i > 0) {
            double sq_i = std::sqrt(det_i);
            IntervalSet bore;
            bore.Add((-pd - sq_i) * inv_dd, (-pd + sq_i) * inv_dd);
            solid = ray_intervals::Subtract(solid, bore);
        }
    }

    // ---- Theta band: {start_theta <= theta(x(t)) <= start_theta + delta} ----
    // theta(x) <= theta0 is z >= cos(theta0)*|x|, split into the exact
    // z-halfspace and the sign set of the quadratic z^2 - c^2*|x|^2 (the
    // squared form covers both cone nappes; the halfspace picks the right
    // one). The {z >= 0} and {z <= 0} sets share their boundary parameter
    // bit-exactly, so unions across the equator plane coalesce exactly.
    if(has_theta_cut_ && solid.n > 0) {
        IntervalSet const h_pos = ray_intervals::HalfLineGEQ(pz, dz);   // {z(t) >= 0}
        IntervalSet const h_neg = ray_intervals::HalfLineGEQ(-pz, -dz); // {z(t) <= 0}

        // Cone quadrics are evaluated about the ray's closest approach to
        // the apex: there the squared apex distance is exactly perp2/dd
        // (free of cancellation) and the linear coefficient collapses,
        // conditioning far-origin rays.
        double t_ca = -pd * inv_dd;
        double qz_ca = pz + t_ca * dz;
        double qd_ca = pd + t_ca * dd; // ~0; kept so the shifted polynomial is exact
        double qq_ca = perp2 * inv_dd;
        auto cone_quadric = [&](double c2, bool geq) {
            double A = dz*dz - c2 * dd;
            double B_half = qz_ca * dz - c2 * qd_ca;
            double C = qz_ca * qz_ca - c2 * qq_ca;
            IntervalSet s = geq ? ray_intervals::QuadraticGEQ(A, B_half, C)
                                : ray_intervals::QuadraticLEQ(A, B_half, C);
            s.Shift(t_ca);
            return s;
        };

        // {theta(x(t)) <= theta0}, i.e. {z >= cos(theta0)*|x|}
        auto theta_below = [&](double theta0) {
            if(std::fabs(theta0 - M_PI / 2.0) < 1e-12) return h_pos;
            double c = std::cos(theta0);
            if(c > 0) return ray_intervals::Intersect(h_pos, cone_quadric(c*c, true));
            return ray_intervals::Union(h_pos,
                ray_intervals::Intersect(h_neg, cone_quadric(c*c, false)));
        };
        // {theta(x(t)) >= theta0}, i.e. {z <= cos(theta0)*|x|}
        auto theta_above = [&](double theta0) {
            if(std::fabs(theta0 - M_PI / 2.0) < 1e-12) return h_neg;
            double c = std::cos(theta0);
            if(c < 0) return ray_intervals::Intersect(h_neg, cone_quadric(c*c, true));
            return ray_intervals::Union(h_neg,
                ray_intervals::Intersect(h_pos, cone_quadric(c*c, false)));
        };

        double theta_max = start_theta_ + delta_theta_;
        if(theta_max < M_PI - 1e-12)
            solid = ray_intervals::Intersect(solid, theta_below(theta_max));
        if(start_theta_ > 1e-12 && solid.n > 0)
            solid = ray_intervals::Intersect(solid, theta_above(start_theta_));
    }
    if(solid.n == 0) return {};

    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    TaggedHit all_hits[2 * IntervalSet::CAP + 2];
    int n_all = 0;

    // Emit interval endpoints, keeping the on-border convention that a
    // boundary within GEOMETRY_PRECISION ahead of the ray origin counts as
    // distance zero; a piece the snap collapses is dropped and pieces it
    // makes touch are merged, so no zero-measure pair is ever emitted.
    {
        double merged_lo = 0, merged_hi = 0;
        bool have_piece = false;
        auto emit_piece = [&](double lo, double hi) {
            all_hits[n_all++] = {lo, siren::math::Vector3D(px + lo*dx, py + lo*dy, pz + lo*dz), true, 0};
            all_hits[n_all++] = {hi, siren::math::Vector3D(px + hi*dx, py + hi*dy, pz + hi*dz), false, 0};
        };
        for(int i = 0; i < solid.n; ++i) {
            double lo = solid.lo[i];
            double hi = solid.hi[i];
            if(lo > 0 && lo < GEOMETRY_PRECISION) lo = 0;
            if(hi > 0 && hi < GEOMETRY_PRECISION) hi = 0;
            if(!(lo < hi)) continue;
            if(have_piece && lo <= merged_hi) {
                if(hi > merged_hi) merged_hi = hi;
                continue;
            }
            if(have_piece) emit_piece(merged_lo, merged_hi);
            merged_lo = lo;
            merged_hi = hi;
            have_piece = true;
        }
        if(have_piece) emit_piece(merged_lo, merged_hi);
    }
    if(n_all == 0) return {};

    if(!has_phi_cut_) {
        std::vector<Intersection> result;
        result.reserve(n_all);
        for(int i = 0; i < n_all; ++i) {
            Intersection isect;
            isect.distance = all_hits[i].distance;
            isect.hierarchy = 0;
            isect.entering = all_hits[i].entering;
            isect.position = all_hits[i].position;
            result.push_back(isect);
        }
        return result;
    }

    // Phi cut: compute infinite wedge intersections (two half-planes from
    // the z-axis). See Polycone.cxx for method description; same pattern in
    // all phi-cut shapes.
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
            entering = ZAxisWedgeEntering(dx, dy, start_phi_, delta_phi_);
        }
        all_hits[n_all] = {t, siren::math::Vector3D(hx, hy, hz), entering, 1};
        n_all++;
    }

    if(n_all == 0) return {};

    std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    bool in_surface = false;
    bool in_wedge = InitialPhiState(px, py, dx, dy, start_phi_, delta_phi_);

    bool was_inside = in_surface && in_wedge;

    std::vector<Intersection> result;
    for(int i = 0; i < n_all; ++i) {
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
AABB Sphere::GetBoundingBox() const {
    return AABB(
        math::Vector3D(-radius_, -radius_, -radius_),
        math::Vector3D( radius_,  radius_,  radius_)
    );
}

} // namespace geometry
} // namespace siren

CEREAL_REGISTER_DYNAMIC_INIT(siren_Sphere);

