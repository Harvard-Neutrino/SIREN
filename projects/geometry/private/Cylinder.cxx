#include "SIREN/geometry/Cylinder.h"

#include <cmath>
#include <tuple>
#include <limits>
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

namespace siren {
namespace geometry {

namespace {
using phi_utils::NormalizePhi;
using phi_utils::PhiInRange;
using phi_utils::InitialPhiState;
using phi_utils::ZAxisWedgeEntering;
using phi_utils::TWO_PI;
} // anonymous namespace

Cylinder::Cylinder() : Geometry("Cylinder"), radius_(0), inner_radius_(0), z_(0), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Cylinder::Cylinder(double radius, double inner_radius, double z) : Geometry("Cylinder"), radius_(radius), inner_radius_(inner_radius), z_(z), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    RecomputeWorldAABB();
}
Cylinder::Cylinder(Placement const & p) : Geometry("Cylinder", p), radius_(0), inner_radius_(0), z_(0), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Cylinder::Cylinder(Placement const & p, double radius, double inner_radius, double z) : Geometry("Cylinder", p), radius_(radius), inner_radius_(inner_radius), z_(z), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    RecomputeWorldAABB();
}
Cylinder::Cylinder(double radius, double inner_radius, double z, double start_phi, double delta_phi) : Geometry("Cylinder"), radius_(radius), inner_radius_(inner_radius), z_(z), start_phi_(start_phi), delta_phi_(delta_phi) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    RecomputeWorldAABB();
}
Cylinder::Cylinder(Placement const & p, double radius, double inner_radius, double z, double start_phi, double delta_phi) : Geometry("Cylinder", p), radius_(radius), inner_radius_(inner_radius), z_(z), start_phi_(start_phi), delta_phi_(delta_phi) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    RecomputeWorldAABB();
}
Cylinder::Cylinder(const Cylinder& o) : Geometry(o), radius_(o.radius_), inner_radius_(o.inner_radius_), z_(o.z_), start_phi_(o.start_phi_), delta_phi_(o.delta_phi_), has_phi_cut_(o.has_phi_cut_) { RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(Cylinder, radius_, inner_radius_, z_, start_phi_, delta_phi_, has_phi_cut_)
SIREN_GEOMETRY_ASSIGN(Cylinder)
SIREN_GEOMETRY_EQUAL(Cylinder, radius_, inner_radius_, z_, start_phi_, delta_phi_, has_phi_cut_)
SIREN_GEOMETRY_LESS(Cylinder, radius_, inner_radius_, z_, start_phi_, delta_phi_, has_phi_cut_)

void Cylinder::print(std::ostream& os) const {
    os << "Cylinder(" << radius_ << ", " << inner_radius_ << ", " << z_;
    if(has_phi_cut_) os << ", " << start_phi_ << ", " << delta_phi_;
    os << ")\n";
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Cylinder::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();
    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();

    double hz = 0.5 * z_;
    double r2_outer = radius_ * radius_;
    double r2_inner = inner_radius_ * inner_radius_;

    // Interval (slab) method. The closed solid is
    //     (infinite outer cylinder INTERSECT z-slab) MINUS infinite inner cylinder,
    // so the ray-parameter set inside it is one or two disjoint intervals.
    // Emitting each nonempty interval as an (enter, exit) pair makes the
    // hit list parity-consistent by construction. Testing each surface
    // independently with tolerance windows (the previous approach) let a
    // crossing near the cap/barrel corner pass one window and fail the
    // other, emitting an unpaired hit that broke every consumer relying
    // on enter/exit alternation (DetectorModel::SectorLoop in particular).

    double const inf = std::numeric_limits<double>::infinity();

    // Ray-parameter interval inside the z-slab |z| <= hz.
    double tz_lo, tz_hi;
    if(dz != 0) {
        double inv = 1.0 / dz;
        tz_lo = (-hz - pz) * inv;
        tz_hi = ( hz - pz) * inv;
        if(tz_lo > tz_hi) std::swap(tz_lo, tz_hi);
    } else {
        if(pz < -hz || pz > hz) return {};
        tz_lo = -inf;
        tz_hi = inf;
    }

    // Ray-parameter interval inside the infinite outer cylinder.
    // det = C*r2 - (px*dy - py*dx)^2 avoids catastrophic cancellation at
    // large distances.
    double C = dx*dx + dy*dy;
    double B_half = px*dx + py*dy;
    double cross_z = px * dy - py * dx;
    double to_lo, to_hi;
    if(C != 0) {
        double det = C * r2_outer - cross_z * cross_z;
        if(det <= 0) return {}; // miss, or tangent line (zero measure)
        double sq = std::sqrt(det);
        double inv_C = 1.0 / C;
        to_lo = (-B_half - sq) * inv_C;
        to_hi = (-B_half + sq) * inv_C;
    } else {
        if(px*px + py*py > r2_outer) return {};
        to_lo = -inf;
        to_hi = inf;
    }

    double a = std::max(tz_lo, to_lo);
    double b = std::min(tz_hi, to_hi);
    if(!(a < b)) return {};

    // Subtract the bore of a hollow cylinder: up to two pieces remain.
    double pieces[2][2];
    int n_pieces = 0;
    bool have_bore = false;
    double ti_lo = 0, ti_hi = 0;
    if(inner_radius_ > 0) {
        if(C != 0) {
            double det_i = C * r2_inner - cross_z * cross_z;
            if(det_i > 0) {
                double sq_i = std::sqrt(det_i);
                double inv_C = 1.0 / C;
                ti_lo = (-B_half - sq_i) * inv_C;
                ti_hi = (-B_half + sq_i) * inv_C;
                have_bore = true;
            }
        } else if(px*px + py*py < r2_inner) {
            return {}; // axis-parallel ray inside the bore
        }
    }
    if(have_bore) {
        double lo1 = a, hi1 = std::min(b, ti_lo);
        double lo2 = std::max(a, ti_hi), hi2 = b;
        if(lo1 < hi1) { pieces[n_pieces][0] = lo1; pieces[n_pieces][1] = hi1; ++n_pieces; }
        if(lo2 < hi2) { pieces[n_pieces][0] = lo2; pieces[n_pieces][1] = hi2; ++n_pieces; }
    } else {
        pieces[0][0] = a;
        pieces[0][1] = b;
        n_pieces = 1;
    }
    if(n_pieces == 0) return {};

    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    TaggedHit all_hits[8]; // max: 2 pieces x 2 endpoints + 2 wedge planes
    int n_all = 0;

    // Emit interval endpoints, keeping the on-border convention that a
    // boundary within GEOMETRY_PRECISION ahead of the ray origin counts as
    // distance zero; a pair the snap collapses to zero measure is dropped.
    for(int i = 0; i < n_pieces; ++i) {
        double lo = pieces[i][0];
        double hi = pieces[i][1];
        if(lo > 0 && lo < GEOMETRY_PRECISION) lo = 0;
        if(hi > 0 && hi < GEOMETRY_PRECISION) hi = 0;
        if(!(lo < hi)) continue;
        all_hits[n_all++] = {lo, siren::math::Vector3D(px + lo*dx, py + lo*dy, pz + lo*dz), true, 0};
        all_hits[n_all++] = {hi, siren::math::Vector3D(px + hi*dx, py + hi*dy, pz + hi*dz), false, 0};
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

    // Phi cut: compute infinite wedge intersections (two half-planes from z-axis).
    // See Polycone.cxx for method description; same pattern in all phi-cut shapes.
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

        double hx = px + t*dx, hy = py + t*dy, h_z = pz + t*dz;
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;

        bool entering = (n_dot_d < 0);
        if(hx*hx + hy*hy < GEOMETRY_PRECISION * 1e3 * GEOMETRY_PRECISION * 1e3) {
            if(z_axis_hit_emitted) continue;
            z_axis_hit_emitted = true;
            entering = ZAxisWedgeEntering(dx, dy, start_phi_, delta_phi_);
        }
        all_hits[n_all] = {t, siren::math::Vector3D(hx, hy, h_z), entering, 1};
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
AABB Cylinder::GetBoundingBox() const {
    double hz = z_ * 0.5;
    return AABB(
        math::Vector3D(-radius_, -radius_, -hz),
        math::Vector3D( radius_,  radius_,  hz)
    );
}

} // namespace geometry
} // namespace siren
