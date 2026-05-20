#include "SIREN/geometry/Cylinder.h"

#include <cmath>
#include <tuple>
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

    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    TaggedHit all_hits[12]; // max: 2 barrel + 2 caps + 2 inner barrel + 2 inner caps + 2 wedge + spare
    int n_all = 0;

    auto add_surface_hit = [&](double t, double ix, double iy, double iz, bool entering) {
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        all_hits[n_all] = {t, siren::math::Vector3D(ix, iy, iz), entering, 0};
        n_all++;
    };

    // Test barrel surface for a given radius squared; invert_entering flips
    // the radial entering logic for inner surfaces
    auto test_barrel = [&](double r2, bool invert_entering) {
        if(dx == 0 && dy == 0) return;
        double C = dx*dx + dy*dy;
        double B_half = px*dx + py*dy;
        // det = C*r2 - (px*dy - py*dx)^2
        // avoids catastrophic cancellation at large distances
        double cross_z = px * dy - py * dx;
        double det = C * r2 - cross_z * cross_z;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double inv_C = 1.0 / C;
        double t1 = (-B_half - sq) * inv_C;
        double t2 = (-B_half + sq) * inv_C;
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? t1 : t2;
            double iz = pz + t * dz;
            if(iz > -hz - GEOMETRY_PRECISION && iz < hz + GEOMETRY_PRECISION) {
                double ix = px + t * dx;
                double iy = py + t * dy;
                bool radial_entering = (ix * dx + iy * dy) < 0;
                add_surface_hit(t, ix, iy, iz, invert_entering ? !radial_entering : radial_entering);
            }
        }
    };

    // Test endcap at z_cap for a given annular ring [r2_lo, r2_hi]
    auto test_cap = [&](double z_cap, double r2_lo, double r2_hi, bool enter) {
        if(dz == 0) return;
        double t = (z_cap - pz) / dz;
        double ix = px + t * dx;
        double iy = py + t * dy;
        double r2_hit = ix*ix + iy*iy;
        if(r2_hit <= r2_hi + GEOMETRY_PRECISION && r2_hit >= r2_lo - GEOMETRY_PRECISION) {
            add_surface_hit(t, ix, iy, pz + t * dz, enter);
        }
    };

    // Outer barrel
    test_barrel(r2_outer, false);
    // Endcaps (annular disk between inner and outer radius)
    test_cap( hz, r2_inner, r2_outer, dz < 0);
    test_cap(-hz, r2_inner, r2_outer, dz > 0);
    // Inner barrel (if hollow)
    if(inner_radius_ > 0) {
        test_barrel(r2_inner, true);
    }

    if(!has_phi_cut_) {
        if(n_all == 0) return {};
        std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
            return a.distance < b.distance;
        });
        std::vector<Intersection> result;
        result.reserve(n_all);
        for(int i = 0; i < n_all; ++i) {
            if(!result.empty() && std::fabs(all_hits[i].distance - result.back().distance) < GEOMETRY_PRECISION) {
                continue;
            }
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
