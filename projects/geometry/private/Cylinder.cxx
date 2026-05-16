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

namespace siren {
namespace geometry {

namespace {

static const double TWO_PI = 2.0 * M_PI;

// Normalize angle to [0, 2*pi)
double NormalizePhi(double phi) {
    phi = std::fmod(phi, TWO_PI);
    if(phi < 0) phi += TWO_PI;
    return phi;
}

// Check if the azimuthal angle of point (x,y) falls within
// [start_phi, start_phi + delta_phi]. Handles wraparound.
bool PhiInRange(double x, double y, double start_phi, double delta_phi) {
    double phi = NormalizePhi(std::atan2(y, x));
    double sp = NormalizePhi(start_phi);
    double ep = sp + delta_phi;
    if(ep <= TWO_PI + 1e-9) {
        return phi >= sp - 1e-9 && phi <= ep + 1e-9;
    } else {
        // Wraps around 2*pi: phi >= sp OR phi <= (ep - 2*pi)
        return phi >= sp - 1e-9 || phi <= NormalizePhi(ep) + 1e-9;
    }
}

} // anonymous namespace

Cylinder::Cylinder()
    : Geometry((std::string)("Cylinder"))
    , radius_(0.0)
    , inner_radius_(0.0)
    , z_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    // Do nothing here
}

Cylinder::Cylinder(double radius, double inner_radius, double z)
    : Geometry((std::string)("Cylinder"))
    , radius_(radius)
    , inner_radius_(inner_radius)
    , z_(z)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(inner_radius_ > radius_) {
        std::swap(inner_radius_, radius_);
    }
    if(inner_radius_ == radius_) {
        //log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
    }
}

Cylinder::Cylinder(Placement const & placement)
    : Geometry((std::string)("Cylinder"), placement)
    , radius_(0.0)
    , inner_radius_(0.0)
    , z_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    // Do nothing here
}

Cylinder::Cylinder(Placement const & placement, double radius, double inner_radius, double z)
    : Geometry((std::string)("Cylinder"), placement)
    , radius_(radius)
    , inner_radius_(inner_radius)
    , z_(z)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(inner_radius_ > radius_) {
        std::swap(inner_radius_, radius_);
    }
    if(inner_radius_ == radius_) {
        //log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
    }
}

Cylinder::Cylinder(double radius, double inner_radius, double z, double start_phi, double delta_phi)
    : Geometry((std::string)("Cylinder"))
    , radius_(radius)
    , inner_radius_(inner_radius)
    , z_(z)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(inner_radius_ > radius_) {
        std::swap(inner_radius_, radius_);
    }
    if(inner_radius_ == radius_) {
        //log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Cylinder delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Cylinder::Cylinder(Placement const & placement, double radius, double inner_radius, double z,
                   double start_phi, double delta_phi)
    : Geometry((std::string)("Cylinder"), placement)
    , radius_(radius)
    , inner_radius_(inner_radius)
    , z_(z)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(inner_radius_ > radius_) {
        std::swap(inner_radius_, radius_);
    }
    if(inner_radius_ == radius_) {
        //log_error("Warning: Inner radius %f == radius %f (Volume is 0)", inner_radius_, radius_);
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Cylinder delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Cylinder::Cylinder(const Cylinder& cylinder)
    : Geometry(cylinder)
    , radius_(cylinder.radius_)
    , inner_radius_(cylinder.inner_radius_)
    , z_(cylinder.z_)
    , start_phi_(cylinder.start_phi_)
    , delta_phi_(cylinder.delta_phi_)
    , has_phi_cut_(cylinder.has_phi_cut_) {
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Cylinder::swap(Geometry& geometry) {
    Cylinder* cylinder = dynamic_cast<Cylinder*>(&geometry);
    if(!cylinder) {
        //log_warn("Cannot swap Cylinder!");
        return;
    }

    Geometry::swap(*cylinder);

    std::swap(inner_radius_, cylinder->inner_radius_);
    std::swap(radius_, cylinder->radius_);
    std::swap(z_, cylinder->z_);
    std::swap(start_phi_, cylinder->start_phi_);
    std::swap(delta_phi_, cylinder->delta_phi_);
    std::swap(has_phi_cut_, cylinder->has_phi_cut_);
}

//------------------------------------------------------------------------- //
Cylinder& Cylinder::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);
        if(!cylinder) {
            //log_warn("Cannot assign Cylinder!");
            return *this;
        }
        Cylinder tmp(*cylinder);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Cylinder::equal(const Geometry& geometry) const {
    const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);

    if(!cylinder)
        return false;
    else if(inner_radius_ != cylinder->inner_radius_)
        return false;
    else if(radius_ != cylinder->radius_)
        return false;
    else if(z_ != cylinder->z_)
        return false;
    else if(start_phi_ != cylinder->start_phi_)
        return false;
    else if(delta_phi_ != cylinder->delta_phi_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Cylinder::less(const Geometry& geometry) const {
    const Cylinder* cylinder = dynamic_cast<const Cylinder*>(&geometry);
    if(!cylinder) return false;

    return
        std::tie(inner_radius_, radius_, z_, start_phi_, delta_phi_)
        <
        std::tie(cylinder->inner_radius_, cylinder->radius_, cylinder->z_, cylinder->start_phi_, cylinder->delta_phi_);
}

void Cylinder::print(std::ostream& os) const {
    os << "Radius: " << radius_ << "\tInnner radius: " << inner_radius_ << " Height: " << z_;
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    os << '\n';
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
        double A = px*px + py*py - r2;
        double det = B_half*B_half - C*A;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double inv_C = 1.0 / C;
        double t1 = (-B_half - sq) * inv_C;
        double t2 = (-B_half + sq) * inv_C;
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? t1 : t2;
            double iz = pz + t * dz;
            if(iz > -hz && iz < hz) {
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
        if(r2_hit <= r2_hi && r2_hit >= r2_lo) {
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
        // No phi cut: surface hits are the final result
        if(n_all == 0) return {};
        std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
            return a.distance < b.distance;
        });
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

    // Phi cut: compute infinite wedge intersections (two half-planes from z-axis)
    // No cross-section filtering -- the CSG walk handles clipping naturally.
    for(int face = 0; face < 2; ++face) {
        double alpha = start_phi_ + face * delta_phi_;
        double ca = std::cos(alpha), sa = std::sin(alpha);
        // Outward-pointing normal (away from phi range interior)
        double nx, ny;
        if(face == 0) { nx = sa; ny = -ca; }
        else { nx = -sa; ny = ca; }
        double n_dot_d = nx*dx + ny*dy;
        if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue;
        double n_dot_p = nx*px + ny*py;
        double t = -n_dot_p / n_dot_d;
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

        double hx = px + t*dx, hy = py + t*dy, h_z = pz + t*dz;
        // Must be on the correct half-plane (outward from z-axis)
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;

        bool entering = (n_dot_d < 0);
        all_hits[n_all] = {t, siren::math::Vector3D(hx, hy, h_z), entering, 1};
        n_all++;
    }

    if(n_all == 0) return {};

    // Sort all hits by distance
    std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // CSG intersection walk: the phi-cut solid is (full cylinder) AND (phi wedge).
    // Walk through sorted hits, tracking in_surface and in_wedge states.
    // Emit an intersection whenever the combined state changes.
    bool in_surface = false;
    bool in_wedge = false;

    // Determine initial in_wedge state
    bool has_wedge_hit = false;
    for(int i = 0; i < n_all; ++i) {
        if(all_hits[i].source == 1) {
            // First wedge hit: if entering, we started outside; if exiting, started inside
            in_wedge = !all_hits[i].entering;
            has_wedge_hit = true;
            break;
        }
    }
    if(!has_wedge_hit) {
        // No wedge crossings: ray is entirely in or entirely out of wedge.
        // Check the ray origin (or any point along the ray).
        in_wedge = PhiInRange(px, py, start_phi_, delta_phi_);
    }

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
