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

namespace siren {
namespace geometry {

namespace {
using phi_utils::NormalizePhi;
using phi_utils::PhiInRange;
using phi_utils::InitialPhiState;
using phi_utils::ZAxisWedgeEntering;
using phi_utils::TWO_PI;

// Check if the polar angle of point (x,y,z) falls within
// [start_theta, start_theta + delta_theta].
bool ThetaInRange(double x, double y, double z, double start_theta, double delta_theta) {
    double r = std::sqrt(x*x + y*y + z*z);
    if(r < 1e-15) return true; // origin is degenerate
    double cos_theta = z / r;
    if(cos_theta > 1.0) cos_theta = 1.0;
    if(cos_theta < -1.0) cos_theta = -1.0;
    double theta = std::acos(cos_theta);
    return theta >= start_theta - 1e-9 && theta <= start_theta + delta_theta + 1e-9;
}
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
    double pp = px*px + py*py + pz*pz;
    double pd = px*dx + py*dy + pz*dz;

    // The sphere with angular cuts is treated as:
    //   (spherical shell) AND (phi wedge) AND (theta band)
    // Each component produces independent enter/exit pairs. A CSG intersection
    // walk combines them without tolerance-dependent filtering or heuristics.
    // Source tags: 0 = shell, 1 = phi wedge, 2 = theta band
    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source;
    };

    TaggedHit all_hits[16]; // max: 4 shell + 4 phi planes + 8 theta cones
    int n_all = 0;

    auto add_hit = [&](double t, bool entering, int source) {
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        all_hits[n_all] = {t, siren::math::Vector3D(px + t*dx, py + t*dy, pz + t*dz), entering, source};
        n_all++;
    };

    // ---- Shell intersections (no angular filtering) ----
    // The discriminant is r^2 - |p x d|^2 (squared perpendicular distance
    // from origin to ray). Computing |p x d|^2 directly avoids catastrophic
    // cancellation that occurs when the ray origin is far from the sphere.
    double cx = py * dz - pz * dy;
    double cy = pz * dx - px * dz;
    double cz = px * dy - py * dx;
    double perp2 = cx*cx + cy*cy + cz*cz;
    auto solve_shell = [&](double r, bool is_outer) {
        double det = r * r - perp2;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double t1 = -pd - sq;
        double t2 = -pd + sq;
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? t1 : t2;
            bool entering = is_outer ? (k == 0) : (k == 1);
            add_hit(t, entering, 0);
        }
    };

    solve_shell(radius_, true);
    if(inner_radius_ > 0) {
        solve_shell(inner_radius_, false);
    }

    // ---- Infinite phi wedge (no radial/theta bounds) ----
    if(has_phi_cut_) {
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
            double hx = px + t*dx, hy = py + t*dy;
            if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;
            bool entering = (n_dot_d < 0);
            if(hx*hx + hy*hy < GEOMETRY_PRECISION * 1e3 * GEOMETRY_PRECISION * 1e3) {
                if(z_axis_hit_emitted) continue;
                z_axis_hit_emitted = true;
                entering = ZAxisWedgeEntering(dx, dy, start_phi_, delta_phi_);
            }
            add_hit(t, entering, 1);
        }
    }

    // ---- Infinite theta band boundaries (no radial/phi bounds) ----
    if(has_theta_cut_) {
        bool theta_apex_emitted = false;
        for(int face = 0; face < 2; ++face) {
            double theta0 = start_theta_ + face * delta_theta_;
            if(theta0 < 1e-12 || std::fabs(theta0 - M_PI) < 1e-12) continue;

            if(std::fabs(theta0 - M_PI / 2.0) < 1e-12) {
                if(std::fabs(dz) < GEOMETRY_PRECISION) continue;
                double t = -pz / dz;
                bool entering = (face == 0) ? (dz < 0) : (dz > 0);
                add_hit(t, entering, 2);
                continue;
            }

            double tan_t = std::tan(theta0);
            double tan2 = tan_t * tan_t;
            double A = dx*dx + dy*dy - dz*dz * tan2;
            double B_half = px*dx + py*dy - pz*dz * tan2;
            double C = px*px + py*py - pz*pz * tan2;

            double det = B_half * B_half - A * C;
            if(det < 0) continue;
            double sq = std::sqrt(std::fmax(det, 0.0));

            for(int root = 0; root < 2; ++root) {
                double t;
                if(std::fabs(A) > GEOMETRY_PRECISION) {
                    t = (root == 0) ? (-B_half - sq) / A : (-B_half + sq) / A;
                } else if(root == 0 && std::fabs(B_half) > GEOMETRY_PRECISION) {
                    t = -C / (2.0 * B_half);
                } else {
                    continue;
                }

                double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
                if(theta0 < M_PI / 2.0 && hz < -GEOMETRY_PRECISION) continue;
                if(theta0 > M_PI / 2.0 && hz > GEOMETRY_PRECISION) continue;

                double h2 = hx*hx + hy*hy + hz*hz;
                if(h2 < GEOMETRY_PRECISION * 1e3 * GEOMETRY_PRECISION * 1e3) {
                    if(theta_apex_emitted) continue;
                    theta_apex_emitted = true;
                    double eps = GEOMETRY_PRECISION * 1e3;
                    bool entering = ThetaInRange(hx + eps*dx, hy + eps*dy, hz + eps*dz,
                                                 start_theta_, delta_theta_);
                    add_hit(t, entering, 2);
                    continue;
                }

                double gnx = hx, gny = hy, gnz = -hz * tan2;
                double g_dot_d = gnx*dx + gny*dy + gnz*dz;
                if(theta0 > M_PI / 2.0) g_dot_d = -g_dot_d;
                bool entering;
                if(face == 0) {
                    entering = (g_dot_d > 0);
                } else {
                    entering = (g_dot_d < 0);
                }
                add_hit(t, entering, 2);
            }
        }
    }

    if(n_all == 0) return {};

    // No angular cuts: shell hits are the final result
    if(!has_phi_cut_ && !has_theta_cut_) {
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

    // Sort all hits by distance
    std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // Determine initial states at t = -infinity.
    // Shell: always starts outside (finite closed surface).
    // Phi wedge: determine from first wedge hit or PhiInRange check.
    // Theta band: evaluate ThetaInRange along -direction (the t=-inf limit).
    bool in_shell = false;
    bool in_phi = !has_phi_cut_; // no phi cut means always "inside" the wedge
    bool in_theta = !has_theta_cut_; // no theta cut means always "inside" the band

    if(has_phi_cut_) {
        in_phi = InitialPhiState(px, py, dx, dy, start_phi_, delta_phi_);
    }

    if(has_theta_cut_) {
        in_theta = ThetaInRange(-dx, -dy, -dz, start_theta_, delta_theta_);
    }

    bool was_inside = in_shell && in_phi && in_theta;

    // CSG intersection walk
    std::vector<Intersection> result;
    for(int i = 0; i < n_all; ++i) {
        switch(all_hits[i].source) {
            case 0: in_shell = all_hits[i].entering; break;
            case 1: in_phi = all_hits[i].entering; break;
            case 2: in_theta = all_hits[i].entering; break;
        }
        bool now_inside = in_shell && in_phi && in_theta;
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

