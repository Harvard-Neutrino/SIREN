#include "SIREN/geometry/Torus.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

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

// =========================================================================
// Cubic and quartic solvers for ray-torus intersection
// =========================================================================

// Solve depressed cubic t^3 + p*t + q = 0 using trigonometric method.
// Returns the number of real roots found (1 or 3).
// Roots are stored in roots[0..n-1], unsorted.
int SolveDepressedCubic(double p, double q, double roots[3]) {
    double disc = -4.0 * p * p * p - 27.0 * q * q;

    if(disc >= 0) {
        // Three real roots (or repeated) -- trigonometric method
        double mp3 = -p / 3.0;
        if(mp3 < 0) mp3 = 0;
        double r = std::sqrt(mp3);
        double cos_arg = 0;
        if(r > 0) {
            cos_arg = -q / (2.0 * r * r * r);
            if(cos_arg > 1.0) cos_arg = 1.0;
            if(cos_arg < -1.0) cos_arg = -1.0;
        }
        double theta = std::acos(cos_arg);
        roots[0] = 2.0 * r * std::cos(theta / 3.0);
        roots[1] = 2.0 * r * std::cos((theta + 2.0 * M_PI) / 3.0);
        roots[2] = 2.0 * r * std::cos((theta + 4.0 * M_PI) / 3.0);
        return 3;
    } else {
        // One real root -- Cardano's formula
        double sq = std::sqrt(q * q / 4.0 + p * p * p / 27.0);
        double u = std::cbrt(-q / 2.0 + sq);
        double v = std::cbrt(-q / 2.0 - sq);
        roots[0] = u + v;
        return 1;
    }
}

// Solve monic quartic x^4 + a*x^3 + b*x^2 + c*x + d = 0
// using Ferrari's method. Returns the number of real roots (0-4).
// Roots are stored in roots[0..n-1], unsorted.
int SolveQuartic(double a, double b, double c, double d, double roots[4]) {
    // Depress the quartic: substitute x = t - a/4
    // Gives t^4 + p*t^2 + q*t + r = 0
    double a2 = a * a;
    double p = b - 3.0 * a2 / 8.0;
    double q = c - a * b / 2.0 + a2 * a / 8.0;
    double r = d - a * c / 4.0 + a2 * b / 16.0 - 3.0 * a2 * a2 / 256.0;

    int n_roots = 0;

    // Biquadratic threshold: scale relative to coefficient magnitudes
    // so that far-field rays (large coefficients) still trigger this path
    // when the quartic is structurally biquadratic (q = 0 analytically)
    double coeff_scale = std::fmax(1.0, std::fmax(std::fabs(p), std::sqrt(std::fmax(std::fabs(r), 0.0))));
    if(std::fabs(q) < 1e-10 * coeff_scale) {
        // Biquadratic: t^4 + p*t^2 + r = 0
        double disc = p * p - 4.0 * r;
        double disc_tol = 1e-10 * std::fmax(p * p, std::fabs(r));
        if(disc < -disc_tol) return 0;
        if(disc < 0) disc = 0;
        double sq = std::sqrt(disc);
        double u1 = (-p + sq) / 2.0;
        double u2 = (-p - sq) / 2.0;
        double u_tol = 1e-10 * std::fmax(std::fabs(p), sq);
        if(u1 >= -u_tol) {
            if(u1 < 0) u1 = 0;
            double s = std::sqrt(u1);
            roots[n_roots++] = s - a / 4.0;
            roots[n_roots++] = -s - a / 4.0;
        }
        if(u2 >= -u_tol && std::fabs(u2 - u1) > u_tol) {
            if(u2 < 0) u2 = 0;
            double s = std::sqrt(u2);
            roots[n_roots++] = s - a / 4.0;
            roots[n_roots++] = -s - a / 4.0;
        }
        return n_roots;
    }

    // Ferrari's resolvent cubic
    double rp = -p / 2.0;
    double rq = -r;
    double rs = p * r / 2.0 - q * q / 8.0;
    double P = rq - rp * rp / 3.0;
    double Q = rs - rp * rq / 3.0 + 2.0 * rp * rp * rp / 27.0;

    double cubic_roots[3];
    int nc = SolveDepressedCubic(P, Q, cubic_roots);

    // Pick the largest real root of the resolvent cubic
    double y = cubic_roots[0] - rp / 3.0;
    for(int i = 1; i < nc; ++i) {
        double yi = cubic_roots[i] - rp / 3.0;
        if(yi > y) y = yi;
    }

    // Factor the depressed quartic into two quadratics using y:
    // t^4 + p*t^2 + q*t + r = (t^2 + s*t + u)(t^2 - s*t + v)
    // where s^2 = 2*y - p, u = y - q/(2*s), v = y + q/(2*s)
    double s2 = 2.0 * y - p;
    if(s2 < 0) s2 = 0;
    double s = std::sqrt(s2);

    double u, v;
    if(s > 1e-14) {
        u = y - q / (2.0 * s);
        v = y + q / (2.0 * s);
    } else {
        u = y;
        v = y;
        s = 0;
    }

    // Solve the two quadratics
    double disc1 = s * s - 4.0 * u;
    if(disc1 >= 0) {
        double sq1 = std::sqrt(disc1);
        roots[n_roots++] = (-s + sq1) / 2.0 - a / 4.0;
        roots[n_roots++] = (-s - sq1) / 2.0 - a / 4.0;
    }
    double disc2 = s * s - 4.0 * v;
    if(disc2 >= 0) {
        double sq2 = std::sqrt(disc2);
        roots[n_roots++] = (s + sq2) / 2.0 - a / 4.0;
        roots[n_roots++] = (s - sq2) / 2.0 - a / 4.0;
    }

    return n_roots;
}

} // anonymous namespace


// =========================================================================
// Constructors
// =========================================================================

Torus::Torus()
    : Geometry((std::string)("Torus"))
    , rtor_(0.0)
    , rmax_(0.0)
    , rmin_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
}

Torus::Torus(double rtor, double rmax, double rmin)
    : Geometry((std::string)("Torus"))
    , rtor_(rtor)
    , rmax_(rmax)
    , rmin_(rmin)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(rtor_ <= 0) {
        throw std::invalid_argument("Torus major radius must be positive!");
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("Torus outer tube radius must be positive!");
    }
    if(rmin_ < 0) {
        throw std::invalid_argument("Torus inner tube radius must be non-negative!");
    }
    if(rmin_ >= rmax_) {
        throw std::invalid_argument("Torus inner tube radius must be less than outer tube radius!");
    }
}

Torus::Torus(Placement const & placement)
    : Geometry((std::string)("Torus"), placement)
    , rtor_(0.0)
    , rmax_(0.0)
    , rmin_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
}

Torus::Torus(Placement const & placement, double rtor, double rmax, double rmin)
    : Geometry((std::string)("Torus"), placement)
    , rtor_(rtor)
    , rmax_(rmax)
    , rmin_(rmin)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(rtor_ <= 0) {
        throw std::invalid_argument("Torus major radius must be positive!");
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("Torus outer tube radius must be positive!");
    }
    if(rmin_ < 0) {
        throw std::invalid_argument("Torus inner tube radius must be non-negative!");
    }
    if(rmin_ >= rmax_) {
        throw std::invalid_argument("Torus inner tube radius must be less than outer tube radius!");
    }
}

Torus::Torus(double rtor, double rmax, double rmin, double start_phi, double delta_phi)
    : Geometry((std::string)("Torus"))
    , rtor_(rtor)
    , rmax_(rmax)
    , rmin_(rmin)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(rtor_ <= 0) {
        throw std::invalid_argument("Torus major radius must be positive!");
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("Torus outer tube radius must be positive!");
    }
    if(rmin_ < 0) {
        throw std::invalid_argument("Torus inner tube radius must be non-negative!");
    }
    if(rmin_ >= rmax_) {
        throw std::invalid_argument("Torus inner tube radius must be less than outer tube radius!");
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Torus delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Torus::Torus(Placement const & placement, double rtor, double rmax, double rmin,
             double start_phi, double delta_phi)
    : Geometry((std::string)("Torus"), placement)
    , rtor_(rtor)
    , rmax_(rmax)
    , rmin_(rmin)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(rtor_ <= 0) {
        throw std::invalid_argument("Torus major radius must be positive!");
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("Torus outer tube radius must be positive!");
    }
    if(rmin_ < 0) {
        throw std::invalid_argument("Torus inner tube radius must be non-negative!");
    }
    if(rmin_ >= rmax_) {
        throw std::invalid_argument("Torus inner tube radius must be less than outer tube radius!");
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Torus delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Torus::Torus(const Torus& torus)
    : Geometry(torus)
    , rtor_(torus.rtor_)
    , rmax_(torus.rmax_)
    , rmin_(torus.rmin_)
    , start_phi_(torus.start_phi_)
    , delta_phi_(torus.delta_phi_)
    , has_phi_cut_(torus.has_phi_cut_) {
}

// =========================================================================
// Geometry interface
// =========================================================================

void Torus::swap(Geometry& geometry) {
    Torus* torus = dynamic_cast<Torus*>(&geometry);
    if(!torus) return;

    Geometry::swap(*torus);
    std::swap(rtor_, torus->rtor_);
    std::swap(rmax_, torus->rmax_);
    std::swap(rmin_, torus->rmin_);
    std::swap(start_phi_, torus->start_phi_);
    std::swap(delta_phi_, torus->delta_phi_);
    std::swap(has_phi_cut_, torus->has_phi_cut_);
}

Torus& Torus::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Torus* torus = dynamic_cast<const Torus*>(&geometry);
        if(!torus) return *this;
        Torus tmp(*torus);
        swap(tmp);
    }
    return *this;
}

bool Torus::equal(const Geometry& geometry) const {
    const Torus* torus = dynamic_cast<const Torus*>(&geometry);
    if(!torus) return false;
    return rtor_ == torus->rtor_ && rmax_ == torus->rmax_ && rmin_ == torus->rmin_
        && start_phi_ == torus->start_phi_ && delta_phi_ == torus->delta_phi_;
}

bool Torus::less(const Geometry& geometry) const {
    const Torus* torus = dynamic_cast<const Torus*>(&geometry);
    if(!torus) return false;
    return std::tie(rtor_, rmax_, rmin_, start_phi_, delta_phi_)
         < std::tie(torus->rtor_, torus->rmax_, torus->rmin_, torus->start_phi_, torus->delta_phi_);
}

void Torus::print(std::ostream& os) const {
    os << "Major radius: " << rtor_
       << "\tOuter tube radius: " << rmax_
       << "\tInner tube radius: " << rmin_;
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    os << '\n';
}

// =========================================================================
// ComputeIntersections
//
// Ray-torus intersection. The implicit torus equation
//   (sqrt(x^2 + y^2) - R)^2 + z^2 = r^2
// after squaring becomes the quartic surface:
//   (x^2 + y^2 + z^2 + R^2 - r^2)^2 = 4*R^2*(x^2 + y^2)
//
// Substituting ray P + t*D and expanding gives a monic quartic in t.
// We solve it with Ferrari's method.
// =========================================================================

std::vector<Geometry::Intersection> Torus::ComputeIntersections(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const {

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    double R = rtor_;

    // At most 4 outer + 4 inner + phi-plane hits = 16
    Intersection hits[16];
    int n_hits = 0;

    auto solve_surface = [&](double r, bool invert_entering) {
        // Shift the ray origin to the closest approach point along the ray.
        // This keeps the quartic coefficients small for far-field rays,
        // avoiding catastrophic cancellation in the depressed quartic.
        double t_shift = -(px*dx + py*dy + pz*dz);
        double qx = px + t_shift * dx;
        double qy = py + t_shift * dy;
        double qz = pz + t_shift * dz;

        double qq = qx*qx + qy*qy + qz*qz;
        double qd = qx*dx + qy*dy + qz*dz; // ~ 0 by construction
        double Qxy2 = qx*qx + qy*qy;
        double Dxy2 = dx*dx + dy*dy;
        double QDxy = qx*dx + qy*dy;

        double alpha = qq + R*R - r*r;

        // Quartic coefficients: t^4 + c3*t^3 + c2*t^2 + c1*t + c0 = 0
        double c3 = 4.0 * qd;
        double c2 = 4.0 * qd*qd + 2.0 * alpha - 4.0 * R*R * Dxy2;
        double c1 = 4.0 * qd * alpha - 8.0 * R*R * QDxy;
        double c0 = alpha * alpha - 4.0 * R*R * Qxy2;

        double raw_roots[4];
        int n_raw = SolveQuartic(c3, c2, c1, c0, raw_roots);
        for(int i = 0; i < n_raw; ++i) raw_roots[i] += t_shift;

        // Validate and deduplicate roots.
        // Ferrari's method can produce duplicate roots (same root from both
        // quadratics) and the quartic comes from squaring the torus equation
        // which can theoretically introduce extraneous roots.
        // For each root: verify it satisfies the torus equation, then deduplicate.
        std::sort(raw_roots, raw_roots + n_raw);
        double roots[4];
        int n = 0;
        for(int i = 0; i < n_raw; ++i) {
            double t = raw_roots[i];
            double ix = px + t * dx;
            double iy = py + t * dy;
            double iz = pz + t * dz;
            // Evaluate torus equation: (sqrt(x^2+y^2) - R)^2 + z^2 - r^2
            double rxy = std::sqrt(ix*ix + iy*iy);
            double residual = (rxy - R) * (rxy - R) + iz * iz - r * r;
            // Scale tolerance by r^2 for relative error
            double tol = 1e-4 * (r * r + R * R);
            if(std::fabs(residual) > tol) continue; // not a true intersection
            // Deduplicate: skip if too close to the previous accepted root
            double dedup_tol = 1e-6 * (1.0 + std::fabs(t));
            if(n > 0 && std::fabs(t - roots[n - 1]) < dedup_tol) continue;
            roots[n++] = t;
        }

        for(int i = 0; i < n; ++i) {
            double t = roots[i];
            if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

            double ix = px + t * dx;
            double iy = py + t * dy;
            double iz = pz + t * dz;

            if(has_phi_cut_ && !PhiInRange(ix, iy, start_phi_, delta_phi_)) continue;

            // Determine entering/exiting using the torus surface normal.
            // The outward normal at a point on the torus surface is:
            //   n = (point - R * (Pxy_hat, 0))
            // where Pxy_hat is the unit vector from the z-axis to the point in the xy-plane.
            double rxy = std::sqrt(ix*ix + iy*iy);
            double nx, ny, nz;
            if(rxy > GEOMETRY_PRECISION) {
                double scale = R / rxy;
                nx = ix - scale * ix;  // = ix * (1 - R/rxy)
                ny = iy - scale * iy;  // = iy * (1 - R/rxy)
                nz = iz;
            } else {
                // Point is on the z-axis (degenerate case)
                nx = ix;
                ny = iy;
                nz = iz;
            }
            double n_dot_d = nx * dx + ny * dy + nz * dz;
            bool entering = n_dot_d < 0;
            if(invert_entering) entering = !entering;

            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = entering;
            hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
            n_hits++;
        }
    };

    // Outer surface
    solve_surface(rmax_, false);
    // Inner surface (hollow torus)
    if(rmin_ > 0) {
        solve_surface(rmin_, true);
    }

    // Phi-boundary plane intersections
    if(has_phi_cut_) {
        for(int face = 0; face < 2; ++face) {
            double alpha = start_phi_ + face * delta_phi_;
            double ca = std::cos(alpha), sa = std::sin(alpha);
            // Outward-pointing normal (away from phi range interior)
            // Face 0 (start_phi): outward = clockwise from radial
            // Face 1 (end_phi): outward = counter-clockwise from radial
            double nx, ny;
            if(face == 0) { nx = sa; ny = -ca; }
            else { nx = -sa; ny = ca; }
            double n_dot_d = nx*dx + ny*dy;
            if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue;
            double n_dot_p = nx*px + ny*py;
            double t = -n_dot_p / n_dot_d;
            if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

            double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
            // Must be on outward side of z-axis
            if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;
            // Must be within torus cross-section
            double rxy = std::sqrt(hx*hx + hy*hy);
            double d2 = (rxy - R)*(rxy - R) + hz*hz;
            if(d2 > rmax_*rmax_ + GEOMETRY_PRECISION) continue;
            if(rmin_ > 0 && d2 < rmin_*rmin_ - GEOMETRY_PRECISION) continue;
            // Entering
            bool entering = (n_dot_d < 0);
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = entering;
            hits[n_hits].position = siren::math::Vector3D(hx, hy, hz);
            n_hits++;
        }
    }

    if(n_hits == 0) return {};
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    // Phi-cut edge fixup: rays grazing the boundary where a phi plane meets
    // the curved surface can produce an odd hit count due to one surface's
    // intersection being just inside tolerance and its counterpart just outside.
    // Discard the unpaired boundary hit to maintain even parity.
    if(has_phi_cut_ && n_hits % 2 != 0) {
        // Find and remove the hit closest to a phi boundary plane
        int best = 0;
        double best_dist = std::numeric_limits<double>::max();
        for(int i = 0; i < n_hits; ++i) {
            double hx = hits[i].position.GetX(), hy = hits[i].position.GetY();
            for(int f = 0; f < 2; ++f) {
                double alpha = start_phi_ + f * delta_phi_;
                double nx = std::sin(alpha), ny = -std::cos(alpha);
                double plane_dist = std::fabs(nx*hx + ny*hy);
                if(plane_dist < best_dist) {
                    best_dist = plane_dist;
                    best = i;
                }
            }
        }
        // Remove the hit at index best by shifting
        for(int i = best; i < n_hits - 1; ++i) hits[i] = hits[i + 1];
        n_hits--;
    }
    return {hits, hits + n_hits};
}

// =========================================================================
// GetBoundingBox
// =========================================================================

AABB Torus::GetBoundingBox() const {
    double extent_xy = rtor_ + rmax_;
    double extent_z = rmax_;
    return AABB(
        math::Vector3D(-extent_xy, -extent_xy, -extent_z),
        math::Vector3D( extent_xy,  extent_xy,  extent_z)
    );
}

} // namespace geometry
} // namespace siren
