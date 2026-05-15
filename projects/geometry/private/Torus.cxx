#include "SIREN/geometry/Torus.h"

#include <cmath>
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

    if(std::fabs(q) < 1e-14) {
        // Biquadratic: t^4 + p*t^2 + r = 0
        double disc = p * p - 4.0 * r;
        if(disc < -1e-14) return 0;
        if(disc < 0) disc = 0;
        double sq = std::sqrt(disc);
        double u1 = (-p + sq) / 2.0;
        double u2 = (-p - sq) / 2.0;
        if(u1 >= 0) {
            double s = std::sqrt(u1);
            roots[n_roots++] = s - a / 4.0;
            roots[n_roots++] = -s - a / 4.0;
        }
        if(u2 >= 0 && std::fabs(u2 - u1) > 1e-14) {
            double s = std::sqrt(u2);
            roots[n_roots++] = s - a / 4.0;
            roots[n_roots++] = -s - a / 4.0;
        }
        return n_roots;
    }

    // Ferrari's resolvent cubic: y^3 - p*y^2/2 - r*y + (p*r/2 - q^2/8) = 0
    // We solve the depressed form. First convert to depressed cubic in y:
    // Let y = m + p/6 to get m^3 + P*m + Q = 0
    double rp = -p / 2.0;          // coefficient of y^2 in resolvent (negated for standard form)
    double rq = -r;                 // coefficient of y^1 in resolvent (negated)
    double rs = p * r / 2.0 - q * q / 8.0; // constant in resolvent
    // Resolvent: y^3 + rp*y^2 + rq*y + rs = 0
    // Depress: y = m - rp/3
    double P = rq - rp * rp / 3.0;
    double Q = rs - rp * rq / 3.0 + 2.0 * rp * rp * rp / 27.0;

    double cubic_roots[3];
    int nc = SolveDepressedCubic(P, Q, cubic_roots);

    // Pick the largest real root of the resolvent cubic (most numerically stable)
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
        // s ~ 0: both quadratics have same linear term
        u = y;
        v = y;
        s = 0;
    }

    // Solve the two quadratics:
    // t^2 + s*t + u = 0
    double disc1 = s * s - 4.0 * u;
    if(disc1 >= 0) {
        double sq1 = std::sqrt(disc1);
        roots[n_roots++] = (-s + sq1) / 2.0 - a / 4.0;
        roots[n_roots++] = (-s - sq1) / 2.0 - a / 4.0;
    }
    // t^2 - s*t + v = 0
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
      , rmin_(0.0) {
}

Torus::Torus(double rtor, double rmax, double rmin)
    : Geometry((std::string)("Torus"))
    , rtor_(rtor)
    , rmax_(rmax)
      , rmin_(rmin) {
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
      , rmin_(0.0) {
}

Torus::Torus(Placement const & placement, double rtor, double rmax, double rmin)
    : Geometry((std::string)("Torus"), placement)
    , rtor_(rtor)
    , rmax_(rmax)
      , rmin_(rmin) {
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

Torus::Torus(const Torus& torus)
    : Geometry(torus)
    , rtor_(torus.rtor_)
    , rmax_(torus.rmax_)
      , rmin_(torus.rmin_) {
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
    return rtor_ == torus->rtor_ && rmax_ == torus->rmax_ && rmin_ == torus->rmin_;
}

bool Torus::less(const Geometry& geometry) const {
    const Torus* torus = dynamic_cast<const Torus*>(&geometry);
    if(!torus) return false;
    return std::tie(rtor_, rmax_, rmin_)
         < std::tie(torus->rtor_, torus->rmax_, torus->rmin_);
}

void Torus::print(std::ostream& os) const {
    os << "Major radius: " << rtor_
       << "\tOuter tube radius: " << rmax_
       << "\tInner tube radius: " << rmin_ << '\n';
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

    // At most 4 outer + 4 inner = 8 intersections for hollow torus
    Intersection hits[8];
    int n_hits = 0;

    auto solve_surface = [&](double r, bool invert_entering) {
        // Derived quantities
        double pp = px*px + py*py + pz*pz;  // P.P
        double pd = px*dx + py*dy + pz*dz;  // P.D
        double Pxy2 = px*px + py*py;        // Px^2 + Py^2
        double Dxy2 = dx*dx + dy*dy;        // Dx^2 + Dy^2
        double PDxy = px*dx + py*dy;        // Px*Dx + Py*Dy

        double alpha = pp + R*R - r*r;      // P.P + R^2 - r^2

        // Quartic coefficients: t^4 + c3*t^3 + c2*t^2 + c1*t + c0 = 0
        double c3 = 4.0 * pd;
        double c2 = 4.0 * pd*pd + 2.0 * alpha - 4.0 * R*R * Dxy2;
        double c1 = 4.0 * pd * alpha - 8.0 * R*R * PDxy;
        double c0 = alpha * alpha - 4.0 * R*R * Pxy2;

        double raw_roots[4];
        int n_raw = SolveQuartic(c3, c2, c1, c0, raw_roots);

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

    if(n_hits == 0) return {};
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
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
