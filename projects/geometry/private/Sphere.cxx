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

namespace siren {
namespace geometry {

Sphere::Sphere()
    : Geometry((std::string)("Sphere"))
    , radius_(0.0)
    , inner_radius_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , start_theta_(0.0)
    , delta_theta_(M_PI) {
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
}

Sphere::Sphere(double radius, double inner_radius)
    : Geometry("Sphere")
    , radius_(radius)
    , inner_radius_(inner_radius)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , start_theta_(0.0)
    , delta_theta_(M_PI) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
}

Sphere::Sphere(Placement const & placement)
    : Geometry((std::string)("Sphere"), placement)
    , radius_(0.0)
    , inner_radius_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , start_theta_(0.0)
    , delta_theta_(M_PI) {
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
}

Sphere::Sphere(Placement const & placement, double radius, double inner_radius)
    : Geometry((std::string)("Sphere"), placement)
    , radius_(radius)
    , inner_radius_(inner_radius)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , start_theta_(0.0)
    , delta_theta_(M_PI) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
}

Sphere::Sphere(double radius, double inner_radius,
               double start_phi, double delta_phi,
               double start_theta, double delta_theta)
    : Geometry("Sphere")
    , radius_(radius)
    , inner_radius_(inner_radius)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi)
    , start_theta_(start_theta)
    , delta_theta_(delta_theta) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9)
        throw std::invalid_argument("Sphere delta_phi must be in (0, 2*pi]!");
    if(start_theta_ < -1e-9 || start_theta_ > M_PI + 1e-9)
        throw std::invalid_argument("Sphere start_theta must be in [0, pi]!");
    if(delta_theta_ <= 0 || start_theta_ + delta_theta_ > M_PI + 1e-9)
        throw std::invalid_argument("Sphere start_theta + delta_theta must be in (0, pi]!");
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
}

Sphere::Sphere(Placement const & placement, double radius, double inner_radius,
               double start_phi, double delta_phi,
               double start_theta, double delta_theta)
    : Geometry((std::string)("Sphere"), placement)
    , radius_(radius)
    , inner_radius_(inner_radius)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi)
    , start_theta_(start_theta)
    , delta_theta_(delta_theta) {
    if(inner_radius_ > radius_) std::swap(inner_radius_, radius_);
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9)
        throw std::invalid_argument("Sphere delta_phi must be in (0, 2*pi]!");
    if(start_theta_ < -1e-9 || start_theta_ > M_PI + 1e-9)
        throw std::invalid_argument("Sphere start_theta must be in [0, pi]!");
    if(delta_theta_ <= 0 || start_theta_ + delta_theta_ > M_PI + 1e-9)
        throw std::invalid_argument("Sphere start_theta + delta_theta must be in (0, pi]!");
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    has_theta_cut_ = (start_theta_ > 1e-9 || delta_theta_ < M_PI - 1e-9);
}

Sphere::Sphere(const Sphere& sphere)
    : Geometry(sphere)
    , radius_(sphere.radius_)
    , inner_radius_(sphere.inner_radius_)
    , start_phi_(sphere.start_phi_)
    , delta_phi_(sphere.delta_phi_)
    , start_theta_(sphere.start_theta_)
    , delta_theta_(sphere.delta_theta_)
    , has_phi_cut_(sphere.has_phi_cut_)
    , has_theta_cut_(sphere.has_theta_cut_) {
}



/*Sphere::Sphere(const nlohmann::json& config)
  : Geometry(config)
  {
  if(not config.is_object()) throw std::invalid_argument("No json object found.");
  if(not config.at("outer_radius").is_number())
  throw std::invalid_argument("Outer radius is not a number.");

  config["outer_radius"].get_to(radius_);
  inner_radius_ = config.value("inner_radius", 0);

  if(inner_radius_ < 0) throw std::logic_error("inner radius must be >= 0");
  if(radius_ < inner_radius_)
  throw std::logic_error("radius must be larger than inner radius");
  }*/

// ------------------------------------------------------------------------- //
void Sphere::swap(Geometry& geometry)
{
    Sphere* sphere = dynamic_cast<Sphere*>(&geometry);
    if (!sphere)
    {
        //log_warn("Cannot swap Sphere!");
        return;
    }

    Geometry::swap(*sphere);

    std::swap(inner_radius_, sphere->inner_radius_);
    std::swap(radius_, sphere->radius_);
    std::swap(start_phi_, sphere->start_phi_);
    std::swap(delta_phi_, sphere->delta_phi_);
    std::swap(start_theta_, sphere->start_theta_);
    std::swap(delta_theta_, sphere->delta_theta_);
    std::swap(has_phi_cut_, sphere->has_phi_cut_);
    std::swap(has_theta_cut_, sphere->has_theta_cut_);
}

//------------------------------------------------------------------------- //
Sphere& Sphere::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);
        if (!sphere)
        {
            //log_warn("Cannot assign Sphere!");
            return *this;
        }

        Sphere tmp(*sphere);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Sphere::equal(const Geometry& geometry) const
{
    const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);

    if(!sphere)
        return false;
    else if(inner_radius_ != sphere->inner_radius_)
        return false;
    else if(radius_ != sphere->radius_)
        return false;
    else if(start_phi_ != sphere->start_phi_)
        return false;
    else if(delta_phi_ != sphere->delta_phi_)
        return false;
    else if(start_theta_ != sphere->start_theta_)
        return false;
    else if(delta_theta_ != sphere->delta_theta_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Sphere::less(const Geometry& geometry) const
{
    const Sphere* sphere = dynamic_cast<const Sphere*>(&geometry);
    if(!sphere) return false;

    return std::tie(inner_radius_, radius_, start_phi_, delta_phi_, start_theta_, delta_theta_)
         < std::tie(sphere->inner_radius_, sphere->radius_, sphere->start_phi_, sphere->delta_phi_, sphere->start_theta_, sphere->delta_theta_);
}

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

    // Max intersections: 4 (sphere surfaces) + 4 (phi planes) + 8 (theta cones) = 16
    Intersection hits[16];
    int n_hits = 0;

    auto add_hit = [&](double t, bool entering) {
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        hits[n_hits].distance = t;
        hits[n_hits].hierarchy = 0;
        hits[n_hits].entering = entering;
        hits[n_hits].position = siren::math::Vector3D(px + t*dx, py + t*dy, pz + t*dz);
        n_hits++;
    };

    // ---- Spherical surface intersections ----
    // For each shell (outer, inner), solve the quadratic and filter by angular range.

    auto solve_shell = [&](double r, bool is_outer) {
        double A = pp - r * r;
        double det = pd * pd - A;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double t1 = -pd - sq;
        double t2 = -pd + sq;
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? t1 : t2;
            double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
            // Angular filter
            if(has_phi_cut_ && !PhiInRange(hx, hy, start_phi_, delta_phi_)) continue;
            if(has_theta_cut_ && !ThetaInRange(hx, hy, hz, start_theta_, delta_theta_)) continue;
            bool entering = is_outer ? (k == 0) : (k == 1);
            add_hit(t, entering);
        }
    };

    solve_shell(radius_, true);
    if(inner_radius_ > 0) {
        solve_shell(inner_radius_, false);
    }

    // ---- Phi boundary planes ----
    if(has_phi_cut_) {
        double r2_outer = radius_ * radius_;
        double r2_inner = inner_radius_ * inner_radius_;
        // Two planes at start_phi and start_phi + delta_phi
        for(int face = 0; face < 2; ++face) {
            double alpha = start_phi_ + face * delta_phi_;
            double ca = std::cos(alpha), sa = std::sin(alpha);
            // Outward-pointing normal (away from the phi range interior).
            // Face 0 (start_phi boundary): outward = clockwise from radial = (sin(alpha), -cos(alpha), 0)
            // Face 1 (end_phi boundary): outward = counter-clockwise from radial = (-sin(alpha), cos(alpha), 0)
            double nx, ny;
            if(face == 0) {
                nx = sa; ny = -ca;
            } else {
                nx = -sa; ny = ca;
            }
            double n_dot_d = nx*dx + ny*dy;
            if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue; // parallel
            double n_dot_p = nx*px + ny*py;
            double t = -n_dot_p / n_dot_d;
            if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

            double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
            // Must be on outward side of z-axis (the half-plane, not the full plane)
            if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;
            // Radial bounds
            double r2_hit = hx*hx + hy*hy + hz*hz;
            if(r2_hit > r2_outer + GEOMETRY_PRECISION) continue;
            if(r2_hit < r2_inner - GEOMETRY_PRECISION) continue;
            // Theta filter
            if(has_theta_cut_ && !ThetaInRange(hx, hy, hz, start_theta_, delta_theta_)) continue;
            // Entering: ray crosses into the phi range
            bool entering = (n_dot_d < 0);
            add_hit(t, entering);
        }
    }

    // ---- Theta boundary cones ----
    if(has_theta_cut_) {
        double r2_outer = radius_ * radius_;
        double r2_inner = inner_radius_ * inner_radius_;
        // Two cones at start_theta and start_theta + delta_theta
        for(int face = 0; face < 2; ++face) {
            double theta0 = start_theta_ + face * delta_theta_;

            // Special case: theta0 ~ 0 or theta0 ~ pi => degenerate cone (z-axis)
            if(theta0 < 1e-12 || std::fabs(theta0 - M_PI) < 1e-12) continue;

            // Special case: theta0 ~ pi/2 => flat disk at z = 0
            if(std::fabs(theta0 - M_PI / 2.0) < 1e-12) {
                if(std::fabs(dz) < GEOMETRY_PRECISION) continue; // parallel to disk
                double t = -pz / dz;
                if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
                double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
                double r2_hit = hx*hx + hy*hy + hz*hz;
                if(r2_hit > r2_outer + GEOMETRY_PRECISION) continue;
                if(r2_hit < r2_inner - GEOMETRY_PRECISION) continue;
                if(has_phi_cut_ && !PhiInRange(hx, hy, start_phi_, delta_phi_)) continue;
                // Entering: for face 0 (upper theta boundary), entering means
                // going from lower theta (above) to higher theta (below)
                // Normal of z=0 plane points in +z direction
                bool entering = (face == 0) ? (dz < 0) : (dz > 0);
                add_hit(t, entering);
                continue;
            }

            // General case: cone at polar angle theta0
            // Equation: x^2 + y^2 - z^2 * tan^2(theta0) = 0
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
                if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

                double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
                // Check z-sign: for theta0 < pi/2, the cone is in the +z hemisphere
                // for theta0 > pi/2, the cone is in the -z hemisphere
                if(theta0 < M_PI / 2.0 && hz < -GEOMETRY_PRECISION) continue;
                if(theta0 > M_PI / 2.0 && hz > GEOMETRY_PRECISION) continue;
                // Radial bounds
                double r2_hit = hx*hx + hy*hy + hz*hz;
                if(r2_hit > r2_outer + GEOMETRY_PRECISION) continue;
                if(r2_hit < r2_inner - GEOMETRY_PRECISION) continue;
                // Phi filter
                if(has_phi_cut_ && !PhiInRange(hx, hy, start_phi_, delta_phi_)) continue;
                // Entering: the cone normal at the hit point.
                // The cone surface x^2+y^2-z^2*tan2=0 has gradient (2x, 2y, -2z*tan2).
                // For face 0 (start_theta), the inward-pointing normal (into the theta range)
                // points toward higher theta (away from z-axis for theta < pi/2).
                // For face 1 (end theta), inward normal points toward lower theta.
                double gnx = hx, gny = hy, gnz = -hz * tan2;
                double g_dot_d = gnx*dx + gny*dy + gnz*dz;
                // face 0: entering when moving INTO the theta range (away from cone on the lower-theta side)
                // face 1: entering when moving INTO the theta range (away from cone on the higher-theta side)
                bool entering;
                if(face == 0) {
                    // Start-theta cone: inward = toward higher theta = outward from z-axis
                    entering = (g_dot_d > 0); // moving outward from cone = into the theta range
                } else {
                    // End-theta cone: inward = toward lower theta = inward toward z-axis
                    entering = (g_dot_d < 0);
                }
                add_hit(t, entering);
            }
        }
    }

    if(n_hits == 0) return {};
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    // Angular-cut edge fixup: rays grazing the boundary where a phi plane
    // or theta cone meets the curved surface can produce an odd hit count
    // due to numerical tolerance mismatches. Discard the unpaired hit
    // closest to a boundary surface to maintain even parity.
    if((has_phi_cut_ || has_theta_cut_) && n_hits % 2 != 0) {
        int best = 0;
        double best_dist = std::numeric_limits<double>::max();
        for(int i = 0; i < n_hits; ++i) {
            double hx = hits[i].position.GetX(), hy = hits[i].position.GetY();
            double hz = hits[i].position.GetZ();
            double min_d = best_dist;
            if(has_phi_cut_) {
                for(int f = 0; f < 2; ++f) {
                    double alpha = start_phi_ + f * delta_phi_;
                    double nx = std::sin(alpha), ny = -std::cos(alpha);
                    double d = std::fabs(nx*hx + ny*hy);
                    if(d < min_d) min_d = d;
                }
            }
            if(has_theta_cut_) {
                // Distance to theta-boundary cones
                double rxy2 = hx*hx + hy*hy;
                for(int f = 0; f < 2; ++f) {
                    double theta0 = start_theta_ + f * delta_theta_;
                    if(theta0 < 1e-12 || std::fabs(theta0 - M_PI) < 1e-12) continue;
                    double tan_t = std::tan(theta0);
                    double d = std::fabs(std::sqrt(rxy2) - std::fabs(hz) * tan_t);
                    if(d < min_d) min_d = d;
                }
            }
            if(min_d < best_dist) {
                best_dist = min_d;
                best = i;
            }
        }
        for(int i = best; i < n_hits - 1; ++i) hits[i] = hits[i + 1];
        n_hits--;
    }
    return {hits, hits + n_hits};
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

