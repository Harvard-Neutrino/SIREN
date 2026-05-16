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
    auto solve_shell = [&](double r, bool is_outer) {
        double A = pp - r * r;
        double det = pd * pd - A;
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
            // Must be on the correct half-plane (outward from z-axis)
            if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;
            bool entering = (n_dot_d < 0);
            add_hit(t, entering, 1);
        }
    }

    // ---- Infinite theta band boundaries (no radial/phi bounds) ----
    if(has_theta_cut_) {
        for(int face = 0; face < 2; ++face) {
            double theta0 = start_theta_ + face * delta_theta_;
            // Degenerate cone at theta=0 or theta=pi (z-axis): skip
            if(theta0 < 1e-12 || std::fabs(theta0 - M_PI) < 1e-12) continue;

            if(std::fabs(theta0 - M_PI / 2.0) < 1e-12) {
                // theta = pi/2 is the z=0 plane
                if(std::fabs(dz) < GEOMETRY_PRECISION) continue;
                double t = -pz / dz;
                bool entering = (face == 0) ? (dz < 0) : (dz > 0);
                add_hit(t, entering, 2);
                continue;
            }

            // General cone at polar angle theta0
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
                // Correct hemisphere check
                if(theta0 < M_PI / 2.0 && hz < -GEOMETRY_PRECISION) continue;
                if(theta0 > M_PI / 2.0 && hz > GEOMETRY_PRECISION) continue;

                // Entering/exiting the theta band
                double gnx = hx, gny = hy, gnz = -hz * tan2;
                double g_dot_d = gnx*dx + gny*dy + gnz*dz;
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
    // Theta band: determine from first theta hit or ThetaInRange check.
    bool in_shell = false;
    bool in_phi = !has_phi_cut_; // no phi cut means always "inside" the wedge
    bool in_theta = !has_theta_cut_; // no theta cut means always "inside" the band

    if(has_phi_cut_) {
        bool found = false;
        for(int i = 0; i < n_all; ++i) {
            if(all_hits[i].source == 1) {
                in_phi = !all_hits[i].entering;
                found = true;
                break;
            }
        }
        if(!found) {
            in_phi = PhiInRange(px, py, start_phi_, delta_phi_);
        }
    }

    if(has_theta_cut_) {
        bool found = false;
        for(int i = 0; i < n_all; ++i) {
            if(all_hits[i].source == 2) {
                in_theta = !all_hits[i].entering;
                found = true;
                break;
            }
        }
        if(!found) {
            in_theta = ThetaInRange(px, py, pz, start_theta_, delta_theta_);
        }
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

