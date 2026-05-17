#include "SIREN/geometry/Polycone.h"

#include <cmath>
#include <tuple>
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

void Polycone::validate() const {
    if(z_planes_.size() < 2) {
        throw std::runtime_error("Polycone requires at least 2 z-planes!");
    }
    if(z_planes_.size() != rmin_.size() || z_planes_.size() != rmax_.size()) {
        throw std::runtime_error("Polycone z_planes, rmin, and rmax vectors must have the same size!");
    }
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        if(rmin_[i] < 0 || rmax_[i] < 0) {
            throw std::runtime_error("Polycone radii must be non-negative!");
        }
        if(rmin_[i] > rmax_[i]) {
            throw std::runtime_error("Polycone inner radius must not exceed outer radius at each z-plane!");
        }
    }
    for(size_t i = 1; i < z_planes_.size(); ++i) {
        if(z_planes_[i] < z_planes_[i - 1]) {
            throw std::runtime_error("Polycone z-planes must be sorted in ascending order!");
        }
    }
}

Polycone::Polycone()
    : Geometry((std::string)("Polycone"))
    , z_planes_()
    , rmin_()
    , rmax_()
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    // Do nothing here
}

Polycone::Polycone(std::vector<double> const & z_planes,
                   std::vector<double> const & rmin,
                   std::vector<double> const & rmax)
    : Geometry((std::string)("Polycone"))
    , z_planes_(z_planes)
    , rmin_(rmin)
    , rmax_(rmax)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    validate();
}

Polycone::Polycone(Placement const & placement)
    : Geometry((std::string)("Polycone"), placement)
    , z_planes_()
    , rmin_()
    , rmax_()
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    // Do nothing here
}

Polycone::Polycone(Placement const & placement,
                   std::vector<double> const & z_planes,
                   std::vector<double> const & rmin,
                   std::vector<double> const & rmax)
    : Geometry((std::string)("Polycone"), placement)
    , z_planes_(z_planes)
    , rmin_(rmin)
    , rmax_(rmax)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    validate();
}

Polycone::Polycone(std::vector<double> const & z_planes,
                   std::vector<double> const & rmin,
                   std::vector<double> const & rmax,
                   double start_phi, double delta_phi)
    : Geometry((std::string)("Polycone"))
    , z_planes_(z_planes)
    , rmin_(rmin)
    , rmax_(rmax)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Polycone delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Polycone::Polycone(Placement const & placement,
                   std::vector<double> const & z_planes,
                   std::vector<double> const & rmin,
                   std::vector<double> const & rmax,
                   double start_phi, double delta_phi)
    : Geometry((std::string)("Polycone"), placement)
    , z_planes_(z_planes)
    , rmin_(rmin)
    , rmax_(rmax)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Polycone delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Polycone::Polycone(const Polycone& polycone)
    : Geometry(polycone)
    , z_planes_(polycone.z_planes_)
    , rmin_(polycone.rmin_)
    , rmax_(polycone.rmax_)
    , start_phi_(polycone.start_phi_)
    , delta_phi_(polycone.delta_phi_)
    , has_phi_cut_(polycone.has_phi_cut_) {
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Polycone::swap(Geometry& geometry) {
    Polycone* polycone = dynamic_cast<Polycone*>(&geometry);
    if(!polycone) {
        return;
    }

    Geometry::swap(*polycone);

    std::swap(z_planes_, polycone->z_planes_);
    std::swap(rmin_, polycone->rmin_);
    std::swap(rmax_, polycone->rmax_);
    std::swap(start_phi_, polycone->start_phi_);
    std::swap(delta_phi_, polycone->delta_phi_);
    std::swap(has_phi_cut_, polycone->has_phi_cut_);
}

// ------------------------------------------------------------------------- //
Polycone& Polycone::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Polycone* polycone = dynamic_cast<const Polycone*>(&geometry);
        if(!polycone) {
            return *this;
        }
        Polycone tmp(*polycone);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Polycone::equal(const Geometry& geometry) const
{
    const Polycone* polycone = dynamic_cast<const Polycone*>(&geometry);

    if(!polycone)
        return false;
    else if(z_planes_ != polycone->z_planes_)
        return false;
    else if(rmin_ != polycone->rmin_)
        return false;
    else if(rmax_ != polycone->rmax_)
        return false;
    else if(start_phi_ != polycone->start_phi_)
        return false;
    else if(delta_phi_ != polycone->delta_phi_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Polycone::less(const Geometry& geometry) const
{
    const Polycone* polycone = dynamic_cast<const Polycone*>(&geometry);
    if(!polycone) return false;

    return
        std::tie(z_planes_, rmin_, rmax_, start_phi_, delta_phi_)
        <
        std::tie(polycone->z_planes_, polycone->rmin_, polycone->rmax_, polycone->start_phi_, polycone->delta_phi_);
}

void Polycone::print(std::ostream& os) const
{
    os << "Polycone with " << z_planes_.size() << " z-planes:";
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        os << " [z=" << z_planes_[i]
           << " rmin=" << rmin_[i]
           << " rmax=" << rmax_[i] << "]";
    }
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    os << '\n';
}

// ------------------------------------------------------------------------- //
// ComputeIntersections
//
// The polycone is decomposed into N-1 conical frustum sections between
// adjacent z-planes. For each section, we compute ray intersections with
// the outer cone surface, inner cone surface (if hollow), and the end caps
// (top/bottom annular disks only for the first and last z-plane, since
// internal z-planes are shared boundaries, not real surfaces).
//
// All coordinates are in local frame. Distance can be negative (full-line
// intersections, same as Cone.cxx).
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

    // Estimate max intersections: each segment can produce up to 4 (2 outer +
    // 2 inner barrel), plus 2 endcaps at top/bottom, plus up to 2 per internal
    // z-plane. Conservative bound: 6 * num_segments.
    size_t max_hits = 6 * (z_planes_.size() - 1) + 4;

    // Fast path: use a stack buffer when the polycone is small enough.
    // Slow path: heap-allocate for large polycones (> ~5 segments).
    static constexpr int STACK_CAP = 32;
    Intersection stack_hits[STACK_CAP];
    std::vector<Intersection> heap_hits;
    bool use_heap = (max_hits > STACK_CAP);
    if(use_heap) heap_hits.reserve(max_hits);
    int n_hits = 0;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    // Write a hit to whichever buffer is active
    auto emit = [&](double t, bool entering) {
        Intersection isect;
        isect.distance = t;
        isect.hierarchy = 0;
        isect.entering = entering;
        isect.position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
        if(use_heap) {
            heap_hits.push_back(isect);
        } else {
            assert(n_hits < STACK_CAP);
            stack_hits[n_hits] = isect;
        }
        n_hits++;
    };

    // Cone surface entering test using the gradient normal.
    // For surface x^2 + y^2 = (a + b*z)^2, outward normal is (x, y, -b*(a+b*z)).
    // Entering when normal . direction < 0.
    auto entering_cone = [&](double a, double b) {
        double r_val = a + b * intersection_z;
        return (intersection_x * dx + intersection_y * dy - b * r_val * dz) < 0;
    };

    size_t n = z_planes_.size();

    // Process each frustum section between adjacent z-planes
    for(size_t seg = 0; seg + 1 < n; ++seg) {
        double z_lo = z_planes_[seg];
        double z_hi = z_planes_[seg + 1];

        // Skip degenerate zero-height sections
        if(z_hi <= z_lo) {
            continue;
        }

        // Horizontal ray cannot hit barrel surfaces outside its z-range.
        // Use half-open interval [z_lo, z_hi) so each boundary belongs to
        // exactly one section and horizontal rays at internal z-planes are
        // not skipped by both adjacent sections.
        if(std::fabs(dz) < GEOMETRY_PRECISION && (pz < z_lo || pz >= z_hi)) {
            continue;
        }

        double dz_sec = z_hi - z_lo;

        // --- Outer conical surface for this section ---
        {
            double rmax_lo = rmax_[seg];
            double rmax_hi = rmax_[seg + 1];

            // r_outer(z) = rmax_lo + (rmax_hi - rmax_lo) * (z - z_lo) / dz_sec
            //            = a + b * z   where:
            //   a = rmax_lo - (rmax_hi - rmax_lo) * z_lo / dz_sec
            //   b = (rmax_hi - rmax_lo) / dz_sec
            double b_outer = (rmax_hi - rmax_lo) / dz_sec;
            double a_outer = rmax_lo - b_outer * z_lo;

            // Surface equation: x^2 + y^2 = (a + b*z)^2
            // Substituting ray: x = px + t*dx, y = py + t*dy, z = pz + t*dz
            // gives:  A*t^2 + B*t + C = 0
            double A = dx * dx + dy * dy - b_outer * b_outer * dz * dz;
            double B = 2.0 * (px * dx + py * dy - b_outer * dz * (a_outer + b_outer * pz));
            double C = px * px + py * py - (a_outer + b_outer * pz) * (a_outer + b_outer * pz);

            if(!(dx == 0 && dy == 0 && b_outer == 0)) {
                double determinant = B * B - 4.0 * A * C;

                if(std::fabs(A) > GEOMETRY_PRECISION) {
                    if(determinant > 0) {
                        double sqrt_det = std::sqrt(determinant);
                        double t1 = (-B + sqrt_det) / (2.0 * A);
                        double t2 = (-B - sqrt_det) / (2.0 * A);

                        if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                            t1 = 0;
                        if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                            t2 = 0;

                        intersection_z = pz + t1 * dz;
                        if(intersection_z >= z_lo && intersection_z < z_hi) {
                            intersection_x = px + t1 * dx;
                            intersection_y = py + t1 * dy;
                            double r_at_z = a_outer + b_outer * intersection_z;
                            if(r_at_z >= 0) {
                                bool entering = entering_cone(a_outer, b_outer);
                                emit(t1, entering);
                            }
                        }

                        intersection_z = pz + t2 * dz;
                        if(intersection_z >= z_lo && intersection_z < z_hi) {
                            intersection_x = px + t2 * dx;
                            intersection_y = py + t2 * dy;
                            double r_at_z = a_outer + b_outer * intersection_z;
                            if(r_at_z >= 0) {
                                bool entering = entering_cone(a_outer, b_outer);
                                emit(t2, entering);
                            }
                        }
                    }
                } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                    // Linear case: ray parallel to cone surface
                    double t1 = -C / B;
                    if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                        t1 = 0;

                    intersection_z = pz + t1 * dz;
                    if(intersection_z >= z_lo && intersection_z < z_hi) {
                        intersection_x = px + t1 * dx;
                        intersection_y = py + t1 * dy;
                        double r_at_z = a_outer + b_outer * intersection_z;
                        if(r_at_z >= 0) {
                            bool entering = entering_cone(a_outer, b_outer);
                            emit(t1, entering);
                        }
                    }
                }
            }
        }

        // --- Inner conical surface for this section (hollow) ---
        if(rmin_[seg] > 0 || rmin_[seg + 1] > 0) {
            double rmin_lo = rmin_[seg];
            double rmin_hi = rmin_[seg + 1];

            double b_inner = (rmin_hi - rmin_lo) / dz_sec;
            double a_inner = rmin_lo - b_inner * z_lo;

            double A = dx * dx + dy * dy - b_inner * b_inner * dz * dz;
            double B = 2.0 * (px * dx + py * dy - b_inner * dz * (a_inner + b_inner * pz));
            double C = px * px + py * py - (a_inner + b_inner * pz) * (a_inner + b_inner * pz);

            if(!(dx == 0 && dy == 0 && b_inner == 0)) {
                double determinant = B * B - 4.0 * A * C;

                if(std::fabs(A) > GEOMETRY_PRECISION) {
                    if(determinant > 0) {
                        double sqrt_det = std::sqrt(determinant);
                        double t1 = (-B + sqrt_det) / (2.0 * A);
                        double t2 = (-B - sqrt_det) / (2.0 * A);

                        if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                            t1 = 0;
                        if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                            t2 = 0;

                        intersection_z = pz + t1 * dz;
                        if(intersection_z >= z_lo && intersection_z < z_hi) {
                            intersection_x = px + t1 * dx;
                            intersection_y = py + t1 * dy;
                            double r_at_z = a_inner + b_inner * intersection_z;
                            if(r_at_z >= 0) {
                                bool entering = !entering_cone(a_inner, b_inner);
                                emit(t1, entering);
                            }
                        }

                        intersection_z = pz + t2 * dz;
                        if(intersection_z >= z_lo && intersection_z < z_hi) {
                            intersection_x = px + t2 * dx;
                            intersection_y = py + t2 * dy;
                            double r_at_z = a_inner + b_inner * intersection_z;
                            if(r_at_z >= 0) {
                                bool entering = !entering_cone(a_inner, b_inner);
                                emit(t2, entering);
                            }
                        }
                    }
                } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                    double t1 = -C / B;
                    if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                        t1 = 0;

                    intersection_z = pz + t1 * dz;
                    if(intersection_z >= z_lo && intersection_z < z_hi) {
                        intersection_x = px + t1 * dx;
                        intersection_y = py + t1 * dy;
                        double r_at_z = a_inner + b_inner * intersection_z;
                        if(r_at_z >= 0) {
                            bool entering = !entering_cone(a_inner, b_inner);
                            emit(t1, entering);
                        }
                    }
                }
            }
        }
    }

    // --- End caps: bottom (z_planes_[0]) and top (z_planes_[n-1]) ---
    if(dz != 0) {
        // Bottom cap
        {
            double z_cap = z_planes_.front();
            double t = (z_cap - pz) / dz;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            intersection_x = px + t * dx;
            intersection_y = py + t * dy;

            double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;
            double rmax_cap = rmax_.front();
            double rmin_cap = rmin_.front();
            if(r2_hit <= rmax_cap * rmax_cap && r2_hit >= rmin_cap * rmin_cap) {
                intersection_z = pz + t * dz;
                emit(t, dz > 0);
            }
        }

        // Top cap
        {
            double z_cap = z_planes_.back();
            double t = (z_cap - pz) / dz;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            intersection_x = px + t * dx;
            intersection_y = py + t * dy;

            double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;
            double rmax_cap = rmax_.back();
            double rmin_cap = rmin_.back();
            if(r2_hit <= rmax_cap * rmax_cap && r2_hit >= rmin_cap * rmin_cap) {
                intersection_z = pz + t * dz;
                emit(t, dz < 0);
            }
        }
    }

    // --- Internal caps for step changes in radius ---
    // These arise at z-planes where consecutive entries share the same z
    // but have different radii (step discontinuity). We group consecutive
    // z-planes at the same z-value and compare radii across the group.
    if(dz != 0) {
        size_t i = 1; // skip first z-plane (endcap handles it)
        while(i + 1 < n) { // skip last z-plane
            double z_here = z_planes_[i];

            // Find the extent of z-planes at the same z
            size_t j = i;
            while(j + 1 < n - 1 && z_planes_[j + 1] == z_here) {
                j++;
            }

            // Radii at the top of the section below this group = values at index i
            // Radii at the bottom of the section above this group = values at index j
            double rmin_below = rmin_[i];
            double rmax_below = rmax_[i];
            double rmin_above = rmin_[j];
            double rmax_above = rmax_[j];

            if(rmax_below == rmax_above && rmin_below == rmin_above) {
                i = j + 1;
                continue;
            }

            double t = (z_here - pz) / dz;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            intersection_x = px + t * dx;
            intersection_y = py + t * dy;
            intersection_z = pz + t * dz;

            double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;

            // Outer annulus: region where rmax changed
            if(rmax_below != rmax_above) {
                double rlo = std::min(rmax_below, rmax_above);
                double rhi = std::max(rmax_below, rmax_above);
                if(r2_hit >= rlo * rlo && r2_hit <= rhi * rhi) {
                    // Material exists on the side with larger rmax
                    // entering = moving into material
                    bool entering = (rmax_above > rmax_below) ? (dz > 0) : (dz < 0);
                    emit(t, entering);
                }
            }

            // Inner annulus: region where rmin changed
            if(rmin_below != rmin_above) {
                double rlo = std::min(rmin_below, rmin_above);
                double rhi = std::max(rmin_below, rmin_above);
                if(r2_hit >= rlo * rlo && r2_hit <= rhi * rhi) {
                    // Hole grows on the side with larger rmin,
                    // so material exists on the side with smaller rmin
                    bool entering = (rmin_above > rmin_below) ? (dz < 0) : (dz > 0);
                    emit(t, entering);
                }
            }

            i = j + 1;
        }
    }

    if(n_hits == 0) return {};

    auto cmp = [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    };

    if(!has_phi_cut_) {
        // No phi cut: surface hits are the final result
        if(use_heap) {
            std::sort(heap_hits.begin(), heap_hits.end(), cmp);
            return heap_hits;
        }
        std::sort(stack_hits, stack_hits + n_hits, cmp);
        return {stack_hits, stack_hits + n_hits};
    }

    // Phi cut: merge surface hits with infinite wedge hits and run CSG walk.
    // Method: intersect the ray with the full-rotation solid AND an infinite
    // wedge (two half-planes from the z-axis). A sorted walk over both hit
    // lists produces the CSG intersection. This pattern is duplicated across
    // Polycone, GenericPolycone, Sphere, Torus, Cylinder, Cone, CutTube.
    // Collect surface hits into a sorted vector of tagged hits.
    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    // Sort surface hits first
    if(use_heap) {
        std::sort(heap_hits.begin(), heap_hits.end(), cmp);
    } else {
        std::sort(stack_hits, stack_hits + n_hits, cmp);
    }

    std::vector<TaggedHit> all_hits;
    all_hits.reserve(n_hits + 2);
    for(int i = 0; i < n_hits; ++i) {
        Intersection const & h = use_heap ? heap_hits[i] : stack_hits[i];
        all_hits.push_back({h.distance, h.position, h.entering, 0});
    }

    // Compute infinite wedge intersections (two half-planes from z-axis)
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

        double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
        // Must be on the correct half-plane (outward from z-axis)
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;

        bool entering = (n_dot_d < 0);
        all_hits.push_back({t, siren::math::Vector3D(hx, hy, hz), entering, 1});
    }

    if(all_hits.empty()) return {};

    // Sort all hits by distance
    std::sort(all_hits.begin(), all_hits.end(), [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // CSG intersection walk: the phi-cut solid is (full surface) AND (phi wedge).
    // Walk through sorted hits, tracking in_surface and in_wedge states.
    bool in_surface = false;
    bool in_wedge = false;

    // Determine initial in_wedge state
    bool has_wedge_hit = false;
    for(size_t i = 0; i < all_hits.size(); ++i) {
        if(all_hits[i].source == 1) {
            // First wedge hit: if entering, we started outside; if exiting, started inside
            in_wedge = !all_hits[i].entering;
            has_wedge_hit = true;
            break;
        }
    }
    if(!has_wedge_hit) {
        // No wedge crossings: ray is entirely in or entirely out of wedge.
        in_wedge = PhiInRange(px, py, start_phi_, delta_phi_);
    }

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
    double sp = NormalizePhi(start_phi_);
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
        if(cn < sp) cn += TWO_PI;
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
