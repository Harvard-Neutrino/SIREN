#include "SIREN/geometry/Polyhedra.h"

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

static const double TWO_PI = 2.0 * M_PI;

double NormalizePhi(double phi) {
    phi = std::fmod(phi, TWO_PI);
    if(phi < 0) phi += TWO_PI;
    return phi;
}

bool PhiInRange(double x, double y, double start_phi, double delta_phi) {
    double phi = NormalizePhi(std::atan2(y, x));
    double sp = NormalizePhi(start_phi);
    double ep = sp + delta_phi;
    if(ep <= TWO_PI + 1e-9) {
        return phi >= sp - 1e-9 && phi <= ep + 1e-9;
    } else {
        return phi >= sp - 1e-9 || phi <= NormalizePhi(ep) + 1e-9;
    }
}

} // anonymous namespace

void Polyhedra::validate() const {
    if(num_sides_ < 3) {
        throw std::runtime_error("Polyhedra requires at least 3 sides!");
    }
    if(num_sides_ > 64) {
        throw std::runtime_error("Polyhedra supports at most 64 sides!");
    }
    if(z_planes_.size() < 2) {
        throw std::runtime_error("Polyhedra requires at least 2 z-planes!");
    }
    if(z_planes_.size() != rmin_.size() || z_planes_.size() != rmax_.size()) {
        throw std::runtime_error("Polyhedra z_planes, rmin, and rmax vectors must have the same size!");
    }
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        if(rmin_[i] < 0 || rmax_[i] < 0) {
            throw std::runtime_error("Polyhedra radii must be non-negative!");
        }
        if(rmin_[i] > rmax_[i]) {
            throw std::runtime_error("Polyhedra inner radius must not exceed outer radius at each z-plane!");
        }
    }
    for(size_t i = 1; i < z_planes_.size(); ++i) {
        if(z_planes_[i] < z_planes_[i - 1]) {
            throw std::runtime_error("Polyhedra z-planes must be sorted in ascending order!");
        }
    }
}

void Polyhedra::precompute_trig() {
    if(num_sides_ < 3) {
        cos_phi_.clear();
        sin_phi_.clear();
        return;
    }
    double dphi = 2.0 * M_PI / num_sides_;
    cos_phi_.resize(num_sides_ + 1);
    sin_phi_.resize(num_sides_ + 1);
    for(int k = 0; k < num_sides_; ++k) {
        double angle = start_phi_ + k * dphi;
        cos_phi_[k] = std::cos(angle);
        sin_phi_[k] = std::sin(angle);
    }
    // Wrap-around entry so index [num_sides_] == index [0]
    cos_phi_[num_sides_] = cos_phi_[0];
    sin_phi_[num_sides_] = sin_phi_[0];
}

Polyhedra::Polyhedra()
    : Geometry((std::string)("Polyhedra"))
    , num_sides_(0)
    , start_phi_(0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false)
    , z_planes_()
    , rmin_()
      , rmax_() {
    // Do nothing here
}

Polyhedra::Polyhedra(int num_sides,
                     double start_phi,
                     std::vector<double> const & z_planes,
                     std::vector<double> const & rmin,
                     std::vector<double> const & rmax)
    : Geometry((std::string)("Polyhedra"))
    , num_sides_(num_sides)
    , start_phi_(start_phi)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false)
    , z_planes_(z_planes)
    , rmin_(rmin)
      , rmax_(rmax) {
    validate();
    precompute_trig();
}

Polyhedra::Polyhedra(Placement const & placement)
    : Geometry((std::string)("Polyhedra"), placement)
    , num_sides_(0)
    , start_phi_(0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false)
    , z_planes_()
    , rmin_()
      , rmax_() {
    // Do nothing here
}

Polyhedra::Polyhedra(Placement const & placement,
                     int num_sides,
                     double start_phi,
                     std::vector<double> const & z_planes,
                     std::vector<double> const & rmin,
                     std::vector<double> const & rmax)
    : Geometry((std::string)("Polyhedra"), placement)
    , num_sides_(num_sides)
    , start_phi_(start_phi)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false)
    , z_planes_(z_planes)
    , rmin_(rmin)
      , rmax_(rmax) {
    validate();
    precompute_trig();
}

Polyhedra::Polyhedra(int num_sides,
                     double start_phi,
                     std::vector<double> const & z_planes,
                     std::vector<double> const & rmin,
                     std::vector<double> const & rmax,
                     double delta_phi)
    : Geometry((std::string)("Polyhedra"))
    , num_sides_(num_sides)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi)
    , has_phi_cut_(std::fabs(delta_phi - 2.0 * M_PI) > 1e-9)
    , z_planes_(z_planes)
    , rmin_(rmin)
      , rmax_(rmax) {
    validate();
    precompute_trig();
}

Polyhedra::Polyhedra(Placement const & placement,
                     int num_sides,
                     double start_phi,
                     std::vector<double> const & z_planes,
                     std::vector<double> const & rmin,
                     std::vector<double> const & rmax,
                     double delta_phi)
    : Geometry((std::string)("Polyhedra"), placement)
    , num_sides_(num_sides)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi)
    , has_phi_cut_(std::fabs(delta_phi - 2.0 * M_PI) > 1e-9)
    , z_planes_(z_planes)
    , rmin_(rmin)
      , rmax_(rmax) {
    validate();
    precompute_trig();
}

Polyhedra::Polyhedra(const Polyhedra& polyhedra)
    : Geometry(polyhedra)
    , num_sides_(polyhedra.num_sides_)
    , start_phi_(polyhedra.start_phi_)
    , delta_phi_(polyhedra.delta_phi_)
    , has_phi_cut_(polyhedra.has_phi_cut_)
    , z_planes_(polyhedra.z_planes_)
    , rmin_(polyhedra.rmin_)
    , rmax_(polyhedra.rmax_)
    , cos_phi_(polyhedra.cos_phi_)
      , sin_phi_(polyhedra.sin_phi_) {
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Polyhedra::swap(Geometry& geometry) {
    Polyhedra* polyhedra = dynamic_cast<Polyhedra*>(&geometry);
    if(!polyhedra) {
        return;
    }

    Geometry::swap(*polyhedra);

    std::swap(num_sides_, polyhedra->num_sides_);
    std::swap(start_phi_, polyhedra->start_phi_);
    std::swap(delta_phi_, polyhedra->delta_phi_);
    std::swap(has_phi_cut_, polyhedra->has_phi_cut_);
    std::swap(z_planes_, polyhedra->z_planes_);
    std::swap(rmin_, polyhedra->rmin_);
    std::swap(rmax_, polyhedra->rmax_);
    std::swap(cos_phi_, polyhedra->cos_phi_);
    std::swap(sin_phi_, polyhedra->sin_phi_);
}

// ------------------------------------------------------------------------- //
Polyhedra& Polyhedra::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Polyhedra* polyhedra = dynamic_cast<const Polyhedra*>(&geometry);
        if(!polyhedra) {
            return *this;
        }
        Polyhedra tmp(*polyhedra);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Polyhedra::equal(const Geometry& geometry) const
{
    const Polyhedra* polyhedra = dynamic_cast<const Polyhedra*>(&geometry);

    if(!polyhedra)
        return false;
    else if(num_sides_ != polyhedra->num_sides_)
        return false;
    else if(start_phi_ != polyhedra->start_phi_)
        return false;
    else if(delta_phi_ != polyhedra->delta_phi_)
        return false;
    else if(z_planes_ != polyhedra->z_planes_)
        return false;
    else if(rmin_ != polyhedra->rmin_)
        return false;
    else if(rmax_ != polyhedra->rmax_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Polyhedra::less(const Geometry& geometry) const
{
    const Polyhedra* polyhedra = dynamic_cast<const Polyhedra*>(&geometry);
    if(!polyhedra) return false;

    return
        std::tie(num_sides_, start_phi_, delta_phi_, z_planes_, rmin_, rmax_)
        <
        std::tie(polyhedra->num_sides_, polyhedra->start_phi_, polyhedra->delta_phi_, polyhedra->z_planes_, polyhedra->rmin_, polyhedra->rmax_);
}

void Polyhedra::print(std::ostream& os) const
{
    os << "Polyhedra with " << num_sides_ << " sides, "
       << z_planes_.size() << " z-planes, start_phi=" << start_phi_
       << ", delta_phi=" << delta_phi_ << ":";
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        os << " [z=" << z_planes_[i]
           << " rmin=" << rmin_[i]
           << " rmax=" << rmax_[i] << "]";
    }
    os << '\n';
}

// ------------------------------------------------------------------------- //
// ComputeIntersections
//
// The polyhedra is decomposed into sections between adjacent z-planes.
// For each section, each polygon side produces a planar quadrilateral face
// (outer and optionally inner). We also check the top/bottom end caps
// as annular polygonal faces.
//
// All coordinates are in local frame. Distance can be negative (full-line
// intersections, same convention as other SIREN geometry classes).
// ------------------------------------------------------------------------- //

// Helper: check if a point (px, py) lies inside a convex polygon defined
// by vertices (vx[], vy[]) with nv vertices, listed in order.
// Uses the cross-product sign test for convex polygons.
// Tolerance eps matches PointInConvexFace3D for consistent boundary handling.
static bool PointInConvexPolygon(double px, double py,
                                 const double* vx, const double* vy, int nv) {
    static constexpr double eps = Geometry::GEOMETRY_PRECISION;
    bool has_pos = false;
    bool has_neg = false;
    for(int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        double ex = vx[j] - vx[i];
        double ey = vy[j] - vy[i];
        double cross = ex * (py - vy[i]) - ey * (px - vx[i]);
        if(cross > eps) has_pos = true;
        if(cross < -eps) has_neg = true;
        if(has_pos && has_neg) return false;
    }
    return true;
}

// Helper: check if a point (px, py) lies inside the annular polygon region
// (between inner polygon at rmin and outer polygon at rmax).
// Uses precomputed cos_phi/sin_phi arrays (num_sides entries) to avoid trig.
static bool PointInAnnularPolygon(double px, double py,
                                   int num_sides,
                                   const double* cos_phi, const double* sin_phi,
                                   double rmin, double rmax) {
    // Outer polygon vertices
    double ovx[64], ovy[64]; // max 64 sides
    for(int k = 0; k < num_sides; ++k) {
        ovx[k] = rmax * cos_phi[k];
        ovy[k] = rmax * sin_phi[k];
    }
    if(!PointInConvexPolygon(px, py, ovx, ovy, num_sides)) {
        return false;
    }

    // Inner polygon check (if hollow)
    if(rmin > 0) {
        double ivx[64], ivy[64];
        for(int k = 0; k < num_sides; ++k) {
            ivx[k] = rmin * cos_phi[k];
            ivy[k] = rmin * sin_phi[k];
        }
        if(PointInConvexPolygon(px, py, ivx, ivy, num_sides)) {
            return false; // inside inner hole
        }
    }

    return true;
}

static bool PointInConvexFace3D(double qx, double qy, double qz,
                                double nx, double ny, double nz,
                                const double* vx, const double* vy, const double* vz,
                                int nv, double eps) {
    bool has_pos = false;
    bool has_neg = false;
    for(int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        double ex = vx[j] - vx[i], ey = vy[j] - vy[i], ez = vz[j] - vz[i];
        double fx = qx - vx[i], fy = qy - vy[i], fz = qz - vz[i];
        double c = nx * (ey*fz - ez*fy) + ny * (ez*fx - ex*fz) + nz * (ex*fy - ey*fx);
        if(c > eps) has_pos = true;
        if(c < -eps) has_neg = true;
        if(has_pos && has_neg) return false;
    }
    return true;
}

std::vector<Geometry::Intersection> Polyhedra::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {

    if(z_planes_.size() < 2 || num_sides_ < 3) {
        return {};
    }

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    // Stack buffer for typical cases; heap fallback for complex polyhedra.
    static constexpr int STACK_CAPACITY = 64;
    Intersection stack_hits[STACK_CAPACITY];
    std::vector<Intersection> heap_hits;
    int n_hits = 0;
    bool using_heap = false;

    auto add_hit = [&](double dist, int hierarchy, bool entering,
                       double hx, double hy, double hz) {
        Intersection isect;
        isect.distance = dist;
        isect.hierarchy = hierarchy;
        isect.entering = entering;
        isect.position = siren::math::Vector3D(hx, hy, hz);
        if(!using_heap) {
            if(n_hits < STACK_CAPACITY) {
                stack_hits[n_hits++] = isect;
            } else {
                using_heap = true;
                heap_hits.assign(stack_hits, stack_hits + STACK_CAPACITY);
                heap_hits.push_back(isect);
                n_hits++;
            }
        } else {
            heap_hits.push_back(isect);
            n_hits++;
        }
    };

    double intersection_x;
    double intersection_y;
    double intersection_z;

    size_t n = z_planes_.size();

    // Precomputed cos/sin arrays: cos_phi_[k], sin_phi_[k] for k = 0..num_sides_
    // with wrap-around at index [num_sides_] == [0].
    const double* cp = cos_phi_.data();
    const double* sp = sin_phi_.data();

    // ---- Lateral faces (outer and inner) for each z-section and each side ----
    for(size_t seg = 0; seg + 1 < n; ++seg) {
        double z_lo = z_planes_[seg];
        double z_hi = z_planes_[seg + 1];

        if(z_hi <= z_lo) {
            continue;
        }

        // Horizontal ray cannot hit lateral faces outside its z-range.
        // Use half-open interval [z_lo, z_hi) so each boundary belongs to
        // exactly one section and horizontal rays at internal z-planes are
        // not skipped by both adjacent sections.
        if(std::fabs(dz) < GEOMETRY_PRECISION && (pz < z_lo || pz >= z_hi)) {
            continue;
        }

        double rmax_lo = rmax_[seg];
        double rmax_hi = rmax_[seg + 1];
        double rmin_lo = rmin_[seg];
        double rmin_hi = rmin_[seg + 1];

        for(int side = 0; side < num_sides_; ++side) {
            double cos_k = cp[side];
            double sin_k = sp[side];
            double cos_k1 = cp[side + 1];
            double sin_k1 = sp[side + 1];

            // ---- Outer lateral face ----
            // The four vertices of this quad face (going around the quad):
            //   v0 = (rmax_lo * cos_k,  rmax_lo * sin_k,  z_lo)
            //   v1 = (rmax_lo * cos_k1, rmax_lo * sin_k1, z_lo)
            //   v2 = (rmax_hi * cos_k1, rmax_hi * sin_k1, z_hi)
            //   v3 = (rmax_hi * cos_k,  rmax_hi * sin_k,  z_hi)
            // When one radius is zero, v0==v1 or v2==v3, making the quad
            // degenerate to a triangle. Use a non-degenerate edge pair
            // for the normal computation.
            {
                double v0x = rmax_lo * cos_k,  v0y = rmax_lo * sin_k,  v0z = z_lo;
                double v1x = rmax_lo * cos_k1, v1y = rmax_lo * sin_k1, v1z = z_lo;
                double v3x = rmax_hi * cos_k,  v3y = rmax_hi * sin_k,  v3z = z_hi;

                double e1x, e1y, e1z, e2x, e2y, e2z;
                if(rmax_lo < GEOMETRY_PRECISION) {
                    double v2x = rmax_hi * cos_k1, v2y = rmax_hi * sin_k1;
                    e1x = v3x - v2x; e1y = v3y - v2y; e1z = 0;
                    e2x = v0x - v3x; e2y = v0y - v3y; e2z = v0z - v3z;
                } else {
                    e1x = v1x - v0x; e1y = v1y - v0y; e1z = v1z - v0z;
                    e2x = v3x - v0x; e2y = v3y - v0y; e2z = v3z - v0z;
                }

                // Normal = e1 x e2 (points outward for outer face)
                double nx = e1y * e2z - e1z * e2y;
                double ny = e1z * e2x - e1x * e2z;
                double nz = e1x * e2y - e1y * e2x;

                double n_dot_d = nx * dx + ny * dy + nz * dz;

                if(std::fabs(n_dot_d) > GEOMETRY_PRECISION) {
                    // Ray-plane intersection: t = n . (v0 - p) / (n . d)
                    double n_dot_pv = nx * (v0x - px) + ny * (v0y - py) + nz * (v0z - pz);
                    double t = n_dot_pv / n_dot_d;

                    if(t > 0 && t < GEOMETRY_PRECISION)
                        t = 0;

                    intersection_x = px + t * dx;
                    intersection_y = py + t * dy;
                    intersection_z = pz + t * dz;

                    // Check if hit point is within the z-range
                    if(intersection_z >= z_lo - GEOMETRY_PRECISION &&
                        intersection_z <= z_hi + GEOMETRY_PRECISION) {

                        double v2x = rmax_hi * cos_k1, v2y = rmax_hi * sin_k1, v2z = z_hi;
                        double qx = intersection_x, qy = intersection_y, qz = intersection_z;

                        bool inside_face;
                        if(rmax_lo < GEOMETRY_PRECISION) {
                            double fvx[] = {v0x, v2x, v3x};
                            double fvy[] = {v0y, v2y, v3y};
                            double fvz[] = {v0z, v2z, v3z};
                            inside_face = PointInConvexFace3D(qx, qy, qz, nx, ny, nz, fvx, fvy, fvz, 3, GEOMETRY_PRECISION);
                        } else if(rmax_hi < GEOMETRY_PRECISION) {
                            double fvx[] = {v0x, v1x, v2x};
                            double fvy[] = {v0y, v1y, v2y};
                            double fvz[] = {v0z, v1z, v2z};
                            inside_face = PointInConvexFace3D(qx, qy, qz, nx, ny, nz, fvx, fvy, fvz, 3, GEOMETRY_PRECISION);
                        } else {
                            double fvx[] = {v0x, v1x, v2x, v3x};
                            double fvy[] = {v0y, v1y, v2y, v3y};
                            double fvz[] = {v0z, v1z, v2z, v3z};
                            inside_face = PointInConvexFace3D(qx, qy, qz, nx, ny, nz, fvx, fvy, fvz, 4, GEOMETRY_PRECISION);
                        }

                        if(inside_face) {
                            add_hit(t, 0, (n_dot_d < 0), intersection_x, intersection_y, intersection_z);
                        }
                    }
                }
            }

            // ---- Inner lateral face (if hollow) ----
            if(rmin_lo > 0 || rmin_hi > 0) {
                double v0x = rmin_lo * cos_k,  v0y = rmin_lo * sin_k,  v0z = z_lo;
                double v1x = rmin_lo * cos_k1, v1y = rmin_lo * sin_k1, v1z = z_lo;
                double v3x = rmin_hi * cos_k,  v3y = rmin_hi * sin_k,  v3z = z_hi;

                double e1x, e1y, e1z, e2x, e2y, e2z;
                if(rmin_lo < GEOMETRY_PRECISION) {
                    double v2x = rmin_hi * cos_k1, v2y = rmin_hi * sin_k1;
                    e1x = v3x - v2x; e1y = v3y - v2y; e1z = 0;
                    e2x = v0x - v3x; e2y = v0y - v3y; e2z = v0z - v3z;
                } else {
                    e1x = v1x - v0x; e1y = v1y - v0y; e1z = v1z - v0z;
                    e2x = v3x - v0x; e2y = v3y - v0y; e2z = v3z - v0z;
                }

                double nx = e1y * e2z - e1z * e2y;
                double ny = e1z * e2x - e1x * e2z;
                double nz = e1x * e2y - e1y * e2x;

                double n_dot_d = nx * dx + ny * dy + nz * dz;

                if(std::fabs(n_dot_d) > GEOMETRY_PRECISION) {
                    double n_dot_pv = nx * (v0x - px) + ny * (v0y - py) + nz * (v0z - pz);
                    double t = n_dot_pv / n_dot_d;

                    if(t > 0 && t < GEOMETRY_PRECISION)
                        t = 0;

                    intersection_x = px + t * dx;
                    intersection_y = py + t * dy;
                    intersection_z = pz + t * dz;

                    if(intersection_z >= z_lo - GEOMETRY_PRECISION &&
                        intersection_z <= z_hi + GEOMETRY_PRECISION) {

                        double v2x = rmin_hi * cos_k1, v2y = rmin_hi * sin_k1, v2z = z_hi;
                        double qx = intersection_x, qy = intersection_y, qz = intersection_z;

                        bool inside_face;
                        if(rmin_lo < GEOMETRY_PRECISION) {
                            double fvx[] = {v0x, v2x, v3x};
                            double fvy[] = {v0y, v2y, v3y};
                            double fvz[] = {v0z, v2z, v3z};
                            inside_face = PointInConvexFace3D(qx, qy, qz, nx, ny, nz, fvx, fvy, fvz, 3, GEOMETRY_PRECISION);
                        } else if(rmin_hi < GEOMETRY_PRECISION) {
                            double fvx[] = {v0x, v1x, v2x};
                            double fvy[] = {v0y, v1y, v2y};
                            double fvz[] = {v0z, v1z, v2z};
                            inside_face = PointInConvexFace3D(qx, qy, qz, nx, ny, nz, fvx, fvy, fvz, 3, GEOMETRY_PRECISION);
                        } else {
                            double fvx[] = {v0x, v1x, v2x, v3x};
                            double fvy[] = {v0y, v1y, v2y, v3y};
                            double fvz[] = {v0z, v1z, v2z, v3z};
                            inside_face = PointInConvexFace3D(qx, qy, qz, nx, ny, nz, fvx, fvy, fvz, 4, GEOMETRY_PRECISION);
                        }

                        if(inside_face) {
                            add_hit(t, 0, (n_dot_d > 0), intersection_x, intersection_y, intersection_z);
                        }
                    }
                }
            }
        }
    }

    // ---- End caps: bottom (z_planes_[0]) and top (z_planes_[n-1]) ----
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
        // Bottom cap
        {
            double z_cap = z_planes_.front();
            double t = (z_cap - pz) / dz;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            intersection_x = px + t * dx;
            intersection_y = py + t * dy;
            intersection_z = pz + t * dz;

            if(rmax_.front() >= GEOMETRY_PRECISION &&
               PointInAnnularPolygon(intersection_x, intersection_y,
                                       num_sides_, cp, sp,
                                       rmin_.front(), rmax_.front())) {
                add_hit(t, 0, (dz > 0), intersection_x, intersection_y, intersection_z);
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
            intersection_z = pz + t * dz;

            if(rmax_.back() >= GEOMETRY_PRECISION &&
               PointInAnnularPolygon(intersection_x, intersection_y,
                                       num_sides_, cp, sp,
                                       rmin_.back(), rmax_.back())) {
                add_hit(t, 0, (dz < 0), intersection_x, intersection_y, intersection_z);
            }
        }
    }

    // ---- Internal caps for step changes in radius ----
    // These arise at z-planes where consecutive entries share the same z
    // but have different radii (step discontinuity). We group consecutive
    // z-planes at the same z-value and compare radii across the group.
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
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

            // Outer annulus: region where rmax changed
            if(rmax_below != rmax_above) {
                double rlo = std::min(rmax_below, rmax_above);
                double rhi = std::max(rmax_below, rmax_above);
                if(rhi >= GEOMETRY_PRECISION &&
                   PointInAnnularPolygon(intersection_x, intersection_y,
                                           num_sides_, cp, sp, rlo, rhi)) {
                    bool entering = (rmax_above > rmax_below) ? (dz > 0) : (dz < 0);
                    add_hit(t, 0, entering, intersection_x, intersection_y, intersection_z);
                }
            }

            // Inner annulus: region where rmin changed
            if(rmin_below != rmin_above) {
                double rlo = std::min(rmin_below, rmin_above);
                double rhi = std::max(rmin_below, rmin_above);
                if(rhi >= GEOMETRY_PRECISION &&
                   PointInAnnularPolygon(intersection_x, intersection_y,
                                           num_sides_, cp, sp, rlo, rhi)) {
                    bool entering = (rmin_above > rmin_below) ? (dz < 0) : (dz > 0);
                    add_hit(t, 0, entering, intersection_x, intersection_y, intersection_z);
                }
            }

            i = j + 1;
        }
    }

    // Sort hits by distance
    auto cmp = [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    };

    if(!has_phi_cut_) {
        if(using_heap) {
            std::sort(heap_hits.begin(), heap_hits.end(), cmp);
            return heap_hits;
        }
        std::sort(stack_hits, stack_hits + n_hits, cmp);
        return {stack_hits, stack_hits + n_hits};
    }

    // Phi cut: CSG intersection of the full-rotation solid with an infinite wedge.
    if(using_heap) {
        std::sort(heap_hits.begin(), heap_hits.end(), cmp);
    } else {
        std::sort(stack_hits, stack_hits + n_hits, cmp);
    }

    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    std::vector<TaggedHit> all_hits;
    all_hits.reserve(n_hits + 2);
    for(int i = 0; i < n_hits; ++i) {
        Intersection const & h = using_heap ? heap_hits[i] : stack_hits[i];
        all_hits.push_back({h.distance, h.position, h.entering, 0});
    }

    // Infinite wedge: two half-planes from the z-axis
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
        all_hits.push_back({t, siren::math::Vector3D(hx, hy, hz), entering, 1});
    }

    if(all_hits.empty()) return {};

    std::sort(all_hits.begin(), all_hits.end(), [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // CSG walk
    bool in_surface = false;
    bool in_wedge = false;

    bool has_wedge_hit = false;
    for(size_t i = 0; i < all_hits.size(); ++i) {
        if(all_hits[i].source == 1) {
            in_wedge = !all_hits[i].entering;
            has_wedge_hit = true;
            break;
        }
    }
    if(!has_wedge_hit) {
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
AABB Polyhedra::GetBoundingBox() const {
    if(z_planes_.empty()) {
        return AABB(math::Vector3D(0, 0, 0), math::Vector3D(0, 0, 0));
    }
    double max_r = *std::max_element(rmax_.begin(), rmax_.end());
    double circum_r = max_r / std::cos(M_PI / num_sides_);
    double z_lo = z_planes_.front();
    double z_hi = z_planes_.back();

    if(!has_phi_cut_) {
        return AABB(
            math::Vector3D(-circum_r, -circum_r, z_lo),
            math::Vector3D( circum_r,  circum_r, z_hi)
        );
    }

    double sp = NormalizePhi(start_phi_);
    double ep = sp + delta_phi_;

    double x_min = 0, x_max = 0, y_min = 0, y_max = 0;

    double angles[2] = {sp, sp + delta_phi_};
    for(int i = 0; i < 2; ++i) {
        double a = angles[i];
        double cx = std::cos(a) * circum_r;
        double cy = std::sin(a) * circum_r;
        if(cx < x_min) x_min = cx;
        if(cx > x_max) x_max = cx;
        if(cy < y_min) y_min = cy;
        if(cy > y_max) y_max = cy;
    }

    double cardinals[4] = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    for(int i = 0; i < 4; ++i) {
        double c = cardinals[i];
        double cn = c;
        if(cn < sp) cn += TWO_PI;
        if(cn <= ep + 1e-9) {
            double cx = std::cos(cardinals[i]) * circum_r;
            double cy = std::sin(cardinals[i]) * circum_r;
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
