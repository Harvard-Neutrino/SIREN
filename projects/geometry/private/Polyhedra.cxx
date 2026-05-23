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
#include "GeometryMacros.h"
#include "PhiUtils.h"

static constexpr int POLYHEDRA_MAX_SIDES = 64;

namespace siren {
namespace geometry {

void Polyhedra::validate() {
    if(delta_phi_ <= 0) {
        throw std::invalid_argument("Polyhedra delta_phi must be positive!");
    }
    if(delta_phi_ > 2.0 * M_PI) {
        delta_phi_ = 2.0 * M_PI;
        has_phi_cut_ = false;
    }
    if(num_sides_ < 3) {
        throw std::runtime_error("Polyhedra requires at least 3 sides!");
    }
    if(num_sides_ > POLYHEDRA_MAX_SIDES) {
        throw std::runtime_error("Polyhedra supports at most 64 sides!");
    }
    if(z_planes_.size() < 2) {
        throw std::runtime_error("Polyhedra requires at least 2 z-planes!");
    }
    if(z_planes_.size() != rmin_.size() || z_planes_.size() != rmax_.size()) {
        throw std::runtime_error("Polyhedra z_planes, rmin, and rmax vectors must have the same size!");
    }
    // Accept ascending or descending z-order (including duplicate z-values).
    // If descending, reverse all three vectors so the geometry math sees
    // ascending order. Non-monotonic sequences are rejected.
    // This must happen before the rmin/rmax check since the reversal
    // reorders the radii arrays.
    bool ascending = true, descending = true;
    for(size_t i = 1; i < z_planes_.size(); ++i) {
        if(z_planes_[i] < z_planes_[i - 1]) ascending = false;
        if(z_planes_[i] > z_planes_[i - 1]) descending = false;
    }
    if(!ascending && !descending) {
        throw std::runtime_error("Polyhedra z-planes must be monotonic (ascending or descending)!");
    }
    if(descending && !ascending) {
        std::reverse(z_planes_.begin(), z_planes_.end());
        std::reverse(rmin_.begin(), rmin_.end());
        std::reverse(rmax_.begin(), rmax_.end());
    }
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        if(rmin_[i] < 0 || rmax_[i] < 0) {
            throw std::runtime_error("Polyhedra radii must be non-negative!");
        }
        if(rmin_[i] > rmax_[i]) {
            throw std::runtime_error("Polyhedra inner radius must not exceed outer radius at each z-plane!");
        }
    }
}

void Polyhedra::precompute_trig() {
    if(num_sides_ < 3) {
        cos_phi_.clear();
        sin_phi_.clear();
        return;
    }
    double dphi = has_phi_cut_ ? delta_phi_ / num_sides_ : 2.0 * M_PI / num_sides_;
    cos_phi_.resize(num_sides_ + 1);
    sin_phi_.resize(num_sides_ + 1);
    for(int k = 0; k < num_sides_; ++k) {
        double angle = start_phi_ + k * dphi;
        cos_phi_[k] = std::cos(angle);
        sin_phi_[k] = std::sin(angle);
    }
    if(has_phi_cut_) {
        double last_angle = start_phi_ + delta_phi_;
        cos_phi_[num_sides_] = std::cos(last_angle);
        sin_phi_[num_sides_] = std::sin(last_angle);
    } else {
        cos_phi_[num_sides_] = cos_phi_[0];
        sin_phi_[num_sides_] = sin_phi_[0];
    }
}

Polyhedra::Polyhedra() : Geometry("Polyhedra"), num_sides_(0), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Polyhedra::Polyhedra(int num_sides, double start_phi, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax) : Geometry("Polyhedra"), num_sides_(num_sides), start_phi_(start_phi), delta_phi_(2.0 * M_PI), has_phi_cut_(false), z_planes_(z_planes), rmin_(rmin), rmax_(rmax) { validate(); precompute_trig(); RecomputeWorldAABB(); }
Polyhedra::Polyhedra(Placement const & p) : Geometry("Polyhedra", p), num_sides_(0), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Polyhedra::Polyhedra(Placement const & p, int num_sides, double start_phi, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax) : Geometry("Polyhedra", p), num_sides_(num_sides), start_phi_(start_phi), delta_phi_(2.0 * M_PI), has_phi_cut_(false), z_planes_(z_planes), rmin_(rmin), rmax_(rmax) { validate(); precompute_trig(); RecomputeWorldAABB(); }
Polyhedra::Polyhedra(int num_sides, double start_phi, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax, double delta_phi) : Geometry("Polyhedra"), num_sides_(num_sides), start_phi_(start_phi), delta_phi_(delta_phi), has_phi_cut_(std::fabs(delta_phi - 2.0 * M_PI) > 1e-9), z_planes_(z_planes), rmin_(rmin), rmax_(rmax) { validate(); precompute_trig(); RecomputeWorldAABB(); }
Polyhedra::Polyhedra(Placement const & p, int num_sides, double start_phi, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax, double delta_phi) : Geometry("Polyhedra", p), num_sides_(num_sides), start_phi_(start_phi), delta_phi_(delta_phi), has_phi_cut_(std::fabs(delta_phi - 2.0 * M_PI) > 1e-9), z_planes_(z_planes), rmin_(rmin), rmax_(rmax) { validate(); precompute_trig(); RecomputeWorldAABB(); }
Polyhedra::Polyhedra(const Polyhedra& o) : Geometry(o), num_sides_(o.num_sides_), start_phi_(o.start_phi_), delta_phi_(o.delta_phi_), has_phi_cut_(o.has_phi_cut_), z_planes_(o.z_planes_), rmin_(o.rmin_), rmax_(o.rmax_), cos_phi_(o.cos_phi_), sin_phi_(o.sin_phi_) { RecomputeWorldAABB(); }

// Swap includes derived fields; equal/less use only primary fields
SIREN_GEOMETRY_SWAP(Polyhedra, num_sides_, start_phi_, delta_phi_, has_phi_cut_, z_planes_, rmin_, rmax_, cos_phi_, sin_phi_)
SIREN_GEOMETRY_ASSIGN(Polyhedra)
SIREN_GEOMETRY_EQUAL(Polyhedra, num_sides_, start_phi_, delta_phi_, z_planes_, rmin_, rmax_)
SIREN_GEOMETRY_LESS(Polyhedra, num_sides_, start_phi_, delta_phi_, z_planes_, rmin_, rmax_)

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
    double ovx[POLYHEDRA_MAX_SIDES], ovy[POLYHEDRA_MAX_SIDES];
    for(int k = 0; k < num_sides; ++k) {
        ovx[k] = rmax * cos_phi[k];
        ovy[k] = rmax * sin_phi[k];
    }
    if(!PointInConvexPolygon(px, py, ovx, ovy, num_sides)) {
        return false;
    }

    // Inner polygon check (if hollow)
    if(rmin > 0) {
        double ivx[POLYHEDRA_MAX_SIDES], ivy[POLYHEDRA_MAX_SIDES];
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

// Ray-casting point-in-polygon for simple (possibly non-convex) polygons.
// Needed for partial-phi annular sector end caps where the closing edge
// makes the polygon non-convex when delta_phi > pi.
static bool PointInSimplePolygon(double px, double py,
                                 const double* vx, const double* vy, int nv) {
    static constexpr double eps = Geometry::GEOMETRY_PRECISION;
    for(int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        double ex = vx[j] - vx[i], ey = vy[j] - vy[i];
        double fx = px - vx[i], fy = py - vy[i];
        double edge_len2 = ex*ex + ey*ey;
        if(edge_len2 < 1e-30) continue;
        double cross = ex*fy - ey*fx;
        if(cross*cross <= eps*eps * edge_len2) {
            double dot = fx*ex + fy*ey;
            if(dot >= -eps && dot <= edge_len2 + eps) return true;
        }
    }
    bool inside = false;
    for(int i = 0, j = nv - 1; i < nv; j = i++) {
        if(((vy[i] > py) != (vy[j] > py)) &&
           (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])) {
            inside = !inside;
        }
    }
    return inside;
}

// Check if (px, py) lies inside an annular sector polygon for partial-phi shapes.
// The sector is bounded by outer arc edges, inner arc edges (reversed), and two
// radial closing edges connecting the arcs at start_phi and end_phi.
// cos_phi/sin_phi must have num_sides+1 entries (including the endpoint).
// rmin/rmax are circumscribed (vertex) radii.
static bool PointInAnnularSectorPolygon(double px, double py,
                                        int num_sides,
                                        const double* cos_phi, const double* sin_phi,
                                        double rmin, double rmax) {
    static constexpr double eps = Geometry::GEOMETRY_PRECISION;
    if(rmin < eps) {
        int nv = num_sides + 2;
        double vx[POLYHEDRA_MAX_SIDES + 2], vy[POLYHEDRA_MAX_SIDES + 2];
        vx[0] = 0.0; vy[0] = 0.0;
        for(int k = 0; k <= num_sides; ++k) {
            vx[k + 1] = rmax * cos_phi[k];
            vy[k + 1] = rmax * sin_phi[k];
        }
        return PointInSimplePolygon(px, py, vx, vy, nv);
    }
    int nv = 2 * (num_sides + 1);
    double vx[2 * (POLYHEDRA_MAX_SIDES + 1)], vy[2 * (POLYHEDRA_MAX_SIDES + 1)];
    for(int k = 0; k <= num_sides; ++k) {
        vx[k] = rmax * cos_phi[k];
        vy[k] = rmax * sin_phi[k];
    }
    for(int k = 0; k <= num_sides; ++k) {
        int idx = num_sides + 1 + k;
        int src = num_sides - k;
        vx[idx] = rmin * cos_phi[src];
        vy[idx] = rmin * sin_phi[src];
    }
    return PointInSimplePolygon(px, py, vx, vy, nv);
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

    // Early-out: check if the ray can possibly intersect the bounding cylinder.
    // The circumscribed radius bounds all vertices; any intersection must be
    // within this radius in the xy-plane.
    double max_cr = *std::max_element(rmax_.begin(), rmax_.end());
    double face_dphi = has_phi_cut_ ? delta_phi_ / num_sides_ : 2.0 * M_PI / num_sides_;
    double bounding_r = max_cr / std::cos(0.5 * face_dphi);
    double bounding_r2 = bounding_r * bounding_r;
    {
        // Closest approach of ray to z-axis in xy-plane
        double C = dx*dx + dy*dy;
        if(C > 0) {
            double cross_z = px*dy - py*dx;
            if(cross_z * cross_z > C * bounding_r2) return {};
        } else {
            // Ray parallel to z-axis
            if(px*px + py*py > bounding_r2) return {};
        }
        // Z-range check
        double z_lo = z_planes_.front(), z_hi = z_planes_.back();
        if(std::fabs(dz) < GEOMETRY_PRECISION) {
            if(pz < z_lo || pz > z_hi) return {};
        } else {
            double t_lo = (z_lo - pz) / dz;
            double t_hi = (z_hi - pz) / dz;
            if(t_lo > t_hi) std::swap(t_lo, t_hi);
            double rxy_at_lo = (px + t_lo*dx)*(px + t_lo*dx) + (py + t_lo*dy)*(py + t_lo*dy);
            double rxy_at_hi = (px + t_hi*dx)*(px + t_hi*dx) + (py + t_hi*dy)*(py + t_hi*dy);
            double rxy_at_closest = px*px + py*py;
            if(C > 0) {
                double t_closest = -(px*dx + py*dy) / C;
                if(t_closest > t_lo && t_closest < t_hi) {
                    rxy_at_closest = (px + t_closest*dx)*(px + t_closest*dx) + (py + t_closest*dy)*(py + t_closest*dy);
                } else {
                    rxy_at_closest = std::min(rxy_at_lo, rxy_at_hi);
                }
            }
            // If closest approach within z-range is still outside bounding radius
            // AND both z-endpoints are outside, we can skip
            if(rxy_at_closest > bounding_r2 && rxy_at_lo > bounding_r2 && rxy_at_hi > bounding_r2) {
                return {};
            }
        }
    }

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

    // rmin_/rmax_ are apothems (perpendicular distance to edge midpoint).
    // Vertices sit at the circumradius = apothem / cos(pi/n).
    // face_dphi was computed in the early-out block above.
    double cr_scale = 1.0 / std::cos(0.5 * face_dphi);

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

        double rmax_lo = rmax_[seg] * cr_scale;
        double rmax_hi = rmax_[seg + 1] * cr_scale;
        double rmin_lo = rmin_[seg] * cr_scale;
        double rmin_hi = rmin_[seg + 1] * cr_scale;

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

            if(rmax_.front() >= GEOMETRY_PRECISION) {
                bool in_cap = has_phi_cut_
                    ? PointInAnnularSectorPolygon(intersection_x, intersection_y,
                                                   num_sides_, cp, sp,
                                                   rmin_.front() * cr_scale, rmax_.front() * cr_scale)
                    : PointInAnnularPolygon(intersection_x, intersection_y,
                                             num_sides_, cp, sp,
                                             rmin_.front() * cr_scale, rmax_.front() * cr_scale);
                if(in_cap) {
                    add_hit(t, 0, (dz > 0), intersection_x, intersection_y, intersection_z);
                }
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

            if(rmax_.back() >= GEOMETRY_PRECISION) {
                bool in_cap = has_phi_cut_
                    ? PointInAnnularSectorPolygon(intersection_x, intersection_y,
                                                   num_sides_, cp, sp,
                                                   rmin_.back() * cr_scale, rmax_.back() * cr_scale)
                    : PointInAnnularPolygon(intersection_x, intersection_y,
                                             num_sides_, cp, sp,
                                             rmin_.back() * cr_scale, rmax_.back() * cr_scale);
                if(in_cap) {
                    add_hit(t, 0, (dz < 0), intersection_x, intersection_y, intersection_z);
                }
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
            double rmin_below = rmin_[i] * cr_scale;
            double rmax_below = rmax_[i] * cr_scale;
            double rmin_above = rmin_[j] * cr_scale;
            double rmax_above = rmax_[j] * cr_scale;

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
                if(rhi >= GEOMETRY_PRECISION) {
                    bool in_cap = has_phi_cut_
                        ? PointInAnnularSectorPolygon(intersection_x, intersection_y,
                                                       num_sides_, cp, sp, rlo, rhi)
                        : PointInAnnularPolygon(intersection_x, intersection_y,
                                                 num_sides_, cp, sp, rlo, rhi);
                    if(in_cap) {
                        bool entering = (rmax_above > rmax_below) ? (dz > 0) : (dz < 0);
                        add_hit(t, 0, entering, intersection_x, intersection_y, intersection_z);
                    }
                }
            }

            // Inner annulus: region where rmin changed
            if(rmin_below != rmin_above) {
                double rlo = std::min(rmin_below, rmin_above);
                double rhi = std::max(rmin_below, rmin_above);
                if(rhi >= GEOMETRY_PRECISION) {
                    bool in_cap = has_phi_cut_
                        ? PointInAnnularSectorPolygon(intersection_x, intersection_y,
                                                       num_sides_, cp, sp, rlo, rhi)
                        : PointInAnnularPolygon(intersection_x, intersection_y,
                                                 num_sides_, cp, sp, rlo, rhi);
                    if(in_cap) {
                        bool entering = (rmin_above > rmin_below) ? (dz < 0) : (dz > 0);
                        add_hit(t, 0, entering, intersection_x, intersection_y, intersection_z);
                    }
                }
            }

            i = j + 1;
        }
    }

    // Sort hits by distance
    auto cmp = [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    };

    // Net-parity dedup: collapse coincident hits from shared face boundaries.
    auto net_parity_dedup = [](Intersection* hits, int count) -> std::vector<Intersection> {
        std::vector<Intersection> result;
        int i = 0;
        while(i < count) {
            int j = i + 1;
            while(j < count && std::fabs(hits[j].distance - hits[i].distance) <= GEOMETRY_PRECISION)
                ++j;
            int net = 0;
            for(int k = i; k < j; ++k)
                net += hits[k].entering ? 1 : -1;
            if(net > 0) {
                hits[i].entering = true;
                result.push_back(hits[i]);
            } else if(net < 0) {
                hits[i].entering = false;
                result.push_back(hits[i]);
            }
            i = j;
        }
        return result;
    };

    if(!has_phi_cut_) {
        if(using_heap) {
            std::sort(heap_hits.begin(), heap_hits.end(), cmp);
            return net_parity_dedup(heap_hits.data(), n_hits);
        }
        std::sort(stack_hits, stack_hits + n_hits, cmp);
        return net_parity_dedup(stack_hits, n_hits);
    }

    // Phi boundary face intersections for partial-phi polyhedra.
    // Each boundary is a half-plane from the z-axis; we intersect the ray
    // with the plane, then check z and radial bounds per z-section.
    bool z_axis_phi_hit = false;
    for(int face = 0; face < 2; ++face) {
        double cos_a = (face == 0) ? cp[0] : cp[num_sides_];
        double sin_a = (face == 0) ? sp[0] : sp[num_sides_];
        double nx, ny;
        if(face == 0) { nx = sin_a; ny = -cos_a; }
        else { nx = -sin_a; ny = cos_a; }

        double n_dot_d = nx * dx + ny * dy;
        if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue;

        double n_dot_p = nx * px + ny * py;
        double t = -n_dot_p / n_dot_d;
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

        double hx = px + t * dx;
        double hy = py + t * dy;
        double hz = pz + t * dz;

        double r_hit = hx * cos_a + hy * sin_a;
        if(r_hit < -GEOMETRY_PRECISION) continue;

        if(r_hit < GEOMETRY_PRECISION * 1e3) {
            if(z_axis_phi_hit) continue;
            z_axis_phi_hit = true;
        }

        for(size_t seg = 0; seg + 1 < n; ++seg) {
            double z_lo = z_planes_[seg];
            double z_hi = z_planes_[seg + 1];
            if(z_hi <= z_lo) continue;
            if(hz < z_lo - GEOMETRY_PRECISION || hz > z_hi + GEOMETRY_PRECISION) continue;

            double frac = (hz - z_lo) / (z_hi - z_lo);
            frac = std::max(0.0, std::min(1.0, frac));
            double rmin_here = rmin_[seg] + frac * (rmin_[seg + 1] - rmin_[seg]);
            double rmax_here = rmax_[seg] + frac * (rmax_[seg + 1] - rmax_[seg]);

            if(r_hit >= rmin_here * cr_scale - GEOMETRY_PRECISION &&
               r_hit <= rmax_here * cr_scale + GEOMETRY_PRECISION) {
                add_hit(t, 0, (n_dot_d < 0), hx, hy, hz);
                break;
            }
        }
    }

    if(using_heap) {
        std::sort(heap_hits.begin(), heap_hits.end(), cmp);
        return net_parity_dedup(heap_hits.data(), n_hits);
    }
    std::sort(stack_hits, stack_hits + n_hits, cmp);
    return net_parity_dedup(stack_hits, n_hits);
}

// ------------------------------------------------------------------------- //
AABB Polyhedra::GetBoundingBox() const {
    if(z_planes_.empty()) {
        return AABB(math::Vector3D(0, 0, 0), math::Vector3D(0, 0, 0));
    }
    double max_r = *std::max_element(rmax_.begin(), rmax_.end());
    double face_dphi = has_phi_cut_ ? delta_phi_ / num_sides_ : 2.0 * M_PI / num_sides_;
    double circum_r = max_r / std::cos(0.5 * face_dphi);
    double z_lo = z_planes_.front();
    double z_hi = z_planes_.back();

    if(!has_phi_cut_) {
        return AABB(
            math::Vector3D(-circum_r, -circum_r, z_lo),
            math::Vector3D( circum_r,  circum_r, z_hi)
        );
    }

    double sp = phi_utils::NormalizePhi(start_phi_);
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
        if(cn < sp) cn += phi_utils::TWO_PI;
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
