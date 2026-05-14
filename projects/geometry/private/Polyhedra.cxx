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

    return
        std::tie(num_sides_, start_phi_, z_planes_, rmin_, rmax_)
        <
        std::tie(polyhedra->num_sides_, polyhedra->start_phi_, polyhedra->z_planes_, polyhedra->rmin_, polyhedra->rmax_);
}

void Polyhedra::print(std::ostream& os) const
{
    os << "Polyhedra with " << num_sides_ << " sides, "
       << z_planes_.size() << " z-planes, start_phi=" << start_phi_ << ":";
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
static bool PointInConvexPolygon(double px, double py,
                                 const double* vx, const double* vy, int nv) {
    bool has_pos = false;
    bool has_neg = false;
    for(int i = 0; i < nv; ++i) {
        int j = (i + 1) % nv;
        double ex = vx[j] - vx[i];
        double ey = vy[j] - vy[i];
        double cross = ex * (py - vy[i]) - ey * (px - vx[i]);
        if(cross > 0) has_pos = true;
        if(cross < 0) has_neg = true;
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

    // Fixed-size hit buffer to avoid heap allocation.
    // Max theoretical hits = 2*(num_sides*num_segments + num_segments+1),
    // but realistic polyhedra produce far fewer.
    static constexpr int MAX_HITS = 64;
    Intersection hits[MAX_HITS];
    int n_hits = 0;

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
            {
                double v0x = rmax_lo * cos_k,  v0y = rmax_lo * sin_k,  v0z = z_lo;
                double v1x = rmax_lo * cos_k1, v1y = rmax_lo * sin_k1, v1z = z_lo;
                double v3x = rmax_hi * cos_k,  v3y = rmax_hi * sin_k,  v3z = z_hi;

                // Two edge vectors from v0
                double e1x = v1x - v0x, e1y = v1y - v0y, e1z = v1z - v0z;
                double e2x = v3x - v0x, e2y = v3y - v0y, e2z = v3z - v0z;

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

                        // Check if hit point is inside the quad using
                        // the cross-product winding test.
                        double v2x = rmax_hi * cos_k1, v2y = rmax_hi * sin_k1, v2z = z_hi;

                        double qx = intersection_x, qy = intersection_y, qz = intersection_z;

                        // Cross product of (v1-v0) x (q-v0) dot n
                        double c0 = nx * ((v1y - v0y) * (qz - v0z) - (v1z - v0z) * (qy - v0y))
                                  + ny * ((v1z - v0z) * (qx - v0x) - (v1x - v0x) * (qz - v0z))
                                  + nz * ((v1x - v0x) * (qy - v0y) - (v1y - v0y) * (qx - v0x));

                        double c1 = nx * ((v2y - v1y) * (qz - v1z) - (v2z - v1z) * (qy - v1y))
                                  + ny * ((v2z - v1z) * (qx - v1x) - (v2x - v1x) * (qz - v1z))
                                  + nz * ((v2x - v1x) * (qy - v1y) - (v2y - v1y) * (qx - v1x));

                        double c2 = nx * ((v3y - v2y) * (qz - v2z) - (v3z - v2z) * (qy - v2y))
                                  + ny * ((v3z - v2z) * (qx - v2x) - (v3x - v2x) * (qz - v2z))
                                  + nz * ((v3x - v2x) * (qy - v2y) - (v3y - v2y) * (qx - v2x));

                        double c3 = nx * ((v0y - v3y) * (qz - v3z) - (v0z - v3z) * (qy - v3y))
                                  + ny * ((v0z - v3z) * (qx - v3x) - (v0x - v3x) * (qz - v3z))
                                  + nz * ((v0x - v3x) * (qy - v3y) - (v0y - v3y) * (qx - v3x));

                        if((c0 >= -GEOMETRY_PRECISION && c1 >= -GEOMETRY_PRECISION &&
                             c2 >= -GEOMETRY_PRECISION && c3 >= -GEOMETRY_PRECISION) ||
                            (c0 <= GEOMETRY_PRECISION && c1 <= GEOMETRY_PRECISION &&
                             c2 <= GEOMETRY_PRECISION && c3 <= GEOMETRY_PRECISION)) {

                            if(n_hits < MAX_HITS) {
                                hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                                hits[n_hits].distance = t;
                                hits[n_hits].hierarchy = 0;
                                hits[n_hits].entering = (n_dot_d < 0);
                                ++n_hits;
                            }
                        }
                    }
                }
            }

            // ---- Inner lateral face (if hollow) ----
            if(rmin_lo > 0 || rmin_hi > 0) {
                double v0x = rmin_lo * cos_k,  v0y = rmin_lo * sin_k,  v0z = z_lo;
                double v1x = rmin_lo * cos_k1, v1y = rmin_lo * sin_k1, v1z = z_lo;
                double v3x = rmin_hi * cos_k,  v3y = rmin_hi * sin_k,  v3z = z_hi;

                // Edge vectors from v0
                double e1x = v1x - v0x, e1y = v1y - v0y, e1z = v1z - v0z;
                double e2x = v3x - v0x, e2y = v3y - v0y, e2z = v3z - v0z;

                // Normal = e1 x e2 (points outward from center, but for the
                // inner surface the "outside" of the solid is inward)
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

                        double c0 = nx * ((v1y - v0y) * (qz - v0z) - (v1z - v0z) * (qy - v0y))
                                  + ny * ((v1z - v0z) * (qx - v0x) - (v1x - v0x) * (qz - v0z))
                                  + nz * ((v1x - v0x) * (qy - v0y) - (v1y - v0y) * (qx - v0x));

                        double c1 = nx * ((v2y - v1y) * (qz - v1z) - (v2z - v1z) * (qy - v1y))
                                  + ny * ((v2z - v1z) * (qx - v1x) - (v2x - v1x) * (qz - v1z))
                                  + nz * ((v2x - v1x) * (qy - v1y) - (v2y - v1y) * (qx - v1x));

                        double c2 = nx * ((v3y - v2y) * (qz - v2z) - (v3z - v2z) * (qy - v2y))
                                  + ny * ((v3z - v2z) * (qx - v2x) - (v3x - v2x) * (qz - v2z))
                                  + nz * ((v3x - v2x) * (qy - v2y) - (v3y - v2y) * (qx - v2x));

                        double c3 = nx * ((v0y - v3y) * (qz - v3z) - (v0z - v3z) * (qy - v3y))
                                  + ny * ((v0z - v3z) * (qx - v3x) - (v0x - v3x) * (qz - v3z))
                                  + nz * ((v0x - v3x) * (qy - v3y) - (v0y - v3y) * (qx - v3x));

                        if((c0 >= -GEOMETRY_PRECISION && c1 >= -GEOMETRY_PRECISION &&
                             c2 >= -GEOMETRY_PRECISION && c3 >= -GEOMETRY_PRECISION) ||
                            (c0 <= GEOMETRY_PRECISION && c1 <= GEOMETRY_PRECISION &&
                             c2 <= GEOMETRY_PRECISION && c3 <= GEOMETRY_PRECISION)) {

                            // Inner surface: entering is inverted relative to
                            // the normal direction.
                            if(n_hits < MAX_HITS) {
                                hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                                hits[n_hits].distance = t;
                                hits[n_hits].hierarchy = 0;
                                hits[n_hits].entering = (n_dot_d > 0);
                                ++n_hits;
                            }
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

            if(PointInAnnularPolygon(intersection_x, intersection_y,
                                       num_sides_, cp, sp,
                                       rmin_.front(), rmax_.front())) {
                if(n_hits < MAX_HITS) {
                    hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = (dz > 0);
                    ++n_hits;
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

            if(PointInAnnularPolygon(intersection_x, intersection_y,
                                       num_sides_, cp, sp,
                                       rmin_.back(), rmax_.back())) {
                if(n_hits < MAX_HITS) {
                    hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = (dz < 0);
                    ++n_hits;
                }
            }
        }
    }

    // ---- Internal annular caps at each intermediate z-plane ----
    // These arise where the radii change discontinuously between sections.
    // Skip caps where adjacent sections have continuous radii.
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
        for(size_t i = 1; i + 1 < n; ++i) {
            // Check if the section below and above are both valid (non-degenerate)
            bool below_valid = (z_planes_[i] > z_planes_[i - 1]);
            bool above_valid = (z_planes_[i + 1] > z_planes_[i]);

            // Get the radii from the section below (at its top) and above (at its bottom)
            double rmin_below = below_valid ? rmin_[i] : 0;
            double rmax_below = below_valid ? rmax_[i] : 0;
            double rmin_above = above_valid ? rmin_[i] : 0;
            double rmax_above = above_valid ? rmax_[i] : 0;

            // For continuous radii these are equal, so skip
            if(rmax_below == rmax_above && rmin_below == rmin_above && below_valid && above_valid) {
                continue;
            }

            double z_cap = z_planes_[i];
            double t = (z_cap - pz) / dz;

            if(t > 0 && t < GEOMETRY_PRECISION)
                t = 0;

            intersection_x = px + t * dx;
            intersection_y = py + t * dy;
            intersection_z = pz + t * dz;

            if(PointInAnnularPolygon(intersection_x, intersection_y,
                                       num_sides_, cp, sp,
                                       rmin_[i], rmax_[i])) {
                // Determine entering based on step direction
                bool step_faces_up = (rmax_above > rmax_below) || (rmin_above < rmin_below);
                bool cap_entering = step_faces_up ? (dz < 0) : (dz > 0);
                if(n_hits < MAX_HITS) {
                    hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                    hits[n_hits].distance = t;
                    hits[n_hits].hierarchy = 0;
                    hits[n_hits].entering = cap_entering;
                    ++n_hits;
                }
            }
        }
    }

    // Sort hits by distance using a plain lambda (no std::function overhead)
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });

    return std::vector<Intersection>(hits, hits + n_hits);
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Polyhedra::ComputeDistanceToBorder(const siren::math::Vector3D& position, const siren::math::Vector3D& direction) const
{
    // Compute the surface intersections
    std::vector<Intersection> intersections = Intersections(position, direction);
    std::vector<double> dists;
    bool first = true;
    for(unsigned int i=0; i<intersections.size(); ++i) {
        Intersection const & obj = intersections[i];
        if(obj.distance > 0) {
            if(first) {
                first = false;
                dists.push_back(obj.distance);
                if(not obj.entering) {
                    break;
                }
            }
            else {
                if(not obj.entering) {
                    dists.push_back(obj.distance);
                    break;
                }
                else {
                    throw(std::runtime_error("There should never be two \"entering\" intersections in a row!"));
                }
            }
        }
    }

    std::pair<double, double> distance;

    if(dists.size() < 1) {
        distance.first  = -1;
        distance.second = -1;
    } else if(dists.size() == 1) {
        distance.first  = dists.at(0);
        distance.second = -1;

    } else if(dists.size() == 2) {
        distance.first  = dists.at(0);
        distance.second = dists.at(1);

        if(distance.second < distance.first) {
            std::swap(distance.first, distance.second);
        }

    } else {
        //log_error("This point should never be reached");
    }

    if(distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if(distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if(distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

// ------------------------------------------------------------------------- //
AABB Polyhedra::GetBoundingBox() const {
    if(z_planes_.empty()) {
        return AABB(math::Vector3D(0, 0, 0), math::Vector3D(0, 0, 0));
    }
    double max_r = *std::max_element(rmax_.begin(), rmax_.end());
    // The inscribed radius is rmax (distance from center to edge midpoint).
    // The circumradius (distance from center to vertex) is rmax / cos(pi/N).
    double circum_r = max_r / std::cos(M_PI / num_sides_);
    double z_lo = z_planes_.front();
    double z_hi = z_planes_.back();
    return AABB(
        math::Vector3D(-circum_r, -circum_r, z_lo),
        math::Vector3D( circum_r,  circum_r, z_hi)
    );
}

} // namespace geometry
} // namespace siren
