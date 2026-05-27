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
#include "GeometryMacros.h"
#include "PhiUtils.h"

namespace siren {
namespace geometry {

void Polycone::validate() {
    if(z_planes_.size() < 2) {
        throw std::runtime_error("Polycone requires at least 2 z-planes!");
    }
    if(z_planes_.size() != rmin_.size() || z_planes_.size() != rmax_.size()) {
        throw std::runtime_error("Polycone z_planes, rmin, and rmax vectors must have the same size!");
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
        throw std::runtime_error("Polycone z-planes must be monotonic (ascending or descending)!");
    }
    if(descending && !ascending) {
        std::reverse(z_planes_.begin(), z_planes_.end());
        std::reverse(rmin_.begin(), rmin_.end());
        std::reverse(rmax_.begin(), rmax_.end());
    }
    for(size_t i = 0; i < z_planes_.size(); ++i) {
        if(rmin_[i] < 0 || rmax_[i] < 0) {
            throw std::runtime_error("Polycone radii must be non-negative!");
        }
        if(rmin_[i] > rmax_[i]) {
            throw std::runtime_error("Polycone inner radius must not exceed outer radius at each z-plane!");
        }
    }
}

Polycone::Polycone() : Geometry("Polycone"), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Polycone::Polycone(std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax) : Geometry("Polycone"), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { validate(); RecomputeWorldAABB(); }
Polycone::Polycone(Placement const & p) : Geometry("Polycone", p), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { RecomputeWorldAABB(); }
Polycone::Polycone(Placement const & p, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax) : Geometry("Polycone", p), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) { validate(); RecomputeWorldAABB(); }
Polycone::Polycone(std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax, double start_phi, double delta_phi) : Geometry("Polycone"), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(start_phi), delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    RecomputeWorldAABB();
}
Polycone::Polycone(Placement const & p, std::vector<double> const & z_planes, std::vector<double> const & rmin, std::vector<double> const & rmax, double start_phi, double delta_phi) : Geometry("Polycone", p), z_planes_(z_planes), rmin_(rmin), rmax_(rmax), start_phi_(start_phi), delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    RecomputeWorldAABB();
}
Polycone::Polycone(const Polycone& o) : Geometry(o), z_planes_(o.z_planes_), rmin_(o.rmin_), rmax_(o.rmax_), start_phi_(o.start_phi_), delta_phi_(o.delta_phi_), has_phi_cut_(o.has_phi_cut_) { RecomputeWorldAABB(); }

// Swap includes derived has_phi_cut_; equal/less use only primary fields
SIREN_GEOMETRY_SWAP(Polycone, z_planes_, rmin_, rmax_, start_phi_, delta_phi_, has_phi_cut_)
SIREN_GEOMETRY_ASSIGN(Polycone)
SIREN_GEOMETRY_EQUAL(Polycone, z_planes_, rmin_, rmax_, start_phi_, delta_phi_)
SIREN_GEOMETRY_LESS(Polycone, z_planes_, rmin_, rmax_, start_phi_, delta_phi_)

void Polycone::print(std::ostream& os) const {
    os << "Polycone(" << z_planes_.size() << " z-planes";
    if(has_phi_cut_) os << ", " << start_phi_ << ", " << delta_phi_;
    os << ")\n";
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

    // Origin shift: move ray origin to closest approach to the coordinate
    // origin. This keeps quadratic coefficients small for far-field rays,
    // avoiding catastrophic cancellation in the discriminant B^2 - 4AC.
    double t_shift = -(px * dx + py * dy + pz * dz);
    double qx = px + t_shift * dx;
    double qy = py + t_shift * dy;
    double qz = pz + t_shift * dz;

    static constexpr int STACK_CAP = 32;
    Intersection stack_hits[STACK_CAP];
    std::vector<Intersection> heap_hits;
    bool using_heap = false;
    int n_hits = 0;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    // Write a hit, promoting from stack to heap if the stack fills up.
    auto emit = [&](double t, bool entering) {
        Intersection isect;
        isect.distance = t;
        isect.hierarchy = 0;
        isect.entering = entering;
        isect.position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
        if(!using_heap) {
            if(n_hits < STACK_CAP) {
                stack_hits[n_hits++] = isect;
            } else {
                using_heap = true;
                heap_hits.assign(stack_hits, stack_hits + STACK_CAP);
                heap_hits.push_back(isect);
                n_hits++;
            }
        } else {
            heap_hits.push_back(isect);
            n_hits++;
        }
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
            // Substituting ray with shifted origin (q) gives: A*s^2 + B*s + C = 0
            // where s = t - t_shift (roots are shifted back to absolute t).
            double A = dx * dx + dy * dy - b_outer * b_outer * dz * dz;
            double B = 2.0 * (qx * dx + qy * dy - b_outer * dz * (a_outer + b_outer * qz));
            double C = qx * qx + qy * qy - (a_outer + b_outer * qz) * (a_outer + b_outer * qz);

            if(!(dx == 0 && dy == 0 && b_outer == 0)) {
                double determinant = B * B - 4.0 * A * C;

                if(std::fabs(A) > GEOMETRY_PRECISION) {
                    if(determinant > 0) {
                        double sqrt_det = std::sqrt(determinant);
                        double t1 = (-B + sqrt_det) / (2.0 * A) + t_shift;
                        double t2 = (-B - sqrt_det) / (2.0 * A) + t_shift;

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
                    double t1 = -C / B + t_shift;
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
            double B = 2.0 * (qx * dx + qy * dy - b_inner * dz * (a_inner + b_inner * qz));
            double C = qx * qx + qy * qy - (a_inner + b_inner * qz) * (a_inner + b_inner * qz);

            if(!(dx == 0 && dy == 0 && b_inner == 0)) {
                double determinant = B * B - 4.0 * A * C;

                if(std::fabs(A) > GEOMETRY_PRECISION) {
                    if(determinant > 0) {
                        double sqrt_det = std::sqrt(determinant);
                        double t1 = (-B + sqrt_det) / (2.0 * A) + t_shift;
                        double t2 = (-B - sqrt_det) / (2.0 * A) + t_shift;

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
                    double t1 = -C / B + t_shift;
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
    if(std::fabs(dz) > GEOMETRY_PRECISION) {
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
            double rmax2 = rmax_cap * rmax_cap;
            double rmin2 = rmin_cap * rmin_cap;
            double r2_tol = GEOMETRY_PRECISION * std::fmax(1.0, rmax2);
            if(r2_hit <= rmax2 + r2_tol && r2_hit >= rmin2 - r2_tol) {
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
            double rmax2 = rmax_cap * rmax_cap;
            double rmin2 = rmin_cap * rmin_cap;
            double r2_tol = GEOMETRY_PRECISION * std::fmax(1.0, rmax2);
            if(r2_hit <= rmax2 + r2_tol && r2_hit >= rmin2 - r2_tol) {
                intersection_z = pz + t * dz;
                emit(t, dz < 0);
            }
        }
    }

    // --- Internal caps for step changes in radius ---
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
        if(using_heap) {
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
    if(using_heap) {
        std::sort(heap_hits.begin(), heap_hits.end(), cmp);
    } else {
        std::sort(stack_hits, stack_hits + n_hits, cmp);
    }

    std::vector<TaggedHit> all_hits;
    all_hits.reserve(n_hits + 2);
    for(int i = 0; i < n_hits; ++i) {
        Intersection const & h = using_heap ? heap_hits[i] : stack_hits[i];
        all_hits.push_back({h.distance, h.position, h.entering, 0});
    }

    // Compute infinite wedge intersections (two half-planes from z-axis)
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
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

        double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;

        bool entering = (n_dot_d < 0);
        if(hx*hx + hy*hy < GEOMETRY_PRECISION * 1e3 * GEOMETRY_PRECISION * 1e3) {
            if(z_axis_hit_emitted) continue;
            z_axis_hit_emitted = true;
            entering = phi_utils::ZAxisWedgeEntering(dx, dy, start_phi_, delta_phi_);
        }
        all_hits.push_back({t, siren::math::Vector3D(hx, hy, hz), entering, 1});
    }

    if(all_hits.empty()) return {};

    std::sort(all_hits.begin(), all_hits.end(), [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    bool in_surface = false;
    bool in_wedge = phi_utils::InitialPhiState(px, py, dx, dy, start_phi_, delta_phi_);

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
    double sp = phi_utils::NormalizePhi(start_phi_);
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
        if(cn < sp) cn += phi_utils::TWO_PI;
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
