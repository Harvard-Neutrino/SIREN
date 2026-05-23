#include "SIREN/geometry/GenericPolycone.h"

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

namespace {
using phi_utils::NormalizePhi;
using phi_utils::PhiInRange;
using phi_utils::InitialPhiState;
using phi_utils::ZAxisWedgeEntering;
using phi_utils::TWO_PI;
} // anonymous namespace

void GenericPolycone::validate() const {
    if(r_.size() < 3) {
        throw std::runtime_error("GenericPolycone requires at least 3 vertices!");
    }
    if(r_.size() != z_.size()) {
        throw std::runtime_error("GenericPolycone r and z vectors must have the same size!");
    }
    for(size_t i = 0; i < r_.size(); ++i) {
        if(r_[i] < 0) {
            throw std::runtime_error("GenericPolycone radii must be non-negative!");
        }
    }
}

void GenericPolycone::ComputeWinding() {
    double signed_area = 0;
    size_t n = r_.size();
    for(size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        signed_area += r_[i] * z_[j] - r_[j] * z_[i];
    }
    winding_ = (signed_area >= 0) ? 1 : -1;
}

GenericPolycone::GenericPolycone() : Geometry("GenericPolycone"), r_(), z_(), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false), winding_(1) { RecomputeWorldAABB(); }
GenericPolycone::GenericPolycone(std::vector<double> const & r, std::vector<double> const & z) : Geometry("GenericPolycone"), r_(r), z_(z), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) {
    validate();
    ComputeWinding();
    RecomputeWorldAABB();
}
GenericPolycone::GenericPolycone(Placement const & p, std::vector<double> const & r, std::vector<double> const & z) : Geometry("GenericPolycone", p), r_(r), z_(z), start_phi_(0), delta_phi_(2.0 * M_PI), has_phi_cut_(false) {
    validate();
    ComputeWinding();
    RecomputeWorldAABB();
}
GenericPolycone::GenericPolycone(std::vector<double> const & r, std::vector<double> const & z, double start_phi, double delta_phi) : Geometry("GenericPolycone"), r_(r), z_(z), start_phi_(start_phi), delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    ComputeWinding();
    RecomputeWorldAABB();
}
GenericPolycone::GenericPolycone(Placement const & p, std::vector<double> const & r, std::vector<double> const & z, double start_phi, double delta_phi) : Geometry("GenericPolycone", p), r_(r), z_(z), start_phi_(start_phi), delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0) throw std::invalid_argument("delta_phi must be positive!"); if(delta_phi_ > 2.0 * M_PI) delta_phi_ = 2.0 * M_PI;
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    ComputeWinding();
    RecomputeWorldAABB();
}
GenericPolycone::GenericPolycone(const GenericPolycone& o) : Geometry(o), r_(o.r_), z_(o.z_), start_phi_(o.start_phi_), delta_phi_(o.delta_phi_), has_phi_cut_(o.has_phi_cut_), winding_(o.winding_) { RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(GenericPolycone, r_, z_, start_phi_, delta_phi_, has_phi_cut_, winding_)
SIREN_GEOMETRY_ASSIGN(GenericPolycone)
SIREN_GEOMETRY_EQUAL(GenericPolycone, r_, z_, start_phi_, delta_phi_)
SIREN_GEOMETRY_LESS(GenericPolycone, r_, z_, start_phi_, delta_phi_)

void GenericPolycone::print(std::ostream& os) const {
    os << "GenericPolycone with " << r_.size() << " vertices:";
    for(size_t i = 0; i < r_.size(); ++i) {
        os << " [r=" << r_[i] << " z=" << z_[i] << "]";
    }
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    os << '\n';
}

// -------------------------------------------------------------------------
// ComputeIntersections
//
// The generic polycone is defined by a closed polygon in (r, z) space that
// is revolved around the z-axis. Each edge of the polygon, when revolved,
// produces a conical frustum (general case), a cylinder (vertical edge),
// or an annular disk (horizontal edge).
//
// For each edge we compute ray-surface intersections, clip to the edge's
// parameter range, and determine entering/exiting from the polygon winding.
// -------------------------------------------------------------------------
std::vector<Geometry::Intersection> GenericPolycone::ComputeIntersections(
        siren::math::Vector3D const & position,
        siren::math::Vector3D const & direction) const {

    size_t n = r_.size();
    if(n < 3) return {};

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    // Origin shift: move ray origin to closest approach to coordinate origin.
    // This keeps quadratic coefficients small for far-field rays, avoiding
    // catastrophic cancellation in the discriminant B^2 - 4AC.
    double t_shift = -(px * dx + py * dy + pz * dz);
    double qx = px + t_shift * dx;
    double qy = py + t_shift * dy;
    double qz = pz + t_shift * dz;

    size_t max_hits = 2 * n + 4;

    static constexpr int STACK_CAP = 32;
    Intersection stack_hits[STACK_CAP];
    std::vector<Intersection> heap_hits;
    bool use_heap = (max_hits > STACK_CAP);
    if(use_heap) heap_hits.reserve(max_hits);
    int n_hits = 0;

    double intersection_x;
    double intersection_y;
    double intersection_z;

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

    auto entering_cone = [&](double a, double b) {
        double r_val = a + b * intersection_z;
        return (intersection_x * dx + intersection_y * dy - b * r_val * dz) < 0;
    };

    for(size_t edge = 0; edge < n; ++edge) {
        size_t i0 = edge;
        size_t i1 = (edge + 1) % n;
        size_t i2 = (edge + 2) % n;

        double r0 = r_[i0], z0 = z_[i0];
        double r1 = r_[i1], z1 = z_[i1];

        double dr_edge = r1 - r0;
        double dz_edge = z1 - z0;

        if(std::fabs(dr_edge) < GEOMETRY_PRECISION && std::fabs(dz_edge) < GEOMETRY_PRECISION) {
            continue;
        }
        if(r0 < GEOMETRY_PRECISION && r1 < GEOMETRY_PRECISION) {
            continue;
        }

        // Include s=1 (claim the endpoint vertex) when the next edge cannot
        // produce a hit there: either the next edge is a horizontal disk, or
        // the shared vertex sits at r=0 on the axis (the next edge's surface
        // is degenerate and will never emit a hit at this vertex).
        double dz_next = z_[i2] - z_[i1];
        bool next_is_disk = std::fabs(dz_next) <= GEOMETRY_PRECISION;
        bool endpoint_on_axis = r1 < GEOMETRY_PRECISION;
        bool has_axis_vertex = r0 < GEOMETRY_PRECISION || r1 < GEOMETRY_PRECISION;
        bool claim_endpoint = next_is_disk || endpoint_on_axis;

        if(std::fabs(dz_edge) > GEOMETRY_PRECISION) {
            // Non-horizontal edge: conical/cylindrical surface
            // r(z) = a + b*z where b = dr/dz, a = r0 - b*z0
            double b = dr_edge / dz_edge;
            double a = r0 - b * z0;

            // Surface: x^2 + y^2 = (a + b*z)^2
            // Use shifted origin (q) for numerical stability at large distances.
            double A = dx*dx + dy*dy - b*b*dz*dz;
            double B = 2.0 * (qx*dx + qy*dy - b*dz*(a + b*qz));
            double C = qx*qx + qy*qy - (a + b*qz)*(a + b*qz);

            bool flip = (winding_ * dz_edge < 0);

            auto try_hit = [&](double t) {
                if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
                intersection_z = pz + t * dz;
                double s = (intersection_z - z0) / dz_edge;
                if(s < 0 || (claim_endpoint ? s > 1 : s >= 1)) return;
                intersection_x = px + t * dx;
                intersection_y = py + t * dy;
                double r_at_z = a + b * intersection_z;
                if(r_at_z < -GEOMETRY_PRECISION) return;
                bool entering;
                if(r_at_z < GEOMETRY_PRECISION) {
                    entering = (b * dz > 0);
                } else {
                    entering = entering_cone(a, b);
                }
                if(flip) entering = !entering;
                emit(t, entering);
            };

            if(std::fabs(A) > GEOMETRY_PRECISION) {
                double det = B*B - 4.0*A*C;
                if(det > 0) {
                    double sqrt_det = std::sqrt(det);
                    try_hit((-B + sqrt_det) / (2.0 * A) + t_shift);
                    try_hit((-B - sqrt_det) / (2.0 * A) + t_shift);
                } else if(has_axis_vertex && det > -std::fabs(A) * GEOMETRY_PRECISION) {
                    try_hit(-B / (2.0 * A) + t_shift);
                }
            } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                try_hit(-C / B + t_shift);
            }
        } else {
            // Horizontal edge: annular disk at z = z0
            if(std::fabs(dz) < GEOMETRY_PRECISION) continue;

            double t = (z0 - pz) / dz;
            if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

            intersection_x = px + t * dx;
            intersection_y = py + t * dy;
            intersection_z = pz + t * dz;

            double r_hit = std::sqrt(intersection_x*intersection_x + intersection_y*intersection_y);

            if(std::fabs(dr_edge) < GEOMETRY_PRECISION) continue;
            double s = (r_hit - r0) / dr_edge;
            bool disk_claim_endpoint = (r1 < GEOMETRY_PRECISION);
            if(s < 0 || (disk_claim_endpoint ? s > 1 : s >= 1)) continue;

            bool entering = (dz * winding_ * dr_edge > 0);
            emit(t, entering);
        }
    }

    if(n_hits == 0) return {};

    auto cmp = [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    };

    // Net-parity dedup: collapse coincident hits from shared edge boundaries.
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
        if(use_heap) {
            std::sort(heap_hits.begin(), heap_hits.end(), cmp);
            return net_parity_dedup(heap_hits.data(), n_hits);
        }
        std::sort(stack_hits, stack_hits + n_hits, cmp);
        return net_parity_dedup(stack_hits, n_hits);
    }

    // Phi cut: merge surface hits with infinite wedge hits and run CSG walk.
    // See Polycone.cxx for method description; same pattern in all phi-cut shapes.
    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

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
            entering = ZAxisWedgeEntering(dx, dy, start_phi_, delta_phi_);
        }
        all_hits.push_back({t, siren::math::Vector3D(hx, hy, hz), entering, 1});
    }

    if(all_hits.empty()) return {};

    std::sort(all_hits.begin(), all_hits.end(), [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    bool in_surface = false;
    bool in_wedge = InitialPhiState(px, py, dx, dy, start_phi_, delta_phi_);

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
    return net_parity_dedup(result.data(), (int)result.size());
}

AABB GenericPolycone::GetBoundingBox() const {
    if(r_.empty()) {
        return AABB(math::Vector3D(0, 0, 0), math::Vector3D(0, 0, 0));
    }
    double max_r = *std::max_element(r_.begin(), r_.end());
    double z_lo = *std::min_element(z_.begin(), z_.end());
    double z_hi = *std::max_element(z_.begin(), z_.end());

    if(!has_phi_cut_) {
        return AABB(
            math::Vector3D(-max_r, -max_r, z_lo),
            math::Vector3D( max_r,  max_r, z_hi)
        );
    }

    double sp = NormalizePhi(start_phi_);
    double ep = sp + delta_phi_;

    double x_min = 0, x_max = 0, y_min = 0, y_max = 0;

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

    double cardinals[4] = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    for(int i = 0; i < 4; ++i) {
        double c = cardinals[i];
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
