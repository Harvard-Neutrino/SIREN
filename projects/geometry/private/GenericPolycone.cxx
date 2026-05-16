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

GenericPolycone::GenericPolycone()
    : Geometry((std::string)("GenericPolycone"))
    , r_()
    , z_()
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false)
    , winding_(1) {
}

GenericPolycone::GenericPolycone(std::vector<double> const & r,
                                 std::vector<double> const & z)
    : Geometry((std::string)("GenericPolycone"))
    , r_(r)
    , z_(z)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    validate();
    ComputeWinding();
}

GenericPolycone::GenericPolycone(Placement const & placement,
                                 std::vector<double> const & r,
                                 std::vector<double> const & z)
    : Geometry((std::string)("GenericPolycone"), placement)
    , r_(r)
    , z_(z)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    validate();
    ComputeWinding();
}

GenericPolycone::GenericPolycone(std::vector<double> const & r,
                                 std::vector<double> const & z,
                                 double start_phi, double delta_phi)
    : Geometry((std::string)("GenericPolycone"))
    , r_(r)
    , z_(z)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("GenericPolycone delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    ComputeWinding();
}

GenericPolycone::GenericPolycone(Placement const & placement,
                                 std::vector<double> const & r,
                                 std::vector<double> const & z,
                                 double start_phi, double delta_phi)
    : Geometry((std::string)("GenericPolycone"), placement)
    , r_(r)
    , z_(z)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    validate();
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("GenericPolycone delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
    ComputeWinding();
}

GenericPolycone::GenericPolycone(const GenericPolycone& other)
    : Geometry(other)
    , r_(other.r_)
    , z_(other.z_)
    , start_phi_(other.start_phi_)
    , delta_phi_(other.delta_phi_)
    , has_phi_cut_(other.has_phi_cut_)
    , winding_(other.winding_) {
}

void GenericPolycone::swap(Geometry& geometry) {
    GenericPolycone* other = dynamic_cast<GenericPolycone*>(&geometry);
    if(!other) return;
    Geometry::swap(*other);
    std::swap(r_, other->r_);
    std::swap(z_, other->z_);
    std::swap(start_phi_, other->start_phi_);
    std::swap(delta_phi_, other->delta_phi_);
    std::swap(has_phi_cut_, other->has_phi_cut_);
    std::swap(winding_, other->winding_);
}

GenericPolycone& GenericPolycone::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const GenericPolycone* other = dynamic_cast<const GenericPolycone*>(&geometry);
        if(!other) return *this;
        GenericPolycone tmp(*other);
        swap(tmp);
    }
    return *this;
}

bool GenericPolycone::equal(const Geometry& geometry) const {
    const GenericPolycone* other = dynamic_cast<const GenericPolycone*>(&geometry);
    if(!other) return false;
    if(r_ != other->r_) return false;
    if(z_ != other->z_) return false;
    if(start_phi_ != other->start_phi_) return false;
    if(delta_phi_ != other->delta_phi_) return false;
    return true;
}

bool GenericPolycone::less(const Geometry& geometry) const {
    const GenericPolycone* other = dynamic_cast<const GenericPolycone*>(&geometry);
    if(!other) return false;
    return std::tie(r_, z_, start_phi_, delta_phi_)
         < std::tie(other->r_, other->z_, other->start_phi_, other->delta_phi_);
}

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

        double r0 = r_[i0], z0 = z_[i0];
        double r1 = r_[i1], z1 = z_[i1];

        double dr_edge = r1 - r0;
        double dz_edge = z1 - z0;

        if(std::fabs(dr_edge) < GEOMETRY_PRECISION && std::fabs(dz_edge) < GEOMETRY_PRECISION) {
            continue;
        }

        if(std::fabs(dz_edge) > GEOMETRY_PRECISION) {
            // Non-horizontal edge: conical/cylindrical surface
            // r(z) = a + b*z where b = dr/dz, a = r0 - b*z0
            double b = dr_edge / dz_edge;
            double a = r0 - b * z0;

            // Surface: x^2 + y^2 = (a + b*z)^2
            double A = dx*dx + dy*dy - b*b*dz*dz;
            double B = 2.0 * (px*dx + py*dy - b*dz*(a + b*pz));
            double C = px*px + py*py - (a + b*pz)*(a + b*pz);

            bool flip = (winding_ * dz_edge < 0);

            auto try_hit = [&](double t) {
                if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
                intersection_z = pz + t * dz;
                // Clip to edge parameter range: s in (0, 1)
                double s = (intersection_z - z0) / dz_edge;
                if(s <= 0 || s >= 1) return;
                intersection_x = px + t * dx;
                intersection_y = py + t * dy;
                double r_at_z = a + b * intersection_z;
                if(r_at_z < -GEOMETRY_PRECISION) return;
                bool entering = flip ? !entering_cone(a, b) : entering_cone(a, b);
                emit(t, entering);
            };

            if(std::fabs(A) > GEOMETRY_PRECISION) {
                double det = B*B - 4.0*A*C;
                if(det > 0) {
                    double sqrt_det = std::sqrt(det);
                    try_hit((-B + sqrt_det) / (2.0 * A));
                    try_hit((-B - sqrt_det) / (2.0 * A));
                }
            } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                try_hit(-C / B);
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

            // Clip to edge parameter range: s in (0, 1)
            if(std::fabs(dr_edge) < GEOMETRY_PRECISION) continue;
            double s = (r_hit - r0) / dr_edge;
            if(s <= 0 || s >= 1) continue;

            // Entering: D . N_outward < 0 where N_outward z-component = -winding * dr_edge
            bool entering = (dz * winding_ * dr_edge > 0);
            emit(t, entering);
        }
    }

    if(n_hits == 0) return {};

    auto cmp = [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    };

    if(!has_phi_cut_) {
        if(use_heap) {
            std::sort(heap_hits.begin(), heap_hits.end(), cmp);
            return heap_hits;
        }
        std::sort(stack_hits, stack_hits + n_hits, cmp);
        return {stack_hits, stack_hits + n_hits};
    }

    // Phi cut: merge surface hits with infinite wedge hits and run CSG walk.
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
