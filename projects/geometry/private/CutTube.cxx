#include "SIREN/geometry/CutTube.h"

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

// ------------------------------------------------------------------------- //
// Helper: normalize a vector and enforce the sign constraint on its z-component.
// For the low normal, z must be < 0; for the high normal, z must be > 0.
// ------------------------------------------------------------------------- //
static siren::math::Vector3D prepare_normal(siren::math::Vector3D n, bool want_positive_z) {
    double mag = n.magnitude();
    if(mag <= 0) {
        // Default to straight cap
        if(want_positive_z)
            return siren::math::Vector3D(0, 0, 1);
        else
            return siren::math::Vector3D(0, 0, -1);
    }
    n = n / mag;
    if(want_positive_z && n.GetZ() < 0) {
        n = -n;
    } else if(!want_positive_z && n.GetZ() > 0) {
        n = -n;
    }
    // If z-component is exactly 0, fall back to default
    if(n.GetZ() == 0) {
        if(want_positive_z)
            return siren::math::Vector3D(0, 0, 1);
        else
            return siren::math::Vector3D(0, 0, -1);
    }
    return n;
}

// ------------------------------------------------------------------------- //
CutTube::CutTube()
    : Geometry((std::string)("CutTube"))
    , rmin_(0.0)
    , rmax_(0.0)
    , dz_(0.0)
    , low_norm_(siren::math::Vector3D(0, 0, -1))
    , high_norm_(siren::math::Vector3D(0, 0, 1))
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    // Do nothing here
}

CutTube::CutTube(double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm)
    : Geometry((std::string)("CutTube"))
    , rmin_(rmin)
    , rmax_(rmax)
    , dz_(dz)
    , low_norm_(low_norm)
    , high_norm_(high_norm)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(rmin_ < 0) {
        throw std::invalid_argument("CutTube inner radius rmin must be non-negative!");
    }
    if(rmax_ < 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be non-negative!");
    }
    if(rmin_ > rmax_) {
        std::swap(rmin_, rmax_);
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be positive!");
    }
    if(dz_ <= 0) {
        throw std::invalid_argument("CutTube half-height dz must be positive!");
    }
    low_norm_ = prepare_normal(low_norm_, false);
    high_norm_ = prepare_normal(high_norm_, true);
}

CutTube::CutTube(Placement const & placement)
    : Geometry((std::string)("CutTube"), placement)
    , rmin_(0.0)
    , rmax_(0.0)
    , dz_(0.0)
    , low_norm_(siren::math::Vector3D(0, 0, -1))
    , high_norm_(siren::math::Vector3D(0, 0, 1))
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    // Do nothing here
}

CutTube::CutTube(Placement const & placement, double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm)
    : Geometry((std::string)("CutTube"), placement)
    , rmin_(rmin)
    , rmax_(rmax)
    , dz_(dz)
    , low_norm_(low_norm)
    , high_norm_(high_norm)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(rmin_ < 0) {
        throw std::invalid_argument("CutTube inner radius rmin must be non-negative!");
    }
    if(rmax_ < 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be non-negative!");
    }
    if(rmin_ > rmax_) {
        std::swap(rmin_, rmax_);
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be positive!");
    }
    if(dz_ <= 0) {
        throw std::invalid_argument("CutTube half-height dz must be positive!");
    }
    low_norm_ = prepare_normal(low_norm_, false);
    high_norm_ = prepare_normal(high_norm_, true);
}

CutTube::CutTube(double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm,
                 double start_phi, double delta_phi)
    : Geometry((std::string)("CutTube"))
    , rmin_(rmin)
    , rmax_(rmax)
    , dz_(dz)
    , low_norm_(low_norm)
    , high_norm_(high_norm)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(rmin_ < 0) {
        throw std::invalid_argument("CutTube inner radius rmin must be non-negative!");
    }
    if(rmax_ < 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be non-negative!");
    }
    if(rmin_ > rmax_) {
        std::swap(rmin_, rmax_);
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be positive!");
    }
    if(dz_ <= 0) {
        throw std::invalid_argument("CutTube half-height dz must be positive!");
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("CutTube delta_phi must be in (0, 2*pi]!");
    }
    low_norm_ = prepare_normal(low_norm_, false);
    high_norm_ = prepare_normal(high_norm_, true);
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

CutTube::CutTube(Placement const & placement, double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm,
                 double start_phi, double delta_phi)
    : Geometry((std::string)("CutTube"), placement)
    , rmin_(rmin)
    , rmax_(rmax)
    , dz_(dz)
    , low_norm_(low_norm)
    , high_norm_(high_norm)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(rmin_ < 0) {
        throw std::invalid_argument("CutTube inner radius rmin must be non-negative!");
    }
    if(rmax_ < 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be non-negative!");
    }
    if(rmin_ > rmax_) {
        std::swap(rmin_, rmax_);
    }
    if(rmax_ <= 0) {
        throw std::invalid_argument("CutTube outer radius rmax must be positive!");
    }
    if(dz_ <= 0) {
        throw std::invalid_argument("CutTube half-height dz must be positive!");
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("CutTube delta_phi must be in (0, 2*pi]!");
    }
    low_norm_ = prepare_normal(low_norm_, false);
    high_norm_ = prepare_normal(high_norm_, true);
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

CutTube::CutTube(const CutTube& other)
    : Geometry(other)
    , rmin_(other.rmin_)
    , rmax_(other.rmax_)
    , dz_(other.dz_)
    , low_norm_(other.low_norm_)
    , high_norm_(other.high_norm_)
    , start_phi_(other.start_phi_)
    , delta_phi_(other.delta_phi_)
    , has_phi_cut_(other.has_phi_cut_) {
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void CutTube::swap(Geometry& geometry) {
    CutTube* ct = dynamic_cast<CutTube*>(&geometry);
    if(!ct) {
        //log_warn("Cannot swap CutTube!");
        return;
    }

    Geometry::swap(*ct);

    std::swap(rmin_, ct->rmin_);
    std::swap(rmax_, ct->rmax_);
    std::swap(dz_, ct->dz_);
    std::swap(low_norm_, ct->low_norm_);
    std::swap(high_norm_, ct->high_norm_);
    std::swap(start_phi_, ct->start_phi_);
    std::swap(delta_phi_, ct->delta_phi_);
    std::swap(has_phi_cut_, ct->has_phi_cut_);
}

//------------------------------------------------------------------------- //
CutTube& CutTube::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const CutTube* ct = dynamic_cast<const CutTube*>(&geometry);
        if(!ct) {
            //log_warn("Cannot assign CutTube!");
            return *this;
        }

        CutTube tmp(*ct);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool CutTube::equal(const Geometry& geometry) const
{
    const CutTube* ct = dynamic_cast<const CutTube*>(&geometry);

    if(!ct)
        return false;
    else if(rmin_ != ct->rmin_)
        return false;
    else if(rmax_ != ct->rmax_)
        return false;
    else if(dz_ != ct->dz_)
        return false;
    else if(low_norm_ != ct->low_norm_)
        return false;
    else if(high_norm_ != ct->high_norm_)
        return false;
    else if(start_phi_ != ct->start_phi_)
        return false;
    else if(delta_phi_ != ct->delta_phi_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool CutTube::less(const Geometry& geometry) const
{
    const CutTube* ct = dynamic_cast<const CutTube*>(&geometry);
    if(!ct) return false;

    return
        std::tie(rmin_, rmax_, dz_, low_norm_, high_norm_, start_phi_, delta_phi_)
        <
        std::tie(ct->rmin_, ct->rmax_, ct->dz_, ct->low_norm_, ct->high_norm_, ct->start_phi_, ct->delta_phi_);
}

// ------------------------------------------------------------------------- //
void CutTube::print(std::ostream& os) const
{
    os << "rmin: " << rmin_ << "\trmax: " << rmax_
       << "\tdz: " << dz_
       << "\tlow_norm: " << low_norm_
       << "\thigh_norm: " << high_norm_;
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    os << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> CutTube::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // CutTube: a tube/cylinder with tilted end-cap planes.
    //
    // Surfaces:
    //   Outer barrel: x^2 + y^2 = rmax^2
    //   Inner barrel: x^2 + y^2 = rmin^2 (if rmin > 0)
    //   Low cut plane: passes through (0,0,-dz) with outward normal low_norm_
    //   High cut plane: passes through (0,0,+dz) with outward normal high_norm_
    //
    // A point is between the cut planes when:
    //   low_norm_ . (point - (0,0,-dz)) <= 0  (inside the low plane)
    //   high_norm_ . (point - (0,0,+dz)) <= 0 (inside the high plane)

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dirx = direction.GetX();
    double diry = direction.GetY();
    double dirz = direction.GetZ();

    // Origin shift: move ray origin to closest approach to coordinate origin.
    // This keeps quadratic barrel coefficients small for far-field rays.
    double t_shift = -(px * dirx + py * diry + pz * dirz);
    double qx = px + t_shift * dirx;
    double qy = py + t_shift * diry;
    double qz = pz + t_shift * dirz;

    double r2_outer = rmax_ * rmax_;
    double r2_inner = rmin_ * rmin_;

    double lnx = low_norm_.GetX();
    double lny = low_norm_.GetY();
    double lnz = low_norm_.GetZ();
    double hnx = high_norm_.GetX();
    double hny = high_norm_.GetY();
    double hnz = high_norm_.GetZ();

    // Max 8 intersections: 2 outer barrel + 2 inner barrel + 2 low cap + 2 high cap
    Intersection hits[8];
    int n_hits = 0;

    auto add_hit = [&](double t, double ix, double iy, double iz, bool entering) {
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        hits[n_hits].distance = t;
        hits[n_hits].hierarchy = 0;
        hits[n_hits].entering = entering;
        hits[n_hits].position = siren::math::Vector3D(ix, iy, iz);
        n_hits++;
    };

    // Check whether a point is between both cut planes (on the interior side).
    // low_norm_ . (P - (0,0,-dz)) <= 0  AND  high_norm_ . (P - (0,0,+dz)) <= 0
    auto between_planes = [&](double x, double y, double z) -> bool {
        double low_dot = lnx * x + lny * y + lnz * (z + dz_);
        double high_dot = hnx * x + hny * y + hnz * (z - dz_);
        return (low_dot <= GEOMETRY_PRECISION) && (high_dot <= GEOMETRY_PRECISION);
    };

    // ------ Barrel surfaces ------ //
    // Quadratic in shifted frame: A*s^2 + 2*B_half*s + C = 0
    // Uses shifted origin (q) for far-field numerical stability.
    auto test_barrel = [&](double r2, bool invert_entering) {
        double A = dirx * dirx + diry * diry;
        if(A == 0) return;
        double B_half = qx * dirx + qy * diry;
        double C = qx * qx + qy * qy - r2;
        double det = B_half * B_half - A * C;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double inv_A = 1.0 / A;
        double t1 = (-B_half - sq) * inv_A + t_shift;
        double t2 = (-B_half + sq) * inv_A + t_shift;
        for(int k = 0; k < 2; ++k) {
            double t = (k == 0) ? t1 : t2;
            double ix = px + t * dirx;
            double iy = py + t * diry;
            double iz = pz + t * dirz;
            if(between_planes(ix, iy, iz)) {
                // Radial entering: ray moves toward the axis
                bool radial_entering = (ix * dirx + iy * diry) < 0;
                add_hit(t, ix, iy, iz, invert_entering ? !radial_entering : radial_entering);
            }
        }
    };

    // ------ Cut plane surfaces ------ //
    // Plane through center_point with outward normal n:
    //   n . (ray(t) - center) = 0
    //   t = -n . (origin - center) / (n . dir)
    // At the hit point, check annular ring: rmin^2 <= ix^2 + iy^2 <= rmax^2
    // Entering: n . dir < 0 (ray going into the solid through the cap)
    auto test_cap = [&](double cx, double cy, double cz,
                        double nx, double ny, double nz) {
        double n_dot_dir = nx * dirx + ny * diry + nz * dirz;
        if(std::fabs(n_dot_dir) < GEOMETRY_PRECISION) return;
        // t = -n . (origin - center) / (n . dir)
        double n_dot_off = nx * (px - cx) + ny * (py - cy) + nz * (pz - cz);
        double t = -n_dot_off / n_dot_dir;
        double ix = px + t * dirx;
        double iy = py + t * diry;
        double iz = pz + t * dirz;
        double r2_hit = ix * ix + iy * iy;
        if(r2_hit >= r2_inner - GEOMETRY_PRECISION && r2_hit <= r2_outer + GEOMETRY_PRECISION) {
            bool entering = n_dot_dir < 0;
            add_hit(t, ix, iy, iz, entering);
        }
    };

    // Outer barrel
    test_barrel(r2_outer, false);
    // Inner barrel (if hollow)
    if(rmin_ > 0) {
        test_barrel(r2_inner, true);
    }
    // Low cut plane: center at (0, 0, -dz_), normal = low_norm_
    test_cap(0, 0, -dz_, lnx, lny, lnz);
    // High cut plane: center at (0, 0, +dz_), normal = high_norm_
    test_cap(0, 0, dz_, hnx, hny, hnz);

    if(n_hits == 0 && !has_phi_cut_) return {};

    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });

    if(!has_phi_cut_) {
        return {hits, hits + n_hits};
    }

    // Phi cut: merge surface hits with infinite wedge hits and run CSG walk.
    // See Polycone.cxx for method description; same pattern in all phi-cut shapes.
    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    // Max 6 tagged hits: 4 surface (2 barrel + 2 cap) + 2 wedge
    TaggedHit all_hits[10];
    int n_all = 0;

    for(int i = 0; i < n_hits; ++i) {
        all_hits[n_all] = {hits[i].distance, hits[i].position, hits[i].entering, 0};
        n_all++;
    }

    // Compute infinite wedge intersections (two half-planes from z-axis)
    for(int face = 0; face < 2; ++face) {
        double alpha = start_phi_ + face * delta_phi_;
        double ca = std::cos(alpha), sa = std::sin(alpha);
        // Outward-pointing normal (away from phi range interior)
        double nx, ny;
        if(face == 0) { nx = sa; ny = -ca; }
        else { nx = -sa; ny = ca; }
        double n_dot_d = nx*dirx + ny*diry;
        if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue;
        double n_dot_p = nx*px + ny*py;
        double t = -n_dot_p / n_dot_d;
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;

        double hx = px + t*dirx, hy = py + t*diry, hz = pz + t*dirz;
        // Must be on the correct half-plane (outward from z-axis)
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;

        bool entering = (n_dot_d < 0);
        all_hits[n_all] = {t, siren::math::Vector3D(hx, hy, hz), entering, 1};
        n_all++;
    }

    if(n_all == 0) return {};

    // Sort all hits by distance
    std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // CSG intersection walk: the phi-cut solid is (full surface) AND (phi wedge).
    bool in_surface = false;
    bool in_wedge = false;

    // Determine initial in_wedge state
    bool has_wedge_hit = false;
    for(int i = 0; i < n_all; ++i) {
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
    for(int i = 0; i < n_all; ++i) {
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
AABB CutTube::GetBoundingBox() const {
    // Conservative bounding box.
    // z extents: the cut planes can extend the z range beyond +/-dz.
    // At the barrel circle (radius rmax), the most extreme z-intercept of
    // the low plane is at the point on the barrel circle where the
    // horizontal component of the normal is aligned:
    //   z_min = -dz - rmax * sqrt(lnx^2 + lny^2) / |lnz|
    //   z_max =  dz + rmax * sqrt(hnx^2 + hny^2) / |hnz|
    double lnx = low_norm_.GetX();
    double lny = low_norm_.GetY();
    double lnz = low_norm_.GetZ();
    double hnx = high_norm_.GetX();
    double hny = high_norm_.GetY();
    double hnz = high_norm_.GetZ();

    double low_tilt = std::sqrt(lnx * lnx + lny * lny) / std::fabs(lnz);
    double high_tilt = std::sqrt(hnx * hnx + hny * hny) / std::fabs(hnz);

    double z_min = -dz_ - rmax_ * low_tilt;
    double z_max =  dz_ + rmax_ * high_tilt;

    if(!has_phi_cut_) {
        return AABB(
            math::Vector3D(-rmax_, -rmax_, z_min),
            math::Vector3D( rmax_,  rmax_, z_max)
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
        double cx = std::cos(a) * rmax_;
        double cy = std::sin(a) * rmax_;
        if(cx < x_min) x_min = cx;
        if(cx > x_max) x_max = cx;
        if(cy < y_min) y_min = cy;
        if(cy > y_max) y_max = cy;
    }

    // Check cardinal directions if they fall in the sector
    double cardinals[4] = {0.0, M_PI / 2.0, M_PI, 3.0 * M_PI / 2.0};
    for(int i = 0; i < 4; ++i) {
        double c = cardinals[i];
        double cn = c;
        if(cn < sp) cn += TWO_PI;
        if(cn <= ep + 1e-9) {
            double cx = std::cos(cardinals[i]) * rmax_;
            double cy = std::sin(cardinals[i]) * rmax_;
            if(cx < x_min) x_min = cx;
            if(cx > x_max) x_max = cx;
            if(cy < y_min) y_min = cy;
            if(cy > y_max) y_max = cy;
        }
    }

    return AABB(
        math::Vector3D(x_min, y_min, z_min),
        math::Vector3D(x_max, y_max, z_max)
    );
}

} // namespace geometry
} // namespace siren
