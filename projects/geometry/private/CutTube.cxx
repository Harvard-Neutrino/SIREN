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
      , high_norm_(siren::math::Vector3D(0, 0, 1)) {
    // Do nothing here
}

CutTube::CutTube(double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm)
    : Geometry((std::string)("CutTube"))
    , rmin_(rmin)
    , rmax_(rmax)
    , dz_(dz)
    , low_norm_(low_norm)
      , high_norm_(high_norm) {
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
      , high_norm_(siren::math::Vector3D(0, 0, 1)) {
    // Do nothing here
}

CutTube::CutTube(Placement const & placement, double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm)
    : Geometry((std::string)("CutTube"), placement)
    , rmin_(rmin)
    , rmax_(rmax)
    , dz_(dz)
    , low_norm_(low_norm)
      , high_norm_(high_norm) {
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

CutTube::CutTube(const CutTube& other)
    : Geometry(other)
    , rmin_(other.rmin_)
    , rmax_(other.rmax_)
    , dz_(other.dz_)
    , low_norm_(other.low_norm_)
      , high_norm_(other.high_norm_) {
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
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool CutTube::less(const Geometry& geometry) const
{
    const CutTube* ct = dynamic_cast<const CutTube*>(&geometry);
    if(!ct) return false;

    return
        std::tie(rmin_, rmax_, dz_, low_norm_, high_norm_)
        <
        std::tie(ct->rmin_, ct->rmax_, ct->dz_, ct->low_norm_, ct->high_norm_);
}

// ------------------------------------------------------------------------- //
void CutTube::print(std::ostream& os) const
{
    os << "rmin: " << rmin_ << "\trmax: " << rmax_
       << "\tdz: " << dz_
       << "\tlow_norm: " << low_norm_
       << "\thigh_norm: " << high_norm_ << '\n';
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
    // Quadratic: A*t^2 + 2*B_half*t + C = 0
    // A = dirx^2 + diry^2
    // B_half = px*dirx + py*diry
    // C = px^2 + py^2 - r^2
    auto test_barrel = [&](double r2, bool invert_entering) {
        double A = dirx * dirx + diry * diry;
        if(A == 0) return;
        double B_half = px * dirx + py * diry;
        double C = px * px + py * py - r2;
        double det = B_half * B_half - A * C;
        if(det <= 0) return;
        double sq = std::sqrt(det);
        double inv_A = 1.0 / A;
        double t1 = (-B_half - sq) * inv_A;
        double t2 = (-B_half + sq) * inv_A;
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

    if(n_hits == 0) return {};
    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    return {hits, hits + n_hits};
}

// ------------------------------------------------------------------------- //
AABB CutTube::GetBoundingBox() const {
    // Conservative bounding box.
    // x, y extents: [-rmax, rmax]
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

    return AABB(
        math::Vector3D(-rmax_, -rmax_, z_min),
        math::Vector3D( rmax_,  rmax_, z_max)
    );
}

} // namespace geometry
} // namespace siren
