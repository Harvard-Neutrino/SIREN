#include "SIREN/geometry/ExtrPoly.h"

#include <cmath>
#include <tuple>
#include <limits>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

ExtrPoly::ExtrPoly()
    : Geometry((std::string)("ExtrPoly"))
      , polygon_({})
, zsections_({})
{
    // Do nothing here
}


ExtrPoly::ExtrPoly(const std::vector<std::vector<double>>& polygon,
        const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry("ExtrPoly")
    , polygon_(polygon)
      , zsections_(zsections)
{
    if (polygon.size() < 3)
    {
        std::cout << "Need 3 polygon vertices at least!! Give it another shot";
        return;
    }
    ComputeLateralPlanes();
}

ExtrPoly::ExtrPoly(Placement const & placement)
    : Geometry((std::string)("ExtrPoly"), placement)
      , polygon_({})
, zsections_({})
{
    ComputeLateralPlanes();
}


ExtrPoly::ExtrPoly(Placement const & placement, const std::vector<std::vector<double>>& polygon,
        const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry((std::string)("ExtrPoly"), placement)
    , polygon_(polygon)
      , zsections_(zsections)
{
    if (polygon.size() < 3)
    {
        std::cout << "Need 3 polygon vertices at least!! Give it another shot";
        return;
    }
    ComputeLateralPlanes();
}

ExtrPoly::ExtrPoly(const ExtrPoly& extr)
    : Geometry(extr)
    , polygon_(extr.polygon_)
      , zsections_(extr.zsections_)
{
    ComputeLateralPlanes();
}


// ------------------------------------------------------------------------- //
void ExtrPoly::swap(Geometry& geometry)
{
    ExtrPoly* extr = dynamic_cast<ExtrPoly*>(&geometry);
    if (!extr)
    {
        //log_warn("Cannot swap ExtrPoly!");
        return;
    }

    Geometry::swap(*extr);

    std::swap(polygon_, extr->polygon_);
    std::swap(zsections_, extr->zsections_);
}

//------------------------------------------------------------------------- //
ExtrPoly& ExtrPoly::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);
        if (!extr)
        {
            //log_warn("Cannot assign ExtrPoly!");
            return *this;
        }

        ExtrPoly tmp(*extr);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::equal(const Geometry& geometry) const
{
    const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);

    if (!extr)
        return false;
    else if (polygon_ != extr->polygon_)
        return false;
    else if (!(zsections_ == extr->zsections_))
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::less(const Geometry& geometry) const
{
    const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);
    if(!extr) return false;

    return
        std::tie(polygon_, zsections_)
        <
        std::tie(extr->polygon_, extr->zsections_);
}

// ------------------------------------------------------------------------- //
void ExtrPoly::print(std::ostream& os) const
{
    os << "Not implemented yet\n";
    //os << "Polygon: " << polygon_ << "\tZsections: " << zsections_ << '\n';
}

// ------------------------------------------------------------------------- //
void ExtrPoly::ComputeLateralPlanes()
{
    int Nv = polygon_.size();
    if(Nv < 3) {
        planes_.clear();
        return;
    }
    planes_.resize(Nv);

    // Compute signed area to determine winding direction.
    // Positive = CCW, negative = CW.
    double signed_area = 0;
    for(int i = 0, k = Nv - 1; i < Nv; k = i++) {
        signed_area += polygon_[k][0] * polygon_[i][1]
                     - polygon_[i][0] * polygon_[k][1];
    }
    // sign = +1 for CW (outward normals from cross product),
    //       -1 for CCW (flip to make outward)
    double sign = (signed_area < 0) ? 1.0 : -1.0;

    for(int i = 0, k = Nv - 1; i < Nv; k = i++) {
        double ex = polygon_[i][0] - polygon_[k][0];
        double ey = polygon_[i][1] - polygon_[k][1];
        double norm = std::sqrt(ex * ex + ey * ey);
        if(norm < 1e-15) norm = 1e-15;
        double inv_norm = sign / norm;
        planes_[i].a = -ey * inv_norm;
        planes_[i].b =  ex * inv_norm;
        planes_[i].c = 0;
        planes_[i].d = (ey * polygon_[i][0] - ex * polygon_[i][1]) * inv_norm;
    }
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> ExtrPoly::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Full-line slab intersection with an extruded convex polygon.
    //
    // Between each pair of adjacent z-sections, the cross-section is a convex
    // polygon whose vertices are scaled and offset: v' = v * scale(z) + offset(z).
    // The lateral face plane coefficients vary linearly with z (and hence with
    // ray parameter t), giving a linear equation for each face intersection.

    int Nz = zsections_.size();
    if(Nz < 2 || planes_.empty()) return {};

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    int np = planes_.size();

    Intersection hits[64];
    int n_hits = 0;

    for(int k = 0; k + 1 < Nz; ++k) {
        double zk = zsections_[k].zpos;
        double zk1 = zsections_[k + 1].zpos;
        double dzk = zk1 - zk;
        if(dzk <= 0) continue;

        double t_enter = -std::numeric_limits<double>::infinity();
        double t_exit = std::numeric_limits<double>::infinity();

        // Z-slab for this section
        if(std::fabs(dz) > GEOMETRY_PRECISION) {
            double t_z0 = (zk - pz) / dz;
            double t_z1 = (zk1 - pz) / dz;
            if(t_z0 > t_z1) std::swap(t_z0, t_z1);
            t_enter = std::max(t_enter, t_z0);
            t_exit = std::min(t_exit, t_z1);
        } else {
            if(pz < zk || pz > zk1) continue;
        }

        // Scale and offset at the ray origin's z, interpolated within this
        // section. The effective plane constant at parameter t is:
        //   D(t) = scale(z(t)) * d - a * offset_x(z(t)) - b * offset_y(z(t))
        // which is linear in t.
        double inv_dzk = 1.0 / dzk;
        double frac_pz = (pz - zk) * inv_dzk;

        double sk = zsections_[k].scale;
        double sk1 = zsections_[k + 1].scale;
        double oxk = zsections_[k].offset[0];
        double oxk1 = zsections_[k + 1].offset[0];
        double oyk = zsections_[k].offset[1];
        double oyk1 = zsections_[k + 1].offset[1];

        // Linear coefficients: value(t) = val0 + val1 * t
        double s0 = sk + (sk1 - sk) * frac_pz;
        double s1 = (sk1 - sk) * inv_dzk * dz;
        double ox0 = oxk + (oxk1 - oxk) * frac_pz;
        double ox1 = (oxk1 - oxk) * inv_dzk * dz;
        double oy0 = oyk + (oyk1 - oyk) * frac_pz;
        double oy1 = (oyk1 - oyk) * inv_dzk * dz;

        bool missed = false;
        for(int i = 0; i < np; ++i) {
            double a = planes_[i].a;
            double b = planes_[i].b;
            double d = planes_[i].d;

            // Effective plane equation along the ray:
            //   a*(px + t*dx - ox(t)) + b*(py + t*dy - oy(t)) + d*s(t) = 0
            // => cosa*t + dist = 0
            // Normals point OUTWARD, so interior is dist + cosa*t < 0.
            double cosa = a * (dx - ox1) + b * (dy - oy1) + d * s1;
            double dist = a * (px - ox0) + b * (py - oy0) + d * s0;

            if(std::fabs(cosa) > GEOMETRY_PRECISION) {
                double t = -dist / cosa;
                if(cosa < 0) {
                    // Interior for t > -dist/cosa (entering)
                    t_enter = std::max(t_enter, t);
                } else {
                    // Interior for t < -dist/cosa (exiting)
                    t_exit = std::min(t_exit, t);
                }
            } else {
                if(dist >= 0) {
                    // Entirely outside this half-plane
                    missed = true;
                    break;
                }
            }
        }

        if(missed || t_enter >= t_exit - GEOMETRY_PRECISION) continue;

        if(n_hits + 2 <= 64) {
            hits[n_hits].distance = t_enter;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = true;
            hits[n_hits].position = siren::math::Vector3D(
                px + dx * t_enter, py + dy * t_enter, pz + dz * t_enter);
            n_hits++;

            hits[n_hits].distance = t_exit;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = false;
            hits[n_hits].position = siren::math::Vector3D(
                px + dx * t_exit, py + dy * t_exit, pz + dz * t_exit);
            n_hits++;
        }
    }

    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    return {hits, hits + n_hits};
}

// ------------------------------------------------------------------------- //
AABB ExtrPoly::GetBoundingBox() const {
    AABB box;
    // Iterate over all z-sections, applying offset and scale to each polygon vertex
    for(auto const & zsec : zsections_) {
        for(auto const & vert : polygon_) {
            if(vert.size() >= 2) {
                double x = vert[0] * zsec.scale + zsec.offset[0];
                double y = vert[1] * zsec.scale + zsec.offset[1];
                box.ExpandToInclude(math::Vector3D(x, y, zsec.zpos));
            }
        }
    }
    return box;
}

} // namespace geometry
} // namespace siren
