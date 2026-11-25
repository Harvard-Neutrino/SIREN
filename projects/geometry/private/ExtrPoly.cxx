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
#include "GeometryMacros.h"

namespace siren {
namespace geometry {

namespace {

// A unit ray moving sideways slower than this is along the extrusion.
constexpr double PARALLEL_PRECISION = 1e-14;
// Below this rate of approach a lateral plane fixes the crossing poorly.
constexpr double CROSSING_PRECISION = 1e-6;

} // namespace

ExtrPoly::ExtrPoly() : Geometry("ExtrPoly"), polygon_({}), zsections_({}) { RecomputeWorldAABB(); }
ExtrPoly::ExtrPoly(const std::vector<std::vector<double>>& polygon, const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry("ExtrPoly"), polygon_(polygon), zsections_(zsections) {
    if(polygon.size() < 3) { throw std::invalid_argument("ExtrPoly requires at least 3 polygon vertices"); }
    ComputeLateralPlanes();
    RecomputeWorldAABB();
}
ExtrPoly::ExtrPoly(Placement const & p) : Geometry("ExtrPoly", p), polygon_({}), zsections_({}) { ComputeLateralPlanes(); RecomputeWorldAABB(); }
ExtrPoly::ExtrPoly(Placement const & p, const std::vector<std::vector<double>>& polygon, const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry("ExtrPoly", p), polygon_(polygon), zsections_(zsections) {
    if(polygon.size() < 3) { throw std::invalid_argument("ExtrPoly requires at least 3 polygon vertices"); }
    ComputeLateralPlanes();
    RecomputeWorldAABB();
}
ExtrPoly::ExtrPoly(const ExtrPoly& o) : Geometry(o), polygon_(o.polygon_), zsections_(o.zsections_) { ComputeLateralPlanes(); RecomputeWorldAABB(); }

SIREN_GEOMETRY_SWAP(ExtrPoly, polygon_, zsections_, planes_)
SIREN_GEOMETRY_ASSIGN(ExtrPoly)
SIREN_GEOMETRY_EQUAL(ExtrPoly, polygon_, zsections_)
SIREN_GEOMETRY_LESS(ExtrPoly, polygon_, zsections_)

void ExtrPoly::print(std::ostream& os) const {
    os << "ExtrPoly(nvert=" << polygon_.size() << ", nsec=" << zsections_.size() << ")\n";
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
    // Full-line intersection with an extruded polygon, convex or not.
    //
    // Between adjacent z-sections the cross-section is the polygon scaled and
    // offset: v' = v * scale(z) + offset(z). Undoing that map takes the ray to
    // the line U(t) = (x - offset_x, y - offset_y, scale) in the homogeneous
    // coordinates of the polygon plane, with every component linear in t.
    //
    // Each vertex is assigned to one side of that line, and an outline edge is
    // crossed exactly when its end points are on different sides. The side is
    // a function of the vertex alone, so a closed outline always has an even
    // number of crossings, including for rays through vertices, and
    // consecutive crossings bound the stretches inside the outline.

    int Nz = zsections_.size();
    int Nv = planes_.size();
    if(Nz < 2 || Nv < 3) return {};

    double const inf = std::numeric_limits<double>::infinity();

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    // The polygon plane is taken relative to the first vertex, so an outline far
    // from its origin loses nothing to cancellation.
    double ref_x = polygon_[0][0];
    double ref_y = polygon_[0][1];

    // Pieces of the ray inside the solid, over all z-sections
    std::vector<std::pair<double, double>> pieces;
    std::vector<double> crossings;

    for(int k = 0; k + 1 < Nz; ++k) {
        double zk = zsections_[k].zpos;
        double zk1 = zsections_[k + 1].zpos;
        double dzk = zk1 - zk;
        if(dzk <= 0) continue;

        // Z-slab for this section
        double t_min = -inf;
        double t_max = inf;
        double t_mid = 0;
        if(std::fabs(dz) > GEOMETRY_PRECISION) {
            t_min = (zk - pz) / dz;
            t_max = (zk1 - pz) / dz;
            if(t_min > t_max) std::swap(t_min, t_max);
            t_mid = 0.5 * (t_min + t_max);
        } else {
            if(pz < zk || pz > zk1) continue;
        }

        double inv_dzk = 1.0 / dzk;
        double frac_pz = (pz - zk) * inv_dzk;

        double sk = zsections_[k].scale;
        double sk1 = zsections_[k + 1].scale;
        double oxk = zsections_[k].offset[0] + sk * ref_x;
        double oxk1 = zsections_[k + 1].offset[0] + sk1 * ref_x;
        double oyk = zsections_[k].offset[1] + sk * ref_y;
        double oyk1 = zsections_[k + 1].offset[1] + sk1 * ref_y;

        // U(t) = U0 + U1 * t
        double w0 = sk + (sk1 - sk) * frac_pz;
        double w1 = (sk1 - sk) * inv_dzk * dz;
        double ux0 = px - (oxk + (oxk1 - oxk) * frac_pz);
        double ux1 = dx - (oxk1 - oxk) * inv_dzk * dz;
        double uy0 = py - (oyk + (oyk1 - oyk) * frac_pz);
        double uy1 = dy - (oyk1 - oyk) * inv_dzk * dz;

        double w_mid = w0 + w1 * t_mid;
        if(!(w_mid > 0)) continue;

        // The line through U0 and U1 is nx * x + ny * y + nw = 0.
        double nx = uy0 * w1 - w0 * uy1;
        double ny = w0 * ux1 - ux0 * w1;
        double nw = ux0 * uy1 - uy0 * ux1;

        // The ray moves sideways at |n| / w_mid in the middle of the slab. Like n,
        // the bounds below scale with w, so the outline and the sections may be
        // rescaled reciprocally. The second terms are the rounding error of n.
        double eps = 4 * std::numeric_limits<double>::epsilon();
        double tol = PARALLEL_PRECISION * w_mid;
        bool parallel =
            std::fabs(nx) <= tol + eps * (std::fabs(uy0 * w1) + std::fabs(w0 * uy1)) &&
            std::fabs(ny) <= tol + eps * (std::fabs(w0 * ux1) + std::fabs(ux0 * w1));
        if(parallel && std::isfinite(t_min)) {
            // Along the extrusion: a single point in the polygon plane.
            if(InsidePolygon((ux0 + ux1 * t_mid) / w_mid, (uy0 + uy1 * t_mid) / w_mid)) {
                pieces.emplace_back(t_min, t_max);
            }
            continue;
        }

        // Fix the sign of the line so that the side of a vertex exactly on it
        // does not depend on the direction of travel.
        if(nx < 0 || (nx == 0 && ny < 0)) {
            nx = -nx;
            ny = -ny;
            nw = -nw;
        }

        crossings.clear();
        double xb = polygon_[Nv - 1][0] - ref_x;
        double yb = polygon_[Nv - 1][1] - ref_y;
        double side_b = nx * xb + ny * yb + nw;
        for(int i = 0; i < Nv; ++i) {
            double xa = xb, ya = yb, side_a = side_b;
            xb = polygon_[i][0] - ref_x;
            yb = polygon_[i][1] - ref_y;
            side_b = nx * xb + ny * yb + nw;
            if((side_a >= 0) == (side_b >= 0)) continue;

            // Lateral plane along the ray: dist + cosa * t = 0
            plane const & pl = planes_[i];
            double cosa = pl.a * (ux1 - xb * w1) + pl.b * (uy1 - yb * w1);
            double t;
            if(std::fabs(cosa) > CROSSING_PRECISION) {
                t = -(pl.a * (ux0 - xb * w0) + pl.b * (uy0 - yb * w0)) / cosa;
            } else {
                // Nearly in the plane of the face, which then fixes t poorly:
                // take the point where the edge meets the line instead.
                double f = side_a / (side_a - side_b);
                double x = xa + f * (xb - xa);
                double y = ya + f * (yb - ya);
                double ax = ux1 - x * w1;
                double ay = uy1 - y * w1;
                t = (ax * (x * w0 - ux0) + ay * (y * w0 - uy0)) / (ax * ax + ay * ay);
            }
            // The section only exists where the scale is positive. Other
            // crossings lie past the apex that the scale extrapolates to.
            if(w0 + w1 * t > 0) crossings.push_back(t);
        }

        // Pair the crossings starting from the apex, which is outside. On the
        // far end the line may stay inside forever; the z-slab bounds it.
        std::sort(crossings.begin(), crossings.end());
        if(crossings.size() % 2 == 1) {
            if(w1 < 0) crossings.insert(crossings.begin(), -inf);
            else crossings.push_back(inf);
        }
        for(size_t i = 0; i + 1 < crossings.size(); i += 2) {
            double lo = std::max(crossings[i], t_min);
            double hi = std::min(crossings[i + 1], t_max);
            if(lo < hi) pieces.emplace_back(lo, hi);
        }
    }

    // Join pieces that overlap or touch, which happens where adjacent
    // sections meet, and drop what is too short to resolve.
    std::sort(pieces.begin(), pieces.end());
    std::vector<Intersection> hits;
    auto add_hit = [&](double t, bool entering) {
        Intersection isect;
        isect.distance = t;
        isect.hierarchy = 0;
        isect.entering = entering;
        isect.matID = 0;
        isect.position = siren::math::Vector3D(px + dx * t, py + dy * t, pz + dz * t);
        hits.push_back(isect);
    };
    for(size_t i = 0; i < pieces.size(); ) {
        double lo = pieces[i].first;
        double hi = pieces[i].second;
        for(++i; i < pieces.size() && pieces[i].first <= hi + GEOMETRY_PRECISION; ++i) {
            hi = std::max(hi, pieces[i].second);
        }
        if(hi - lo <= GEOMETRY_PRECISION) continue;
        add_hit(lo, true);
        add_hit(hi, false);
    }
    return hits;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::InsidePolygon(double x, double y) const {
    // Even-odd rule on half-open edges, so that a vertex is counted once.
    // The point is relative to the first vertex.
    bool inside = false;
    int Nv = polygon_.size();
    for(int i = 0, k = Nv - 1; i < Nv; k = i++) {
        double xi = polygon_[i][0] - polygon_[0][0], yi = polygon_[i][1] - polygon_[0][1];
        double xk = polygon_[k][0] - polygon_[0][0], yk = polygon_[k][1] - polygon_[0][1];
        if((yi > y) == (yk > y)) continue;
        if(x < xi + (y - yi) * (xk - xi) / (yk - yi)) inside = !inside;
    }
    return inside;
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
