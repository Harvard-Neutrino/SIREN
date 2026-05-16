#include "SIREN/geometry/Trap.h"

#include <cmath>
#include <tuple>
#include <math.h>
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

Trap::Trap()
    : Geometry((std::string)("Trap"))
    , dz_(0.0)
    , theta_(0.0)
    , phi_(0.0)
    , dy1_(0.0)
    , dx1_(0.0)
    , dx2_(0.0)
    , alpha1_(0.0)
    , dy2_(0.0)
    , dx3_(0.0)
    , dx4_(0.0)
      , alpha2_(0.0) {
    ComputeVertices();
    ComputePlanes();
}

Trap::Trap(double dz, double theta, double phi,
           double dy1, double dx1, double dx2, double alpha1,
           double dy2, double dx3, double dx4, double alpha2)
    : Geometry((std::string)("Trap"))
    , dz_(dz)
    , theta_(theta)
    , phi_(phi)
    , dy1_(dy1)
    , dx1_(dx1)
    , dx2_(dx2)
    , alpha1_(alpha1)
    , dy2_(dy2)
    , dx3_(dx3)
    , dx4_(dx4)
      , alpha2_(alpha2) {
    if(dz_ <= 0) {
        throw std::invalid_argument("Trap half-height dz must be positive!");
    }
    if(dx1_ < 0 || dx2_ < 0 || dx3_ < 0 || dx4_ < 0) {
        throw std::invalid_argument("Trap half-widths (dx1, dx2, dx3, dx4) must be non-negative!");
    }
    if(dy1_ < 0 || dy2_ < 0) {
        throw std::invalid_argument("Trap half-heights (dy1, dy2) must be non-negative!");
    }
    ComputeVertices();
    ComputePlanes();
}

Trap::Trap(Placement const & placement)
    : Geometry((std::string)("Trap"), placement)
    , dz_(0.0)
    , theta_(0.0)
    , phi_(0.0)
    , dy1_(0.0)
    , dx1_(0.0)
    , dx2_(0.0)
    , alpha1_(0.0)
    , dy2_(0.0)
    , dx3_(0.0)
    , dx4_(0.0)
      , alpha2_(0.0) {
    ComputeVertices();
    ComputePlanes();
}

Trap::Trap(Placement const & placement,
           double dz, double theta, double phi,
           double dy1, double dx1, double dx2, double alpha1,
           double dy2, double dx3, double dx4, double alpha2)
    : Geometry((std::string)("Trap"), placement)
    , dz_(dz)
    , theta_(theta)
    , phi_(phi)
    , dy1_(dy1)
    , dx1_(dx1)
    , dx2_(dx2)
    , alpha1_(alpha1)
    , dy2_(dy2)
    , dx3_(dx3)
    , dx4_(dx4)
      , alpha2_(alpha2) {
    if(dz_ <= 0) {
        throw std::invalid_argument("Trap half-height dz must be positive!");
    }
    if(dx1_ < 0 || dx2_ < 0 || dx3_ < 0 || dx4_ < 0) {
        throw std::invalid_argument("Trap half-widths (dx1, dx2, dx3, dx4) must be non-negative!");
    }
    if(dy1_ < 0 || dy2_ < 0) {
        throw std::invalid_argument("Trap half-heights (dy1, dy2) must be non-negative!");
    }
    ComputeVertices();
    ComputePlanes();
}

Trap::Trap(const Trap& trap)
    : Geometry(trap)
    , dz_(trap.dz_)
    , theta_(trap.theta_)
    , phi_(trap.phi_)
    , dy1_(trap.dy1_)
    , dx1_(trap.dx1_)
    , dx2_(trap.dx2_)
    , alpha1_(trap.alpha1_)
    , dy2_(trap.dy2_)
    , dx3_(trap.dx3_)
    , dx4_(trap.dx4_)
      , alpha2_(trap.alpha2_) {
    for(int i = 0; i < 8; ++i) {
        for(int j = 0; j < 3; ++j) {
            vertices_[i][j] = trap.vertices_[i][j];
        }
    }
    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j < 4; ++j) {
            planes_[i][j] = trap.planes_[i][j];
        }
    }
}

// ------------------------------------------------------------------------- //
void Trap::ComputeVertices() {
    double tan_theta = std::tan(theta_);
    double cos_phi = std::cos(phi_);
    double sin_phi = std::sin(phi_);
    double tan_alpha1 = std::tan(alpha1_);
    double tan_alpha2 = std::tan(alpha2_);

    // Bottom face at z = -dz
    double xc_lo = -dz_ * tan_theta * cos_phi;
    double yc_lo = -dz_ * tan_theta * sin_phi;

    // v0: bottom-left (y = -dy1)
    vertices_[0][0] = xc_lo - dy1_ * tan_alpha1 - dx1_;
    vertices_[0][1] = yc_lo - dy1_;
    vertices_[0][2] = -dz_;

    // v1: bottom-right (y = -dy1)
    vertices_[1][0] = xc_lo - dy1_ * tan_alpha1 + dx1_;
    vertices_[1][1] = yc_lo - dy1_;
    vertices_[1][2] = -dz_;

    // v2: top-right (y = +dy1)
    vertices_[2][0] = xc_lo + dy1_ * tan_alpha1 + dx2_;
    vertices_[2][1] = yc_lo + dy1_;
    vertices_[2][2] = -dz_;

    // v3: top-left (y = +dy1)
    vertices_[3][0] = xc_lo + dy1_ * tan_alpha1 - dx2_;
    vertices_[3][1] = yc_lo + dy1_;
    vertices_[3][2] = -dz_;

    // Top face at z = +dz
    double xc_hi = dz_ * tan_theta * cos_phi;
    double yc_hi = dz_ * tan_theta * sin_phi;

    // v4: bottom-left (y = -dy2)
    vertices_[4][0] = xc_hi - dy2_ * tan_alpha2 - dx3_;
    vertices_[4][1] = yc_hi - dy2_;
    vertices_[4][2] = dz_;

    // v5: bottom-right (y = -dy2)
    vertices_[5][0] = xc_hi - dy2_ * tan_alpha2 + dx3_;
    vertices_[5][1] = yc_hi - dy2_;
    vertices_[5][2] = dz_;

    // v6: top-right (y = +dy2)
    vertices_[6][0] = xc_hi + dy2_ * tan_alpha2 + dx4_;
    vertices_[6][1] = yc_hi + dy2_;
    vertices_[6][2] = dz_;

    // v7: top-left (y = +dy2)
    vertices_[7][0] = xc_hi + dy2_ * tan_alpha2 - dx4_;
    vertices_[7][1] = yc_hi + dy2_;
    vertices_[7][2] = dz_;
}

// ------------------------------------------------------------------------- //
void Trap::ComputePlane(const double v0[3], const double v1[3], const double v2[3], double plane[4]) const {
    // edge1 = v1 - v0
    double e1x = v1[0] - v0[0];
    double e1y = v1[1] - v0[1];
    double e1z = v1[2] - v0[2];

    // edge2 = v2 - v0
    double e2x = v2[0] - v0[0];
    double e2y = v2[1] - v0[1];
    double e2z = v2[2] - v0[2];

    // normal = cross(edge1, edge2)
    double nx = e1y * e2z - e1z * e2y;
    double ny = e1z * e2x - e1x * e2z;
    double nz = e1x * e2y - e1y * e2x;

    // normalize
    double len = std::sqrt(nx * nx + ny * ny + nz * nz);
    if(len > 0.0) {
        nx /= len;
        ny /= len;
        nz /= len;
    }

    plane[0] = nx;
    plane[1] = ny;
    plane[2] = nz;
    plane[3] = -(nx * v0[0] + ny * v0[1] + nz * v0[2]);
}

// ------------------------------------------------------------------------- //
void Trap::ComputePlanes() {
    // Compute the centroid of the solid for outward-normal orientation
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for(int i = 0; i < 8; ++i) {
        cx += vertices_[i][0];
        cy += vertices_[i][1];
        cz += vertices_[i][2];
    }
    cx /= 8.0;
    cy /= 8.0;
    cz /= 8.0;

    // Face 0 (bottom, z=-dz): vertices 0,1,2,3
    ComputePlane(vertices_[0], vertices_[1], vertices_[2], planes_[0]);

    // Face 1 (top, z=+dz): vertices 7,6,5,4 (reversed winding for outward normal toward +z)
    ComputePlane(vertices_[7], vertices_[6], vertices_[5], planes_[1]);

    // Face 2 (front, y-low): vertices 0,1,5,4
    ComputePlane(vertices_[0], vertices_[1], vertices_[5], planes_[2]);

    // Face 3 (back, y-high): vertices 3,2,6,7 -- use 3,2,6
    ComputePlane(vertices_[3], vertices_[2], vertices_[6], planes_[3]);

    // Face 4 (left, x-low): vertices 0,3,7,4
    ComputePlane(vertices_[0], vertices_[3], vertices_[7], planes_[4]);

    // Face 5 (right, x-high): vertices 1,2,6,5
    ComputePlane(vertices_[1], vertices_[2], vertices_[6], planes_[5]);

    // Ensure all normals point outward (away from centroid)
    for(int f = 0; f < 6; ++f) {
        // Pick a vertex on this face and check that the centroid is on the
        // negative side of the plane (n.centroid + d < 0)
        int vi;
        switch(f) {
            case 0: vi = 0; break;
            case 1: vi = 7; break;
            case 2: vi = 0; break;
            case 3: vi = 3; break;
            case 4: vi = 0; break;
            case 5: vi = 1; break;
            default: vi = 0; break;
        }
        double dot = planes_[f][0] * cx + planes_[f][1] * cy + planes_[f][2] * cz + planes_[f][3];
        if(dot > 0.0) {
            // Normal points inward -- flip it
            planes_[f][0] = -planes_[f][0];
            planes_[f][1] = -planes_[f][1];
            planes_[f][2] = -planes_[f][2];
            planes_[f][3] = -planes_[f][3];
        }
    }
}

// ------------------------------------------------------------------------- //
void Trap::swap(Geometry& geometry) {
    Trap* trap = dynamic_cast<Trap*>(&geometry);
    if(!trap) {
        //log_warn("Cannot swap Trap!");
        return;
    }

    Geometry::swap(*trap);

    std::swap(dz_, trap->dz_);
    std::swap(theta_, trap->theta_);
    std::swap(phi_, trap->phi_);
    std::swap(dy1_, trap->dy1_);
    std::swap(dx1_, trap->dx1_);
    std::swap(dx2_, trap->dx2_);
    std::swap(alpha1_, trap->alpha1_);
    std::swap(dy2_, trap->dy2_);
    std::swap(dx3_, trap->dx3_);
    std::swap(dx4_, trap->dx4_);
    std::swap(alpha2_, trap->alpha2_);

    for(int i = 0; i < 8; ++i) {
        for(int j = 0; j < 3; ++j) {
            std::swap(vertices_[i][j], trap->vertices_[i][j]);
        }
    }
    for(int i = 0; i < 6; ++i) {
        for(int j = 0; j < 4; ++j) {
            std::swap(planes_[i][j], trap->planes_[i][j]);
        }
    }
}

//------------------------------------------------------------------------- //
Trap& Trap::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Trap* trap = dynamic_cast<const Trap*>(&geometry);
        if(!trap) {
            //log_warn("Cannot assign Trap!");
            return *this;
        }

        Trap tmp(*trap);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Trap::equal(const Geometry& geometry) const
{
    const Trap* trap = dynamic_cast<const Trap*>(&geometry);

    if(!trap)
        return false;
    else if(dz_ != trap->dz_)
        return false;
    else if(theta_ != trap->theta_)
        return false;
    else if(phi_ != trap->phi_)
        return false;
    else if(dy1_ != trap->dy1_)
        return false;
    else if(dx1_ != trap->dx1_)
        return false;
    else if(dx2_ != trap->dx2_)
        return false;
    else if(alpha1_ != trap->alpha1_)
        return false;
    else if(dy2_ != trap->dy2_)
        return false;
    else if(dx3_ != trap->dx3_)
        return false;
    else if(dx4_ != trap->dx4_)
        return false;
    else if(alpha2_ != trap->alpha2_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Trap::less(const Geometry& geometry) const
{
    const Trap* trap = dynamic_cast<const Trap*>(&geometry);
    if(!trap) return false;

    return
        std::tie(dz_, theta_, phi_, dy1_, dx1_, dx2_, alpha1_, dy2_, dx3_, dx4_, alpha2_)
        <
        std::tie(trap->dz_, trap->theta_, trap->phi_, trap->dy1_, trap->dx1_, trap->dx2_, trap->alpha1_, trap->dy2_, trap->dx3_, trap->dx4_, trap->alpha2_);
}

// ------------------------------------------------------------------------- //
void Trap::print(std::ostream& os) const
{
    os << "dz: " << dz_ << "\ttheta: " << theta_ << "\tphi: " << phi_
       << "\tdy1: " << dy1_ << "\tdx1: " << dx1_ << "\tdx2: " << dx2_
       << "\talpha1: " << alpha1_
       << "\tdy2: " << dy2_ << "\tdx3: " << dx3_ << "\tdx4: " << dx4_
       << "\talpha2: " << alpha2_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Trap::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Generalized slab method for a convex polyhedron with 6 planar faces.
    // For each face plane n.x + d = 0, compute the ray parameter t where
    // the ray hits that plane. Track the latest entry and earliest exit.

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dirx = direction.GetX();
    double diry = direction.GetY();
    double dirz = direction.GetZ();

    double t_enter = -1e30;
    double t_exit = 1e30;
    int enter_face = -1;
    int exit_face = -1;

    for(int f = 0; f < 6; ++f) {
        double ndotd = planes_[f][0] * dirx + planes_[f][1] * diry + planes_[f][2] * dirz;
        double ndoto = planes_[f][0] * px + planes_[f][1] * py + planes_[f][2] * pz + planes_[f][3];

        if(std::fabs(ndotd) < GEOMETRY_PRECISION) {
            // Ray parallel to this face
            if(ndoto > 0) return {}; // outside and parallel = miss
            continue;
        }

        double t = -ndoto / ndotd;

        if(ndotd < 0) {
            // Entering through this face
            if(t > t_enter) {
                t_enter = t;
                enter_face = f;
            }
        } else {
            // Exiting through this face
            if(t < t_exit) {
                t_exit = t;
                exit_face = f;
            }
        }
    }

    if(t_enter >= t_exit) return {}; // miss

    // Two intersections
    Intersection hits[2];
    siren::math::Vector3D enter_pos(px + t_enter * dirx, py + t_enter * diry, pz + t_enter * dirz);
    siren::math::Vector3D exit_pos(px + t_exit * dirx, py + t_exit * diry, pz + t_exit * dirz);
    hits[0].distance = t_enter;
    hits[0].hierarchy = 0;
    hits[0].entering = true;
    hits[0].position = enter_pos;
    hits[1].distance = t_exit;
    hits[1].hierarchy = 0;
    hits[1].entering = false;
    hits[1].position = exit_pos;
    return {hits[0], hits[1]};
}

// ------------------------------------------------------------------------- //
AABB Trap::GetBoundingBox() const {
    AABB box;
    for(int i = 0; i < 8; ++i) {
        box.ExpandToInclude(math::Vector3D(vertices_[i][0], vertices_[i][1], vertices_[i][2]));
    }
    return box;
}

} // namespace geometry
} // namespace siren
