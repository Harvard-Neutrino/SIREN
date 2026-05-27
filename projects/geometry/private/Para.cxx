#include "SIREN/geometry/Para.h"

#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <limits>
#include <ostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"
#include "GeometryMacros.h"

namespace siren {
namespace geometry {

void Para::ComputeFaceData() {
    // Edge vectors of the parallelepiped
    // e1 = (dx, 0, 0)
    // e2 = (dy*tan(alpha), dy, 0)
    // e3 = (dz*sin(theta)*cos(phi), dz*sin(theta)*sin(phi), dz*cos(theta))
    double ta = std::tan(alpha_);
    double tt = std::tan(theta_);
    double sp = std::sin(phi_);
    double cp = std::cos(phi_);

    double e1x = dx_, e1y = 0.0, e1z = 0.0;
    double e2x = dy_ * ta, e2y = dy_, e2z = 0.0;
    double e3x = dz_ * tt * cp, e3y = dz_ * tt * sp, e3z = dz_;

    // Face normals via cross products of edge pairs, then normalize
    // n1 = normalize(cross(e2, e3)) -- normal to face pair perpendicular to e1
    double n1x = e2y * e3z - e2z * e3y;
    double n1y = e2z * e3x - e2x * e3z;
    double n1z = e2x * e3y - e2y * e3x;
    double mag1 = std::sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
    if(mag1 > 0) { n1x /= mag1; n1y /= mag1; n1z /= mag1; }

    // n2 = normalize(cross(e3, e1)) -- normal to face pair perpendicular to e2
    double n2x = e3y * e1z - e3z * e1y;
    double n2y = e3z * e1x - e3x * e1z;
    double n2z = e3x * e1y - e3y * e1x;
    double mag2 = std::sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
    if(mag2 > 0) { n2x /= mag2; n2y /= mag2; n2z /= mag2; }

    // n3 = normalize(cross(e1, e2)) -- normal to face pair perpendicular to e3
    double n3x = e1y * e2z - e1z * e2y;
    double n3y = e1z * e2x - e1x * e2z;
    double n3z = e1x * e2y - e1y * e2x;
    double mag3 = std::sqrt(n3x * n3x + n3y * n3y + n3z * n3z);
    if(mag3 > 0) { n3x /= mag3; n3y /= mag3; n3z /= mag3; }

    face_normals_[0][0] = n1x; face_normals_[0][1] = n1y; face_normals_[0][2] = n1z;
    face_normals_[1][0] = n2x; face_normals_[1][1] = n2y; face_normals_[1][2] = n2z;
    face_normals_[2][0] = n3x; face_normals_[2][1] = n3y; face_normals_[2][2] = n3z;

    // Half-slab distances: d_i = |n_i . e_i|
    face_dists_[0] = std::fabs(n1x * e1x + n1y * e1y + n1z * e1z);
    face_dists_[1] = std::fabs(n2x * e2x + n2y * e2y + n2z * e2z);
    face_dists_[2] = std::fabs(n3x * e3x + n3y * e3y + n3z * e3z);
}

Para::Para() : Geometry("Para"), dx_(0), dy_(0), dz_(0), alpha_(0), theta_(0), phi_(0) {
    face_normals_[0][0] = 0; face_normals_[0][1] = 0; face_normals_[0][2] = 0;
    face_normals_[1][0] = 0; face_normals_[1][1] = 0; face_normals_[1][2] = 0;
    face_normals_[2][0] = 0; face_normals_[2][1] = 0; face_normals_[2][2] = 0;
    face_dists_[0] = 0; face_dists_[1] = 0; face_dists_[2] = 0;
    RecomputeWorldAABB();
}
Para::Para(double dx, double dy, double dz, double alpha, double theta, double phi)
    : Geometry("Para"), dx_(dx), dy_(dy), dz_(dz), alpha_(alpha), theta_(theta), phi_(phi) {
    if(dx_ <= 0 || dy_ <= 0 || dz_ <= 0) throw std::invalid_argument("Para half-lengths must be positive!");
    ComputeFaceData();
    RecomputeWorldAABB();
}
Para::Para(Placement const & p) : Geometry("Para", p), dx_(0), dy_(0), dz_(0), alpha_(0), theta_(0), phi_(0) {
    face_normals_[0][0] = 0; face_normals_[0][1] = 0; face_normals_[0][2] = 0;
    face_normals_[1][0] = 0; face_normals_[1][1] = 0; face_normals_[1][2] = 0;
    face_normals_[2][0] = 0; face_normals_[2][1] = 0; face_normals_[2][2] = 0;
    face_dists_[0] = 0; face_dists_[1] = 0; face_dists_[2] = 0;
    RecomputeWorldAABB();
}
Para::Para(Placement const & p, double dx, double dy, double dz, double alpha, double theta, double phi)
    : Geometry("Para", p), dx_(dx), dy_(dy), dz_(dz), alpha_(alpha), theta_(theta), phi_(phi) {
    if(dx_ <= 0 || dy_ <= 0 || dz_ <= 0) throw std::invalid_argument("Para half-lengths must be positive!");
    ComputeFaceData();
    RecomputeWorldAABB();
}
Para::Para(const Para& o) : Geometry(o), dx_(o.dx_), dy_(o.dy_), dz_(o.dz_), alpha_(o.alpha_), theta_(o.theta_), phi_(o.phi_) {
    for(int i = 0; i < 3; ++i) {
        face_normals_[i][0] = o.face_normals_[i][0];
        face_normals_[i][1] = o.face_normals_[i][1];
        face_normals_[i][2] = o.face_normals_[i][2];
        face_dists_[i] = o.face_dists_[i];
    }
    RecomputeWorldAABB();
}

void Para::swap(Geometry& geometry) {
    Para* o_ = dynamic_cast<Para*>(&geometry);
    if(!o_) return;
    Geometry::swap(*o_);
    std::swap(dx_, o_->dx_); std::swap(dy_, o_->dy_); std::swap(dz_, o_->dz_);
    std::swap(alpha_, o_->alpha_); std::swap(theta_, o_->theta_); std::swap(phi_, o_->phi_);
    for(int i = 0; i < 3; ++i) {
        std::swap(face_normals_[i][0], o_->face_normals_[i][0]);
        std::swap(face_normals_[i][1], o_->face_normals_[i][1]);
        std::swap(face_normals_[i][2], o_->face_normals_[i][2]);
        std::swap(face_dists_[i], o_->face_dists_[i]);
    }
}

SIREN_GEOMETRY_ASSIGN(Para)
SIREN_GEOMETRY_EQUAL(Para, dx_, dy_, dz_, alpha_, theta_, phi_)
SIREN_GEOMETRY_LESS(Para, dx_, dy_, dz_, alpha_, theta_, phi_)

void Para::print(std::ostream& os) const {
    os << "Para(" << dx_ << ", " << dy_ << ", " << dz_ << ", " << alpha_ << ", " << theta_ << ", " << phi_ << ")\n";
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Para::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Slab intersection method using 3 precomputed face normal pairs.
    // Same algorithm as Box but with oblique slabs defined by the
    // parallelepiped's face normals and distances.

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dirx = direction.GetX();
    double diry = direction.GetY();
    double dirz = direction.GetZ();

    double t_enter = -std::numeric_limits<double>::infinity();
    double t_exit  =  std::numeric_limits<double>::infinity();

    for(int i = 0; i < 3; ++i) {
        double nx = face_normals_[i][0];
        double ny = face_normals_[i][1];
        double nz = face_normals_[i][2];
        double fd = face_dists_[i];

        double ndotd = nx * dirx + ny * diry + nz * dirz;
        double ndoto = nx * px + ny * py + nz * pz;

        if(std::fabs(ndotd) < GEOMETRY_PRECISION) {
            // Ray is parallel to this slab pair
            if(std::fabs(ndoto) > fd) return {};
            continue;
        }

        double t1 = (-fd - ndoto) / ndotd;
        double t2 = ( fd - ndoto) / ndotd;

        if(ndotd < 0) { std::swap(t1, t2); }

        if(t1 > t_enter) t_enter = t1;
        if(t2 < t_exit)  t_exit  = t2;

        if(t_enter > t_exit) return {};
    }

    // Precision control
    if(t_enter > 0 && t_enter < GEOMETRY_PRECISION) t_enter = 0;
    if(t_exit > 0 && t_exit < GEOMETRY_PRECISION) t_exit = 0;

    Intersection hit_enter, hit_exit;
    hit_enter.distance = t_enter;
    hit_enter.hierarchy = 0;
    hit_enter.entering = true;
    hit_enter.position = siren::math::Vector3D(px + t_enter * dirx, py + t_enter * diry, pz + t_enter * dirz);

    hit_exit.distance = t_exit;
    hit_exit.hierarchy = 0;
    hit_exit.entering = false;
    hit_exit.position = siren::math::Vector3D(px + t_exit * dirx, py + t_exit * diry, pz + t_exit * dirz);

    return {hit_enter, hit_exit};
}

// ------------------------------------------------------------------------- //
AABB Para::GetBoundingBox() const {
    // Compute all 8 vertices and take the axis-aligned envelope
    double ta = std::tan(alpha_);
    double tt = std::tan(theta_);
    double sp = std::sin(phi_);
    double cp = std::cos(phi_);

    double e1x = dx_, e1y = 0.0, e1z = 0.0;
    double e2x = dy_ * ta, e2y = dy_, e2z = 0.0;
    double e3x = dz_ * tt * cp, e3y = dz_ * tt * sp, e3z = dz_;

    AABB box;
    for(int s1 = -1; s1 <= 1; s1 += 2) {
        for(int s2 = -1; s2 <= 1; s2 += 2) {
            for(int s3 = -1; s3 <= 1; s3 += 2) {
                double vx = s1 * e1x + s2 * e2x + s3 * e3x;
                double vy = s1 * e1y + s2 * e2y + s3 * e3y;
                double vz = s1 * e1z + s2 * e2z + s3 * e3z;
                box.ExpandToInclude(math::Vector3D(vx, vy, vz));
            }
        }
    }
    return box;
}

} // namespace geometry
} // namespace siren
