#include "SIREN/geometry/GeometryMesh.h"

#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/MeshBuilder.h"

using namespace siren::geometry;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace siren {
namespace geometry {

TriangularMesh::TriangularMesh()
    : Geometry((std::string)("TriangularMesh"))
{
    // Do nothing here
}

TriangularMesh::TriangularMesh(Mesh::TMesh const & mesh)
    : Geometry("TriangularMesh")
    , mesh(mesh)
{
    BuildAccelerationStructure();
}

TriangularMesh::TriangularMesh(Placement const & placement)
    : Geometry((std::string)("TriangularMesh"), placement)
{
    // Do nothing here
}

TriangularMesh::TriangularMesh(Placement const & placement, Mesh::TMesh const & mesh)
    : Geometry((std::string)("TriangularMesh"), placement)
    , mesh(mesh)
{
    BuildAccelerationStructure();
}

TriangularMesh::TriangularMesh(const TriangularMesh& tmesh)
    : Geometry(tmesh)
    , mesh(tmesh.mesh)
    , triangle_data_(tmesh.triangle_data_)
    , kd_root_(tmesh.kd_root_)
{
    // Share the KD-tree via shared_ptr (immutable after construction)
}

Mesh::VAttribute & TriangularMesh::GetVertex(Mesh::Vertex v) {
    return mesh.vmap[v];
}

Mesh::EAttribute & TriangularMesh::GetEdge(Mesh::Edge e) {
    return mesh.emap[e];
}

Mesh::TAttribute & TriangularMesh::GetTriangle(Mesh::Triangle t) {
    return mesh.tmap[t];
}

Mesh::VAttribute const & TriangularMesh::GetVertex(Mesh::Vertex v) const {
    return mesh.vmap[v];
}

Mesh::EAttribute const & TriangularMesh::GetEdge(Mesh::Edge e) const {
    return mesh.emap.at(e);
}

Mesh::TAttribute const & TriangularMesh::GetTriangle(Mesh::Triangle t) const {
    return mesh.tmap.at(t);
}

// ------------------------------------------------------------------------- //
void TriangularMesh::swap(Geometry& geometry)
{
    TriangularMesh* tmesh = dynamic_cast<TriangularMesh*>(&geometry);
    if (!tmesh)
    {
        //log_warn("Cannot swap TriangularMesh!");
        return;
    }

    Geometry::swap(*tmesh);

    std::swap(mesh, tmesh->mesh);
    std::swap(triangle_data_, tmesh->triangle_data_);
    std::swap(kd_root_, tmesh->kd_root_);
}

//------------------------------------------------------------------------- //
TriangularMesh& TriangularMesh::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const TriangularMesh* tmesh = dynamic_cast<const TriangularMesh*>(&geometry);
        if (!tmesh)
        {
            //log_warn("Cannot assign Sphere!");
            return *this;
        }

        TriangularMesh tmp(*tmesh);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool TriangularMesh::equal(const Geometry& geometry) const
{
    const TriangularMesh* tmesh = dynamic_cast<const TriangularMesh*>(&geometry);

    if (!tmesh)
        return false;
    else if (mesh != tmesh->mesh)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool TriangularMesh::less(const Geometry& geometry) const
{
    const TriangularMesh* tmesh = dynamic_cast<const TriangularMesh*>(&geometry);

    return mesh < tmesh->mesh;
}

// ------------------------------------------------------------------------- //
void TriangularMesh::print(std::ostream& os) const
{
}

// ------------------------------------------------------------------------- //
void TriangularMesh::BuildAccelerationStructure() {
    triangle_data_.clear();
    kd_root_.reset();

    if(mesh.tmap.empty()) return;

    // Extract triangle vertex data into a flat vector for the KD-tree builder
    triangle_data_.reserve(mesh.tmap.size());
    for(auto const & tp : mesh.tmap) {
        triangle_data_.push_back(tp.second.data);
    }

    // Build KD-tree with default parameters
    double traversal_cost = 1.0;
    double intersection_cost = 1.5;
    int max_depth = 20;
    kd_root_ = Mesh::BuildKDTree(triangle_data_, traversal_cost, intersection_cost, max_depth);
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> TriangularMesh::ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const {
    std::vector<Geometry::Intersection> result;

    if(!kd_root_ || triangle_data_.empty()) return result;

    Mesh::Point origin = {{position.GetX(), position.GetY(), position.GetZ()}};
    Mesh::Point dir = {{direction.GetX(), direction.GetY(), direction.GetZ()}};
    Mesh::Point inv_dir = {{1.0 / dir[0], 1.0 / dir[1], 1.0 / dir[2]}};

    std::vector<Mesh::RayTriangleHit> hits;
    std::vector<Mesh::TriangleID> hit_tri_ids;
    Mesh::TraverseKDTree(kd_root_, triangle_data_, origin, dir, inv_dir, hits, hit_tri_ids);

    result.reserve(hits.size());
    for(size_t i = 0; i < hits.size(); ++i) {
        Intersection isect;
        isect.distance = hits[i].t;
        isect.hierarchy = 0;
        isect.entering = hits[i].front_face;
        isect.matID = 0;
        isect.position = position + direction * hits[i].t;
        result.push_back(isect);
    }

    // Sort by distance
    std::sort(result.begin(), result.end(),
        [](Intersection const & a, Intersection const & b) {
            return a.distance < b.distance;
        });

    return result;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> TriangularMesh::ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const {
    std::vector<Intersection> intersections = ComputeIntersections(position, direction);

    if(intersections.empty()) {
        return std::make_pair(-1.0, -1.0);
    }

    // Find the first positive-distance intersection
    double dist_1 = -1.0;
    double dist_2 = -1.0;

    for(auto const & isect : intersections) {
        if(isect.distance > GEOMETRY_PRECISION) {
            if(dist_1 < 0) {
                dist_1 = isect.distance;
            } else {
                dist_2 = isect.distance;
                break;
            }
        }
    }

    return std::make_pair(dist_1, dist_2);
}

// ------------------------------------------------------------------------- //
AABB TriangularMesh::GetBoundingBox() const {
    AABB box;
    for(auto const & vattr : mesh.vmap) {
        box.ExpandToInclude(math::Vector3D(vattr.data[0], vattr.data[1], vattr.data[2]));
    }
    return box;
}

} // namespace geometry
} // namespace siren
