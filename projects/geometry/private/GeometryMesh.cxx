#include "SIREN/geometry/GeometryMesh.h"

#include <set>
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
}

TriangularMesh::TriangularMesh(std::vector<std::array<math::Vector3D, 3>> const & triangles)
    : Geometry("TriangularMesh")
    , mesh(BuildMeshFromTriangles(triangles))
{
    BuildAccelerationStructure();
}

TriangularMesh::TriangularMesh(Placement const & placement, std::vector<std::array<math::Vector3D, 3>> const & triangles)
    : Geometry("TriangularMesh", placement)
    , mesh(BuildMeshFromTriangles(triangles))
{
    BuildAccelerationStructure();
}

Mesh::TMesh TriangularMesh::BuildMeshFromTriangles(std::vector<std::array<math::Vector3D, 3>> const & triangles) {
    Mesh::TMesh m;
    int vertex_id = 0;
    for(auto const & tri : triangles) {
        int v0 = vertex_id++;
        int v1 = vertex_id++;
        int v2 = vertex_id++;

        Mesh::VData d0 = {{tri[0].GetX(), tri[0].GetY(), tri[0].GetZ()}};
        Mesh::VData d1 = {{tri[1].GetX(), tri[1].GetY(), tri[1].GetZ()}};
        Mesh::VData d2 = {{tri[2].GetX(), tri[2].GetY(), tri[2].GetZ()}};

        Mesh::VAttribute va0; va0.data = d0;
        Mesh::VAttribute va1; va1.data = d1;
        Mesh::VAttribute va2; va2.data = d2;

        Mesh::Triangle t = {{v0, v1, v2}};
        Mesh::TData tdata = {{d0, d1, d2}};
        Mesh::TAttribute ta; ta.data = tdata;

        va0.tset.insert(t);
        va1.tset.insert(t);
        va2.tset.insert(t);

        Mesh::Edge e01 = {{std::min(v0, v1), std::max(v0, v1)}};
        Mesh::Edge e12 = {{std::min(v1, v2), std::max(v1, v2)}};
        Mesh::Edge e20 = {{std::min(v2, v0), std::max(v2, v0)}};

        va0.eset.insert(e01); va0.eset.insert(e20);
        va1.eset.insert(e01); va1.eset.insert(e12);
        va2.eset.insert(e12); va2.eset.insert(e20);

        m.vmap.push_back(va0);
        m.vmap.push_back(va1);
        m.vmap.push_back(va2);

        m.tmap[t] = ta;

        Mesh::EData ed01 = {{d0, d1}};
        Mesh::EData ed12 = {{d1, d2}};
        Mesh::EData ed20 = {{d2, d0}};

        if(m.emap.find(e01) == m.emap.end()) {
            Mesh::EAttribute ea; ea.data = ed01;
            ea.tset.insert(t);
            m.emap[e01] = ea;
        } else {
            m.emap[e01].tset.insert(t);
        }
        if(m.emap.find(e12) == m.emap.end()) {
            Mesh::EAttribute ea; ea.data = ed12;
            ea.tset.insert(t);
            m.emap[e12] = ea;
        } else {
            m.emap[e12].tset.insert(t);
        }
        if(m.emap.find(e20) == m.emap.end()) {
            Mesh::EAttribute ea; ea.data = ed20;
            ea.tset.insert(t);
            m.emap[e20] = ea;
        } else {
            m.emap[e20].tset.insert(t);
        }
    }
    return m;
}

void TriangularMesh::RebuildFromTriangleData() {
    std::vector<std::array<math::Vector3D, 3>> triangles;
    triangles.reserve(triangle_data_.size());
    for(auto const & td : triangle_data_) {
        std::array<math::Vector3D, 3> tri = {{
            math::Vector3D(td[0][0], td[0][1], td[0][2]),
            math::Vector3D(td[1][0], td[1][1], td[1][2]),
            math::Vector3D(td[2][0], td[2][1], td[2][2])
        }};
        triangles.push_back(tri);
    }
    mesh = BuildMeshFromTriangles(triangles);
    kd_root_ = Mesh::BuildKDTree(triangle_data_, 1.0, 1.5, 20);
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
    if(!tmesh) return false;

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

    // Build KD-tree. For closed-surface meshes (tessellated solids, arb8),
    // triangles overlap heavily in bounding-box terms and the SAH builder
    // degenerates (duplicating all triangles into both children). Use a flat
    // leaf (brute-force) for meshes up to 1024 triangles. For larger meshes,
    // cap depth at log2(n) to bound degenerate cases.
    double traversal_cost = 1.0;
    double intersection_cost = 1.5;
    int n = (int)triangle_data_.size();
    int max_depth;
    if(n <= 1024) {
        max_depth = 0;
    } else {
        max_depth = 0;
        int tmp = n;
        while(tmp > 1) { ++max_depth; tmp >>= 1; }
        if(max_depth > 16) max_depth = 16;
    }
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

    // Deduplicate hits by triangle ID.
    // A triangle can appear in multiple KD-tree leaves and produce
    // duplicate intersections. Build an index sorted by distance,
    // then keep only the first occurrence of each triangle ID.
    std::vector<size_t> order(hits.size());
    for(size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(),
        [&hits](size_t a, size_t b) {
            return hits[a].t < hits[b].t;
        });

    std::set<Mesh::TriangleID> seen_tris;
    std::vector<Intersection> unique_tri_hits;
    unique_tri_hits.reserve(hits.size());
    for(size_t idx : order) {
        if(!seen_tris.insert(hit_tri_ids[idx]).second) {
            continue;
        }
        Intersection isect;
        isect.distance = hits[idx].t;
        isect.hierarchy = 0;
        isect.entering = hits[idx].front_face;
        isect.matID = 0;
        isect.position = position + direction * hits[idx].t;
        unique_tri_hits.push_back(isect);
    }

    // Merge coincident hits from adjacent triangles sharing an edge/vertex.
    // When a ray hits a shared edge, both triangles report a hit at the same
    // distance. Keep only one hit per coincident group.
    result.reserve(unique_tri_hits.size());
    for(auto const & isect : unique_tri_hits) {
        if(!result.empty()) {
            double scale = std::fmax(std::fabs(isect.distance), std::fabs(result.back().distance));
            double tol = std::fmax(1e-9, scale * 1e-9);
            if(std::fabs(isect.distance - result.back().distance) < tol) {
                continue;
            }
        }
        result.push_back(isect);
    }

    return result;
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
