#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/Placement.h"
#include "LeptonInjector/geometry/MeshBuilder.h"
#include "LeptonInjector/geometry/GeometryMesh.h"

using namespace LI::geometry;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace LI {
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
    // Do nothing here
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
    // Do nothing here
}

TriangularMesh::TriangularMesh(const TriangularMesh& tmesh)
    : Geometry()
    , mesh(tmesh.mesh)
{
    // Nothing to do here
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
std::vector<Geometry::Intersection> TriangularMesh::ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const {
    return std::vector<Geometry::Intersection>();
}

// ------------------------------------------------------------------------- //
std::pair<double, double> TriangularMesh::ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const
{
    return std::pair<double, double>();
}

} // namespace geometry
} // namespace LI
