#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/Placement.h"
#include "earthmodel-service/MeshBuilder.h"
#include "earthmodel-service/GeometryMesh.h"

using namespace earthmodel;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace earthmodel {

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
std::vector<Geometry::Intersection> TriangularMesh::ComputeIntersections(Vector3D const & position, Vector3D const & direction) const {
}

// ------------------------------------------------------------------------- //
std::pair<double, double> TriangularMesh::ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const
{
}

} // namespace earthmodel
