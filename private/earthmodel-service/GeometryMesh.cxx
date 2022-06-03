#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/Placement.h"
#include "earthmodel-service/GeometryMesh.h"

using namespace earthmodel;

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace earthmodel {

bool TriangularMesh::TMesh::operator<(TriangularMesh::TMesh other) {
    return false;
}

TriangularMesh::TriangularMesh()
    : Geometry((std::string)("TriangularMesh"))
{
    // Do nothing here
}

TriangularMesh::TriangularMesh(TriangularMesh::TMesh const & mesh)
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

TriangularMesh::TriangularMesh(Placement const & placement, TriangularMesh const & mesh)
    : Geometry((std::string)("TriangularMesh"), placement)
    , mesh(mesh)
{
    // Do nothing here
}

TriangularMesh::TriangularMesh(const TriangularMesh& tmesh)
    : Geometry(box)
    , mesh(tmesh.mesh)
{
    // Nothing to do here
}

TriangularMesh::VAttribute & TriangularMesh::GetVertex(TriangularMesh::Vertex v) {
    return mesh.vmap[v];
}

TriangularMesh::EAttribute & TriangularMesh::GetEdge(TriangularMesh::Edge e) {
    return mesh.emap[e];
}

TriangularMesh::TAttribute & TriangularMesh::GetTriangle(TriangularMesh::Triangle t) {
    return mesh.tmap[t];
}

TriangularMesh::VAttribute const & TriangularMesh::GetVertex(TriangularMesh::Vertex v) const {
    return mesh.vmap[v];
}

TriangularMesh::EAttribute const & TriangularMesh::GetEdge(TriangularMesh::Edge e) const {
    return mesh.emap[e];
}

TriangularMesh::TAttribute const & TriangularMesh::GetTriangle(TriangularMesh::Triangle t) const {
    return mesh.tmap[t];
}

// ------------------------------------------------------------------------- //
void TriangularMesh::swap(Geometry& geometry)
{
    TriangularMesh* tmesh = dynamic_cast<TriangularMesh*>(&geometry);
    if (!mesh)
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
        if (!box)
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
