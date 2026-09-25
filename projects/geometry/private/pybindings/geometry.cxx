
#include <array>
#include <vector>

#include "../../public/SIREN/geometry/Placement.h"
#include "../../public/SIREN/geometry/Geometry.h"
#include "../../public/SIREN/geometry/ExtrPoly.h"
#include "../../public/SIREN/geometry/Cylinder.h"
#include "../../public/SIREN/geometry/Box.h"
#include "../../public/SIREN/geometry/Cone.h"
#include "../../public/SIREN/geometry/Polycone.h"
#include "../../public/SIREN/geometry/Polyhedra.h"
#include "../../public/SIREN/geometry/GenericPolycone.h"
#include "../../public/SIREN/geometry/Trd.h"
#include "../../public/SIREN/geometry/Torus.h"
#include "../../public/SIREN/geometry/Sphere.h"
#include "../../public/SIREN/geometry/BooleanGeometry.h"
#include "../../public/SIREN/geometry/AABB.h"
#include "../../public/SIREN/geometry/EllipticalTube.h"
#include "../../public/SIREN/geometry/CutTube.h"
#include "../../public/SIREN/geometry/Trap.h"
#include "../../public/SIREN/geometry/Ellipsoid.h"
#include "../../public/SIREN/geometry/Para.h"
#include "../../public/SIREN/geometry/GeometryMesh.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>



using namespace pybind11;

PYBIND11_MODULE(geometry,m) {
    using namespace siren::geometry;

    // geometry

    class_<Geometry, std::shared_ptr<Geometry>>(m, "Geometry")
        .def("IsInside", (bool (Geometry::*)(siren::math::Vector3D const &, siren::math::Vector3D const &) const) &Geometry::IsInside)
        .def("IsInside", (bool (Geometry::*)(siren::math::Vector3D const &) const) &Geometry::IsInside)
        .def("Intersections", (std::vector<Geometry::Intersection> (Geometry::*)(siren::math::Vector3D const &, siren::math::Vector3D const &) const) &Geometry::Intersections)
        .def("DistanceToClosestApproach",&Geometry::DistanceToClosestApproach)
        .def("GetBoundingBox",&Geometry::GetBoundingBox)
        .def("GetWorldBoundingBox",&Geometry::GetWorldBoundingBox)
        .def_property_readonly("name",&Geometry::GetName)
        .def_property("placement",&Geometry::GetPlacement, &Geometry::SetPlacement)
        .def("ComputeIntersections",&Geometry::ComputeIntersections)
        .def("create",&Geometry::create);

    // data structs in Geometry class

    class_<Geometry::Intersection>(m, "Intersection")
        .def_readwrite("distance", &Geometry::Intersection::distance)
        .def_readwrite("hierarchy", &Geometry::Intersection::hierarchy)
        .def_readwrite("entering", &Geometry::Intersection::entering)
        .def_readwrite("matID", &Geometry::Intersection::matID)
        .def_readwrite("position", &Geometry::Intersection::position);

    class_<Geometry::IntersectionList>(m, "IntersectionList")
        .def_readwrite("position", &Geometry::IntersectionList::position)
        .def_readwrite("direction", &Geometry::IntersectionList::direction)
        .def_readwrite("intersections", &Geometry::IntersectionList::intersections);

    // ExtrPoly

    class_<ExtrPoly::ZSection>(m, "ZSection")
        .def(init<>())
        .def(init<double, double*, double>())
        .def_readwrite("zpos", &ExtrPoly::ZSection::zpos)
        .def_readwrite("scale", &ExtrPoly::ZSection::scale)
        .def_readwrite("offset", &ExtrPoly::ZSection::offset);

    class_<ExtrPoly::plane>(m, "plane")
        .def_readwrite("a", &ExtrPoly::plane::a)
        .def_readwrite("b", &ExtrPoly::plane::b)
        .def_readwrite("c", &ExtrPoly::plane::c)
        .def_readwrite("d", &ExtrPoly::plane::d);

    class_<ExtrPoly, std::shared_ptr<ExtrPoly>, Geometry>(m, "ExtrPoly")
        .def(init<>())
        .def(init<const std::vector<std::vector<double>>&,
                const std::vector<ExtrPoly::ZSection>&>())
        .def(init<Placement const &,
                const std::vector<std::vector<double>>&,
                const std::vector<ExtrPoly::ZSection>&>())
        .def(init<Placement const &>())
        .def(init<const ExtrPoly&>())
        .def_property("polygon",&ExtrPoly::GetPolygon, &ExtrPoly::SetPolygon)
        .def_property("zsections",&ExtrPoly::GetZSections, &ExtrPoly::SetZSections)
        .def("ComputeLateralPlanes",&ExtrPoly::ComputeLateralPlanes);

    // Box

    class_<Box, std::shared_ptr<Box>, Geometry>(m, "Box")
        .def(init<>())
        // Positional Box(x, y, z) places the box at the ORIGIN, which reads as
        // a size-only constructor and silently drops the intended center. Warn
        // and keep the origin placement rather than raising, so positional
        // callers are not broken immediately.
        .def(init([](double x, double y, double z) {
            PyErr_WarnEx(PyExc_DeprecationWarning,
                "Box(x, y, z) places the box at the ORIGIN; pass "
                "Box(widths=(x, y, z), center=(cx, cy, cz))",
                1);
            return std::make_shared<Box>(x, y, z);
        }))
        .def(init([](std::array<double, 3> widths, std::array<double, 3> center) {
            Placement placement(siren::math::Vector3D(center[0], center[1], center[2]));
            return std::make_shared<Box>(placement, widths[0], widths[1], widths[2]);
        }), arg("widths"), arg("center") = std::array<double, 3>{0.0, 0.0, 0.0})
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<const Box&>())
        .def_property("X",&Box::GetX, &Box::SetX)
        .def_property("Y",&Box::GetY, &Box::SetY)
        .def_property("Z",&Box::GetZ, &Box::SetZ);

    // TriangularMesh

    class_<TriangularMesh, std::shared_ptr<TriangularMesh>, Geometry>(m, "TriangularMesh")
        .def(init<>())
        .def(init<std::vector<std::array<siren::math::Vector3D, 3>> const &>())
        .def(init<Placement const &, std::vector<std::array<siren::math::Vector3D, 3>> const &>())
        .def(init<Placement const &>())
        .def(init<const TriangularMesh&>())
        .def("TriangleCount", &TriangularMesh::TriangleCount)
        .def("GetTriangles", &TriangularMesh::GetTriangles)
        .def("ValidateClosed", &TriangularMesh::ValidateClosed);

    // Cone

    class_<Cone, std::shared_ptr<Cone>, Geometry>(m, "Cone")
        .def(init<>())
        .def(init<double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double, double, double>())
        .def(init<double, double, double, double, double, double, double>())
        .def(init<Placement const &, double, double, double, double, double, double, double>())
        .def(init<const Cone&>())
        .def_property_readonly("Rmin1",&Cone::GetRmin1)
        .def_property_readonly("Rmax1",&Cone::GetRmax1)
        .def_property_readonly("Rmin2",&Cone::GetRmin2)
        .def_property_readonly("Rmax2",&Cone::GetRmax2)
        .def_property_readonly("Z",&Cone::GetZ)
        .def_property_readonly("StartPhi",&Cone::GetStartPhi)
        .def_property_readonly("DeltaPhi",&Cone::GetDeltaPhi);

    // Torus

    class_<Torus, std::shared_ptr<Torus>, Geometry>(m, "Torus")
        .def(init<>())
        .def(init<double, double, double>())
        .def(init<double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<Placement const &, double, double, double, double, double>())
        .def(init<const Torus&>())
        .def_property_readonly("MajorRadius",&Torus::GetMajorRadius)
        .def_property_readonly("MinorRadius",&Torus::GetMinorRadius)
        .def_property_readonly("InnerRadius",&Torus::GetInnerRadius)
        .def_property_readonly("StartPhi",&Torus::GetStartPhi)
        .def_property_readonly("DeltaPhi",&Torus::GetDeltaPhi);

    // Polycone

    class_<Polycone, std::shared_ptr<Polycone>, Geometry>(m, "Polycone")
        .def(init<>())
        .def(init<std::vector<double>, std::vector<double>, std::vector<double>>())
        .def(init<Placement const &>())
        .def(init<Placement const &, std::vector<double>, std::vector<double>, std::vector<double>>())
        .def(init<std::vector<double>, std::vector<double>, std::vector<double>, double, double>())
        .def(init<Placement const &, std::vector<double>, std::vector<double>, std::vector<double>, double, double>())
        .def(init<const Polycone&>())
        .def_property_readonly("ZPlanes",&Polycone::GetZPlanes)
        .def_property_readonly("Rmin",&Polycone::GetRmin)
        .def_property_readonly("Rmax",&Polycone::GetRmax)
        .def_property_readonly("StartPhi",&Polycone::GetStartPhi)
        .def_property_readonly("DeltaPhi",&Polycone::GetDeltaPhi);

    // Polyhedra

    class_<Polyhedra, std::shared_ptr<Polyhedra>, Geometry>(m, "Polyhedra")
        .def(init<>())
        .def(init<int, double, std::vector<double>, std::vector<double>, std::vector<double>>())
        .def(init<Placement const &>())
        .def(init<Placement const &, int, double, std::vector<double>, std::vector<double>, std::vector<double>>())
        .def(init<int, double, std::vector<double>, std::vector<double>, std::vector<double>, double>())
        .def(init<Placement const &, int, double, std::vector<double>, std::vector<double>, std::vector<double>, double>())
        .def(init<const Polyhedra&>())
        .def_property_readonly("NumSides",&Polyhedra::GetNumSides)
        .def_property_readonly("StartPhi",&Polyhedra::GetStartPhi)
        .def_property_readonly("DeltaPhi",&Polyhedra::GetDeltaPhi)
        .def_property_readonly("ZPlanes",&Polyhedra::GetZPlanes)
        .def_property_readonly("Rmin",&Polyhedra::GetRmin)
        .def_property_readonly("Rmax",&Polyhedra::GetRmax);

    // GenericPolycone

    class_<GenericPolycone, std::shared_ptr<GenericPolycone>, Geometry>(m, "GenericPolycone")
        .def(init<>())
        .def(init<std::vector<double>, std::vector<double>>())
        .def(init<Placement const &, std::vector<double>, std::vector<double>>())
        .def(init<std::vector<double>, std::vector<double>, double, double>())
        .def(init<Placement const &, std::vector<double>, std::vector<double>, double, double>())
        .def(init<const GenericPolycone&>())
        .def_property_readonly("R",&GenericPolycone::GetR)
        .def_property_readonly("Z",&GenericPolycone::GetZ)
        .def_property_readonly("StartPhi",&GenericPolycone::GetStartPhi)
        .def_property_readonly("DeltaPhi",&GenericPolycone::GetDeltaPhi);

    // Trd

    class_<Trd, std::shared_ptr<Trd>, Geometry>(m, "Trd")
        .def(init<>())
        .def(init<double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double, double, double>())
        .def(init<const Trd&>())
        .def_property_readonly("Dx1",&Trd::GetDx1)
        .def_property_readonly("Dx2",&Trd::GetDx2)
        .def_property_readonly("Dy1",&Trd::GetDy1)
        .def_property_readonly("Dy2",&Trd::GetDy2)
        .def_property_readonly("Dz",&Trd::GetDz);

    // Cylinder

    class_<Cylinder, std::shared_ptr<Cylinder>, Geometry>(m, "Cylinder")
        .def(init<>())
        .def(init<double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<double, double, double, double, double>())
        .def(init<Placement const &, double, double, double, double, double>())
        .def(init<const Cylinder&>())
        .def_property("InnerRadius",&Cylinder::GetInnerRadius, &Cylinder::SetInnerRadius)
        .def_property("Radius",&Cylinder::GetRadius, &Cylinder::SetRadius)
        .def_property("Z",&Cylinder::GetZ, &Cylinder::SetZ)
        .def_property_readonly("StartPhi",&Cylinder::GetStartPhi)
        .def_property_readonly("DeltaPhi",&Cylinder::GetDeltaPhi);

    // Sphere

    class_<Sphere, std::shared_ptr<Sphere>, Geometry>(m, "Sphere")
        .def(init<>())
        .def(init<double, double>())
        .def(init<double, double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double>())
        .def(init<Placement const &, double, double, double, double, double, double>())
        .def(init<const Sphere&>())
        .def_property_readonly("InnerRadius",&Sphere::GetInnerRadius)
        .def_property_readonly("Radius",&Sphere::GetRadius)
        .def_property_readonly("StartPhi",&Sphere::GetStartPhi)
        .def_property_readonly("DeltaPhi",&Sphere::GetDeltaPhi)
        .def_property_readonly("StartTheta",&Sphere::GetStartTheta)
        .def_property_readonly("DeltaTheta",&Sphere::GetDeltaTheta);

    // BooleanGeometry

    enum_<BooleanOperation>(m, "BooleanOperation")
        .value("UNION", BooleanOperation::UNION)
        .value("SUBTRACTION", BooleanOperation::SUBTRACTION)
        .value("INTERSECTION", BooleanOperation::INTERSECTION);

    class_<BooleanGeometry, std::shared_ptr<BooleanGeometry>, Geometry>(m, "BooleanGeometry")
        .def(init<>())
        .def(init<BooleanOperation, std::shared_ptr<const Geometry>, std::shared_ptr<const Geometry>>())
        .def(init<Placement const &, BooleanOperation, std::shared_ptr<const Geometry>, std::shared_ptr<const Geometry>>())
        .def(init<const BooleanGeometry&>())
        .def_property_readonly("Operation",&BooleanGeometry::GetOperation)
        .def_property_readonly("Left",&BooleanGeometry::GetLeft)
        .def_property_readonly("Right",&BooleanGeometry::GetRight);

    // EllipticalTube

    class_<EllipticalTube, std::shared_ptr<EllipticalTube>, Geometry>(m, "EllipticalTube")
        .def(init<>())
        .def(init<double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<const EllipticalTube&>())
        .def_property_readonly("Dx",&EllipticalTube::GetDx)
        .def_property_readonly("Dy",&EllipticalTube::GetDy)
        .def_property_readonly("Dz",&EllipticalTube::GetDz);

    // CutTube

    class_<CutTube, std::shared_ptr<CutTube>, Geometry>(m, "CutTube")
        .def(init<>())
        .def(init<double, double, double, siren::math::Vector3D, siren::math::Vector3D>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double, siren::math::Vector3D, siren::math::Vector3D>())
        .def(init<double, double, double, siren::math::Vector3D, siren::math::Vector3D, double, double>())
        .def(init<Placement const &, double, double, double, siren::math::Vector3D, siren::math::Vector3D, double, double>())
        .def(init<const CutTube&>())
        .def_property_readonly("Rmin",&CutTube::GetRmin)
        .def_property_readonly("Rmax",&CutTube::GetRmax)
        .def_property_readonly("Dz",&CutTube::GetDz)
        .def_property_readonly("LowNorm",&CutTube::GetLowNorm)
        .def_property_readonly("HighNorm",&CutTube::GetHighNorm)
        .def_property_readonly("StartPhi",&CutTube::GetStartPhi)
        .def_property_readonly("DeltaPhi",&CutTube::GetDeltaPhi);

    // Trap

    class_<Trap, std::shared_ptr<Trap>, Geometry>(m, "Trap")
        .def(init<>())
        .def(init<double, double, double, double, double, double, double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double, double, double, double, double, double, double, double, double>())
        .def(init<const Trap&>())
        .def_property_readonly("Dz",&Trap::GetDz)
        .def_property_readonly("Theta",&Trap::GetTheta)
        .def_property_readonly("Phi",&Trap::GetPhi)
        .def_property_readonly("Dy1",&Trap::GetDy1)
        .def_property_readonly("Dx1",&Trap::GetDx1)
        .def_property_readonly("Dx2",&Trap::GetDx2)
        .def_property_readonly("Alpha1",&Trap::GetAlpha1)
        .def_property_readonly("Dy2",&Trap::GetDy2)
        .def_property_readonly("Dx3",&Trap::GetDx3)
        .def_property_readonly("Dx4",&Trap::GetDx4)
        .def_property_readonly("Alpha2",&Trap::GetAlpha2);

    // Ellipsoid

    class_<Ellipsoid, std::shared_ptr<Ellipsoid>, Geometry>(m, "Ellipsoid")
        .def(init<>())
        .def(init<double, double, double>())
        .def(init<double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<Placement const &, double, double, double, double, double>())
        .def(init<const Ellipsoid&>())
        .def_property_readonly("Ax",&Ellipsoid::GetAx)
        .def_property_readonly("By",&Ellipsoid::GetBy)
        .def_property_readonly("Cz",&Ellipsoid::GetCz)
        .def_property_readonly("Zcut1",&Ellipsoid::GetZcut1)
        .def_property_readonly("Zcut2",&Ellipsoid::GetZcut2);

    // Para

    class_<Para, std::shared_ptr<Para>, Geometry>(m, "Para")
        .def(init<>())
        .def(init<double, double, double, double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double, double, double, double>())
        .def(init<const Para&>())
        .def_property_readonly("Dx",&Para::GetDx)
        .def_property_readonly("Dy",&Para::GetDy)
        .def_property_readonly("Dz",&Para::GetDz)
        .def_property_readonly("Alpha",&Para::GetAlpha)
        .def_property_readonly("Theta",&Para::GetTheta)
        .def_property_readonly("Phi",&Para::GetPhi);

    // AABB

    class_<AABB>(m, "AABB")
        .def(init<>())
        .def(init<siren::math::Vector3D const &, siren::math::Vector3D const &>())
        .def_readwrite("min_corner", &AABB::min_corner)
        .def_readwrite("max_corner", &AABB::max_corner)
        .def("ExpandToInclude", (void (AABB::*)(siren::math::Vector3D const &)) &AABB::ExpandToInclude)
        .def("Centroid", &AABB::Centroid)
        .def("SurfaceArea", &AABB::SurfaceArea)
        .def("LargestAxis", &AABB::LargestAxis)
        .def("IsValid", &AABB::IsValid)
        .def("Contains", &AABB::Contains);

    class_<Placement, std::shared_ptr<Placement>>(m, "Placement")
        .def(init<>())
        .def(init<siren::math::Vector3D const &>())
        .def(init<siren::math::Quaternion const &>())
        .def(init<siren::math::Vector3D const &, siren::math::Quaternion const &>())
        .def(init<Placement const &>())
        .def_property("Position", &Placement::GetPosition, &Placement::SetPosition)
        .def_property("Quaternion", &Placement::GetQuaternion, &Placement::SetQuaternion);
}
