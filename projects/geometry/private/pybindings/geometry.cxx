
#include <vector>

#include "../../public/LeptonInjector/geometry/Placement.h"
#include "../../public/LeptonInjector/geometry/Geometry.h"
#include "../../public/LeptonInjector/geometry/ExtrPoly.h"
#include "../../public/LeptonInjector/geometry/Cylinder.h"
#include "../../public/LeptonInjector/geometry/Box.h"
#include "../../public/LeptonInjector/geometry/Sphere.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>



using namespace pybind11;

PYBIND11_MODULE(geometry,m) {
    using namespace LI::geometry;

    // geometry

    class_<Geometry, std::shared_ptr<Geometry>>(m, "Geometry")
        .def("IsInside",&Geometry::IsInside)
        .def("IsInfront",&Geometry::IsInfront)
        .def("IsBehind",&Geometry::IsBehind)
        .def("DistanceToBorder",&Geometry::DistanceToBorder)
        .def("Intersections",&Geometry::Intersections)
        .def("DistanceToClosestApproach",&Geometry::DistanceToClosestApproach)
        .def("GetLocation",&Geometry::GetLocation)
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
        .def(init<double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<const Box&>())
        .def_property("X",&Box::GetX, &Box::SetX)
        .def_property("Y",&Box::GetY, &Box::SetY)
        .def_property("Z",&Box::GetZ, &Box::SetZ);

    // Cylinder

    class_<Cylinder, std::shared_ptr<Cylinder>, Geometry>(m, "Cylinder")
        .def(init<>())
        .def(init<double, double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double, double>())
        .def(init<const Cylinder&>())
        .def_property("InnerRadius",&Cylinder::GetInnerRadius, &Cylinder::SetInnerRadius)
        .def_property("Radius",&Cylinder::GetRadius, &Cylinder::SetRadius)
        .def_property("Z",&Cylinder::GetZ, &Cylinder::SetZ);

    // Sphere

    class_<Sphere, std::shared_ptr<Sphere>, Geometry>(m, "Sphere")
        .def(init<>())
        .def(init<double, double>())
        .def(init<Placement const &>())
        .def(init<Placement const &, double, double>())
        .def(init<const Sphere&>())
        .def_property("InnerRadius",&Sphere::GetInnerRadius, &Sphere::SetInnerRadius)
        .def_property("Radius",&Sphere::GetRadius, &Sphere::SetRadius);

    class_<Placement, std::shared_ptr<Placement>>(m, "Placement")
        .def(init<>())
        .def(init<LI::math::Vector3D const &>())
        .def(init<LI::math::Quaternion const &>())
        .def(init<LI::math::Vector3D const &, LI::math::Quaternion const &>())
        .def(init<Placement const &>())
        .def_property("Position", &Placement::GetPosition, &Placement::SetPosition)
        .def_property("Quaternion", &Placement::GetQuaternion, &Placement::SetQuaternion);
}
