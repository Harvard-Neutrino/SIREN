
#include <vector>

#include "../../public/LeptonInjector/geometry/Placement.h"
#include "../../public/LeptonInjector/geometry/Geometry.h"
#include "../../public/LeptonInjector/geometry/ExtrPoly.h"
#include "../../public/LeptonInjector/geometry/Cylinder.h"
#include "../../public/LeptonInjector/geometry/Box.h"
#include "../../public/LeptonInjector/geometry/Sphere.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(Geometry,m) {
  using namespace LI::geometry;

  class_<Geometry>(m, "Geometry")
    //.def(init<>())
    //.def(init<const std::string>())
    //.def(init<const std::string, Placement const &>())
    //.def(init<Placement const &>())
    //.def(init<const Geometry&>())
    .def("IsInside",&Geometry::IsInside)
    .def("IsInfront",&Geometry::IsInfront)
    .def("IsBehind",&Geometry::IsBehind)
    .def("DistanceToBorder",&Geometry::DistanceToBorder)
    .def("Intersections",&Geometry::Intersections)
    .def("DistanceToClosestApproach",&Geometry::DistanceToClosestApproach)
    .def("GetLocation",&Geometry::GetLocation)
    .def_property_readonly("name",&Geometry::GetName)
    .def_property("placement",&Geometry::GetPlacement, &Geometry::SetPlacement)
    //.def("ComputeDistanceToBorder",&Geometry::ComputeDistanceToBorder)
    .def("ComputeIntersections",&Geometry::ComputeIntersections)
    .def("create",&Geometry::create);

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
  
  class_<ExtrPoly, Geometry>(m, "ExtrPoly")
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

  class_<Box, Geometry>(m, "Box")
    .def(init<>())
    .def(init<double, double, double>())
    .def(init<Placement const &>())
    .def(init<Placement const &, double, double, double>())
    .def(init<const Box&>())
    .def_property("X",&Box::GetX, &Box::SetX)
    .def_property("Y",&Box::GetY, &Box::SetY)
    .def_property("Z",&Box::GetZ, &Box::SetZ);
  
  // Cylinder

  class_<Cylinder, Geometry>(m, "Cylinder")
    .def(init<>())
    .def(init<double, double, double>())
    .def(init<Placement const &>())
    .def(init<Placement const &, double, double, double>())
    .def(init<const Cylinder&>())
    .def_property("InnerRadius",&Cylinder::GetInnerRadius, &Cylinder::SetInnerRadius)
    .def_property("Radius",&Cylinder::GetRadius, &Cylinder::SetRadius)
    .def_property("Z",&Cylinder::GetZ, &Cylinder::SetZ);
  
  // Sphere

  class_<Sphere, Geometry>(m, "Sphere")
    .def(init<>())
    .def(init<double, double>())
    .def(init<Placement const &>())
    .def(init<Placement const &, double, double>())
    .def(init<const Sphere&>())
    .def_property("InnerRadius",&Sphere::GetInnerRadius, &Sphere::SetInnerRadius)
    .def_property("Radius",&Sphere::GetRadius, &Sphere::SetRadius);


}
