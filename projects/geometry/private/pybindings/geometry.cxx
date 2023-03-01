
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <tuple>

#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Coordinates.h>
#include <LeptonInjector/Random.h>

#include <LeptonInjector/geometry/Geometry.h>
#include <LeptonInjector/geometry/ExtrPoly.h>

#include <pybind11/pybind11.h>


#include "container_conversions.h"

using namespace pybind11;

PYBIND11_MODULE(Geometry,m) {
  using namespace LI::geometry;

  class_<Geometry>(m, "Geometry")
    .def(init<const std::string>())
    .def(init<const std::string, Placement const &>())
    .def(init<Placement const &>())
    .def(init<const Geometry&>())
    .def("IsInside",&Geometry::IsInside)
    .def("IsInfront",&Geometry::IsInfront)
    .def("IsBehind",&Geometry::IsBehind)
    .def("DistanceToBorder",&Geometry::DistanceToBorder)
    .def("Intersections",&Geometry::Intersections)
    .def("DistanceToClosestApproach",&Geometry::DistanceToClosestApproach)
    .def("GetLocation",&Geometry::GetLocation)
    .def_property_readonly("name",&Geometry::GetName)
    .def_property("placement",&Geometry::GetPlacement, &Geometry::SetPlacement)
    .def("ComputeDistanceToBorder",&Geometry::ComputeDistanceToBorder);
    .def("ComputeIntersections",&Geometry::ComputeIntersections);

  class_<ExtrPoly, Geometry>(m, "ExtrPoly")
    .def(init<const std::vector<std::vector<double>>&,
              const std::vector<ExtrPoly::ZSection>&>>())
    .def(init<Placement const &,
              const std::vector<std::vector<double>>&,
              const std::vector<ExtrPoly::ZSection>&>>())
}
