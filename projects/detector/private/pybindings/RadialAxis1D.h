#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/detector/DetectorModel.h"
#include "../../public/SIREN/detector/Axis1D.h"
#include "../../public/SIREN/detector/RadialAxis1D.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"

void register_RadialAxis1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::detector;

    class_<RadialAxis1D, std::shared_ptr<RadialAxis1D>>(m, "RadialAxis1D")
        .def(init<>())
        .def(init<siren::math::Vector3D>())
        .def(init<siren::math::Vector3D, siren::math::Vector3D>())
        .def(self == self)
        .def(self != self)
        .def("_compare", &RadialAxis1D::compare)
        .def("_clone", &RadialAxis1D::clone)
        .def("_create", &RadialAxis1D::create)
        .def("GetX", &RadialAxis1D::GetX)
        .def("GetdX", &RadialAxis1D::GetdX)
        ;
}

