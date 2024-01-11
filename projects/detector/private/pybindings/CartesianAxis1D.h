#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/CartesianAxis1D.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

using namespace pybind11;
using namespace LI::detector;

void register_CartesianAxis1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<CartesianAxis1D, Axis1D, std::shared_ptr<CartesianAxis1D>>(m, "CartesianAxis1D")
        .def(init<>())
        .def(init<LI::math::Vector3D const &, LI::math::Vector3D const &>())
        .def("_compare", &CartesianAxis1D::compare)
        .def("_clone", &CartesianAxis1D::clone)
        .def("_create", &CartesianAxis1D::create)
        .def("GetX", &CartesianAxis1D::GetX)
        .def("GetdX", &CartesianAxis1D::GetdX)
        ;
}
