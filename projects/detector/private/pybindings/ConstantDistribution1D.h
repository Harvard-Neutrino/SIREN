#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/detector/DetectorModel.h"
#include "../../public/SIREN/detector/Distribution1D.h"
#include "../../public/SIREN/detector/ConstantDistribution1D.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"

void register_ConstantDistribution1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::detector;

    class_<ConstantDistribution1D, std::shared_ptr<ConstantDistribution1D>>(m, "ConstantDistribution1D")
        .def(init<>())
        .def(self == self)
        .def(self != self)
        .def("_compare", &ConstantDistribution1D::compare)
        .def("_clone", &ConstantDistribution1D::clone)
        .def("_create", &ConstantDistribution1D::create)
        .def("Derivative", &ConstantDistribution1D::Derivative)
        .def("AntiDerivative", &ConstantDistribution1D::AntiDerivative)
        .def("Evaluate", &ConstantDistribution1D::Evaluate)
        ;
}

