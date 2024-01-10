#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/Distribution1D.h"
#include "../../public/LeptonInjector/detector/ExponentialDistribution1D.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

void register_ExponentialDistribution1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<ExponentialDistribution1D, std::shared_ptr<ExponentialDistribution1D>>(m, "ExponentialDistribution1D")
        .def(init<>())
        .def(init<ExponentialDistribution1D>())
        .def(init<double>())
        .def("_compare", &ExponentialDistribution1D::compare)
        .def("_clone", &ExponentialDistribution1D::clone)
        .def("_create", &ExponentialDistribution1D::create)
        .def("Derivative", &ExponentialDistribution1D::Derivative)
        .def("AntiDerivative", &ExponentialDistribution1D::AntiDerivative)
        .def("Evaluate", &ExponentialDistribution1D::Evaluate)
        ;
}

