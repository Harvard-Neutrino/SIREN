#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/Distribution1D.h"
#include "../../public/LeptonInjector/detector/PolynomialDistribution1D.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

void register_PolynomialDistribution1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<PolynomialDistribution1D, std::shared_ptr<PolynomialDistribution1D>>(m, "PolynomialDistribution1D")
        .def(init<>())
        .def(init<PolynomialDistribution1D>())
        .def(init<LI::math::Polynom>())
        .def(init<std::vector<double>>())
        .def("_compare", &PolynomialDistribution1D::compare)
        .def("_clone", &PolynomialDistribution1D::clone)
        .def("_create", &PolynomialDistribution1D::create)
        .def("Derivative", &PolynomialDistribution1D::Derivative)
        .def("AntiDerivative", &PolynomialDistribution1D::AntiDerivative)
        .def("Evaluate", &PolynomialDistribution1D::Evaluate)
        ;
}

