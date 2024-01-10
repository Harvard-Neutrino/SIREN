#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/Distribution1D.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

using namespace pybind11;
using namespace LI::detector;
class PyDistribution1D : public LI::detector::Distribution1D {
public:
    using Distribution1D::Distribution1D;

    bool compare(const Distribution1D& density_distr) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            bool,
            Distribution1D,
            "_compare",
            compare,
            density_distr
        );
    }

    Distribution1D * clone() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            Distribution1D *,
            Distribution1D,
            "_clone",
            clone
        );
    }
    std::shared_ptr<const Distribution1D> create() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::shared_ptr<const Distribution1D>,
            Distribution1D,
            "_create",
            create
        );
    }

    double Derivative(double x) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            Distribution1D,
            Derivative,
            x
        );
    }
    double AntiDerivative(double x) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            Distribution1D,
            AntiDerivative,
            x
        );
    }
    double Evaluate(double x) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            Distribution1D,
            Evaluate,
            x
        );
    }
};

void register_Distribution1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<Distribution1D, std::shared_ptr<Distribution1D>, PyDistribution1D>(m, "Distribution1D")
        .def(init<>())
        .def("__eq__", [](const Distribution1D &self, const Distribution1D &other){ return self == other; })
        .def("__ne__", [](const Distribution1D &self, const Distribution1D &other){ return self != other; }) 
        .def("_compare", &Distribution1D::compare)
        .def("_clone", &Distribution1D::clone)
        .def("_create", &Distribution1D::create)
        .def("Derivative", &Distribution1D::Derivative)
        .def("AntiDerivative", &Distribution1D::AntiDerivative)
        .def("Evaluate", &Distribution1D::Evaluate)
        ;
}
