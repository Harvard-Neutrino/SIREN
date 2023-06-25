#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/EarthModel.h"
#include "../../public/LeptonInjector/detector/DensityDistribution.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

using namespace pybind11;
using namespace LI::detector;
class PyDensityDistribution : public LI::detector::DensityDistribution {
public:
    using DensityDistribution::DensityDistribution;

    bool compare(const DensityDistribution& density_distr) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            bool,
            DensityDistribution,
            "_compare",
            compare,
            density_distr
        );
    }

    DensityDistribution * clone() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            DensityDistribution *,
            DensityDistribution,
            "_clone",
            clone
        );
    }
    std::shared_ptr<const DensityDistribution> create() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::shared_ptr<const DensityDistribution>,
            DensityDistribution,
            "_create",
            create
        );
    }

    double Derivative(LI::math::Vector3D const & xi,
            LI::math::Vector3D const & direction) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            Derivative,
            xi,
            direction
        );
    }
    double AntiDerivative(LI::math::Vector3D const & xi,
            LI::math::Vector3D const & direction) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            AntiDerivative,
            xi,
            direction
        );
    }
    double Integral(LI::math::Vector3D const & xi,
            LI::math::Vector3D const & direction,
            double distance) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            Integral,
            xi,
            direction,
            distance
        );
    }
    double Integral(LI::math::Vector3D const & xi,
            LI::math::Vector3D const & xj) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            Integral,
            xi,
            xj
        );
    }
    double InverseIntegral(LI::math::Vector3D const & xi,
            LI::math::Vector3D const & direction,
            double integral,
            double max_distance) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            InverseIntegral,
            xi,
            direction,
            integral,
            max_distance
        );
    }
    double InverseIntegral(LI::math::Vector3D const & xi,
            LI::math::Vector3D const & direction,
            double constant,
            double integral,
            double max_distance) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            InverseIntegral,
            xi,
            direction,
            constant,
            integral,
            max_distance
        );
    }
    double Evaluate(const LI::math::Vector3D& xi) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            Evaluate,
            xi
        );
    }
};

void register_DensityDistribution(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<DensityDistribution, std::shared_ptr<DensityDistribution>, PyDensityDistribution>(m, "DensityDistribution")
        .def(init<>())
        .def(self == self)
        .def(self != self)
        .def("_compare", &DensityDistribution::compare)
        .def("_clone", &DensityDistribution::clone)
        .def("_create", &DensityDistribution::create)
        .def("Derivative", &DensityDistribution::Derivative)
        .def("AntiDerivative", &DensityDistribution::AntiDerivative)
        .def("Integral", (
                    double (DensityDistribution::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double) const)(&DensityDistribution::Integral)
                    )
        .def("Integral", (
                    double (DensityDistribution::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &) const)(&DensityDistribution::Integral)
                    )
        .def("InverseIntegral", (
                    double (DensityDistribution::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        double) const)(&DensityDistribution::InverseIntegral)
                    )
        .def("InverseIntegral", (
                    double (DensityDistribution::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        double,
                        double) const)(&DensityDistribution::InverseIntegral)
                    )
        .def("Evaluate", &DensityDistribution::Evaluate)
        ;
}
