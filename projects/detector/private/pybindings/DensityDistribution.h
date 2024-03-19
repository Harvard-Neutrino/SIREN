#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/detector/DetectorModel.h"
#include "../../public/SIREN/detector/DensityDistribution.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"

using namespace pybind11;
using namespace siren::detector;
class PyDensityDistribution : public siren::detector::DensityDistribution {
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

    virtual DensityDistribution * clone() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            DensityDistribution *,
            DensityDistribution,
            "_clone",
            clone
        );
    }
    virtual std::shared_ptr<DensityDistribution> create() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::shared_ptr<DensityDistribution>,
            DensityDistribution,
            "_create",
            create
        );
    }

    double Derivative(siren::math::Vector3D const & xi,
            siren::math::Vector3D const & direction) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            Derivative,
            xi,
            direction
        );
    }
    double AntiDerivative(siren::math::Vector3D const & xi,
            siren::math::Vector3D const & direction) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            AntiDerivative,
            xi,
            direction
        );
    }
    double Integral(siren::math::Vector3D const & xi,
            siren::math::Vector3D const & direction,
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
    double Integral(siren::math::Vector3D const & xi,
            siren::math::Vector3D const & xj) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            DensityDistribution,
            Integral,
            xi,
            xj
        );
    }
    double InverseIntegral(siren::math::Vector3D const & xi,
            siren::math::Vector3D const & direction,
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
    double InverseIntegral(siren::math::Vector3D const & xi,
            siren::math::Vector3D const & direction,
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
    double Evaluate(const siren::math::Vector3D& xi) const override {
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
    using namespace siren::detector;

    class_<DensityDistribution, std::shared_ptr<DensityDistribution>, PyDensityDistribution>(m, "DensityDistribution")
        .def(init<>())
        .def("__eq__", [](const DensityDistribution &self, const DensityDistribution &other){ return self == other; })
        .def("__ne__", [](const DensityDistribution &self, const DensityDistribution &other){ return self != other; })
        .def("_compare", &DensityDistribution::compare)
        .def("_clone", &DensityDistribution::clone)
        .def("_create", &DensityDistribution::create)
        .def("Derivative", &DensityDistribution::Derivative)
        .def("AntiDerivative", &DensityDistribution::AntiDerivative)
        .def("Integral", (
                    double (DensityDistribution::*)(
                        siren::math::Vector3D const &,
                        siren::math::Vector3D const &,
                        double) const)(&DensityDistribution::Integral)
                    )
        .def("Integral", (
                    double (DensityDistribution::*)(
                        siren::math::Vector3D const &,
                        siren::math::Vector3D const &) const)(&DensityDistribution::Integral)
                    )
        .def("InverseIntegral", (
                    double (DensityDistribution::*)(
                        siren::math::Vector3D const &,
                        siren::math::Vector3D const &,
                        double,
                        double) const)(&DensityDistribution::InverseIntegral)
                    )
        .def("InverseIntegral", (
                    double (DensityDistribution::*)(
                        siren::math::Vector3D const &,
                        siren::math::Vector3D const &,
                        double,
                        double,
                        double) const)(&DensityDistribution::InverseIntegral)
                    )
        .def("Evaluate", &DensityDistribution::Evaluate)
        ;
}
