#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/EarthModel.h"
#include "../../public/LeptonInjector/detector/Axis1D.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

using namespace pybind11;
using namespace LI::detector;
class PyAxis1D : public LI::detector::Axis1D {
public:
    using Axis1D::Axis1D;

    bool compare(const Axis1D& density_distr) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            bool,
            Axis1D,
            "_compare",
            compare,
            density_distr
        );
    }

    Axis1D * clone() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            Axis1D *,
            Axis1D,
            "_clone",
            clone
        );
    }
    std::shared_ptr<const Axis1D> create() const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            std::shared_ptr<const Axis1D>,
            Axis1D,
            "_create",
            create
        );
    }

    double GetX(const LI::math::Vector3D& xi) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            Axis1D,
            GetX,
            xi
        );
    }

    double GetdX(const LI::math::Vector3D& xi, const LI::math::Vector3D& direction) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            Axis1D,
            GetdX,
            xi,
            direction
        );
    }
};

void register_Axis1D(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<Axis1D, std::shared_ptr<Axis1D>, PyAxis1D>(m, "Axis1D")
        .def(init<>())
        .def(self == self)
        .def(self != self)
        .def("_compare", &Axis1D::compare)
        .def("_clone", &Axis1D::clone)
        .def("_create", &Axis1D::create)
        .def("GetX", &Axis1D::GetX)
        .def("GetdX", &Axis1D::GetdX)
        .def("GetAxis", &Axis1D::GetAxis)
        .def("GetFp0", &Axis1D::GetFp0)
        ;
}
