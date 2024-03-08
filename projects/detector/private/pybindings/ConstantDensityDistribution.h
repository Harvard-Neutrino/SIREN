#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/detector/DetectorModel.h"
#include "../../public/SIREN/detector/ConstantDensityDistribution.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"

void register_ConstantDensityDistribution(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace SI::detector;

    typedef CartesianAxis1D AxisT;
    typedef ConstantDistribution1D DistributionT;
    typedef DensityDistribution1D<AxisT,DistributionT> DDist1DT;

    std::string name = "ConstantDensityDistribution";

    class_<DDist1DT, std::shared_ptr<DDist1DT>, DensityDistribution>(m, name.c_str())
        .def(init<>())
        .def(init<AxisT, DistributionT>())
        .def(init<AxisT, double>())
        .def(init<double>())
        .def(init<DDist1DT>())
        .def("_compare", &DDist1DT::compare)
        .def("_clone", &DDist1DT::clone)
        .def("_create", &DDist1DT::create)
        .def("Derivative", &DDist1DT::Derivative)
        .def("AntiDerivative", &DDist1DT::AntiDerivative)
        .def("Integral", (
                    double (DDist1DT::*)(
                        SI::math::Vector3D const &,
                        SI::math::Vector3D const &,
                        double) const)(&DDist1DT::Integral)
                    )
        .def("Integral", (
                    double (DDist1DT::*)(
                        SI::math::Vector3D const &,
                        SI::math::Vector3D const &) const)(&DDist1DT::Integral)
                    )
        .def("InverseIntegral", (
                    double (DDist1DT::*)(
                        SI::math::Vector3D const &,
                        SI::math::Vector3D const &,
                        double,
                        double) const)(&DDist1DT::InverseIntegral)
                    )
        .def("InverseIntegral", (
                    double (DDist1DT::*)(
                        SI::math::Vector3D const &,
                        SI::math::Vector3D const &,
                        double,
                        double,
                        double) const)(&DDist1DT::InverseIntegral)
                    )
        .def("Evaluate", &DDist1DT::Evaluate)
        ;
}


