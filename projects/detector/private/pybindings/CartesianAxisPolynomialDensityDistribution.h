#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/EarthModel.h"
#include "../../public/LeptonInjector/detector/CartesianAxisPolynomialDensityDistribution.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

void register_CartesianAxisPolynomialDensityDistribution(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    typedef CartesianAxis1D AxisT;
    typedef PolynomialDistribution1D DistributionT;
    typedef DensityDistribution1D<AxisT,DistributionT> DDist1DT;


    std::string name = "CartesianAxisPolynomialDensityDistribution";

    class_<DDist1DT, std::shared_ptr<DDist1DT>, DensityDistribution>(m, name.c_str())
        .def(init<>())
        .def(init<AxisT, DistributionT>())
        .def(init<DDist1DT>())
        .def("_compare", &DDist1DT::compare)
        .def("_clone", &DDist1DT::clone)
        .def("_create", &DDist1DT::create)
        .def("Derivative", &DDist1DT::Derivative)
        .def("AntiDerivative", &DDist1DT::AntiDerivative)
        .def("Integral", (
                    double (DDist1DT::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double) const)(&DDist1DT::Integral)
                    )
        .def("Integral", (
                    double (DDist1DT::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &) const)(&DDist1DT::Integral)
                    )
        .def("InverseIntegral", (
                    double (DDist1DT::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        double) const)(&DDist1DT::InverseIntegral)
                    )
        .def("InverseIntegral", (
                    double (DDist1DT::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        double,
                        double) const)(&DDist1DT::InverseIntegral)
                    )
        .def("Evaluate", &DDist1DT::Evaluate)
        ;
}

