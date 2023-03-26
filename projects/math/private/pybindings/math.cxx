#include <vector>

#include "../../public/LeptonInjector/math/Vector3D.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

using namespace pybind11;

PYBIND11_MODULE(Math,m) {
  using namespace LI::math;

  class_<Vector3D, std::shared_ptr<Vector3D>>(m, "Vector3D")
    .def(init<>())
    .def(init<const double, const double, const double>())
    .def(init<const Vector3D&>())
    .def(init<std::array<double, 3> const &>())
    .def(self + self)
    .def(self += self)
    .def(self - self)
    .def(self -= self)
    .def(self * double())
    .def(self *= double())
    .def(double() * self)
    .def(self / double())
    .def(self /= double())
    .def("magnitude",&Vector3D::magnitude)
    .def("normalize",&Vector3D::normalize)
    .def("normalized",&Vector3D::normalized)
    .def("deflect",&Vector3D::deflect)
    .def("invert",&Vector3D::invert)
    .def("inverted",&Vector3D::inverted)
    .def("CalculateCartesianFromSpherical",&Vector3D::CalculateCartesianFromSpherical)
    .def("CalculateSphericalCoordinates",&Vector3D::CalculateSphericalCoordinates)
    .def("GetX",&Vector3D::GetX)
    .def("GetY",&Vector3D::GetY)
    .def("GetZ",&Vector3D::GetZ)
    .def("GetRadius",&Vector3D::GetRadius)
    .def("GetPhi",&Vector3D::GetPhi)
    .def("GetTheta",&Vector3D::GetTheta);
}
