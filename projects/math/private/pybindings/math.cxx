#include <vector>

#include "../../public/LeptonInjector/math/Vector3D.h"
#include "../../public/LeptonInjector/math/Quaternion.h"
#include "../../public/LeptonInjector/math/Matrix3D.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

using namespace pybind11;

PYBIND11_MODULE(math,m) {
    using namespace LI::math;

    class_<Vector3D, std::shared_ptr<Vector3D>>(m, "Vector3D")
        .def(init<>())
        .def(init<const double, const double, const double>())
        .def(init<const Vector3D&>())
        .def(init<std::array<double, 3> const &>())
        .def(self == self)
        .def(self != self)
        .def(self < self)
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

    class_<Quaternion, std::shared_ptr<Quaternion>>(m, "Quaternion")
        .def(init<>())
        .def(init<const double, const double, const double, const double>())
        .def(init<const Quaternion &>())
        .def(init<const Vector3D &>())
        //.def(init<geom3::Rotation3::Quaternion const &>())
        .def(self == self)
        .def(self != self)
        .def(self < self)
        .def(self * self)
        .def(self *= self)
        .def(self * double())
        .def(self *= double())
        .def(self + self)
        .def(self += self)
        .def(self + double())
        .def(self += double())
        .def("__invert__", [](Quaternion const & q)->Quaternion{return ~q;})
        .def("__not__", [](Quaternion const & q)->Quaternion{return !q;})
        .def("SetPosition", &Quaternion::SetPosition)
        .def("GetMatrix", (void     (Quaternion::*)(Matrix3D &) const)(&Quaternion::GetMatrix))
        .def("GetMatrix", (Matrix3D (Quaternion::*)(          ) const)(&Quaternion::GetMatrix))
        .def("SetMatrix", &Quaternion::SetMatrix)
        .def("invert", &Quaternion::invert)
        .def("inverted", &Quaternion::inverted)
        .def("conjugate", &Quaternion::conjugate)
        .def("conjugated", &Quaternion::conjugated)
        .def("normalize", &Quaternion::normalize)
        .def("normalized", &Quaternion::normalized)
        .def("magnitude", &Quaternion::magnitude)
        .def("magnitudesq", &Quaternion::magnitudesq)
        .def("DotProduct", &Quaternion::DotProduct)
        .def("lerp", &Quaternion::lerp)
        .def("slerp", &Quaternion::slerp)
        .def("SetAxisAngle", &Quaternion::SetAxisAngle)
        .def("GetAxisAngle", (void (Quaternion::*)(Vector3D &, double &) const)(&Quaternion::GetAxisAngle))
        .def("GetAxisAngle", (std::tuple<Vector3D, double> (Quaternion::*)()const)(&Quaternion::GetAxisAngle))
        .def("GetEulerAngles", &Quaternion::GetEulerAngles)
        .def("GetEulerAnglesZXZr", &Quaternion::GetEulerAnglesZXZr)
        .def("GetEulerAnglesXYZs", &Quaternion::GetEulerAnglesXYZs)
        .def("SetEulerAngles", &Quaternion::SetEulerAngles)
        .def("SetEulerAnglesZXZr", &Quaternion::SetEulerAnglesZXZr)
        .def("SetEulerAnglesXYZs", &Quaternion::SetEulerAnglesXYZs)
        .def("rotate", (Quaternion (Quaternion::*)(Quaternion const &, bool) const)(&Quaternion::rotate))
        .def("rotate", (Vector3D (Quaternion::*)(Vector3D const &, bool) const)(&Quaternion::rotate))
        .def_property("X", &Quaternion::GetX, &Quaternion::SetX)
        .def_property("Y", &Quaternion::GetY, &Quaternion::SetY)
        .def_property("Z", &Quaternion::GetZ, &Quaternion::SetZ)
        .def_property("W", &Quaternion::GetW, &Quaternion::SetW)
        .def_static("rotation_between", [](object, Vector3D const & a, Vector3D const & b)->Quaternion{return rotation_between(a, b);});
}
