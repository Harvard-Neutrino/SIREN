#include <vector>
#include <optional>
#include "SIREN/math/FourMomentum.h"

#include "../../public/SIREN/math/Vector3D.h"
#include "../../public/SIREN/math/Quaternion.h"
#include "../../public/SIREN/math/Matrix3D.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

using namespace pybind11;

PYBIND11_MODULE(math,m) {
    using namespace siren::math;

    class_<Vector3D, std::shared_ptr<Vector3D>>(m, "Vector3D")
        .def(init<>())
        .def(init<const double, const double, const double>())
        .def(init<const Vector3D&>())
        .def(init<std::array<double, 3> const &>())
        .def("__str__", [](Vector3D const & v) { std::stringstream ss; ss << v; return ss.str(); })
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
        .def("GetTheta",&Vector3D::GetTheta)
        .def_static("scalar_product", [](Vector3D const & a, Vector3D const & b)->double{return scalar_product(a, b);})
        .def_static("vector_product", [](Vector3D const & a, Vector3D const & b)->Vector3D{return vector_product(a, b);})
        .def_static("cross_product", [](Vector3D const & a, Vector3D const & b)->Vector3D{return cross_product(a, b);});

    class_<FourMomentum>(m, "FourMomentum")
        .def_static("from_components", &FourMomentum::FromComponents, arg("components"))
        .def_static("from_mass", &FourMomentum::FromMass,
                    arg("mass"), arg("momentum"), arg("negative_energy") = false)
        .def("components", &FourMomentum::Components)
        .def("energy", &FourMomentum::Energy)
        .def("momentum", &FourMomentum::Momentum)
        .def("mass", &FourMomentum::Mass)
        .def("mass_squared", &FourMomentum::MassSquared)
        .def("has_known_mass", &FourMomentum::HasKnownMass)
        .def("boost_to_lab", &FourMomentum::BoostToLab, arg("frame"))
        .def("boost_to_rest", &FourMomentum::BoostToRest, arg("frame"))
        .def_static("invariant_mass", &FourMomentum::InvariantMass, arg("a"), arg("b"))
        .def(self + self)
        .def(self - self)
        .def(pickle(
            [](FourMomentum const & p) {
                std::optional<double> mass;
                if (p.HasKnownMass()) mass = p.Mass();
                return pybind11::make_tuple(p.Components(), mass);
            },
            [](tuple state) {
                if (state.size() != 2)
                    throw std::runtime_error("Invalid FourMomentum pickle state");
                auto components = state[0].cast<std::array<double, 4>>();
                auto mass = state[1].cast<std::optional<double>>();
                if (!mass) return FourMomentum::FromComponents(components);
                auto p = FourMomentum::FromMass(*mass,
                    Vector3D(components[1], components[2], components[3]),
                    std::signbit(components[0]));
                if (p.Components() != components)
                    throw std::runtime_error("FourMomentum pickle energy disagrees with its mass");
                return p;
            }));

    class_<Matrix3D, std::shared_ptr<Matrix3D>>(m, "Matrix3D")
        .def(init<>())
        .def(init<double,double,double,double,double,double,double,double,double>(),
             arg("xx"), arg("xy"), arg("xz"),
             arg("yx"), arg("yy"), arg("yz"),
             arg("zx"), arg("zy"), arg("zz"))
        .def(init<const Matrix3D &>())
        .def("__str__", [](Matrix3D const & m) { std::stringstream ss; ss << m; return ss.str(); })
        .def(self == self)
        .def(self != self);

    class_<Quaternion, std::shared_ptr<Quaternion>>(m, "Quaternion")
        .def(init<>())
        .def(init<const double, const double, const double, const double>())
        .def(init<const Quaternion &>())
        .def(init<const Vector3D &>())
        //.def(init<geom3::Rotation3::Quaternion const &>())
        .def("__str__", [](Quaternion const & q) { std::stringstream ss; ss << q; return ss.str(); })
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
        .def("SetAxisAngle", &Quaternion::SetAxisAngle)
        .def("GetAxisAngle", (void (Quaternion::*)(Vector3D &, double &) const)(&Quaternion::GetAxisAngle))
        .def("GetAxisAngle", (std::tuple<Vector3D, double> (Quaternion::*)()const)(&Quaternion::GetAxisAngle))
        .def("GetEulerAngles", &Quaternion::GetEulerAngles)
        .def("GetEulerAnglesZXZr", [](Quaternion const & q) {
            double a, b, g; q.GetEulerAnglesZXZr(a, b, g);
            return std::make_tuple(a, b, g);
        })
        .def("GetEulerAnglesXYZs", [](Quaternion const & q) {
            double a, b, g; q.GetEulerAnglesXYZs(a, b, g);
            return std::make_tuple(a, b, g);
        })
        .def("SetEulerAngles", &Quaternion::SetEulerAngles)
        .def("SetEulerAnglesZXZr", &Quaternion::SetEulerAnglesZXZr)
        .def("SetEulerAnglesXYZs", &Quaternion::SetEulerAnglesXYZs)
        .def("rotate", (Quaternion (Quaternion::*)(Quaternion const &, bool) const)(&Quaternion::rotate))
        .def("rotate", (Vector3D (Quaternion::*)(Vector3D const &, bool) const)(&Quaternion::rotate))
        .def_property("X", &Quaternion::GetX, &Quaternion::SetX)
        .def_property("Y", &Quaternion::GetY, &Quaternion::SetY)
        .def_property("Z", &Quaternion::GetZ, &Quaternion::SetZ)
        .def_property("W", &Quaternion::GetW, &Quaternion::SetW)
        .def_static("lerp", &Quaternion::lerp)
        .def_static("slerp", &Quaternion::slerp)
        .def_static("rotation_between", [](Vector3D const & a, Vector3D const & b)->Quaternion{return rotation_between(a, b);});
}
