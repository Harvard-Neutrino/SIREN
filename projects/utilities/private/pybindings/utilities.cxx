
#include <vector>

#include "../../public/SIREN/utilities/Random.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost
#include <rk/LT.hh>                                           // for LT

using namespace pybind11;

PYBIND11_MODULE(utilities,m) {
  using namespace siren::utilities;

  class_<SIREN_random, std::shared_ptr<SIREN_random>>(m, "SIREN_random")
    .def(init<>())
    .def(init<unsigned int>())
    .def("Uniform",&SIREN_random::Uniform)
    .def("set_seed",&SIREN_random::set_seed)
    .def("get_seed",&SIREN_random::get_seed)
    .def_static("generate_seed",&SIREN_random::generate_seed);

    auto py_rk = m.def_submodule("rk", "Bindings for rk types");

    py_rk.doc() = "Python bindings for the rk relativistic kinematics library";

    // --- Minimal bindings for the geom3 types ---
    // These are only placeholders. In practice, you should provide full bindings.
    auto py_geom3 = py_rk.def_submodule("geom3", "Bindings for rk/geom3.hh types");

    py_geom3.doc() = "Python bindings for the geom3 library";

    // =============================
    // Bindings for UnitVector3
    // =============================
    class_<geom3::UnitVector3>(py_geom3, "UnitVector3")
        .def(init<double, double>(), arg("theta"), arg("phi"))
        .def(init<double, double, double, bool>(),
             arg("x"), arg("y"), arg("z"), arg("mustNormalize") = true)
        .def(init<double>(), arg("phi"),
             "Fast constructor for directions in the XY plane (z is set to 0).")
        .def("setTheta", &geom3::UnitVector3::setTheta, "Set the polar angle theta.")
        .def("setPhi", &geom3::UnitVector3::setPhi, "Set the azimuthal angle phi.")
        .def("setEta", &geom3::UnitVector3::setEta, "Set the pseudorapidity eta.")
        .def("x", &geom3::UnitVector3::x)
        .def("y", &geom3::UnitVector3::y)
        .def("z", &geom3::UnitVector3::z)
        .def("theta", &geom3::UnitVector3::theta)
        .def("cosTheta", &geom3::UnitVector3::cosTheta)
        .def("phi", &geom3::UnitVector3::phi)
        .def("eta", &geom3::UnitVector3::eta)
        // Overloads for dot, cross and angle methods:
        .def("dot", (double (geom3::UnitVector3::*)(const geom3::UnitVector3&) const)
             &geom3::UnitVector3::dot, "Dot product with another UnitVector3.")
        .def("dot_vec", (double (geom3::UnitVector3::*)(const geom3::Vector3&) const)
             &geom3::UnitVector3::dot, "Dot product with a Vector3.")
        .def("cross", (geom3::Vector3 (geom3::UnitVector3::*)(const geom3::UnitVector3&) const)
             &geom3::UnitVector3::cross, "Cross product with another UnitVector3.")
        .def("cross_vec", (geom3::Vector3 (geom3::UnitVector3::*)(const geom3::Vector3&) const)
             &geom3::UnitVector3::cross, "Cross product with a Vector3.")
        .def("angle", (double (geom3::UnitVector3::*)(const geom3::UnitVector3&) const)
             &geom3::UnitVector3::angle, "Angle between two UnitVector3 objects.")
        .def("angle_vec", (double (geom3::UnitVector3::*)(const geom3::Vector3&) const)
             &geom3::UnitVector3::angle, "Angle between a UnitVector3 and a Vector3.")
        // Unary operators
        .def(-self)
        .def(+self)
        // Comparison and scalar operators
        .def(self == self)
        .def(self != self)
        .def(self * double())
        .def(double() * self)
        .def(self / double())
        // Static member functions
        .def_static("xAxis", &geom3::UnitVector3::xAxis, "Return the unit vector along the x-axis.")
        .def_static("yAxis", &geom3::UnitVector3::yAxis, "Return the unit vector along the y-axis.")
        .def_static("zAxis", &geom3::UnitVector3::zAxis, "Return the unit vector along the z-axis.")
        .def_static("random", &geom3::UnitVector3::random, arg("rnd0"), arg("rnd1"),
                    "Generate a random direction uniformly distributed in 4pi solid angle.")
        ;

    // =============================
    // Bindings for Vector3
    // =============================
    class_<geom3::Vector3>(py_geom3, "Vector3")
        .def(init<>())
        .def(init<double, double, double>(), arg("x"), arg("y"), arg("z"))
        .def(init<const std::array<double, 3>&>(), arg("vec"))
        .def(init<double, const geom3::UnitVector3&>(),
             arg("length"), arg("direction"),
             "Construct a Vector3 given a length and a UnitVector3 for the direction.")
        .def(init<const geom3::Point3&, const geom3::Point3&>(),
             arg("p"), arg("origin"),
             "Construct a vector as the difference between two Point3 objects.")
        .def("setLength", &geom3::Vector3::setLength, "Set the vector length (absolute value is used).")
        .def("set", &geom3::Vector3::set, "Set a coordinate value by index.")
        .def("x", &geom3::Vector3::x)
        .def("y", &geom3::Vector3::y)
        .def("z", &geom3::Vector3::z)
        .def("length", &geom3::Vector3::length)
        .def("lengthSquared", &geom3::Vector3::lengthSquared)
        .def("direction", &geom3::Vector3::direction, "Return the unit vector in the direction of this vector.")
        .def("theta", &geom3::Vector3::theta)
        .def("cosTheta", &geom3::Vector3::cosTheta)
        .def("phi", &geom3::Vector3::phi)
        .def("eta", &geom3::Vector3::eta)
        // Overloads for dot, cross and angle methods:
        .def("dot", (double (geom3::Vector3::*)(const geom3::Vector3&) const)
             &geom3::Vector3::dot, "Dot product with another Vector3.")
        .def("dot_uv", (double (geom3::Vector3::*)(const geom3::UnitVector3&) const)
             &geom3::Vector3::dot, "Dot product with a UnitVector3.")
        .def("cross", (geom3::Vector3 (geom3::Vector3::*)(const geom3::Vector3&) const)
             &geom3::Vector3::cross, "Cross product with another Vector3.")
        .def("cross_uv", (geom3::Vector3 (geom3::Vector3::*)(const geom3::UnitVector3&) const)
             &geom3::Vector3::cross, "Cross product with a UnitVector3.")
        .def("angle", (double (geom3::Vector3::*)(const geom3::Vector3&) const)
             &geom3::Vector3::angle, "Angle between two vectors.")
        .def("angle_uv", (double (geom3::Vector3::*)(const geom3::UnitVector3&) const)
             &geom3::Vector3::angle, "Angle between a Vector3 and a UnitVector3.")
        // Operators
        .def(-self)
        .def(+self)
        .def(self == self)
        .def(self != self)
        .def(self * double())
        .def(double() * self)
        .def(self / double())
        .def(self + self)
        .def(self - self)
        .def(self += self)
        .def(self -= self)
        .def(self *= double())
        .def(self /= double())
        // Provide index access (read-only)
        .def("__getitem__", [](const geom3::Vector3 &v, unsigned i) { return v[i]; })
        ;

    // =============================
    // Bindings for Point3
    // =============================
    class_<geom3::Point3>(py_geom3, "Point3")
        .def(init<>())
        .def(init<double, double, double>(), arg("x"), arg("y"), arg("z"))
        .def(init<const std::array<double, 3>&>(), arg("p"))
        .def("set", &geom3::Point3::set, "Set a coordinate value by index.")
        .def("x", &geom3::Point3::x)
        .def("y", &geom3::Point3::y)
        .def("z", &geom3::Point3::z)
        .def(self == self)
        .def(self != self)
        // The following operators allow addition/subtraction with Vector3 objects.
        .def(self - self)
        //.def(self + self)
        .def("__getitem__", [](const geom3::Point3 &p, unsigned i) { return p[i]; })
        ;

    // =============================
    // Bindings for Matrix3x3
    // =============================
    class_<geom3::Matrix3x3>(py_geom3, "Matrix3x3")
        .def(init<>(), "Default constructor makes a unit matrix.")
        .def(init<const double*>(), arg("data"),
             "Construct from a C-style array (row-by-row).")
        .def(init<double, double, double,
                      double, double, double,
                      double, double, double>(),
             arg("m00"), arg("m01"), arg("m02"),
             arg("m10"), arg("m11"), arg("m12"),
             arg("m20"), arg("m21"), arg("m22"))
        .def(init<const geom3::Vector3&, const geom3::Vector3&, const geom3::Vector3&>(),
             arg("row0"), arg("row1"), arg("row2"))
        .def(init<const geom3::UnitVector3&, const geom3::UnitVector3&, const geom3::UnitVector3&>(),
             arg("row0"), arg("row1"), arg("row2"))
        .def("set", &geom3::Matrix3x3::set, "Set an element of the matrix by indices.")
        .def("__getitem__", [](const geom3::Matrix3x3 &m, unsigned i) -> const geom3::Vector3& {
            return m[i];
        }, return_value_policy::reference_internal)
        .def("fillCArray", [](const geom3::Matrix3x3 &m) {
            std::array<double, 9> arr;
            m.fillCArray(arr.data());
            return arr;
        }, "Fill a C-style array (row-by-row) with matrix elements.")
        .def("fillFArray", [](const geom3::Matrix3x3 &m) {
            std::array<double, 9> arr;
            m.fillFArray(arr.data());
            return arr;
        }, "Fill a Fortran-style array (column-by-column) with matrix elements.")
        .def("T", &geom3::Matrix3x3::T, "Return the transposed matrix.")
        .def("inverse", &geom3::Matrix3x3::inverse, "Return the inverse matrix.")
        .def("det", &geom3::Matrix3x3::det, "Return the determinant.")
        .def("tr", &geom3::Matrix3x3::tr, "Return the trace of the matrix.")
        // Operators
        .def(-self)
        .def(+self)
        .def(self == self)
        .def(self != self)
        .def(self * double())
        .def(double() * self)
        .def(self * self)
        .def(self / double())
        .def(self / self)
        .def(self + self)
        .def(self - self)
        .def(self += self)
        .def(self -= self)
        .def(self *= double())
        .def(self /= double())
        ;

    // =============================
    // Bindings for Rotation3
    // =============================
    class_<geom3::Rotation3>(py_geom3, "Rotation3")
        .def(init<>(), "Identity rotation.")
        .def(init<const geom3::UnitVector3&, double>(),
             arg("axis"), arg("angle"),
             "Construct a rotation about a given axis by an angle (in radians).")
        .def(init<const geom3::Matrix3x3&>(), arg("matrix"),
             "Construct a rotation from a 3x3 matrix.")
        .def("axis", &geom3::Rotation3::axis,
             return_value_policy::reference_internal,
             "Return the axis of rotation.")
        .def("angle", &geom3::Rotation3::angle, "Return the rotation angle (in radians).")
        // Overload operator* for rotation of Vector3, UnitVector3, and for composing rotations.
        .def("__mul__", overload_cast<const geom3::Vector3&>(&geom3::Rotation3::operator*, const_))
        .def("__mul__", overload_cast<const geom3::UnitVector3&>(&geom3::Rotation3::operator*, const_))
        .def("__mul__", overload_cast<const geom3::Rotation3&>(&geom3::Rotation3::operator*, const_))
        .def("rotate", &geom3::Rotation3::rotate, "Rotate a Vector3 using a fast rotation method.")
        .def("__imul__", &geom3::Rotation3::operator*=, "Compose another rotation (in-place).")
        .def("inverse", &geom3::Rotation3::inverse, "Return the inverse rotation.")
        .def("distance", &geom3::Rotation3::distance, "Return the angular distance to another rotation.")
        .def("matrix", &geom3::Rotation3::matrix, "Return the corresponding 3x3 matrix.")
        // Static functions for interpolation and random rotation
        .def_static("interpolate", overload_cast<double, double, const geom3::Rotation3&, const geom3::Rotation3&, double>
             (&geom3::Rotation3::interpolate),
             arg("t0"), arg("t1"), arg("r0"), arg("r1"), arg("t"),
             "Linear rotation interpolation.")
        .def_static("interpolate_cubic", overload_cast<double, double, const geom3::Rotation3&,
             const geom3::Rotation3&, const geom3::Rotation3&, const geom3::Rotation3&, double>
             (&geom3::Rotation3::interpolate),
             arg("t0"), arg("t1"), arg("r0"), arg("r1_3"), arg("r2_3"), arg("r1"), arg("t"),
             "Cubic rotation interpolation.")
        .def_static("random", &geom3::Rotation3::random, arg("rnd0"), arg("rnd1"), arg("rnd2"),
             "Generate a random rotation given three random numbers in [0, 1).")
        ;

    // =============================
    // Bind free functions
    // =============================
    py_geom3.def("deg2rad", &geom3::deg2rad, arg("deg"),
          "Convert degrees to radians.");
    py_geom3.def("rad2deg", &geom3::rad2deg, arg("rad"),
          "Convert radians to degrees.");
    py_geom3.def("angularVelocity", &geom3::angularVelocity,
          arg("r1"), arg("r2"), arg("delta_t"),
          "Compute the average angular velocity between two orientations.");

    // --- Binding for rk::P4 ---
    class_<rk::P4>(py_rk, "P4")
        .def(init<>())
        .def(init<const geom3::Vector3&, double, bool>(),
             arg("p"), arg("m"), arg("isEnergyNegative") = false)
        .def(init<double, const geom3::Vector3&>(),
             arg("e"), arg("p"))
        .def("momentum", &rk::P4::momentum, return_value_policy::reference_internal)
        .def("e", &rk::P4::e)
        .def("m", &rk::P4::m)
        .def("px", &rk::P4::px)
        .def("py", &rk::P4::py)
        .def("pz", &rk::P4::pz)
        .def("p", &rk::P4::p)
        .def("pt", &rk::P4::pt)
        .def("transverse", &rk::P4::transverse)
        .def("et", &rk::P4::et)
        .def("rapidity", &rk::P4::rapidity)
        .def("eta", &rk::P4::eta)
        .def("theta", &rk::P4::theta)
        .def("cosTheta", &rk::P4::cosTheta)
        .def("phi", &rk::P4::phi)
        .def("beta", &rk::P4::beta)
        .def("gamma", &rk::P4::gamma)
        .def("betaGamma", &rk::P4::betaGamma)
        .def("velocity", &rk::P4::velocity)
        .def("fourVelocity", &rk::P4::fourVelocity)
        .def("dot", &rk::P4::dot)
        .def("squared", &rk::P4::squared)
        .def("boost", &rk::P4::boost)
        .def("rotate", &rk::P4::rotate)
        .def("transform", &rk::P4::transform)
        .def("restBoost", &rk::P4::restBoost)
        .def("labBoost", &rk::P4::labBoost)
        .def("reverseE", &rk::P4::reverseE)
        .def("reverseP", &rk::P4::reverseP)
        // Operators
        .def(self == self)
        .def(self != self)
        .def(self + self)
        .def(self - self)
        .def(self * double())
        .def(double() * self)
        .def(self / double())
        .def(self += self)
        .def(self -= self)
        .def(self *= double())
        .def(self /= double())
        .def(-self)
        .def(+self)
        ;

    // --- Binding for rk::Boost ---
    class_<rk::Boost>(py_rk, "Boost")
        .def(init<>())
        .def(init<const geom3::Vector3&>(),
             arg("velocity_in_units_of_c"))
        .def(init<const geom3::UnitVector3&, double>(),
             arg("direction"), arg("rapidity"))
        .def("direction", &rk::Boost::direction, return_value_policy::reference_internal)
        .def("rapidity", &rk::Boost::rapidity)
        .def("beta", &rk::Boost::beta)
        .def("gamma", &rk::Boost::gamma)
        .def("betaGamma", &rk::Boost::betaGamma)
        .def("velocity", &rk::Boost::velocity)
        .def("inverse", &rk::Boost::inverse)
        .def(self == self)
        .def(self != self)
        ;

    // --- Bind operator* for Boost acting on P4 ---
    py_rk.def("boost_p4", [](const rk::Boost &b, const rk::P4 &p) {
        return b * p;
    }, arg("b"), arg("p"));

    // --- Binding for rk::Point4 ---
    class_<rk::Point4>(py_rk, "Point4")
        .def(init<>())
        .def(init<double, const geom3::Point3&>(),
             arg("t"), arg("location"))
        .def("location", &rk::Point4::location, return_value_policy::reference_internal)
        .def("t", &rk::Point4::t)
        .def("x", &rk::Point4::x)
        .def("y", &rk::Point4::y)
        .def("z", &rk::Point4::z)
        .def(self == self)
        .def(self != self)
        // If needed, you can bind the + and - operators that mix Point4 and P4.
        .def("__iadd__", &rk::Point4::operator+=)
        .def("__isub__", &rk::Point4::operator-=)
        ;

    // --- Binding for free functions ---
    py_rk.def("transformVelocity", &rk::transformVelocity,
          arg("b"), arg("v"));
    py_rk.def("invMass", (double(*)(const rk::P4&, const rk::P4&)) &rk::invMass,
          arg("p1"), arg("p2"));
    py_rk.def("invMass3", (double(*)(const rk::P4&, const rk::P4&, const rk::P4&)) &rk::invMass,
          arg("p1"), arg("p2"), arg("p3"));
    py_rk.def("lambda", &rk::lambda,
          arg("x"), arg("y"), arg("z"));
    py_rk.def("phaseSpaceDecay", (std::pair<rk::P4, rk::P4>(*)(const rk::P4&, double, double, double, double))
          &rk::phaseSpaceDecay,
          arg("parent"), arg("m1"), arg("m2"),
          arg("rnd0"), arg("rnd1"));
    // For the void version that writes into pointers, we bind a helper lambda that returns a tuple.
    py_rk.def("phaseSpaceDecay_inplace",
          [](const rk::P4 &parent, double m1, double m2, double rnd0, double rnd1) {
              rk::P4 dau1, dau2;
              rk::phaseSpaceDecay(parent, m1, m2, rnd0, rnd1, &dau1, &dau2);
              return std::make_tuple(dau1, dau2);
          },
          arg("parent"), arg("m1"), arg("m2"),
          arg("rnd0"), arg("rnd1"));
    py_rk.def("peripheralSplit",
          [](const rk::P4 &pa, const rk::P4 &pb, double m1, double m2, double rnd0, double rnd1) {
              rk::P4 dau1, dau2;
              double t0, t1, t;
              rk::peripheralSplit(pa, pb, m1, m2, rnd0, rnd1, &dau1, &dau2, &t0, &t1, &t);
              return std::make_tuple(dau1, dau2, t0, t1, t);
          },
          arg("pa"), arg("pb"), arg("m1"), arg("m2"),
          arg("rnd0"), arg("rnd1"));
}

