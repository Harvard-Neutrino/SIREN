#ifndef LI_STGA3_Utilities_H
#define LI_STGA3_Utilities_H

#include <iostream>

#include "stga3/STGA3.h"
#include "stga3/Typedefs.h"

namespace stga3 {

template<typename T>
inline Boost<T> beta_to_boost(Beta<T> const & beta) {
    Boost<T> boost = exp(beta / 2.0);
    return boost;
}

template<typename T>
inline Boost<T> boost_from_beta(Beta<T> const & beta) {
    return exp(beta / 2.0);
}

template<typename T, template <typename> class U>
decltype(std::declval<Boost<T>>().conjugate(std::declval<U<T>>())) apply_boost(Boost<T> const & boost, U<T> const & element) {
    return boost.conjugate(element);
}

template<typename T, template <typename> class U>
decltype(std::declval<Boost<T>>().conjugate(std::declval<U<T>>())) apply_boost(Beta<T> const & beta, U<T> const & element) {
    Boost<T> boost = beta_to_boost(beta);
    return boost.conjugate(element);
}

template<typename T>
inline Plane<T> tangent_plane(ThreeVector<T> const & vec) {
    return Plane<T>{vec.e1(), vec.e2(), vec.e3()};
}

template<typename T>
inline Plane<T> unit_tangent_plane(ThreeVector<T> const & vec) {
    T const & x1 = vec.e1();
    T const & x2 = vec.e2();
    T const & x3 = vec.e3();
    T const norm = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    return Plane<T>{x1/norm, x2/norm, x3/norm};
}

template<typename T>
inline Plane<T> tangent_plane(FourVector<T> const & vec) {
    return Plane<T>{vec.e1(), vec.e2(), vec.e3()};
}

template<typename T>
inline Plane<T> unit_tangent_plane(FourVector<T> const & vec) {
    T const & x1 = vec.e1();
    T const & x2 = vec.e2();
    T const & x3 = vec.e3();
    T const norm = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    return Plane<T>{x1/norm, x2/norm, x3/norm};
}

template<typename T>
inline Rotation<T> rotation_about(Plane<T> const & plane) {
    return exp(plane / 2.0);
}

template<typename T>
inline Rotation<T> rotation_about(Plane<T> const & plane, T phi) {
    return exp(plane * (plane.invnorm() * phi / 2.0));
}

template<typename T>
inline Rotation<T> rotation_about(ThreeVector<T> const & vec, T phi) {
    Plane<T> plane = unit_tangent_plane(vec);
    return exp(plane * (phi / 2.0));
}

template<typename T>
inline Rotation<T> rotation_about(FourVector<T> const & vec, T phi) {
    Plane<T> plane = unit_tangent_plane(vec);
    return exp(plane * (phi / 2.0));
}

template<typename T, template <typename> class U>
decltype(std::declval<Rotation<T>>().conjugate(std::declval<U<T>>())) apply_rotation(Rotation<T> const & rot, U<T> const & element) {
    return rot.conjugate(element);
}

template<typename T, template <typename> class U>
decltype(std::declval<Rotation<T>>().conjugate(std::declval<U<T>>())) apply_rotation(Plane<T> const & plane, U<T> const & element) {
    return rotation_about(plane).conjugate(element);
}

template<typename T, template <typename> class U>
decltype(std::declval<Rotation<T>>().conjugate(std::declval<U<T>>())) apply_rotation(Plane<T> const & plane, T phi, U<T> const & element) {
    return rotation_about(plane, phi).conjugate(element);
}

template<typename T, template <typename> class U>
decltype(std::declval<Rotation<T>>().conjugate(std::declval<U<T>>())) apply_rotation(FourVector<T> const & vec, T phi, U<T> const & element) {
    return rotation_about(vec, phi).conjugate(element);
}

template<typename T>
inline Beta<T> beta_to_rest_frame_of(FourVector<T> const & vec) {
    R130B1Sp1<T> gamma0{1};
    return (vec ^ gamma0) / (vec | gamma0);
}

template<typename T>
inline Beta<T> beta_between(FourVector<T> const & vec0, FourVector<T> const & vec1) {
    return (vec0 ^ vec1) / (vec0 | vec1);
}

template<typename T>
inline Boost<T> boost_to_rest_frame_of(FourVector<T> const & vec) {
    return boost_from_beta(beta_to_rest_frame_of(vec));
}

//template<typename T>
//inline Boost<T> boost_between(FourVector<T> const & vec0, FourVector<T> const & vec1) {
//    return boost_from_beta(beta_between(vec0, vec1));
//}

template<typename T>
Rotation<T> rotation_between(ThreeVector<T> const & vec0, ThreeVector<T> const & vec1) {
    stga3::ThreeVector<T> dir0 = vec0 * vec0.invnorm();
    stga3::ThreeVector<T> dir1 = vec1 * vec1.invnorm();
    stga3::Rotation<T> rot = dir0 * dir1;
    rot.scalar() -= 1.0;
    rot = rot * rot.invnorm();
    return rot;
}

template<typename T>
Rotation<T> rotation_between(FourVector<T> const & vec0, FourVector<T> const & vec1) {
    T const & x1 = vec0.e1();
    T const & x2 = vec0.e2();
    T const & x3 = vec0.e3();
    T const x_norm = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    stga3::ThreeVector<T> dir0{x1 / x_norm, x2 / x_norm, x3 / x_norm};
    T const & y1 = vec1.e1();
    T const & y2 = vec1.e2();
    T const & y3 = vec1.e3();
    T const y_norm = std::sqrt(y1*y1 + y2*y2 + y3*y3);
    stga3::ThreeVector<T> dir1{y1 / y_norm, y2 / y_norm, y3 / y_norm};
    stga3::Rotation<T> rot = dir0 * dir1;
    rot.scalar() -= 1.0;
    rot = rot * rot.invnorm();
    return rot;
}

template<typename T>
Rotation<T> rotation_between(ThreeVector<T> const & vec0, FourVector<T> const & vec1) {
    stga3::ThreeVector<T> dir0 = vec0 * vec0.invnorm();
    T const & y1 = vec1.e1();
    T const & y2 = vec1.e2();
    T const & y3 = vec1.e3();
    T const y_norm = std::sqrt(y1*y1 + y2*y2 + y3*y3);
    stga3::ThreeVector<T> dir1{y1 / y_norm, y2 / y_norm, y3 / y_norm};
    stga3::Rotation<T> rot = dir0 * dir1;
    rot.scalar() -= 1.0;
    rot = rot * rot.invnorm();
    return rot;
}

template<typename T>
Rotation<T> rotation_between(FourVector<T> const & vec0, ThreeVector<T> const & vec1) {
    T const & x1 = vec0.e1();
    T const & x2 = vec0.e2();
    T const & x3 = vec0.e3();
    T const x_norm = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    stga3::ThreeVector<T> dir0{x1 / x_norm, x2 / x_norm, x3 / x_norm};
    stga3::ThreeVector<T> dir1 = vec1 * vec1.invnorm();
    stga3::Rotation<T> rot = dir0 * dir1;
    rot.scalar() -= 1.0;
    rot = rot * rot.invnorm();
    return rot;
}

} // namespace stga3

#endif // LI_STGA3_Utilities_H
