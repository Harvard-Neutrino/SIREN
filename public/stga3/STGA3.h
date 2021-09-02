#ifndef LI_STGA3_STGA3_H
#define LI_STGA3_STGA3_H

#include <array>
#include <cmath>
#include <algorithm>
#include <initializer_list>

namespace stga3 {
template <typename T>
class Boost;
template <typename T>
class R130B0;
template <typename T>
class R130B1;
template <typename T>
class R130B1Sm1;
template <typename T>
class R130B1Sp1;
template <typename T>
class R130B2;
template <typename T>
class R130B2Sm1;
template <typename T>
class R130B2Sp1;
template <typename T>
class R130B3;
template <typename T>
class R130B3Sm1;
template <typename T>
class R130B3Sp1;
template <typename T>
class R130B4;
template <typename T>
class R130MV;
template <typename T>
class Rotation;
template <typename T>
class Rotor;

template <typename T>
using Scalar = R130B0<T>;
template <typename T>
using Vector = R130B1<T>;
template <typename T>
using Bivector = R130B2<T>;
template <typename T>
using Trivector = R130B3<T>;
template <typename T>
using Pseudoscalar = R130B4<T>;
template <typename T>
using Multivector = R130MV<T>;

template<typename T>
class Boost {
private:
    std::array<T, 4> mvec;
public:
    Boost() {mvec.fill(0);}
    Boost(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & scalar() {return (*this)[0];}
    T const & scalar() const {return (*this)[0];}
    T & e10() {return (*this)[1];}
    T const & e10() const {return (*this)[1];}
    T & e20() {return (*this)[2];}
    T const & e20() const {return (*this)[2];}
    T & e30() {return (*this)[3];}
    T const & e30() const {return (*this)[3];}
    inline operator R130MV<T>() const;
    inline operator Rotor<T>() const;
    inline Boost<T> negation() const;
    inline Boost<T> involution() const;
    inline Boost<T> reversion() const;
    inline Boost<T> conjugate() const;
    inline R130MV<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline R130MV<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline R130MV<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B0 {
private:
    std::array<T, 1> mvec;
public:
    R130B0() {mvec.fill(0);}
    R130B0(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & scalar() {return (*this)[0];}
    T const & scalar() const {return (*this)[0];}
    inline operator T() const;
    inline operator Boost<T>() const;
    inline operator R130MV<T>() const;
    inline operator Rotation<T>() const;
    inline operator Rotor<T>() const;
    inline R130B0<T> negation() const;
    inline R130B0<T> involution() const;
    inline R130B0<T> reversion() const;
    inline R130B0<T> conjugate() const;
    inline R130B4<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B1 {
private:
    std::array<T, 4> mvec;
public:
    R130B1() {mvec.fill(0);}
    R130B1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e0() {return (*this)[0];}
    T const & e0() const {return (*this)[0];}
    T & e1() {return (*this)[1];}
    T const & e1() const {return (*this)[1];}
    T & e2() {return (*this)[2];}
    T const & e2() const {return (*this)[2];}
    T & e3() {return (*this)[3];}
    T const & e3() const {return (*this)[3];}
    inline operator R130MV<T>() const;
    inline R130B1<T> negation() const;
    inline R130B1<T> involution() const;
    inline R130B1<T> reversion() const;
    inline R130B1<T> conjugate() const;
    inline R130B3<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline R130MV<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline R130MV<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B1Sm1 {
private:
    std::array<T, 3> mvec;
public:
    R130B1Sm1() {mvec.fill(0);}
    R130B1Sm1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e1() {return (*this)[0];}
    T const & e1() const {return (*this)[0];}
    T & e2() {return (*this)[1];}
    T const & e2() const {return (*this)[1];}
    T & e3() {return (*this)[2];}
    T const & e3() const {return (*this)[2];}
    inline operator R130B1<T>() const;
    inline operator R130MV<T>() const;
    inline R130B1Sm1<T> negation() const;
    inline R130B1Sm1<T> involution() const;
    inline R130B1Sm1<T> reversion() const;
    inline R130B1Sm1<T> conjugate() const;
    inline R130B3Sm1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B1Sp1 {
private:
    std::array<T, 1> mvec;
public:
    R130B1Sp1() {mvec.fill(0);}
    R130B1Sp1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e0() {return (*this)[0];}
    T const & e0() const {return (*this)[0];}
    inline operator R130B1<T>() const;
    inline operator R130MV<T>() const;
    inline R130B1Sp1<T> negation() const;
    inline R130B1Sp1<T> involution() const;
    inline R130B1Sp1<T> reversion() const;
    inline R130B1Sp1<T> conjugate() const;
    inline R130B3Sp1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B2 {
private:
    std::array<T, 6> mvec;
public:
    R130B2() {mvec.fill(0);}
    R130B2(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e10() {return (*this)[0];}
    T const & e10() const {return (*this)[0];}
    T & e20() {return (*this)[1];}
    T const & e20() const {return (*this)[1];}
    T & e30() {return (*this)[2];}
    T const & e30() const {return (*this)[2];}
    T & e32() {return (*this)[3];}
    T const & e32() const {return (*this)[3];}
    T & e13() {return (*this)[4];}
    T const & e13() const {return (*this)[4];}
    T & e21() {return (*this)[5];}
    T const & e21() const {return (*this)[5];}
    inline operator R130MV<T>() const;
    inline operator Rotor<T>() const;
    inline R130B2<T> negation() const;
    inline R130B2<T> involution() const;
    inline R130B2<T> reversion() const;
    inline R130B2<T> conjugate() const;
    inline R130B2<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Rotor<T> conjugate(const Boost<T> &A) const;
    inline R130MV<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130MV<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotor<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B2Sm1 {
private:
    std::array<T, 3> mvec;
public:
    R130B2Sm1() {mvec.fill(0);}
    R130B2Sm1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e32() {return (*this)[0];}
    T const & e32() const {return (*this)[0];}
    T & e13() {return (*this)[1];}
    T const & e13() const {return (*this)[1];}
    T & e21() {return (*this)[2];}
    T const & e21() const {return (*this)[2];}
    inline operator R130B2<T>() const;
    inline operator R130MV<T>() const;
    inline operator Rotation<T>() const;
    inline operator Rotor<T>() const;
    inline R130B2Sm1<T> negation() const;
    inline R130B2Sm1<T> involution() const;
    inline R130B2Sm1<T> reversion() const;
    inline R130B2Sm1<T> conjugate() const;
    inline R130B2Sp1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B2Sp1 {
private:
    std::array<T, 3> mvec;
public:
    R130B2Sp1() {mvec.fill(0);}
    R130B2Sp1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e10() {return (*this)[0];}
    T const & e10() const {return (*this)[0];}
    T & e20() {return (*this)[1];}
    T const & e20() const {return (*this)[1];}
    T & e30() {return (*this)[2];}
    T const & e30() const {return (*this)[2];}
    inline operator Boost<T>() const;
    inline operator R130B2<T>() const;
    inline operator R130MV<T>() const;
    inline operator Rotor<T>() const;
    inline R130B2Sp1<T> negation() const;
    inline R130B2Sp1<T> involution() const;
    inline R130B2Sp1<T> reversion() const;
    inline R130B2Sp1<T> conjugate() const;
    inline R130B2Sm1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B3 {
private:
    std::array<T, 4> mvec;
public:
    R130B3() {mvec.fill(0);}
    R130B3(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e123() {return (*this)[0];}
    T const & e123() const {return (*this)[0];}
    T & e320() {return (*this)[1];}
    T const & e320() const {return (*this)[1];}
    T & e130() {return (*this)[2];}
    T const & e130() const {return (*this)[2];}
    T & e210() {return (*this)[3];}
    T const & e210() const {return (*this)[3];}
    inline operator R130MV<T>() const;
    inline R130B3<T> negation() const;
    inline R130B3<T> involution() const;
    inline R130B3<T> reversion() const;
    inline R130B3<T> conjugate() const;
    inline R130B1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline R130MV<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline R130MV<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B3Sm1 {
private:
    std::array<T, 3> mvec;
public:
    R130B3Sm1() {mvec.fill(0);}
    R130B3Sm1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e320() {return (*this)[0];}
    T const & e320() const {return (*this)[0];}
    T & e130() {return (*this)[1];}
    T const & e130() const {return (*this)[1];}
    T & e210() {return (*this)[2];}
    T const & e210() const {return (*this)[2];}
    inline operator R130B3<T>() const;
    inline operator R130MV<T>() const;
    inline R130B3Sm1<T> negation() const;
    inline R130B3Sm1<T> involution() const;
    inline R130B3Sm1<T> reversion() const;
    inline R130B3Sm1<T> conjugate() const;
    inline R130B1Sm1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B3Sp1 {
private:
    std::array<T, 1> mvec;
public:
    R130B3Sp1() {mvec.fill(0);}
    R130B3Sp1(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e123() {return (*this)[0];}
    T const & e123() const {return (*this)[0];}
    inline operator R130B3<T>() const;
    inline operator R130MV<T>() const;
    inline R130B3Sp1<T> negation() const;
    inline R130B3Sp1<T> involution() const;
    inline R130B3Sp1<T> reversion() const;
    inline R130B3Sp1<T> conjugate() const;
    inline R130B1Sp1<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130B4 {
private:
    std::array<T, 1> mvec;
public:
    R130B4() {mvec.fill(0);}
    R130B4(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & e0123() {return (*this)[0];}
    T const & e0123() const {return (*this)[0];}
    T & pseudoscalar() {return e0123();}
    T const & pseudoscalar() const {return e0123();}
    inline operator R130MV<T>() const;
    inline operator Rotor<T>() const;
    inline R130B4<T> negation() const;
    inline R130B4<T> involution() const;
    inline R130B4<T> reversion() const;
    inline R130B4<T> conjugate() const;
    inline R130B0<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class R130MV {
private:
    std::array<T, 16> mvec;
public:
    R130MV() {mvec.fill(0);}
    R130MV(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & scalar() {return (*this)[0];}
    T const & scalar() const {return (*this)[0];}
    T & e0() {return (*this)[1];}
    T const & e0() const {return (*this)[1];}
    T & e1() {return (*this)[2];}
    T const & e1() const {return (*this)[2];}
    T & e2() {return (*this)[3];}
    T const & e2() const {return (*this)[3];}
    T & e3() {return (*this)[4];}
    T const & e3() const {return (*this)[4];}
    T & e10() {return (*this)[5];}
    T const & e10() const {return (*this)[5];}
    T & e20() {return (*this)[6];}
    T const & e20() const {return (*this)[6];}
    T & e30() {return (*this)[7];}
    T const & e30() const {return (*this)[7];}
    T & e32() {return (*this)[8];}
    T const & e32() const {return (*this)[8];}
    T & e13() {return (*this)[9];}
    T const & e13() const {return (*this)[9];}
    T & e21() {return (*this)[10];}
    T const & e21() const {return (*this)[10];}
    T & e123() {return (*this)[11];}
    T const & e123() const {return (*this)[11];}
    T & e320() {return (*this)[12];}
    T const & e320() const {return (*this)[12];}
    T & e130() {return (*this)[13];}
    T const & e130() const {return (*this)[13];}
    T & e210() {return (*this)[14];}
    T const & e210() const {return (*this)[14];}
    T & e0123() {return (*this)[15];}
    T const & e0123() const {return (*this)[15];}
    T & pseudoscalar() {return e0123();}
    T const & pseudoscalar() const {return e0123();}
    inline R130MV<T> negation() const;
    inline R130MV<T> involution() const;
    inline R130MV<T> reversion() const;
    inline R130MV<T> conjugate() const;
    inline R130MV<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline R130MV<T> conjugate(const Boost<T> &A) const;
    inline R130MV<T> conjugate(const R130B0<T> &A) const;
    inline R130MV<T> conjugate(const R130B1<T> &A) const;
    inline R130MV<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130MV<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130MV<T> conjugate(const R130B2<T> &A) const;
    inline R130MV<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130MV<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130MV<T> conjugate(const R130B3<T> &A) const;
    inline R130MV<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130MV<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130MV<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline R130MV<T> conjugate(const Rotation<T> &A) const;
    inline R130MV<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class Rotation {
private:
    std::array<T, 4> mvec;
public:
    Rotation() {mvec.fill(0);}
    Rotation(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & scalar() {return (*this)[0];}
    T const & scalar() const {return (*this)[0];}
    T & e32() {return (*this)[1];}
    T const & e32() const {return (*this)[1];}
    T & e13() {return (*this)[2];}
    T const & e13() const {return (*this)[2];}
    T & e21() {return (*this)[3];}
    T const & e21() const {return (*this)[3];}
    inline operator R130MV<T>() const;
    inline operator Rotor<T>() const;
    inline Rotation<T> negation() const;
    inline Rotation<T> involution() const;
    inline Rotation<T> reversion() const;
    inline Rotation<T> conjugate() const;
    inline R130MV<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Boost<T> conjugate(const Boost<T> &A) const;
    inline R130B0<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1Sm1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1Sp1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2Sm1<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2Sp1<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3Sm1<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3Sp1<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130B4<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotation<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};

template<typename T>
class Rotor {
private:
    std::array<T, 8> mvec;
public:
    Rotor() {mvec.fill(0);}
    Rotor(std::initializer_list<T> const & v) {std::copy(v.begin(), v.end(), mvec.begin());}
    T & operator [] (size_t idx) {return mvec[idx];}
    T const & operator [] (size_t idx) const {return mvec[idx];}
    T & scalar() {return (*this)[0];}
    T const & scalar() const {return (*this)[0];}
    T & e10() {return (*this)[1];}
    T const & e10() const {return (*this)[1];}
    T & e20() {return (*this)[2];}
    T const & e20() const {return (*this)[2];}
    T & e30() {return (*this)[3];}
    T const & e30() const {return (*this)[3];}
    T & e32() {return (*this)[4];}
    T const & e32() const {return (*this)[4];}
    T & e13() {return (*this)[5];}
    T const & e13() const {return (*this)[5];}
    T & e21() {return (*this)[6];}
    T const & e21() const {return (*this)[6];}
    T & e0123() {return (*this)[7];}
    T const & e0123() const {return (*this)[7];}
    T & pseudoscalar() {return e0123();}
    T const & pseudoscalar() const {return e0123();}
    inline operator R130MV<T>() const;
    inline Rotor<T> negation() const;
    inline Rotor<T> involution() const;
    inline Rotor<T> reversion() const;
    inline Rotor<T> conjugate() const;
    inline Rotor<T> dual() const;
    inline R130B0<T> norm() const;
    inline R130B0<T> invnorm() const;
    inline Rotor<T> conjugate(const Boost<T> &A) const;
    inline R130MV<T> conjugate(const R130B0<T> &A) const;
    inline R130B1<T> conjugate(const R130B1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sm1<T> &A) const;
    inline R130B1<T> conjugate(const R130B1Sp1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sm1<T> &A) const;
    inline R130B2<T> conjugate(const R130B2Sp1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sm1<T> &A) const;
    inline R130B3<T> conjugate(const R130B3Sp1<T> &A) const;
    inline R130MV<T> conjugate(const R130B4<T> &A) const;
    inline R130MV<T> conjugate(const R130MV<T> &A) const;
    inline Rotor<T> conjugate(const Rotation<T> &A) const;
    inline Rotor<T> conjugate(const Rotor<T> &A) const;
};


} // namespace stga3

#include "stga3/UnaryOperators.h"
#include "stga3/BinaryOperators.h"

#endif // LI_STGA3_STGA3_H
