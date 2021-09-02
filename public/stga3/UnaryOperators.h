#ifndef LI_STGA3_UnaryOperators_H
#define LI_STGA3_UnaryOperators_H

#include <array>
#include <cmath>

namespace stga3 {

//-----------------------
// Boost unary operations
//-----------------------

template<typename T>
inline Boost<T> Boost<T>::negation() const {
    Boost<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline Boost<T> operator-(const Boost<T> &A) {
    return A.negation();
};

template<typename T>
inline Boost<T> Boost<T>::involution() const {
    Boost<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    return res;
};

template<typename T>
inline Boost<T> Boost<T>::reversion() const {
    Boost<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline Boost<T> Boost<T>::conjugate() const {
    Boost<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130MV<T> Boost<T>::dual() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=(*this)[1];
    res[9]=(*this)[2];
    res[10]=(*this)[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> Boost<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2)));
    return res;
};

template<typename T>
inline R130B0<T> Boost<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2)));
    return res;
};



//------------------------
// R130B0 unary operations
//------------------------

template<typename T>
inline R130B0<T> R130B0<T>::negation() const {
    R130B0<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> operator-(const R130B0<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B0<T> R130B0<T>::involution() const {
    R130B0<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B0<T>::reversion() const {
    R130B0<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B0<T>::conjugate() const {
    R130B0<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B0<T>::dual() const {
    R130B4<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B0<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B0<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130B0<T> exp(const R130B0<T> &A) {
    R130B0<T> res;
    res[0]=std::exp(A[0]);
    return res;
};

//------------------------
// R130B1 unary operations
//------------------------

template<typename T>
inline R130B1<T> R130B1<T>::negation() const {
    R130B1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B1<T> R130B1<T>::involution() const {
    R130B1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B1<T>::reversion() const {
    R130B1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B1<T>::conjugate() const {
    R130B1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B1<T>::dual() const {
    R130B3<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    return res;
};

template<typename T>
inline R130B0<T> R130B1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(-std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(-std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2)));
    return res;
};



//---------------------------
// R130B1Sm1 unary operations
//---------------------------

template<typename T>
inline R130B1Sm1<T> R130B1Sm1<T>::negation() const {
    R130B1Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator-(const R130B1Sm1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B1Sm1<T> R130B1Sm1<T>::involution() const {
    R130B1Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B1Sm1<T>::reversion() const {
    R130B1Sm1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B1Sm1<T>::conjugate() const {
    R130B1Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B1Sm1<T>::dual() const {
    R130B3Sm1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B0<T> R130B1Sm1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B1Sm1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline R130MV<T> exp(const R130B1Sm1<T> &A) {
    R130MV<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2) + std::pow(A[1], 2) + std::pow(A[2], 2));
    T x1 = std::sin(x0)/x0;
    res[0]=std::cos(x0);
    res[1]=0;
    res[2]=A[0]*x1;
    res[3]=A[1]*x1;
    res[4]=A[2]*x1;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//---------------------------
// R130B1Sp1 unary operations
//---------------------------

template<typename T>
inline R130B1Sp1<T> R130B1Sp1<T>::negation() const {
    R130B1Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator-(const R130B1Sp1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B1Sp1<T> R130B1Sp1<T>::involution() const {
    R130B1Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B1Sp1<T>::reversion() const {
    R130B1Sp1<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B1Sp1<T>::conjugate() const {
    R130B1Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B1Sp1<T>::dual() const {
    R130B3Sp1<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B1Sp1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(-std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B1Sp1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(-std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130MV<T> exp(const R130B1Sp1<T> &A) {
    R130MV<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2));
    res[0]=std::cosh(x0);
    res[1]=A[0]*std::sinh(x0)/x0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//------------------------
// R130B2 unary operations
//------------------------

template<typename T>
inline R130B2<T> R130B2<T>::negation() const {
    R130B2<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B2<T> R130B2<T>::involution() const {
    R130B2<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    res[4]=(*this)[4];
    res[5]=(*this)[5];
    return res;
};

template<typename T>
inline R130B2<T> R130B2<T>::reversion() const {
    R130B2<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    return res;
};

template<typename T>
inline R130B2<T> R130B2<T>::conjugate() const {
    R130B2<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    return res;
};

template<typename T>
inline R130B2<T> R130B2<T>::dual() const {
    R130B2<T> res;
    res[0]=(*this)[3];
    res[1]=(*this)[4];
    res[2]=(*this)[5];
    res[3]=(*this)[0];
    res[4]=(*this)[1];
    res[5]=(*this)[2];
    return res;
};

template<typename T>
inline R130B0<T> R130B2<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B2<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2)));
    return res;
};



//---------------------------
// R130B2Sm1 unary operations
//---------------------------

template<typename T>
inline R130B2Sm1<T> R130B2Sm1<T>::negation() const {
    R130B2Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator-(const R130B2Sm1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B2Sm1<T> R130B2Sm1<T>::involution() const {
    R130B2Sm1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B2Sm1<T>::reversion() const {
    R130B2Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B2Sm1<T>::conjugate() const {
    R130B2Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B2Sm1<T>::dual() const {
    R130B2Sp1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B0<T> R130B2Sm1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B2Sm1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline Rotation<T> exp(const R130B2Sm1<T> &A) {
    Rotation<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2) + std::pow(A[1], 2) + std::pow(A[2], 2));
    T x1 = std::sin(x0)/x0;
    res[0]=std::cos(x0);
    res[1]=A[0]*x1;
    res[2]=A[1]*x1;
    res[3]=A[2]*x1;
    return res;
};

//---------------------------
// R130B2Sp1 unary operations
//---------------------------

template<typename T>
inline R130B2Sp1<T> R130B2Sp1<T>::negation() const {
    R130B2Sp1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator-(const R130B2Sp1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B2Sp1<T> R130B2Sp1<T>::involution() const {
    R130B2Sp1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B2Sp1<T>::reversion() const {
    R130B2Sp1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B2Sp1<T>::conjugate() const {
    R130B2Sp1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B2Sp1<T>::dual() const {
    R130B2Sm1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B0<T> R130B2Sp1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B2Sp1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline Boost<T> exp(const R130B2Sp1<T> &A) {
    Boost<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2) + std::pow(A[1], 2) + std::pow(A[2], 2));
    T x1 = std::sinh(x0)/x0;
    res[0]=std::cosh(x0);
    res[1]=A[0]*x1;
    res[2]=A[1]*x1;
    res[3]=A[2]*x1;
    return res;
};

//------------------------
// R130B3 unary operations
//------------------------

template<typename T>
inline R130B3<T> R130B3<T>::negation() const {
    R130B3<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B3<T> R130B3<T>::involution() const {
    R130B3<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B3<T>::reversion() const {
    R130B3<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B3<T>::conjugate() const {
    R130B3<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B3<T>::dual() const {
    R130B1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130B0<T> R130B3<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B3<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2)));
    return res;
};



//---------------------------
// R130B3Sm1 unary operations
//---------------------------

template<typename T>
inline R130B3Sm1<T> R130B3Sm1<T>::negation() const {
    R130B3Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator-(const R130B3Sm1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B3Sm1<T> R130B3Sm1<T>::involution() const {
    R130B3Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B3Sm1<T>::reversion() const {
    R130B3Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B3Sm1<T>::conjugate() const {
    R130B3Sm1<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B3Sm1<T>::dual() const {
    R130B1Sm1<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    return res;
};

template<typename T>
inline R130B0<T> R130B3Sm1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B3Sm1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2)));
    return res;
};

template<typename T>
inline R130MV<T> exp(const R130B3Sm1<T> &A) {
    R130MV<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2) + std::pow(A[1], 2) + std::pow(A[2], 2));
    T x1 = std::sin(x0)/x0;
    res[0]=std::cos(x0);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*x1;
    res[13]=A[1]*x1;
    res[14]=A[2]*x1;
    res[15]=0;
    return res;
};

//---------------------------
// R130B3Sp1 unary operations
//---------------------------

template<typename T>
inline R130B3Sp1<T> R130B3Sp1<T>::negation() const {
    R130B3Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator-(const R130B3Sp1<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B3Sp1<T> R130B3Sp1<T>::involution() const {
    R130B3Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B3Sp1<T>::reversion() const {
    R130B3Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B3Sp1<T>::conjugate() const {
    R130B3Sp1<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B3Sp1<T>::dual() const {
    R130B1Sp1<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B3Sp1<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B3Sp1<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130MV<T> exp(const R130B3Sp1<T> &A) {
    R130MV<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2));
    res[0]=std::cosh(x0);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*std::sinh(x0)/x0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//------------------------
// R130B4 unary operations
//------------------------

template<typename T>
inline R130B4<T> R130B4<T>::negation() const {
    R130B4<T> res;
    res[0]=-(*this)[0];
    return res;
};

template<typename T>
inline R130B4<T> operator-(const R130B4<T> &A) {
    return A.negation();
};

template<typename T>
inline R130B4<T> R130B4<T>::involution() const {
    R130B4<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B4<T>::reversion() const {
    R130B4<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B4<T>::conjugate() const {
    R130B4<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B4<T>::dual() const {
    R130B0<T> res;
    res[0]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130B4<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(-std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130B4<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(-std::pow((*this)[0], 2)));
    return res;
};

template<typename T>
inline R130MV<T> exp(const R130B4<T> &A) {
    R130MV<T> res;
    T x0 = std::sqrt(std::pow(A[0], 2));
    res[0]=std::cos(x0);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*std::sin(x0)/x0;
    return res;
};

//------------------------
// R130MV unary operations
//------------------------

template<typename T>
inline R130MV<T> R130MV<T>::negation() const {
    R130MV<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    res[6]=-(*this)[6];
    res[7]=-(*this)[7];
    res[8]=-(*this)[8];
    res[9]=-(*this)[9];
    res[10]=-(*this)[10];
    res[11]=-(*this)[11];
    res[12]=-(*this)[12];
    res[13]=-(*this)[13];
    res[14]=-(*this)[14];
    res[15]=-(*this)[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A) {
    return A.negation();
};

template<typename T>
inline R130MV<T> R130MV<T>::involution() const {
    R130MV<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=(*this)[5];
    res[6]=(*this)[6];
    res[7]=(*this)[7];
    res[8]=(*this)[8];
    res[9]=(*this)[9];
    res[10]=(*this)[10];
    res[11]=-(*this)[11];
    res[12]=-(*this)[12];
    res[13]=-(*this)[13];
    res[14]=-(*this)[14];
    res[15]=(*this)[15];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::reversion() const {
    R130MV<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    res[4]=(*this)[4];
    res[5]=-(*this)[5];
    res[6]=-(*this)[6];
    res[7]=-(*this)[7];
    res[8]=-(*this)[8];
    res[9]=-(*this)[9];
    res[10]=-(*this)[10];
    res[11]=-(*this)[11];
    res[12]=-(*this)[12];
    res[13]=-(*this)[13];
    res[14]=-(*this)[14];
    res[15]=(*this)[15];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate() const {
    R130MV<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    res[6]=-(*this)[6];
    res[7]=-(*this)[7];
    res[8]=-(*this)[8];
    res[9]=-(*this)[9];
    res[10]=-(*this)[10];
    res[11]=(*this)[11];
    res[12]=(*this)[12];
    res[13]=(*this)[13];
    res[14]=(*this)[14];
    res[15]=(*this)[15];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::dual() const {
    R130MV<T> res;
    res[0]=(*this)[15];
    res[1]=-(*this)[11];
    res[2]=-(*this)[12];
    res[3]=-(*this)[13];
    res[4]=-(*this)[14];
    res[5]=(*this)[8];
    res[6]=(*this)[9];
    res[7]=(*this)[10];
    res[8]=(*this)[5];
    res[9]=(*this)[6];
    res[10]=(*this)[7];
    res[11]=(*this)[1];
    res[12]=(*this)[2];
    res[13]=(*this)[3];
    res[14]=(*this)[4];
    res[15]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> R130MV<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[10], 2) + std::pow((*this)[11], 2) - std::pow((*this)[12], 2) - std::pow((*this)[13], 2) - std::pow((*this)[14], 2) - std::pow((*this)[15], 2) - std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) - std::pow((*this)[5], 2) - std::pow((*this)[6], 2) - std::pow((*this)[7], 2) + std::pow((*this)[8], 2) + std::pow((*this)[9], 2)));
    return res;
};

template<typename T>
inline R130B0<T> R130MV<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[10], 2) + std::pow((*this)[11], 2) - std::pow((*this)[12], 2) - std::pow((*this)[13], 2) - std::pow((*this)[14], 2) - std::pow((*this)[15], 2) - std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) - std::pow((*this)[5], 2) - std::pow((*this)[6], 2) - std::pow((*this)[7], 2) + std::pow((*this)[8], 2) + std::pow((*this)[9], 2)));
    return res;
};



//--------------------------
// Rotation unary operations
//--------------------------

template<typename T>
inline Rotation<T> Rotation<T>::negation() const {
    Rotation<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const Rotation<T> &A) {
    return A.negation();
};

template<typename T>
inline Rotation<T> Rotation<T>::involution() const {
    Rotation<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    return res;
};

template<typename T>
inline Rotation<T> Rotation<T>::reversion() const {
    Rotation<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline Rotation<T> Rotation<T>::conjugate() const {
    Rotation<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    return res;
};

template<typename T>
inline R130MV<T> Rotation<T>::dual() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[1];
    res[6]=(*this)[2];
    res[7]=(*this)[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> Rotation<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2)));
    return res;
};

template<typename T>
inline R130B0<T> Rotation<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2)));
    return res;
};



//-----------------------
// Rotor unary operations
//-----------------------

template<typename T>
inline Rotor<T> Rotor<T>::negation() const {
    Rotor<T> res;
    res[0]=-(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    res[6]=-(*this)[6];
    res[7]=-(*this)[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A) {
    return A.negation();
};

template<typename T>
inline Rotor<T> Rotor<T>::involution() const {
    Rotor<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    res[4]=(*this)[4];
    res[5]=(*this)[5];
    res[6]=(*this)[6];
    res[7]=(*this)[7];
    return res;
};

template<typename T>
inline Rotor<T> Rotor<T>::reversion() const {
    Rotor<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    res[6]=-(*this)[6];
    res[7]=(*this)[7];
    return res;
};

template<typename T>
inline Rotor<T> Rotor<T>::conjugate() const {
    Rotor<T> res;
    res[0]=(*this)[0];
    res[1]=-(*this)[1];
    res[2]=-(*this)[2];
    res[3]=-(*this)[3];
    res[4]=-(*this)[4];
    res[5]=-(*this)[5];
    res[6]=-(*this)[6];
    res[7]=(*this)[7];
    return res;
};

template<typename T>
inline Rotor<T> Rotor<T>::dual() const {
    Rotor<T> res;
    res[0]=(*this)[7];
    res[1]=(*this)[4];
    res[2]=(*this)[5];
    res[3]=(*this)[6];
    res[4]=(*this)[1];
    res[5]=(*this)[2];
    res[6]=(*this)[3];
    res[7]=(*this)[0];
    return res;
};

template<typename T>
inline R130B0<T> Rotor<T>::norm() const {
    R130B0<T> res;
    res[0]=std::fabs(std::sqrt(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) - std::pow((*this)[7], 2)));
    return res;
};

template<typename T>
inline R130B0<T> Rotor<T>::invnorm() const {
    R130B0<T> res;
    res[0]=1.0/std::fabs(std::sqrt(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) - std::pow((*this)[7], 2)));
    return res;
};




} // namespace stga3

#endif // LI_STGA3_UnaryOperators_H
