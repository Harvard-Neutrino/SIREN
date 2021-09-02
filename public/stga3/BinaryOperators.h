#ifndef LI_STGA3_BinaryOperators_H
#define LI_STGA3_BinaryOperators_H

#include <array>
#include <cmath>

namespace stga3 {

//---------------------------------
// (Boost, Boost) binary operations
//---------------------------------

template<typename T>
inline Boost<T> operator+(const Boost<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    return res;
};

template<typename T>
inline Boost<T> operator-(const Boost<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[1] + A[1]*B[0];
    res[6]=A[0]*B[2] + A[2]*B[0];
    res[7]=A[0]*B[3] + A[3]*B[0];
    res[8]=A[2]*B[3] - A[3]*B[2];
    res[9]=-A[1]*B[3] + A[3]*B[1];
    res[10]=A[1]*B[2] - A[2]*B[1];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Boost<T> operator|(const Boost<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0];
    res[2]=A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[3] + A[3]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const Boost<T> &A, const Boost<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[2]*B[3] - A[3]*B[2];
    res[1]=-A[1]*B[3] + A[3]*B[1];
    res[2]=A[1]*B[2] - A[2]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> Boost<T>::conjugate(const Boost<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[1];
    T x5 = (*this)[2]*x4;
    T x6 = (*this)[3]*A[3];
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[1];
    T x9 = (*this)[3]*A[2];
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[0]*A[3];
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x5 - x4*x6;
    res[6]=-A[1]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - x6*x7;
    res[7]=A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 - x4*x8 - x7*x9;
    res[8]=x10*x9 - x11*x7;
    res[9]=-x10*x8 + x11*x4;
    res[10]=(*this)[0]*(A[1]*x7 - A[2]*x4);
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//----------------------------------
// (Boost, R130B0) binary operations
//----------------------------------

template<typename T>
inline Boost<T> operator+(const Boost<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline Boost<T> operator-(const Boost<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline Boost<T> operator*(const Boost<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline Boost<T> operator/(const Boost<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

template<typename T>
inline Boost<T> operator|(const Boost<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Boost<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> Boost<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2));
    return res;
};

//----------------------------------
// (Boost, R130B1) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
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

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
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

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[2]=A[0]*B[1] + A[1]*B[0];
    res[3]=A[0]*B[2] + A[2]*B[0];
    res[4]=A[0]*B[3] + A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[2]*B[3] - A[3]*B[2];
    res[13]=-A[1]*B[3] + A[3]*B[1];
    res[14]=A[1]*B[2] - A[2]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[2]*B[3] - A[3]*B[2];
    res[13]=-A[1]*B[3] + A[3]*B[1];
    res[14]=A[1]*B[2] - A[2]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator^(const Boost<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> Boost<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[2];
    T x7 = (*this)[3]*A[3];
    T x8 = 2.0*(*this)[1];
    T x9 = A[1]*x8;
    T x10 = A[0]*x4;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[1]*x5 - x4*x6 - x4*x7;
    res[1]=-A[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + x6*x8 + x7*x8;
    res[2]=-(*this)[2]*x10 + 2.0*(*this)[2]*x7 + (*this)[2]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3;
    res[3]=-(*this)[3]*x10 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//-------------------------------------
// (Boost, R130B1Sm1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
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

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
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

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[2]=A[0]*B[0];
    res[3]=A[0]*B[1];
    res[4]=A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[2]*B[2] - A[3]*B[1];
    res[13]=-A[1]*B[2] + A[3]*B[0];
    res[14]=A[1]*B[1] - A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[0]*B[1];
    res[4]=A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[2]*B[2] - A[3]*B[1];
    res[13]=-A[1]*B[2] + A[3]*B[0];
    res[14]=A[1]*B[1] - A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator^(const Boost<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B1<T> Boost<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*A[0];
    T x2 = (*this)[2]*A[1];
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=-x0*(x1 + x2 + x3);
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x6 - A[0]*x7 + x2*x8 + x3*x8;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + x1*x10 + x10*x2;
    return res;
};

//-------------------------------------
// (Boost, R130B1Sp1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
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

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
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

template<typename T>
inline R130B1<T> operator*(const Boost<T> &A, const R130B1Sp1<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator|(const Boost<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const Boost<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[1]*B[0];
    res[1]=A[2]*B[0];
    res[2]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> Boost<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0]*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    res[1]=-(*this)[1]*x0;
    res[2]=-(*this)[2]*x0;
    res[3]=-(*this)[3]*x0;
    return res;
};

//----------------------------------
// (Boost, R130B2) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1] + B[0];
    res[6]=A[2] + B[1];
    res[7]=A[3] + B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1] - B[0];
    res[6]=A[2] - B[1];
    res[7]=A[3] - B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Boost<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0] - A[2]*B[5] + A[3]*B[4];
    res[2]=A[0]*B[1] + A[1]*B[5] - A[3]*B[3];
    res[3]=A[0]*B[2] - A[1]*B[4] + A[2]*B[3];
    res[4]=A[0]*B[3] + A[2]*B[2] - A[3]*B[1];
    res[5]=A[0]*B[4] - A[1]*B[2] + A[3]*B[0];
    res[6]=A[0]*B[5] + A[1]*B[1] - A[2]*B[0];
    res[7]=A[1]*B[3] + A[2]*B[4] + A[3]*B[5];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Boost<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=A[0]*B[4];
    res[6]=A[0]*B[5];
    res[7]=A[1]*B[3] + A[2]*B[4] + A[3]*B[5];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Boost<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[5] + A[3]*B[4];
    res[1]=A[1]*B[5] - A[3]*B[3];
    res[2]=-A[1]*B[4] + A[2]*B[3];
    res[3]=A[2]*B[2] - A[3]*B[1];
    res[4]=-A[1]*B[2] + A[3]*B[0];
    res[5]=A[1]*B[1] - A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> Boost<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[1], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[5]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = 2.0*(*this)[1];
    T x8 = (*this)[2]*x7;
    T x9 = (*this)[3]*A[2];
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[1]*x4;
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x7;
    T x14 = (*this)[3]*x10;
    res[0]=(*this)[2]*x5 + A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[1]*x8 - A[4]*x6 - x7*x9;
    res[1]=-(*this)[1]*x5 - A[0]*x8 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 + A[3]*x6 - x10*x9;
    res[2]=-A[0]*x13 - A[1]*x14 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[3]*x12 + A[4]*x11;
    res[3]=A[1]*x6 - A[2]*x12 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 - A[4]*x8 - A[5]*x13;
    res[4]=-A[0]*x6 + A[2]*x11 - A[3]*x8 + A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3 - A[5]*x14;
    res[5]=A[0]*x12 - A[1]*x11 - A[3]*x13 - A[4]*x14 + A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x3;
    return res;
};

//-------------------------------------
// (Boost, R130B2Sm1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[2]*B[2] + A[3]*B[1];
    res[6]=A[1]*B[2] - A[3]*B[0];
    res[7]=-A[1]*B[1] + A[2]*B[0];
    res[8]=A[0]*B[0];
    res[9]=A[0]*B[1];
    res[10]=A[0]*B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[0];
    res[9]=A[0]*B[1];
    res[10]=A[0]*B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const Boost<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> Boost<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[2]*x0;
    T x2 = (*this)[3]*x0;
    T x3 = (*this)[1]*A[1];
    T x4 = (*this)[2]*A[0];
    T x5 = std::pow((*this)[0], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = std::pow((*this)[1], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*A[2];
    res[0]=(*this)[2]*x1 - A[1]*x2;
    res[1]=-(*this)[1]*x1 + A[0]*x2;
    res[2]=x0*(x3 - x4);
    res[3]=A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - x10*x11 - x3*x9;
    res[4]=A[1]*x5 - A[1]*x6 + A[1]*x7 + A[1]*x8 - x10*x4 - x11*x9;
    res[5]=-(*this)[3]*A[0]*x10 - (*this)[3]*A[1]*x9 + A[2]*x5 + A[2]*x6 - A[2]*x7 + A[2]*x8;
    return res;
};

//-------------------------------------
// (Boost, R130B2Sp1) binary operations
//-------------------------------------

template<typename T>
inline Boost<T> operator+(const Boost<T> &A, const R130B2Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    return res;
};

template<typename T>
inline Boost<T> operator-(const Boost<T> &A, const R130B2Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0];
    res[6]=A[0]*B[1];
    res[7]=A[0]*B[2];
    res[8]=A[2]*B[2] - A[3]*B[1];
    res[9]=-A[1]*B[2] + A[3]*B[0];
    res[10]=A[1]*B[1] - A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Boost<T> operator|(const Boost<T> &A, const R130B2Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const Boost<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[2]*B[2] - A[3]*B[1];
    res[1]=-A[1]*B[2] + A[3]*B[0];
    res[2]=A[1]*B[1] - A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> Boost<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[1], 2);
    T x4 = 2.0*(*this)[1];
    T x5 = (*this)[2]*x4;
    T x6 = (*this)[3]*A[2];
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[0];
    T x9 = (*this)[3]*A[1];
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[0]*A[2];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[1]*x5 - x4*x6;
    res[1]=-A[0]*x5 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - x6*x7;
    res[2]=A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - x4*x8 - x7*x9;
    res[3]=x10*x9 - x11*x7;
    res[4]=-x10*x8 + x11*x4;
    res[5]=(*this)[0]*(A[0]*x7 - A[1]*x4);
    return res;
};

//----------------------------------
// (Boost, R130B3) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[2]*B[3] + A[3]*B[2];
    res[3]=A[1]*B[3] - A[3]*B[1];
    res[4]=-A[1]*B[2] + A[2]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=A[0]*B[1] - A[1]*B[0];
    res[13]=A[0]*B[2] - A[2]*B[0];
    res[14]=A[0]*B[3] - A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[2]*B[3] + A[3]*B[2];
    res[3]=A[1]*B[3] - A[3]*B[1];
    res[4]=-A[1]*B[2] + A[2]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator^(const Boost<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    res[3]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> Boost<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[2];
    T x7 = (*this)[3]*A[3];
    T x8 = 2.0*(*this)[1];
    T x9 = A[0]*x4;
    T x10 = A[1]*x8;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[1]*x5 + x4*x6 + x4*x7;
    res[1]=A[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + x6*x8 + x7*x8;
    res[2]=(*this)[2]*x10 + 2.0*(*this)[2]*x7 + (*this)[2]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3;
    res[3]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//-------------------------------------
// (Boost, R130B3Sm1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[2]*B[2] + A[3]*B[1];
    res[3]=A[1]*B[2] - A[3]*B[0];
    res[4]=-A[1]*B[1] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[2]*B[2] + A[3]*B[1];
    res[3]=A[1]*B[2] - A[3]*B[0];
    res[4]=-A[1]*B[1] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator^(const Boost<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    return res;
};

template<typename T>
inline R130B3<T> Boost<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*A[0];
    T x2 = (*this)[2]*A[1];
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=x0*(x1 + x2 + x3);
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x6 - A[0]*x7 + x2*x8 + x3*x8;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + x1*x10 + x10*x2;
    return res;
};

//-------------------------------------
// (Boost, R130B3Sp1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator*(const Boost<T> &A, const R130B3Sp1<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    res[3]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator|(const Boost<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const Boost<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[1]*B[0];
    res[1]=-A[2]*B[0];
    res[2]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> Boost<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0]*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    res[1]=(*this)[1]*x0;
    res[2]=(*this)[2]*x0;
    res[3]=(*this)[3]*x0;
    return res;
};

//----------------------------------
// (Boost, R130B4) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1]*B[0];
    res[9]=A[2]*B[0];
    res[10]=A[3]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1]*B[0];
    res[9]=A[2]*B[0];
    res[10]=A[3]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Boost<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> Boost<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2));
    return res;
};

//----------------------------------
// (Boost, R130MV) binary operations
//----------------------------------

template<typename T>
inline Boost<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=(*this)[0];
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=A[1] + B[5];
    res[6]=A[2] + B[6];
    res[7]=A[3] + B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=A[1] - B[5];
    res[6]=A[2] - B[6];
    res[7]=A[3] - B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Boost<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7];
    res[1]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3] + A[3]*B[4];
    res[2]=A[0]*B[2] + A[1]*B[1] - A[2]*B[14] + A[3]*B[13];
    res[3]=A[0]*B[3] + A[1]*B[14] + A[2]*B[1] - A[3]*B[12];
    res[4]=A[0]*B[4] - A[1]*B[13] + A[2]*B[12] + A[3]*B[1];
    res[5]=A[0]*B[5] + A[1]*B[0] - A[2]*B[10] + A[3]*B[9];
    res[6]=A[0]*B[6] + A[1]*B[10] + A[2]*B[0] - A[3]*B[8];
    res[7]=A[0]*B[7] - A[1]*B[9] + A[2]*B[8] + A[3]*B[0];
    res[8]=A[0]*B[8] + A[1]*B[15] + A[2]*B[7] - A[3]*B[6];
    res[9]=A[0]*B[9] - A[1]*B[7] + A[2]*B[15] + A[3]*B[5];
    res[10]=A[0]*B[10] + A[1]*B[6] - A[2]*B[5] + A[3]*B[15];
    res[11]=A[0]*B[11] - A[1]*B[12] - A[2]*B[13] - A[3]*B[14];
    res[12]=A[0]*B[12] - A[1]*B[11] + A[2]*B[4] - A[3]*B[3];
    res[13]=A[0]*B[13] - A[1]*B[4] - A[2]*B[11] + A[3]*B[2];
    res[14]=A[0]*B[14] + A[1]*B[3] - A[2]*B[2] - A[3]*B[11];
    res[15]=A[0]*B[15] + A[1]*B[8] + A[2]*B[9] + A[3]*B[10];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Boost<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[14] + A[3]*B[13];
    res[3]=A[0]*B[3] + A[1]*B[14] - A[3]*B[12];
    res[4]=A[0]*B[4] - A[1]*B[13] + A[2]*B[12];
    res[5]=A[0]*B[5] + A[1]*B[0];
    res[6]=A[0]*B[6] + A[2]*B[0];
    res[7]=A[0]*B[7] + A[3]*B[0];
    res[8]=A[0]*B[8] + A[1]*B[15];
    res[9]=A[0]*B[9] + A[2]*B[15];
    res[10]=A[0]*B[10] + A[3]*B[15];
    res[11]=A[0]*B[11];
    res[12]=A[0]*B[12] + A[2]*B[4] - A[3]*B[3];
    res[13]=A[0]*B[13] - A[1]*B[4] + A[3]*B[2];
    res[14]=A[0]*B[14] + A[1]*B[3] - A[2]*B[2];
    res[15]=A[0]*B[15] + A[1]*B[8] + A[2]*B[9] + A[3]*B[10];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Boost<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[1]*B[2] + A[2]*B[3] + A[3]*B[4];
    res[2]=A[1]*B[1];
    res[3]=A[2]*B[1];
    res[4]=A[3]*B[1];
    res[5]=-A[2]*B[10] + A[3]*B[9];
    res[6]=A[1]*B[10] - A[3]*B[8];
    res[7]=-A[1]*B[9] + A[2]*B[8];
    res[8]=A[2]*B[7] - A[3]*B[6];
    res[9]=-A[1]*B[7] + A[3]*B[5];
    res[10]=A[1]*B[6] - A[2]*B[5];
    res[11]=-A[1]*B[12] - A[2]*B[13] - A[3]*B[14];
    res[12]=-A[1]*B[11];
    res[13]=-A[2]*B[11];
    res[14]=-A[3]*B[11];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> Boost<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[3];
    T x7 = (*this)[3]*A[4];
    T x8 = 2.0*(*this)[1];
    T x9 = A[2]*x8;
    T x10 = 2.0*(*this)[2];
    T x11 = A[1]*x4;
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[2]*x8;
    T x15 = (*this)[3]*A[7];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*x10;
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x5 - x4*x6 - x4*x7;
    res[2]=-A[1]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + x6*x8 + x7*x8;
    res[3]=-(*this)[2]*x11 + (*this)[2]*x9 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 + x10*x7;
    res[4]=-(*this)[3]*x11 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3;
    res[5]=A[10]*x12 + A[5]*x0 - A[5]*x1 + A[5]*x2 + A[5]*x3 - A[6]*x14 - A[9]*x13 - x15*x8;
    res[6]=-A[10]*x5 - A[5]*x14 + A[6]*x0 + A[6]*x1 - A[6]*x2 + A[6]*x3 + A[8]*x13 - x10*x15;
    res[7]=-A[5]*x16 - A[6]*x17 + A[7]*x0 + A[7]*x1 + A[7]*x2 - A[7]*x3 - A[8]*x12 + A[9]*x5;
    res[8]=-A[10]*x16 + A[6]*x13 - A[7]*x12 + A[8]*x0 - A[8]*x1 + A[8]*x2 + A[8]*x3 - A[9]*x14;
    res[9]=-A[10]*x17 - A[5]*x13 + A[7]*x5 - A[8]*x14 + A[9]*x0 + A[9]*x1 - A[9]*x2 + A[9]*x3;
    res[10]=A[10]*x0 + A[10]*x1 + A[10]*x2 - A[10]*x3 + A[5]*x12 - A[6]*x5 - A[8]*x16 - A[9]*x17;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[12]*x5 + A[13]*x12 + A[14]*x13;
    res[12]=A[11]*x5 + A[12]*x0 + A[12]*x1 - A[12]*x2 - A[12]*x3 + A[13]*x14 + A[14]*x16;
    res[13]=A[11]*x12 + A[12]*x14 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 + A[14]*x17;
    res[14]=A[11]*x13 + A[12]*x16 + A[13]*x17 + A[14]*x0 - A[14]*x1 - A[14]*x2 + A[14]*x3;
    res[15]=A[15]*(x0 - x1 - x2 - x3);
    return res;
};

//------------------------------------
// (Boost, Rotation) binary operations
//------------------------------------

template<typename T>
inline R130MV<T> operator+(const Boost<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Boost<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Boost<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0] - A[2]*B[3] + A[3]*B[2];
    res[2]=A[1]*B[3] + A[2]*B[0] - A[3]*B[1];
    res[3]=-A[1]*B[2] + A[2]*B[1] + A[3]*B[0];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    res[6]=A[0]*B[3];
    res[7]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Boost<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    res[6]=A[0]*B[3];
    res[7]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const Boost<T> &A, const Rotation<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> Boost<T>::conjugate(const Rotation<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[3]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = (*this)[1]*A[2];
    T x8 = (*this)[2]*A[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*A[3];
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[2]*x5 - A[2]*x6;
    res[6]=-(*this)[1]*x5 + A[1]*x6;
    res[7]=x4*(x7 - x8);
    res[8]=A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - x10*x11 - x7*x9;
    res[9]=A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - x10*x8 - x11*x9;
    res[10]=-(*this)[3]*A[1]*x10 - (*this)[3]*A[2]*x9 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//---------------------------------
// (Boost, Rotor) binary operations
//---------------------------------

template<typename T>
inline Boost<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=(*this)[3];
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator+(const Boost<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Boost<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Boost<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0] - A[2]*B[6] + A[3]*B[5];
    res[2]=A[0]*B[2] + A[1]*B[6] + A[2]*B[0] - A[3]*B[4];
    res[3]=A[0]*B[3] - A[1]*B[5] + A[2]*B[4] + A[3]*B[0];
    res[4]=A[0]*B[4] + A[1]*B[7] + A[2]*B[3] - A[3]*B[2];
    res[5]=A[0]*B[5] - A[1]*B[3] + A[2]*B[7] + A[3]*B[1];
    res[6]=A[0]*B[6] + A[1]*B[2] - A[2]*B[1] + A[3]*B[7];
    res[7]=A[0]*B[7] + A[1]*B[4] + A[2]*B[5] + A[3]*B[6];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Boost<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0];
    res[2]=A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[3] + A[3]*B[0];
    res[4]=A[0]*B[4] + A[1]*B[7];
    res[5]=A[0]*B[5] + A[2]*B[7];
    res[6]=A[0]*B[6] + A[3]*B[7];
    res[7]=A[0]*B[7] + A[1]*B[4] + A[2]*B[5] + A[3]*B[6];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Boost<T> &A, const Rotor<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[6] + A[3]*B[5];
    res[1]=A[1]*B[6] - A[3]*B[4];
    res[2]=-A[1]*B[5] + A[2]*B[4];
    res[3]=A[2]*B[3] - A[3]*B[2];
    res[4]=-A[1]*B[3] + A[3]*B[1];
    res[5]=A[1]*B[2] - A[2]*B[1];
    return res;
};

template<typename T>
inline Rotor<T> Boost<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[6]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = 2.0*(*this)[1];
    T x8 = (*this)[2]*x7;
    T x9 = (*this)[3]*A[3];
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[1]*x4;
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x7;
    T x14 = (*this)[3]*x10;
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=(*this)[2]*x5 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x8 - A[5]*x6 - x7*x9;
    res[2]=-(*this)[1]*x5 - A[1]*x8 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 + A[4]*x6 - x10*x9;
    res[3]=-A[1]*x13 - A[2]*x14 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 - A[4]*x12 + A[5]*x11;
    res[4]=A[2]*x6 - A[3]*x12 + A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3 - A[5]*x8 - A[6]*x13;
    res[5]=-A[1]*x6 + A[3]*x11 - A[4]*x8 + A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x3 - A[6]*x14;
    res[6]=A[1]*x12 - A[2]*x11 - A[4]*x13 - A[5]*x14 + A[6]*x0 + A[6]*x1 + A[6]*x2 - A[6]*x3;
    res[7]=A[7]*(x0 - x1 - x2 - x3);
    return res;
};

//----------------------------------
// (Boost, scalar) binary operations
//----------------------------------

template<typename T>
inline Boost<T> operator*(const Boost<T> &A, const T &B) {
    Boost<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    return res;
};
template<typename T>
inline Boost<T> operator*(const T &A, const Boost<T> &B) {
    return B*A;
};
template<typename T>
inline Boost<T> operator/(const Boost<T> &A, const T &B) {
    return A * (1.0 / B);
};

//----------------------------------
// (R130B0, Boost) binary operations
//----------------------------------

template<typename T>
inline R130B0<T>::operator Boost<T>() const {
    Boost<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    return res;
};

template<typename T>
inline Boost<T> operator+(const R130B0<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    return res;
};

template<typename T>
inline Boost<T> operator-(const R130B0<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    return res;
};

template<typename T>
inline Boost<T> operator*(const R130B0<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline Boost<T> operator|(const R130B0<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Boost<T> R130B0<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//-----------------------------------
// (R130B0, R130B0) binary operations
//-----------------------------------

template<typename T>
inline R130B0<T> operator+(const R130B0<T> &A, const R130B0<T> &B) {
    R130B0<T> res;
    res[0]=A[0] + B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator-(const R130B0<T> &A, const R130B0<T> &B) {
    R130B0<T> res;
    res[0]=A[0] - B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator*(const R130B0<T> &A, const R130B0<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator/(const R130B0<T> &A, const R130B0<T> &B) {
    R130B0<T> res;
    res[0]=A[0]/B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B0<T> &A, const R130B0<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B0<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B0, R130B1) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
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

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
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

template<typename T>
inline R130B1<T> operator*(const R130B0<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B0<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> R130B0<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//--------------------------------------
// (R130B0, R130B1Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
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

template<typename T>
inline R130B1Sm1<T> operator*(const R130B0<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B0<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B0<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B0, R130B1Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
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

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
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

template<typename T>
inline R130B1Sp1<T> operator*(const R130B0<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator|(const R130B0<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B0<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B0, R130B2) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2<T> operator*(const R130B0<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[0]*B[4];
    res[5]=A[0]*B[5];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B0<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[0]*B[4];
    res[5]=A[0]*B[5];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2<T> R130B0<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    return res;
};

//--------------------------------------
// (R130B0, R130B2Sm1) binary operations
//--------------------------------------

template<typename T>
inline Rotation<T> operator+(const R130B0<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const R130B0<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B0<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B0<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B0<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B0, R130B2Sp1) binary operations
//--------------------------------------

template<typename T>
inline Boost<T> operator+(const R130B0<T> &A, const R130B2Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    return res;
};

template<typename T>
inline Boost<T> operator-(const R130B0<T> &A, const R130B2Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B0<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator|(const R130B0<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B0<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//-----------------------------------
// (R130B0, R130B3) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator*(const R130B0<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B0<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> R130B0<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//--------------------------------------
// (R130B0, R130B3Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B0<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B0<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B0<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B0, R130B3Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator*(const R130B0<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator|(const R130B0<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B0<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B0, R130B4) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
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
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator*(const R130B0<T> &A, const R130B4<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator|(const R130B0<T> &A, const R130B4<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> R130B0<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B0, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B0<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=(*this)[0];
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B0<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B0<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B0<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[0]*B[4];
    res[5]=A[0]*B[5];
    res[6]=A[0]*B[6];
    res[7]=A[0]*B[7];
    res[8]=A[0]*B[8];
    res[9]=A[0]*B[9];
    res[10]=A[0]*B[10];
    res[11]=A[0]*B[11];
    res[12]=A[0]*B[12];
    res[13]=A[0]*B[13];
    res[14]=A[0]*B[14];
    res[15]=A[0]*B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B0<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[0]*B[4];
    res[5]=A[0]*B[5];
    res[6]=A[0]*B[6];
    res[7]=A[0]*B[7];
    res[8]=A[0]*B[8];
    res[9]=A[0]*B[9];
    res[10]=A[0]*B[10];
    res[11]=A[0]*B[11];
    res[12]=A[0]*B[12];
    res[13]=A[0]*B[13];
    res[14]=A[0]*B[14];
    res[15]=A[0]*B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B0<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    res[6]=A[6]*x0;
    res[7]=A[7]*x0;
    res[8]=A[8]*x0;
    res[9]=A[9]*x0;
    res[10]=A[10]*x0;
    res[11]=A[11]*x0;
    res[12]=A[12]*x0;
    res[13]=A[13]*x0;
    res[14]=A[14]*x0;
    res[15]=A[15]*x0;
    return res;
};

//-------------------------------------
// (R130B0, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130B0<T>::operator Rotation<T>() const {
    Rotation<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    return res;
};

template<typename T>
inline Rotation<T> operator+(const R130B0<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const R130B0<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const R130B0<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator|(const R130B0<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Rotation<T> R130B0<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//----------------------------------
// (R130B0, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130B0<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator+(const R130B0<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const R130B0<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B0<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[0]*B[4];
    res[5]=A[0]*B[5];
    res[6]=A[0]*B[6];
    res[7]=A[0]*B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B0<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[0]*B[4];
    res[5]=A[0]*B[5];
    res[6]=A[0]*B[6];
    res[7]=A[0]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B0<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B0<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    res[6]=A[6]*x0;
    res[7]=A[7]*x0;
    return res;
};

//-----------------------------------
// (R130B0, scalar) binary operations
//-----------------------------------

template<typename T>
inline R130B0<T> operator*(const R130B0<T> &A, const T &B) {
    R130B0<T> res;
    res[0]=A[0]*B;
    return res;
};
template<typename T>
inline R130B0<T> operator*(const T &A, const R130B0<T> &B) {
    return B*A;
};
template<typename T>
inline R130B0<T> operator/(const R130B0<T> &A, const T &B) {
    return A * (1.0 / B);
};

//----------------------------------
// (R130B1, Boost) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
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

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
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

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    res[3]=-A[0]*B[2] + A[2]*B[0];
    res[4]=-A[0]*B[3] + A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[2]*B[3] + A[3]*B[2];
    res[13]=A[1]*B[3] - A[3]*B[1];
    res[14]=-A[1]*B[2] + A[2]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[2]*B[3] + A[3]*B[2];
    res[13]=A[1]*B[3] - A[3]*B[1];
    res[14]=-A[1]*B[2] + A[2]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B1<T> &A, const Boost<T> &B) {
    R130B1<T> res;
    res[0]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    res[3]=-A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> R130B1<T>::conjugate(const Boost<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = 2.0*(*this)[1];
    T x5 = (*this)[2]*x4;
    T x6 = (*this)[3]*A[3];
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[1];
    T x9 = (*this)[3]*A[2];
    T x10 = (*this)[0]*A[3];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*(x0 + x1 + x2 - x3);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x5 - x4*x6;
    res[6]=-A[1]*x5 + A[2]*x0 - A[2]*x1 + A[2]*x2 + A[2]*x3 - x6*x7;
    res[7]=A[3]*x0 + A[3]*x1 - A[3]*x2 + A[3]*x3 - x4*x8 - x7*x9;
    res[8]=x10*x7 - x11*x9;
    res[9]=-x10*x4 + x11*x8;
    res[10]=(*this)[0]*(-A[1]*x7 + A[2]*x4);
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//-----------------------------------
// (R130B1, R130B0) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
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

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
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

template<typename T>
inline R130B1<T> operator*(const R130B1<T> &A, const R130B0<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator/(const R130B1<T> &A, const R130B0<T> &B) {
    R130B1<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B1<T> &A, const R130B0<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=A[0]*(-std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    return res;
};

//-----------------------------------
// (R130B1, R130B1) binary operations
//-----------------------------------

template<typename T>
inline R130B1<T> operator+(const R130B1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[0]*B[1] + A[1]*B[0];
    res[6]=-A[0]*B[2] + A[2]*B[0];
    res[7]=-A[0]*B[3] + A[3]*B[0];
    res[8]=-A[2]*B[3] + A[3]*B[2];
    res[9]=A[1]*B[3] - A[3]*B[1];
    res[10]=-A[1]*B[2] + A[2]*B[1];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1<T> &A, const R130B1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B1<T> &A, const R130B1<T> &B) {
    R130B2<T> res;
    res[0]=-A[0]*B[1] + A[1]*B[0];
    res[1]=-A[0]*B[2] + A[2]*B[0];
    res[2]=-A[0]*B[3] + A[3]*B[0];
    res[3]=-A[2]*B[3] + A[3]*B[2];
    res[4]=A[1]*B[3] - A[3]*B[1];
    res[5]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130B1<T> R130B1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[2];
    T x7 = (*this)[3]*A[3];
    T x8 = 2.0*(*this)[1];
    T x9 = A[1]*x8;
    T x10 = A[0]*x4;
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x5 + x4*x6 + x4*x7;
    res[1]=-A[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + x6*x8 + x7*x8;
    res[2]=-(*this)[2]*x10 + 2.0*(*this)[2]*x7 + (*this)[2]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3;
    res[3]=-(*this)[3]*x10 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//--------------------------------------
// (R130B1, R130B1Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130B1<T> operator+(const R130B1<T> &A, const R130B1Sm1<T> &B) {
    R130B1<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1<T> &A, const R130B1Sm1<T> &B) {
    R130B1<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[0]*B[0];
    res[6]=-A[0]*B[1];
    res[7]=-A[0]*B[2];
    res[8]=-A[2]*B[2] + A[3]*B[1];
    res[9]=A[1]*B[2] - A[3]*B[0];
    res[10]=-A[1]*B[1] + A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1<T> &A, const R130B1Sm1<T> &B) {
    R130B0<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B1<T> &A, const R130B1Sm1<T> &B) {
    R130B2<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    res[3]=-A[2]*B[2] + A[3]*B[1];
    res[4]=A[1]*B[2] - A[3]*B[0];
    res[5]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*A[0];
    T x2 = (*this)[2]*A[1];
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=x0*(x1 + x2 + x3);
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x6 - A[0]*x7 + x2*x8 + x3*x8;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + x1*x10 + x10*x2;
    return res;
};

//--------------------------------------
// (R130B1, R130B1Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130B1<T> operator+(const R130B1<T> &A, const R130B1Sp1<T> &B) {
    R130B1<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1<T> &A, const R130B1Sp1<T> &B) {
    R130B1<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline Boost<T> operator*(const R130B1<T> &A, const R130B1Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1<T> &A, const R130B1Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B1<T> &A, const R130B1Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[1]*B[0];
    res[1]=A[2]*B[0];
    res[2]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0]*A[0];
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    res[1]=-(*this)[1]*x0;
    res[2]=-(*this)[2]*x0;
    res[3]=-(*this)[3]*x0;
    return res;
};

//-----------------------------------
// (R130B1, R130B2) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[2]=-A[0]*B[0] - A[2]*B[5] + A[3]*B[4];
    res[3]=-A[0]*B[1] + A[1]*B[5] - A[3]*B[3];
    res[4]=-A[0]*B[2] - A[1]*B[4] + A[2]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[3] - A[2]*B[4] - A[3]*B[5];
    res[12]=A[0]*B[3] - A[2]*B[2] + A[3]*B[1];
    res[13]=A[0]*B[4] + A[1]*B[2] - A[3]*B[0];
    res[14]=A[0]*B[5] - A[1]*B[1] + A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B1<T> &A, const R130B2<T> &B) {
    R130B3<T> res;
    res[0]=-A[1]*B[3] - A[2]*B[4] - A[3]*B[5];
    res[1]=A[0]*B[3] - A[2]*B[2] + A[3]*B[1];
    res[2]=A[0]*B[4] + A[1]*B[2] - A[3]*B[0];
    res[3]=A[0]*B[5] - A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B1<T> &A, const R130B2<T> &B) {
    R130B1<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=-A[0]*B[0] - A[2]*B[5] + A[3]*B[4];
    res[2]=-A[0]*B[1] + A[1]*B[5] - A[3]*B[3];
    res[3]=-A[0]*B[2] - A[1]*B[4] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> R130B1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[1], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[5]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = 2.0*(*this)[1];
    T x8 = (*this)[2]*x7;
    T x9 = (*this)[3]*A[2];
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[1]*x4;
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x7;
    T x14 = (*this)[3]*x10;
    res[0]=(*this)[2]*x5 + A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[1]*x8 - A[4]*x6 - x7*x9;
    res[1]=-(*this)[1]*x5 - A[0]*x8 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 + A[3]*x6 - x10*x9;
    res[2]=-A[0]*x13 - A[1]*x14 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[3]*x12 + A[4]*x11;
    res[3]=-A[1]*x6 + A[2]*x12 - A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + A[4]*x8 + A[5]*x13;
    res[4]=A[0]*x6 - A[2]*x11 + A[3]*x8 - A[4]*x0 + A[4]*x1 - A[4]*x2 - A[4]*x3 + A[5]*x14;
    res[5]=-A[0]*x12 + A[1]*x11 + A[3]*x13 + A[4]*x14 - A[5]*x0 - A[5]*x1 + A[5]*x2 - A[5]*x3;
    return res;
};

//--------------------------------------
// (R130B1, R130B2Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[2]*B[2] + A[3]*B[1];
    res[3]=A[1]*B[2] - A[3]*B[0];
    res[4]=-A[1]*B[1] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B1<T> &A, const R130B2Sm1<T> &B) {
    R130B3<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1<T> &A, const R130B2Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> R130B1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[2]*x0;
    T x2 = (*this)[3]*x0;
    T x3 = (*this)[1]*A[1];
    T x4 = (*this)[2]*A[0];
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[0], 2);
    T x7 = std::pow((*this)[2], 2);
    T x8 = std::pow((*this)[3], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*A[2];
    res[0]=(*this)[2]*x1 - A[1]*x2;
    res[1]=-(*this)[1]*x1 + A[0]*x2;
    res[2]=x0*(x3 - x4);
    res[3]=A[0]*x5 - A[0]*x6 - A[0]*x7 - A[0]*x8 + x10*x11 + x3*x9;
    res[4]=-A[1]*x5 - A[1]*x6 + A[1]*x7 - A[1]*x8 + x10*x4 + x11*x9;
    res[5]=(*this)[3]*A[0]*x10 + (*this)[3]*A[1]*x9 - A[2]*x5 - A[2]*x6 - A[2]*x7 + A[2]*x8;
    return res;
};

//--------------------------------------
// (R130B1, R130B2Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
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

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[2]=-A[0]*B[0];
    res[3]=-A[0]*B[1];
    res[4]=-A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[2]*B[2] + A[3]*B[1];
    res[13]=A[1]*B[2] - A[3]*B[0];
    res[14]=-A[1]*B[1] + A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B1<T> &A, const R130B2Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B1<T> &A, const R130B2Sp1<T> &B) {
    R130B1<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=-A[0]*B[0];
    res[2]=-A[0]*B[1];
    res[3]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> R130B1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[1], 2);
    T x4 = 2.0*(*this)[1];
    T x5 = (*this)[2]*x4;
    T x6 = (*this)[3]*A[2];
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[0];
    T x9 = (*this)[3]*A[1];
    T x10 = (*this)[0]*A[2];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[1]*x5 - x4*x6;
    res[1]=-A[0]*x5 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - x6*x7;
    res[2]=A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - x4*x8 - x7*x9;
    res[3]=x10*x7 - x11*x9;
    res[4]=-x10*x4 + x11*x8;
    res[5]=(*this)[0]*(-A[0]*x7 + A[1]*x4);
    return res;
};

//-----------------------------------
// (R130B1, R130B3) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[2]*B[3] + A[3]*B[2];
    res[6]=A[1]*B[3] - A[3]*B[1];
    res[7]=-A[1]*B[2] + A[2]*B[1];
    res[8]=A[0]*B[1] + A[1]*B[0];
    res[9]=A[0]*B[2] + A[2]*B[0];
    res[10]=A[0]*B[3] + A[3]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B1<T> &A, const R130B3<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=A[0]*B[1] + A[1]*B[0];
    res[4]=A[0]*B[2] + A[2]*B[0];
    res[5]=A[0]*B[3] + A[3]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1<T> &A, const R130B3<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[2];
    T x7 = (*this)[3]*A[3];
    T x8 = 2.0*(*this)[1];
    T x9 = A[0]*x4;
    T x10 = A[1]*x8;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[1]*x5 + x4*x6 + x4*x7;
    res[1]=-A[0]*x5 - A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - x6*x8 - x7*x8;
    res[2]=-(*this)[2]*x10 - 2.0*(*this)[2]*x7 - (*this)[2]*x9 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3;
    res[3]=-(*this)[3]*x10 - 2.0*(*this)[3]*x6 - (*this)[3]*x9 - A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3;
    return res;
};

//--------------------------------------
// (R130B1, R130B3Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[2]*B[2] + A[3]*B[1];
    res[6]=A[1]*B[2] - A[3]*B[0];
    res[7]=-A[1]*B[1] + A[2]*B[0];
    res[8]=A[0]*B[0];
    res[9]=A[0]*B[1];
    res[10]=A[0]*B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B1<T> &A, const R130B3Sm1<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    res[3]=A[0]*B[0];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1<T> &A, const R130B3Sm1<T> &B) {
    R130B4<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B3<T> R130B1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*A[0];
    T x2 = (*this)[2]*A[1];
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[2], 2);
    T x5 = std::pow((*this)[3], 2);
    T x6 = std::pow((*this)[0], 2);
    T x7 = std::pow((*this)[1], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=x0*(x1 + x2 + x3);
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x6 - A[0]*x7 - x2*x8 - x3*x8;
    res[2]=-A[1]*x4 + A[1]*x5 - A[1]*x6 + A[1]*x7 - x1*x9 - x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 - x1*x10 - x10*x2;
    return res;
};

//--------------------------------------
// (R130B1, R130B3Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1]*B[0];
    res[9]=A[2]*B[0];
    res[10]=A[3]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B1<T> &A, const R130B3Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[1]*B[0];
    res[1]=A[2]*B[0];
    res[2]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1<T> &A, const R130B3Sp1<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0]*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    res[1]=-(*this)[1]*x0;
    res[2]=-(*this)[2]*x0;
    res[3]=-(*this)[3]*x0;
    return res;
};

//-----------------------------------
// (R130B1, R130B4) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
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
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
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
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator*(const R130B1<T> &A, const R130B4<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    res[3]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B1<T> &A, const R130B4<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    res[3]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2));
    return res;
};

//-----------------------------------
// (R130B1, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    res[4]=(*this)[3];
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

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    res[4]=A[3] + B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    res[4]=A[3] - B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[1] - A[1]*B[2] - A[2]*B[3] - A[3]*B[4];
    res[1]=A[0]*B[0] - A[1]*B[5] - A[2]*B[6] - A[3]*B[7];
    res[2]=-A[0]*B[5] + A[1]*B[0] - A[2]*B[10] + A[3]*B[9];
    res[3]=-A[0]*B[6] + A[1]*B[10] + A[2]*B[0] - A[3]*B[8];
    res[4]=-A[0]*B[7] - A[1]*B[9] + A[2]*B[8] + A[3]*B[0];
    res[5]=-A[0]*B[2] + A[1]*B[1] - A[2]*B[14] + A[3]*B[13];
    res[6]=-A[0]*B[3] + A[1]*B[14] + A[2]*B[1] - A[3]*B[12];
    res[7]=-A[0]*B[4] - A[1]*B[13] + A[2]*B[12] + A[3]*B[1];
    res[8]=A[0]*B[12] + A[1]*B[11] - A[2]*B[4] + A[3]*B[3];
    res[9]=A[0]*B[13] + A[1]*B[4] + A[2]*B[11] - A[3]*B[2];
    res[10]=A[0]*B[14] - A[1]*B[3] + A[2]*B[2] + A[3]*B[11];
    res[11]=A[0]*B[15] - A[1]*B[8] - A[2]*B[9] - A[3]*B[10];
    res[12]=A[0]*B[8] - A[1]*B[15] - A[2]*B[7] + A[3]*B[6];
    res[13]=A[0]*B[9] + A[1]*B[7] - A[2]*B[15] - A[3]*B[5];
    res[14]=A[0]*B[10] - A[1]*B[6] + A[2]*B[5] - A[3]*B[15];
    res[15]=A[0]*B[11] + A[1]*B[12] + A[2]*B[13] + A[3]*B[14];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[1] - A[1]*B[2] - A[2]*B[3] - A[3]*B[4];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=-A[2]*B[14] + A[3]*B[13];
    res[6]=A[1]*B[14] - A[3]*B[12];
    res[7]=-A[1]*B[13] + A[2]*B[12];
    res[8]=A[0]*B[12] + A[1]*B[11];
    res[9]=A[0]*B[13] + A[2]*B[11];
    res[10]=A[0]*B[14] + A[3]*B[11];
    res[11]=-A[1]*B[8] - A[2]*B[9] - A[3]*B[10];
    res[12]=A[0]*B[8] - A[2]*B[7] + A[3]*B[6];
    res[13]=A[0]*B[9] + A[1]*B[7] - A[3]*B[5];
    res[14]=A[0]*B[10] - A[1]*B[6] + A[2]*B[5];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[5] - A[2]*B[6] - A[3]*B[7];
    res[2]=-A[0]*B[5] - A[2]*B[10] + A[3]*B[9];
    res[3]=-A[0]*B[6] + A[1]*B[10] - A[3]*B[8];
    res[4]=-A[0]*B[7] - A[1]*B[9] + A[2]*B[8];
    res[5]=-A[0]*B[2] + A[1]*B[1];
    res[6]=-A[0]*B[3] + A[2]*B[1];
    res[7]=-A[0]*B[4] + A[3]*B[1];
    res[8]=-A[2]*B[4] + A[3]*B[3];
    res[9]=A[1]*B[4] - A[3]*B[2];
    res[10]=-A[1]*B[3] + A[2]*B[2];
    res[11]=A[0]*B[15];
    res[12]=-A[1]*B[15];
    res[13]=-A[2]*B[15];
    res[14]=-A[3]*B[15];
    res[15]=A[0]*B[11] + A[1]*B[12] + A[2]*B[13] + A[3]*B[14];
    return res;
};

template<typename T>
inline R130MV<T> R130B1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[3];
    T x7 = (*this)[3]*A[4];
    T x8 = 2.0*(*this)[1];
    T x9 = A[2]*x8;
    T x10 = 2.0*(*this)[2];
    T x11 = A[1]*x4;
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[2]*x8;
    T x15 = (*this)[3]*A[7];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*x10;
    res[0]=A[0]*(x0 + x1 + x2 - x3);
    res[1]=-A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x5 + x4*x6 + x4*x7;
    res[2]=-A[1]*x5 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 + x6*x8 + x7*x8;
    res[3]=-(*this)[2]*x11 + (*this)[2]*x9 - A[3]*x0 + A[3]*x1 - A[3]*x2 + A[3]*x3 + x10*x7;
    res[4]=-(*this)[3]*x11 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 - A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3;
    res[5]=A[10]*x12 - A[5]*x0 + A[5]*x1 + A[5]*x2 + A[5]*x3 - A[6]*x14 - A[9]*x13 - x15*x8;
    res[6]=-A[10]*x5 - A[5]*x14 + A[6]*x0 - A[6]*x1 + A[6]*x2 + A[6]*x3 + A[8]*x13 - x10*x15;
    res[7]=-A[5]*x16 - A[6]*x17 + A[7]*x0 + A[7]*x1 - A[7]*x2 + A[7]*x3 - A[8]*x12 + A[9]*x5;
    res[8]=A[10]*x16 - A[6]*x13 + A[7]*x12 + A[8]*x0 - A[8]*x1 - A[8]*x2 - A[8]*x3 + A[9]*x14;
    res[9]=A[10]*x17 + A[5]*x13 - A[7]*x5 + A[8]*x14 - A[9]*x0 + A[9]*x1 - A[9]*x2 - A[9]*x3;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 - A[10]*x3 - A[5]*x12 + A[6]*x5 + A[8]*x16 + A[9]*x17;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[12]*x5 + A[13]*x12 + A[14]*x13;
    res[12]=-A[11]*x5 - A[12]*x0 + A[12]*x1 + A[12]*x2 - A[12]*x3 - A[13]*x14 - A[14]*x16;
    res[13]=-A[11]*x12 - A[12]*x14 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 - A[14]*x17;
    res[14]=-A[11]*x13 - A[12]*x16 - A[13]*x17 + A[14]*x0 + A[14]*x1 - A[14]*x2 - A[14]*x3;
    res[15]=A[15]*(-x0 - x1 - x2 + x3);
    return res;
};

//-------------------------------------
// (R130B1, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0] - A[2]*B[3] + A[3]*B[2];
    res[3]=A[1]*B[3] + A[2]*B[0] - A[3]*B[1];
    res[4]=-A[1]*B[2] + A[2]*B[1] + A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1<T> &A, const Rotation<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> R130B1<T>::conjugate(const Rotation<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[3]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = (*this)[1]*A[2];
    T x8 = (*this)[2]*A[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*A[3];
    res[0]=A[0]*(x0 + x1 + x2 - x3);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[2]*x5 - A[2]*x6;
    res[6]=-(*this)[1]*x5 + A[1]*x6;
    res[7]=x4*(x7 - x8);
    res[8]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + x10*x11 + x7*x9;
    res[9]=-A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + x10*x8 + x11*x9;
    res[10]=(*this)[3]*A[1]*x10 + (*this)[3]*A[2]*x9 - A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//----------------------------------
// (R130B1, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=B[4];
    res[9]=B[5];
    res[10]=B[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    res[4]=A[3];
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=-B[4];
    res[9]=-B[5];
    res[10]=-B[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=-A[0]*B[1] + A[1]*B[0] - A[2]*B[6] + A[3]*B[5];
    res[3]=-A[0]*B[2] + A[1]*B[6] + A[2]*B[0] - A[3]*B[4];
    res[4]=-A[0]*B[3] - A[1]*B[5] + A[2]*B[4] + A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[7] - A[1]*B[4] - A[2]*B[5] - A[3]*B[6];
    res[12]=A[0]*B[4] - A[1]*B[7] - A[2]*B[3] + A[3]*B[2];
    res[13]=A[0]*B[5] + A[1]*B[3] - A[2]*B[7] - A[3]*B[1];
    res[14]=A[0]*B[6] - A[1]*B[2] + A[2]*B[1] - A[3]*B[7];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[4] - A[2]*B[5] - A[3]*B[6];
    res[12]=A[0]*B[4] - A[2]*B[3] + A[3]*B[2];
    res[13]=A[0]*B[5] + A[1]*B[3] - A[3]*B[1];
    res[14]=A[0]*B[6] - A[1]*B[2] + A[2]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=-A[0]*B[1] - A[2]*B[6] + A[3]*B[5];
    res[3]=-A[0]*B[2] + A[1]*B[6] - A[3]*B[4];
    res[4]=-A[0]*B[3] - A[1]*B[5] + A[2]*B[4];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[7];
    res[12]=-A[1]*B[7];
    res[13]=-A[2]*B[7];
    res[14]=-A[3]*B[7];
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[6]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = 2.0*(*this)[1];
    T x8 = (*this)[2]*x7;
    T x9 = (*this)[3]*A[3];
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[1]*x4;
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x7;
    T x14 = (*this)[3]*x10;
    res[0]=A[0]*(x0 + x1 + x2 - x3);
    res[1]=(*this)[2]*x5 - A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x8 - A[5]*x6 - x7*x9;
    res[2]=-(*this)[1]*x5 - A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 + A[2]*x3 + A[4]*x6 - x10*x9;
    res[3]=-A[1]*x13 - A[2]*x14 + A[3]*x0 + A[3]*x1 - A[3]*x2 + A[3]*x3 - A[4]*x12 + A[5]*x11;
    res[4]=-A[2]*x6 + A[3]*x12 + A[4]*x0 - A[4]*x1 - A[4]*x2 - A[4]*x3 + A[5]*x8 + A[6]*x13;
    res[5]=A[1]*x6 - A[3]*x11 + A[4]*x8 - A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[6]*x14;
    res[6]=-A[1]*x12 + A[2]*x11 + A[4]*x13 + A[5]*x14 - A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3;
    res[7]=A[7]*(-x0 - x1 - x2 + x3);
    return res;
};

//-----------------------------------
// (R130B1, scalar) binary operations
//-----------------------------------

template<typename T>
inline R130B1<T> operator*(const R130B1<T> &A, const T &B) {
    R130B1<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    return res;
};
template<typename T>
inline R130B1<T> operator*(const T &A, const R130B1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B1<T> operator/(const R130B1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//-------------------------------------
// (R130B1Sm1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
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

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
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

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[1]*B[3] + A[2]*B[2];
    res[13]=A[0]*B[3] - A[2]*B[1];
    res[14]=-A[0]*B[2] + A[1]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[1]*B[3] + A[2]*B[2];
    res[13]=A[0]*B[3] - A[2]*B[1];
    res[14]=-A[0]*B[2] + A[1]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator^(const R130B1Sm1<T> &A, const Boost<T> &B) {
    R130B1Sp1<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    return res;
};

template<typename T>
inline Boost<T> R130B1Sm1<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-(*this)[2]*A[1]*x3 - (*this)[2]*A[2]*x6 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    return res;
};

//--------------------------------------
// (R130B1Sm1, R130B0) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
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

template<typename T>
inline R130B1Sm1<T> operator*(const R130B1Sm1<T> &A, const R130B0<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator/(const R130B1Sm1<T> &A, const R130B0<T> &B) {
    R130B1Sm1<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B1Sm1<T> &A, const R130B0<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B1Sm1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B1Sm1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sm1<T>::operator R130B1<T>() const {
    R130B1<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    return res;
};

template<typename T>
inline R130B1<T> operator+(const R130B1Sm1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1Sm1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0];
    res[6]=A[1]*B[0];
    res[7]=A[2]*B[0];
    res[8]=-A[1]*B[3] + A[2]*B[2];
    res[9]=A[0]*B[3] - A[2]*B[1];
    res[10]=-A[0]*B[2] + A[1]*B[1];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1Sm1<T> &A, const R130B1<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B1Sm1<T> &A, const R130B1<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=-A[1]*B[3] + A[2]*B[2];
    res[4]=A[0]*B[3] - A[2]*B[1];
    res[5]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline R130B1<T> R130B1Sm1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B1Sm1, R130B1Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130B1Sm1<T> operator+(const R130B1Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator-(const R130B1Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const R130B1Sm1<T> &A, const R130B1Sm1<T> &B) {
    Rotation<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const R130B1Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B1Sm1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B1Sm1, R130B1Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130B1<T> operator+(const R130B1Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B1<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B1<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B1Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B1Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B1Sm1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B1Sm1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[2]=-A[1]*B[5] + A[2]*B[4];
    res[3]=A[0]*B[5] - A[2]*B[3];
    res[4]=-A[0]*B[4] + A[1]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[3] - A[1]*B[4] - A[2]*B[5];
    res[12]=-A[1]*B[2] + A[2]*B[1];
    res[13]=A[0]*B[2] - A[2]*B[0];
    res[14]=-A[0]*B[1] + A[1]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B1Sm1<T> &A, const R130B2<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[3] - A[1]*B[4] - A[2]*B[5];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B1Sm1<T> &A, const R130B2<T> &B) {
    R130B1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[1]=-A[1]*B[5] + A[2]*B[4];
    res[2]=A[0]*B[5] - A[2]*B[3];
    res[3]=-A[0]*B[4] + A[1]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> R130B1Sm1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-A[0]*x7 - A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    res[3]=-A[3]*x0 - A[3]*x1 + A[3]*x2 + A[4]*x4 + A[5]*x7;
    res[4]=A[3]*x4 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[5]*x8;
    res[5]=A[3]*x7 + A[4]*x8 - A[5]*x0 + A[5]*x1 - A[5]*x2;
    return res;
};

//-----------------------------------------
// (R130B1Sm1, R130B2Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=A[0]*B[2] - A[2]*B[0];
    res[4]=-A[0]*B[1] + A[1]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator|(const R130B1Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B1Sm1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B1Sm1, R130B2Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
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

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
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
    res[12]=-A[1]*B[2] + A[2]*B[1];
    res[13]=A[0]*B[2] - A[2]*B[0];
    res[14]=-A[0]*B[1] + A[1]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B1Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator^(const R130B1Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B1Sm1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-(*this)[2]*A[0]*x3 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//--------------------------------------
// (R130B1Sm1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[3] + A[2]*B[2];
    res[6]=A[0]*B[3] - A[2]*B[1];
    res[7]=-A[0]*B[2] + A[1]*B[1];
    res[8]=A[0]*B[0];
    res[9]=A[1]*B[0];
    res[10]=A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B1Sm1<T> &A, const R130B3<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    res[3]=A[0]*B[0];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1Sm1<T> &A, const R130B3<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B1Sm1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-(*this)[2]*A[1]*x3 - (*this)[2]*A[2]*x6 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B1Sm1, R130B3Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[2] + A[2]*B[1];
    res[6]=A[0]*B[2] - A[2]*B[0];
    res[7]=-A[0]*B[1] + A[1]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator|(const R130B1Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B1Sm1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-(*this)[2]*A[0]*x3 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B1Sm1, R130B3Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B1Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B1Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B1Sm1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B1Sm1, R130B4) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
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
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
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
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B1Sm1<T> &A, const R130B4<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B1Sm1<T> &A, const R130B4<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B1Sm1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B1Sm1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sm1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=(*this)[0];
    res[3]=(*this)[1];
    res[4]=(*this)[2];
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

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=A[0] + B[2];
    res[3]=A[1] + B[3];
    res[4]=A[2] + B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=A[0] - B[2];
    res[3]=A[1] - B[3];
    res[4]=A[2] - B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[2] - A[1]*B[3] - A[2]*B[4];
    res[1]=-A[0]*B[5] - A[1]*B[6] - A[2]*B[7];
    res[2]=A[0]*B[0] - A[1]*B[10] + A[2]*B[9];
    res[3]=A[0]*B[10] + A[1]*B[0] - A[2]*B[8];
    res[4]=-A[0]*B[9] + A[1]*B[8] + A[2]*B[0];
    res[5]=A[0]*B[1] - A[1]*B[14] + A[2]*B[13];
    res[6]=A[0]*B[14] + A[1]*B[1] - A[2]*B[12];
    res[7]=-A[0]*B[13] + A[1]*B[12] + A[2]*B[1];
    res[8]=A[0]*B[11] - A[1]*B[4] + A[2]*B[3];
    res[9]=A[0]*B[4] + A[1]*B[11] - A[2]*B[2];
    res[10]=-A[0]*B[3] + A[1]*B[2] + A[2]*B[11];
    res[11]=-A[0]*B[8] - A[1]*B[9] - A[2]*B[10];
    res[12]=-A[0]*B[15] - A[1]*B[7] + A[2]*B[6];
    res[13]=A[0]*B[7] - A[1]*B[15] - A[2]*B[5];
    res[14]=-A[0]*B[6] + A[1]*B[5] - A[2]*B[15];
    res[15]=A[0]*B[12] + A[1]*B[13] + A[2]*B[14];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[2] - A[1]*B[3] - A[2]*B[4];
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=-A[1]*B[14] + A[2]*B[13];
    res[6]=A[0]*B[14] - A[2]*B[12];
    res[7]=-A[0]*B[13] + A[1]*B[12];
    res[8]=A[0]*B[11];
    res[9]=A[1]*B[11];
    res[10]=A[2]*B[11];
    res[11]=-A[0]*B[8] - A[1]*B[9] - A[2]*B[10];
    res[12]=-A[1]*B[7] + A[2]*B[6];
    res[13]=A[0]*B[7] - A[2]*B[5];
    res[14]=-A[0]*B[6] + A[1]*B[5];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[5] - A[1]*B[6] - A[2]*B[7];
    res[2]=-A[1]*B[10] + A[2]*B[9];
    res[3]=A[0]*B[10] - A[2]*B[8];
    res[4]=-A[0]*B[9] + A[1]*B[8];
    res[5]=A[0]*B[1];
    res[6]=A[1]*B[1];
    res[7]=A[2]*B[1];
    res[8]=-A[1]*B[4] + A[2]*B[3];
    res[9]=A[0]*B[4] - A[2]*B[2];
    res[10]=-A[0]*B[3] + A[1]*B[2];
    res[11]=0;
    res[12]=-A[0]*B[15];
    res[13]=-A[1]*B[15];
    res[14]=-A[2]*B[15];
    res[15]=A[0]*B[12] + A[1]*B[13] + A[2]*B[14];
    return res;
};

template<typename T>
inline R130MV<T> R130B1Sm1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[4];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*(x0 + x1 + x2);
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=-A[5]*x0 + A[5]*x1 + A[5]*x2 - A[6]*x4 - A[7]*x7;
    res[6]=-A[5]*x4 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[7]*x8;
    res[7]=-A[5]*x7 - A[6]*x8 + A[7]*x0 + A[7]*x1 - A[7]*x2;
    res[8]=A[10]*x7 + A[8]*x0 - A[8]*x1 - A[8]*x2 + A[9]*x4;
    res[9]=A[10]*x8 + A[8]*x4 - A[9]*x0 + A[9]*x1 - A[9]*x2;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 + A[8]*x7 + A[9]*x8;
    res[11]=A[11]*(x0 + x1 + x2);
    res[12]=-A[12]*x0 + A[12]*x1 + A[12]*x2 - A[13]*x4 - A[14]*x7;
    res[13]=-A[12]*x4 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[14]*x8;
    res[14]=-A[12]*x7 - A[13]*x8 + A[14]*x0 + A[14]*x1 - A[14]*x2;
    res[15]=-A[15]*(x0 + x1 + x2);
    return res;
};

//----------------------------------------
// (R130B1Sm1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0] - A[1]*B[3] + A[2]*B[2];
    res[3]=A[0]*B[3] + A[1]*B[0] - A[2]*B[1];
    res[4]=-A[0]*B[2] + A[1]*B[1] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1Sm1<T> &A, const Rotation<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline Rotation<T> R130B1Sm1<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-------------------------------------
// (R130B1Sm1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=B[4];
    res[9]=B[5];
    res[10]=B[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=A[0];
    res[3]=A[1];
    res[4]=A[2];
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=-B[4];
    res[9]=-B[5];
    res[10]=-B[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[2]=A[0]*B[0] - A[1]*B[6] + A[2]*B[5];
    res[3]=A[0]*B[6] + A[1]*B[0] - A[2]*B[4];
    res[4]=-A[0]*B[5] + A[1]*B[4] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[4] - A[1]*B[5] - A[2]*B[6];
    res[12]=-A[0]*B[7] - A[1]*B[3] + A[2]*B[2];
    res[13]=A[0]*B[3] - A[1]*B[7] - A[2]*B[1];
    res[14]=-A[0]*B[2] + A[1]*B[1] - A[2]*B[7];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[4] - A[1]*B[5] - A[2]*B[6];
    res[12]=-A[1]*B[3] + A[2]*B[2];
    res[13]=A[0]*B[3] - A[2]*B[1];
    res[14]=-A[0]*B[2] + A[1]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[2]=-A[1]*B[6] + A[2]*B[5];
    res[3]=A[0]*B[6] - A[2]*B[4];
    res[4]=-A[0]*B[5] + A[1]*B[4];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[0]*B[7];
    res[13]=-A[1]*B[7];
    res[14]=-A[2]*B[7];
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B1Sm1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-A[1]*x7 - A[2]*x8 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    res[4]=A[4]*x0 - A[4]*x1 - A[4]*x2 + A[5]*x4 + A[6]*x7;
    res[5]=A[4]*x4 - A[5]*x0 + A[5]*x1 - A[5]*x2 + A[6]*x8;
    res[6]=A[4]*x7 + A[5]*x8 - A[6]*x0 - A[6]*x1 + A[6]*x2;
    res[7]=-A[7]*(x0 + x1 + x2);
    return res;
};

//--------------------------------------
// (R130B1Sm1, scalar) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sm1<T> operator*(const R130B1Sm1<T> &A, const T &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    return res;
};
template<typename T>
inline R130B1Sm1<T> operator*(const T &A, const R130B1Sm1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B1Sm1<T> operator/(const R130B1Sm1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//-------------------------------------
// (R130B1Sp1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
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

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
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

template<typename T>
inline R130B1<T> operator*(const R130B1Sp1<T> &A, const Boost<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    res[3]=-A[0]*B[3];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator|(const R130B1Sp1<T> &A, const Boost<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1Sp1<T> &A, const Boost<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[0]*B[1];
    res[1]=-A[0]*B[2];
    res[2]=-A[0]*B[3];
    return res;
};

template<typename T>
inline Boost<T> R130B1Sp1<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//--------------------------------------
// (R130B1Sp1, R130B0) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
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

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
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

template<typename T>
inline R130B1Sp1<T> operator*(const R130B1Sp1<T> &A, const R130B0<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator/(const R130B1Sp1<T> &A, const R130B0<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]/B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator|(const R130B1Sp1<T> &A, const R130B0<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B1Sp1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=-std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B1Sp1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sp1<T>::operator R130B1<T>() const {
    R130B1<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator+(const R130B1Sp1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1Sp1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    return res;
};

template<typename T>
inline Boost<T> operator*(const R130B1Sp1<T> &A, const R130B1<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    res[3]=-A[0]*B[3];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1Sp1<T> &A, const R130B1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B1Sp1<T> &A, const R130B1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[1];
    res[1]=-A[0]*B[2];
    res[2]=-A[0]*B[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B1Sp1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//-----------------------------------------
// (R130B1Sp1, R130B1Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130B1<T> operator+(const R130B1Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B1<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    return res;
};

template<typename T>
inline R130B1<T> operator-(const R130B1Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B1<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B1Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B1Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B1Sp1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//-----------------------------------------
// (R130B1Sp1, R130B1Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130B1Sp1<T> operator+(const R130B1Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0] + B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator-(const R130B1Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0] - B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator*(const R130B1Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B1Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B1Sp1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=-std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B1Sp1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[0]*B[0];
    res[3]=-A[0]*B[1];
    res[4]=-A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[3];
    res[13]=A[0]*B[4];
    res[14]=A[0]*B[5];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B1Sp1<T> &A, const R130B2<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[3];
    res[1]=A[0]*B[4];
    res[2]=A[0]*B[5];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1Sp1<T> &A, const R130B2<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> R130B1Sp1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=-A[3]*x0;
    res[4]=-A[4]*x0;
    res[5]=-A[5]*x0;
    return res;
};

//-----------------------------------------
// (R130B1Sp1, R130B2Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B1Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B1Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B1Sp1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    return res;
};

//-----------------------------------------
// (R130B1Sp1, R130B2Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
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

template<typename T>
inline R130B1Sm1<T> operator*(const R130B1Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B1Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B1Sp1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B1Sp1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sp1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[1];
    res[9]=A[0]*B[2];
    res[10]=A[0]*B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B1Sp1<T> &A, const R130B3<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[1];
    res[1]=A[0]*B[2];
    res[2]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1Sp1<T> &A, const R130B3<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B1Sp1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    return res;
};

//-----------------------------------------
// (R130B1Sp1, R130B3Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
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
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
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
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B1Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B1Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B1Sp1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    return res;
};

//-----------------------------------------
// (R130B1Sp1, R130B3Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> operator*(const R130B1Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B1Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B1Sp1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B1Sp1, R130B4) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
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
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0];
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
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator*(const R130B1Sp1<T> &A, const R130B4<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator^(const R130B1Sp1<T> &A, const R130B4<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B1Sp1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B1Sp1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sp1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=(*this)[0];
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

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[1];
    res[1]=A[0]*B[0];
    res[2]=-A[0]*B[5];
    res[3]=-A[0]*B[6];
    res[4]=-A[0]*B[7];
    res[5]=-A[0]*B[2];
    res[6]=-A[0]*B[3];
    res[7]=-A[0]*B[4];
    res[8]=A[0]*B[12];
    res[9]=A[0]*B[13];
    res[10]=A[0]*B[14];
    res[11]=A[0]*B[15];
    res[12]=A[0]*B[8];
    res[13]=A[0]*B[9];
    res[14]=A[0]*B[10];
    res[15]=A[0]*B[11];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[1];
    res[1]=A[0]*B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[12];
    res[9]=A[0]*B[13];
    res[10]=A[0]*B[14];
    res[11]=0;
    res[12]=A[0]*B[8];
    res[13]=A[0]*B[9];
    res[14]=A[0]*B[10];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[0]*B[5];
    res[3]=-A[0]*B[6];
    res[4]=-A[0]*B[7];
    res[5]=-A[0]*B[2];
    res[6]=-A[0]*B[3];
    res[7]=-A[0]*B[4];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[15];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[11];
    return res;
};

template<typename T>
inline R130MV<T> R130B1Sp1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    res[6]=A[6]*x0;
    res[7]=A[7]*x0;
    res[8]=-A[8]*x0;
    res[9]=-A[9]*x0;
    res[10]=-A[10]*x0;
    res[11]=A[11]*x0;
    res[12]=-A[12]*x0;
    res[13]=-A[13]*x0;
    res[14]=-A[14]*x0;
    res[15]=A[15]*x0;
    return res;
};

//----------------------------------------
// (R130B1Sp1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
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
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
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
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Rotation<T> R130B1Sp1<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    return res;
};

//-------------------------------------
// (R130B1Sp1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B1Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=B[4];
    res[9]=B[5];
    res[10]=B[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B1Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=-B[4];
    res[9]=-B[5];
    res[10]=-B[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B1Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=-A[0]*B[1];
    res[3]=-A[0]*B[2];
    res[4]=-A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[7];
    res[12]=A[0]*B[4];
    res[13]=A[0]*B[5];
    res[14]=A[0]*B[6];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B1Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
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
    res[12]=A[0]*B[4];
    res[13]=A[0]*B[5];
    res[14]=A[0]*B[6];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B1Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[0]*B[1];
    res[3]=-A[0]*B[2];
    res[4]=-A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[7];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B1Sp1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=-A[4]*x0;
    res[5]=-A[5]*x0;
    res[6]=-A[6]*x0;
    res[7]=A[7]*x0;
    return res;
};

//--------------------------------------
// (R130B1Sp1, scalar) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sp1<T> operator*(const R130B1Sp1<T> &A, const T &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B;
    return res;
};
template<typename T>
inline R130B1Sp1<T> operator*(const T &A, const R130B1Sp1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B1Sp1<T> operator/(const R130B1Sp1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//----------------------------------
// (R130B2, Boost) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0] + B[1];
    res[6]=A[1] + B[2];
    res[7]=A[2] + B[3];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0] - B[1];
    res[6]=A[1] - B[2];
    res[7]=A[2] - B[3];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0] - A[4]*B[3] + A[5]*B[2];
    res[2]=A[1]*B[0] + A[3]*B[3] - A[5]*B[1];
    res[3]=A[2]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[4]=A[1]*B[3] - A[2]*B[2] + A[3]*B[0];
    res[5]=-A[0]*B[3] + A[2]*B[1] + A[4]*B[0];
    res[6]=A[0]*B[2] - A[1]*B[1] + A[5]*B[0];
    res[7]=A[3]*B[1] + A[4]*B[2] + A[5]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B2<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=A[4]*B[0];
    res[6]=A[5]*B[0];
    res[7]=A[3]*B[1] + A[4]*B[2] + A[5]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2<T> &A, const Boost<T> &B) {
    R130B2<T> res;
    res[0]=-A[4]*B[3] + A[5]*B[2];
    res[1]=A[3]*B[3] - A[5]*B[1];
    res[2]=-A[3]*B[2] + A[4]*B[1];
    res[3]=A[1]*B[3] - A[2]*B[2];
    res[4]=-A[0]*B[3] + A[2]*B[1];
    res[5]=A[0]*B[2] - A[1]*B[1];
    return res;
};

template<typename T>
inline Rotor<T> R130B2<T>::conjugate(const Boost<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[3], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[4]*x6;
    T x8 = (*this)[5]*A[3];
    T x9 = 2.0*(*this)[0];
    T x10 = (*this)[1]*A[2];
    T x11 = (*this)[2]*A[3];
    T x12 = 2.0*(*this)[4];
    T x13 = A[1]*x9;
    T x14 = 2.0*(*this)[1];
    T x15 = (*this)[5]*A[1];
    T x16 = (*this)[5]*A[2];
    T x17 = 2.0*(*this)[2];
    T x18 = (*this)[1]*A[1];
    T x19 = (*this)[0]*x6;
    T x20 = (*this)[1]*x12;
    res[0]=A[0]*(x0 + x1 + x2 - x3 - x4 - x5);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[2]*x7 - x10*x9 - x11*x9 + x6*x8;
    res[2]=-(*this)[1]*x13 + A[1]*x7 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - x11*x14 + x12*x8;
    res[3]=-(*this)[2]*x13 - A[3]*x0 - A[3]*x1 + A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 - x10*x17 + x12*x16 + x15*x6;
    res[4]=-(*this)[4]*A[2]*x9 - A[1]*x19 - x10*x6 - x11*x6 + x12*x18 + x15*x17 - x8*x9;
    res[5]=-(*this)[4]*x13 + A[2]*x19 - x10*x12 - x11*x12 - x14*x8 + x16*x17 - x18*x6;
    res[6]=-(*this)[2]*A[1]*x6 - (*this)[2]*A[2]*x12 - 2.0*(*this)[5]*x10 - (*this)[5]*x13 + A[3]*x19 + A[3]*x20 - x17*x8;
    res[7]=-A[0]*((*this)[5]*x17 + x19 + x20);
    return res;
};

//-----------------------------------
// (R130B2, R130B0) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2<T> operator*(const R130B2<T> &A, const R130B0<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[4]*B[0];
    res[5]=A[5]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> operator/(const R130B2<T> &A, const R130B0<T> &B) {
    R130B2<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B2<T> &A, const R130B0<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[4]*B[0];
    res[5]=A[5]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B2<T>::conjugate(const R130B0<T> &A) const {
    R130MV<T> res;
    res[0]=A[0]*(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2));
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
    res[15]=-2.0*A[0]*((*this)[0]*(*this)[3] + (*this)[1]*(*this)[4] + (*this)[2]*(*this)[5]);
    return res;
};

//-----------------------------------
// (R130B2, R130B1) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[2]=A[0]*B[0] - A[4]*B[3] + A[5]*B[2];
    res[3]=A[1]*B[0] + A[3]*B[3] - A[5]*B[1];
    res[4]=A[2]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[3]*B[1] - A[4]*B[2] - A[5]*B[3];
    res[12]=A[1]*B[3] - A[2]*B[2] + A[3]*B[0];
    res[13]=-A[0]*B[3] + A[2]*B[1] + A[4]*B[0];
    res[14]=A[0]*B[2] - A[1]*B[1] + A[5]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B2<T> &A, const R130B1<T> &B) {
    R130B3<T> res;
    res[0]=-A[3]*B[1] - A[4]*B[2] - A[5]*B[3];
    res[1]=A[1]*B[3] - A[2]*B[2] + A[3]*B[0];
    res[2]=-A[0]*B[3] + A[2]*B[1] + A[4]*B[0];
    res[3]=A[0]*B[2] - A[1]*B[1] + A[5]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B2<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0] - A[4]*B[3] + A[5]*B[2];
    res[2]=A[1]*B[0] + A[3]*B[3] - A[5]*B[1];
    res[3]=A[2]*B[0] - A[3]*B[2] + A[4]*B[1];
    return res;
};

template<typename T>
inline R130B1<T> R130B2<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = std::pow((*this)[4], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = 2.0*(*this)[0];
    T x7 = A[2]*x6;
    T x8 = 2.0*(*this)[1];
    T x9 = (*this)[3]*A[3];
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[4]*x10;
    T x12 = A[3]*x6;
    T x13 = (*this)[5]*x8;
    T x14 = (*this)[3]*A[2];
    T x15 = 2.0*(*this)[4];
    T x16 = 2.0*(*this)[5];
    T x17 = A[1]*x6;
    T x18 = (*this)[2]*x8;
    T x19 = (*this)[3]*A[0];
    T x20 = (*this)[3]*A[1];
    T x21 = (*this)[5]*x15;
    T x22 = A[0]*x6;
    res[0]=-(*this)[4]*x12 + (*this)[5]*x7 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[1]*x11 - A[1]*x13 - x10*x14 + x8*x9;
    res[1]=(*this)[1]*x7 + (*this)[2]*x12 - A[0]*x11 + A[0]*x13 + A[1]*x0 - A[1]*x1 - A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + x14*x15 + x16*x9;
    res[2]=(*this)[1]*x17 - (*this)[5]*x22 - A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[3]*x18 + A[3]*x21 + x10*x19 + x15*x20;
    res[3]=(*this)[2]*x17 + (*this)[4]*x22 + A[2]*x18 + A[2]*x21 - A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 + x16*x20 - x19*x8;
    return res;
};

//--------------------------------------
// (R130B2, R130B1Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[2]=-A[4]*B[2] + A[5]*B[1];
    res[3]=A[3]*B[2] - A[5]*B[0];
    res[4]=-A[3]*B[1] + A[4]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[3]*B[0] - A[4]*B[1] - A[5]*B[2];
    res[12]=A[1]*B[2] - A[2]*B[1];
    res[13]=-A[0]*B[2] + A[2]*B[0];
    res[14]=A[0]*B[1] - A[1]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B2<T> &A, const R130B1Sm1<T> &B) {
    R130B3<T> res;
    res[0]=-A[3]*B[0] - A[4]*B[1] - A[5]*B[2];
    res[1]=A[1]*B[2] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B2<T> &A, const R130B1Sm1<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[1]=-A[4]*B[2] + A[5]*B[1];
    res[2]=A[3]*B[2] - A[5]*B[0];
    res[3]=-A[3]*B[1] + A[4]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B2<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[1]*x0;
    T x2 = 2.0*(*this)[1];
    T x3 = (*this)[3]*A[2];
    T x4 = 2.0*(*this)[2];
    T x5 = (*this)[4]*A[0];
    T x6 = A[2]*x0;
    T x7 = (*this)[5]*A[0];
    T x8 = (*this)[3]*A[1];
    T x9 = std::pow((*this)[0], 2);
    T x10 = std::pow((*this)[3], 2);
    T x11 = std::pow((*this)[1], 2);
    T x12 = std::pow((*this)[2], 2);
    T x13 = std::pow((*this)[4], 2);
    T x14 = std::pow((*this)[5], 2);
    T x15 = 2.0*(*this)[4];
    T x16 = A[0]*x0;
    T x17 = (*this)[2]*x2;
    T x18 = 2.0*(*this)[3];
    T x19 = (*this)[5]*x15;
    res[0]=-(*this)[4]*x6 + (*this)[5]*x1 + x2*x3 - x2*x7 + x4*x5 - x4*x8;
    res[1]=(*this)[1]*x1 + (*this)[2]*x6 + 2.0*(*this)[5]*x3 + A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 + A[0]*x9 + x15*x8;
    res[2]=(*this)[1]*x16 - A[1]*x10 + A[1]*x11 - A[1]*x12 + A[1]*x13 - A[1]*x14 - A[1]*x9 + A[2]*x17 + A[2]*x19 + x18*x5;
    res[3]=(*this)[2]*x16 + A[1]*x17 + A[1]*x19 - A[2]*x10 - A[2]*x11 + A[2]*x12 - A[2]*x13 + A[2]*x14 - A[2]*x9 + x18*x7;
    return res;
};

//--------------------------------------
// (R130B2, R130B1Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[3]*B[0];
    res[13]=A[4]*B[0];
    res[14]=A[5]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B2<T> &A, const R130B1Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[3]*B[0];
    res[1]=A[4]*B[0];
    res[2]=A[5]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B2<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B2<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2));
    res[1]=x0*((*this)[1]*(*this)[5] - (*this)[2]*(*this)[4]);
    res[2]=x0*(-(*this)[0]*(*this)[5] + (*this)[2]*(*this)[3]);
    res[3]=x0*((*this)[0]*(*this)[4] - (*this)[1]*(*this)[3]);
    return res;
};

//-----------------------------------
// (R130B2, R130B2) binary operations
//-----------------------------------

template<typename T>
inline R130B2<T> operator+(const R130B2<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    res[4]=A[4] + B[4];
    res[5]=A[5] + B[5];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    res[4]=A[4] - B[4];
    res[5]=A[5] - B[5];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] - A[3]*B[3] - A[4]*B[4] - A[5]*B[5];
    res[1]=-A[1]*B[5] + A[2]*B[4] - A[4]*B[2] + A[5]*B[1];
    res[2]=A[0]*B[5] - A[2]*B[3] + A[3]*B[2] - A[5]*B[0];
    res[3]=-A[0]*B[4] + A[1]*B[3] - A[3]*B[1] + A[4]*B[0];
    res[4]=A[1]*B[2] - A[2]*B[1] - A[4]*B[5] + A[5]*B[4];
    res[5]=-A[0]*B[2] + A[2]*B[0] + A[3]*B[5] - A[5]*B[3];
    res[6]=A[0]*B[1] - A[1]*B[0] - A[3]*B[4] + A[4]*B[3];
    res[7]=A[0]*B[3] + A[1]*B[4] + A[2]*B[5] + A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] - A[3]*B[3] - A[4]*B[4] - A[5]*B[5];
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
    res[15]=A[0]*B[3] + A[1]*B[4] + A[2]*B[5] + A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[5] + A[2]*B[4] - A[4]*B[2] + A[5]*B[1];
    res[1]=A[0]*B[5] - A[2]*B[3] + A[3]*B[2] - A[5]*B[0];
    res[2]=-A[0]*B[4] + A[1]*B[3] - A[3]*B[1] + A[4]*B[0];
    res[3]=A[1]*B[2] - A[2]*B[1] - A[4]*B[5] + A[5]*B[4];
    res[4]=-A[0]*B[2] + A[2]*B[0] + A[3]*B[5] - A[5]*B[3];
    res[5]=A[0]*B[1] - A[1]*B[0] - A[3]*B[4] + A[4]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> R130B2<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = std::pow((*this)[4], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = 2.0*(*this)[0];
    T x7 = A[3]*x6;
    T x8 = (*this)[4]*A[4];
    T x9 = (*this)[5]*A[5];
    T x10 = 2.0*(*this)[3];
    T x11 = (*this)[1]*A[4];
    T x12 = (*this)[2]*A[5];
    T x13 = (*this)[4]*x10;
    T x14 = (*this)[5]*A[2];
    T x15 = (*this)[1]*x6;
    T x16 = (*this)[2]*A[2];
    T x17 = 2.0*A[3];
    T x18 = (*this)[1]*(*this)[4];
    T x19 = (*this)[2]*(*this)[5];
    T x20 = A[3]*x10;
    T x21 = 2.0*(*this)[1];
    T x22 = 2.0*(*this)[4];
    T x23 = (*this)[3]*x6;
    T x24 = 2.0*x19;
    T x25 = 2.0*(*this)[5];
    T x26 = 2.0*(*this)[2];
    T x27 = A[0]*x10;
    T x28 = (*this)[5]*A[1];
    T x29 = A[0]*x6;
    T x30 = (*this)[2]*A[1];
    T x31 = 2.0*x18;
    res[0]=(*this)[3]*x7 + A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[0]*x4 - A[0]*x5 + A[1]*x13 - A[1]*x15 + x10*x11 + x10*x12 + x10*x14 - x16*x6 - x17*x18 - x17*x19 + x6*x8 + x6*x9;
    res[1]=(*this)[1]*x20 + (*this)[4]*x7 + A[0]*x13 - A[0]*x15 - A[1]*x0 + A[1]*x1 - A[1]*x2 + A[1]*x3 + A[1]*x4 - A[1]*x5 - A[4]*x23 - A[4]*x24 + x12*x22 + x14*x22 - x16*x21 + x21*x8 + x21*x9;
    res[2]=(*this)[2]*x20 - (*this)[2]*x29 + (*this)[5]*x27 + (*this)[5]*x7 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - A[5]*x23 - A[5]*x31 + x11*x25 - x21*x30 + x22*x28 + x26*x8 + x26*x9;
    res[3]=-(*this)[1]*A[1]*x10 - (*this)[4]*A[1]*x6 - A[0]*x23 + A[0]*x24 + A[0]*x31 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 - A[3]*x4 - A[3]*x5 - x10*x16 + x10*x8 + x10*x9 - x11*x6 - x12*x6 - x14*x6;
    res[4]=-(*this)[1]*x27 - (*this)[1]*x7 - (*this)[4]*x29 + A[1]*x23 + A[1]*x24 - A[1]*x31 + A[3]*x13 - A[4]*x0 + A[4]*x1 - A[4]*x2 + A[4]*x3 + A[4]*x4 - A[4]*x5 - x12*x21 - x14*x21 - x16*x22 + x22*x9;
    res[5]=-(*this)[2]*x27 - (*this)[2]*x7 + (*this)[5]*x20 - (*this)[5]*x29 + A[2]*x23 + A[2]*x31 + A[5]*x0 - A[5]*x1 - A[5]*x2 + A[5]*x3 - A[5]*x4 + A[5]*x5 - x11*x26 - x14*x26 - x21*x28 - x22*x30 + x25*x8;
    return res;
};

//--------------------------------------
// (R130B2, R130B2Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130B2<T> operator+(const R130B2<T> &A, const R130B2Sm1<T> &B) {
    R130B2<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3] + B[0];
    res[4]=A[4] + B[1];
    res[5]=A[5] + B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2<T> &A, const R130B2Sm1<T> &B) {
    R130B2<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3] - B[0];
    res[4]=A[4] - B[1];
    res[5]=A[5] - B[2];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2<T> &A, const R130B2Sm1<T> &B) {
    Rotor<T> res;
    res[0]=-A[3]*B[0] - A[4]*B[1] - A[5]*B[2];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    res[4]=-A[4]*B[2] + A[5]*B[1];
    res[5]=A[3]*B[2] - A[5]*B[0];
    res[6]=-A[3]*B[1] + A[4]*B[0];
    res[7]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[3]*B[0] - A[4]*B[1] - A[5]*B[2];
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
    res[15]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2<T> &A, const R130B2Sm1<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    res[3]=-A[4]*B[2] + A[5]*B[1];
    res[4]=A[3]*B[2] - A[5]*B[0];
    res[5]=-A[3]*B[1] + A[4]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> R130B2<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[0]*x0;
    T x2 = (*this)[4]*A[1];
    T x3 = (*this)[5]*A[2];
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[1]*A[1];
    T x6 = (*this)[2]*A[2];
    T x7 = 2.0*A[0];
    T x8 = (*this)[1]*(*this)[4];
    T x9 = (*this)[2]*(*this)[5];
    T x10 = A[0]*x4;
    T x11 = 2.0*(*this)[1];
    T x12 = 2.0*(*this)[4];
    T x13 = (*this)[3]*x0;
    T x14 = 2.0*(*this)[5];
    T x15 = 2.0*(*this)[2];
    T x16 = std::pow((*this)[1], 2);
    T x17 = std::pow((*this)[2], 2);
    T x18 = std::pow((*this)[3], 2);
    T x19 = std::pow((*this)[0], 2);
    T x20 = std::pow((*this)[4], 2);
    T x21 = std::pow((*this)[5], 2);
    res[0]=(*this)[3]*x1 + x0*x2 + x0*x3 + x4*x5 + x4*x6 - x7*x8 - x7*x9;
    res[1]=(*this)[1]*x10 + (*this)[4]*x1 - A[1]*x13 - 2.0*A[1]*x9 + x11*x2 + x11*x3 + x12*x6;
    res[2]=(*this)[2]*x10 + (*this)[5]*x1 - A[2]*x13 - 2.0*A[2]*x8 + x14*x5 + x15*x2 + x15*x3;
    res[3]=A[0]*x16 + A[0]*x17 + A[0]*x18 - A[0]*x19 - A[0]*x20 - A[0]*x21 - x0*x5 - x0*x6 + x2*x4 + x3*x4;
    res[4]=-(*this)[1]*x1 + (*this)[4]*x10 - A[1]*x16 + A[1]*x17 - A[1]*x18 + A[1]*x19 + A[1]*x20 - A[1]*x21 - x11*x6 + x12*x3;
    res[5]=-(*this)[2]*x1 + (*this)[5]*x10 + A[2]*x16 - A[2]*x17 - A[2]*x18 + A[2]*x19 - A[2]*x20 + A[2]*x21 + x14*x2 - x15*x5;
    return res;
};

//--------------------------------------
// (R130B2, R130B2Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130B2<T> operator+(const R130B2<T> &A, const R130B2Sp1<T> &B) {
    R130B2<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2<T> &A, const R130B2Sp1<T> &B) {
    R130B2<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2<T> &A, const R130B2Sp1<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[1]=-A[4]*B[2] + A[5]*B[1];
    res[2]=A[3]*B[2] - A[5]*B[0];
    res[3]=-A[3]*B[1] + A[4]*B[0];
    res[4]=A[1]*B[2] - A[2]*B[1];
    res[5]=-A[0]*B[2] + A[2]*B[0];
    res[6]=A[0]*B[1] - A[1]*B[0];
    res[7]=A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
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
    res[15]=A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2<T> &A, const R130B2Sp1<T> &B) {
    R130B2<T> res;
    res[0]=-A[4]*B[2] + A[5]*B[1];
    res[1]=A[3]*B[2] - A[5]*B[0];
    res[2]=-A[3]*B[1] + A[4]*B[0];
    res[3]=A[1]*B[2] - A[2]*B[1];
    res[4]=-A[0]*B[2] + A[2]*B[0];
    res[5]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> R130B2<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = std::pow((*this)[4], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[4]*x6;
    T x8 = (*this)[5]*A[2];
    T x9 = 2.0*(*this)[0];
    T x10 = (*this)[1]*A[1];
    T x11 = (*this)[2]*A[2];
    T x12 = 2.0*(*this)[4];
    T x13 = A[0]*x9;
    T x14 = 2.0*(*this)[1];
    T x15 = (*this)[5]*A[0];
    T x16 = (*this)[5]*A[1];
    T x17 = 2.0*(*this)[2];
    T x18 = (*this)[1]*A[0];
    T x19 = (*this)[0]*x6;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[0]*x4 - A[0]*x5 + A[1]*x7 - x10*x9 - x11*x9 + x6*x8;
    res[1]=-(*this)[1]*x13 + A[0]*x7 - A[1]*x0 + A[1]*x1 - A[1]*x2 + A[1]*x3 + A[1]*x4 - A[1]*x5 - x11*x14 + x12*x8;
    res[2]=-(*this)[2]*x13 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - x10*x17 + x12*x16 + x15*x6;
    res[3]=-(*this)[4]*A[1]*x9 - A[0]*x19 - x10*x6 - x11*x6 + x12*x18 + x15*x17 - x8*x9;
    res[4]=-(*this)[4]*x13 + A[1]*x19 - x10*x12 - x11*x12 - x14*x8 + x16*x17 - x18*x6;
    res[5]=(*this)[1]*A[2]*x12 - (*this)[2]*A[0]*x6 - (*this)[2]*A[1]*x12 - 2.0*(*this)[5]*x10 - (*this)[5]*x13 + A[2]*x19 - x17*x8;
    return res;
};

//-----------------------------------
// (R130B2, R130B3) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[3]*B[1] - A[4]*B[2] - A[5]*B[3];
    res[2]=-A[1]*B[3] + A[2]*B[2] + A[3]*B[0];
    res[3]=A[0]*B[3] - A[2]*B[1] + A[4]*B[0];
    res[4]=-A[0]*B[2] + A[1]*B[1] + A[5]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[12]=-A[0]*B[0] - A[4]*B[3] + A[5]*B[2];
    res[13]=-A[1]*B[0] + A[3]*B[3] - A[5]*B[1];
    res[14]=-A[2]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B2<T> &A, const R130B3<T> &B) {
    R130B1<T> res;
    res[0]=-A[3]*B[1] - A[4]*B[2] - A[5]*B[3];
    res[1]=-A[1]*B[3] + A[2]*B[2] + A[3]*B[0];
    res[2]=A[0]*B[3] - A[2]*B[1] + A[4]*B[0];
    res[3]=-A[0]*B[2] + A[1]*B[1] + A[5]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B2<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=-A[0]*B[0] - A[4]*B[3] + A[5]*B[2];
    res[2]=-A[1]*B[0] + A[3]*B[3] - A[5]*B[1];
    res[3]=-A[2]*B[0] - A[3]*B[2] + A[4]*B[1];
    return res;
};

template<typename T>
inline R130B3<T> R130B2<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = std::pow((*this)[4], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = 2.0*(*this)[0];
    T x7 = A[3]*x6;
    T x8 = 2.0*(*this)[1];
    T x9 = (*this)[5]*x8;
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[3]*A[2];
    T x12 = A[2]*x6;
    T x13 = (*this)[3]*A[3];
    T x14 = (*this)[4]*x10;
    T x15 = 2.0*(*this)[4];
    T x16 = 2.0*(*this)[5];
    T x17 = A[1]*x6;
    T x18 = A[0]*x6;
    T x19 = (*this)[2]*x8;
    T x20 = (*this)[3]*A[1];
    T x21 = (*this)[5]*x15;
    T x22 = (*this)[3]*A[0];
    res[0]=(*this)[4]*x7 - (*this)[5]*x12 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 - A[1]*x14 + A[1]*x9 + x10*x11 - x13*x8;
    res[1]=(*this)[1]*x12 + (*this)[2]*x7 + A[0]*x14 - A[0]*x9 + A[1]*x0 - A[1]*x1 - A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + x11*x15 + x13*x16;
    res[2]=(*this)[1]*x17 + (*this)[5]*x18 - A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[3]*x19 + A[3]*x21 - x10*x22 + x15*x20;
    res[3]=(*this)[2]*x17 - (*this)[4]*x18 + A[2]*x19 + A[2]*x21 - A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 + x16*x20 + x22*x8;
    return res;
};

//--------------------------------------
// (R130B2, R130B3Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[3]*B[0] - A[4]*B[1] - A[5]*B[2];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=A[0]*B[2] - A[2]*B[0];
    res[4]=-A[0]*B[1] + A[1]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[12]=-A[4]*B[2] + A[5]*B[1];
    res[13]=A[3]*B[2] - A[5]*B[0];
    res[14]=-A[3]*B[1] + A[4]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B2<T> &A, const R130B3Sm1<T> &B) {
    R130B1<T> res;
    res[0]=-A[3]*B[0] - A[4]*B[1] - A[5]*B[2];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B2<T> &A, const R130B3Sm1<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[1]=-A[4]*B[2] + A[5]*B[1];
    res[2]=A[3]*B[2] - A[5]*B[0];
    res[3]=-A[3]*B[1] + A[4]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B2<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[2]*x0;
    T x2 = 2.0*(*this)[1];
    T x3 = (*this)[5]*A[0];
    T x4 = 2.0*(*this)[2];
    T x5 = (*this)[3]*A[1];
    T x6 = A[1]*x0;
    T x7 = (*this)[3]*A[2];
    T x8 = (*this)[4]*A[0];
    T x9 = std::pow((*this)[0], 2);
    T x10 = std::pow((*this)[3], 2);
    T x11 = std::pow((*this)[1], 2);
    T x12 = std::pow((*this)[2], 2);
    T x13 = std::pow((*this)[4], 2);
    T x14 = std::pow((*this)[5], 2);
    T x15 = 2.0*(*this)[4];
    T x16 = A[0]*x0;
    T x17 = (*this)[2]*x2;
    T x18 = 2.0*(*this)[3];
    T x19 = (*this)[5]*x15;
    res[0]=(*this)[4]*x1 - (*this)[5]*x6 + x2*x3 - x2*x7 + x4*x5 - x4*x8;
    res[1]=(*this)[1]*x6 + (*this)[2]*x1 + 2.0*(*this)[5]*x7 + A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 + A[0]*x9 + x15*x5;
    res[2]=(*this)[1]*x16 - A[1]*x10 + A[1]*x11 - A[1]*x12 + A[1]*x13 - A[1]*x14 - A[1]*x9 + A[2]*x17 + A[2]*x19 + x18*x8;
    res[3]=(*this)[2]*x16 + A[1]*x17 + A[1]*x19 - A[2]*x10 - A[2]*x11 + A[2]*x12 - A[2]*x13 + A[2]*x14 - A[2]*x9 + x18*x3;
    return res;
};

//--------------------------------------
// (R130B2, R130B3Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[3]*B[0];
    res[3]=A[4]*B[0];
    res[4]=A[5]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[0]*B[0];
    res[13]=-A[1]*B[0];
    res[14]=-A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B2<T> &A, const R130B3Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[3]*B[0];
    res[1]=A[4]*B[0];
    res[2]=A[5]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B2<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B2<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2));
    res[1]=x0*(-(*this)[1]*(*this)[5] + (*this)[2]*(*this)[4]);
    res[2]=x0*((*this)[0]*(*this)[5] - (*this)[2]*(*this)[3]);
    res[3]=x0*(-(*this)[0]*(*this)[4] + (*this)[1]*(*this)[3]);
    return res;
};

//-----------------------------------
// (R130B2, R130B4) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3];
    res[9]=A[4];
    res[10]=A[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B2<T> operator*(const R130B2<T> &A, const R130B4<T> &B) {
    R130B2<T> res;
    res[0]=-A[3]*B[0];
    res[1]=-A[4]*B[0];
    res[2]=-A[5]*B[0];
    res[3]=A[0]*B[0];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B2<T> &A, const R130B4<T> &B) {
    R130B2<T> res;
    res[0]=-A[3]*B[0];
    res[1]=-A[4]*B[0];
    res[2]=-A[5]*B[0];
    res[3]=A[0]*B[0];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B2<T>::conjugate(const R130B4<T> &A) const {
    R130MV<T> res;
    res[0]=2.0*A[0]*((*this)[0]*(*this)[3] + (*this)[1]*(*this)[4] + (*this)[2]*(*this)[5]);
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
    res[15]=A[0]*(-std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2));
    return res;
};

//-----------------------------------
// (R130B2, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B2<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[0];
    res[6]=(*this)[1];
    res[7]=(*this)[2];
    res[8]=(*this)[3];
    res[9]=(*this)[4];
    res[10]=(*this)[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=A[0] + B[5];
    res[6]=A[1] + B[6];
    res[7]=A[2] + B[7];
    res[8]=A[3] + B[8];
    res[9]=A[4] + B[9];
    res[10]=A[5] + B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=A[0] - B[5];
    res[6]=A[1] - B[6];
    res[7]=A[2] - B[7];
    res[8]=A[3] - B[8];
    res[9]=A[4] - B[9];
    res[10]=A[5] - B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7] - A[3]*B[8] - A[4]*B[9] - A[5]*B[10];
    res[1]=A[0]*B[2] + A[1]*B[3] + A[2]*B[4] - A[3]*B[12] - A[4]*B[13] - A[5]*B[14];
    res[2]=A[0]*B[1] - A[1]*B[14] + A[2]*B[13] + A[3]*B[11] - A[4]*B[4] + A[5]*B[3];
    res[3]=A[0]*B[14] + A[1]*B[1] - A[2]*B[12] + A[3]*B[4] + A[4]*B[11] - A[5]*B[2];
    res[4]=-A[0]*B[13] + A[1]*B[12] + A[2]*B[1] - A[3]*B[3] + A[4]*B[2] + A[5]*B[11];
    res[5]=A[0]*B[0] - A[1]*B[10] + A[2]*B[9] - A[3]*B[15] - A[4]*B[7] + A[5]*B[6];
    res[6]=A[0]*B[10] + A[1]*B[0] - A[2]*B[8] + A[3]*B[7] - A[4]*B[15] - A[5]*B[5];
    res[7]=-A[0]*B[9] + A[1]*B[8] + A[2]*B[0] - A[3]*B[6] + A[4]*B[5] - A[5]*B[15];
    res[8]=A[0]*B[15] + A[1]*B[7] - A[2]*B[6] + A[3]*B[0] - A[4]*B[10] + A[5]*B[9];
    res[9]=-A[0]*B[7] + A[1]*B[15] + A[2]*B[5] + A[3]*B[10] + A[4]*B[0] - A[5]*B[8];
    res[10]=A[0]*B[6] - A[1]*B[5] + A[2]*B[15] - A[3]*B[9] + A[4]*B[8] + A[5]*B[0];
    res[11]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14] - A[3]*B[2] - A[4]*B[3] - A[5]*B[4];
    res[12]=-A[0]*B[11] + A[1]*B[4] - A[2]*B[3] + A[3]*B[1] - A[4]*B[14] + A[5]*B[13];
    res[13]=-A[0]*B[4] - A[1]*B[11] + A[2]*B[2] + A[3]*B[14] + A[4]*B[1] - A[5]*B[12];
    res[14]=A[0]*B[3] - A[1]*B[2] - A[2]*B[11] - A[3]*B[13] + A[4]*B[12] + A[5]*B[1];
    res[15]=A[0]*B[8] + A[1]*B[9] + A[2]*B[10] + A[3]*B[5] + A[4]*B[6] + A[5]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7] - A[3]*B[8] - A[4]*B[9] - A[5]*B[10];
    res[1]=-A[3]*B[12] - A[4]*B[13] - A[5]*B[14];
    res[2]=-A[1]*B[14] + A[2]*B[13] + A[3]*B[11];
    res[3]=A[0]*B[14] - A[2]*B[12] + A[4]*B[11];
    res[4]=-A[0]*B[13] + A[1]*B[12] + A[5]*B[11];
    res[5]=A[0]*B[0] - A[3]*B[15];
    res[6]=A[1]*B[0] - A[4]*B[15];
    res[7]=A[2]*B[0] - A[5]*B[15];
    res[8]=A[0]*B[15] + A[3]*B[0];
    res[9]=A[1]*B[15] + A[4]*B[0];
    res[10]=A[2]*B[15] + A[5]*B[0];
    res[11]=-A[3]*B[2] - A[4]*B[3] - A[5]*B[4];
    res[12]=A[1]*B[4] - A[2]*B[3] + A[3]*B[1];
    res[13]=-A[0]*B[4] + A[2]*B[2] + A[4]*B[1];
    res[14]=A[0]*B[3] - A[1]*B[2] + A[5]*B[1];
    res[15]=A[0]*B[8] + A[1]*B[9] + A[2]*B[10] + A[3]*B[5] + A[4]*B[6] + A[5]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[2] + A[1]*B[3] + A[2]*B[4];
    res[2]=A[0]*B[1] - A[4]*B[4] + A[5]*B[3];
    res[3]=A[1]*B[1] + A[3]*B[4] - A[5]*B[2];
    res[4]=A[2]*B[1] - A[3]*B[3] + A[4]*B[2];
    res[5]=-A[1]*B[10] + A[2]*B[9] - A[4]*B[7] + A[5]*B[6];
    res[6]=A[0]*B[10] - A[2]*B[8] + A[3]*B[7] - A[5]*B[5];
    res[7]=-A[0]*B[9] + A[1]*B[8] - A[3]*B[6] + A[4]*B[5];
    res[8]=A[1]*B[7] - A[2]*B[6] - A[4]*B[10] + A[5]*B[9];
    res[9]=-A[0]*B[7] + A[2]*B[5] + A[3]*B[10] - A[5]*B[8];
    res[10]=A[0]*B[6] - A[1]*B[5] - A[3]*B[9] + A[4]*B[8];
    res[11]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[12]=-A[0]*B[11] - A[4]*B[14] + A[5]*B[13];
    res[13]=-A[1]*B[11] + A[3]*B[14] - A[5]*B[12];
    res[14]=-A[2]*B[11] - A[3]*B[13] + A[4]*B[12];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B2<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[3], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*A[15];
    T x7 = (*this)[0]*(*this)[3];
    T x8 = (*this)[1]*(*this)[4];
    T x9 = (*this)[2]*(*this)[5];
    T x10 = 2.0*(*this)[0];
    T x11 = A[3]*x10;
    T x12 = 2.0*(*this)[1];
    T x13 = (*this)[3]*A[4];
    T x14 = 2.0*(*this)[2];
    T x15 = (*this)[4]*x14;
    T x16 = A[4]*x10;
    T x17 = (*this)[5]*x12;
    T x18 = (*this)[3]*A[3];
    T x19 = 2.0*(*this)[4];
    T x20 = 2.0*(*this)[5];
    T x21 = A[2]*x10;
    T x22 = (*this)[2]*x12;
    T x23 = (*this)[3]*A[1];
    T x24 = (*this)[3]*A[2];
    T x25 = (*this)[5]*x19;
    T x26 = A[1]*x10;
    T x27 = 2.0*A[8];
    T x28 = (*this)[4]*x10;
    T x29 = (*this)[5]*x10;
    T x30 = (*this)[3]*x12;
    T x31 = (*this)[3]*x14;
    T x32 = (*this)[3]*x19;
    T x33 = (*this)[3]*x20;
    T x34 = (*this)[1]*x10;
    T x35 = (*this)[2]*x10;
    T x36 = 2.0*A[9];
    T x37 = 2.0*A[10];
    T x38 = 2.0*A[5];
    T x39 = 2.0*A[6];
    T x40 = 2.0*A[7];
    T x41 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[0]*x4 - A[0]*x5 + x6*x7 + x6*x8 + x6*x9;
    res[1]=-(*this)[4]*x16 + (*this)[5]*x11 + A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 + A[1]*x4 + A[1]*x5 + A[2]*x15 - A[2]*x17 + x12*x13 - x14*x18;
    res[2]=(*this)[1]*x11 + (*this)[2]*x16 - A[1]*x15 + A[1]*x17 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 - A[2]*x5 + x13*x20 + x18*x19;
    res[3]=(*this)[1]*x21 - (*this)[5]*x26 - A[3]*x0 + A[3]*x1 - A[3]*x2 - A[3]*x3 + A[3]*x4 - A[3]*x5 + A[4]*x22 + A[4]*x25 + x14*x23 + x19*x24;
    res[4]=(*this)[2]*x21 + (*this)[4]*x26 + A[3]*x22 + A[3]*x25 - A[4]*x0 - A[4]*x1 + A[4]*x2 - A[4]*x3 - A[4]*x4 + A[4]*x5 - x12*x23 + x20*x24;
    res[5]=A[10]*x29 + A[10]*x31 + A[5]*x0 - A[5]*x1 - A[5]*x2 - A[5]*x3 + A[5]*x4 + A[5]*x5 + A[6]*x32 - A[6]*x34 + A[7]*x33 - A[7]*x35 + A[9]*x28 + A[9]*x30 + x27*x7 - x27*x8 - x27*x9;
    res[6]=A[10]*x15 + A[10]*x17 + A[5]*x32 - A[5]*x34 - A[6]*x0 + A[6]*x1 - A[6]*x2 + A[6]*x3 - A[6]*x4 + A[6]*x5 - A[7]*x22 + A[7]*x25 + A[8]*x28 + A[8]*x30 - x36*x7 + x36*x8 - x36*x9;
    res[7]=A[5]*x33 - A[5]*x35 - A[6]*x22 + A[6]*x25 - A[7]*x0 - A[7]*x1 + A[7]*x2 + A[7]*x3 + A[7]*x4 - A[7]*x5 + A[8]*x29 + A[8]*x31 + A[9]*x15 + A[9]*x17 - x37*x7 - x37*x8 + x37*x9;
    res[8]=A[10]*x33 - A[10]*x35 - A[6]*x28 - A[6]*x30 - A[7]*x29 - A[7]*x31 + A[8]*x0 - A[8]*x1 - A[8]*x2 - A[8]*x3 + A[8]*x4 + A[8]*x5 + A[9]*x32 - A[9]*x34 - x38*x7 + x38*x8 + x38*x9;
    res[9]=-A[10]*x22 + A[10]*x25 - A[5]*x28 - A[5]*x30 - A[7]*x15 - A[7]*x17 + A[8]*x32 - A[8]*x34 - A[9]*x0 + A[9]*x1 - A[9]*x2 + A[9]*x3 - A[9]*x4 + A[9]*x5 + x39*x7 - x39*x8 + x39*x9;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 + A[10]*x3 + A[10]*x4 - A[10]*x5 - A[5]*x29 - A[5]*x31 - A[6]*x15 - A[6]*x17 + A[8]*x33 - A[8]*x35 - A[9]*x22 + A[9]*x25 + x40*x7 + x40*x8 - x40*x9;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[11]*x4 + A[11]*x5 - A[12]*x15 + A[12]*x17 - A[13]*x29 + A[13]*x31 + A[14]*x28 - A[14]*x30;
    res[12]=A[11]*x15 - A[11]*x17 + A[12]*x0 - A[12]*x1 - A[12]*x2 + A[12]*x3 - A[12]*x4 - A[12]*x5 + A[13]*x32 + A[13]*x34 + A[14]*x33 + A[14]*x35;
    res[13]=A[11]*x29 - A[11]*x31 + A[12]*x32 + A[12]*x34 - A[13]*x0 + A[13]*x1 - A[13]*x2 - A[13]*x3 + A[13]*x4 - A[13]*x5 + A[14]*x22 + A[14]*x25;
    res[14]=-A[11]*x28 + A[11]*x30 + A[12]*x33 + A[12]*x35 + A[13]*x22 + A[13]*x25 - A[14]*x0 - A[14]*x1 + A[14]*x2 - A[14]*x3 - A[14]*x4 + A[14]*x5;
    res[15]=A[15]*x0 + A[15]*x1 + A[15]*x2 - A[15]*x3 - A[15]*x4 - A[15]*x5 - x41*x7 - x41*x8 - x41*x9;
    return res;
};

//-------------------------------------
// (R130B2, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3] + B[1];
    res[9]=A[4] + B[2];
    res[10]=A[5] + B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=A[3] - B[1];
    res[9]=A[4] - B[2];
    res[10]=A[5] - B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=-A[3]*B[1] - A[4]*B[2] - A[5]*B[3];
    res[1]=A[0]*B[0] - A[1]*B[3] + A[2]*B[2];
    res[2]=A[0]*B[3] + A[1]*B[0] - A[2]*B[1];
    res[3]=-A[0]*B[2] + A[1]*B[1] + A[2]*B[0];
    res[4]=A[3]*B[0] - A[4]*B[3] + A[5]*B[2];
    res[5]=A[3]*B[3] + A[4]*B[0] - A[5]*B[1];
    res[6]=-A[3]*B[2] + A[4]*B[1] + A[5]*B[0];
    res[7]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B2<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=-A[3]*B[1] - A[4]*B[2] - A[5]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=A[4]*B[0];
    res[6]=A[5]*B[0];
    res[7]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2<T> &A, const Rotation<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    res[3]=-A[4]*B[3] + A[5]*B[2];
    res[4]=A[3]*B[3] - A[5]*B[1];
    res[5]=-A[3]*B[2] + A[4]*B[1];
    return res;
};

template<typename T>
inline Rotor<T> R130B2<T>::conjugate(const Rotation<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[3], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[0];
    T x7 = A[1]*x6;
    T x8 = (*this)[4]*A[2];
    T x9 = (*this)[5]*A[3];
    T x10 = 2.0*(*this)[3];
    T x11 = (*this)[1]*A[2];
    T x12 = (*this)[2]*A[3];
    T x13 = 2.0*A[1];
    T x14 = (*this)[1]*(*this)[4];
    T x15 = (*this)[2]*(*this)[5];
    T x16 = A[1]*x10;
    T x17 = 2.0*(*this)[1];
    T x18 = 2.0*(*this)[4];
    T x19 = (*this)[3]*x6;
    T x20 = 2.0*x15;
    T x21 = 2.0*(*this)[5];
    T x22 = 2.0*(*this)[2];
    T x23 = 2.0*x14;
    res[0]=A[0]*(x0 + x1 + x2 - x3 - x4 - x5);
    res[1]=(*this)[3]*x7 + x10*x11 + x10*x12 - x13*x14 - x13*x15 + x6*x8 + x6*x9;
    res[2]=(*this)[1]*x16 + (*this)[4]*x7 - A[2]*x19 - A[2]*x20 + x12*x18 + x17*x8 + x17*x9;
    res[3]=(*this)[2]*x16 + (*this)[5]*x7 - A[3]*x19 - A[3]*x23 + x11*x21 + x22*x8 + x22*x9;
    res[4]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + x10*x8 + x10*x9 - x11*x6 - x12*x6;
    res[5]=-(*this)[1]*x7 + (*this)[4]*x16 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - x12*x17 + x18*x9;
    res[6]=-(*this)[2]*x7 + (*this)[5]*x16 - A[3]*x0 - A[3]*x1 + A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 - x11*x22 + x21*x8;
    res[7]=-A[0]*(x19 + x20 + x23);
    return res;
};

//----------------------------------
// (R130B2, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130B2<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    res[4]=(*this)[3];
    res[5]=(*this)[4];
    res[6]=(*this)[5];
    res[7]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator+(const R130B2<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    res[4]=A[3] + B[4];
    res[5]=A[4] + B[5];
    res[6]=A[5] + B[6];
    res[7]=B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const R130B2<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    res[4]=A[3] - B[4];
    res[5]=A[4] - B[5];
    res[6]=A[5] - B[6];
    res[7]=-B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3] - A[3]*B[4] - A[4]*B[5] - A[5]*B[6];
    res[1]=A[0]*B[0] - A[1]*B[6] + A[2]*B[5] - A[3]*B[7] - A[4]*B[3] + A[5]*B[2];
    res[2]=A[0]*B[6] + A[1]*B[0] - A[2]*B[4] + A[3]*B[3] - A[4]*B[7] - A[5]*B[1];
    res[3]=-A[0]*B[5] + A[1]*B[4] + A[2]*B[0] - A[3]*B[2] + A[4]*B[1] - A[5]*B[7];
    res[4]=A[0]*B[7] + A[1]*B[3] - A[2]*B[2] + A[3]*B[0] - A[4]*B[6] + A[5]*B[5];
    res[5]=-A[0]*B[3] + A[1]*B[7] + A[2]*B[1] + A[3]*B[6] + A[4]*B[0] - A[5]*B[4];
    res[6]=A[0]*B[2] - A[1]*B[1] + A[2]*B[7] - A[3]*B[5] + A[4]*B[4] + A[5]*B[0];
    res[7]=A[0]*B[4] + A[1]*B[5] + A[2]*B[6] + A[3]*B[1] + A[4]*B[2] + A[5]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B2<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3] - A[3]*B[4] - A[4]*B[5] - A[5]*B[6];
    res[1]=A[0]*B[0] - A[3]*B[7];
    res[2]=A[1]*B[0] - A[4]*B[7];
    res[3]=A[2]*B[0] - A[5]*B[7];
    res[4]=A[0]*B[7] + A[3]*B[0];
    res[5]=A[1]*B[7] + A[4]*B[0];
    res[6]=A[2]*B[7] + A[5]*B[0];
    res[7]=A[0]*B[4] + A[1]*B[5] + A[2]*B[6] + A[3]*B[1] + A[4]*B[2] + A[5]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2<T> &A, const Rotor<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[6] + A[2]*B[5] - A[4]*B[3] + A[5]*B[2];
    res[1]=A[0]*B[6] - A[2]*B[4] + A[3]*B[3] - A[5]*B[1];
    res[2]=-A[0]*B[5] + A[1]*B[4] - A[3]*B[2] + A[4]*B[1];
    res[3]=A[1]*B[3] - A[2]*B[2] - A[4]*B[6] + A[5]*B[5];
    res[4]=-A[0]*B[3] + A[2]*B[1] + A[3]*B[6] - A[5]*B[4];
    res[5]=A[0]*B[2] - A[1]*B[1] - A[3]*B[5] + A[4]*B[4];
    return res;
};

template<typename T>
inline Rotor<T> R130B2<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[3], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[0], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*A[7];
    T x7 = (*this)[0]*(*this)[3];
    T x8 = (*this)[1]*(*this)[4];
    T x9 = (*this)[2]*(*this)[5];
    T x10 = 2.0*A[4];
    T x11 = 2.0*(*this)[0];
    T x12 = (*this)[4]*A[5];
    T x13 = (*this)[5]*A[6];
    T x14 = 2.0*(*this)[3];
    T x15 = (*this)[1]*A[5];
    T x16 = (*this)[2]*A[6];
    T x17 = (*this)[4]*x14;
    T x18 = (*this)[5]*A[3];
    T x19 = (*this)[1]*x11;
    T x20 = (*this)[2]*A[3];
    T x21 = (*this)[0]*x10;
    T x22 = (*this)[3]*x10;
    T x23 = 2.0*A[5];
    T x24 = 2.0*(*this)[1];
    T x25 = 2.0*(*this)[4];
    T x26 = 2.0*(*this)[5];
    T x27 = 2.0*(*this)[2];
    T x28 = 2.0*A[6];
    T x29 = A[1]*x14;
    T x30 = (*this)[5]*A[2];
    T x31 = A[1]*x11;
    T x32 = (*this)[2]*A[2];
    T x33 = 2.0*A[1];
    T x34 = 2.0*A[2];
    T x35 = 2.0*A[3];
    T x36 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[0]*x4 - A[0]*x5 + x6*x7 + x6*x8 + x6*x9;
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[2]*x17 - A[2]*x19 + x10*x7 - x10*x8 - x10*x9 + x11*x12 + x11*x13 - x11*x20 + x14*x15 + x14*x16 + x14*x18;
    res[2]=(*this)[1]*x22 + (*this)[4]*x21 + A[1]*x17 - A[1]*x19 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 + x13*x24 + x16*x25 + x18*x25 - x20*x24 - x23*x7 + x23*x8 - x23*x9;
    res[3]=(*this)[2]*x22 - (*this)[2]*x31 + (*this)[5]*x21 + (*this)[5]*x29 - A[3]*x0 - A[3]*x1 + A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 + x12*x27 + x15*x26 - x24*x32 + x25*x30 - x28*x7 - x28*x8 + x28*x9;
    res[4]=-(*this)[1]*A[2]*x14 - (*this)[4]*A[2]*x11 + A[4]*x0 - A[4]*x1 - A[4]*x2 - A[4]*x3 + A[4]*x4 + A[4]*x5 - x11*x15 - x11*x16 - x11*x18 + x12*x14 + x13*x14 - x14*x20 - x33*x7 + x33*x8 + x33*x9;
    res[5]=-(*this)[1]*x21 - (*this)[1]*x29 + (*this)[4]*x22 - (*this)[4]*x31 - A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x3 - A[5]*x4 + A[5]*x5 + x13*x25 - x16*x24 - x18*x24 - x20*x25 + x34*x7 - x34*x8 + x34*x9;
    res[6]=-(*this)[2]*x21 - (*this)[2]*x29 + (*this)[5]*x22 - (*this)[5]*x31 - A[6]*x0 - A[6]*x1 + A[6]*x2 + A[6]*x3 + A[6]*x4 - A[6]*x5 + x12*x26 - x15*x27 - x24*x30 - x25*x32 + x35*x7 + x35*x8 - x35*x9;
    res[7]=A[7]*x0 + A[7]*x1 + A[7]*x2 - A[7]*x3 - A[7]*x4 - A[7]*x5 - x36*x7 - x36*x8 - x36*x9;
    return res;
};

//-----------------------------------
// (R130B2, scalar) binary operations
//-----------------------------------

template<typename T>
inline R130B2<T> operator*(const R130B2<T> &A, const T &B) {
    R130B2<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    res[4]=A[4]*B;
    res[5]=A[5]*B;
    return res;
};
template<typename T>
inline R130B2<T> operator*(const T &A, const R130B2<T> &B) {
    return B*A;
};
template<typename T>
inline R130B2<T> operator/(const R130B2<T> &A, const T &B) {
    return A * (1.0 / B);
};

//-------------------------------------
// (R130B2Sm1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[3] + A[2]*B[2];
    res[6]=A[0]*B[3] - A[2]*B[1];
    res[7]=-A[0]*B[2] + A[1]*B[1];
    res[8]=A[0]*B[0];
    res[9]=A[1]*B[0];
    res[10]=A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[0];
    res[9]=A[1]*B[0];
    res[10]=A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B2Sm1<T> &A, const Boost<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline Boost<T> R130B2Sm1<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//--------------------------------------
// (R130B2Sm1, R130B0) binary operations
//--------------------------------------

template<typename T>
inline Rotation<T> operator+(const R130B2Sm1<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const R130B2Sm1<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B2Sm1<T> &A, const R130B0<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator/(const R130B2Sm1<T> &A, const R130B0<T> &B) {
    R130B2Sm1<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B2Sm1<T> &A, const R130B0<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B2Sm1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sm1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[1]*B[3] + A[2]*B[2];
    res[3]=A[0]*B[3] - A[2]*B[1];
    res[4]=-A[0]*B[2] + A[1]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[12]=A[0]*B[0];
    res[13]=A[1]*B[0];
    res[14]=A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B2Sm1<T> &A, const R130B1<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B2Sm1<T> &A, const R130B1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline R130B1<T> R130B2Sm1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sm1, R130B1Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=A[0]*B[2] - A[2]*B[0];
    res[4]=-A[0]*B[1] + A[1]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator|(const R130B2Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B2Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B2Sm1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sm1, R130B1Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B2Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B2Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B2Sm1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sm1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sm1<T>::operator R130B2<T>() const {
    R130B2<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=(*this)[0];
    res[4]=(*this)[1];
    res[5]=(*this)[2];
    return res;
};

template<typename T>
inline R130B2<T> operator+(const R130B2Sm1<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=A[0] + B[3];
    res[4]=A[1] + B[4];
    res[5]=A[2] + B[5];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2Sm1<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=A[0] - B[3];
    res[4]=A[1] - B[4];
    res[5]=A[2] - B[5];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2Sm1<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=-A[0]*B[3] - A[1]*B[4] - A[2]*B[5];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    res[4]=-A[1]*B[5] + A[2]*B[4];
    res[5]=A[0]*B[5] - A[2]*B[3];
    res[6]=-A[0]*B[4] + A[1]*B[3];
    res[7]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[3] - A[1]*B[4] - A[2]*B[5];
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
    res[15]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2Sm1<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    res[3]=-A[1]*B[5] + A[2]*B[4];
    res[4]=A[0]*B[5] - A[2]*B[3];
    res[5]=-A[0]*B[4] + A[1]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> R130B2Sm1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=A[0]*x7 + A[1]*x8 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    res[3]=A[3]*x0 - A[3]*x1 - A[3]*x2 + A[4]*x4 + A[5]*x7;
    res[4]=A[3]*x4 - A[4]*x0 + A[4]*x1 - A[4]*x2 + A[5]*x8;
    res[5]=A[3]*x7 + A[4]*x8 - A[5]*x0 - A[5]*x1 + A[5]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sm1, R130B2Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130B2Sm1<T> operator+(const R130B2Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator-(const R130B2Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const R130B2Sm1<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B2Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const R130B2Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B2Sm1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sm1, R130B2Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130B2<T> operator+(const R130B2Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B2<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=A[0];
    res[4]=A[1];
    res[5]=A[2];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B2<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=A[0];
    res[4]=A[1];
    res[5]=A[2];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[2] + A[2]*B[1];
    res[6]=A[0]*B[2] - A[2]*B[0];
    res[7]=-A[0]*B[1] + A[1]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B4<T> operator|(const R130B2Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B2Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B2Sm1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//--------------------------------------
// (R130B2Sm1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[1]*B[3] + A[2]*B[2];
    res[13]=A[0]*B[3] - A[2]*B[1];
    res[14]=-A[0]*B[2] + A[1]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B2Sm1<T> &A, const R130B3<T> &B) {
    R130B1<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B2Sm1<T> &A, const R130B3<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline R130B3<T> R130B2Sm1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sm1, R130B3Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
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
    res[12]=-A[1]*B[2] + A[2]*B[1];
    res[13]=A[0]*B[2] - A[2]*B[0];
    res[14]=-A[0]*B[1] + A[1]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator|(const R130B2Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B2Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B2Sm1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sm1, R130B3Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator*(const R130B2Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B2Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B2Sm1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sm1, R130B4) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0];
    res[9]=A[1];
    res[10]=A[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B2Sm1<T> &A, const R130B4<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator|(const R130B2Sm1<T> &A, const R130B4<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> R130B2Sm1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sm1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sm1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=(*this)[0];
    res[9]=(*this)[1];
    res[10]=(*this)[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B2Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=A[0] + B[8];
    res[9]=A[1] + B[9];
    res[10]=A[2] + B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=A[0] - B[8];
    res[9]=A[1] - B[9];
    res[10]=A[2] - B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[8] - A[1]*B[9] - A[2]*B[10];
    res[1]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[2]=A[0]*B[11] - A[1]*B[4] + A[2]*B[3];
    res[3]=A[0]*B[4] + A[1]*B[11] - A[2]*B[2];
    res[4]=-A[0]*B[3] + A[1]*B[2] + A[2]*B[11];
    res[5]=-A[0]*B[15] - A[1]*B[7] + A[2]*B[6];
    res[6]=A[0]*B[7] - A[1]*B[15] - A[2]*B[5];
    res[7]=-A[0]*B[6] + A[1]*B[5] - A[2]*B[15];
    res[8]=A[0]*B[0] - A[1]*B[10] + A[2]*B[9];
    res[9]=A[0]*B[10] + A[1]*B[0] - A[2]*B[8];
    res[10]=-A[0]*B[9] + A[1]*B[8] + A[2]*B[0];
    res[11]=-A[0]*B[2] - A[1]*B[3] - A[2]*B[4];
    res[12]=A[0]*B[1] - A[1]*B[14] + A[2]*B[13];
    res[13]=A[0]*B[14] + A[1]*B[1] - A[2]*B[12];
    res[14]=-A[0]*B[13] + A[1]*B[12] + A[2]*B[1];
    res[15]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[8] - A[1]*B[9] - A[2]*B[10];
    res[1]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[2]=A[0]*B[11];
    res[3]=A[1]*B[11];
    res[4]=A[2]*B[11];
    res[5]=-A[0]*B[15];
    res[6]=-A[1]*B[15];
    res[7]=-A[2]*B[15];
    res[8]=A[0]*B[0];
    res[9]=A[1]*B[0];
    res[10]=A[2]*B[0];
    res[11]=-A[0]*B[2] - A[1]*B[3] - A[2]*B[4];
    res[12]=A[0]*B[1];
    res[13]=A[1]*B[1];
    res[14]=A[2]*B[1];
    res[15]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[1]*B[4] + A[2]*B[3];
    res[3]=A[0]*B[4] - A[2]*B[2];
    res[4]=-A[0]*B[3] + A[1]*B[2];
    res[5]=-A[1]*B[7] + A[2]*B[6];
    res[6]=A[0]*B[7] - A[2]*B[5];
    res[7]=-A[0]*B[6] + A[1]*B[5];
    res[8]=-A[1]*B[10] + A[2]*B[9];
    res[9]=A[0]*B[10] - A[2]*B[8];
    res[10]=-A[0]*B[9] + A[1]*B[8];
    res[11]=0;
    res[12]=-A[1]*B[14] + A[2]*B[13];
    res[13]=A[0]*B[14] - A[2]*B[12];
    res[14]=-A[0]*B[13] + A[1]*B[12];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B2Sm1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[4];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*(x0 + x1 + x2);
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=A[5]*x0 - A[5]*x1 - A[5]*x2 + A[6]*x4 + A[7]*x7;
    res[6]=A[5]*x4 - A[6]*x0 + A[6]*x1 - A[6]*x2 + A[7]*x8;
    res[7]=A[5]*x7 + A[6]*x8 - A[7]*x0 - A[7]*x1 + A[7]*x2;
    res[8]=A[10]*x7 + A[8]*x0 - A[8]*x1 - A[8]*x2 + A[9]*x4;
    res[9]=A[10]*x8 + A[8]*x4 - A[9]*x0 + A[9]*x1 - A[9]*x2;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 + A[8]*x7 + A[9]*x8;
    res[11]=A[11]*(x0 + x1 + x2);
    res[12]=A[12]*x0 - A[12]*x1 - A[12]*x2 + A[13]*x4 + A[14]*x7;
    res[13]=A[12]*x4 - A[13]*x0 + A[13]*x1 - A[13]*x2 + A[14]*x8;
    res[14]=A[12]*x7 + A[13]*x8 - A[14]*x0 - A[14]*x1 + A[14]*x2;
    res[15]=A[15]*(x0 + x1 + x2);
    return res;
};

//----------------------------------------
// (R130B2Sm1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130B2Sm1<T>::operator Rotation<T>() const {
    Rotation<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    return res;
};

template<typename T>
inline Rotation<T> operator+(const R130B2Sm1<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const R130B2Sm1<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const R130B2Sm1<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=A[0]*B[0] - A[1]*B[3] + A[2]*B[2];
    res[2]=A[0]*B[3] + A[1]*B[0] - A[2]*B[1];
    res[3]=-A[0]*B[2] + A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline Rotation<T> operator|(const R130B2Sm1<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const R130B2Sm1<T> &A, const Rotation<T> &B) {
    R130B2Sm1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline Rotation<T> R130B2Sm1<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-------------------------------------
// (R130B2Sm1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130B2Sm1<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=(*this)[0];
    res[5]=(*this)[1];
    res[6]=(*this)[2];
    res[7]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator+(const R130B2Sm1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=A[0] + B[4];
    res[5]=A[1] + B[5];
    res[6]=A[2] + B[6];
    res[7]=B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const R130B2Sm1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=A[0] - B[4];
    res[5]=A[1] - B[5];
    res[6]=A[2] - B[6];
    res[7]=-B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2Sm1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-A[0]*B[4] - A[1]*B[5] - A[2]*B[6];
    res[1]=-A[0]*B[7] - A[1]*B[3] + A[2]*B[2];
    res[2]=A[0]*B[3] - A[1]*B[7] - A[2]*B[1];
    res[3]=-A[0]*B[2] + A[1]*B[1] - A[2]*B[7];
    res[4]=A[0]*B[0] - A[1]*B[6] + A[2]*B[5];
    res[5]=A[0]*B[6] + A[1]*B[0] - A[2]*B[4];
    res[6]=-A[0]*B[5] + A[1]*B[4] + A[2]*B[0];
    res[7]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B2Sm1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-A[0]*B[4] - A[1]*B[5] - A[2]*B[6];
    res[1]=-A[0]*B[7];
    res[2]=-A[1]*B[7];
    res[3]=-A[2]*B[7];
    res[4]=A[0]*B[0];
    res[5]=A[1]*B[0];
    res[6]=A[2]*B[0];
    res[7]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2Sm1<T> &A, const Rotor<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    res[3]=-A[1]*B[6] + A[2]*B[5];
    res[4]=A[0]*B[6] - A[2]*B[4];
    res[5]=-A[0]*B[5] + A[1]*B[4];
    return res;
};

template<typename T>
inline Rotor<T> R130B2Sm1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=A[1]*x7 + A[2]*x8 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    res[4]=A[4]*x0 - A[4]*x1 - A[4]*x2 + A[5]*x4 + A[6]*x7;
    res[5]=A[4]*x4 - A[5]*x0 + A[5]*x1 - A[5]*x2 + A[6]*x8;
    res[6]=A[4]*x7 + A[5]*x8 - A[6]*x0 - A[6]*x1 + A[6]*x2;
    res[7]=A[7]*(x0 + x1 + x2);
    return res;
};

//--------------------------------------
// (R130B2Sm1, scalar) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sm1<T> operator*(const R130B2Sm1<T> &A, const T &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    return res;
};
template<typename T>
inline R130B2Sm1<T> operator*(const T &A, const R130B2Sm1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B2Sm1<T> operator/(const R130B2Sm1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//-------------------------------------
// (R130B2Sp1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130B2Sp1<T>::operator Boost<T>() const {
    Boost<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    return res;
};

template<typename T>
inline Boost<T> operator+(const R130B2Sp1<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    return res;
};

template<typename T>
inline Boost<T> operator-(const R130B2Sp1<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0];
    res[6]=A[1]*B[0];
    res[7]=A[2]*B[0];
    res[8]=A[1]*B[3] - A[2]*B[2];
    res[9]=-A[0]*B[3] + A[2]*B[1];
    res[10]=A[0]*B[2] - A[1]*B[1];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Boost<T> operator|(const R130B2Sp1<T> &A, const Boost<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const R130B2Sp1<T> &A, const Boost<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[1]*B[3] - A[2]*B[2];
    res[1]=-A[0]*B[3] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[1]*B[1];
    return res;
};

template<typename T>
inline Boost<T> R130B2Sp1<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-(*this)[2]*A[1]*x3 - (*this)[2]*A[2]*x6 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    return res;
};

//--------------------------------------
// (R130B2Sp1, R130B0) binary operations
//--------------------------------------

template<typename T>
inline Boost<T> operator+(const R130B2Sp1<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline Boost<T> operator-(const R130B2Sp1<T> &A, const R130B0<T> &B) {
    Boost<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B2Sp1<T> &A, const R130B0<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator/(const R130B2Sp1<T> &A, const R130B0<T> &B) {
    R130B2Sp1<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator|(const R130B2Sp1<T> &A, const R130B0<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B2Sp1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sp1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
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

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[2]=A[0]*B[0];
    res[3]=A[1]*B[0];
    res[4]=A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[1]*B[3] - A[2]*B[2];
    res[13]=-A[0]*B[3] + A[2]*B[1];
    res[14]=A[0]*B[2] - A[1]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B2Sp1<T> &A, const R130B1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[1]*B[3] - A[2]*B[2];
    res[1]=-A[0]*B[3] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[1]*B[1];
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B2Sp1<T> &A, const R130B1<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B2Sp1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sp1, R130B1Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
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

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
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
    res[12]=A[1]*B[2] - A[2]*B[1];
    res[13]=-A[0]*B[2] + A[2]*B[0];
    res[14]=A[0]*B[1] - A[1]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B2Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[1]*B[2] - A[2]*B[1];
    res[1]=-A[0]*B[2] + A[2]*B[0];
    res[2]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator^(const R130B2Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B2Sp1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sp1, R130B1Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
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

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
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

template<typename T>
inline R130B1Sm1<T> operator*(const R130B2Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B2Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B2Sp1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sp1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sp1<T>::operator R130B2<T>() const {
    R130B2<T> res;
    res[0]=(*this)[0];
    res[1]=(*this)[1];
    res[2]=(*this)[2];
    res[3]=0;
    res[4]=0;
    res[5]=0;
    return res;
};

template<typename T>
inline R130B2<T> operator+(const R130B2Sp1<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2Sp1<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2Sp1<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[1]=-A[1]*B[5] + A[2]*B[4];
    res[2]=A[0]*B[5] - A[2]*B[3];
    res[3]=-A[0]*B[4] + A[1]*B[3];
    res[4]=A[1]*B[2] - A[2]*B[1];
    res[5]=-A[0]*B[2] + A[2]*B[0];
    res[6]=A[0]*B[1] - A[1]*B[0];
    res[7]=A[0]*B[3] + A[1]*B[4] + A[2]*B[5];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
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
    res[15]=A[0]*B[3] + A[1]*B[4] + A[2]*B[5];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2Sp1<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[5] + A[2]*B[4];
    res[1]=A[0]*B[5] - A[2]*B[3];
    res[2]=-A[0]*B[4] + A[1]*B[3];
    res[3]=A[1]*B[2] - A[2]*B[1];
    res[4]=-A[0]*B[2] + A[2]*B[0];
    res[5]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> R130B2Sp1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-A[0]*x7 - A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    res[3]=A[3]*x0 + A[3]*x1 - A[3]*x2 - A[4]*x4 - A[5]*x7;
    res[4]=-A[3]*x4 - A[4]*x0 + A[4]*x1 + A[4]*x2 - A[5]*x8;
    res[5]=-A[3]*x7 - A[4]*x8 + A[5]*x0 - A[5]*x1 + A[5]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sp1, R130B2Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130B2<T> operator+(const R130B2Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B2<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=B[0];
    res[4]=B[1];
    res[5]=B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator-(const R130B2Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B2<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=-B[0];
    res[4]=-B[1];
    res[5]=-B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[2] + A[2]*B[1];
    res[6]=A[0]*B[2] - A[2]*B[0];
    res[7]=-A[0]*B[1] + A[1]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B4<T> operator|(const R130B2Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B2Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B2Sp1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-(*this)[2]*A[0]*x3 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sp1, R130B2Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130B2Sp1<T> operator+(const R130B2Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator-(const R130B2Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const R130B2Sp1<T> &A, const R130B2Sp1<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[1]=A[1]*B[2] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B2Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const R130B2Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[1]*B[2] - A[2]*B[1];
    res[1]=-A[0]*B[2] + A[2]*B[0];
    res[2]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B2Sp1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-(*this)[2]*A[0]*x3 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//--------------------------------------
// (R130B2Sp1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[1]*B[3] + A[2]*B[2];
    res[3]=A[0]*B[3] - A[2]*B[1];
    res[4]=-A[0]*B[2] + A[1]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[12]=-A[0]*B[0];
    res[13]=-A[1]*B[0];
    res[14]=-A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B2Sp1<T> &A, const R130B3<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B2Sp1<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=-A[0]*B[0];
    res[2]=-A[1]*B[0];
    res[3]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B2Sp1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sp1, R130B3Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=A[0]*B[2] - A[2]*B[0];
    res[4]=-A[0]*B[1] + A[1]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B2Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator^(const R130B2Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B2Sp1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B2Sp1, R130B3Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B2Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B2Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B2Sp1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sp1, R130B4) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B2Sp1<T> &A, const R130B4<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B2Sp1<T> &A, const R130B4<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> R130B2Sp1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B2Sp1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sp1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[0];
    res[6]=(*this)[1];
    res[7]=(*this)[2];
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

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=A[0] + B[5];
    res[6]=A[1] + B[6];
    res[7]=A[2] + B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=A[0] - B[5];
    res[6]=A[1] - B[6];
    res[7]=A[2] - B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7];
    res[1]=A[0]*B[2] + A[1]*B[3] + A[2]*B[4];
    res[2]=A[0]*B[1] - A[1]*B[14] + A[2]*B[13];
    res[3]=A[0]*B[14] + A[1]*B[1] - A[2]*B[12];
    res[4]=-A[0]*B[13] + A[1]*B[12] + A[2]*B[1];
    res[5]=A[0]*B[0] - A[1]*B[10] + A[2]*B[9];
    res[6]=A[0]*B[10] + A[1]*B[0] - A[2]*B[8];
    res[7]=-A[0]*B[9] + A[1]*B[8] + A[2]*B[0];
    res[8]=A[0]*B[15] + A[1]*B[7] - A[2]*B[6];
    res[9]=-A[0]*B[7] + A[1]*B[15] + A[2]*B[5];
    res[10]=A[0]*B[6] - A[1]*B[5] + A[2]*B[15];
    res[11]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[12]=-A[0]*B[11] + A[1]*B[4] - A[2]*B[3];
    res[13]=-A[0]*B[4] - A[1]*B[11] + A[2]*B[2];
    res[14]=A[0]*B[3] - A[1]*B[2] - A[2]*B[11];
    res[15]=A[0]*B[8] + A[1]*B[9] + A[2]*B[10];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7];
    res[1]=0;
    res[2]=-A[1]*B[14] + A[2]*B[13];
    res[3]=A[0]*B[14] - A[2]*B[12];
    res[4]=-A[0]*B[13] + A[1]*B[12];
    res[5]=A[0]*B[0];
    res[6]=A[1]*B[0];
    res[7]=A[2]*B[0];
    res[8]=A[0]*B[15];
    res[9]=A[1]*B[15];
    res[10]=A[2]*B[15];
    res[11]=0;
    res[12]=A[1]*B[4] - A[2]*B[3];
    res[13]=-A[0]*B[4] + A[2]*B[2];
    res[14]=A[0]*B[3] - A[1]*B[2];
    res[15]=A[0]*B[8] + A[1]*B[9] + A[2]*B[10];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B2Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[2] + A[1]*B[3] + A[2]*B[4];
    res[2]=A[0]*B[1];
    res[3]=A[1]*B[1];
    res[4]=A[2]*B[1];
    res[5]=-A[1]*B[10] + A[2]*B[9];
    res[6]=A[0]*B[10] - A[2]*B[8];
    res[7]=-A[0]*B[9] + A[1]*B[8];
    res[8]=A[1]*B[7] - A[2]*B[6];
    res[9]=-A[0]*B[7] + A[2]*B[5];
    res[10]=A[0]*B[6] - A[1]*B[5];
    res[11]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[12]=-A[0]*B[11];
    res[13]=-A[1]*B[11];
    res[14]=-A[2]*B[11];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B2Sp1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[4];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=A[1]*(x0 + x1 + x2);
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=-A[5]*x0 + A[5]*x1 + A[5]*x2 - A[6]*x4 - A[7]*x7;
    res[6]=-A[5]*x4 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[7]*x8;
    res[7]=-A[5]*x7 - A[6]*x8 + A[7]*x0 + A[7]*x1 - A[7]*x2;
    res[8]=-A[10]*x7 - A[8]*x0 + A[8]*x1 + A[8]*x2 - A[9]*x4;
    res[9]=-A[10]*x8 - A[8]*x4 + A[9]*x0 - A[9]*x1 + A[9]*x2;
    res[10]=A[10]*x0 + A[10]*x1 - A[10]*x2 - A[8]*x7 - A[9]*x8;
    res[11]=A[11]*(x0 + x1 + x2);
    res[12]=A[12]*x0 - A[12]*x1 - A[12]*x2 + A[13]*x4 + A[14]*x7;
    res[13]=A[12]*x4 - A[13]*x0 + A[13]*x1 - A[13]*x2 + A[14]*x8;
    res[14]=A[12]*x7 + A[13]*x8 - A[14]*x0 - A[14]*x1 + A[14]*x2;
    res[15]=-A[15]*(x0 + x1 + x2);
    return res;
};

//----------------------------------------
// (R130B2Sp1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B2Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B2Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0];
    res[6]=A[1];
    res[7]=A[2];
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B2Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0] - A[1]*B[3] + A[2]*B[2];
    res[6]=A[0]*B[3] + A[1]*B[0] - A[2]*B[1];
    res[7]=-A[0]*B[2] + A[1]*B[1] + A[2]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B2Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0];
    res[6]=A[1]*B[0];
    res[7]=A[2]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B2Sp1<T> &A, const Rotation<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline Rotation<T> R130B2Sp1<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-(*this)[2]*A[1]*x3 - (*this)[2]*A[2]*x6 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    return res;
};

//-------------------------------------
// (R130B2Sp1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130B2Sp1<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator+(const R130B2Sp1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const R130B2Sp1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B2Sp1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0] - A[1]*B[6] + A[2]*B[5];
    res[2]=A[0]*B[6] + A[1]*B[0] - A[2]*B[4];
    res[3]=-A[0]*B[5] + A[1]*B[4] + A[2]*B[0];
    res[4]=A[0]*B[7] + A[1]*B[3] - A[2]*B[2];
    res[5]=-A[0]*B[3] + A[1]*B[7] + A[2]*B[1];
    res[6]=A[0]*B[2] - A[1]*B[1] + A[2]*B[7];
    res[7]=A[0]*B[4] + A[1]*B[5] + A[2]*B[6];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B2Sp1<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[0]*B[7];
    res[5]=A[1]*B[7];
    res[6]=A[2]*B[7];
    res[7]=A[0]*B[4] + A[1]*B[5] + A[2]*B[6];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B2Sp1<T> &A, const Rotor<T> &B) {
    R130B2<T> res;
    res[0]=-A[1]*B[6] + A[2]*B[5];
    res[1]=A[0]*B[6] - A[2]*B[4];
    res[2]=-A[0]*B[5] + A[1]*B[4];
    res[3]=A[1]*B[3] - A[2]*B[2];
    res[4]=-A[0]*B[3] + A[2]*B[1];
    res[5]=A[0]*B[2] - A[1]*B[1];
    return res;
};

template<typename T>
inline Rotor<T> R130B2Sp1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-A[1]*x7 - A[2]*x8 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    res[4]=-A[4]*x0 + A[4]*x1 + A[4]*x2 - A[5]*x4 - A[6]*x7;
    res[5]=-A[4]*x4 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[6]*x8;
    res[6]=-A[4]*x7 - A[5]*x8 + A[6]*x0 + A[6]*x1 - A[6]*x2;
    res[7]=-A[7]*(x0 + x1 + x2);
    return res;
};

//--------------------------------------
// (R130B2Sp1, scalar) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sp1<T> operator*(const R130B2Sp1<T> &A, const T &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    return res;
};
template<typename T>
inline R130B2Sp1<T> operator*(const T &A, const R130B2Sp1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B2Sp1<T> operator/(const R130B2Sp1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//----------------------------------
// (R130B3, Boost) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[2]*B[3] - A[3]*B[2];
    res[3]=-A[1]*B[3] + A[3]*B[1];
    res[4]=A[1]*B[2] - A[2]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[12]=A[0]*B[1] + A[1]*B[0];
    res[13]=A[0]*B[2] + A[2]*B[0];
    res[14]=A[0]*B[3] + A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[2]*B[3] - A[3]*B[2];
    res[3]=-A[1]*B[3] + A[3]*B[1];
    res[4]=A[1]*B[2] - A[2]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B3<T> &A, const Boost<T> &B) {
    R130B3<T> res;
    res[0]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> R130B3<T>::conjugate(const Boost<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[1];
    T x5 = (*this)[2]*x4;
    T x6 = (*this)[3]*A[3];
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[1];
    T x9 = (*this)[3]*A[2];
    T x10 = (*this)[0]*A[3];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x5 + x4*x6;
    res[6]=A[1]*x5 - A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + x6*x7;
    res[7]=-A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + x4*x8 + x7*x9;
    res[8]=x10*x7 - x11*x9;
    res[9]=-x10*x4 + x11*x8;
    res[10]=(*this)[0]*(-A[1]*x7 + A[2]*x4);
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//-----------------------------------
// (R130B3, R130B0) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
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
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
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
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator*(const R130B3<T> &A, const R130B0<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator/(const R130B3<T> &A, const R130B0<T> &B) {
    R130B3<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

template<typename T>
inline R130B3<T> operator|(const R130B3<T> &A, const R130B0<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B3<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2));
    return res;
};

//-----------------------------------
// (R130B3, R130B1) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[2]*B[3] - A[3]*B[2];
    res[6]=-A[1]*B[3] + A[3]*B[1];
    res[7]=A[1]*B[2] - A[2]*B[1];
    res[8]=A[0]*B[1] + A[1]*B[0];
    res[9]=A[0]*B[2] + A[2]*B[0];
    res[10]=A[0]*B[3] + A[3]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B3<T> &A, const R130B1<T> &B) {
    R130B2<T> res;
    res[0]=A[2]*B[3] - A[3]*B[2];
    res[1]=-A[1]*B[3] + A[3]*B[1];
    res[2]=A[1]*B[2] - A[2]*B[1];
    res[3]=A[0]*B[1] + A[1]*B[0];
    res[4]=A[0]*B[2] + A[2]*B[0];
    res[5]=A[0]*B[3] + A[3]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3<T> &A, const R130B1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B3<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[2];
    T x7 = (*this)[3]*A[3];
    T x8 = 2.0*(*this)[1];
    T x9 = A[0]*x4;
    T x10 = A[1]*x8;
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3 - A[1]*x5 - x4*x6 - x4*x7;
    res[1]=A[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + x6*x8 + x7*x8;
    res[2]=(*this)[2]*x10 + 2.0*(*this)[2]*x7 + (*this)[2]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3;
    res[3]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//--------------------------------------
// (R130B3, R130B1Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[2]*B[2] - A[3]*B[1];
    res[6]=-A[1]*B[2] + A[3]*B[0];
    res[7]=A[1]*B[1] - A[2]*B[0];
    res[8]=A[0]*B[0];
    res[9]=A[0]*B[1];
    res[10]=A[0]*B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B3<T> &A, const R130B1Sm1<T> &B) {
    R130B2<T> res;
    res[0]=A[2]*B[2] - A[3]*B[1];
    res[1]=-A[1]*B[2] + A[3]*B[0];
    res[2]=A[1]*B[1] - A[2]*B[0];
    res[3]=A[0]*B[0];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3<T> &A, const R130B1Sm1<T> &B) {
    R130B4<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    return res;
};

template<typename T>
inline R130B1<T> R130B3<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*A[0];
    T x2 = (*this)[2]*A[1];
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=-x0*(x1 + x2 + x3);
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x6 - A[0]*x7 + x2*x8 + x3*x8;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + x1*x10 + x10*x2;
    return res;
};

//--------------------------------------
// (R130B3, R130B1Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1]*B[0];
    res[9]=A[2]*B[0];
    res[10]=A[3]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B3<T> &A, const R130B1Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[1]*B[0];
    res[1]=A[2]*B[0];
    res[2]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3<T> &A, const R130B1Sp1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B3<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[0]*A[0];
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    res[1]=(*this)[1]*x0;
    res[2]=(*this)[2]*x0;
    res[3]=(*this)[3]*x0;
    return res;
};

//-----------------------------------
// (R130B3, R130B2) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[3] - A[2]*B[4] - A[3]*B[5];
    res[2]=A[0]*B[3] + A[2]*B[2] - A[3]*B[1];
    res[3]=A[0]*B[4] - A[1]*B[2] + A[3]*B[0];
    res[4]=A[0]*B[5] + A[1]*B[1] - A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[12]=A[0]*B[0] - A[2]*B[5] + A[3]*B[4];
    res[13]=A[0]*B[1] + A[1]*B[5] - A[3]*B[3];
    res[14]=A[0]*B[2] - A[1]*B[4] + A[2]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B3<T> &A, const R130B2<T> &B) {
    R130B1<T> res;
    res[0]=-A[1]*B[3] - A[2]*B[4] - A[3]*B[5];
    res[1]=A[0]*B[3] + A[2]*B[2] - A[3]*B[1];
    res[2]=A[0]*B[4] - A[1]*B[2] + A[3]*B[0];
    res[3]=A[0]*B[5] + A[1]*B[1] - A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B3<T> &A, const R130B2<T> &B) {
    R130B3<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0] - A[2]*B[5] + A[3]*B[4];
    res[2]=A[0]*B[1] + A[1]*B[5] - A[3]*B[3];
    res[3]=A[0]*B[2] - A[1]*B[4] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> R130B3<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[0], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = (*this)[0]*A[5];
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*x7;
    T x9 = (*this)[0]*x7;
    T x10 = (*this)[3]*x4;
    T x11 = 2.0*(*this)[1];
    T x12 = (*this)[0]*x11;
    T x13 = (*this)[0]*x4;
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x6 + A[2]*x8 - A[4]*x9 + x4*x5;
    res[1]=A[0]*x6 - A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 + A[2]*x10 + A[3]*x9 - x11*x5;
    res[2]=A[0]*x8 + A[1]*x10 - A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 - A[3]*x13 + A[4]*x12;
    res[3]=-A[1]*x9 + A[2]*x13 - A[3]*x0 + A[3]*x1 + A[3]*x2 + A[3]*x3 - A[4]*x6 - A[5]*x8;
    res[4]=A[0]*x9 - A[2]*x12 - A[3]*x6 + A[4]*x0 + A[4]*x1 - A[4]*x2 + A[4]*x3 - A[5]*x10;
    res[5]=-A[0]*x13 + A[1]*x12 - A[3]*x8 - A[4]*x10 + A[5]*x0 + A[5]*x1 + A[5]*x2 - A[5]*x3;
    return res;
};

//--------------------------------------
// (R130B3, R130B2Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[2]=A[0]*B[0];
    res[3]=A[0]*B[1];
    res[4]=A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[2]*B[2] + A[3]*B[1];
    res[13]=A[1]*B[2] - A[3]*B[0];
    res[14]=-A[1]*B[1] + A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B3<T> &A, const R130B2Sm1<T> &B) {
    R130B1<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3<T> &A, const R130B2Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> R130B3<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[2]*x0;
    T x2 = (*this)[3]*x0;
    T x3 = (*this)[1]*A[1];
    T x4 = (*this)[2]*A[0];
    T x5 = std::pow((*this)[0], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = std::pow((*this)[1], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*A[2];
    res[0]=(*this)[2]*x1 - A[1]*x2;
    res[1]=-(*this)[1]*x1 + A[0]*x2;
    res[2]=x0*(x3 - x4);
    res[3]=A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - x10*x11 - x3*x9;
    res[4]=A[1]*x5 - A[1]*x6 + A[1]*x7 + A[1]*x8 - x10*x4 - x11*x9;
    res[5]=-(*this)[3]*A[0]*x10 - (*this)[3]*A[1]*x9 + A[2]*x5 + A[2]*x6 - A[2]*x7 + A[2]*x8;
    return res;
};

//--------------------------------------
// (R130B3, R130B2Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[2]*B[2] - A[3]*B[1];
    res[3]=-A[1]*B[2] + A[3]*B[0];
    res[4]=A[1]*B[1] - A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B3<T> &A, const R130B2Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[2]*B[2] - A[3]*B[1];
    res[1]=-A[1]*B[2] + A[3]*B[0];
    res[2]=A[1]*B[1] - A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B3<T> &A, const R130B2Sp1<T> &B) {
    R130B3<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> R130B3<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[0], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[1];
    T x5 = (*this)[2]*x4;
    T x6 = (*this)[3]*A[2];
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[0];
    T x9 = (*this)[3]*A[1];
    T x10 = (*this)[0]*A[2];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x5 + x4*x6;
    res[1]=A[0]*x5 - A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 + x6*x7;
    res[2]=-A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 + x4*x8 + x7*x9;
    res[3]=x10*x7 - x11*x9;
    res[4]=-x10*x4 + x11*x8;
    res[5]=(*this)[0]*(-A[0]*x7 + A[1]*x4);
    return res;
};

//-----------------------------------
// (R130B3, R130B3) binary operations
//-----------------------------------

template<typename T>
inline R130B3<T> operator+(const R130B3<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[1] - A[1]*B[0];
    res[6]=A[0]*B[2] - A[2]*B[0];
    res[7]=A[0]*B[3] - A[3]*B[0];
    res[8]=-A[2]*B[3] + A[3]*B[2];
    res[9]=A[1]*B[3] - A[3]*B[1];
    res[10]=-A[1]*B[2] + A[2]*B[1];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3<T> &A, const R130B3<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B3<T> &A, const R130B3<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[1] - A[1]*B[0];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=A[0]*B[3] - A[3]*B[0];
    res[3]=-A[2]*B[3] + A[3]*B[2];
    res[4]=A[1]*B[3] - A[3]*B[1];
    res[5]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130B3<T> R130B3<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[2];
    T x7 = (*this)[3]*A[3];
    T x8 = 2.0*(*this)[1];
    T x9 = A[0]*x4;
    T x10 = A[1]*x8;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[1]*x5 - x4*x6 - x4*x7;
    res[1]=A[0]*x5 - A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - x6*x8 - x7*x8;
    res[2]=-(*this)[2]*x10 - 2.0*(*this)[2]*x7 + (*this)[2]*x9 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3;
    res[3]=-(*this)[3]*x10 - 2.0*(*this)[3]*x6 + (*this)[3]*x9 - A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3;
    return res;
};

//--------------------------------------
// (R130B3, R130B3Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130B3<T> operator+(const R130B3<T> &A, const R130B3Sm1<T> &B) {
    R130B3<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3<T> &A, const R130B3Sm1<T> &B) {
    R130B3<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0];
    res[6]=A[0]*B[1];
    res[7]=A[0]*B[2];
    res[8]=-A[2]*B[2] + A[3]*B[1];
    res[9]=A[1]*B[2] - A[3]*B[0];
    res[10]=-A[1]*B[1] + A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3<T> &A, const R130B3Sm1<T> &B) {
    R130B0<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B3<T> &A, const R130B3Sm1<T> &B) {
    R130B2<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=-A[2]*B[2] + A[3]*B[1];
    res[4]=A[1]*B[2] - A[3]*B[0];
    res[5]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B3<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*A[0];
    T x2 = (*this)[2]*A[1];
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[2], 2);
    T x5 = std::pow((*this)[3], 2);
    T x6 = std::pow((*this)[0], 2);
    T x7 = std::pow((*this)[1], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=-x0*(x1 + x2 + x3);
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x6 - A[0]*x7 - x2*x8 - x3*x8;
    res[2]=-A[1]*x4 + A[1]*x5 - A[1]*x6 + A[1]*x7 - x1*x9 - x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 - x1*x10 - x10*x2;
    return res;
};

//--------------------------------------
// (R130B3, R130B3Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130B3<T> operator+(const R130B3<T> &A, const R130B3Sp1<T> &B) {
    R130B3<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3<T> &A, const R130B3Sp1<T> &B) {
    R130B3<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline Boost<T> operator*(const R130B3<T> &A, const R130B3Sp1<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    res[3]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3<T> &A, const R130B3Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B3<T> &A, const R130B3Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[1]*B[0];
    res[1]=-A[2]*B[0];
    res[2]=-A[3]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> R130B3<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0]*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    res[1]=(*this)[1]*x0;
    res[2]=(*this)[2]*x0;
    res[3]=(*this)[3]*x0;
    return res;
};

//-----------------------------------
// (R130B3, R130B4) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B1<T> operator*(const R130B3<T> &A, const R130B4<T> &B) {
    R130B1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B3<T> &A, const R130B4<T> &B) {
    R130B1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B3<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=A[0]*(-std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    return res;
};

//-----------------------------------
// (R130B3, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B3<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=(*this)[0];
    res[12]=(*this)[1];
    res[13]=(*this)[2];
    res[14]=(*this)[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=A[0] + B[11];
    res[12]=A[1] + B[12];
    res[13]=A[2] + B[13];
    res[14]=A[3] + B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=A[0] - B[11];
    res[12]=A[1] - B[12];
    res[13]=A[2] - B[13];
    res[14]=A[3] - B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[11] - A[1]*B[12] - A[2]*B[13] - A[3]*B[14];
    res[1]=-A[0]*B[15] - A[1]*B[8] - A[2]*B[9] - A[3]*B[10];
    res[2]=A[0]*B[8] + A[1]*B[15] + A[2]*B[7] - A[3]*B[6];
    res[3]=A[0]*B[9] - A[1]*B[7] + A[2]*B[15] + A[3]*B[5];
    res[4]=A[0]*B[10] + A[1]*B[6] - A[2]*B[5] + A[3]*B[15];
    res[5]=A[0]*B[12] - A[1]*B[11] + A[2]*B[4] - A[3]*B[3];
    res[6]=A[0]*B[13] - A[1]*B[4] - A[2]*B[11] + A[3]*B[2];
    res[7]=A[0]*B[14] + A[1]*B[3] - A[2]*B[2] - A[3]*B[11];
    res[8]=A[0]*B[2] + A[1]*B[1] - A[2]*B[14] + A[3]*B[13];
    res[9]=A[0]*B[3] + A[1]*B[14] + A[2]*B[1] - A[3]*B[12];
    res[10]=A[0]*B[4] - A[1]*B[13] + A[2]*B[12] + A[3]*B[1];
    res[11]=A[0]*B[0] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7];
    res[12]=A[0]*B[5] + A[1]*B[0] - A[2]*B[10] + A[3]*B[9];
    res[13]=A[0]*B[6] + A[1]*B[10] + A[2]*B[0] - A[3]*B[8];
    res[14]=A[0]*B[7] - A[1]*B[9] + A[2]*B[8] + A[3]*B[0];
    res[15]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3] - A[3]*B[4];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[11] - A[1]*B[12] - A[2]*B[13] - A[3]*B[14];
    res[1]=-A[1]*B[8] - A[2]*B[9] - A[3]*B[10];
    res[2]=A[0]*B[8] + A[2]*B[7] - A[3]*B[6];
    res[3]=A[0]*B[9] - A[1]*B[7] + A[3]*B[5];
    res[4]=A[0]*B[10] + A[1]*B[6] - A[2]*B[5];
    res[5]=A[2]*B[4] - A[3]*B[3];
    res[6]=-A[1]*B[4] + A[3]*B[2];
    res[7]=A[1]*B[3] - A[2]*B[2];
    res[8]=A[0]*B[2] + A[1]*B[1];
    res[9]=A[0]*B[3] + A[2]*B[1];
    res[10]=A[0]*B[4] + A[3]*B[1];
    res[11]=A[0]*B[0];
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[15];
    res[2]=A[1]*B[15];
    res[3]=A[2]*B[15];
    res[4]=A[3]*B[15];
    res[5]=A[0]*B[12] - A[1]*B[11];
    res[6]=A[0]*B[13] - A[2]*B[11];
    res[7]=A[0]*B[14] - A[3]*B[11];
    res[8]=-A[2]*B[14] + A[3]*B[13];
    res[9]=A[1]*B[14] - A[3]*B[12];
    res[10]=-A[1]*B[13] + A[2]*B[12];
    res[11]=A[1]*B[5] + A[2]*B[6] + A[3]*B[7];
    res[12]=A[0]*B[5] - A[2]*B[10] + A[3]*B[9];
    res[13]=A[0]*B[6] + A[1]*B[10] - A[3]*B[8];
    res[14]=A[0]*B[7] - A[1]*B[9] + A[2]*B[8];
    res[15]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3] - A[3]*B[4];
    return res;
};

template<typename T>
inline R130MV<T> R130B3<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[1]*x4;
    T x6 = (*this)[2]*A[3];
    T x7 = (*this)[3]*A[4];
    T x8 = 2.0*(*this)[1];
    T x9 = A[1]*x4;
    T x10 = A[2]*x8;
    T x11 = 2.0*(*this)[2];
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[2]*x8;
    T x14 = (*this)[3]*A[7];
    T x15 = (*this)[3]*x4;
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*x11;
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=-A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x5 - x4*x6 - x4*x7;
    res[2]=A[1]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + x6*x8 + x7*x8;
    res[3]=(*this)[2]*x10 + (*this)[2]*x9 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 + x11*x7;
    res[4]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3;
    res[5]=A[10]*x12 - A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[6]*x13 - A[9]*x15 + x14*x8;
    res[6]=-A[10]*x5 + A[5]*x13 - A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3 + A[8]*x15 + x11*x14;
    res[7]=A[5]*x16 + A[6]*x17 - A[7]*x0 - A[7]*x1 - A[7]*x2 + A[7]*x3 - A[8]*x12 + A[9]*x5;
    res[8]=-A[10]*x16 - A[6]*x15 + A[7]*x12 + A[8]*x0 - A[8]*x1 + A[8]*x2 + A[8]*x3 - A[9]*x13;
    res[9]=-A[10]*x17 + A[5]*x15 - A[7]*x5 - A[8]*x13 + A[9]*x0 + A[9]*x1 - A[9]*x2 + A[9]*x3;
    res[10]=A[10]*x0 + A[10]*x1 + A[10]*x2 - A[10]*x3 - A[5]*x12 + A[6]*x5 - A[8]*x16 - A[9]*x17;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 - A[12]*x5 - A[13]*x12 - A[14]*x15;
    res[12]=A[11]*x5 - A[12]*x0 - A[12]*x1 + A[12]*x2 + A[12]*x3 - A[13]*x13 - A[14]*x16;
    res[13]=A[11]*x12 - A[12]*x13 - A[13]*x0 + A[13]*x1 - A[13]*x2 + A[13]*x3 - A[14]*x17;
    res[14]=A[11]*x15 - A[12]*x16 - A[13]*x17 - A[14]*x0 + A[14]*x1 + A[14]*x2 - A[14]*x3;
    res[15]=A[15]*(-x0 + x1 + x2 + x3);
    return res;
};

//-------------------------------------
// (R130B3, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[1]*B[0] - A[2]*B[3] + A[3]*B[2];
    res[13]=A[1]*B[3] + A[2]*B[0] - A[3]*B[1];
    res[14]=-A[1]*B[2] + A[2]*B[1] + A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3<T> &A, const Rotation<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> R130B3<T>::conjugate(const Rotation<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[0];
    T x5 = A[3]*x4;
    T x6 = (*this)[3]*x4;
    T x7 = (*this)[1]*A[2];
    T x8 = (*this)[2]*A[1];
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*A[3];
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[2]*x5 - A[2]*x6;
    res[6]=-(*this)[1]*x5 + A[1]*x6;
    res[7]=x4*(x7 - x8);
    res[8]=A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - x10*x11 - x7*x9;
    res[9]=A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - x10*x8 - x11*x9;
    res[10]=-(*this)[3]*A[1]*x10 - (*this)[3]*A[2]*x9 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//----------------------------------
// (R130B3, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=B[4];
    res[9]=B[5];
    res[10]=B[6];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=-B[4];
    res[9]=-B[5];
    res[10]=-B[6];
    res[11]=A[0];
    res[12]=A[1];
    res[13]=A[2];
    res[14]=A[3];
    res[15]=-B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[7] - A[1]*B[4] - A[2]*B[5] - A[3]*B[6];
    res[2]=A[0]*B[4] + A[1]*B[7] + A[2]*B[3] - A[3]*B[2];
    res[3]=A[0]*B[5] - A[1]*B[3] + A[2]*B[7] + A[3]*B[1];
    res[4]=A[0]*B[6] + A[1]*B[2] - A[2]*B[1] + A[3]*B[7];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[12]=A[0]*B[1] + A[1]*B[0] - A[2]*B[6] + A[3]*B[5];
    res[13]=A[0]*B[2] + A[1]*B[6] + A[2]*B[0] - A[3]*B[4];
    res[14]=A[0]*B[3] - A[1]*B[5] + A[2]*B[4] + A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[4] - A[2]*B[5] - A[3]*B[6];
    res[2]=A[0]*B[4] + A[2]*B[3] - A[3]*B[2];
    res[3]=A[0]*B[5] - A[1]*B[3] + A[3]*B[1];
    res[4]=A[0]*B[6] + A[1]*B[2] - A[2]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[7];
    res[2]=A[1]*B[7];
    res[3]=A[2]*B[7];
    res[4]=A[3]*B[7];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[12]=A[0]*B[1] - A[2]*B[6] + A[3]*B[5];
    res[13]=A[0]*B[2] + A[1]*B[6] - A[3]*B[4];
    res[14]=A[0]*B[3] - A[1]*B[5] + A[2]*B[4];
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B3<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = (*this)[0]*A[6];
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*x7;
    T x9 = (*this)[0]*x7;
    T x10 = (*this)[3]*x4;
    T x11 = 2.0*(*this)[1];
    T x12 = (*this)[0]*x11;
    T x13 = (*this)[0]*x4;
    res[0]=A[0]*(x0 - x1 - x2 - x3);
    res[1]=-A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x6 + A[3]*x8 - A[5]*x9 + x4*x5;
    res[2]=A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[3]*x10 + A[4]*x9 - x11*x5;
    res[3]=A[1]*x8 + A[2]*x10 - A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 - A[4]*x13 + A[5]*x12;
    res[4]=-A[2]*x9 + A[3]*x13 + A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3 - A[5]*x6 - A[6]*x8;
    res[5]=A[1]*x9 - A[3]*x12 - A[4]*x6 + A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x3 - A[6]*x10;
    res[6]=-A[1]*x13 + A[2]*x12 - A[4]*x8 - A[5]*x10 + A[6]*x0 + A[6]*x1 + A[6]*x2 - A[6]*x3;
    res[7]=A[7]*(-x0 + x1 + x2 + x3);
    return res;
};

//-----------------------------------
// (R130B3, scalar) binary operations
//-----------------------------------

template<typename T>
inline R130B3<T> operator*(const R130B3<T> &A, const T &B) {
    R130B3<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    return res;
};
template<typename T>
inline R130B3<T> operator*(const T &A, const R130B3<T> &B) {
    return B*A;
};
template<typename T>
inline R130B3<T> operator/(const R130B3<T> &A, const T &B) {
    return A * (1.0 / B);
};

//-------------------------------------
// (R130B3Sm1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[1]*B[3] - A[2]*B[2];
    res[3]=-A[0]*B[3] + A[2]*B[1];
    res[4]=A[0]*B[2] - A[1]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[12]=A[0]*B[0];
    res[13]=A[1]*B[0];
    res[14]=A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sm1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[1]*B[3] - A[2]*B[2];
    res[3]=-A[0]*B[3] + A[2]*B[1];
    res[4]=A[0]*B[2] - A[1]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[1]*B[0];
    res[14]=A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator^(const R130B3Sm1<T> &A, const Boost<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    return res;
};

template<typename T>
inline Boost<T> R130B3Sm1<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//--------------------------------------
// (R130B3Sm1, R130B0) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
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
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
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
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B3Sm1<T> &A, const R130B0<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator/(const R130B3Sm1<T> &A, const R130B0<T> &B) {
    R130B3Sm1<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator|(const R130B3Sm1<T> &A, const R130B0<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sm1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B3Sm1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B3Sm1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1]*B[3] - A[2]*B[2];
    res[6]=-A[0]*B[3] + A[2]*B[1];
    res[7]=A[0]*B[2] - A[1]*B[1];
    res[8]=A[0]*B[0];
    res[9]=A[1]*B[0];
    res[10]=A[2]*B[0];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B3Sm1<T> &A, const R130B1<T> &B) {
    R130B2<T> res;
    res[0]=A[1]*B[3] - A[2]*B[2];
    res[1]=-A[0]*B[3] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[1]*B[1];
    res[3]=A[0]*B[0];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3Sm1<T> &A, const R130B1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B3Sm1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=(*this)[2]*A[1]*x3 + (*this)[2]*A[2]*x6 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B3Sm1, R130B1Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1]*B[2] - A[2]*B[1];
    res[6]=-A[0]*B[2] + A[2]*B[0];
    res[7]=A[0]*B[1] - A[1]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator|(const R130B3Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[1]*B[2] - A[2]*B[1];
    res[1]=-A[0]*B[2] + A[2]*B[0];
    res[2]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3Sm1<T> &A, const R130B1Sm1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B3Sm1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B3Sm1, R130B1Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
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
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
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
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B3Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B3Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sm1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B3Sm1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B3Sm1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[3] - A[1]*B[4] - A[2]*B[5];
    res[2]=A[1]*B[2] - A[2]*B[1];
    res[3]=-A[0]*B[2] + A[2]*B[0];
    res[4]=A[0]*B[1] - A[1]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[12]=-A[1]*B[5] + A[2]*B[4];
    res[13]=A[0]*B[5] - A[2]*B[3];
    res[14]=-A[0]*B[4] + A[1]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator|(const R130B3Sm1<T> &A, const R130B2<T> &B) {
    R130B1<T> res;
    res[0]=-A[0]*B[3] - A[1]*B[4] - A[2]*B[5];
    res[1]=A[1]*B[2] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B3Sm1<T> &A, const R130B2<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[1]=-A[1]*B[5] + A[2]*B[4];
    res[2]=A[0]*B[5] - A[2]*B[3];
    res[3]=-A[0]*B[4] + A[1]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> R130B3Sm1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=A[0]*x7 + A[1]*x8 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    res[3]=-A[3]*x0 + A[3]*x1 + A[3]*x2 - A[4]*x4 - A[5]*x7;
    res[4]=-A[3]*x4 + A[4]*x0 - A[4]*x1 + A[4]*x2 - A[5]*x8;
    res[5]=-A[3]*x7 - A[4]*x8 + A[5]*x0 + A[5]*x1 - A[5]*x2;
    return res;
};

//-----------------------------------------
// (R130B3Sm1, R130B2Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
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
    res[12]=-A[1]*B[2] + A[2]*B[1];
    res[13]=A[0]*B[2] - A[2]*B[0];
    res[14]=-A[0]*B[1] + A[1]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator|(const R130B3Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3Sm1<T> &A, const R130B2Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B3Sm1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-(*this)[2]*A[0]*x3 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B3Sm1, R130B2Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[1]*B[2] - A[2]*B[1];
    res[3]=-A[0]*B[2] + A[2]*B[0];
    res[4]=A[0]*B[1] - A[1]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B3Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[1]*B[2] - A[2]*B[1];
    res[1]=-A[0]*B[2] + A[2]*B[0];
    res[2]=A[0]*B[1] - A[1]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator^(const R130B3Sm1<T> &A, const R130B2Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B3Sm1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 + A[1]*x4 + x3*x5;
    res[1]=A[0]*x4 - A[1]*x0 + A[1]*x1 - A[1]*x2 + x5*x6;
    res[2]=(*this)[2]*A[0]*x3 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//--------------------------------------
// (R130B3Sm1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sm1<T>::operator R130B3<T>() const {
    R130B3<T> res;
    res[0]=0;
    res[1]=(*this)[0];
    res[2]=(*this)[1];
    res[3]=(*this)[2];
    return res;
};

template<typename T>
inline R130B3<T> operator+(const R130B3Sm1<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=B[0];
    res[1]=A[0] + B[1];
    res[2]=A[1] + B[2];
    res[3]=A[2] + B[3];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3Sm1<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=-B[0];
    res[1]=A[0] - B[1];
    res[2]=A[1] - B[2];
    res[3]=A[2] - B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[0]*B[0];
    res[6]=-A[1]*B[0];
    res[7]=-A[2]*B[0];
    res[8]=-A[1]*B[3] + A[2]*B[2];
    res[9]=A[0]*B[3] - A[2]*B[1];
    res[10]=-A[0]*B[2] + A[1]*B[1];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3Sm1<T> &A, const R130B3<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const R130B3Sm1<T> &A, const R130B3<T> &B) {
    R130B2<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    res[3]=-A[1]*B[3] + A[2]*B[2];
    res[4]=A[0]*B[3] - A[2]*B[1];
    res[5]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline R130B3<T> R130B3Sm1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-(*this)[2]*A[1]*x3 - (*this)[2]*A[2]*x6 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    return res;
};

//-----------------------------------------
// (R130B3Sm1, R130B3Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130B3Sm1<T> operator+(const R130B3Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator-(const R130B3Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const R130B3Sm1<T> &A, const R130B3Sm1<T> &B) {
    Rotation<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    res[1]=-A[1]*B[2] + A[2]*B[1];
    res[2]=A[0]*B[2] - A[2]*B[0];
    res[3]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[0] - A[1]*B[1] - A[2]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const R130B3Sm1<T> &A, const R130B3Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=-A[1]*B[2] + A[2]*B[1];
    res[1]=A[0]*B[2] - A[2]*B[0];
    res[2]=-A[0]*B[1] + A[1]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B3Sm1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[1], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[0], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x2 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x2 - x5*x6;
    res[2]=-(*this)[2]*A[0]*x3 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x2;
    return res;
};

//-----------------------------------------
// (R130B3Sm1, R130B3Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130B3<T> operator+(const R130B3Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B3<T> res;
    res[0]=B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B3<T> res;
    res[0]=-B[0];
    res[1]=A[0];
    res[2]=A[1];
    res[3]=A[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B3Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B3Sm1<T> &A, const R130B3Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[1]*B[0];
    res[2]=-A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B3Sm1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B3Sm1, R130B4) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator*(const R130B3Sm1<T> &A, const R130B4<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sm1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B3Sm1<T> &A, const R130B4<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B3Sm1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2));
    return res;
};

//--------------------------------------
// (R130B3Sm1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sm1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
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
    res[12]=(*this)[0];
    res[13]=(*this)[1];
    res[14]=(*this)[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=A[0] + B[12];
    res[13]=A[1] + B[13];
    res[14]=A[2] + B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=A[0] - B[12];
    res[13]=A[1] - B[13];
    res[14]=A[2] - B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[1]=-A[0]*B[8] - A[1]*B[9] - A[2]*B[10];
    res[2]=A[0]*B[15] + A[1]*B[7] - A[2]*B[6];
    res[3]=-A[0]*B[7] + A[1]*B[15] + A[2]*B[5];
    res[4]=A[0]*B[6] - A[1]*B[5] + A[2]*B[15];
    res[5]=-A[0]*B[11] + A[1]*B[4] - A[2]*B[3];
    res[6]=-A[0]*B[4] - A[1]*B[11] + A[2]*B[2];
    res[7]=A[0]*B[3] - A[1]*B[2] - A[2]*B[11];
    res[8]=A[0]*B[1] - A[1]*B[14] + A[2]*B[13];
    res[9]=A[0]*B[14] + A[1]*B[1] - A[2]*B[12];
    res[10]=-A[0]*B[13] + A[1]*B[12] + A[2]*B[1];
    res[11]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7];
    res[12]=A[0]*B[0] - A[1]*B[10] + A[2]*B[9];
    res[13]=A[0]*B[10] + A[1]*B[0] - A[2]*B[8];
    res[14]=-A[0]*B[9] + A[1]*B[8] + A[2]*B[0];
    res[15]=-A[0]*B[2] - A[1]*B[3] - A[2]*B[4];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[12] - A[1]*B[13] - A[2]*B[14];
    res[1]=-A[0]*B[8] - A[1]*B[9] - A[2]*B[10];
    res[2]=A[1]*B[7] - A[2]*B[6];
    res[3]=-A[0]*B[7] + A[2]*B[5];
    res[4]=A[0]*B[6] - A[1]*B[5];
    res[5]=A[1]*B[4] - A[2]*B[3];
    res[6]=-A[0]*B[4] + A[2]*B[2];
    res[7]=A[0]*B[3] - A[1]*B[2];
    res[8]=A[0]*B[1];
    res[9]=A[1]*B[1];
    res[10]=A[2]*B[1];
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[1]*B[0];
    res[14]=A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sm1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[15];
    res[3]=A[1]*B[15];
    res[4]=A[2]*B[15];
    res[5]=-A[0]*B[11];
    res[6]=-A[1]*B[11];
    res[7]=-A[2]*B[11];
    res[8]=-A[1]*B[14] + A[2]*B[13];
    res[9]=A[0]*B[14] - A[2]*B[12];
    res[10]=-A[0]*B[13] + A[1]*B[12];
    res[11]=A[0]*B[5] + A[1]*B[6] + A[2]*B[7];
    res[12]=-A[1]*B[10] + A[2]*B[9];
    res[13]=A[0]*B[10] - A[2]*B[8];
    res[14]=-A[0]*B[9] + A[1]*B[8];
    res[15]=-A[0]*B[2] - A[1]*B[3] - A[2]*B[4];
    return res;
};

template<typename T>
inline R130MV<T> R130B3Sm1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[4];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*(x0 + x1 + x2);
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=A[5]*x0 - A[5]*x1 - A[5]*x2 + A[6]*x4 + A[7]*x7;
    res[6]=A[5]*x4 - A[6]*x0 + A[6]*x1 - A[6]*x2 + A[7]*x8;
    res[7]=A[5]*x7 + A[6]*x8 - A[7]*x0 - A[7]*x1 + A[7]*x2;
    res[8]=-A[10]*x7 - A[8]*x0 + A[8]*x1 + A[8]*x2 - A[9]*x4;
    res[9]=-A[10]*x8 - A[8]*x4 + A[9]*x0 - A[9]*x1 + A[9]*x2;
    res[10]=A[10]*x0 + A[10]*x1 - A[10]*x2 - A[8]*x7 - A[9]*x8;
    res[11]=A[11]*(x0 + x1 + x2);
    res[12]=-A[12]*x0 + A[12]*x1 + A[12]*x2 - A[13]*x4 - A[14]*x7;
    res[13]=-A[12]*x4 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[14]*x8;
    res[14]=-A[12]*x7 - A[13]*x8 + A[14]*x0 + A[14]*x1 - A[14]*x2;
    res[15]=A[15]*(x0 + x1 + x2);
    return res;
};

//----------------------------------------
// (R130B3Sm1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
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
    res[12]=A[0]*B[0] - A[1]*B[3] + A[2]*B[2];
    res[13]=A[0]*B[3] + A[1]*B[0] - A[2]*B[1];
    res[14]=-A[0]*B[2] + A[1]*B[1] + A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sm1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[1] - A[1]*B[2] - A[2]*B[3];
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
    res[12]=A[0]*B[0];
    res[13]=A[1]*B[0];
    res[14]=A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3Sm1<T> &A, const Rotation<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[1]*B[3] + A[2]*B[2];
    res[1]=A[0]*B[3] - A[2]*B[1];
    res[2]=-A[0]*B[2] + A[1]*B[1];
    return res;
};

template<typename T>
inline Rotation<T> R130B3Sm1<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-(*this)[2]*A[1]*x3 - (*this)[2]*A[2]*x6 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    return res;
};

//-------------------------------------
// (R130B3Sm1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=B[4];
    res[9]=B[5];
    res[10]=B[6];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=-B[4];
    res[9]=-B[5];
    res[10]=-B[6];
    res[11]=0;
    res[12]=A[0];
    res[13]=A[1];
    res[14]=A[2];
    res[15]=-B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[4] - A[1]*B[5] - A[2]*B[6];
    res[2]=A[0]*B[7] + A[1]*B[3] - A[2]*B[2];
    res[3]=-A[0]*B[3] + A[1]*B[7] + A[2]*B[1];
    res[4]=A[0]*B[2] - A[1]*B[1] + A[2]*B[7];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[12]=A[0]*B[0] - A[1]*B[6] + A[2]*B[5];
    res[13]=A[0]*B[6] + A[1]*B[0] - A[2]*B[4];
    res[14]=-A[0]*B[5] + A[1]*B[4] + A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[4] - A[1]*B[5] - A[2]*B[6];
    res[2]=A[1]*B[3] - A[2]*B[2];
    res[3]=-A[0]*B[3] + A[2]*B[1];
    res[4]=A[0]*B[2] - A[1]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[1]*B[0];
    res[14]=A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sm1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[7];
    res[3]=A[1]*B[7];
    res[4]=A[2]*B[7];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3];
    res[12]=-A[1]*B[6] + A[2]*B[5];
    res[13]=A[0]*B[6] - A[2]*B[4];
    res[14]=-A[0]*B[5] + A[1]*B[4];
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B3Sm1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[3];
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x3;
    T x8 = (*this)[2]*x6;
    res[0]=-A[0]*(x0 + x1 + x2);
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=A[1]*x7 + A[2]*x8 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    res[4]=-A[4]*x0 + A[4]*x1 + A[4]*x2 - A[5]*x4 - A[6]*x7;
    res[5]=-A[4]*x4 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[6]*x8;
    res[6]=-A[4]*x7 - A[5]*x8 + A[6]*x0 + A[6]*x1 - A[6]*x2;
    res[7]=A[7]*(x0 + x1 + x2);
    return res;
};

//--------------------------------------
// (R130B3Sm1, scalar) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sm1<T> operator*(const R130B3Sm1<T> &A, const T &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    return res;
};
template<typename T>
inline R130B3Sm1<T> operator*(const T &A, const R130B3Sm1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B3Sm1<T> operator/(const R130B3Sm1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//-------------------------------------
// (R130B3Sp1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator*(const R130B3Sp1<T> &A, const Boost<T> &B) {
    R130B3<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator|(const R130B3Sp1<T> &A, const Boost<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3Sp1<T> &A, const Boost<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[1];
    res[1]=A[0]*B[2];
    res[2]=A[0]*B[3];
    return res;
};

template<typename T>
inline Boost<T> R130B3Sp1<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    return res;
};

//--------------------------------------
// (R130B3Sp1, R130B0) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
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
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
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
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator*(const R130B3Sp1<T> &A, const R130B0<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator/(const R130B3Sp1<T> &A, const R130B0<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]/B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator|(const R130B3Sp1<T> &A, const R130B0<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B3Sp1<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B3Sp1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sp1<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[1];
    res[9]=A[0]*B[2];
    res[10]=A[0]*B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B3Sp1<T> &A, const R130B1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[1];
    res[1]=A[0]*B[2];
    res[2]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3Sp1<T> &A, const R130B1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1<T> R130B3Sp1<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//-----------------------------------------
// (R130B3Sp1, R130B1Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B3Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B3Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B3Sp1<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//-----------------------------------------
// (R130B3Sp1, R130B1Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> operator*(const R130B3Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> operator^(const R130B3Sp1<T> &A, const R130B1Sp1<T> &B) {
    R130B4<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B3Sp1<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=-std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B3Sp1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sp1<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[3];
    res[3]=A[0]*B[4];
    res[4]=A[0]*B[5];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B3Sp1<T> &A, const R130B2<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[3];
    res[1]=A[0]*B[4];
    res[2]=A[0]*B[5];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3Sp1<T> &A, const R130B2<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> R130B3Sp1<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    return res;
};

//-----------------------------------------
// (R130B3Sp1, R130B2Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator*(const R130B3Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator|(const R130B3Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B3Sp1<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//-----------------------------------------
// (R130B3Sp1, R130B2Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B3Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B3Sp1<T> &A, const R130B2Sp1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B3Sp1<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B3Sp1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sp1<T>::operator R130B3<T>() const {
    R130B3<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator+(const R130B3Sp1<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3Sp1<T> &A, const R130B3<T> &B) {
    R130B3<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    return res;
};

template<typename T>
inline Boost<T> operator*(const R130B3Sp1<T> &A, const R130B3<T> &B) {
    Boost<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3Sp1<T> &A, const R130B3<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B3Sp1<T> &A, const R130B3<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[1];
    res[1]=A[0]*B[2];
    res[2]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B3Sp1<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    return res;
};

//-----------------------------------------
// (R130B3Sp1, R130B3Sm1) binary operations
//-----------------------------------------

template<typename T>
inline R130B3<T> operator+(const R130B3Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B3<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    return res;
};

template<typename T>
inline R130B3<T> operator-(const R130B3Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B3<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B3Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const R130B3Sp1<T> &A, const R130B3Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B3Sp1<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    return res;
};

//-----------------------------------------
// (R130B3Sp1, R130B3Sp1) binary operations
//-----------------------------------------

template<typename T>
inline R130B3Sp1<T> operator+(const R130B3Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0] + B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator-(const R130B3Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=A[0] - B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator*(const R130B3Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B3Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130B0<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B3Sp1<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B3Sp1, R130B4) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator*(const R130B3Sp1<T> &A, const R130B4<T> &B) {
    R130B1Sp1<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator^(const R130B3Sp1<T> &A, const R130B4<T> &B) {
    R130B1Sp1<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> R130B3Sp1<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=-std::pow((*this)[0], 2)*A[0];
    return res;
};

//--------------------------------------
// (R130B3Sp1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sp1<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=(*this)[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=A[0] + B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=A[0] - B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[11];
    res[1]=-A[0]*B[15];
    res[2]=A[0]*B[8];
    res[3]=A[0]*B[9];
    res[4]=A[0]*B[10];
    res[5]=A[0]*B[12];
    res[6]=A[0]*B[13];
    res[7]=A[0]*B[14];
    res[8]=A[0]*B[2];
    res[9]=A[0]*B[3];
    res[10]=A[0]*B[4];
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[5];
    res[13]=A[0]*B[6];
    res[14]=A[0]*B[7];
    res[15]=-A[0]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[11];
    res[1]=0;
    res[2]=A[0]*B[8];
    res[3]=A[0]*B[9];
    res[4]=A[0]*B[10];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[2];
    res[9]=A[0]*B[3];
    res[10]=A[0]*B[4];
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[15];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[12];
    res[6]=A[0]*B[13];
    res[7]=A[0]*B[14];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[5];
    res[13]=A[0]*B[6];
    res[14]=A[0]*B[7];
    res[15]=-A[0]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> R130B3Sp1<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=-A[5]*x0;
    res[6]=-A[6]*x0;
    res[7]=-A[7]*x0;
    res[8]=A[8]*x0;
    res[9]=A[9]*x0;
    res[10]=A[10]*x0;
    res[11]=A[11]*x0;
    res[12]=-A[12]*x0;
    res[13]=-A[13]*x0;
    res[14]=-A[14]*x0;
    res[15]=-A[15]*x0;
    return res;
};

//----------------------------------------
// (R130B3Sp1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Rotation<T> R130B3Sp1<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//-------------------------------------
// (R130B3Sp1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B3Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=B[4];
    res[9]=B[5];
    res[10]=B[6];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B3Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=-B[4];
    res[9]=-B[5];
    res[10]=-B[6];
    res[11]=A[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B3Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[7];
    res[2]=A[0]*B[4];
    res[3]=A[0]*B[5];
    res[4]=A[0]*B[6];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B3Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[4];
    res[3]=A[0]*B[5];
    res[4]=A[0]*B[6];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B3Sp1<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[0]*B[7];
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
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B3Sp1<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    res[6]=A[6]*x0;
    res[7]=-A[7]*x0;
    return res;
};

//--------------------------------------
// (R130B3Sp1, scalar) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sp1<T> operator*(const R130B3Sp1<T> &A, const T &B) {
    R130B3Sp1<T> res;
    res[0]=A[0]*B;
    return res;
};
template<typename T>
inline R130B3Sp1<T> operator*(const T &A, const R130B3Sp1<T> &B) {
    return B*A;
};
template<typename T>
inline R130B3Sp1<T> operator/(const R130B3Sp1<T> &A, const T &B) {
    return A * (1.0 / B);
};

//----------------------------------
// (R130B4, Boost) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B4<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[1];
    res[9]=A[0]*B[2];
    res[10]=A[0]*B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[0]*B[1];
    res[9]=A[0]*B[2];
    res[10]=A[0]*B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Boost<T> R130B4<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    return res;
};

//-----------------------------------
// (R130B4, R130B0) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B4<T> operator*(const R130B4<T> &A, const R130B0<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator/(const R130B4<T> &A, const R130B0<T> &B) {
    R130B4<T> res;
    res[0]=A[0]/B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator|(const R130B4<T> &A, const R130B0<T> &B) {
    R130B4<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> R130B4<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=-std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B4, R130B1) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B3<T> operator*(const R130B4<T> &A, const R130B1<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> operator^(const R130B4<T> &A, const R130B1<T> &B) {
    R130B3<T> res;
    res[0]=-A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    return res;
};

template<typename T>
inline R130B1<T> R130B4<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//--------------------------------------
// (R130B4, R130B1Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator*(const R130B4<T> &A, const R130B1Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const R130B4<T> &A, const R130B1Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B1Sm1<T> R130B4<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B4, R130B1Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=B[0];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-B[0];
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
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator*(const R130B4<T> &A, const R130B1Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> operator^(const R130B4<T> &A, const R130B1Sp1<T> &B) {
    R130B3Sp1<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> R130B4<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B4, R130B2) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=B[3];
    res[9]=B[4];
    res[10]=B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=-B[3];
    res[9]=-B[4];
    res[10]=-B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B2<T> operator*(const R130B4<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[0]*B[3];
    res[1]=-A[0]*B[4];
    res[2]=-A[0]*B[5];
    res[3]=A[0]*B[0];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator|(const R130B4<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[0]*B[3];
    res[1]=-A[0]*B[4];
    res[2]=-A[0]*B[5];
    res[3]=A[0]*B[0];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2<T> R130B4<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    res[4]=-A[4]*x0;
    res[5]=-A[5]*x0;
    return res;
};

//--------------------------------------
// (R130B4, R130B2Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[0];
    res[9]=B[1];
    res[10]=B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[0];
    res[9]=-B[1];
    res[10]=-B[2];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator*(const R130B4<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator|(const R130B4<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sm1<T> R130B4<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B4, R130B2Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator*(const R130B4<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator|(const R130B4<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B2Sp1<T> R130B4<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    return res;
};

//-----------------------------------
// (R130B4, R130B3) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B1<T> operator*(const R130B4<T> &A, const R130B3<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    res[3]=-A[0]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> operator^(const R130B4<T> &A, const R130B3<T> &B) {
    R130B1<T> res;
    res[0]=A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    res[3]=-A[0]*B[3];
    return res;
};

template<typename T>
inline R130B3<T> R130B4<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

//--------------------------------------
// (R130B4, R130B3Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator*(const R130B4<T> &A, const R130B3Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const R130B4<T> &A, const R130B3Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[0]*B[0];
    res[1]=-A[0]*B[1];
    res[2]=-A[0]*B[2];
    return res;
};

template<typename T>
inline R130B3Sm1<T> R130B4<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    return res;
};

//--------------------------------------
// (R130B4, R130B3Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator*(const R130B4<T> &A, const R130B3Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> operator^(const R130B4<T> &A, const R130B3Sp1<T> &B) {
    R130B1Sp1<T> res;
    res[0]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130B3Sp1<T> R130B4<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B4, R130B4) binary operations
//-----------------------------------

template<typename T>
inline R130B4<T> operator+(const R130B4<T> &A, const R130B4<T> &B) {
    R130B4<T> res;
    res[0]=A[0] + B[0];
    return res;
};

template<typename T>
inline R130B4<T> operator-(const R130B4<T> &A, const R130B4<T> &B) {
    R130B4<T> res;
    res[0]=A[0] - B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator*(const R130B4<T> &A, const R130B4<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130B0<T> operator|(const R130B4<T> &A, const R130B4<T> &B) {
    R130B0<T> res;
    res[0]=-A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> R130B4<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=-std::pow((*this)[0], 2)*A[0];
    return res;
};

//-----------------------------------
// (R130B4, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B4<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=(*this)[0];
    return res;
};

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=B[8];
    res[9]=B[9];
    res[10]=B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=A[0] + B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=-B[8];
    res[9]=-B[9];
    res[10]=-B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=A[0] - B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B4<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=-A[0]*B[15];
    res[1]=A[0]*B[11];
    res[2]=-A[0]*B[12];
    res[3]=-A[0]*B[13];
    res[4]=-A[0]*B[14];
    res[5]=-A[0]*B[8];
    res[6]=-A[0]*B[9];
    res[7]=-A[0]*B[10];
    res[8]=A[0]*B[5];
    res[9]=A[0]*B[6];
    res[10]=A[0]*B[7];
    res[11]=-A[0]*B[1];
    res[12]=A[0]*B[2];
    res[13]=A[0]*B[3];
    res[14]=A[0]*B[4];
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B4<T> &A, const R130MV<T> &B) {
    Rotor<T> res;
    res[0]=-A[0]*B[15];
    res[1]=-A[0]*B[8];
    res[2]=-A[0]*B[9];
    res[3]=-A[0]*B[10];
    res[4]=A[0]*B[5];
    res[5]=A[0]*B[6];
    res[6]=A[0]*B[7];
    res[7]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[11];
    res[2]=-A[0]*B[12];
    res[3]=-A[0]*B[13];
    res[4]=-A[0]*B[14];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[0]*B[1];
    res[12]=A[0]*B[2];
    res[13]=A[0]*B[3];
    res[14]=A[0]*B[4];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130B4<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=-A[5]*x0;
    res[6]=-A[6]*x0;
    res[7]=-A[7]*x0;
    res[8]=-A[8]*x0;
    res[9]=-A[9]*x0;
    res[10]=-A[10]*x0;
    res[11]=A[11]*x0;
    res[12]=A[12]*x0;
    res[13]=A[13]*x0;
    res[14]=A[14]*x0;
    res[15]=-A[15]*x0;
    return res;
};

//-------------------------------------
// (R130B4, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130B4<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=B[1];
    res[9]=B[2];
    res[10]=B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130B4<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=-B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=-B[1];
    res[9]=-B[2];
    res[10]=-B[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130B4<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[0]*B[1];
    res[6]=-A[0]*B[2];
    res[7]=-A[0]*B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130B4<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[0]*B[1];
    res[6]=-A[0]*B[2];
    res[7]=-A[0]*B[3];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Rotation<T> R130B4<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    return res;
};

//----------------------------------
// (R130B4, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130B4<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=(*this)[0];
    return res;
};

template<typename T>
inline Rotor<T> operator+(const R130B4<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=A[0] + B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const R130B4<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=A[0] - B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const R130B4<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-A[0]*B[7];
    res[1]=-A[0]*B[4];
    res[2]=-A[0]*B[5];
    res[3]=-A[0]*B[6];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    res[6]=A[0]*B[3];
    res[7]=A[0]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130B4<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=-A[0]*B[7];
    res[1]=-A[0]*B[4];
    res[2]=-A[0]*B[5];
    res[3]=-A[0]*B[6];
    res[4]=A[0]*B[1];
    res[5]=A[0]*B[2];
    res[6]=A[0]*B[3];
    res[7]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130B4<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> R130B4<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    res[0]=-A[0]*x0;
    res[1]=-A[1]*x0;
    res[2]=-A[2]*x0;
    res[3]=-A[3]*x0;
    res[4]=-A[4]*x0;
    res[5]=-A[5]*x0;
    res[6]=-A[6]*x0;
    res[7]=-A[7]*x0;
    return res;
};

//-----------------------------------
// (R130B4, scalar) binary operations
//-----------------------------------

template<typename T>
inline R130B4<T> operator*(const R130B4<T> &A, const T &B) {
    R130B4<T> res;
    res[0]=A[0]*B;
    return res;
};
template<typename T>
inline R130B4<T> operator*(const T &A, const R130B4<T> &B) {
    return B*A;
};
template<typename T>
inline R130B4<T> operator/(const R130B4<T> &A, const T &B) {
    return A * (1.0 / B);
};

//----------------------------------
// (R130MV, Boost) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] + B[1];
    res[6]=A[6] + B[2];
    res[7]=A[7] + B[3];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] - B[1];
    res[6]=A[6] - B[2];
    res[7]=A[7] - B[3];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3];
    res[1]=A[1]*B[0] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[2]=A[13]*B[3] - A[14]*B[2] - A[1]*B[1] + A[2]*B[0];
    res[3]=-A[12]*B[3] + A[14]*B[1] - A[1]*B[2] + A[3]*B[0];
    res[4]=A[12]*B[2] - A[13]*B[1] - A[1]*B[3] + A[4]*B[0];
    res[5]=A[0]*B[1] + A[10]*B[2] + A[5]*B[0] - A[9]*B[3];
    res[6]=A[0]*B[2] - A[10]*B[1] + A[6]*B[0] + A[8]*B[3];
    res[7]=A[0]*B[3] + A[7]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[8]=A[15]*B[1] + A[6]*B[3] - A[7]*B[2] + A[8]*B[0];
    res[9]=A[15]*B[2] - A[5]*B[3] + A[7]*B[1] + A[9]*B[0];
    res[10]=A[10]*B[0] + A[15]*B[3] + A[5]*B[2] - A[6]*B[1];
    res[11]=A[11]*B[0] + A[12]*B[1] + A[13]*B[2] + A[14]*B[3];
    res[12]=A[11]*B[1] + A[12]*B[0] - A[3]*B[3] + A[4]*B[2];
    res[13]=A[11]*B[2] + A[13]*B[0] + A[2]*B[3] - A[4]*B[1];
    res[14]=A[11]*B[3] + A[14]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[15]=A[10]*B[3] + A[15]*B[0] + A[8]*B[1] + A[9]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3];
    res[1]=A[1]*B[0];
    res[2]=A[13]*B[3] - A[14]*B[2] + A[2]*B[0];
    res[3]=-A[12]*B[3] + A[14]*B[1] + A[3]*B[0];
    res[4]=A[12]*B[2] - A[13]*B[1] + A[4]*B[0];
    res[5]=A[0]*B[1] + A[5]*B[0];
    res[6]=A[0]*B[2] + A[6]*B[0];
    res[7]=A[0]*B[3] + A[7]*B[0];
    res[8]=A[15]*B[1] + A[8]*B[0];
    res[9]=A[15]*B[2] + A[9]*B[0];
    res[10]=A[10]*B[0] + A[15]*B[3];
    res[11]=A[11]*B[0];
    res[12]=A[12]*B[0] - A[3]*B[3] + A[4]*B[2];
    res[13]=A[13]*B[0] + A[2]*B[3] - A[4]*B[1];
    res[14]=A[14]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[15]=A[10]*B[3] + A[15]*B[0] + A[8]*B[1] + A[9]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[2]=-A[1]*B[1];
    res[3]=-A[1]*B[2];
    res[4]=-A[1]*B[3];
    res[5]=A[10]*B[2] - A[9]*B[3];
    res[6]=-A[10]*B[1] + A[8]*B[3];
    res[7]=-A[8]*B[2] + A[9]*B[1];
    res[8]=A[6]*B[3] - A[7]*B[2];
    res[9]=-A[5]*B[3] + A[7]*B[1];
    res[10]=A[5]*B[2] - A[6]*B[1];
    res[11]=A[12]*B[1] + A[13]*B[2] + A[14]*B[3];
    res[12]=A[11]*B[1];
    res[13]=A[11]*B[2];
    res[14]=A[11]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const Boost<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[10], 2);
    T x2 = std::pow((*this)[11], 2);
    T x3 = std::pow((*this)[2], 2);
    T x4 = std::pow((*this)[3], 2);
    T x5 = std::pow((*this)[4], 2);
    T x6 = std::pow((*this)[8], 2);
    T x7 = std::pow((*this)[9], 2);
    T x8 = std::pow((*this)[12], 2);
    T x9 = std::pow((*this)[13], 2);
    T x10 = std::pow((*this)[14], 2);
    T x11 = std::pow((*this)[15], 2);
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[5], 2);
    T x14 = std::pow((*this)[6], 2);
    T x15 = std::pow((*this)[7], 2);
    T x16 = 2.0*(*this)[0];
    T x17 = A[1]*x16;
    T x18 = A[2]*x16;
    T x19 = A[3]*x16;
    T x20 = 2.0*A[1];
    T x21 = (*this)[3]*x20;
    T x22 = 2.0*A[3];
    T x23 = (*this)[12]*x22;
    T x24 = (*this)[7]*x20;
    T x25 = 2.0*A[2];
    T x26 = (*this)[5]*x25;
    T x27 = (*this)[9]*x22;
    T x28 = (*this)[8]*x25;
    T x29 = (*this)[11]*x22;
    T x30 = (*this)[10]*x25;
    T x31 = (*this)[11]*x20;
    T x32 = (*this)[11]*x25;
    T x33 = (*this)[12]*x20;
    T x34 = (*this)[12]*x25;
    T x35 = (*this)[15]*x25;
    T x36 = (*this)[5]*x22;
    T x37 = (*this)[14]*x22;
    T x38 = (*this)[14]*x20;
    T x39 = (*this)[5]*x20;
    T x40 = (*this)[6]*x25;
    T x41 = (*this)[7]*x22;
    T x42 = (*this)[3]*x22;
    T x43 = (*this)[9]*x20;
    T x44 = (*this)[10]*x22;
    T x45 = (*this)[10]*x20;
    T x46 = (*this)[13]*x20;
    T x47 = (*this)[9]*x25;
    T x48 = (*this)[15]*x20;
    T x49 = (*this)[7]*x25;
    T x50 = (*this)[15]*x22;
    T x51 = (*this)[1]*x22;
    T x52 = (*this)[2]*x20;
    T x53 = (*this)[4]*x22;
    T x54 = (*this)[14]*x25;
    T x55 = (*this)[3]*x25;
    T x56 = (*this)[13]*x25;
    T x57 = (*this)[1]*(*this)[4];
    T x58 = A[0]*x16;
    T x59 = 2.0*A[0];
    T x60 = (*this)[10]*x59;
    T x61 = (*this)[5]*x59;
    T x62 = (*this)[6]*x59;
    T x63 = (*this)[14]*x59;
    T x64 = (*this)[2]*x59;
    T x65 = (*this)[3]*x59;
    T x66 = (*this)[1]*x59;
    res[0]=A[0]*(x0 + x1 - x10 - x11 - x12 - x13 - x14 - x15 + x2 + x3 + x4 + x5 + x6 + x7 - x8 - x9);
    res[1]=(*this)[10]*x21 - (*this)[10]*x29 + (*this)[13]*x24 - (*this)[13]*x35 - (*this)[13]*x36 + (*this)[14]*x26 - (*this)[15]*x33 - (*this)[15]*x37 - (*this)[1]*x39 - (*this)[1]*x40 - (*this)[1]*x41 + (*this)[2]*x17 + (*this)[2]*x27 - (*this)[2]*x30 + (*this)[3]*x18 + (*this)[4]*x19 + (*this)[4]*x28 - (*this)[4]*x43 + (*this)[6]*x23 - (*this)[6]*x38 - (*this)[7]*x34 - (*this)[8]*x31 - (*this)[8]*x42 - (*this)[9]*x32;
    res[2]=(*this)[10]*x23 - (*this)[10]*x38 + (*this)[13]*x19 + (*this)[13]*x28 - (*this)[13]*x43 - (*this)[14]*x18 + (*this)[15]*x31 + (*this)[15]*x42 + (*this)[1]*x17 + (*this)[1]*x27 - (*this)[1]*x30 - (*this)[2]*x39 - (*this)[2]*x40 - (*this)[2]*x41 - (*this)[3]*x26 + (*this)[4]*x24 - (*this)[4]*x35 - (*this)[4]*x36 + (*this)[6]*x21 - (*this)[6]*x29 + (*this)[7]*x32 + (*this)[8]*x33 + (*this)[8]*x37 + (*this)[9]*x34;
    res[3]=-(*this)[11]*x24 - (*this)[12]*x19 - (*this)[12]*x28 + (*this)[13]*x44 + (*this)[13]*x47 + (*this)[14]*x17 + (*this)[14]*x27 - (*this)[14]*x30 + (*this)[15]*x32 + (*this)[1]*x18 + (*this)[1]*x45 + (*this)[2]*x26 - (*this)[2]*x50 - (*this)[3]*x40 - (*this)[3]*x41 + (*this)[4]*x48 + (*this)[4]*x49 - (*this)[5]*x21 + (*this)[5]*x29 - (*this)[6]*x52 - (*this)[6]*x53 + (*this)[8]*x46 - (*this)[8]*x51 + (*this)[9]*x33;
    res[4]=(*this)[10]*x33 + (*this)[10]*x37 - (*this)[11]*x26 + (*this)[12]*x18 - (*this)[13]*x17 - (*this)[13]*x27 + (*this)[13]*x30 + (*this)[14]*x47 - (*this)[15]*x21 + (*this)[15]*x29 + (*this)[1]*x19 + (*this)[1]*x28 - (*this)[1]*x43 - (*this)[2]*x24 + (*this)[2]*x35 + (*this)[2]*x36 - (*this)[3]*x49 - (*this)[4]*x39 - (*this)[4]*x40 - (*this)[4]*x41 + (*this)[6]*x31 + (*this)[6]*x42 - (*this)[8]*x23 + (*this)[8]*x38;
    res[5]=-(*this)[10]*x18 + (*this)[13]*x34 + (*this)[13]*x51 + (*this)[14]*x23 - (*this)[1]*x54 - (*this)[2]*x53 - (*this)[2]*x55 - (*this)[3]*x29 + (*this)[4]*x32 - (*this)[6]*x26 + (*this)[6]*x50 - (*this)[7]*x35 - (*this)[7]*x36 + (*this)[8]*x44 + (*this)[9]*x19 + (*this)[9]*x28 + A[1]*x0 - A[1]*x1 - A[1]*x10 - A[1]*x11 + A[1]*x12 - A[1]*x13 + A[1]*x14 + A[1]*x15 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 + A[1]*x8 - A[1]*x9;
    res[6]=(*this)[10]*x17 + (*this)[10]*x27 + (*this)[13]*x33 + (*this)[13]*x37 + (*this)[15]*x24 - (*this)[15]*x36 - (*this)[1]*x23 + (*this)[1]*x38 - (*this)[2]*x21 + (*this)[2]*x29 - (*this)[4]*x31 - (*this)[4]*x42 - (*this)[6]*x39 - (*this)[6]*x41 - (*this)[8]*x19 + (*this)[8]*x43 + A[2]*x0 - A[2]*x1 - A[2]*x10 - A[2]*x11 + A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 - A[2]*x8 + A[2]*x9;
    res[7]=(*this)[11]*x21 + (*this)[13]*x54 + (*this)[14]*x33 + (*this)[15]*x26 + (*this)[1]*x34 - (*this)[1]*x46 - (*this)[2]*x32 - (*this)[4]*x52 - (*this)[4]*x55 - (*this)[5]*x24 - (*this)[6]*x48 - (*this)[7]*x40 + (*this)[8]*x18 + (*this)[8]*x45 - (*this)[9]*x17 + (*this)[9]*x30 + A[3]*x0 + A[3]*x1 + A[3]*x10 - A[3]*x11 + A[3]*x12 + A[3]*x13 + A[3]*x14 - A[3]*x15 - A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 - A[3]*x6 - A[3]*x7 - A[3]*x8 - A[3]*x9;
    res[8]=(*this)[10]*x24 - (*this)[10]*x36 - (*this)[13]*x21 + (*this)[13]*x29 - (*this)[14]*x32 + (*this)[15]*x17 + (*this)[15]*x27 - (*this)[15]*x30 + (*this)[1]*x31 + (*this)[1]*x42 + (*this)[2]*x33 + (*this)[2]*x37 + (*this)[2]*x56 + (*this)[3]*x34 + (*this)[4]*x23 - (*this)[4]*x38 - (*this)[6]*x19 - (*this)[6]*x28 + (*this)[6]*x43 + (*this)[7]*x18 - (*this)[8]*x39 - (*this)[8]*x41 - (*this)[9]*x26 - x25*x57;
    res[9]=-(*this)[11]*x23 + (*this)[12]*x21 + (*this)[13]*x53 + (*this)[13]*x55 + (*this)[14]*x31 + (*this)[15]*x18 + (*this)[15]*x45 + (*this)[1]*x32 - (*this)[2]*x34 + (*this)[2]*x46 - (*this)[2]*x51 + (*this)[3]*x37 - (*this)[4]*x54 + (*this)[5]*x19 - (*this)[6]*(*this)[8]*x20 - (*this)[6]*x44 - (*this)[7]*x17 - (*this)[7]*x27 + (*this)[7]*x30 + (*this)[8]*x26 - (*this)[8]*x50 - (*this)[9]*x39 - (*this)[9]*x40 + x20*x57;
    res[10]=-(*this)[10]*x39 - (*this)[10]*x41 + (*this)[12]*x32 - (*this)[13]*x31 - (*this)[13]*x42 + (*this)[15]*x19 + (*this)[15]*x28 - (*this)[15]*x43 + (*this)[1]*(*this)[2]*x25 - (*this)[1]*x21 + (*this)[1]*x29 - (*this)[2]*x23 + (*this)[2]*x38 + (*this)[3]*x54 + (*this)[4]*x33 + (*this)[4]*x37 + (*this)[4]*x56 - (*this)[5]*x18 + (*this)[6]*x17 + (*this)[6]*x27 - (*this)[6]*x30 - (*this)[7]*x47 - (*this)[8]*x24 + (*this)[8]*x36;
    res[11]=(*this)[11]*x58 + (*this)[12]*x61 + (*this)[13]*x62 - (*this)[15]*x66 + (*this)[4]*x60 + (*this)[7]*x63 + (*this)[8]*x64 + (*this)[9]*x65;
    res[12]=(*this)[11]*x61 + (*this)[12]*x58 - (*this)[13]*x60 + (*this)[15]*x64 - (*this)[4]*x62 + (*this)[7]*x65 - (*this)[8]*x66 + (*this)[9]*x63;
    res[13]=(*this)[11]*x62 + (*this)[12]*x60 + (*this)[13]*x58 + (*this)[15]*x65 + (*this)[4]*x61 - (*this)[7]*x64 - (*this)[8]*x63 - (*this)[9]*x66;
    res[14]=(*this)[11]*(*this)[7]*x59 - (*this)[12]*(*this)[9]*x59 + (*this)[13]*(*this)[8]*x59 + (*this)[14]*x58 + (*this)[15]*(*this)[4]*x59 - (*this)[1]*x60 + (*this)[2]*x62 - (*this)[3]*x61;
    res[15]=-(*this)[11]*x66 - (*this)[12]*x64 - (*this)[13]*x65 + (*this)[15]*x58 - (*this)[4]*x63 - (*this)[7]*x60 - (*this)[8]*x61 - (*this)[9]*x62;
    return res;
};

//-----------------------------------
// (R130MV, R130B0) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[4]*B[0];
    res[5]=A[5]*B[0];
    res[6]=A[6]*B[0];
    res[7]=A[7]*B[0];
    res[8]=A[8]*B[0];
    res[9]=A[9]*B[0];
    res[10]=A[10]*B[0];
    res[11]=A[11]*B[0];
    res[12]=A[12]*B[0];
    res[13]=A[13]*B[0];
    res[14]=A[14]*B[0];
    res[15]=A[15]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator/(const R130MV<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    res[6]=A[6]*x0;
    res[7]=A[7]*x0;
    res[8]=A[8]*x0;
    res[9]=A[9]*x0;
    res[10]=A[10]*x0;
    res[11]=A[11]*x0;
    res[12]=A[12]*x0;
    res[13]=A[13]*x0;
    res[14]=A[14]*x0;
    res[15]=A[15]*x0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[4]*B[0];
    res[5]=A[5]*B[0];
    res[6]=A[6]*B[0];
    res[7]=A[7]*B[0];
    res[8]=A[8]*B[0];
    res[9]=A[9]*B[0];
    res[10]=A[10]*B[0];
    res[11]=A[11]*B[0];
    res[12]=A[12]*B[0];
    res[13]=A[13]*B[0];
    res[14]=A[14]*B[0];
    res[15]=A[15]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B0<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[10], 2) + std::pow((*this)[11], 2) - std::pow((*this)[12], 2) - std::pow((*this)[13], 2) - std::pow((*this)[14], 2) - std::pow((*this)[15], 2) - std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) - std::pow((*this)[5], 2) - std::pow((*this)[6], 2) - std::pow((*this)[7], 2) + std::pow((*this)[8], 2) + std::pow((*this)[9], 2));
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
    res[11]=x0*((*this)[0]*(*this)[11] + (*this)[10]*(*this)[4] + (*this)[12]*(*this)[5] + (*this)[13]*(*this)[6] + (*this)[14]*(*this)[7] - (*this)[15]*(*this)[1] + (*this)[2]*(*this)[8] + (*this)[3]*(*this)[9]);
    res[12]=x0*((*this)[0]*(*this)[12] - (*this)[10]*(*this)[13] + (*this)[11]*(*this)[5] + (*this)[14]*(*this)[9] + (*this)[15]*(*this)[2] - (*this)[1]*(*this)[8] + (*this)[3]*(*this)[7] - (*this)[4]*(*this)[6]);
    res[13]=x0*((*this)[0]*(*this)[13] + (*this)[10]*(*this)[12] + (*this)[11]*(*this)[6] - (*this)[14]*(*this)[8] + (*this)[15]*(*this)[3] - (*this)[1]*(*this)[9] - (*this)[2]*(*this)[7] + (*this)[4]*(*this)[5]);
    res[14]=x0*((*this)[0]*(*this)[14] - (*this)[10]*(*this)[1] + (*this)[11]*(*this)[7] - (*this)[12]*(*this)[9] + (*this)[13]*(*this)[8] + (*this)[15]*(*this)[4] + (*this)[2]*(*this)[6] - (*this)[3]*(*this)[5]);
    res[15]=x0*((*this)[0]*(*this)[15] - (*this)[10]*(*this)[7] - (*this)[11]*(*this)[1] - (*this)[12]*(*this)[2] - (*this)[13]*(*this)[3] - (*this)[14]*(*this)[4] - (*this)[5]*(*this)[8] - (*this)[6]*(*this)[9]);
    return res;
};

//-----------------------------------
// (R130MV, R130B1) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    res[4]=A[4] + B[3];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    res[4]=A[4] - B[3];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[1]*B[0] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[1]=A[0]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3];
    res[2]=A[0]*B[1] + A[10]*B[2] + A[5]*B[0] - A[9]*B[3];
    res[3]=A[0]*B[2] - A[10]*B[1] + A[6]*B[0] + A[8]*B[3];
    res[4]=A[0]*B[3] + A[7]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[5]=A[13]*B[3] - A[14]*B[2] - A[1]*B[1] + A[2]*B[0];
    res[6]=-A[12]*B[3] + A[14]*B[1] - A[1]*B[2] + A[3]*B[0];
    res[7]=A[12]*B[2] - A[13]*B[1] - A[1]*B[3] + A[4]*B[0];
    res[8]=A[11]*B[1] + A[12]*B[0] - A[3]*B[3] + A[4]*B[2];
    res[9]=A[11]*B[2] + A[13]*B[0] + A[2]*B[3] - A[4]*B[1];
    res[10]=A[11]*B[3] + A[14]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[11]=-A[10]*B[3] - A[15]*B[0] - A[8]*B[1] - A[9]*B[2];
    res[12]=A[15]*B[1] + A[6]*B[3] - A[7]*B[2] + A[8]*B[0];
    res[13]=A[15]*B[2] - A[5]*B[3] + A[7]*B[1] + A[9]*B[0];
    res[14]=A[10]*B[0] + A[15]*B[3] + A[5]*B[2] - A[6]*B[1];
    res[15]=-A[11]*B[0] - A[12]*B[1] - A[13]*B[2] - A[14]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[1]*B[0] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=A[13]*B[3] - A[14]*B[2];
    res[6]=-A[12]*B[3] + A[14]*B[1];
    res[7]=A[12]*B[2] - A[13]*B[1];
    res[8]=A[11]*B[1] + A[12]*B[0];
    res[9]=A[11]*B[2] + A[13]*B[0];
    res[10]=A[11]*B[3] + A[14]*B[0];
    res[11]=-A[10]*B[3] - A[8]*B[1] - A[9]*B[2];
    res[12]=A[6]*B[3] - A[7]*B[2] + A[8]*B[0];
    res[13]=-A[5]*B[3] + A[7]*B[1] + A[9]*B[0];
    res[14]=A[10]*B[0] + A[5]*B[2] - A[6]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[5]*B[1] + A[6]*B[2] + A[7]*B[3];
    res[2]=A[10]*B[2] + A[5]*B[0] - A[9]*B[3];
    res[3]=-A[10]*B[1] + A[6]*B[0] + A[8]*B[3];
    res[4]=A[7]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[5]=-A[1]*B[1] + A[2]*B[0];
    res[6]=-A[1]*B[2] + A[3]*B[0];
    res[7]=-A[1]*B[3] + A[4]*B[0];
    res[8]=-A[3]*B[3] + A[4]*B[2];
    res[9]=A[2]*B[3] - A[4]*B[1];
    res[10]=-A[2]*B[2] + A[3]*B[1];
    res[11]=-A[15]*B[0];
    res[12]=A[15]*B[1];
    res[13]=A[15]*B[2];
    res[14]=A[15]*B[3];
    res[15]=-A[11]*B[0] - A[12]*B[1] - A[13]*B[2] - A[14]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B1<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[10], 2);
    T x2 = std::pow((*this)[15], 2);
    T x3 = std::pow((*this)[5], 2);
    T x4 = std::pow((*this)[6], 2);
    T x5 = std::pow((*this)[7], 2);
    T x6 = std::pow((*this)[8], 2);
    T x7 = std::pow((*this)[9], 2);
    T x8 = std::pow((*this)[11], 2);
    T x9 = std::pow((*this)[12], 2);
    T x10 = std::pow((*this)[13], 2);
    T x11 = std::pow((*this)[14], 2);
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[2], 2);
    T x14 = std::pow((*this)[3], 2);
    T x15 = std::pow((*this)[4], 2);
    T x16 = 2.0*A[2];
    T x17 = (*this)[5]*x16;
    T x18 = (*this)[12]*x16;
    T x19 = 2.0*(*this)[2];
    T x20 = (*this)[13]*A[3];
    T x21 = 2.0*A[1];
    T x22 = (*this)[14]*(*this)[3];
    T x23 = (*this)[1]*x19;
    T x24 = (*this)[3]*x16;
    T x25 = 2.0*A[3];
    T x26 = (*this)[4]*x25;
    T x27 = (*this)[8]*x25;
    T x28 = (*this)[7]*(*this)[9];
    T x29 = (*this)[0]*(*this)[5];
    T x30 = (*this)[0]*x16;
    T x31 = (*this)[0]*x25;
    T x32 = (*this)[15]*x25;
    T x33 = (*this)[10]*(*this)[6];
    T x34 = (*this)[11]*(*this)[12];
    T x35 = (*this)[11]*x16;
    T x36 = (*this)[11]*x25;
    T x37 = (*this)[12]*x25;
    T x38 = (*this)[13]*(*this)[4];
    T x39 = (*this)[14]*x16;
    T x40 = (*this)[15]*(*this)[8];
    T x41 = (*this)[15]*x16;
    T x42 = (*this)[5]*x25;
    T x43 = (*this)[8]*x16;
    T x44 = 2.0*A[0];
    T x45 = 2.0*x20;
    T x46 = A[3]*x19;
    T x47 = (*this)[10]*x21;
    T x48 = (*this)[10]*(*this)[9];
    T x49 = (*this)[11]*x44;
    T x50 = (*this)[11]*x21;
    T x51 = (*this)[12]*x21;
    T x52 = (*this)[12]*x44;
    T x53 = (*this)[1]*x21;
    T x54 = A[1]*x19;
    T x55 = (*this)[6]*x21;
    T x56 = (*this)[6]*(*this)[7];
    T x57 = (*this)[7]*x44;
    T x58 = (*this)[9]*x21;
    T x59 = (*this)[6]*x44;
    T x60 = (*this)[5]*x44;
    T x61 = A[0]*x19;
    T x62 = (*this)[7]*x21;
    T x63 = (*this)[15]*x44;
    T x64 = (*this)[1]*x44;
    T x65 = (*this)[3]*x44;
    T x66 = (*this)[9]*x25;
    T x67 = (*this)[6]*x16;
    T x68 = (*this)[10]*x16;
    T x69 = (*this)[4]*x44;
    T x70 = (*this)[14]*x21;
    T x71 = (*this)[13]*x21;
    T x72 = (*this)[9]*x16;
    T x73 = (*this)[3]*x21;
    T x74 = (*this)[3]*x25;
    T x75 = (*this)[4]*x21;
    T x76 = (*this)[7]*x16;
    T x77 = (*this)[14]*x25;
    T x78 = (*this)[14]*x44;
    T x79 = (*this)[13]*x44;
    res[0]=0;
    res[1]=(*this)[10]*x17 - (*this)[10]*x32 - (*this)[13]*x35 - (*this)[14]*x36 + (*this)[1]*x24 + (*this)[1]*x26 - (*this)[2]*x39 - (*this)[3]*x37 + (*this)[4]*x18 + (*this)[6]*x27 - (*this)[6]*x30 - (*this)[7]*x31 - (*this)[7]*x43 - (*this)[9]*x41 - (*this)[9]*x42 + A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - A[0]*x9 + A[1]*x23 + x19*x20 + x21*x22 + x21*x28 - x21*x29 - x21*x33 - x21*x34 - x21*x38 - x21*x40;
    res[2]=(*this)[10]*x27 - (*this)[10]*x30 + (*this)[13]*x18 + (*this)[14]*x37 - (*this)[1]*x39 + (*this)[1]*x45 + (*this)[2]*x24 + (*this)[3]*x36 - (*this)[4]*x35 + (*this)[4]*x46 + (*this)[6]*x17 - (*this)[6]*x32 + (*this)[7]*x41 + (*this)[7]*x42 + (*this)[9]*x31 + (*this)[9]*x43 - A[0]*x23 + A[1]*x0 - A[1]*x1 - A[1]*x10 - A[1]*x11 + A[1]*x12 + A[1]*x13 - A[1]*x14 - A[1]*x15 + A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 + A[1]*x8 + A[1]*x9 + x22*x44 - x28*x44 - x29*x44 + x33*x44 + x34*x44 - x38*x44 - x40*x44;
    res[3]=-(*this)[0]*x27 + (*this)[0]*x47 - (*this)[0]*x59 - (*this)[10]*x60 - (*this)[11]*x46 + (*this)[13]*x49 + (*this)[13]*x51 + (*this)[14]*x45 + (*this)[14]*x53 - (*this)[14]*x61 - (*this)[15]*x62 - (*this)[1]*x37 + (*this)[3]*x26 + (*this)[3]*x54 - (*this)[3]*x64 + (*this)[4]*x50 + (*this)[4]*x52 + (*this)[5]*x32 + (*this)[5]*x55 + (*this)[8]*x57 + (*this)[8]*x58 - (*this)[9]*x63 + A[2]*x0 - A[2]*x1 + A[2]*x10 - A[2]*x11 + A[2]*x12 - A[2]*x13 + A[2]*x14 - A[2]*x15 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + A[2]*x8 - A[2]*x9 + x25*x48 + x25*x56;
    res[4]=-(*this)[0]*x57 - (*this)[0]*x58 - (*this)[10]*x63 + (*this)[13]*x39 - (*this)[13]*x53 + (*this)[13]*x61 + (*this)[14]*x49 + (*this)[14]*x51 - (*this)[15]*x17 + (*this)[15]*x55 + (*this)[1]*x18 + (*this)[2]*x35 - (*this)[3]*x50 - (*this)[3]*x52 + (*this)[4]*x24 + (*this)[4]*x54 - (*this)[4]*x64 + (*this)[5]*x62 + (*this)[8]*x30 + (*this)[8]*x47 - (*this)[8]*x59 + (*this)[9]*x60 + A[3]*x0 + A[3]*x1 - A[3]*x10 + A[3]*x11 + A[3]*x12 - A[3]*x13 - A[3]*x14 + A[3]*x15 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7 + A[3]*x8 - A[3]*x9 + x16*x48 + x16*x56;
    res[5]=(*this)[0]*x45 + (*this)[0]*x53 - (*this)[0]*x61 + (*this)[10]*x37 + (*this)[10]*x65 + (*this)[13]*x43 - (*this)[13]*x57 - (*this)[13]*x58 + (*this)[14]*x27 - (*this)[14]*x30 - (*this)[14]*x47 + (*this)[14]*x59 - (*this)[15]*x50 - (*this)[15]*x52 - (*this)[1]*x60 + (*this)[1]*x66 - (*this)[1]*x68 + (*this)[2]*x67 + (*this)[3]*x17 - (*this)[3]*x32 - (*this)[3]*x55 + (*this)[4]*x41 - (*this)[4]*x62 + (*this)[5]*x26 + (*this)[5]*x54 + (*this)[6]*x36 - (*this)[7]*x35 + (*this)[7]*x46 + (*this)[8]*x49 + (*this)[8]*x51 + (*this)[9]*x18 - (*this)[9]*x69;
    res[6]=-(*this)[0]*x65 + (*this)[0]*x70 - (*this)[10]*x39 + (*this)[10]*x45 - (*this)[10]*x61 - (*this)[12]*x31 - (*this)[13]*x63 + (*this)[13]*x72 - (*this)[14]*x60 + (*this)[14]*x66 - (*this)[15]*x35 + (*this)[15]*x46 - (*this)[15]*x75 - (*this)[1]*x27 + (*this)[1]*x30 + (*this)[1]*x47 - (*this)[1]*x59 - (*this)[2]*x17 - (*this)[4]*x76 - (*this)[5]*x36 + (*this)[5]*x73 + (*this)[6]*x24 + (*this)[6]*x26 + (*this)[6]*x54 + (*this)[7]*x50 + (*this)[7]*x52 + (*this)[7]*x74 - (*this)[8]*x18 + (*this)[8]*x69 + (*this)[8]*x71 + (*this)[9]*x49 + (*this)[9]*x51;
    res[7]=(*this)[0]*x18 - (*this)[0]*x69 - (*this)[0]*x71 + (*this)[10]*x49 + (*this)[10]*x77 + (*this)[11]*x17 - (*this)[11]*x32 - (*this)[12]*x27 + (*this)[12]*x47 + (*this)[13]*x60 + (*this)[13]*x68 - (*this)[14]*x63 + (*this)[15]*x73 + (*this)[1]*x31 + (*this)[1]*x43 - (*this)[1]*x57 - (*this)[2]*x41 + (*this)[4]*x67 - (*this)[5]*x46 + (*this)[5]*x75 - (*this)[6]*x50 - (*this)[6]*x52 - (*this)[6]*x74 + (*this)[7]*x24 + (*this)[7]*x26 + (*this)[7]*x54 - (*this)[8]*x65 + (*this)[8]*x70 + (*this)[9]*x39 - (*this)[9]*x45 - (*this)[9]*x53 + (*this)[9]*x61;
    res[8]=(*this)[0]*x50 + (*this)[0]*x52 - (*this)[10]*x35 + (*this)[10]*x46 - (*this)[10]*x79 - (*this)[13]*x17 + (*this)[13]*x55 - (*this)[14]*x42 + (*this)[14]*x62 - (*this)[15]*x39 + (*this)[15]*x45 + (*this)[15]*x53 - (*this)[15]*x61 - (*this)[1]*(*this)[6]*x25 + (*this)[1]*x76 + (*this)[2]*x72 + (*this)[3]*x31 - (*this)[3]*x57 - (*this)[3]*x58 - (*this)[4]*x30 - (*this)[4]*x47 + (*this)[4]*x59 - (*this)[5]*x49 - (*this)[5]*x51 - (*this)[6]*x18 - (*this)[7]*x37 + (*this)[8]*x24 + (*this)[8]*x26 + (*this)[8]*x54 - (*this)[8]*x64 + (*this)[9]*x36 + (*this)[9]*x78;
    res[9]=-(*this)[0]*x46 + (*this)[0]*x75 + (*this)[0]*x79 + (*this)[10]*x52 + (*this)[10]*x74 - (*this)[11]*x27 + (*this)[11]*x30 + (*this)[11]*x47 + (*this)[12]*x17 - (*this)[12]*x32 - (*this)[13]*x67 + (*this)[15]*x70 + (*this)[1]*x41 + (*this)[1]*x42 - (*this)[2]*x43 - (*this)[3]*x63 - (*this)[4]*x60 - (*this)[4]*x68 - (*this)[5]*x71 - (*this)[6]*x49 - (*this)[6]*x51 - (*this)[6]*x77 + (*this)[7]*x39 - (*this)[7]*x45 - (*this)[7]*x53 + (*this)[7]*x61 + (*this)[8]*x73 - (*this)[8]*x78 + (*this)[9]*x24 + (*this)[9]*x26 + (*this)[9]*x54 - (*this)[9]*x64;
    res[10]=-(*this)[0]*x73 + (*this)[0]*x78 + (*this)[10]*x24 + (*this)[10]*x26 + (*this)[10]*x54 - (*this)[10]*x64 + (*this)[11]*x31 - (*this)[13]*x76 + (*this)[15]*x18 - (*this)[15]*x71 - (*this)[1]*x17 + (*this)[1]*x32 + (*this)[2]*x30 + (*this)[3]*x60 - (*this)[3]*x66 - (*this)[4]*x63 + (*this)[4]*x72 + (*this)[5]*x37 - (*this)[5]*x70 - (*this)[6]*x39 + (*this)[6]*x45 + (*this)[6]*x53 - (*this)[6]*x61 - (*this)[7]*x49 - (*this)[7]*x51 - (*this)[7]*x77 + (*this)[8]*x35 - (*this)[8]*x46 + (*this)[8]*x75 + (*this)[8]*x79 - (*this)[9]*x50 - (*this)[9]*x52;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//--------------------------------------
// (R130MV, R130B1Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2] + B[0];
    res[3]=A[3] + B[1];
    res[4]=A[4] + B[2];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2] - B[0];
    res[3]=A[3] - B[1];
    res[4]=A[4] - B[2];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[1]=A[5]*B[0] + A[6]*B[1] + A[7]*B[2];
    res[2]=A[0]*B[0] + A[10]*B[1] - A[9]*B[2];
    res[3]=A[0]*B[1] - A[10]*B[0] + A[8]*B[2];
    res[4]=A[0]*B[2] - A[8]*B[1] + A[9]*B[0];
    res[5]=A[13]*B[2] - A[14]*B[1] - A[1]*B[0];
    res[6]=-A[12]*B[2] + A[14]*B[0] - A[1]*B[1];
    res[7]=A[12]*B[1] - A[13]*B[0] - A[1]*B[2];
    res[8]=A[11]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[9]=A[11]*B[1] + A[2]*B[2] - A[4]*B[0];
    res[10]=A[11]*B[2] - A[2]*B[1] + A[3]*B[0];
    res[11]=-A[10]*B[2] - A[8]*B[0] - A[9]*B[1];
    res[12]=A[15]*B[0] + A[6]*B[2] - A[7]*B[1];
    res[13]=A[15]*B[1] - A[5]*B[2] + A[7]*B[0];
    res[14]=A[15]*B[2] + A[5]*B[1] - A[6]*B[0];
    res[15]=-A[12]*B[0] - A[13]*B[1] - A[14]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[0]*B[1];
    res[4]=A[0]*B[2];
    res[5]=A[13]*B[2] - A[14]*B[1];
    res[6]=-A[12]*B[2] + A[14]*B[0];
    res[7]=A[12]*B[1] - A[13]*B[0];
    res[8]=A[11]*B[0];
    res[9]=A[11]*B[1];
    res[10]=A[11]*B[2];
    res[11]=-A[10]*B[2] - A[8]*B[0] - A[9]*B[1];
    res[12]=A[6]*B[2] - A[7]*B[1];
    res[13]=-A[5]*B[2] + A[7]*B[0];
    res[14]=A[5]*B[1] - A[6]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[5]*B[0] + A[6]*B[1] + A[7]*B[2];
    res[2]=A[10]*B[1] - A[9]*B[2];
    res[3]=-A[10]*B[0] + A[8]*B[2];
    res[4]=-A[8]*B[1] + A[9]*B[0];
    res[5]=-A[1]*B[0];
    res[6]=-A[1]*B[1];
    res[7]=-A[1]*B[2];
    res[8]=-A[3]*B[2] + A[4]*B[1];
    res[9]=A[2]*B[2] - A[4]*B[0];
    res[10]=-A[2]*B[1] + A[3]*B[0];
    res[11]=0;
    res[12]=A[15]*B[0];
    res[13]=A[15]*B[1];
    res[14]=A[15]*B[2];
    res[15]=-A[12]*B[0] - A[13]*B[1] - A[14]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*A[1];
    T x1 = (*this)[5]*x0;
    T x2 = (*this)[12]*x0;
    T x3 = 2.0*(*this)[2];
    T x4 = (*this)[13]*A[2];
    T x5 = 2.0*A[0];
    T x6 = (*this)[14]*x5;
    T x7 = A[0]*x3;
    T x8 = (*this)[3]*x0;
    T x9 = 2.0*A[2];
    T x10 = (*this)[4]*x9;
    T x11 = (*this)[8]*x9;
    T x12 = (*this)[9]*x5;
    T x13 = (*this)[0]*x5;
    T x14 = (*this)[0]*x0;
    T x15 = (*this)[0]*x9;
    T x16 = (*this)[15]*x9;
    T x17 = (*this)[6]*x5;
    T x18 = (*this)[11]*x5;
    T x19 = (*this)[11]*x0;
    T x20 = (*this)[11]*x9;
    T x21 = (*this)[12]*x9;
    T x22 = (*this)[13]*x5;
    T x23 = (*this)[14]*x0;
    T x24 = (*this)[15]*x5;
    T x25 = (*this)[15]*x0;
    T x26 = (*this)[5]*x9;
    T x27 = (*this)[8]*x0;
    T x28 = std::pow((*this)[0], 2);
    T x29 = std::pow((*this)[11], 2);
    T x30 = std::pow((*this)[12], 2);
    T x31 = std::pow((*this)[15], 2);
    T x32 = std::pow((*this)[1], 2);
    T x33 = std::pow((*this)[2], 2);
    T x34 = std::pow((*this)[5], 2);
    T x35 = std::pow((*this)[8], 2);
    T x36 = std::pow((*this)[10], 2);
    T x37 = std::pow((*this)[13], 2);
    T x38 = std::pow((*this)[14], 2);
    T x39 = std::pow((*this)[3], 2);
    T x40 = std::pow((*this)[4], 2);
    T x41 = std::pow((*this)[6], 2);
    T x42 = std::pow((*this)[7], 2);
    T x43 = std::pow((*this)[9], 2);
    T x44 = 2.0*x4;
    T x45 = A[2]*x3;
    T x46 = (*this)[10]*(*this)[9];
    T x47 = (*this)[6]*(*this)[7];
    T x48 = (*this)[8]*x5;
    T x49 = (*this)[7]*x5;
    T x50 = (*this)[9]*x9;
    T x51 = (*this)[6]*x0;
    T x52 = (*this)[10]*(*this)[1];
    T x53 = (*this)[13]*x0;
    T x54 = (*this)[5]*x5;
    T x55 = (*this)[3]*x9;
    T x56 = (*this)[7]*x0;
    T x57 = (*this)[10]*x5;
    T x58 = (*this)[14]*x9;
    T x59 = (*this)[9]*x0;
    res[0]=0;
    res[1]=(*this)[10]*x1 - (*this)[10]*x16 - (*this)[10]*x17 - (*this)[12]*x18 - (*this)[13]*x19 - (*this)[14]*x20 + (*this)[1]*x10 + (*this)[1]*x7 + (*this)[1]*x8 - (*this)[2]*x23 - (*this)[3]*x21 + (*this)[3]*x6 + (*this)[4]*x2 - (*this)[4]*x22 - (*this)[5]*x13 + (*this)[6]*x11 - (*this)[6]*x14 + (*this)[7]*x12 - (*this)[7]*x15 - (*this)[7]*x27 - (*this)[8]*x24 - (*this)[9]*x25 - (*this)[9]*x26 + x3*x4;
    res[2]=(*this)[10]*x11 - (*this)[10]*x14 + (*this)[13]*x2 + (*this)[14]*x21 - (*this)[1]*x23 + (*this)[1]*x44 + (*this)[2]*x8 + (*this)[3]*x20 - (*this)[4]*x19 + (*this)[4]*x45 + (*this)[6]*x1 - (*this)[6]*x16 + (*this)[7]*x25 + (*this)[7]*x26 + (*this)[9]*x15 + (*this)[9]*x27 + A[0]*x28 + A[0]*x29 + A[0]*x30 + A[0]*x31 + A[0]*x32 + A[0]*x33 + A[0]*x34 + A[0]*x35 - A[0]*x36 - A[0]*x37 - A[0]*x38 - A[0]*x39 - A[0]*x40 - A[0]*x41 - A[0]*x42 - A[0]*x43;
    res[3]=-(*this)[0]*x11 + (*this)[10]*x13 - (*this)[11]*x45 + (*this)[12]*x22 + (*this)[14]*x44 - (*this)[1]*x21 + (*this)[1]*x6 + (*this)[3]*x10 + (*this)[3]*x7 + (*this)[4]*x18 + (*this)[5]*x16 + (*this)[5]*x17 - (*this)[7]*x24 + (*this)[8]*x12 + A[1]*x28 + A[1]*x29 - A[1]*x30 + A[1]*x31 + A[1]*x32 - A[1]*x33 - A[1]*x34 - A[1]*x35 - A[1]*x36 + A[1]*x37 - A[1]*x38 + A[1]*x39 - A[1]*x40 + A[1]*x41 - A[1]*x42 + A[1]*x43 + x46*x9 + x47*x9;
    res[4]=-(*this)[0]*x12 + (*this)[10]*x48 + (*this)[12]*x6 + (*this)[13]*x23 - (*this)[15]*x1 + (*this)[15]*x17 + (*this)[1]*x2 - (*this)[1]*x22 + (*this)[2]*x19 - (*this)[3]*x18 + (*this)[4]*x7 + (*this)[4]*x8 + (*this)[5]*x49 + (*this)[8]*x14 + A[2]*x28 + A[2]*x29 - A[2]*x30 + A[2]*x31 + A[2]*x32 - A[2]*x33 - A[2]*x34 - A[2]*x35 + A[2]*x36 - A[2]*x37 + A[2]*x38 - A[2]*x39 + A[2]*x40 - A[2]*x41 + A[2]*x42 - A[2]*x43 + x0*x46 + x0*x47;
    res[5]=(*this)[0]*x44 + (*this)[10]*x21 - (*this)[10]*x6 + (*this)[12]*x48 - (*this)[13]*x12 + (*this)[13]*x27 + (*this)[14]*x11 - (*this)[14]*x14 - (*this)[15]*x18 + (*this)[1]*x13 + (*this)[1]*x50 + (*this)[2]*x51 + (*this)[3]*x1 - (*this)[3]*x16 - (*this)[3]*x17 + (*this)[4]*x25 - (*this)[4]*x49 + (*this)[5]*x10 + (*this)[5]*x7 + (*this)[6]*x20 - (*this)[7]*x19 + (*this)[7]*x45 + (*this)[9]*x2 - x0*x52;
    res[6]=(*this)[0]*x6 - (*this)[10]*x23 + (*this)[10]*x44 + (*this)[12]*x12 - (*this)[12]*x15 + (*this)[14]*x50 - (*this)[15]*x19 + (*this)[15]*x45 - (*this)[1]*x11 + (*this)[1]*x14 - (*this)[2]*x1 + (*this)[3]*x54 - (*this)[4]*x24 - (*this)[4]*x56 - (*this)[5]*x20 + (*this)[6]*x10 + (*this)[6]*x7 + (*this)[6]*x8 + (*this)[7]*x18 + (*this)[7]*x55 - (*this)[8]*x2 + (*this)[8]*x22 + (*this)[9]*x53 + x5*x52;
    res[7]=(*this)[0]*x2 + (*this)[10]*x53 + (*this)[10]*x58 + (*this)[11]*x1 - (*this)[11]*x16 - (*this)[11]*x17 - (*this)[12]*x11 + (*this)[12]*x57 - (*this)[13]*x13 - (*this)[1]*x12 + (*this)[1]*x15 + (*this)[1]*x27 - (*this)[2]*x25 + (*this)[3]*x24 + (*this)[4]*x51 + (*this)[4]*x54 - (*this)[5]*x45 - (*this)[6]*x55 + (*this)[7]*x10 + (*this)[7]*x7 + (*this)[7]*x8 + (*this)[8]*x6 + (*this)[9]*x23 - (*this)[9]*x44;
    res[8]=-(*this)[10]*x19 + (*this)[10]*x45 + (*this)[11]*x13 - (*this)[12]*x54 - (*this)[13]*x1 + (*this)[13]*x17 - (*this)[14]*x26 - (*this)[15]*x23 + (*this)[15]*x44 - (*this)[1]*(*this)[6]*x9 + (*this)[1]*x24 + (*this)[1]*x56 + (*this)[2]*x59 - (*this)[3]*x12 + (*this)[3]*x15 - (*this)[4]*x14 - (*this)[4]*x57 - (*this)[6]*x2 - (*this)[7]*x21 + (*this)[7]*x6 + (*this)[8]*x10 + (*this)[8]*x7 + (*this)[8]*x8 + (*this)[9]*x20;
    res[9]=-(*this)[0]*x45 - (*this)[10]*(*this)[4]*x0 + (*this)[10]*x18 + (*this)[10]*x55 - (*this)[11]*x11 + (*this)[11]*x14 + (*this)[12]*x1 - (*this)[12]*x16 - (*this)[12]*x17 - (*this)[13]*x51 + (*this)[15]*x6 + (*this)[1]*x25 + (*this)[1]*x26 - (*this)[1]*x49 - (*this)[2]*x27 + (*this)[3]*x48 + (*this)[4]*x13 - (*this)[5]*x22 - (*this)[6]*x58 + (*this)[7]*x23 - (*this)[7]*x44 + (*this)[9]*x10 + (*this)[9]*x7 + (*this)[9]*x8;
    res[10]=(*this)[10]*x10 + (*this)[10]*x7 + (*this)[10]*x8 - (*this)[11]*x12 + (*this)[11]*x15 - (*this)[12]*x49 + (*this)[15]*x2 - (*this)[15]*x22 - (*this)[1]*x1 + (*this)[1]*x16 + (*this)[1]*x17 + (*this)[2]*x14 - (*this)[3]*x13 - (*this)[3]*x50 + (*this)[4]*x48 + (*this)[4]*x59 + (*this)[5]*x21 - (*this)[5]*x6 - (*this)[6]*x23 + (*this)[6]*x44 - (*this)[7]*x53 - (*this)[7]*x58 + (*this)[8]*x19 - (*this)[8]*x45;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//--------------------------------------
// (R130MV, R130B1Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[1]*B[0];
    res[1]=A[0]*B[0];
    res[2]=A[5]*B[0];
    res[3]=A[6]*B[0];
    res[4]=A[7]*B[0];
    res[5]=A[2]*B[0];
    res[6]=A[3]*B[0];
    res[7]=A[4]*B[0];
    res[8]=A[12]*B[0];
    res[9]=A[13]*B[0];
    res[10]=A[14]*B[0];
    res[11]=-A[15]*B[0];
    res[12]=A[8]*B[0];
    res[13]=A[9]*B[0];
    res[14]=A[10]*B[0];
    res[15]=-A[11]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[1]*B[0];
    res[1]=A[0]*B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[12]*B[0];
    res[9]=A[13]*B[0];
    res[10]=A[14]*B[0];
    res[11]=0;
    res[12]=A[8]*B[0];
    res[13]=A[9]*B[0];
    res[14]=A[10]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[5]*B[0];
    res[3]=A[6]*B[0];
    res[4]=A[7]*B[0];
    res[5]=A[2]*B[0];
    res[6]=A[3]*B[0];
    res[7]=A[4]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[15]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-A[11]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*A[0];
    res[0]=0;
    res[1]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[10], 2) - std::pow((*this)[11], 2) - std::pow((*this)[12], 2) - std::pow((*this)[13], 2) - std::pow((*this)[14], 2) + std::pow((*this)[15], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2) - std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) + std::pow((*this)[7], 2) + std::pow((*this)[8], 2) + std::pow((*this)[9], 2));
    res[2]=x0*(-(*this)[0]*(*this)[5] + (*this)[10]*(*this)[6] + (*this)[11]*(*this)[12] - (*this)[13]*(*this)[4] + (*this)[14]*(*this)[3] - (*this)[15]*(*this)[8] - (*this)[1]*(*this)[2] - (*this)[7]*(*this)[9]);
    res[3]=x0*(-(*this)[0]*(*this)[6] - (*this)[10]*(*this)[5] + (*this)[11]*(*this)[13] + (*this)[12]*(*this)[4] - (*this)[14]*(*this)[2] - (*this)[15]*(*this)[9] - (*this)[1]*(*this)[3] + (*this)[7]*(*this)[8]);
    res[4]=x0*(-(*this)[0]*(*this)[7] - (*this)[10]*(*this)[15] + (*this)[11]*(*this)[14] - (*this)[12]*(*this)[3] + (*this)[13]*(*this)[2] - (*this)[1]*(*this)[4] + (*this)[5]*(*this)[9] - (*this)[6]*(*this)[8]);
    res[5]=x0*(-(*this)[0]*(*this)[2] + (*this)[10]*(*this)[3] + (*this)[11]*(*this)[8] - (*this)[12]*(*this)[15] - (*this)[13]*(*this)[7] + (*this)[14]*(*this)[6] - (*this)[1]*(*this)[5] - (*this)[4]*(*this)[9]);
    res[6]=x0*(-(*this)[0]*(*this)[3] - (*this)[10]*(*this)[2] + (*this)[11]*(*this)[9] + (*this)[12]*(*this)[7] - (*this)[13]*(*this)[15] - (*this)[14]*(*this)[5] - (*this)[1]*(*this)[6] + (*this)[4]*(*this)[8]);
    res[7]=x0*(-(*this)[0]*(*this)[4] + (*this)[10]*(*this)[11] - (*this)[12]*(*this)[6] + (*this)[13]*(*this)[5] - (*this)[14]*(*this)[15] - (*this)[1]*(*this)[7] + (*this)[2]*(*this)[9] - (*this)[3]*(*this)[8]);
    res[8]=x0*((*this)[0]*(*this)[12] - (*this)[10]*(*this)[13] - (*this)[11]*(*this)[5] + (*this)[14]*(*this)[9] - (*this)[15]*(*this)[2] - (*this)[1]*(*this)[8] - (*this)[3]*(*this)[7] + (*this)[4]*(*this)[6]);
    res[9]=x0*((*this)[0]*(*this)[13] + (*this)[10]*(*this)[12] - (*this)[11]*(*this)[6] - (*this)[14]*(*this)[8] - (*this)[15]*(*this)[3] - (*this)[1]*(*this)[9] + (*this)[2]*(*this)[7] - (*this)[4]*(*this)[5]);
    res[10]=x0*((*this)[0]*(*this)[14] - (*this)[10]*(*this)[1] - (*this)[11]*(*this)[7] - (*this)[12]*(*this)[9] + (*this)[13]*(*this)[8] - (*this)[15]*(*this)[4] - (*this)[2]*(*this)[6] + (*this)[3]*(*this)[5]);
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//-----------------------------------
// (R130MV, R130B2) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] + B[0];
    res[6]=A[6] + B[1];
    res[7]=A[7] + B[2];
    res[8]=A[8] + B[3];
    res[9]=A[9] + B[4];
    res[10]=A[10] + B[5];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] - B[0];
    res[6]=A[6] - B[1];
    res[7]=A[7] - B[2];
    res[8]=A[8] - B[3];
    res[9]=A[9] - B[4];
    res[10]=A[10] - B[5];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=-A[10]*B[5] + A[5]*B[0] + A[6]*B[1] + A[7]*B[2] - A[8]*B[3] - A[9]*B[4];
    res[1]=-A[12]*B[3] - A[13]*B[4] - A[14]*B[5] - A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[2]=A[11]*B[3] + A[13]*B[2] - A[14]*B[1] - A[1]*B[0] - A[3]*B[5] + A[4]*B[4];
    res[3]=A[11]*B[4] - A[12]*B[2] + A[14]*B[0] - A[1]*B[1] + A[2]*B[5] - A[4]*B[3];
    res[4]=A[11]*B[5] + A[12]*B[1] - A[13]*B[0] - A[1]*B[2] - A[2]*B[4] + A[3]*B[3];
    res[5]=A[0]*B[0] + A[10]*B[1] - A[15]*B[3] - A[6]*B[5] + A[7]*B[4] - A[9]*B[2];
    res[6]=A[0]*B[1] - A[10]*B[0] - A[15]*B[4] + A[5]*B[5] - A[7]*B[3] + A[8]*B[2];
    res[7]=A[0]*B[2] - A[15]*B[5] - A[5]*B[4] + A[6]*B[3] - A[8]*B[1] + A[9]*B[0];
    res[8]=A[0]*B[3] + A[10]*B[4] + A[15]*B[0] + A[6]*B[2] - A[7]*B[1] - A[9]*B[5];
    res[9]=A[0]*B[4] - A[10]*B[3] + A[15]*B[1] - A[5]*B[2] + A[7]*B[0] + A[8]*B[5];
    res[10]=A[0]*B[5] + A[15]*B[2] + A[5]*B[1] - A[6]*B[0] - A[8]*B[4] + A[9]*B[3];
    res[11]=A[12]*B[0] + A[13]*B[1] + A[14]*B[2] - A[2]*B[3] - A[3]*B[4] - A[4]*B[5];
    res[12]=A[11]*B[0] - A[13]*B[5] + A[14]*B[4] + A[1]*B[3] - A[3]*B[2] + A[4]*B[1];
    res[13]=A[11]*B[1] + A[12]*B[5] - A[14]*B[3] + A[1]*B[4] + A[2]*B[2] - A[4]*B[0];
    res[14]=A[11]*B[2] - A[12]*B[4] + A[13]*B[3] + A[1]*B[5] - A[2]*B[1] + A[3]*B[0];
    res[15]=A[10]*B[2] + A[5]*B[3] + A[6]*B[4] + A[7]*B[5] + A[8]*B[0] + A[9]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=-A[10]*B[5] + A[5]*B[0] + A[6]*B[1] + A[7]*B[2] - A[8]*B[3] - A[9]*B[4];
    res[1]=-A[12]*B[3] - A[13]*B[4] - A[14]*B[5];
    res[2]=A[11]*B[3] + A[13]*B[2] - A[14]*B[1];
    res[3]=A[11]*B[4] - A[12]*B[2] + A[14]*B[0];
    res[4]=A[11]*B[5] + A[12]*B[1] - A[13]*B[0];
    res[5]=A[0]*B[0] - A[15]*B[3];
    res[6]=A[0]*B[1] - A[15]*B[4];
    res[7]=A[0]*B[2] - A[15]*B[5];
    res[8]=A[0]*B[3] + A[15]*B[0];
    res[9]=A[0]*B[4] + A[15]*B[1];
    res[10]=A[0]*B[5] + A[15]*B[2];
    res[11]=-A[2]*B[3] - A[3]*B[4] - A[4]*B[5];
    res[12]=A[1]*B[3] - A[3]*B[2] + A[4]*B[1];
    res[13]=A[1]*B[4] + A[2]*B[2] - A[4]*B[0];
    res[14]=A[1]*B[5] - A[2]*B[1] + A[3]*B[0];
    res[15]=A[10]*B[2] + A[5]*B[3] + A[6]*B[4] + A[7]*B[5] + A[8]*B[0] + A[9]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[2]=-A[1]*B[0] - A[3]*B[5] + A[4]*B[4];
    res[3]=-A[1]*B[1] + A[2]*B[5] - A[4]*B[3];
    res[4]=-A[1]*B[2] - A[2]*B[4] + A[3]*B[3];
    res[5]=A[10]*B[1] - A[6]*B[5] + A[7]*B[4] - A[9]*B[2];
    res[6]=-A[10]*B[0] + A[5]*B[5] - A[7]*B[3] + A[8]*B[2];
    res[7]=-A[5]*B[4] + A[6]*B[3] - A[8]*B[1] + A[9]*B[0];
    res[8]=A[10]*B[4] + A[6]*B[2] - A[7]*B[1] - A[9]*B[5];
    res[9]=-A[10]*B[3] - A[5]*B[2] + A[7]*B[0] + A[8]*B[5];
    res[10]=A[5]*B[1] - A[6]*B[0] - A[8]*B[4] + A[9]*B[3];
    res[11]=A[12]*B[0] + A[13]*B[1] + A[14]*B[2];
    res[12]=A[11]*B[0] - A[13]*B[5] + A[14]*B[4];
    res[13]=A[11]*B[1] + A[12]*B[5] - A[14]*B[3];
    res[14]=A[11]*B[2] - A[12]*B[4] + A[13]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B2<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[0]*x0;
    T x2 = (*this)[3]*x0;
    T x3 = A[2]*x0;
    T x4 = 2.0*(*this)[10];
    T x5 = (*this)[12]*x4;
    T x6 = A[5]*x4;
    T x7 = (*this)[3]*A[0];
    T x8 = 2.0*(*this)[6];
    T x9 = (*this)[12]*x8;
    T x10 = 2.0*(*this)[13];
    T x11 = (*this)[7]*A[0];
    T x12 = (*this)[8]*x10;
    T x13 = 2.0*(*this)[14];
    T x14 = (*this)[5]*A[1];
    T x15 = A[3]*x13;
    T x16 = 2.0*(*this)[1];
    T x17 = (*this)[8]*A[3];
    T x18 = (*this)[9]*x16;
    T x19 = A[5]*x8;
    T x20 = 2.0*(*this)[9];
    T x21 = (*this)[2]*A[2];
    T x22 = 2.0*A[3];
    T x23 = (*this)[3]*(*this)[7];
    T x24 = 2.0*(*this)[4];
    T x25 = (*this)[5]*A[4];
    T x26 = (*this)[8]*x24;
    T x27 = A[3]*x0;
    T x28 = A[4]*x0;
    T x29 = (*this)[14]*x0;
    T x30 = (*this)[11]*x4;
    T x31 = A[3]*x4;
    T x32 = A[1]*x4;
    T x33 = (*this)[5]*x22;
    T x34 = (*this)[11]*x8;
    T x35 = 2.0*(*this)[11];
    T x36 = (*this)[7]*A[5];
    T x37 = A[0]*x35;
    T x38 = (*this)[11]*x20;
    T x39 = 2.0*(*this)[12];
    T x40 = (*this)[15]*A[0];
    T x41 = (*this)[7]*A[1];
    T x42 = (*this)[12]*x20;
    T x43 = (*this)[15]*x10;
    T x44 = (*this)[5]*A[2];
    T x45 = (*this)[15]*x13;
    T x46 = (*this)[14]*A[0];
    T x47 = (*this)[8]*x13;
    T x48 = (*this)[15]*(*this)[2];
    T x49 = 2.0*A[4];
    T x50 = (*this)[15]*(*this)[3];
    T x51 = (*this)[15]*x24;
    T x52 = (*this)[5]*A[0];
    T x53 = A[1]*x8;
    T x54 = (*this)[7]*x16;
    T x55 = (*this)[2]*x49;
    T x56 = 2.0*(*this)[3];
    T x57 = (*this)[5]*A[5];
    T x58 = (*this)[8]*A[2];
    T x59 = A[3]*x8;
    T x60 = (*this)[4]*x20;
    T x61 = (*this)[8]*x39;
    T x62 = 2.0*A[2];
    T x63 = 2.0*(*this)[2];
    T x64 = (*this)[2]*x20;
    T x65 = (*this)[8]*x49;
    T x66 = (*this)[9]*x10;
    T x67 = (*this)[15]*x16;
    T x68 = 2.0*x21;
    T x69 = (*this)[3]*x20;
    T x70 = A[1]*x0;
    T x71 = A[2]*x4;
    T x72 = (*this)[1]*A[0];
    T x73 = (*this)[15]*x35;
    T x74 = (*this)[15]*x39;
    T x75 = A[3]*x10;
    T x76 = A[4]*x8;
    T x77 = (*this)[9]*x13;
    T x78 = A[5]*x0;
    T x79 = A[4]*x4;
    T x80 = (*this)[8]*x35;
    T x81 = (*this)[7]*A[4];
    T x82 = (*this)[2]*A[0];
    T x83 = 2.0*x7;
    T x84 = A[2]*x8;
    T x85 = (*this)[12]*x22;
    T x86 = 2.0*A[1];
    T x87 = (*this)[8]*A[1];
    T x88 = (*this)[8]*A[5];
    T x89 = A[2]*x24;
    T x90 = std::pow((*this)[0], 2);
    T x91 = std::pow((*this)[12], 2);
    T x92 = std::pow((*this)[1], 2);
    T x93 = std::pow((*this)[3], 2);
    T x94 = std::pow((*this)[4], 2);
    T x95 = std::pow((*this)[6], 2);
    T x96 = std::pow((*this)[7], 2);
    T x97 = std::pow((*this)[8], 2);
    T x98 = std::pow((*this)[10], 2);
    T x99 = std::pow((*this)[11], 2);
    T x100 = std::pow((*this)[13], 2);
    T x101 = std::pow((*this)[14], 2);
    T x102 = std::pow((*this)[15], 2);
    T x103 = std::pow((*this)[2], 2);
    T x104 = std::pow((*this)[5], 2);
    T x105 = std::pow((*this)[9], 2);
    T x106 = A[5]*x10;
    T x107 = (*this)[11]*x16;
    T x108 = (*this)[11]*x24;
    T x109 = (*this)[12]*x10;
    T x110 = (*this)[12]*A[2];
    T x111 = (*this)[3]*A[4];
    T x112 = (*this)[12]*A[5];
    T x113 = A[2]*x10;
    T x114 = (*this)[2]*A[4];
    T x115 = A[5]*x13;
    T x116 = A[5]*x16;
    T x117 = 2.0*x17;
    T x118 = 2.0*(*this)[8];
    T x119 = A[4]*x13;
    T x120 = (*this)[3]*x35;
    T x121 = A[1]*x13;
    T x122 = 2.0*(*this)[15];
    T x123 = (*this)[15]*x20;
    T x124 = (*this)[4]*x16;
    T x125 = (*this)[2]*A[1];
    T x126 = 2.0*(*this)[7];
    T x127 = (*this)[8]*A[0];
    T x128 = (*this)[12]*x35;
    T x129 = A[0]*x13;
    T x130 = (*this)[12]*x16;
    T x131 = A[1]*x10;
    T x132 = A[4]*x10;
    T x133 = (*this)[2]*A[5];
    T x134 = (*this)[3]*x16;
    T x135 = (*this)[3]*A[1];
    T x136 = 2.0*(*this)[5];
    T x137 = (*this)[15]*(*this)[7];
    T x138 = (*this)[4]*x22;
    T x139 = A[2]*x13;
    T x140 = (*this)[3]*x22;
    res[0]=0;
    res[1]=-(*this)[11]*x33 - (*this)[12]*x27 - (*this)[13]*x28 - (*this)[13]*x31 - (*this)[1]*x53 + (*this)[1]*x6 + (*this)[2]*x1 + (*this)[2]*x19 - (*this)[2]*x32 + (*this)[4]*x3 - (*this)[4]*x59 - (*this)[7]*x55 - (*this)[8]*x37 + (*this)[9]*x15 - A[0]*x60 + A[1]*x2 + A[1]*x26 - A[1]*x38 - A[1]*x43 - A[2]*x30 - A[2]*x45 - A[2]*x54 + A[2]*x9 + A[4]*x18 - A[4]*x34 - A[4]*x47 + A[4]*x5 + A[5]*x12 - A[5]*x29 - A[5]*x42 - A[5]*x51 + x10*x11 - x10*x44 + x13*x14 + x16*x17 - x16*x52 + x20*x21 + x22*x23 - x22*x48 + x24*x25 - x35*x36 - x39*x40 - x39*x41 + x4*x7 - x46*x8 - x49*x50 - x56*x57 - x56*x58;
    res[2]=(*this)[11]*x27 + (*this)[12]*x33 + (*this)[13]*x3 - (*this)[13]*x59 + (*this)[15]*x37 + (*this)[1]*x1 + (*this)[1]*x19 - (*this)[1]*x32 - (*this)[2]*x53 + (*this)[2]*x6 + (*this)[3]*x65 - (*this)[4]*x28 - (*this)[4]*x31 - (*this)[7]*x15 - (*this)[7]*x68 + A[0]*x61 - A[0]*x66 + A[1]*x12 - A[1]*x29 + A[1]*x42 - A[1]*x51 + A[2]*x18 - A[2]*x34 + A[2]*x47 + A[2]*x5 - A[3]*x67 - A[3]*x69 - A[4]*x30 + A[4]*x45 - A[4]*x54 + A[4]*x64 + A[4]*x9 + A[5]*x2 + A[5]*x26 + A[5]*x38 - A[5]*x43 + x10*x25 + x11*x24 + x13*x57 - x14*x56 + x17*x63 - x24*x44 + x35*x41 + x36*x39 - x4*x46 + x50*x62 - x52*x63 + x7*x8;
    res[3]=(*this)[11]*x28 - (*this)[12]*x3 + (*this)[13]*x71 + (*this)[13]*x76 + (*this)[14]*x1 + (*this)[14]*x19 - (*this)[14]*x32 - (*this)[15]*x15 - (*this)[15]*x68 + (*this)[1]*x70 - (*this)[2]*x78 - (*this)[3]*x53 + (*this)[3]*x6 + (*this)[4]*x27 - (*this)[4]*x79 - (*this)[4]*x84 + (*this)[5]*x75 - (*this)[5]*x83 - (*this)[8]*x55 + A[0]*x12 + A[0]*x42 - A[1]*x61 + A[1]*x66 + A[1]*x73 + A[2]*x77 + A[3]*x30 + A[3]*x54 + A[3]*x64 + A[3]*x9 - A[4]*x67 + A[4]*x69 + A[5]*x60 + A[5]*x74 - A[5]*x80 + x10*x36 - x11*x35 - x13*x81 + x14*x63 - x16*x57 - x16*x58 + x17*x56 - x23*x62 + x24*x40 + x24*x41 - x25*x39 + x35*x44 + x4*x72 - x8*x82;
    res[4]=(*this)[11]*x78 + (*this)[12]*x70 - (*this)[13]*x1 - (*this)[13]*x19 + (*this)[13]*x32 + (*this)[14]*x71 + (*this)[14]*x76 - (*this)[15]*x83 + (*this)[1]*x3 - (*this)[1]*x59 + (*this)[2]*x28 + (*this)[2]*x31 + (*this)[3]*x79 + (*this)[3]*x84 - (*this)[4]*x53 + (*this)[4]*x6 + (*this)[5]*x15 + (*this)[5]*x68 + (*this)[7]*x85 - (*this)[7]*x89 - A[0]*x18 + A[0]*x34 + A[0]*x47 + A[0]*x5 + A[1]*x77 - A[2]*x66 + A[2]*x73 - A[3]*x2 - A[3]*x38 + A[3]*x43 + A[4]*x60 - A[4]*x74 + A[4]*x80 - A[5]*x67 - A[5]*x69 + x10*x81 - x11*x63 + x13*x36 - x14*x35 + x16*x25 + x16*x87 + x17*x24 - x23*x86 - x24*x52 - x39*x57 - x39*x58 + x48*x86 - x63*x88;
    res[5]=-(*this)[10]*x70 + (*this)[11]*x106 - (*this)[11]*x119 - (*this)[15]*x27 + (*this)[15]*x79 + (*this)[15]*x84 + (*this)[1]*x113 - (*this)[1]*x121 + (*this)[2]*x115 + (*this)[2]*x85 + (*this)[3]*x116 - (*this)[3]*x75 - (*this)[4]*x15 + (*this)[5]*x117 + (*this)[5]*x6 + (*this)[6]*x78 - (*this)[7]*x28 - (*this)[7]*x31 + (*this)[8]*x76 + (*this)[9]*x3 - (*this)[9]*x59 - A[0]*x100 - A[0]*x101 - A[0]*x102 - A[0]*x103 - A[0]*x104 - A[0]*x105 + A[0]*x90 + A[0]*x91 + A[0]*x92 + A[0]*x93 + A[0]*x94 + A[0]*x95 + A[0]*x96 + A[0]*x97 - A[0]*x98 - A[0]*x99 + A[1]*x108 + A[1]*x109 - A[2]*x120 + A[3]*x107 - A[4]*x124 - A[5]*x123 + x10*x114 + x110*x13 + x111*x39 + x112*x24 + x118*x36 - x122*x41 - x125*x56 - x126*x44 - x14*x8 + x20*x25 + x20*x87 - x21*x24 + x4*x58;
    res[6]=(*this)[10]*x1 + (*this)[11]*x15 + (*this)[14]*x113 - (*this)[15]*x28 - (*this)[15]*x31 - (*this)[2]*x116 + (*this)[2]*x75 + (*this)[3]*x115 + (*this)[3]*x85 - (*this)[3]*x89 + (*this)[4]*x106 - (*this)[4]*x119 + (*this)[5]*A[3]*x20 + (*this)[6]*x6 + (*this)[7]*x27 - (*this)[7]*x79 - (*this)[7]*x84 - (*this)[8]*x3 + (*this)[9]*x71 + (*this)[9]*x76 - A[0]*x108 + A[0]*x109 + A[1]*x100 - A[1]*x101 - A[1]*x102 + A[1]*x103 + A[1]*x104 + A[1]*x105 + A[1]*x90 - A[1]*x91 + A[1]*x92 - A[1]*x93 + A[1]*x94 - A[1]*x95 + A[1]*x96 - A[1]*x97 - A[1]*x98 - A[1]*x99 + A[3]*x124 + A[4]*x107 - x0*x57 + x10*x111 + x11*x122 - x110*x16 - x112*x35 - x114*x39 - x118*x25 - x122*x44 + x122*x88 + x127*x20 + x13*x72 + x17*x8 + x20*x36 + x21*x35 - x52*x8 - x63*x7;
    res[7]=-(*this)[11]*x75 + (*this)[12]*x129 + (*this)[14]*x131 - (*this)[15]*x65 - (*this)[15]*x78 + (*this)[2]*x15 - (*this)[3]*x106 + (*this)[4]*x115 + (*this)[4]*x132 + (*this)[4]*x85 + (*this)[5]*x31 - (*this)[6]*x27 + (*this)[6]*x79 + (*this)[7]*x117 + (*this)[7]*x6 + (*this)[8]*x70 - (*this)[9]*x1 - (*this)[9]*x19 + (*this)[9]*x32 + A[1]*x130 - A[2]*x100 + A[2]*x101 - A[2]*x102 + A[2]*x103 + A[2]*x104 - A[2]*x105 + A[2]*x90 - A[2]*x91 + A[2]*x92 + A[2]*x93 - A[2]*x94 + A[2]*x95 - A[2]*x96 - A[2]*x97 + A[2]*x98 - A[2]*x99 + A[3]*x123 - A[3]*x134 + A[4]*x128 + A[5]*x107 + x0*x25 - x10*x72 - x11*x136 + x111*x13 + x114*x16 - x118*x57 + x122*x14 - x125*x35 + x127*x4 - x133*x39 - x135*x24 + x20*x81 - x24*x82 + x35*x7 - x40*x8 - x41*x8;
    res[8]=-(*this)[10]*x28 + (*this)[11]*x113 - (*this)[11]*x121 + (*this)[12]*x89 + (*this)[15]*x1 + (*this)[15]*x19 - (*this)[15]*x32 - (*this)[1]*x106 + (*this)[1]*x119 + (*this)[3]*x55 - (*this)[4]*x129 - (*this)[6]*x3 + (*this)[8]*A[4]*x20 - (*this)[8]*x53 + (*this)[8]*x6 + (*this)[9]*A[0]*x8 + (*this)[9]*x78 + A[0]*x107 - A[1]*x124 + A[2]*x123 + A[2]*x134 + A[3]*x100 + A[3]*x101 - A[3]*x102 + A[3]*x103 - A[3]*x104 - A[3]*x105 + A[3]*x90 - A[3]*x91 - A[3]*x92 - A[3]*x93 - A[3]*x94 + A[3]*x95 + A[3]*x96 + A[3]*x97 - A[3]*x98 + A[3]*x99 - A[4]*x108 - A[4]*x109 + A[5]*x120 + x0*x41 + x10*x125 - x10*x7 + x11*x4 - x112*x13 - x118*x52 - x126*x58 + x13*x21 + x133*x24 + x135*x39 - x136*x36 - x137*x49 - x14*x20 - x25*x8 + x39*x82 - x4*x44;
    res[9]=(*this)[10]*x27 + (*this)[11]*x129 + (*this)[11]*x138 - (*this)[12]*x75 - (*this)[14]*x106 + (*this)[15]*x70 - (*this)[1]*x15 + (*this)[2]*x140 + (*this)[3]*A[5]*x24 + (*this)[3]*x131 + (*this)[3]*x139 + (*this)[4]*x113 - (*this)[4]*x121 + (*this)[5]*x3 - (*this)[5]*x59 - (*this)[6]*x71 - (*this)[7]*A[2]*x20 - (*this)[7]*x1 - (*this)[7]*x19 + (*this)[7]*x32 - (*this)[8]*x78 - (*this)[9]*x53 + (*this)[9]*x6 + A[0]*x124 + A[1]*x107 - A[4]*x100 + A[4]*x101 - A[4]*x102 - A[4]*x103 + A[4]*x104 + A[4]*x105 + A[4]*x90 + A[4]*x91 - A[4]*x92 + A[4]*x93 - A[4]*x94 - A[4]*x95 + A[4]*x96 - A[4]*x97 - A[4]*x98 + A[4]*x99 + x10*x82 - x110*x35 + x112*x16 + x118*x14 - x122*x57 - x122*x58 - x125*x39 - x127*x8 - x133*x35 + x137*x22 - x16*x21 + x17*x20 - x20*x52 + x39*x7 + x4*x40;
    res[10]=-(*this)[11]*A[0]*x10 - (*this)[11]*x140 + (*this)[12]*A[0]*x24 - (*this)[12]*x15 - (*this)[14]*x132 + (*this)[15]*(*this)[8]*x86 + (*this)[15]*x3 - (*this)[15]*x59 + (*this)[1]*x75 + (*this)[2]*x138 - (*this)[3]*x113 + (*this)[3]*x121 + (*this)[4]*x131 + (*this)[4]*x139 + (*this)[6]*x1 - (*this)[6]*x32 - (*this)[7]*x33 - (*this)[7]*x71 - (*this)[7]*x76 + (*this)[8]*x28 - (*this)[9]*x27 + (*this)[9]*x79 + (*this)[9]*x84 + A[1]*x128 + A[2]*x107 - A[4]*x130 + A[5]*x100 - A[5]*x101 - A[5]*x102 - A[5]*x103 + A[5]*x104 - A[5]*x105 + A[5]*x90 + A[5]*x91 - A[5]*x92 - A[5]*x93 + A[5]*x94 + A[5]*x95 - A[5]*x96 - A[5]*x97 + A[5]*x98 + A[5]*x99 - x0*x14 - x11*x118 + x111*x24 + x114*x35 + x118*x44 + x122*x25 + x125*x16 + x13*x82 - x16*x7 + x17*x4 - x20*x40 - x20*x41 - x21*x39 - x4*x52;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//--------------------------------------
// (R130MV, R130B2Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8] + B[0];
    res[9]=A[9] + B[1];
    res[10]=A[10] + B[2];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8] - B[0];
    res[9]=A[9] - B[1];
    res[10]=A[10] - B[2];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[10]*B[2] - A[8]*B[0] - A[9]*B[1];
    res[1]=-A[12]*B[0] - A[13]*B[1] - A[14]*B[2];
    res[2]=A[11]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[3]=A[11]*B[1] + A[2]*B[2] - A[4]*B[0];
    res[4]=A[11]*B[2] - A[2]*B[1] + A[3]*B[0];
    res[5]=-A[15]*B[0] - A[6]*B[2] + A[7]*B[1];
    res[6]=-A[15]*B[1] + A[5]*B[2] - A[7]*B[0];
    res[7]=-A[15]*B[2] - A[5]*B[1] + A[6]*B[0];
    res[8]=A[0]*B[0] + A[10]*B[1] - A[9]*B[2];
    res[9]=A[0]*B[1] - A[10]*B[0] + A[8]*B[2];
    res[10]=A[0]*B[2] - A[8]*B[1] + A[9]*B[0];
    res[11]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[12]=-A[13]*B[2] + A[14]*B[1] + A[1]*B[0];
    res[13]=A[12]*B[2] - A[14]*B[0] + A[1]*B[1];
    res[14]=-A[12]*B[1] + A[13]*B[0] + A[1]*B[2];
    res[15]=A[5]*B[0] + A[6]*B[1] + A[7]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[10]*B[2] - A[8]*B[0] - A[9]*B[1];
    res[1]=-A[12]*B[0] - A[13]*B[1] - A[14]*B[2];
    res[2]=A[11]*B[0];
    res[3]=A[11]*B[1];
    res[4]=A[11]*B[2];
    res[5]=-A[15]*B[0];
    res[6]=-A[15]*B[1];
    res[7]=-A[15]*B[2];
    res[8]=A[0]*B[0];
    res[9]=A[0]*B[1];
    res[10]=A[0]*B[2];
    res[11]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[12]=A[1]*B[0];
    res[13]=A[1]*B[1];
    res[14]=A[1]*B[2];
    res[15]=A[5]*B[0] + A[6]*B[1] + A[7]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B2Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[3]*B[2] + A[4]*B[1];
    res[3]=A[2]*B[2] - A[4]*B[0];
    res[4]=-A[2]*B[1] + A[3]*B[0];
    res[5]=-A[6]*B[2] + A[7]*B[1];
    res[6]=A[5]*B[2] - A[7]*B[0];
    res[7]=-A[5]*B[1] + A[6]*B[0];
    res[8]=A[10]*B[1] - A[9]*B[2];
    res[9]=-A[10]*B[0] + A[8]*B[2];
    res[10]=-A[8]*B[1] + A[9]*B[0];
    res[11]=0;
    res[12]=-A[13]*B[2] + A[14]*B[1];
    res[13]=A[12]*B[2] - A[14]*B[0];
    res[14]=-A[12]*B[1] + A[13]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*(*this)[10];
    T x1 = (*this)[12]*A[1];
    T x2 = A[2]*x0;
    T x3 = 2.0*(*this)[8];
    T x4 = A[2]*x3;
    T x5 = 2.0*(*this)[9];
    T x6 = (*this)[14]*A[0];
    T x7 = A[0]*x3;
    T x8 = A[1]*x5;
    T x9 = 2.0*A[2];
    T x10 = (*this)[6]*x9;
    T x11 = 2.0*A[0];
    T x12 = (*this)[7]*x11;
    T x13 = 2.0*A[1];
    T x14 = (*this)[5]*x13;
    T x15 = (*this)[0]*x11;
    T x16 = (*this)[0]*x13;
    T x17 = (*this)[0]*x9;
    T x18 = A[0]*x0;
    T x19 = (*this)[5]*x11;
    T x20 = (*this)[6]*x13;
    T x21 = (*this)[7]*x9;
    T x22 = A[2]*x5;
    T x23 = A[1]*x3;
    T x24 = (*this)[15]*x11;
    T x25 = (*this)[15]*x13;
    T x26 = (*this)[15]*x9;
    T x27 = (*this)[7]*x13;
    T x28 = (*this)[5]*x9;
    T x29 = (*this)[6]*x11;
    T x30 = 2.0*x1;
    T x31 = A[1]*x0;
    T x32 = 2.0*x6;
    T x33 = A[0]*x5;
    T x34 = (*this)[13]*x9;
    T x35 = (*this)[11]*(*this)[1];
    T x36 = (*this)[12]*x11;
    T x37 = (*this)[12]*x9;
    T x38 = (*this)[13]*(*this)[2];
    T x39 = (*this)[14]*x9;
    T x40 = (*this)[1]*x9;
    T x41 = (*this)[14]*x13;
    T x42 = (*this)[13]*(*this)[3];
    T x43 = (*this)[1]*(*this)[4];
    T x44 = (*this)[4]*x13;
    T x45 = (*this)[2]*x13;
    T x46 = (*this)[11]*x11;
    T x47 = (*this)[3]*x11;
    T x48 = std::pow((*this)[0], 2);
    T x49 = std::pow((*this)[11], 2);
    T x50 = std::pow((*this)[13], 2);
    T x51 = std::pow((*this)[14], 2);
    T x52 = std::pow((*this)[2], 2);
    T x53 = std::pow((*this)[6], 2);
    T x54 = std::pow((*this)[7], 2);
    T x55 = std::pow((*this)[8], 2);
    T x56 = std::pow((*this)[10], 2);
    T x57 = std::pow((*this)[12], 2);
    T x58 = std::pow((*this)[15], 2);
    T x59 = std::pow((*this)[1], 2);
    T x60 = std::pow((*this)[3], 2);
    T x61 = std::pow((*this)[4], 2);
    T x62 = std::pow((*this)[5], 2);
    T x63 = std::pow((*this)[9], 2);
    T x64 = (*this)[3]*x9;
    T x65 = (*this)[2]*x9;
    res[0]=0;
    res[1]=-(*this)[11]*x19 - (*this)[11]*x20 - (*this)[11]*x21 - (*this)[12]*x15 - (*this)[12]*x22 - (*this)[13]*x16 - (*this)[13]*x18 + (*this)[13]*x4 - (*this)[14]*x17 - (*this)[14]*x23 + (*this)[1]*x2 + (*this)[1]*x7 + (*this)[1]*x8 + (*this)[2]*x10 - (*this)[2]*x24 - (*this)[2]*x27 + (*this)[3]*x12 - (*this)[3]*x25 - (*this)[3]*x28 + (*this)[4]*x14 - (*this)[4]*x26 - (*this)[4]*x29 + x0*x1 + x5*x6;
    res[2]=(*this)[11]*x15 + (*this)[11]*x22 - (*this)[11]*x31 + (*this)[12]*x19 + (*this)[12]*x21 + (*this)[13]*x14 - (*this)[13]*x26 - (*this)[13]*x29 + (*this)[14]*x25 + (*this)[14]*x28 + (*this)[1]*x10 - (*this)[1]*x24 - (*this)[1]*x27 + (*this)[2]*x2 + (*this)[2]*x7 + (*this)[2]*x8 + (*this)[3]*x17 + (*this)[3]*x23 - (*this)[3]*x33 - (*this)[4]*x16 - (*this)[4]*x18 + (*this)[4]*x4 + (*this)[6]*x30 - (*this)[7]*x32;
    res[3]=(*this)[11]*x16 + (*this)[11]*x18 - (*this)[11]*x4 + (*this)[12]*x26 + (*this)[12]*x29 + (*this)[13]*x19 + (*this)[13]*x20 + (*this)[13]*x21 + (*this)[14]*x10 - (*this)[14]*x27 - (*this)[15]*x32 + (*this)[1]*x12 - (*this)[1]*x25 - (*this)[1]*x28 - (*this)[2]*x17 - (*this)[2]*x23 + (*this)[2]*x33 + (*this)[3]*x2 + (*this)[3]*x7 + (*this)[3]*x8 + (*this)[4]*x15 + (*this)[4]*x22 - (*this)[4]*x31 - (*this)[5]*x30;
    res[4]=(*this)[11]*x17 + (*this)[11]*x23 - (*this)[11]*x33 + (*this)[12]*x12 - (*this)[12]*x28 - (*this)[13]*x10 + (*this)[13]*x24 + (*this)[13]*x27 + (*this)[14]*x20 + (*this)[14]*x21 - (*this)[15]*x30 + (*this)[1]*x14 - (*this)[1]*x26 - (*this)[1]*x29 + (*this)[2]*x16 + (*this)[2]*x18 - (*this)[2]*x4 - (*this)[3]*x15 - (*this)[3]*x22 + (*this)[3]*x31 + (*this)[4]*x2 + (*this)[4]*x7 + (*this)[4]*x8 + (*this)[5]*x32;
    res[5]=(*this)[0]*x10 + (*this)[11]*x34 - (*this)[11]*x41 - (*this)[15]*x15 - (*this)[15]*x22 + (*this)[15]*x31 + (*this)[2]*x36 + (*this)[2]*x39 + (*this)[3]*x30 + (*this)[3]*x40 - (*this)[4]*x32 + (*this)[4]*x37 + (*this)[5]*x2 + (*this)[5]*x7 + (*this)[5]*x8 + (*this)[6]*x23 - (*this)[6]*x33 - (*this)[7]*x16 - (*this)[7]*x18 + (*this)[7]*x4 + x11*x35 - x11*x42 + x13*x38 - x13*x43;
    res[6]=(*this)[0]*x12 + (*this)[11]*x32 - (*this)[11]*x37 - (*this)[15]*x16 - (*this)[15]*x18 + (*this)[15]*x4 - (*this)[2]*x30 - (*this)[2]*x40 + (*this)[3]*x36 + (*this)[3]*x39 + (*this)[4]*x34 - (*this)[4]*x41 - (*this)[5]*x17 - (*this)[5]*x23 + (*this)[5]*x33 + (*this)[6]*x2 + (*this)[6]*x7 + (*this)[6]*x8 + (*this)[7]*x22 - (*this)[7]*x31 + x11*x38 + x11*x43 + x13*x35 + x13*x42;
    res[7]=(*this)[0]*x14 + (*this)[11]*x30 + (*this)[13]*x44 - (*this)[13]*x46 - (*this)[15]*x17 - (*this)[15]*x23 + (*this)[15]*x33 + (*this)[1]*x45 - (*this)[1]*x47 + (*this)[2]*x32 - (*this)[2]*x37 - (*this)[3]*x34 + (*this)[3]*x41 + (*this)[4]*x36 + (*this)[4]*x39 + (*this)[5]*x18 - (*this)[5]*x4 - (*this)[6]*x15 - (*this)[6]*x22 + (*this)[6]*x31 + (*this)[7]*x2 + (*this)[7]*x7 + (*this)[7]*x8 + x35*x9;
    res[8]=(*this)[0]*x22 - (*this)[0]*x31 - (*this)[11]*x44 + (*this)[11]*x64 - (*this)[13]*x30 - (*this)[14]*x37 + (*this)[15]*x10 - (*this)[1]*x34 + (*this)[1]*x41 + (*this)[3]*x45 + (*this)[4]*x65 - (*this)[5]*x21 - (*this)[6]*x14 - (*this)[7]*x25 + (*this)[8]*x2 + (*this)[9]*x23 + A[0]*x48 + A[0]*x49 + A[0]*x50 + A[0]*x51 + A[0]*x52 + A[0]*x53 + A[0]*x54 + A[0]*x55 - A[0]*x56 - A[0]*x57 - A[0]*x58 - A[0]*x59 - A[0]*x60 - A[0]*x61 - A[0]*x62 - A[0]*x63;
    res[9]=(*this)[0]*x18 - (*this)[0]*x4 - (*this)[11]*x65 - (*this)[13]*x36 - (*this)[14]*x34 + (*this)[15]*x12 - (*this)[1]*x32 + (*this)[1]*x37 + (*this)[2]*x47 + (*this)[4]*x46 + (*this)[4]*x64 - (*this)[5]*x26 - (*this)[6]*x19 - (*this)[7]*x10 + (*this)[9]*x2 + (*this)[9]*x7 + A[1]*x48 + A[1]*x49 - A[1]*x50 + A[1]*x51 - A[1]*x52 - A[1]*x53 + A[1]*x54 - A[1]*x55 - A[1]*x56 + A[1]*x57 - A[1]*x58 - A[1]*x59 + A[1]*x60 - A[1]*x61 + A[1]*x62 + A[1]*x63;
    res[10]=(*this)[0]*x23 - (*this)[0]*x33 + (*this)[11]*x45 - (*this)[12]*x32 + (*this)[13]*(*this)[1]*x11 - (*this)[13]*x41 + (*this)[15]*x14 - (*this)[1]*x30 + (*this)[2]*(*this)[4]*x11 + (*this)[3]*x44 - (*this)[3]*x46 - (*this)[5]*x12 - (*this)[6]*x24 - (*this)[7]*x20 + (*this)[8]*x18 + (*this)[9]*x31 + A[2]*x48 + A[2]*x49 + A[2]*x50 - A[2]*x51 - A[2]*x52 + A[2]*x53 - A[2]*x54 - A[2]*x55 + A[2]*x56 + A[2]*x57 - A[2]*x58 - A[2]*x59 - A[2]*x60 + A[2]*x61 + A[2]*x62 - A[2]*x63;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//--------------------------------------
// (R130MV, R130B2Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] + B[0];
    res[6]=A[6] + B[1];
    res[7]=A[7] + B[2];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] - B[0];
    res[6]=A[6] - B[1];
    res[7]=A[7] - B[2];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[5]*B[0] + A[6]*B[1] + A[7]*B[2];
    res[1]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[2]=A[13]*B[2] - A[14]*B[1] - A[1]*B[0];
    res[3]=-A[12]*B[2] + A[14]*B[0] - A[1]*B[1];
    res[4]=A[12]*B[1] - A[13]*B[0] - A[1]*B[2];
    res[5]=A[0]*B[0] + A[10]*B[1] - A[9]*B[2];
    res[6]=A[0]*B[1] - A[10]*B[0] + A[8]*B[2];
    res[7]=A[0]*B[2] - A[8]*B[1] + A[9]*B[0];
    res[8]=A[15]*B[0] + A[6]*B[2] - A[7]*B[1];
    res[9]=A[15]*B[1] - A[5]*B[2] + A[7]*B[0];
    res[10]=A[15]*B[2] + A[5]*B[1] - A[6]*B[0];
    res[11]=A[12]*B[0] + A[13]*B[1] + A[14]*B[2];
    res[12]=A[11]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[13]=A[11]*B[1] + A[2]*B[2] - A[4]*B[0];
    res[14]=A[11]*B[2] - A[2]*B[1] + A[3]*B[0];
    res[15]=A[10]*B[2] + A[8]*B[0] + A[9]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[5]*B[0] + A[6]*B[1] + A[7]*B[2];
    res[1]=0;
    res[2]=A[13]*B[2] - A[14]*B[1];
    res[3]=-A[12]*B[2] + A[14]*B[0];
    res[4]=A[12]*B[1] - A[13]*B[0];
    res[5]=A[0]*B[0];
    res[6]=A[0]*B[1];
    res[7]=A[0]*B[2];
    res[8]=A[15]*B[0];
    res[9]=A[15]*B[1];
    res[10]=A[15]*B[2];
    res[11]=0;
    res[12]=-A[3]*B[2] + A[4]*B[1];
    res[13]=A[2]*B[2] - A[4]*B[0];
    res[14]=-A[2]*B[1] + A[3]*B[0];
    res[15]=A[10]*B[2] + A[8]*B[0] + A[9]*B[1];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[2]*B[0] - A[3]*B[1] - A[4]*B[2];
    res[2]=-A[1]*B[0];
    res[3]=-A[1]*B[1];
    res[4]=-A[1]*B[2];
    res[5]=A[10]*B[1] - A[9]*B[2];
    res[6]=-A[10]*B[0] + A[8]*B[2];
    res[7]=-A[8]*B[1] + A[9]*B[0];
    res[8]=A[6]*B[2] - A[7]*B[1];
    res[9]=-A[5]*B[2] + A[7]*B[0];
    res[10]=A[5]*B[1] - A[6]*B[0];
    res[11]=A[12]*B[0] + A[13]*B[1] + A[14]*B[2];
    res[12]=A[11]*B[0];
    res[13]=A[11]*B[1];
    res[14]=A[11]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[0]*x0;
    T x2 = A[1]*x0;
    T x3 = A[2]*x0;
    T x4 = 2.0*A[0];
    T x5 = (*this)[3]*x4;
    T x6 = 2.0*A[2];
    T x7 = (*this)[12]*x6;
    T x8 = (*this)[7]*x4;
    T x9 = 2.0*A[1];
    T x10 = (*this)[5]*x9;
    T x11 = (*this)[9]*x6;
    T x12 = (*this)[8]*x9;
    T x13 = (*this)[11]*x6;
    T x14 = (*this)[10]*x9;
    T x15 = (*this)[11]*x4;
    T x16 = (*this)[11]*x9;
    T x17 = (*this)[12]*x4;
    T x18 = (*this)[12]*x9;
    T x19 = (*this)[15]*x9;
    T x20 = (*this)[5]*x6;
    T x21 = (*this)[14]*x6;
    T x22 = (*this)[14]*x4;
    T x23 = (*this)[5]*x4;
    T x24 = (*this)[6]*x9;
    T x25 = (*this)[7]*x6;
    T x26 = (*this)[3]*x6;
    T x27 = (*this)[9]*x4;
    T x28 = (*this)[10]*x6;
    T x29 = (*this)[10]*x4;
    T x30 = (*this)[13]*x4;
    T x31 = (*this)[9]*x9;
    T x32 = (*this)[15]*x4;
    T x33 = (*this)[7]*x9;
    T x34 = (*this)[15]*x6;
    T x35 = (*this)[1]*x6;
    T x36 = (*this)[2]*x4;
    T x37 = (*this)[4]*x6;
    T x38 = std::pow((*this)[0], 2);
    T x39 = std::pow((*this)[12], 2);
    T x40 = std::pow((*this)[1], 2);
    T x41 = std::pow((*this)[3], 2);
    T x42 = std::pow((*this)[4], 2);
    T x43 = std::pow((*this)[6], 2);
    T x44 = std::pow((*this)[7], 2);
    T x45 = std::pow((*this)[8], 2);
    T x46 = std::pow((*this)[10], 2);
    T x47 = std::pow((*this)[11], 2);
    T x48 = std::pow((*this)[13], 2);
    T x49 = std::pow((*this)[14], 2);
    T x50 = std::pow((*this)[15], 2);
    T x51 = std::pow((*this)[2], 2);
    T x52 = std::pow((*this)[5], 2);
    T x53 = std::pow((*this)[9], 2);
    T x54 = (*this)[14]*x9;
    T x55 = (*this)[3]*x9;
    T x56 = (*this)[13]*x9;
    T x57 = (*this)[1]*(*this)[4];
    res[0]=0;
    res[1]=-(*this)[10]*x13 + (*this)[10]*x5 - (*this)[13]*x19 - (*this)[13]*x20 + (*this)[13]*x8 + (*this)[14]*x10 - (*this)[15]*x17 - (*this)[15]*x21 - (*this)[1]*x23 - (*this)[1]*x24 - (*this)[1]*x25 + (*this)[2]*x1 + (*this)[2]*x11 - (*this)[2]*x14 + (*this)[3]*x2 + (*this)[4]*x12 - (*this)[4]*x27 + (*this)[4]*x3 - (*this)[6]*x22 + (*this)[6]*x7 - (*this)[7]*x18 - (*this)[8]*x15 - (*this)[8]*x26 - (*this)[9]*x16;
    res[2]=-(*this)[10]*x22 + (*this)[10]*x7 + (*this)[13]*x12 - (*this)[13]*x27 + (*this)[13]*x3 - (*this)[14]*x2 + (*this)[15]*x15 + (*this)[15]*x26 + (*this)[1]*x1 + (*this)[1]*x11 - (*this)[1]*x14 - (*this)[2]*x23 - (*this)[2]*x24 - (*this)[2]*x25 - (*this)[3]*x10 - (*this)[4]*x19 - (*this)[4]*x20 + (*this)[4]*x8 - (*this)[6]*x13 + (*this)[6]*x5 + (*this)[7]*x16 + (*this)[8]*x17 + (*this)[8]*x21 + (*this)[9]*x18;
    res[3]=-(*this)[11]*x8 - (*this)[12]*x12 - (*this)[12]*x3 + (*this)[13]*x28 + (*this)[13]*x31 + (*this)[14]*x1 + (*this)[14]*x11 - (*this)[14]*x14 + (*this)[15]*x16 + (*this)[1]*x2 + (*this)[1]*x29 + (*this)[2]*x10 - (*this)[2]*x34 - (*this)[3]*x24 - (*this)[3]*x25 + (*this)[4]*x32 + (*this)[4]*x33 + (*this)[5]*x13 - (*this)[5]*x5 - (*this)[6]*x36 - (*this)[6]*x37 + (*this)[8]*x30 - (*this)[8]*x35 + (*this)[9]*x17;
    res[4]=(*this)[10]*x17 + (*this)[10]*x21 - (*this)[11]*x10 + (*this)[12]*x2 - (*this)[13]*x1 - (*this)[13]*x11 + (*this)[13]*x14 + (*this)[14]*x31 + (*this)[15]*x13 - (*this)[15]*x5 + (*this)[1]*x12 - (*this)[1]*x27 + (*this)[1]*x3 + (*this)[2]*x19 + (*this)[2]*x20 - (*this)[2]*x8 - (*this)[3]*x33 - (*this)[4]*x23 - (*this)[4]*x24 - (*this)[4]*x25 + (*this)[6]*x15 + (*this)[6]*x26 + (*this)[8]*x22 - (*this)[8]*x7;
    res[5]=-(*this)[10]*x2 + (*this)[13]*x18 + (*this)[13]*x35 + (*this)[14]*x7 - (*this)[1]*x54 - (*this)[2]*x37 - (*this)[2]*x55 - (*this)[3]*x13 + (*this)[4]*x16 - (*this)[6]*x10 + (*this)[6]*x34 - (*this)[7]*x19 - (*this)[7]*x20 + (*this)[8]*x28 + (*this)[9]*x12 + (*this)[9]*x3 + A[0]*x38 + A[0]*x39 + A[0]*x40 + A[0]*x41 + A[0]*x42 + A[0]*x43 + A[0]*x44 + A[0]*x45 - A[0]*x46 - A[0]*x47 - A[0]*x48 - A[0]*x49 - A[0]*x50 - A[0]*x51 - A[0]*x52 - A[0]*x53;
    res[6]=(*this)[10]*x1 + (*this)[10]*x11 + (*this)[13]*x17 + (*this)[13]*x21 - (*this)[15]*x20 + (*this)[15]*x8 + (*this)[1]*x22 - (*this)[1]*x7 + (*this)[2]*x13 - (*this)[2]*x5 - (*this)[4]*x15 - (*this)[4]*x26 - (*this)[6]*x23 - (*this)[6]*x25 + (*this)[8]*x27 - (*this)[8]*x3 + A[1]*x38 - A[1]*x39 + A[1]*x40 - A[1]*x41 + A[1]*x42 - A[1]*x43 + A[1]*x44 - A[1]*x45 - A[1]*x46 - A[1]*x47 + A[1]*x48 - A[1]*x49 - A[1]*x50 + A[1]*x51 + A[1]*x52 + A[1]*x53;
    res[7]=(*this)[11]*x5 + (*this)[13]*x54 + (*this)[14]*x17 + (*this)[15]*x10 + (*this)[1]*x18 - (*this)[1]*x30 - (*this)[2]*x16 - (*this)[4]*x36 - (*this)[4]*x55 - (*this)[5]*x8 - (*this)[6]*x32 - (*this)[7]*x24 + (*this)[8]*x2 + (*this)[8]*x29 - (*this)[9]*x1 + (*this)[9]*x14 + A[2]*x38 - A[2]*x39 + A[2]*x40 + A[2]*x41 - A[2]*x42 + A[2]*x43 - A[2]*x44 - A[2]*x45 + A[2]*x46 - A[2]*x47 - A[2]*x48 + A[2]*x49 - A[2]*x50 + A[2]*x51 + A[2]*x52 - A[2]*x53;
    res[8]=-(*this)[10]*x20 + (*this)[10]*x8 + (*this)[13]*x13 - (*this)[13]*x5 - (*this)[14]*x16 + (*this)[15]*x1 + (*this)[15]*x11 - (*this)[15]*x14 + (*this)[1]*x15 + (*this)[1]*x26 + (*this)[2]*x17 + (*this)[2]*x21 + (*this)[2]*x56 + (*this)[3]*x18 - (*this)[4]*x22 + (*this)[4]*x7 - (*this)[6]*x12 + (*this)[6]*x27 - (*this)[6]*x3 + (*this)[7]*x2 - (*this)[8]*x23 - (*this)[8]*x25 - (*this)[9]*x10 - x57*x9;
    res[9]=-(*this)[11]*x7 + (*this)[12]*x5 + (*this)[13]*x37 + (*this)[13]*x55 + (*this)[14]*x15 + (*this)[15]*x2 + (*this)[15]*x29 + (*this)[1]*x16 - (*this)[2]*x18 + (*this)[2]*x30 - (*this)[2]*x35 + (*this)[3]*x21 - (*this)[4]*x54 + (*this)[5]*x3 - (*this)[6]*(*this)[8]*x4 - (*this)[6]*x28 - (*this)[7]*x1 - (*this)[7]*x11 + (*this)[7]*x14 + (*this)[8]*x10 - (*this)[8]*x34 - (*this)[9]*x23 - (*this)[9]*x24 + x4*x57;
    res[10]=-(*this)[10]*x23 - (*this)[10]*x25 + (*this)[12]*x16 - (*this)[13]*x15 - (*this)[13]*x26 + (*this)[15]*x12 - (*this)[15]*x27 + (*this)[15]*x3 + (*this)[1]*(*this)[2]*x9 + (*this)[1]*x13 - (*this)[1]*x5 + (*this)[2]*x22 - (*this)[2]*x7 + (*this)[3]*x54 + (*this)[4]*x17 + (*this)[4]*x21 + (*this)[4]*x56 - (*this)[5]*x2 + (*this)[6]*x1 + (*this)[6]*x11 - (*this)[6]*x14 - (*this)[7]*x31 + (*this)[8]*x20 - (*this)[8]*x8;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

//-----------------------------------
// (R130MV, R130B3) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11] + B[0];
    res[12]=A[12] + B[1];
    res[13]=A[13] + B[2];
    res[14]=A[14] + B[3];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11] - B[0];
    res[12]=A[12] - B[1];
    res[13]=A[13] - B[2];
    res[14]=A[14] - B[3];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[11]*B[0] - A[12]*B[1] - A[13]*B[2] - A[14]*B[3];
    res[1]=-A[10]*B[3] + A[15]*B[0] - A[8]*B[1] - A[9]*B[2];
    res[2]=-A[15]*B[1] - A[6]*B[3] + A[7]*B[2] + A[8]*B[0];
    res[3]=-A[15]*B[2] + A[5]*B[3] - A[7]*B[1] + A[9]*B[0];
    res[4]=A[10]*B[0] - A[15]*B[3] - A[5]*B[2] + A[6]*B[1];
    res[5]=A[11]*B[1] - A[12]*B[0] - A[3]*B[3] + A[4]*B[2];
    res[6]=A[11]*B[2] - A[13]*B[0] + A[2]*B[3] - A[4]*B[1];
    res[7]=A[11]*B[3] - A[14]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[8]=-A[13]*B[3] + A[14]*B[2] + A[1]*B[1] + A[2]*B[0];
    res[9]=A[12]*B[3] - A[14]*B[1] + A[1]*B[2] + A[3]*B[0];
    res[10]=-A[12]*B[2] + A[13]*B[1] + A[1]*B[3] + A[4]*B[0];
    res[11]=A[0]*B[0] - A[5]*B[1] - A[6]*B[2] - A[7]*B[3];
    res[12]=A[0]*B[1] + A[10]*B[2] - A[5]*B[0] - A[9]*B[3];
    res[13]=A[0]*B[2] - A[10]*B[1] - A[6]*B[0] + A[8]*B[3];
    res[14]=A[0]*B[3] - A[7]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2] + A[4]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[11]*B[0] - A[12]*B[1] - A[13]*B[2] - A[14]*B[3];
    res[1]=-A[10]*B[3] - A[8]*B[1] - A[9]*B[2];
    res[2]=-A[6]*B[3] + A[7]*B[2] + A[8]*B[0];
    res[3]=A[5]*B[3] - A[7]*B[1] + A[9]*B[0];
    res[4]=A[10]*B[0] - A[5]*B[2] + A[6]*B[1];
    res[5]=-A[3]*B[3] + A[4]*B[2];
    res[6]=A[2]*B[3] - A[4]*B[1];
    res[7]=-A[2]*B[2] + A[3]*B[1];
    res[8]=A[1]*B[1] + A[2]*B[0];
    res[9]=A[1]*B[2] + A[3]*B[0];
    res[10]=A[1]*B[3] + A[4]*B[0];
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[15]*B[0];
    res[2]=-A[15]*B[1];
    res[3]=-A[15]*B[2];
    res[4]=-A[15]*B[3];
    res[5]=A[11]*B[1] - A[12]*B[0];
    res[6]=A[11]*B[2] - A[13]*B[0];
    res[7]=A[11]*B[3] - A[14]*B[0];
    res[8]=-A[13]*B[3] + A[14]*B[2];
    res[9]=A[12]*B[3] - A[14]*B[1];
    res[10]=-A[12]*B[2] + A[13]*B[1];
    res[11]=-A[5]*B[1] - A[6]*B[2] - A[7]*B[3];
    res[12]=A[10]*B[2] - A[5]*B[0] - A[9]*B[3];
    res[13]=-A[10]*B[1] - A[6]*B[0] + A[8]*B[3];
    res[14]=-A[7]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2] + A[4]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B3<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*A[0];
    T x1 = (*this)[0]*x0;
    T x2 = 2.0*(*this)[10];
    T x3 = (*this)[12]*A[2];
    T x4 = A[3]*x2;
    T x5 = (*this)[10]*x0;
    T x6 = 2.0*(*this)[11];
    T x7 = (*this)[5]*A[1];
    T x8 = (*this)[6]*A[2];
    T x9 = (*this)[7]*A[3];
    T x10 = 2.0*A[3];
    T x11 = (*this)[13]*x10;
    T x12 = 2.0*A[1];
    T x13 = (*this)[14]*x12;
    T x14 = (*this)[15]*x0;
    T x15 = (*this)[15]*x12;
    T x16 = 2.0*A[2];
    T x17 = (*this)[15]*x16;
    T x18 = (*this)[4]*x10;
    T x19 = (*this)[1]*x12;
    T x20 = (*this)[1]*x16;
    T x21 = (*this)[7]*x16;
    T x22 = (*this)[2]*x0;
    T x23 = (*this)[5]*x10;
    T x24 = (*this)[3]*x0;
    T x25 = (*this)[4]*x12;
    T x26 = (*this)[12]*x12;
    T x27 = (*this)[0]*x16;
    T x28 = (*this)[0]*x10;
    T x29 = A[1]*x2;
    T x30 = (*this)[12]*x0;
    T x31 = (*this)[12]*x10;
    T x32 = (*this)[13]*x0;
    T x33 = (*this)[7]*x0;
    T x34 = (*this)[14]*x16;
    T x35 = (*this)[6]*x10;
    T x36 = (*this)[7]*x12;
    T x37 = (*this)[4]*x16;
    T x38 = std::pow((*this)[0], 2);
    T x39 = std::pow((*this)[10], 2);
    T x40 = std::pow((*this)[11], 2);
    T x41 = std::pow((*this)[12], 2);
    T x42 = std::pow((*this)[13], 2);
    T x43 = std::pow((*this)[14], 2);
    T x44 = std::pow((*this)[15], 2);
    T x45 = std::pow((*this)[1], 2);
    T x46 = std::pow((*this)[2], 2);
    T x47 = std::pow((*this)[3], 2);
    T x48 = std::pow((*this)[4], 2);
    T x49 = std::pow((*this)[5], 2);
    T x50 = std::pow((*this)[6], 2);
    T x51 = std::pow((*this)[7], 2);
    T x52 = std::pow((*this)[8], 2);
    T x53 = std::pow((*this)[9], 2);
    T x54 = 2.0*(*this)[0];
    T x55 = 2.0*x3;
    T x56 = A[2]*x2;
    T x57 = A[1]*x6;
    T x58 = A[2]*x6;
    T x59 = A[3]*x6;
    T x60 = 2.0*(*this)[5];
    T x61 = (*this)[8]*(*this)[9];
    T x62 = (*this)[2]*(*this)[3];
    T x63 = 2.0*(*this)[6];
    T x64 = 2.0*(*this)[7];
    T x65 = (*this)[0]*x12;
    T x66 = 2.0*(*this)[1];
    res[0]=-(*this)[0]*x26 + (*this)[11]*x1 - (*this)[13]*x27 - (*this)[13]*x29 - (*this)[14]*x28 - (*this)[14]*x33 + (*this)[15]*x18 + (*this)[1]*x14 + (*this)[1]*x4 + (*this)[2]*x15 + (*this)[2]*x21 - (*this)[2]*x35 + (*this)[3]*x17 + (*this)[3]*x23 - (*this)[3]*x36 + (*this)[4]*x5 - (*this)[5]*x30 - (*this)[5]*x37 + (*this)[6]*x25 - (*this)[6]*x32 + (*this)[8]*x11 + (*this)[8]*x19 + (*this)[8]*x22 - (*this)[8]*x34 + (*this)[9]*x13 + (*this)[9]*x20 + (*this)[9]*x24 - (*this)[9]*x31 + x2*x3 + x6*x7 + x6*x8 + x6*x9;
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
    res[11]=-(*this)[12]*x57 - (*this)[13]*x25 - (*this)[13]*x58 - (*this)[14]*x59 + (*this)[15]*x4 + (*this)[1]*x18 + (*this)[2]*x11 + (*this)[2]*x19 - (*this)[2]*x34 + (*this)[3]*x13 + (*this)[3]*x20 - (*this)[3]*x31 + (*this)[4]*x55 - (*this)[5]*x56 + (*this)[6]*x29 + (*this)[8]*x15 + (*this)[8]*x21 - (*this)[8]*x35 + (*this)[9]*x17 + (*this)[9]*x23 - (*this)[9]*x36 + A[0]*x38 + A[0]*x39 + A[0]*x40 + A[0]*x41 + A[0]*x42 + A[0]*x43 + A[0]*x44 + A[0]*x45 + A[0]*x46 + A[0]*x47 + A[0]*x48 + A[0]*x49 + A[0]*x50 + A[0]*x51 + A[0]*x52 + A[0]*x53 + x54*x7 + x54*x8 + x54*x9;
    res[12]=-(*this)[0]*x56 + (*this)[11]*x30 - (*this)[13]*x55 + (*this)[14]*x20 + (*this)[14]*x24 - (*this)[14]*x31 - (*this)[15]*x35 - (*this)[1]*x11 - (*this)[1]*x22 - (*this)[2]*x18 - (*this)[3]*x59 - (*this)[4]*x32 + (*this)[4]*x58 + (*this)[5]*x1 - (*this)[6]*x5 + (*this)[7]*x17 + (*this)[8]*x14 + (*this)[8]*x4 + (*this)[9]*x28 + (*this)[9]*x33 + A[1]*x38 - A[1]*x39 - A[1]*x40 - A[1]*x41 + A[1]*x42 + A[1]*x43 + A[1]*x44 - A[1]*x45 - A[1]*x46 + A[1]*x47 + A[1]*x48 + A[1]*x49 - A[1]*x50 - A[1]*x51 + A[1]*x52 - A[1]*x53 + x16*x61 - x16*x62 + x60*x8 + x60*x9;
    res[13]=(*this)[0]*x29 + (*this)[11]*x32 - (*this)[13]*x26 - (*this)[14]*x11 - (*this)[14]*x22 + (*this)[15]*x23 - (*this)[1]*x13 - (*this)[1]*x24 + (*this)[1]*x31 + (*this)[2]*x59 - (*this)[3]*x18 + (*this)[4]*x30 - (*this)[4]*x57 + (*this)[5]*x5 + (*this)[6]*x1 - (*this)[7]*x15 - (*this)[8]*x28 - (*this)[8]*x33 + (*this)[9]*x14 + (*this)[9]*x4 + A[2]*x38 - A[2]*x39 - A[2]*x40 + A[2]*x41 - A[2]*x42 + A[2]*x43 + A[2]*x44 - A[2]*x45 + A[2]*x46 - A[2]*x47 + A[2]*x48 - A[2]*x49 + A[2]*x50 - A[2]*x51 - A[2]*x52 + A[2]*x53 + x12*x61 - x12*x62 + x63*x7 + x63*x9;
    res[14]=(*this)[11]*(*this)[14]*x0 - (*this)[12]*x13 - (*this)[12]*x24 + (*this)[13]*x19 + (*this)[13]*x22 - (*this)[13]*x34 + (*this)[15]*x5 - (*this)[1]*(*this)[4]*x0 - (*this)[1]*x55 - (*this)[2]*x25 - (*this)[2]*x58 - (*this)[3]*x37 + (*this)[3]*x57 - (*this)[5]*(*this)[9]*x0 - (*this)[5]*x17 + (*this)[6]*(*this)[8]*x0 + (*this)[6]*x15 + (*this)[7]*x1 + (*this)[8]*x27 + (*this)[8]*x29 + (*this)[9]*x56 - (*this)[9]*x65 + A[3]*x38 + A[3]*x39 - A[3]*x40 + A[3]*x41 + A[3]*x42 - A[3]*x43 + A[3]*x44 - A[3]*x45 + A[3]*x46 + A[3]*x47 - A[3]*x48 - A[3]*x49 - A[3]*x50 + A[3]*x51 - A[3]*x52 - A[3]*x53 + x64*x7 + x64*x8;
    res[15]=-(*this)[0]*x18 + (*this)[11]*x14 + (*this)[11]*x4 - (*this)[12]*x15 - (*this)[13]*x17 + (*this)[13]*x36 - (*this)[14]*(*this)[15]*x10 - (*this)[14]*x5 - (*this)[1]*x1 - (*this)[2]*(*this)[9]*x10 + (*this)[2]*x56 - (*this)[2]*x65 + (*this)[3]*(*this)[8]*x10 - (*this)[3]*x27 - (*this)[3]*x29 - (*this)[4]*x33 - (*this)[5]*x11 - (*this)[5]*x22 + (*this)[5]*x34 - (*this)[6]*x13 - (*this)[6]*x24 + (*this)[6]*x31 - (*this)[7]*x55 - (*this)[8]*x30 - (*this)[8]*x37 + (*this)[8]*x57 + (*this)[9]*x25 - (*this)[9]*x32 + (*this)[9]*x58 - x66*x7 - x66*x8 - x66*x9;
    return res;
};

//--------------------------------------
// (R130MV, R130B3Sm1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12] + B[0];
    res[13]=A[13] + B[1];
    res[14]=A[14] + B[2];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12] - B[0];
    res[13]=A[13] - B[1];
    res[14]=A[14] - B[2];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[12]*B[0] - A[13]*B[1] - A[14]*B[2];
    res[1]=-A[10]*B[2] - A[8]*B[0] - A[9]*B[1];
    res[2]=-A[15]*B[0] - A[6]*B[2] + A[7]*B[1];
    res[3]=-A[15]*B[1] + A[5]*B[2] - A[7]*B[0];
    res[4]=-A[15]*B[2] - A[5]*B[1] + A[6]*B[0];
    res[5]=A[11]*B[0] - A[3]*B[2] + A[4]*B[1];
    res[6]=A[11]*B[1] + A[2]*B[2] - A[4]*B[0];
    res[7]=A[11]*B[2] - A[2]*B[1] + A[3]*B[0];
    res[8]=-A[13]*B[2] + A[14]*B[1] + A[1]*B[0];
    res[9]=A[12]*B[2] - A[14]*B[0] + A[1]*B[1];
    res[10]=-A[12]*B[1] + A[13]*B[0] + A[1]*B[2];
    res[11]=-A[5]*B[0] - A[6]*B[1] - A[7]*B[2];
    res[12]=A[0]*B[0] + A[10]*B[1] - A[9]*B[2];
    res[13]=A[0]*B[1] - A[10]*B[0] + A[8]*B[2];
    res[14]=A[0]*B[2] - A[8]*B[1] + A[9]*B[0];
    res[15]=A[2]*B[0] + A[3]*B[1] + A[4]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=-A[12]*B[0] - A[13]*B[1] - A[14]*B[2];
    res[1]=-A[10]*B[2] - A[8]*B[0] - A[9]*B[1];
    res[2]=-A[6]*B[2] + A[7]*B[1];
    res[3]=A[5]*B[2] - A[7]*B[0];
    res[4]=-A[5]*B[1] + A[6]*B[0];
    res[5]=-A[3]*B[2] + A[4]*B[1];
    res[6]=A[2]*B[2] - A[4]*B[0];
    res[7]=-A[2]*B[1] + A[3]*B[0];
    res[8]=A[1]*B[0];
    res[9]=A[1]*B[1];
    res[10]=A[1]*B[2];
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[15]*B[0];
    res[3]=-A[15]*B[1];
    res[4]=-A[15]*B[2];
    res[5]=A[11]*B[0];
    res[6]=A[11]*B[1];
    res[7]=A[11]*B[2];
    res[8]=-A[13]*B[2] + A[14]*B[1];
    res[9]=A[12]*B[2] - A[14]*B[0];
    res[10]=-A[12]*B[1] + A[13]*B[0];
    res[11]=-A[5]*B[0] - A[6]*B[1] - A[7]*B[2];
    res[12]=A[10]*B[1] - A[9]*B[2];
    res[13]=-A[10]*B[0] + A[8]*B[2];
    res[14]=-A[8]*B[1] + A[9]*B[0];
    res[15]=A[2]*B[0] + A[3]*B[1] + A[4]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*(*this)[10];
    T x1 = (*this)[12]*A[1];
    T x2 = A[2]*x0;
    T x3 = 2.0*(*this)[11];
    T x4 = (*this)[5]*A[0];
    T x5 = (*this)[6]*A[1];
    T x6 = (*this)[7]*A[2];
    T x7 = 2.0*A[2];
    T x8 = (*this)[13]*x7;
    T x9 = 2.0*A[0];
    T x10 = (*this)[14]*x9;
    T x11 = (*this)[15]*x9;
    T x12 = 2.0*A[1];
    T x13 = (*this)[15]*x12;
    T x14 = (*this)[4]*x7;
    T x15 = (*this)[1]*x9;
    T x16 = (*this)[1]*x12;
    T x17 = (*this)[7]*x12;
    T x18 = (*this)[5]*x7;
    T x19 = (*this)[4]*x9;
    T x20 = (*this)[12]*x9;
    T x21 = (*this)[0]*x12;
    T x22 = (*this)[0]*x7;
    T x23 = A[0]*x0;
    T x24 = (*this)[12]*x7;
    T x25 = (*this)[14]*x12;
    T x26 = (*this)[6]*x7;
    T x27 = (*this)[7]*x9;
    T x28 = (*this)[4]*x12;
    T x29 = 2.0*(*this)[0];
    T x30 = 2.0*x1;
    T x31 = A[1]*x0;
    T x32 = A[0]*x3;
    T x33 = A[1]*x3;
    T x34 = A[2]*x3;
    T x35 = std::pow((*this)[0], 2);
    T x36 = std::pow((*this)[13], 2);
    T x37 = std::pow((*this)[14], 2);
    T x38 = std::pow((*this)[15], 2);
    T x39 = std::pow((*this)[3], 2);
    T x40 = std::pow((*this)[4], 2);
    T x41 = std::pow((*this)[5], 2);
    T x42 = std::pow((*this)[8], 2);
    T x43 = std::pow((*this)[10], 2);
    T x44 = std::pow((*this)[11], 2);
    T x45 = std::pow((*this)[12], 2);
    T x46 = std::pow((*this)[1], 2);
    T x47 = std::pow((*this)[2], 2);
    T x48 = std::pow((*this)[6], 2);
    T x49 = std::pow((*this)[7], 2);
    T x50 = std::pow((*this)[9], 2);
    T x51 = 2.0*(*this)[5];
    T x52 = (*this)[8]*(*this)[9];
    T x53 = (*this)[2]*(*this)[3];
    T x54 = 2.0*(*this)[6];
    T x55 = 2.0*(*this)[7];
    T x56 = (*this)[0]*x9;
    T x57 = 2.0*(*this)[1];
    res[0]=-(*this)[0]*x20 - (*this)[13]*x21 - (*this)[13]*x23 - (*this)[14]*x22 + (*this)[15]*x14 + (*this)[1]*x2 + (*this)[2]*x11 + (*this)[2]*x17 - (*this)[2]*x26 + (*this)[3]*x13 + (*this)[3]*x18 - (*this)[3]*x27 - (*this)[5]*x28 + (*this)[6]*x19 + (*this)[8]*x15 - (*this)[8]*x25 + (*this)[8]*x8 + (*this)[9]*x10 + (*this)[9]*x16 - (*this)[9]*x24 + x0*x1 + x3*x4 + x3*x5 + x3*x6;
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
    res[11]=-(*this)[12]*x32 - (*this)[13]*x19 - (*this)[13]*x33 - (*this)[14]*x34 + (*this)[15]*x2 + (*this)[1]*x14 + (*this)[2]*x15 - (*this)[2]*x25 + (*this)[2]*x8 + (*this)[3]*x10 + (*this)[3]*x16 - (*this)[3]*x24 + (*this)[4]*x30 - (*this)[5]*x31 + (*this)[6]*x23 + (*this)[8]*x11 + (*this)[8]*x17 - (*this)[8]*x26 + (*this)[9]*x13 + (*this)[9]*x18 - (*this)[9]*x27 + x29*x4 + x29*x5 + x29*x6;
    res[12]=-(*this)[0]*x31 - (*this)[13]*x30 + (*this)[14]*x16 - (*this)[14]*x24 - (*this)[15]*x26 - (*this)[1]*x8 - (*this)[2]*x14 - (*this)[3]*x34 + (*this)[4]*x33 + (*this)[7]*x13 + (*this)[8]*x2 + (*this)[9]*x22 + A[0]*x35 + A[0]*x36 + A[0]*x37 + A[0]*x38 + A[0]*x39 + A[0]*x40 + A[0]*x41 + A[0]*x42 - A[0]*x43 - A[0]*x44 - A[0]*x45 - A[0]*x46 - A[0]*x47 - A[0]*x48 - A[0]*x49 - A[0]*x50 + x12*x52 - x12*x53 + x5*x51 + x51*x6;
    res[13]=(*this)[0]*x23 - (*this)[13]*x20 - (*this)[14]*x8 + (*this)[15]*x18 - (*this)[1]*x10 + (*this)[1]*x24 + (*this)[2]*x34 - (*this)[3]*x14 - (*this)[4]*x32 - (*this)[7]*x11 - (*this)[8]*x22 + (*this)[9]*x2 + A[1]*x35 - A[1]*x36 + A[1]*x37 + A[1]*x38 - A[1]*x39 + A[1]*x40 - A[1]*x41 - A[1]*x42 - A[1]*x43 - A[1]*x44 + A[1]*x45 - A[1]*x46 + A[1]*x47 + A[1]*x48 - A[1]*x49 + A[1]*x50 + x4*x54 + x52*x9 - x53*x9 + x54*x6;
    res[14]=-(*this)[12]*x10 + (*this)[13]*x15 - (*this)[13]*x25 - (*this)[1]*x30 - (*this)[2]*x19 - (*this)[2]*x33 - (*this)[3]*x28 + (*this)[3]*x32 - (*this)[5]*x13 + (*this)[6]*x11 + (*this)[8]*x21 + (*this)[8]*x23 + (*this)[9]*x31 - (*this)[9]*x56 + A[2]*x35 + A[2]*x36 - A[2]*x37 + A[2]*x38 + A[2]*x39 - A[2]*x40 - A[2]*x41 - A[2]*x42 + A[2]*x43 - A[2]*x44 + A[2]*x45 - A[2]*x46 + A[2]*x47 - A[2]*x48 + A[2]*x49 - A[2]*x50 + x4*x55 + x5*x55;
    res[15]=-(*this)[0]*x14 + (*this)[11]*x2 - (*this)[12]*x11 - (*this)[13]*x13 + (*this)[13]*x27 - (*this)[14]*(*this)[15]*x7 - (*this)[2]*(*this)[9]*x7 + (*this)[2]*x31 - (*this)[2]*x56 + (*this)[3]*(*this)[8]*x7 - (*this)[3]*x21 - (*this)[3]*x23 + (*this)[5]*x25 - (*this)[5]*x8 - (*this)[6]*x10 + (*this)[6]*x24 - (*this)[7]*x30 - (*this)[8]*x28 + (*this)[8]*x32 + (*this)[9]*x19 + (*this)[9]*x33 - x4*x57 - x5*x57 - x57*x6;
    return res;
};

//--------------------------------------
// (R130MV, R130B3Sp1) binary operations
//--------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11] + B[0];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11] - B[0];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[11]*B[0];
    res[1]=A[15]*B[0];
    res[2]=A[8]*B[0];
    res[3]=A[9]*B[0];
    res[4]=A[10]*B[0];
    res[5]=-A[12]*B[0];
    res[6]=-A[13]*B[0];
    res[7]=-A[14]*B[0];
    res[8]=A[2]*B[0];
    res[9]=A[3]*B[0];
    res[10]=A[4]*B[0];
    res[11]=A[0]*B[0];
    res[12]=-A[5]*B[0];
    res[13]=-A[6]*B[0];
    res[14]=-A[7]*B[0];
    res[15]=A[1]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[11]*B[0];
    res[1]=0;
    res[2]=A[8]*B[0];
    res[3]=A[9]*B[0];
    res[4]=A[10]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[2]*B[0];
    res[9]=A[3]*B[0];
    res[10]=A[4]*B[0];
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[15]*B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[12]*B[0];
    res[6]=-A[13]*B[0];
    res[7]=-A[14]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=-A[5]*B[0];
    res[13]=-A[6]*B[0];
    res[14]=-A[7]*B[0];
    res[15]=A[1]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*A[0];
    res[0]=x0*((*this)[0]*(*this)[11] + (*this)[10]*(*this)[4] - (*this)[12]*(*this)[5] - (*this)[13]*(*this)[6] - (*this)[14]*(*this)[7] + (*this)[15]*(*this)[1] + (*this)[2]*(*this)[8] + (*this)[3]*(*this)[9]);
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
    res[11]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[10], 2) + std::pow((*this)[11], 2) + std::pow((*this)[12], 2) + std::pow((*this)[13], 2) + std::pow((*this)[14], 2) + std::pow((*this)[15], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) + std::pow((*this)[7], 2) + std::pow((*this)[8], 2) + std::pow((*this)[9], 2));
    res[12]=x0*((*this)[0]*(*this)[5] - (*this)[10]*(*this)[6] + (*this)[11]*(*this)[12] - (*this)[13]*(*this)[4] + (*this)[14]*(*this)[3] + (*this)[15]*(*this)[8] - (*this)[1]*(*this)[2] + (*this)[7]*(*this)[9]);
    res[13]=x0*((*this)[0]*(*this)[6] + (*this)[10]*(*this)[5] + (*this)[11]*(*this)[13] + (*this)[12]*(*this)[4] - (*this)[14]*(*this)[2] + (*this)[15]*(*this)[9] - (*this)[1]*(*this)[3] - (*this)[7]*(*this)[8]);
    res[14]=x0*((*this)[0]*(*this)[7] + (*this)[10]*(*this)[15] + (*this)[11]*(*this)[14] - (*this)[12]*(*this)[3] + (*this)[13]*(*this)[2] - (*this)[1]*(*this)[4] - (*this)[5]*(*this)[9] + (*this)[6]*(*this)[8]);
    res[15]=x0*(-(*this)[0]*(*this)[1] - (*this)[10]*(*this)[14] + (*this)[11]*(*this)[15] - (*this)[12]*(*this)[8] - (*this)[13]*(*this)[9] - (*this)[2]*(*this)[5] - (*this)[3]*(*this)[6] - (*this)[4]*(*this)[7]);
    return res;
};

//-----------------------------------
// (R130MV, R130B4) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15] + B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8];
    res[9]=A[9];
    res[10]=A[10];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15] - B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=-A[15]*B[0];
    res[1]=-A[11]*B[0];
    res[2]=A[12]*B[0];
    res[3]=A[13]*B[0];
    res[4]=A[14]*B[0];
    res[5]=-A[8]*B[0];
    res[6]=-A[9]*B[0];
    res[7]=-A[10]*B[0];
    res[8]=A[5]*B[0];
    res[9]=A[6]*B[0];
    res[10]=A[7]*B[0];
    res[11]=A[1]*B[0];
    res[12]=-A[2]*B[0];
    res[13]=-A[3]*B[0];
    res[14]=-A[4]*B[0];
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const R130MV<T> &A, const R130B4<T> &B) {
    Rotor<T> res;
    res[0]=-A[15]*B[0];
    res[1]=-A[8]*B[0];
    res[2]=-A[9]*B[0];
    res[3]=-A[10]*B[0];
    res[4]=A[5]*B[0];
    res[5]=A[6]*B[0];
    res[6]=A[7]*B[0];
    res[7]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[11]*B[0];
    res[2]=A[12]*B[0];
    res[3]=A[13]*B[0];
    res[4]=A[14]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[1]*B[0];
    res[12]=-A[2]*B[0];
    res[13]=-A[3]*B[0];
    res[14]=-A[4]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130B4<T> &A) const {
    R130MV<T> res;
    T x0 = 2.0*A[0];
    res[0]=x0*(-(*this)[0]*(*this)[15] + (*this)[10]*(*this)[7] - (*this)[11]*(*this)[1] - (*this)[12]*(*this)[2] - (*this)[13]*(*this)[3] - (*this)[14]*(*this)[4] + (*this)[5]*(*this)[8] + (*this)[6]*(*this)[9]);
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
    res[11]=x0*(-(*this)[0]*(*this)[1] - (*this)[10]*(*this)[14] - (*this)[11]*(*this)[15] - (*this)[12]*(*this)[8] - (*this)[13]*(*this)[9] + (*this)[2]*(*this)[5] + (*this)[3]*(*this)[6] + (*this)[4]*(*this)[7]);
    res[12]=x0*((*this)[0]*(*this)[2] - (*this)[10]*(*this)[3] - (*this)[11]*(*this)[8] - (*this)[12]*(*this)[15] - (*this)[13]*(*this)[7] + (*this)[14]*(*this)[6] - (*this)[1]*(*this)[5] + (*this)[4]*(*this)[9]);
    res[13]=x0*((*this)[0]*(*this)[3] + (*this)[10]*(*this)[2] - (*this)[11]*(*this)[9] + (*this)[12]*(*this)[7] - (*this)[13]*(*this)[15] - (*this)[14]*(*this)[5] - (*this)[1]*(*this)[6] - (*this)[4]*(*this)[8]);
    res[14]=x0*((*this)[0]*(*this)[4] - (*this)[10]*(*this)[11] - (*this)[12]*(*this)[6] + (*this)[13]*(*this)[5] - (*this)[14]*(*this)[15] - (*this)[1]*(*this)[7] - (*this)[2]*(*this)[9] + (*this)[3]*(*this)[8]);
    res[15]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[10], 2) - std::pow((*this)[11], 2) + std::pow((*this)[12], 2) + std::pow((*this)[13], 2) + std::pow((*this)[14], 2) - std::pow((*this)[15], 2) + std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2) - std::pow((*this)[4], 2) - std::pow((*this)[5], 2) - std::pow((*this)[6], 2) - std::pow((*this)[7], 2) + std::pow((*this)[8], 2) + std::pow((*this)[9], 2));
    return res;
};

//-----------------------------------
// (R130MV, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    res[4]=A[4] + B[4];
    res[5]=A[5] + B[5];
    res[6]=A[6] + B[6];
    res[7]=A[7] + B[7];
    res[8]=A[8] + B[8];
    res[9]=A[9] + B[9];
    res[10]=A[10] + B[10];
    res[11]=A[11] + B[11];
    res[12]=A[12] + B[12];
    res[13]=A[13] + B[13];
    res[14]=A[14] + B[14];
    res[15]=A[15] + B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    res[4]=A[4] - B[4];
    res[5]=A[5] - B[5];
    res[6]=A[6] - B[6];
    res[7]=A[7] - B[7];
    res[8]=A[8] - B[8];
    res[9]=A[9] - B[9];
    res[10]=A[10] - B[10];
    res[11]=A[11] - B[11];
    res[12]=A[12] - B[12];
    res[13]=A[13] - B[13];
    res[14]=A[14] - B[14];
    res[15]=A[15] - B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[10]*B[10] + A[11]*B[11] - A[12]*B[12] - A[13]*B[13] - A[14]*B[14] - A[15]*B[15] + A[1]*B[1] - A[2]*B[2] - A[3]*B[3] - A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7] - A[8]*B[8] - A[9]*B[9];
    res[1]=A[0]*B[1] - A[10]*B[14] - A[11]*B[15] - A[12]*B[8] - A[13]*B[9] - A[14]*B[10] + A[15]*B[11] + A[1]*B[0] - A[2]*B[5] - A[3]*B[6] - A[4]*B[7] + A[5]*B[2] + A[6]*B[3] + A[7]*B[4] - A[8]*B[12] - A[9]*B[13];
    res[2]=A[0]*B[2] + A[10]*B[3] + A[11]*B[8] + A[12]*B[15] + A[13]*B[7] - A[14]*B[6] - A[15]*B[12] - A[1]*B[5] + A[2]*B[0] - A[3]*B[10] + A[4]*B[9] + A[5]*B[1] - A[6]*B[14] + A[7]*B[13] + A[8]*B[11] - A[9]*B[4];
    res[3]=A[0]*B[3] - A[10]*B[2] + A[11]*B[9] - A[12]*B[7] + A[13]*B[15] + A[14]*B[5] - A[15]*B[13] - A[1]*B[6] + A[2]*B[10] + A[3]*B[0] - A[4]*B[8] + A[5]*B[14] + A[6]*B[1] - A[7]*B[12] + A[8]*B[4] + A[9]*B[11];
    res[4]=A[0]*B[4] + A[10]*B[11] + A[11]*B[10] + A[12]*B[6] - A[13]*B[5] + A[14]*B[15] - A[15]*B[14] - A[1]*B[7] - A[2]*B[9] + A[3]*B[8] + A[4]*B[0] - A[5]*B[13] + A[6]*B[12] + A[7]*B[1] - A[8]*B[3] + A[9]*B[2];
    res[5]=A[0]*B[5] + A[10]*B[6] + A[11]*B[12] - A[12]*B[11] + A[13]*B[4] - A[14]*B[3] - A[15]*B[8] - A[1]*B[2] + A[2]*B[1] - A[3]*B[14] + A[4]*B[13] + A[5]*B[0] - A[6]*B[10] + A[7]*B[9] - A[8]*B[15] - A[9]*B[7];
    res[6]=A[0]*B[6] - A[10]*B[5] + A[11]*B[13] - A[12]*B[4] - A[13]*B[11] + A[14]*B[2] - A[15]*B[9] - A[1]*B[3] + A[2]*B[14] + A[3]*B[1] - A[4]*B[12] + A[5]*B[10] + A[6]*B[0] - A[7]*B[8] + A[8]*B[7] - A[9]*B[15];
    res[7]=A[0]*B[7] - A[10]*B[15] + A[11]*B[14] + A[12]*B[3] - A[13]*B[2] - A[14]*B[11] - A[15]*B[10] - A[1]*B[4] - A[2]*B[13] + A[3]*B[12] + A[4]*B[1] - A[5]*B[9] + A[6]*B[8] + A[7]*B[0] - A[8]*B[6] + A[9]*B[5];
    res[8]=A[0]*B[8] + A[10]*B[9] + A[11]*B[2] + A[12]*B[1] - A[13]*B[14] + A[14]*B[13] + A[15]*B[5] + A[1]*B[12] + A[2]*B[11] - A[3]*B[4] + A[4]*B[3] + A[5]*B[15] + A[6]*B[7] - A[7]*B[6] + A[8]*B[0] - A[9]*B[10];
    res[9]=A[0]*B[9] - A[10]*B[8] + A[11]*B[3] + A[12]*B[14] + A[13]*B[1] - A[14]*B[12] + A[15]*B[6] + A[1]*B[13] + A[2]*B[4] + A[3]*B[11] - A[4]*B[2] - A[5]*B[7] + A[6]*B[15] + A[7]*B[5] + A[8]*B[10] + A[9]*B[0];
    res[10]=A[0]*B[10] + A[10]*B[0] + A[11]*B[4] - A[12]*B[13] + A[13]*B[12] + A[14]*B[1] + A[15]*B[7] + A[1]*B[14] - A[2]*B[3] + A[3]*B[2] + A[4]*B[11] + A[5]*B[6] - A[6]*B[5] + A[7]*B[15] - A[8]*B[9] + A[9]*B[8];
    res[11]=A[0]*B[11] - A[10]*B[4] + A[11]*B[0] + A[12]*B[5] + A[13]*B[6] + A[14]*B[7] - A[15]*B[1] + A[1]*B[15] - A[2]*B[8] - A[3]*B[9] - A[4]*B[10] - A[5]*B[12] - A[6]*B[13] - A[7]*B[14] - A[8]*B[2] - A[9]*B[3];
    res[12]=A[0]*B[12] + A[10]*B[13] + A[11]*B[5] + A[12]*B[0] - A[13]*B[10] + A[14]*B[9] + A[15]*B[2] + A[1]*B[8] - A[2]*B[15] - A[3]*B[7] + A[4]*B[6] - A[5]*B[11] + A[6]*B[4] - A[7]*B[3] + A[8]*B[1] - A[9]*B[14];
    res[13]=A[0]*B[13] - A[10]*B[12] + A[11]*B[6] + A[12]*B[10] + A[13]*B[0] - A[14]*B[8] + A[15]*B[3] + A[1]*B[9] + A[2]*B[7] - A[3]*B[15] - A[4]*B[5] - A[5]*B[4] - A[6]*B[11] + A[7]*B[2] + A[8]*B[14] + A[9]*B[1];
    res[14]=A[0]*B[14] + A[10]*B[1] + A[11]*B[7] - A[12]*B[9] + A[13]*B[8] + A[14]*B[0] + A[15]*B[4] + A[1]*B[10] - A[2]*B[6] + A[3]*B[5] - A[4]*B[15] + A[5]*B[3] - A[6]*B[2] - A[7]*B[11] - A[8]*B[13] + A[9]*B[12];
    res[15]=A[0]*B[15] + A[10]*B[7] - A[11]*B[1] - A[12]*B[2] - A[13]*B[3] - A[14]*B[4] + A[15]*B[0] + A[1]*B[11] + A[2]*B[12] + A[3]*B[13] + A[4]*B[14] + A[5]*B[8] + A[6]*B[9] + A[7]*B[10] + A[8]*B[5] + A[9]*B[6];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[10]*B[10] + A[11]*B[11] - A[12]*B[12] - A[13]*B[13] - A[14]*B[14] - A[15]*B[15] + A[1]*B[1] - A[2]*B[2] - A[3]*B[3] - A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7] - A[8]*B[8] - A[9]*B[9];
    res[1]=A[0]*B[1] - A[10]*B[14] - A[12]*B[8] - A[13]*B[9] - A[14]*B[10] + A[1]*B[0] - A[8]*B[12] - A[9]*B[13];
    res[2]=A[0]*B[2] + A[11]*B[8] + A[13]*B[7] - A[14]*B[6] + A[2]*B[0] - A[6]*B[14] + A[7]*B[13] + A[8]*B[11];
    res[3]=A[0]*B[3] + A[11]*B[9] - A[12]*B[7] + A[14]*B[5] + A[3]*B[0] + A[5]*B[14] - A[7]*B[12] + A[9]*B[11];
    res[4]=A[0]*B[4] + A[10]*B[11] + A[11]*B[10] + A[12]*B[6] - A[13]*B[5] + A[4]*B[0] - A[5]*B[13] + A[6]*B[12];
    res[5]=A[0]*B[5] + A[13]*B[4] - A[14]*B[3] - A[15]*B[8] - A[3]*B[14] + A[4]*B[13] + A[5]*B[0] - A[8]*B[15];
    res[6]=A[0]*B[6] - A[12]*B[4] + A[14]*B[2] - A[15]*B[9] + A[2]*B[14] - A[4]*B[12] + A[6]*B[0] - A[9]*B[15];
    res[7]=A[0]*B[7] - A[10]*B[15] + A[12]*B[3] - A[13]*B[2] - A[15]*B[10] - A[2]*B[13] + A[3]*B[12] + A[7]*B[0];
    res[8]=A[0]*B[8] + A[11]*B[2] + A[12]*B[1] + A[15]*B[5] + A[1]*B[12] + A[2]*B[11] + A[5]*B[15] + A[8]*B[0];
    res[9]=A[0]*B[9] + A[11]*B[3] + A[13]*B[1] + A[15]*B[6] + A[1]*B[13] + A[3]*B[11] + A[6]*B[15] + A[9]*B[0];
    res[10]=A[0]*B[10] + A[10]*B[0] + A[11]*B[4] + A[14]*B[1] + A[15]*B[7] + A[1]*B[14] + A[4]*B[11] + A[7]*B[15];
    res[11]=A[0]*B[11] - A[10]*B[4] + A[11]*B[0] - A[2]*B[8] - A[3]*B[9] - A[4]*B[10] - A[8]*B[2] - A[9]*B[3];
    res[12]=A[0]*B[12] + A[12]*B[0] + A[1]*B[8] - A[3]*B[7] + A[4]*B[6] + A[6]*B[4] - A[7]*B[3] + A[8]*B[1];
    res[13]=A[0]*B[13] + A[13]*B[0] + A[1]*B[9] + A[2]*B[7] - A[4]*B[5] - A[5]*B[4] + A[7]*B[2] + A[9]*B[1];
    res[14]=A[0]*B[14] + A[10]*B[1] + A[14]*B[0] + A[1]*B[10] - A[2]*B[6] + A[3]*B[5] + A[5]*B[3] - A[6]*B[2];
    res[15]=A[0]*B[15] + A[10]*B[7] + A[15]*B[0] + A[5]*B[8] + A[6]*B[9] + A[7]*B[10] + A[8]*B[5] + A[9]*B[6];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[11]*B[15] + A[15]*B[11] - A[2]*B[5] - A[3]*B[6] - A[4]*B[7] + A[5]*B[2] + A[6]*B[3] + A[7]*B[4];
    res[2]=A[10]*B[3] + A[12]*B[15] - A[15]*B[12] - A[1]*B[5] - A[3]*B[10] + A[4]*B[9] + A[5]*B[1] - A[9]*B[4];
    res[3]=-A[10]*B[2] + A[13]*B[15] - A[15]*B[13] - A[1]*B[6] + A[2]*B[10] - A[4]*B[8] + A[6]*B[1] + A[8]*B[4];
    res[4]=A[14]*B[15] - A[15]*B[14] - A[1]*B[7] - A[2]*B[9] + A[3]*B[8] + A[7]*B[1] - A[8]*B[3] + A[9]*B[2];
    res[5]=A[10]*B[6] + A[11]*B[12] - A[12]*B[11] - A[1]*B[2] + A[2]*B[1] - A[6]*B[10] + A[7]*B[9] - A[9]*B[7];
    res[6]=-A[10]*B[5] + A[11]*B[13] - A[13]*B[11] - A[1]*B[3] + A[3]*B[1] + A[5]*B[10] - A[7]*B[8] + A[8]*B[7];
    res[7]=A[11]*B[14] - A[14]*B[11] - A[1]*B[4] + A[4]*B[1] - A[5]*B[9] + A[6]*B[8] - A[8]*B[6] + A[9]*B[5];
    res[8]=A[10]*B[9] - A[13]*B[14] + A[14]*B[13] - A[3]*B[4] + A[4]*B[3] + A[6]*B[7] - A[7]*B[6] - A[9]*B[10];
    res[9]=-A[10]*B[8] + A[12]*B[14] - A[14]*B[12] + A[2]*B[4] - A[4]*B[2] - A[5]*B[7] + A[7]*B[5] + A[8]*B[10];
    res[10]=-A[12]*B[13] + A[13]*B[12] - A[2]*B[3] + A[3]*B[2] + A[5]*B[6] - A[6]*B[5] - A[8]*B[9] + A[9]*B[8];
    res[11]=A[12]*B[5] + A[13]*B[6] + A[14]*B[7] - A[15]*B[1] + A[1]*B[15] - A[5]*B[12] - A[6]*B[13] - A[7]*B[14];
    res[12]=A[10]*B[13] + A[11]*B[5] - A[13]*B[10] + A[14]*B[9] + A[15]*B[2] - A[2]*B[15] - A[5]*B[11] - A[9]*B[14];
    res[13]=-A[10]*B[12] + A[11]*B[6] + A[12]*B[10] - A[14]*B[8] + A[15]*B[3] - A[3]*B[15] - A[6]*B[11] + A[8]*B[14];
    res[14]=A[11]*B[7] - A[12]*B[9] + A[13]*B[8] + A[15]*B[4] - A[4]*B[15] - A[7]*B[11] - A[8]*B[13] + A[9]*B[12];
    res[15]=-A[11]*B[1] - A[12]*B[2] - A[13]*B[3] - A[14]*B[4] + A[1]*B[11] + A[2]*B[12] + A[3]*B[13] + A[4]*B[14];
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[10], 2);
    T x2 = std::pow((*this)[11], 2);
    T x3 = std::pow((*this)[2], 2);
    T x4 = std::pow((*this)[3], 2);
    T x5 = std::pow((*this)[4], 2);
    T x6 = std::pow((*this)[8], 2);
    T x7 = std::pow((*this)[9], 2);
    T x8 = std::pow((*this)[12], 2);
    T x9 = std::pow((*this)[13], 2);
    T x10 = std::pow((*this)[14], 2);
    T x11 = std::pow((*this)[15], 2);
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[5], 2);
    T x14 = std::pow((*this)[6], 2);
    T x15 = std::pow((*this)[7], 2);
    T x16 = 2.0*A[11];
    T x17 = (*this)[0]*(*this)[11];
    T x18 = 2.0*(*this)[10];
    T x19 = (*this)[12]*x18;
    T x20 = (*this)[1]*x18;
    T x21 = (*this)[10]*x16;
    T x22 = (*this)[7]*x18;
    T x23 = 2.0*(*this)[11];
    T x24 = (*this)[5]*x23;
    T x25 = (*this)[6]*x23;
    T x26 = (*this)[7]*x23;
    T x27 = 2.0*A[14];
    T x28 = (*this)[13]*(*this)[8];
    T x29 = 2.0*A[12];
    T x30 = (*this)[14]*(*this)[9];
    T x31 = (*this)[15]*(*this)[1];
    T x32 = (*this)[15]*(*this)[2];
    T x33 = 2.0*A[13];
    T x34 = (*this)[15]*(*this)[3];
    T x35 = (*this)[15]*(*this)[4];
    T x36 = (*this)[1]*(*this)[8];
    T x37 = (*this)[1]*(*this)[9];
    T x38 = (*this)[2]*(*this)[7];
    T x39 = (*this)[2]*(*this)[8];
    T x40 = (*this)[3]*(*this)[5];
    T x41 = (*this)[3]*(*this)[9];
    T x42 = (*this)[4]*(*this)[6];
    T x43 = 2.0*A[15];
    T x44 = (*this)[5]*(*this)[8];
    T x45 = (*this)[6]*(*this)[9];
    T x46 = (*this)[0]*(*this)[12];
    T x47 = (*this)[0]*(*this)[13];
    T x48 = (*this)[0]*(*this)[14];
    T x49 = (*this)[0]*(*this)[15];
    T x50 = (*this)[13]*x18;
    T x51 = (*this)[1]*x23;
    T x52 = (*this)[12]*(*this)[2];
    T x53 = (*this)[12]*(*this)[5];
    T x54 = (*this)[12]*(*this)[9];
    T x55 = (*this)[13]*(*this)[3];
    T x56 = (*this)[13]*(*this)[6];
    T x57 = (*this)[14]*(*this)[4];
    T x58 = (*this)[14]*(*this)[7];
    T x59 = (*this)[14]*(*this)[8];
    T x60 = (*this)[2]*(*this)[6];
    T x61 = (*this)[3]*(*this)[7];
    T x62 = (*this)[4]*(*this)[5];
    T x63 = 2.0*(*this)[0];
    T x64 = A[5]*x63;
    T x65 = (*this)[3]*x63;
    T x66 = (*this)[4]*x63;
    T x67 = (*this)[3]*A[5];
    T x68 = (*this)[5]*A[3];
    T x69 = 2.0*(*this)[12];
    T x70 = (*this)[4]*A[3];
    T x71 = (*this)[6]*x69;
    T x72 = 2.0*(*this)[13];
    T x73 = (*this)[2]*A[4];
    T x74 = (*this)[7]*A[5];
    T x75 = 2.0*A[10];
    T x76 = 2.0*(*this)[14];
    T x77 = (*this)[3]*A[2];
    T x78 = (*this)[5]*A[6];
    T x79 = 2.0*A[8];
    T x80 = 2.0*(*this)[1];
    T x81 = (*this)[2]*A[2];
    T x82 = (*this)[3]*A[3];
    T x83 = (*this)[4]*A[4];
    T x84 = 2.0*A[9];
    T x85 = 2.0*(*this)[9];
    T x86 = (*this)[2]*A[7];
    T x87 = 2.0*(*this)[8];
    T x88 = (*this)[4]*A[6];
    T x89 = (*this)[6]*A[4];
    T x90 = (*this)[7]*x85;
    T x91 = (*this)[5]*x63;
    T x92 = (*this)[6]*x63;
    T x93 = A[4]*x63;
    T x94 = (*this)[11]*x18;
    T x95 = A[4]*x18;
    T x96 = (*this)[2]*x18;
    T x97 = (*this)[6]*x18;
    T x98 = (*this)[12]*x23;
    T x99 = (*this)[13]*A[3];
    T x100 = A[4]*x23;
    T x101 = A[5]*x23;
    T x102 = (*this)[9]*x23;
    T x103 = A[5]*x69;
    T x104 = A[4]*x69;
    T x105 = (*this)[7]*x69;
    T x106 = (*this)[15]*x72;
    T x107 = (*this)[4]*x72;
    T x108 = (*this)[5]*x72;
    T x109 = (*this)[15]*x76;
    T x110 = A[3]*x76;
    T x111 = (*this)[6]*x76;
    T x112 = (*this)[15]*x87;
    T x113 = (*this)[15]*A[3];
    T x114 = (*this)[5]*A[5];
    T x115 = (*this)[6]*A[6];
    T x116 = (*this)[7]*x80;
    T x117 = (*this)[3]*A[7];
    T x118 = (*this)[4]*x85;
    T x119 = (*this)[5]*A[4];
    T x120 = (*this)[7]*x87;
    T x121 = 2.0*A[7];
    T x122 = 2.0*A[6];
    T x123 = A[4]*x72;
    T x124 = A[1]*x76;
    T x125 = (*this)[14]*x75;
    T x126 = 2.0*(*this)[7];
    T x127 = (*this)[1]*x75;
    T x128 = 2.0*(*this)[2];
    T x129 = 2.0*(*this)[4];
    T x130 = (*this)[2]*(*this)[9];
    T x131 = 2.0*(*this)[6];
    T x132 = (*this)[3]*(*this)[8];
    T x133 = (*this)[4]*x75;
    T x134 = (*this)[8]*x85;
    T x135 = (*this)[0]*x18;
    T x136 = (*this)[14]*x18;
    T x137 = (*this)[4]*x18;
    T x138 = (*this)[9]*x72;
    T x139 = 2.0*(*this)[15];
    T x140 = A[1]*x80;
    T x141 = 2.0*A[5];
    T x142 = (*this)[1]*x63;
    T x143 = (*this)[3]*x18;
    T x144 = A[1]*x23;
    T x145 = (*this)[15]*x23;
    T x146 = (*this)[4]*A[2];
    T x147 = A[2]*x69;
    T x148 = (*this)[15]*x69;
    T x149 = A[1]*x69;
    T x150 = (*this)[7]*x72;
    T x151 = (*this)[1]*A[2];
    T x152 = (*this)[1]*x79;
    T x153 = 2.0*(*this)[3];
    T x154 = A[2]*x131;
    T x155 = (*this)[2]*x63;
    T x156 = (*this)[5]*A[1];
    T x157 = (*this)[8]*x23;
    T x158 = (*this)[8]*x69;
    T x159 = A[2]*x126;
    T x160 = (*this)[15]*A[1];
    T x161 = (*this)[8]*x63;
    T x162 = (*this)[8]*x18;
    T x163 = (*this)[9]*x18;
    T x164 = (*this)[2]*x23;
    T x165 = (*this)[1]*x69;
    T x166 = (*this)[14]*x72;
    T x167 = (*this)[2]*x72;
    T x168 = (*this)[5]*A[8];
    T x169 = A[9]*x80;
    T x170 = 2.0*(*this)[5];
    T x171 = (*this)[8]*x79;
    T x172 = (*this)[9]*x84;
    T x173 = (*this)[7]*x63;
    T x174 = (*this)[9]*x63;
    T x175 = (*this)[7]*x121;
    T x176 = A[1]*x87;
    T x177 = 2.0*A[4];
    T x178 = (*this)[15]*x18;
    T x179 = (*this)[5]*x18;
    T x180 = (*this)[13]*x23;
    T x181 = A[6]*x69;
    T x182 = A[7]*x69;
    T x183 = (*this)[3]*x69;
    T x184 = (*this)[4]*x69;
    T x185 = 2.0*A[3];
    T x186 = (*this)[1]*x72;
    T x187 = (*this)[15]*x121;
    T x188 = (*this)[8]*x84;
    T x189 = (*this)[8]*x75;
    T x190 = (*this)[14]*x23;
    T x191 = (*this)[1]*x76;
    T x192 = (*this)[15]*(*this)[7];
    T x193 = (*this)[9]*x75;
    T x194 = (*this)[2]*(*this)[3];
    T x195 = 2.0*A[2];
    T x196 = (*this)[9]*x79;
    T x197 = (*this)[2]*x76;
    T x198 = (*this)[3]*x76;
    T x199 = (*this)[2]*x85;
    T x200 = (*this)[15]*A[5];
    T x201 = (*this)[4]*A[5];
    T x202 = 2.0*A[1];
    T x203 = (*this)[3]*x23;
    T x204 = (*this)[15]*x75;
    T x205 = (*this)[4]*x23;
    T x206 = (*this)[13]*x69;
    T x207 = (*this)[14]*x69;
    T x208 = (*this)[5]*(*this)[6];
    T x209 = (*this)[7]*x75;
    T x210 = (*this)[15]*(*this)[5];
    T x211 = (*this)[2]*(*this)[4];
    T x212 = (*this)[3]*(*this)[4];
    T x213 = (*this)[15]*(*this)[6];
    T x214 = (*this)[5]*(*this)[7];
    T x215 = (*this)[6]*(*this)[7];
    T x216 = 2.0*A[0];
    T x217 = (*this)[0]*(*this)[5];
    T x218 = (*this)[0]*(*this)[6];
    T x219 = (*this)[0]*x27;
    T x220 = (*this)[12]*x33;
    T x221 = (*this)[13]*x27;
    T x222 = (*this)[14]*(*this)[3];
    T x223 = (*this)[15]*(*this)[8];
    T x224 = (*this)[9]*x33;
    T x225 = (*this)[1]*(*this)[2];
    T x226 = (*this)[1]*x33;
    T x227 = (*this)[1]*x27;
    T x228 = (*this)[2]*x43;
    T x229 = (*this)[6]*x43;
    T x230 = (*this)[4]*x43;
    T x231 = (*this)[5]*(*this)[9];
    T x232 = (*this)[7]*(*this)[8];
    T x233 = (*this)[1]*x43;
    T x234 = (*this)[12]*x27;
    T x235 = (*this)[12]*x43;
    T x236 = (*this)[13]*(*this)[4];
    T x237 = (*this)[13]*x43;
    T x238 = (*this)[14]*(*this)[2];
    T x239 = (*this)[6]*(*this)[8];
    T x240 = (*this)[7]*(*this)[9];
    T x241 = (*this)[11]*x16;
    T x242 = (*this)[0]*(*this)[3];
    T x243 = (*this)[12]*x16;
    T x244 = (*this)[9]*x16;
    T x245 = (*this)[9]*x29;
    T x246 = (*this)[13]*x29;
    T x247 = (*this)[14]*x29;
    T x248 = (*this)[14]*x43;
    T x249 = (*this)[1]*x16;
    T x250 = (*this)[7]*x16;
    T x251 = (*this)[8]*x33;
    T x252 = (*this)[2]*x16;
    T x253 = (*this)[14]*x33;
    res[0]=(*this)[4]*x21 + A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - A[0]*x9 + A[12]*x24 - A[12]*x50 + A[13]*x19 + A[13]*x25 + A[14]*x20 + A[14]*x26 + A[15]*x22 - A[15]*x51 + x16*x17 + x16*x31 + x16*x39 + x16*x41 - x16*x53 - x16*x56 - x16*x58 + x27*x28 + x27*x35 + x27*x40 - x27*x48 - x27*x54 - x27*x60 + x29*x30 + x29*x32 + x29*x36 + x29*x42 - x29*x46 - x29*x61 + x33*x34 + x33*x37 + x33*x38 - x33*x47 - x33*x59 - x33*x62 + x43*x44 + x43*x45 - x43*x49 - x43*x52 - x43*x55 - x43*x57;
    res[1]=-(*this)[14]*x100 - (*this)[15]*x103 - (*this)[15]*x95 - (*this)[2]*x110 + (*this)[2]*x64 - (*this)[3]*x104 - (*this)[7]*x93 - (*this)[8]*x101 + A[10]*x20 - A[10]*x26 + A[1]*x0 + A[1]*x1 - A[1]*x10 + A[1]*x11 - A[1]*x12 + A[1]*x13 + A[1]*x14 + A[1]*x15 - A[1]*x2 - A[1]*x3 - A[1]*x4 - A[1]*x5 + A[1]*x6 + A[1]*x7 - A[1]*x8 - A[1]*x9 - A[2]*x107 - A[2]*x112 + A[2]*x90 - A[2]*x91 - A[2]*x97 - A[2]*x98 - A[3]*x120 - A[3]*x92 - A[5]*x111 - A[5]*x118 - A[6]*x102 - A[6]*x105 - A[6]*x106 + A[6]*x65 - A[6]*x96 - A[7]*x108 - A[7]*x109 - A[7]*x116 + A[7]*x66 + A[7]*x71 - A[7]*x94 - A[8]*x24 - A[8]*x50 + A[9]*x19 - A[9]*x25 - x113*x85 - x114*x80 - x115*x80 - x117*x87 - x119*x85 + x18*x67 + x18*x68 - x23*x99 + x28*x75 + x30*x79 - x32*x79 - x34*x84 - x35*x75 + x36*x79 + x37*x84 - x38*x84 - x40*x75 - x42*x79 - x46*x79 - x47*x84 - x48*x75 - x54*x75 - x59*x84 + x60*x75 + x61*x79 + x62*x84 + x69*x70 + x72*x73 + x72*x74 + x76*x77 + x76*x78 + x80*x81 + x80*x82 + x80*x83 + x85*x86 + x87*x88 + x87*x89;
    res[2]=(*this)[14]*x104 + (*this)[15]*x101 - (*this)[1]*x110 + (*this)[1]*x123 + (*this)[1]*x64 - (*this)[2]*x140 + (*this)[3]*x100 + (*this)[3]*x124 + (*this)[5]*x125 + (*this)[6]*x127 + (*this)[8]*x103 + (*this)[8]*x133 + (*this)[8]*x95 + (*this)[9]*x93 + A[10]*x102 + A[10]*x105 - A[10]*x106 + A[10]*x65 + A[10]*x96 - A[1]*x107 - A[1]*x112 - A[1]*x90 - A[1]*x91 + A[1]*x97 + A[1]*x98 + A[2]*x0 - A[2]*x1 - A[2]*x10 + A[2]*x11 + A[2]*x12 + A[2]*x13 - A[2]*x14 - A[2]*x15 + A[2]*x2 + A[2]*x3 - A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7 + A[2]*x8 - A[2]*x9 + A[3]*x134 - A[3]*x135 - A[5]*x136 - A[5]*x138 - A[6]*x20 + A[6]*x26 + A[7]*x19 - A[7]*x25 - A[8]*x137 + A[9]*x108 + A[9]*x109 - A[9]*x116 - A[9]*x66 + A[9]*x71 - A[9]*x94 + x113*x126 - x114*x128 + x119*x126 + x121*x34 + x121*x37 - x121*x38 + x121*x47 + x121*x59 - x121*x62 + x122*x28 - x122*x35 - x122*x40 - x122*x48 + x122*x54 - x122*x60 + x128*x82 + x129*x73 + x129*x74 + x130*x84 + x131*x67 + x131*x68 + x132*x84 - x139*x89 + x17*x79 - x23*x70 - x31*x79 + x39*x79 - x41*x79 + x53*x79 - x56*x79 - x58*x79 + x69*x99;
    res[3]=(*this)[13]*x144 + (*this)[13]*x147 + (*this)[14]*x123 - (*this)[15]*x159 - (*this)[1]*x104 - (*this)[2]*x124 - (*this)[3]*x140 + (*this)[4]*x149 - (*this)[5]*x127 + (*this)[5]*x154 + (*this)[6]*x125 + (*this)[7]*x152 - (*this)[8]*x93 + (*this)[9]*x133 + (*this)[9]*x95 + A[10]*x143 + A[10]*x148 + A[10]*x150 - A[10]*x155 - A[10]*x157 + A[1]*x120 - A[1]*x92 + A[2]*x134 + A[2]*x135 + A[3]*x0 - A[3]*x1 - A[3]*x10 + A[3]*x11 + A[3]*x12 - A[3]*x13 + A[3]*x14 - A[3]*x15 + A[3]*x2 - A[3]*x3 + A[3]*x4 - A[3]*x5 - A[3]*x6 + A[3]*x7 - A[3]*x8 + A[3]*x9 + A[5]*x20 - A[5]*x26 - A[6]*x136 + A[6]*x138 + A[6]*x142 + A[6]*x145 - A[6]*x158 + A[7]*x24 + A[7]*x50 + A[8]*x108 - A[8]*x109 + A[8]*x66 + A[8]*x71 + A[8]*x94 - A[9]*x137 - x115*x153 + x119*x139 + x121*x30 - x121*x32 - x121*x36 - x121*x42 - x121*x46 - x121*x61 + x126*x88 + x126*x89 + x128*x77 + x128*x78 + x130*x79 + x132*x79 + x141*x28 + x141*x35 - x141*x40 + x141*x48 + x141*x54 - x141*x60 + x146*x23 + x151*x76 + x153*x83 - x156*x18 - x160*x85 + x17*x84 - x23*x73 - x31*x84 - x39*x84 + x41*x84 - x53*x84 + x56*x84 - x58*x84;
    res[4]=(*this)[14]*x144 + (*this)[14]*x147 + (*this)[15]*x154 - (*this)[3]*x149 - (*this)[4]*x140 + (*this)[4]*x171 + (*this)[4]*x172 - (*this)[4]*x175 + (*this)[5]*x159 + (*this)[5]*x169 + (*this)[6]*A[3]*x126 - (*this)[6]*x152 - (*this)[6]*x176 + A[10]*x137 + A[1]*x167 - A[1]*x173 + A[2]*x162 - A[2]*x174 + A[3]*x161 + A[3]*x163 + A[3]*x164 + A[3]*x165 + A[3]*x166 + A[4]*x0 + A[4]*x1 + A[4]*x10 + A[4]*x11 + A[4]*x12 - A[4]*x13 - A[4]*x14 + A[4]*x15 + A[4]*x2 - A[4]*x3 - A[4]*x4 + A[4]*x5 - A[4]*x6 - A[4]*x7 - A[4]*x8 - A[4]*x9 + A[5]*x19 + A[5]*x25 - A[6]*x24 + A[6]*x50 + A[7]*x136 - A[7]*x138 + A[7]*x142 + A[7]*x145 - A[7]*x158 - A[8]*x102 + A[8]*x105 + A[8]*x106 - A[8]*x65 + A[8]*x96 + A[9]*x111 + A[9]*x143 - A[9]*x148 + A[9]*x150 + A[9]*x155 + A[9]*x157 + x117*x131 + x122*x30 + x122*x32 + x122*x36 - x122*x42 + x122*x46 - x122*x61 + x129*x81 - x139*x68 - x141*x34 - x141*x37 - x141*x38 - x141*x47 + x141*x59 - x141*x62 - x151*x72 + x153*x70 + x156*x85 - x160*x18 + x168*x76 + x17*x75 + x170*x86 - x23*x77 - x31*x75 - x39*x75 - x41*x75 - x53*x75 - x56*x75 + x58*x75;
    res[5]=(*this)[13]*x181 + (*this)[14]*x182 - (*this)[15]*x193 + (*this)[2]*x125 + (*this)[3]*x127 - (*this)[4]*x169 - (*this)[5]*x140 + (*this)[5]*x172 - (*this)[5]*x175 + (*this)[6]*x187 + (*this)[6]*x188 + (*this)[7]*x189 + (*this)[8]*x144 + (*this)[8]*x147 + A[10]*x179 + A[10]*x180 + A[10]*x184 + A[10]*x92 + A[1]*x111 - A[1]*x118 + A[1]*x143 - A[1]*x148 - A[1]*x150 - A[1]*x155 - A[2]*x136 - A[2]*x138 + A[2]*x142 - A[2]*x145 - A[3]*x20 - A[3]*x26 + A[4]*x19 + A[4]*x25 + A[5]*x0 - A[5]*x1 - A[5]*x10 - A[5]*x11 + A[5]*x12 - A[5]*x13 + A[5]*x14 + A[5]*x15 - A[5]*x2 - A[5]*x3 + A[5]*x4 + A[5]*x5 + A[5]*x6 - A[5]*x7 + A[5]*x8 - A[5]*x9 + A[6]*x134 - A[6]*x135 - A[6]*x191 + A[7]*x162 + A[7]*x174 + A[7]*x186 - A[8]*x22 + A[8]*x51 + A[9]*x167 - A[9]*x173 + A[9]*x178 + A[9]*x183 - A[9]*x190 - x117*x23 - x122*x192 - x122*x194 - x126*x146 - x129*x86 - x131*x77 - x131*x78 + x170*x81 - x177*x34 + x177*x37 + x177*x38 + x177*x47 + x177*x59 + x177*x62 + x185*x28 + x185*x35 + x185*x40 - x185*x48 + x185*x54 + x185*x60 + x23*x88 + x44*x79 - x45*x79 - x49*x79 + x52*x79 - x55*x79 - x57*x79;
    res[6]=(*this)[13]*x103 + (*this)[15]*x189 - (*this)[2]*x127 + (*this)[3]*x125 - (*this)[4]*x101 + (*this)[4]*x152 + (*this)[4]*x176 - (*this)[5]*x124 - (*this)[5]*x187 + (*this)[5]*x196 - (*this)[6]*x140 + (*this)[6]*x171 - (*this)[6]*x175 + (*this)[7]*x193 + A[10]*x107 - A[10]*x91 + A[10]*x97 - A[10]*x98 + A[1]*x102 + A[1]*x105 - A[1]*x106 - A[1]*x65 - A[1]*x96 + A[2]*x20 + A[2]*x26 - A[3]*x136 + A[3]*x138 + A[3]*x142 - A[3]*x158 - A[4]*x24 + A[4]*x50 + A[5]*x134 + A[5]*x135 + A[5]*x191 + A[6]*x0 - A[6]*x1 - A[6]*x10 - A[6]*x11 + A[6]*x12 + A[6]*x13 - A[6]*x14 + A[6]*x15 - A[6]*x2 + A[6]*x3 - A[6]*x4 + A[6]*x5 - A[6]*x6 + A[6]*x7 - A[6]*x8 + A[6]*x9 - A[7]*x161 + A[7]*x163 - A[7]*x165 + A[7]*x166 + A[8]*x167 + A[8]*x173 - A[8]*x178 + A[8]*x183 + A[8]*x190 - A[9]*x22 + A[9]*x51 - x113*x23 - x114*x131 - x117*x129 - x126*x70 - x128*x67 - x128*x68 + x131*x82 + x139*x74 + x177*x30 + x177*x32 - x177*x36 + x177*x42 - x177*x46 + x177*x61 + x195*x28 - x195*x35 + x195*x40 + x195*x48 + x195*x54 + x195*x60 + x23*x86 - x44*x84 + x45*x84 - x49*x84 - x52*x84 + x55*x84 - x57*x84;
    res[7]=(*this)[14]*x103 + (*this)[14]*x95 - (*this)[15]*x100 - (*this)[15]*x188 + (*this)[15]*x196 + (*this)[1]*x93 + (*this)[2]*x169 - (*this)[3]*x152 - (*this)[3]*x176 + (*this)[7]*x171 + (*this)[7]*x172 - (*this)[8]*x104 - (*this)[9]*x123 - (*this)[9]*x64 + A[10]*x22 + A[10]*x51 + A[1]*x108 - A[1]*x109 - A[1]*x116 + A[1]*x199 - A[1]*x66 - A[1]*x71 + A[1]*x94 + A[2]*x19 - A[2]*x25 + A[3]*x24 + A[3]*x50 + A[5]*x162 - A[5]*x186 + A[6]*x161 + A[6]*x163 - A[6]*x164 + A[6]*x165 + A[6]*x166 + A[7]*x0 + A[7]*x1 + A[7]*x10 - A[7]*x11 + A[7]*x12 + A[7]*x13 + A[7]*x14 - A[7]*x15 - A[7]*x2 + A[7]*x3 + A[7]*x4 - A[7]*x5 - A[7]*x6 - A[7]*x7 - A[7]*x8 - A[7]*x9 - A[8]*x180 + A[8]*x184 + A[8]*x197 - A[8]*x92 + A[9]*x107 + A[9]*x198 + A[9]*x91 + A[9]*x97 + A[9]*x98 - x115*x126 + x126*x83 - x128*x201 - x131*x200 + x139*x78 - x153*x88 - x153*x89 + x168*x18 - x170*x73 - x170*x74 + x185*x30 - x185*x32 + x185*x36 + x185*x42 + x185*x46 + x185*x61 + x195*x34 - x195*x37 + x195*x38 - x195*x47 + x195*x59 + x195*x62 + x23*x67 - x44*x75 - x45*x75 - x49*x75 - x52*x75 - x55*x75 + x57*x75;
    res[8]=(*this)[15]*A[7]*x85 + (*this)[2]*x133 + (*this)[3]*x181 + (*this)[4]*x182 - (*this)[5]*x209 + (*this)[6]*x204 - (*this)[7]*x104 + (*this)[8]*x172 + (*this)[9]*x100 + A[10]*x162 + A[10]*x174 - A[10]*x186 + A[10]*x203 - A[10]*x207 - A[1]*x24 - A[1]*x50 - A[2]*x137 - A[3]*x109 + A[3]*x116 + A[3]*x199 - A[3]*x66 - A[3]*x71 - A[3]*x94 + A[4]*x106 + A[4]*x65 + A[5]*x22 + A[5]*x51 + A[6]*x167 + A[6]*x173 - A[6]*x178 - A[6]*x190 - A[7]*x120 - A[7]*x179 + A[7]*x180 - A[7]*x92 + A[8]*x0 - A[8]*x1 + A[8]*x10 - A[8]*x11 - A[8]*x12 - A[8]*x13 + A[8]*x14 + A[8]*x15 + A[8]*x2 + A[8]*x3 - A[8]*x4 - A[8]*x5 + A[8]*x6 - A[8]*x7 - A[8]*x8 + A[8]*x9 - A[9]*x135 + A[9]*x191 - A[9]*x205 - A[9]*x206 - x115*x87 + x117*x80 - x119*x76 - x141*x44 + x141*x45 + x141*x49 + x141*x52 - x141*x55 - x141*x57 + x17*x195 + x18*x73 - x192*x84 + x194*x84 + x195*x31 + x195*x39 - x195*x41 - x195*x53 + x195*x56 + x195*x58 + x202*x30 - x202*x32 - x202*x36 + x202*x42 + x202*x46 - x202*x61 - x208*x84 - x68*x72 + x76*x86 - x78*x85 - x80*x88 - x80*x89 + x82*x87 + x83*x87;
    res[9]=(*this)[14]*x101 - (*this)[15]*x104 + (*this)[3]*x133 + (*this)[3]*x95 - (*this)[5]*x204 - (*this)[6]*A[5]*x87 - (*this)[6]*x209 - (*this)[7]*x123 - (*this)[7]*x64 - (*this)[8]*x100 + (*this)[9]*x171 - A[10]*x161 + A[10]*x163 - A[10]*x164 + A[10]*x165 - A[10]*x166 + A[1]*x19 - A[1]*x25 - A[2]*x108 + A[2]*x109 - A[2]*x116 + A[2]*x66 - A[2]*x71 + A[2]*x94 + A[5]*x167 + A[5]*x178 + A[6]*x22 + A[6]*x51 + A[7]*x107 - A[7]*x112 - A[7]*x90 + A[7]*x91 - A[7]*x97 - A[7]*x98 + A[8]*x135 - A[8]*x191 + A[8]*x205 - A[8]*x206 + A[9]*x0 - A[9]*x1 + A[9]*x10 - A[9]*x11 - A[9]*x12 + A[9]*x13 - A[9]*x14 + A[9]*x15 + A[9]*x2 - A[9]*x3 + A[9]*x4 - A[9]*x5 - A[9]*x6 + A[9]*x7 + A[9]*x8 - A[9]*x9 - x114*x85 + x117*x76 + x119*x80 + x122*x44 - x122*x45 + x122*x49 - x122*x52 + x122*x55 - x122*x57 + x17*x185 - x18*x70 + x185*x31 - x185*x39 + x185*x41 + x185*x53 - x185*x56 + x185*x58 + x192*x79 + x194*x79 + x201*x80 - x202*x34 - x202*x37 + x202*x38 + x202*x47 - x202*x59 - x202*x62 - x208*x79 - x63*x73 + x67*x69 - x76*x89 + x77*x87 - x80*x86 + x81*x85 + x83*x85;
    res[10]=-(*this)[13]*x101 + (*this)[2]*A[6]*x80 + (*this)[4]*x103 - (*this)[5]*A[2]*x76 + (*this)[6]*A[2]*x80 - (*this)[6]*x110 + (*this)[6]*x64 + A[10]*x0 + A[10]*x1 - A[10]*x10 - A[10]*x11 - A[10]*x12 + A[10]*x13 + A[10]*x14 - A[10]*x15 + A[10]*x2 - A[10]*x3 - A[10]*x4 + A[10]*x5 - A[10]*x6 - A[10]*x7 + A[10]*x8 + A[10]*x9 - A[1]*x20 - A[1]*x26 - A[2]*x102 - A[2]*x105 - A[2]*x106 - A[2]*x65 - A[3]*x150 + A[3]*x155 + A[3]*x157 + A[5]*x197 + A[6]*x112 + A[6]*x198 - A[6]*x90 - A[6]*x97 + A[6]*x98 - A[7]*x22 + A[7]*x51 + A[8]*x162 - A[8]*x174 + A[8]*x186 - A[8]*x203 - A[8]*x207 + A[9]*x161 + A[9]*x163 + A[9]*x164 - A[9]*x165 - A[9]*x166 + x113*x69 - x114*x18 + x121*x44 + x121*x45 + x121*x49 - x121*x52 - x121*x55 + x121*x57 + x146*x87 + x17*x177 + x177*x31 - x177*x39 - x177*x41 + x177*x53 + x177*x56 - x177*x58 + x18*x81 + x18*x82 + x18*x83 - x200*x85 + x202*x28 - x202*x35 + x202*x40 + x202*x48 - x202*x54 - x202*x60 + x210*x84 + x211*x79 + x212*x84 - x213*x79 - x214*x79 - x215*x84 - x63*x78 - x67*x80 - x68*x80 + x70*x85 + x72*x88 - x74*x87;
    res[11]=-(*this)[0]*x233 + (*this)[15]*x224 + (*this)[2]*x221 + (*this)[3]*x226 + (*this)[3]*x229 - (*this)[3]*x234 + (*this)[4]*x220 + (*this)[4]*x227 + (*this)[5]*x228 + (*this)[7]*x219 + (*this)[7]*x230 - (*this)[8]*x235 - (*this)[9]*x237 + A[0]*x137 + A[11]*x0 + A[11]*x1 + A[11]*x10 + A[11]*x11 + A[11]*x12 + A[11]*x13 + A[11]*x14 + A[11]*x15 + A[11]*x2 + A[11]*x3 + A[11]*x4 + A[11]*x5 + A[11]*x6 + A[11]*x7 + A[11]*x8 + A[11]*x9 + A[12]*x97 - A[12]*x98 - A[13]*x179 - A[13]*x180 + A[14]*x178 - A[14]*x190 - A[15]*x136 - A[15]*x145 + x17*x216 - x216*x31 + x216*x39 + x216*x41 + x216*x53 + x216*x56 + x216*x58 + x217*x29 + x218*x33 + x222*x29 + x223*x29 + x225*x29 + x231*x27 + x232*x33 - x236*x29 - x238*x33 - x239*x27 - x240*x29;
    res[12]=(*this)[0]*x228 + (*this)[12]*x241 - (*this)[13]*x220 + (*this)[14]*x226 + (*this)[14]*x229 - (*this)[14]*x234 - (*this)[15]*x235 - (*this)[1]*x221 - (*this)[5]*x233 - (*this)[6]*x21 - (*this)[7]*x237 + (*this)[8]*x224 + (*this)[9]*x219 + (*this)[9]*x230 + A[0]*x24 - A[0]*x50 + A[12]*x0 - A[12]*x1 + A[12]*x10 + A[12]*x11 - A[12]*x12 + A[12]*x13 - A[12]*x14 - A[12]*x15 - A[12]*x2 - A[12]*x3 + A[12]*x4 + A[12]*x5 + A[12]*x6 - A[12]*x7 - A[12]*x8 + A[12]*x9 - A[13]*x135 + A[13]*x205 + A[14]*x162 - A[14]*x203 - A[15]*x143 - A[15]*x157 + x16*x217 + x16*x222 + x16*x223 - x16*x225 - x16*x236 + x16*x240 + x192*x33 - x194*x33 + x208*x33 - x211*x27 - x213*x27 + x214*x27 + x216*x30 + x216*x32 - x216*x36 - x216*x42 + x216*x46 + x216*x61;
    res[13]=(*this)[12]*x227 - (*this)[12]*x246 + (*this)[13]*x241 - (*this)[14]*x221 - (*this)[15]*x237 + (*this)[15]*x244 - (*this)[1]*x229 - (*this)[1]*x247 - (*this)[3]*x249 + (*this)[4]*x243 + (*this)[5]*x21 - (*this)[5]*x248 + (*this)[7]*x235 - (*this)[8]*x219 - (*this)[8]*x230 + (*this)[8]*x245 + A[0]*x19 + A[0]*x25 + A[12]*x135 - A[12]*x205 + A[13]*x0 - A[13]*x1 + A[13]*x10 + A[13]*x11 - A[13]*x12 - A[13]*x13 + A[13]*x14 - A[13]*x15 - A[13]*x2 + A[13]*x3 - A[13]*x4 + A[13]*x5 - A[13]*x6 + A[13]*x7 + A[13]*x8 - A[13]*x9 + A[14]*x163 + A[14]*x164 - A[15]*x102 + A[15]*x96 + x16*x218 - x16*x232 - x16*x238 - x192*x29 - x194*x29 + x208*x29 + x210*x27 - x212*x27 + x215*x27 + x216*x34 - x216*x37 - x216*x38 + x216*x47 - x216*x59 + x216*x62 + x242*x43;
    res[14]=(*this)[0]*x230 - (*this)[0]*x245 + (*this)[0]*x250 + (*this)[0]*x251 - (*this)[12]*x229 - (*this)[12]*x247 + (*this)[13]*x252 - (*this)[13]*x253 + (*this)[14]*x241 + (*this)[15]*x21 - (*this)[15]*x248 - (*this)[1]*x220 + (*this)[1]*x246 - (*this)[3]*x243 - (*this)[4]*x249 + (*this)[5]*x237 - (*this)[7]*x233 - A[0]*x20 + A[0]*x26 + A[12]*x162 + A[12]*x203 + A[13]*x163 - A[13]*x164 + A[14]*x0 + A[14]*x1 - A[14]*x10 + A[14]*x11 - A[14]*x12 - A[14]*x13 - A[14]*x14 + A[14]*x15 - A[14]*x2 + A[14]*x3 + A[14]*x4 - A[14]*x5 - A[14]*x6 - A[14]*x7 + A[14]*x8 + A[14]*x9 - A[15]*x94 - x130*x43 + x132*x43 - x16*x231 + x16*x239 - x210*x33 - x211*x29 - x212*x33 + x213*x29 + x214*x29 + x215*x33 + x216*x28 + x216*x35 - x216*x40 + x216*x48 - x216*x54 + x216*x60;
    res[15]=-(*this)[0]*(*this)[2]*x29 - (*this)[0]*x249 - (*this)[12]*(*this)[15]*x29 - (*this)[13]*(*this)[15]*x33 - (*this)[13]*x244 - (*this)[14]*(*this)[15]*x27 - (*this)[14]*x21 + (*this)[15]*x241 - (*this)[1]*(*this)[5]*x29 - (*this)[3]*(*this)[6]*x16 - (*this)[4]*x219 + (*this)[4]*x245 - (*this)[4]*x250 - (*this)[4]*x251 - (*this)[5]*x221 - (*this)[5]*x252 + (*this)[5]*x253 - (*this)[6]*x226 + (*this)[6]*x234 - (*this)[6]*x247 - (*this)[7]*x220 - (*this)[7]*x227 + (*this)[7]*x246 - (*this)[8]*x243 - A[0]*x22 - A[0]*x51 - A[12]*x143 + A[12]*x157 + A[13]*x102 + A[13]*x96 + A[14]*x94 + A[15]*x0 + A[15]*x1 + A[15]*x10 - A[15]*x11 + A[15]*x12 - A[15]*x13 - A[15]*x14 - A[15]*x15 - A[15]*x2 - A[15]*x3 - A[15]*x4 - A[15]*x5 + A[15]*x6 + A[15]*x7 + A[15]*x8 + A[15]*x9 - x130*x27 + x132*x27 - x216*x44 - x216*x45 + x216*x49 - x216*x52 - x216*x55 - x216*x57 - x242*x33;
    return res;
};

//-------------------------------------
// (R130MV, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8] + B[1];
    res[9]=A[9] + B[2];
    res[10]=A[10] + B[3];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    res[8]=A[8] - B[1];
    res[9]=A[9] - B[2];
    res[10]=A[10] - B[3];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[10]*B[3] - A[8]*B[1] - A[9]*B[2];
    res[1]=-A[12]*B[1] - A[13]*B[2] - A[14]*B[3] + A[1]*B[0];
    res[2]=A[11]*B[1] + A[2]*B[0] - A[3]*B[3] + A[4]*B[2];
    res[3]=A[11]*B[2] + A[2]*B[3] + A[3]*B[0] - A[4]*B[1];
    res[4]=A[11]*B[3] - A[2]*B[2] + A[3]*B[1] + A[4]*B[0];
    res[5]=-A[15]*B[1] + A[5]*B[0] - A[6]*B[3] + A[7]*B[2];
    res[6]=-A[15]*B[2] + A[5]*B[3] + A[6]*B[0] - A[7]*B[1];
    res[7]=-A[15]*B[3] - A[5]*B[2] + A[6]*B[1] + A[7]*B[0];
    res[8]=A[0]*B[1] + A[10]*B[2] + A[8]*B[0] - A[9]*B[3];
    res[9]=A[0]*B[2] - A[10]*B[1] + A[8]*B[3] + A[9]*B[0];
    res[10]=A[0]*B[3] + A[10]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[11]=A[11]*B[0] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[12]=A[12]*B[0] - A[13]*B[3] + A[14]*B[2] + A[1]*B[1];
    res[13]=A[12]*B[3] + A[13]*B[0] - A[14]*B[1] + A[1]*B[2];
    res[14]=-A[12]*B[2] + A[13]*B[1] + A[14]*B[0] + A[1]*B[3];
    res[15]=A[15]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[10]*B[3] - A[8]*B[1] - A[9]*B[2];
    res[1]=-A[12]*B[1] - A[13]*B[2] - A[14]*B[3] + A[1]*B[0];
    res[2]=A[11]*B[1] + A[2]*B[0];
    res[3]=A[11]*B[2] + A[3]*B[0];
    res[4]=A[11]*B[3] + A[4]*B[0];
    res[5]=-A[15]*B[1] + A[5]*B[0];
    res[6]=-A[15]*B[2] + A[6]*B[0];
    res[7]=-A[15]*B[3] + A[7]*B[0];
    res[8]=A[0]*B[1] + A[8]*B[0];
    res[9]=A[0]*B[2] + A[9]*B[0];
    res[10]=A[0]*B[3] + A[10]*B[0];
    res[11]=A[11]*B[0] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[12]=A[12]*B[0] + A[1]*B[1];
    res[13]=A[13]*B[0] + A[1]*B[2];
    res[14]=A[14]*B[0] + A[1]*B[3];
    res[15]=A[15]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const Rotation<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[3]*B[3] + A[4]*B[2];
    res[3]=A[2]*B[3] - A[4]*B[1];
    res[4]=-A[2]*B[2] + A[3]*B[1];
    res[5]=-A[6]*B[3] + A[7]*B[2];
    res[6]=A[5]*B[3] - A[7]*B[1];
    res[7]=-A[5]*B[2] + A[6]*B[1];
    res[8]=A[10]*B[2] - A[9]*B[3];
    res[9]=-A[10]*B[1] + A[8]*B[3];
    res[10]=-A[8]*B[2] + A[9]*B[1];
    res[11]=0;
    res[12]=-A[13]*B[3] + A[14]*B[2];
    res[13]=A[12]*B[3] - A[14]*B[1];
    res[14]=-A[12]*B[2] + A[13]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const Rotation<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[10], 2);
    T x2 = std::pow((*this)[11], 2);
    T x3 = std::pow((*this)[2], 2);
    T x4 = std::pow((*this)[3], 2);
    T x5 = std::pow((*this)[4], 2);
    T x6 = std::pow((*this)[8], 2);
    T x7 = std::pow((*this)[9], 2);
    T x8 = std::pow((*this)[12], 2);
    T x9 = std::pow((*this)[13], 2);
    T x10 = std::pow((*this)[14], 2);
    T x11 = std::pow((*this)[15], 2);
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[5], 2);
    T x14 = std::pow((*this)[6], 2);
    T x15 = std::pow((*this)[7], 2);
    T x16 = 2.0*(*this)[10];
    T x17 = (*this)[12]*A[2];
    T x18 = A[3]*x16;
    T x19 = 2.0*(*this)[8];
    T x20 = A[3]*x19;
    T x21 = 2.0*(*this)[9];
    T x22 = (*this)[14]*A[1];
    T x23 = A[1]*x19;
    T x24 = A[2]*x21;
    T x25 = 2.0*A[3];
    T x26 = (*this)[6]*x25;
    T x27 = 2.0*A[1];
    T x28 = (*this)[7]*x27;
    T x29 = 2.0*A[2];
    T x30 = (*this)[5]*x29;
    T x31 = (*this)[0]*x27;
    T x32 = (*this)[0]*x29;
    T x33 = (*this)[0]*x25;
    T x34 = A[1]*x16;
    T x35 = (*this)[5]*x27;
    T x36 = (*this)[6]*x29;
    T x37 = (*this)[7]*x25;
    T x38 = A[3]*x21;
    T x39 = A[2]*x19;
    T x40 = (*this)[15]*x27;
    T x41 = (*this)[15]*x29;
    T x42 = (*this)[15]*x25;
    T x43 = (*this)[7]*x29;
    T x44 = (*this)[5]*x25;
    T x45 = (*this)[6]*x27;
    T x46 = 2.0*x17;
    T x47 = A[2]*x16;
    T x48 = 2.0*x22;
    T x49 = A[1]*x21;
    T x50 = (*this)[13]*x25;
    T x51 = (*this)[11]*(*this)[1];
    T x52 = (*this)[12]*x27;
    T x53 = (*this)[12]*x25;
    T x54 = (*this)[13]*(*this)[2];
    T x55 = (*this)[14]*x25;
    T x56 = (*this)[1]*x25;
    T x57 = (*this)[14]*x29;
    T x58 = (*this)[13]*(*this)[3];
    T x59 = (*this)[1]*(*this)[4];
    T x60 = (*this)[4]*x29;
    T x61 = (*this)[2]*x29;
    T x62 = (*this)[11]*x27;
    T x63 = (*this)[3]*x27;
    T x64 = (*this)[3]*x25;
    T x65 = (*this)[2]*x25;
    T x66 = 2.0*A[0];
    T x67 = (*this)[0]*x66;
    T x68 = A[0]*x16;
    T x69 = (*this)[5]*x66;
    T x70 = (*this)[6]*x66;
    T x71 = (*this)[7]*x66;
    T x72 = A[0]*x19;
    T x73 = A[0]*x21;
    T x74 = (*this)[15]*x66;
    res[0]=A[0]*(x0 + x1 - x10 - x11 - x12 - x13 - x14 - x15 + x2 + x3 + x4 + x5 + x6 + x7 - x8 - x9);
    res[1]=-(*this)[11]*x35 - (*this)[11]*x36 - (*this)[11]*x37 - (*this)[12]*x31 - (*this)[12]*x38 + (*this)[13]*x20 - (*this)[13]*x32 - (*this)[13]*x34 - (*this)[14]*x33 - (*this)[14]*x39 + (*this)[1]*x18 + (*this)[1]*x23 + (*this)[1]*x24 + (*this)[2]*x26 - (*this)[2]*x40 - (*this)[2]*x43 + (*this)[3]*x28 - (*this)[3]*x41 - (*this)[3]*x44 + (*this)[4]*x30 - (*this)[4]*x42 - (*this)[4]*x45 + x16*x17 + x21*x22;
    res[2]=(*this)[11]*x31 + (*this)[11]*x38 - (*this)[11]*x47 + (*this)[12]*x35 + (*this)[12]*x37 + (*this)[13]*x30 - (*this)[13]*x42 - (*this)[13]*x45 + (*this)[14]*x41 + (*this)[14]*x44 + (*this)[1]*x26 - (*this)[1]*x40 - (*this)[1]*x43 + (*this)[2]*x18 + (*this)[2]*x23 + (*this)[2]*x24 + (*this)[3]*x33 + (*this)[3]*x39 - (*this)[3]*x49 + (*this)[4]*x20 - (*this)[4]*x32 - (*this)[4]*x34 + (*this)[6]*x46 - (*this)[7]*x48;
    res[3]=-(*this)[11]*x20 + (*this)[11]*x32 + (*this)[11]*x34 + (*this)[12]*x42 + (*this)[12]*x45 + (*this)[13]*x35 + (*this)[13]*x36 + (*this)[13]*x37 + (*this)[14]*x26 - (*this)[14]*x43 - (*this)[15]*x48 + (*this)[1]*x28 - (*this)[1]*x41 - (*this)[1]*x44 - (*this)[2]*x33 - (*this)[2]*x39 + (*this)[2]*x49 + (*this)[3]*x18 + (*this)[3]*x23 + (*this)[3]*x24 + (*this)[4]*x31 + (*this)[4]*x38 - (*this)[4]*x47 - (*this)[5]*x46;
    res[4]=(*this)[11]*x33 + (*this)[11]*x39 - (*this)[11]*x49 + (*this)[12]*x28 - (*this)[12]*x44 - (*this)[13]*x26 + (*this)[13]*x40 + (*this)[13]*x43 + (*this)[14]*x36 + (*this)[14]*x37 - (*this)[15]*x46 + (*this)[1]*x30 - (*this)[1]*x42 - (*this)[1]*x45 - (*this)[2]*x20 + (*this)[2]*x32 + (*this)[2]*x34 - (*this)[3]*x31 - (*this)[3]*x38 + (*this)[3]*x47 + (*this)[4]*x18 + (*this)[4]*x23 + (*this)[4]*x24 + (*this)[5]*x48;
    res[5]=(*this)[0]*x26 + (*this)[11]*x50 - (*this)[11]*x57 - (*this)[15]*x31 - (*this)[15]*x38 + (*this)[15]*x47 + (*this)[2]*x52 + (*this)[2]*x55 + (*this)[3]*x46 + (*this)[3]*x56 - (*this)[4]*x48 + (*this)[4]*x53 + (*this)[5]*x18 + (*this)[5]*x23 + (*this)[5]*x24 + (*this)[6]*x39 - (*this)[6]*x49 + (*this)[7]*x20 - (*this)[7]*x32 - (*this)[7]*x34 + x27*x51 - x27*x58 + x29*x54 - x29*x59;
    res[6]=(*this)[0]*x28 + (*this)[11]*x48 - (*this)[11]*x53 + (*this)[15]*x20 - (*this)[15]*x32 - (*this)[15]*x34 - (*this)[2]*x46 - (*this)[2]*x56 + (*this)[3]*x52 + (*this)[3]*x55 + (*this)[4]*x50 - (*this)[4]*x57 - (*this)[5]*x33 - (*this)[5]*x39 + (*this)[5]*x49 + (*this)[6]*x18 + (*this)[6]*x23 + (*this)[6]*x24 + (*this)[7]*x38 - (*this)[7]*x47 + x27*x54 + x27*x59 + x29*x51 + x29*x58;
    res[7]=(*this)[0]*x30 + (*this)[11]*x46 + (*this)[13]*x60 - (*this)[13]*x62 - (*this)[15]*x33 - (*this)[15]*x39 + (*this)[15]*x49 + (*this)[1]*x61 - (*this)[1]*x63 + (*this)[2]*x48 - (*this)[2]*x53 - (*this)[3]*x50 + (*this)[3]*x57 + (*this)[4]*x52 + (*this)[4]*x55 - (*this)[5]*x20 + (*this)[5]*x34 - (*this)[6]*x31 - (*this)[6]*x38 + (*this)[6]*x47 + (*this)[7]*x18 + (*this)[7]*x23 + (*this)[7]*x24 + x25*x51;
    res[8]=(*this)[0]*x38 - (*this)[0]*x47 - (*this)[11]*x60 + (*this)[11]*x64 - (*this)[13]*x46 - (*this)[14]*x53 + (*this)[15]*x26 - (*this)[1]*x50 + (*this)[1]*x57 + (*this)[3]*x61 + (*this)[4]*x65 - (*this)[5]*x37 - (*this)[6]*x30 - (*this)[7]*x41 + (*this)[8]*x18 + (*this)[9]*x39 + A[1]*x0 - A[1]*x1 + A[1]*x10 - A[1]*x11 - A[1]*x12 - A[1]*x13 + A[1]*x14 + A[1]*x15 + A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 - A[1]*x8 + A[1]*x9;
    res[9]=-(*this)[0]*x20 + (*this)[0]*x34 - (*this)[11]*x65 - (*this)[13]*x52 - (*this)[14]*x50 + (*this)[15]*x28 - (*this)[1]*x48 + (*this)[1]*x53 + (*this)[2]*x63 + (*this)[4]*x62 + (*this)[4]*x64 - (*this)[5]*x42 - (*this)[6]*x35 - (*this)[7]*x26 + (*this)[9]*x18 + (*this)[9]*x23 + A[2]*x0 - A[2]*x1 + A[2]*x10 - A[2]*x11 - A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + A[2]*x8 - A[2]*x9;
    res[10]=(*this)[0]*x39 - (*this)[0]*x49 + (*this)[11]*x61 - (*this)[12]*x48 + (*this)[13]*(*this)[1]*x27 - (*this)[13]*x57 + (*this)[15]*x30 - (*this)[1]*x46 + (*this)[2]*(*this)[4]*x27 + (*this)[3]*x60 - (*this)[3]*x62 - (*this)[5]*x28 - (*this)[6]*x40 - (*this)[7]*x36 + (*this)[8]*x34 + (*this)[9]*x47 + A[3]*x0 + A[3]*x1 - A[3]*x10 - A[3]*x11 - A[3]*x12 + A[3]*x13 + A[3]*x14 - A[3]*x15 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7 + A[3]*x8 + A[3]*x9;
    res[11]=(*this)[11]*x67 + (*this)[12]*x69 + (*this)[13]*x70 + (*this)[14]*x71 - (*this)[1]*x74 + (*this)[2]*x72 + (*this)[3]*x73 + (*this)[4]*x68;
    res[12]=(*this)[11]*x69 + (*this)[12]*x67 - (*this)[13]*x68 + (*this)[14]*x73 - (*this)[1]*x72 + (*this)[2]*x74 + (*this)[3]*x71 - (*this)[4]*x70;
    res[13]=(*this)[11]*x70 + (*this)[12]*x68 + (*this)[13]*x67 - (*this)[14]*x72 - (*this)[1]*x73 - (*this)[2]*x71 + (*this)[3]*x74 + (*this)[4]*x69;
    res[14]=(*this)[11]*x71 - (*this)[12]*x73 + (*this)[13]*x72 + (*this)[14]*x67 - (*this)[1]*x68 + (*this)[2]*x70 - (*this)[3]*x69 + (*this)[4]*x74;
    res[15]=-(*this)[12]*(*this)[2]*x66 - (*this)[14]*(*this)[4]*x66 + (*this)[15]*x67 - (*this)[5]*x72 - (*this)[6]*x73 - (*this)[7]*x68 - x51*x66 - x58*x66;
    return res;
};

//----------------------------------
// (R130MV, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const R130MV<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] + B[1];
    res[6]=A[6] + B[2];
    res[7]=A[7] + B[3];
    res[8]=A[8] + B[4];
    res[9]=A[9] + B[5];
    res[10]=A[10] + B[6];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15] + B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const R130MV<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5] - B[1];
    res[6]=A[6] - B[2];
    res[7]=A[7] - B[3];
    res[8]=A[8] - B[4];
    res[9]=A[9] - B[5];
    res[10]=A[10] - B[6];
    res[11]=A[11];
    res[12]=A[12];
    res[13]=A[13];
    res[14]=A[14];
    res[15]=A[15] - B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[10]*B[6] - A[15]*B[7] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3] - A[8]*B[4] - A[9]*B[5];
    res[1]=-A[11]*B[7] - A[12]*B[4] - A[13]*B[5] - A[14]*B[6] + A[1]*B[0] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[2]=A[11]*B[4] + A[12]*B[7] + A[13]*B[3] - A[14]*B[2] - A[1]*B[1] + A[2]*B[0] - A[3]*B[6] + A[4]*B[5];
    res[3]=A[11]*B[5] - A[12]*B[3] + A[13]*B[7] + A[14]*B[1] - A[1]*B[2] + A[2]*B[6] + A[3]*B[0] - A[4]*B[4];
    res[4]=A[11]*B[6] + A[12]*B[2] - A[13]*B[1] + A[14]*B[7] - A[1]*B[3] - A[2]*B[5] + A[3]*B[4] + A[4]*B[0];
    res[5]=A[0]*B[1] + A[10]*B[2] - A[15]*B[4] + A[5]*B[0] - A[6]*B[6] + A[7]*B[5] - A[8]*B[7] - A[9]*B[3];
    res[6]=A[0]*B[2] - A[10]*B[1] - A[15]*B[5] + A[5]*B[6] + A[6]*B[0] - A[7]*B[4] + A[8]*B[3] - A[9]*B[7];
    res[7]=A[0]*B[3] - A[10]*B[7] - A[15]*B[6] - A[5]*B[5] + A[6]*B[4] + A[7]*B[0] - A[8]*B[2] + A[9]*B[1];
    res[8]=A[0]*B[4] + A[10]*B[5] + A[15]*B[1] + A[5]*B[7] + A[6]*B[3] - A[7]*B[2] + A[8]*B[0] - A[9]*B[6];
    res[9]=A[0]*B[5] - A[10]*B[4] + A[15]*B[2] - A[5]*B[3] + A[6]*B[7] + A[7]*B[1] + A[8]*B[6] + A[9]*B[0];
    res[10]=A[0]*B[6] + A[10]*B[0] + A[15]*B[3] + A[5]*B[2] - A[6]*B[1] + A[7]*B[7] - A[8]*B[5] + A[9]*B[4];
    res[11]=A[11]*B[0] + A[12]*B[1] + A[13]*B[2] + A[14]*B[3] + A[1]*B[7] - A[2]*B[4] - A[3]*B[5] - A[4]*B[6];
    res[12]=A[11]*B[1] + A[12]*B[0] - A[13]*B[6] + A[14]*B[5] + A[1]*B[4] - A[2]*B[7] - A[3]*B[3] + A[4]*B[2];
    res[13]=A[11]*B[2] + A[12]*B[6] + A[13]*B[0] - A[14]*B[4] + A[1]*B[5] + A[2]*B[3] - A[3]*B[7] - A[4]*B[1];
    res[14]=A[11]*B[3] - A[12]*B[5] + A[13]*B[4] + A[14]*B[0] + A[1]*B[6] - A[2]*B[2] + A[3]*B[1] - A[4]*B[7];
    res[15]=A[0]*B[7] + A[10]*B[3] + A[15]*B[0] + A[5]*B[4] + A[6]*B[5] + A[7]*B[6] + A[8]*B[1] + A[9]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const R130MV<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[10]*B[6] - A[15]*B[7] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3] - A[8]*B[4] - A[9]*B[5];
    res[1]=-A[12]*B[4] - A[13]*B[5] - A[14]*B[6] + A[1]*B[0];
    res[2]=A[11]*B[4] + A[13]*B[3] - A[14]*B[2] + A[2]*B[0];
    res[3]=A[11]*B[5] - A[12]*B[3] + A[14]*B[1] + A[3]*B[0];
    res[4]=A[11]*B[6] + A[12]*B[2] - A[13]*B[1] + A[4]*B[0];
    res[5]=A[0]*B[1] - A[15]*B[4] + A[5]*B[0] - A[8]*B[7];
    res[6]=A[0]*B[2] - A[15]*B[5] + A[6]*B[0] - A[9]*B[7];
    res[7]=A[0]*B[3] - A[10]*B[7] - A[15]*B[6] + A[7]*B[0];
    res[8]=A[0]*B[4] + A[15]*B[1] + A[5]*B[7] + A[8]*B[0];
    res[9]=A[0]*B[5] + A[15]*B[2] + A[6]*B[7] + A[9]*B[0];
    res[10]=A[0]*B[6] + A[10]*B[0] + A[15]*B[3] + A[7]*B[7];
    res[11]=A[11]*B[0] - A[2]*B[4] - A[3]*B[5] - A[4]*B[6];
    res[12]=A[12]*B[0] + A[1]*B[4] - A[3]*B[3] + A[4]*B[2];
    res[13]=A[13]*B[0] + A[1]*B[5] + A[2]*B[3] - A[4]*B[1];
    res[14]=A[14]*B[0] + A[1]*B[6] - A[2]*B[2] + A[3]*B[1];
    res[15]=A[0]*B[7] + A[10]*B[3] + A[15]*B[0] + A[5]*B[4] + A[6]*B[5] + A[7]*B[6] + A[8]*B[1] + A[9]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const R130MV<T> &A, const Rotor<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[11]*B[7] - A[2]*B[1] - A[3]*B[2] - A[4]*B[3];
    res[2]=A[12]*B[7] - A[1]*B[1] - A[3]*B[6] + A[4]*B[5];
    res[3]=A[13]*B[7] - A[1]*B[2] + A[2]*B[6] - A[4]*B[4];
    res[4]=A[14]*B[7] - A[1]*B[3] - A[2]*B[5] + A[3]*B[4];
    res[5]=A[10]*B[2] - A[6]*B[6] + A[7]*B[5] - A[9]*B[3];
    res[6]=-A[10]*B[1] + A[5]*B[6] - A[7]*B[4] + A[8]*B[3];
    res[7]=-A[5]*B[5] + A[6]*B[4] - A[8]*B[2] + A[9]*B[1];
    res[8]=A[10]*B[5] + A[6]*B[3] - A[7]*B[2] - A[9]*B[6];
    res[9]=-A[10]*B[4] - A[5]*B[3] + A[7]*B[1] + A[8]*B[6];
    res[10]=A[5]*B[2] - A[6]*B[1] - A[8]*B[5] + A[9]*B[4];
    res[11]=A[12]*B[1] + A[13]*B[2] + A[14]*B[3] + A[1]*B[7];
    res[12]=A[11]*B[1] - A[13]*B[6] + A[14]*B[5] - A[2]*B[7];
    res[13]=A[11]*B[2] + A[12]*B[6] - A[14]*B[4] - A[3]*B[7];
    res[14]=A[11]*B[3] - A[12]*B[5] + A[13]*B[4] - A[4]*B[7];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> R130MV<T>::conjugate(const Rotor<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[10], 2);
    T x2 = std::pow((*this)[11], 2);
    T x3 = std::pow((*this)[2], 2);
    T x4 = std::pow((*this)[3], 2);
    T x5 = std::pow((*this)[4], 2);
    T x6 = std::pow((*this)[8], 2);
    T x7 = std::pow((*this)[9], 2);
    T x8 = std::pow((*this)[12], 2);
    T x9 = std::pow((*this)[13], 2);
    T x10 = std::pow((*this)[14], 2);
    T x11 = std::pow((*this)[15], 2);
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[5], 2);
    T x14 = std::pow((*this)[6], 2);
    T x15 = std::pow((*this)[7], 2);
    T x16 = 2.0*A[7];
    T x17 = (*this)[10]*(*this)[7];
    T x18 = (*this)[5]*(*this)[8];
    T x19 = (*this)[6]*(*this)[9];
    T x20 = (*this)[0]*(*this)[15];
    T x21 = (*this)[11]*(*this)[1];
    T x22 = (*this)[12]*(*this)[2];
    T x23 = (*this)[13]*(*this)[3];
    T x24 = (*this)[14]*(*this)[4];
    T x25 = 2.0*(*this)[0];
    T x26 = A[1]*x25;
    T x27 = (*this)[3]*x25;
    T x28 = A[3]*x25;
    T x29 = 2.0*(*this)[10];
    T x30 = (*this)[12]*x29;
    T x31 = A[6]*x29;
    T x32 = (*this)[3]*A[1];
    T x33 = 2.0*(*this)[6];
    T x34 = (*this)[12]*x33;
    T x35 = 2.0*(*this)[13];
    T x36 = (*this)[7]*A[1];
    T x37 = (*this)[8]*x35;
    T x38 = 2.0*(*this)[14];
    T x39 = (*this)[5]*A[2];
    T x40 = A[4]*x38;
    T x41 = 2.0*(*this)[1];
    T x42 = (*this)[8]*A[4];
    T x43 = (*this)[9]*x41;
    T x44 = A[6]*x33;
    T x45 = 2.0*(*this)[9];
    T x46 = (*this)[2]*A[3];
    T x47 = 2.0*A[4];
    T x48 = (*this)[3]*(*this)[7];
    T x49 = 2.0*(*this)[4];
    T x50 = (*this)[5]*A[5];
    T x51 = (*this)[8]*x49;
    T x52 = A[4]*x25;
    T x53 = A[5]*x25;
    T x54 = (*this)[14]*x25;
    T x55 = (*this)[11]*x29;
    T x56 = A[4]*x29;
    T x57 = A[2]*x29;
    T x58 = (*this)[5]*x47;
    T x59 = (*this)[11]*x33;
    T x60 = 2.0*(*this)[11];
    T x61 = (*this)[7]*A[6];
    T x62 = A[1]*x60;
    T x63 = (*this)[11]*x45;
    T x64 = 2.0*(*this)[12];
    T x65 = (*this)[15]*A[1];
    T x66 = (*this)[7]*A[2];
    T x67 = (*this)[12]*x45;
    T x68 = (*this)[15]*x35;
    T x69 = (*this)[5]*A[3];
    T x70 = (*this)[15]*x38;
    T x71 = (*this)[14]*A[1];
    T x72 = (*this)[8]*x38;
    T x73 = (*this)[15]*(*this)[2];
    T x74 = 2.0*A[5];
    T x75 = (*this)[15]*(*this)[3];
    T x76 = (*this)[15]*x49;
    T x77 = (*this)[5]*A[1];
    T x78 = A[2]*x33;
    T x79 = (*this)[7]*x41;
    T x80 = (*this)[2]*x74;
    T x81 = 2.0*(*this)[3];
    T x82 = (*this)[5]*A[6];
    T x83 = (*this)[8]*A[3];
    T x84 = A[4]*x33;
    T x85 = (*this)[4]*x45;
    T x86 = (*this)[8]*x64;
    T x87 = 2.0*A[3];
    T x88 = 2.0*(*this)[2];
    T x89 = (*this)[2]*x45;
    T x90 = (*this)[8]*x74;
    T x91 = (*this)[9]*x35;
    T x92 = (*this)[15]*x41;
    T x93 = 2.0*x46;
    T x94 = (*this)[3]*x45;
    T x95 = A[2]*x25;
    T x96 = A[3]*x29;
    T x97 = (*this)[1]*A[1];
    T x98 = (*this)[15]*x60;
    T x99 = (*this)[15]*x64;
    T x100 = A[4]*x35;
    T x101 = A[5]*x33;
    T x102 = (*this)[9]*x38;
    T x103 = A[6]*x25;
    T x104 = A[5]*x29;
    T x105 = (*this)[8]*x60;
    T x106 = (*this)[7]*A[5];
    T x107 = (*this)[2]*A[1];
    T x108 = 2.0*x32;
    T x109 = A[3]*x33;
    T x110 = (*this)[12]*x47;
    T x111 = 2.0*A[2];
    T x112 = (*this)[8]*A[2];
    T x113 = (*this)[8]*A[6];
    T x114 = A[3]*x49;
    T x115 = A[6]*x35;
    T x116 = (*this)[11]*x49;
    T x117 = (*this)[12]*x35;
    T x118 = (*this)[12]*A[3];
    T x119 = (*this)[3]*A[5];
    T x120 = (*this)[12]*A[6];
    T x121 = A[3]*x35;
    T x122 = A[5]*x35;
    T x123 = A[6]*x38;
    T x124 = A[6]*x41;
    T x125 = 2.0*x61;
    T x126 = A[5]*x38;
    T x127 = (*this)[3]*x60;
    T x128 = A[2]*x38;
    T x129 = 2.0*(*this)[15];
    T x130 = (*this)[15]*x45;
    T x131 = (*this)[4]*x41;
    T x132 = (*this)[2]*A[2];
    T x133 = 2.0*(*this)[7];
    T x134 = (*this)[8]*A[1];
    T x135 = 2.0*A[6];
    T x136 = (*this)[12]*x60;
    T x137 = A[1]*x38;
    T x138 = (*this)[12]*x41;
    T x139 = A[2]*x35;
    T x140 = (*this)[2]*A[5];
    T x141 = (*this)[3]*x41;
    T x142 = (*this)[3]*A[2];
    T x143 = 2.0*x36;
    T x144 = 2.0*A[1];
    T x145 = A[6]*x49;
    T x146 = (*this)[15]*(*this)[7];
    T x147 = (*this)[4]*x47;
    T x148 = (*this)[3]*x47;
    T x149 = A[0]*x25;
    T x150 = A[0]*x29;
    T x151 = (*this)[5]*A[0];
    T x152 = A[0]*x33;
    T x153 = (*this)[7]*A[0];
    T x154 = (*this)[2]*x16;
    T x155 = (*this)[8]*A[0];
    T x156 = (*this)[6]*x16;
    T x157 = (*this)[4]*x16;
    T x158 = (*this)[1]*x16;
    T x159 = (*this)[10]*x16;
    T x160 = (*this)[11]*x16;
    T x161 = (*this)[12]*x16;
    T x162 = (*this)[13]*x16;
    T x163 = 2.0*A[0];
    T x164 = (*this)[3]*x16;
    T x165 = (*this)[14]*x16;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - A[0]*x9 + x16*x17 + x16*x18 + x16*x19 - x16*x20 - x16*x21 - x16*x22 - x16*x23 - x16*x24;
    res[1]=-(*this)[11]*x58 - (*this)[12]*x52 - (*this)[13]*x53 - (*this)[13]*x56 + (*this)[1]*x31 - (*this)[1]*x78 + (*this)[2]*x26 + (*this)[2]*x44 - (*this)[2]*x57 + (*this)[4]*x28 - (*this)[4]*x84 - (*this)[7]*x80 - (*this)[8]*x62 + (*this)[9]*x40 - A[1]*x85 + A[2]*x27 + A[2]*x51 - A[2]*x63 - A[2]*x68 + A[3]*x34 - A[3]*x55 - A[3]*x70 - A[3]*x79 + A[5]*x30 + A[5]*x43 - A[5]*x59 - A[5]*x72 + A[6]*x37 - A[6]*x54 - A[6]*x67 - A[6]*x76 + x29*x32 - x33*x71 + x35*x36 - x35*x69 + x38*x39 + x41*x42 - x41*x77 + x45*x46 + x47*x48 - x47*x73 + x49*x50 - x60*x61 - x64*x65 - x64*x66 - x74*x75 - x81*x82 - x81*x83;
    res[2]=(*this)[11]*x52 + (*this)[12]*x58 + (*this)[13]*x28 - (*this)[13]*x84 + (*this)[15]*x62 + (*this)[1]*x26 + (*this)[1]*x44 - (*this)[1]*x57 + (*this)[2]*x31 - (*this)[2]*x78 + (*this)[3]*x90 - (*this)[4]*x53 - (*this)[4]*x56 - (*this)[7]*x40 - (*this)[7]*x93 + A[1]*x86 - A[1]*x91 + A[2]*x37 - A[2]*x54 + A[2]*x67 - A[2]*x76 + A[3]*x30 + A[3]*x43 - A[3]*x59 + A[3]*x72 - A[4]*x92 - A[4]*x94 + A[5]*x34 - A[5]*x55 + A[5]*x70 - A[5]*x79 + A[5]*x89 + A[6]*x27 + A[6]*x51 + A[6]*x63 - A[6]*x68 - x29*x71 + x32*x33 + x35*x50 + x36*x49 + x38*x82 - x39*x81 + x42*x88 - x49*x69 + x60*x66 + x61*x64 + x75*x87 - x77*x88;
    res[3]=(*this)[11]*x53 - (*this)[12]*x28 + (*this)[13]*x101 + (*this)[13]*x96 + (*this)[14]*x26 + (*this)[14]*x44 - (*this)[14]*x57 - (*this)[15]*x40 - (*this)[15]*x93 + (*this)[1]*x95 - (*this)[2]*x103 + (*this)[3]*x31 - (*this)[3]*x78 - (*this)[4]*x104 - (*this)[4]*x109 + (*this)[4]*x52 + (*this)[5]*x100 - (*this)[5]*x108 - (*this)[8]*x80 + A[1]*x37 + A[1]*x67 - A[2]*x86 + A[2]*x91 + A[2]*x98 + A[3]*x102 + A[4]*x34 + A[4]*x55 + A[4]*x79 + A[4]*x89 - A[5]*x92 + A[5]*x94 - A[6]*x105 + A[6]*x85 + A[6]*x99 - x106*x38 - x107*x33 + x29*x97 + x35*x61 - x36*x60 + x39*x88 - x41*x82 - x41*x83 + x42*x81 - x48*x87 + x49*x65 + x49*x66 - x50*x64 + x60*x69;
    res[4]=(*this)[11]*x103 + (*this)[12]*x95 - (*this)[13]*x26 - (*this)[13]*x44 + (*this)[13]*x57 + (*this)[14]*x101 + (*this)[14]*x96 - (*this)[15]*x108 + (*this)[1]*x28 - (*this)[1]*x84 + (*this)[2]*x53 + (*this)[2]*x56 + (*this)[3]*x104 + (*this)[3]*x109 + (*this)[4]*x31 - (*this)[4]*x78 + (*this)[5]*x40 + (*this)[5]*x93 + (*this)[7]*x110 - (*this)[7]*x114 + A[1]*x30 - A[1]*x43 + A[1]*x59 + A[1]*x72 + A[2]*x102 - A[3]*x91 + A[3]*x98 - A[4]*x27 - A[4]*x63 + A[4]*x68 + A[5]*x105 + A[5]*x85 - A[5]*x99 - A[6]*x92 - A[6]*x94 + x106*x35 - x111*x48 + x111*x73 + x112*x41 - x113*x88 - x36*x88 + x38*x61 - x39*x60 + x41*x50 + x42*x49 - x49*x77 - x64*x82 - x64*x83;
    res[5]=-(*this)[10]*x95 + (*this)[11]*x115 - (*this)[11]*x126 + (*this)[15]*x104 + (*this)[15]*x109 + (*this)[1]*x121 - (*this)[1]*x128 + (*this)[2]*x122 + (*this)[2]*x123 + (*this)[3]*x124 + (*this)[5]*x31 + (*this)[6]*x103 - (*this)[7]*x53 + (*this)[8]*x101 + (*this)[8]*x125 + (*this)[9]*x28 + A[1]*x0 - A[1]*x1 - A[1]*x10 - A[1]*x11 + A[1]*x12 - A[1]*x13 + A[1]*x14 + A[1]*x15 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 + A[1]*x8 - A[1]*x9 + A[2]*x116 + A[2]*x117 - A[3]*x127 - A[5]*x131 - A[6]*x130 + x112*x45 + x118*x38 + x119*x64 + x120*x49 - x129*x66 - x132*x81 - x133*x69 - x17*x47 + x18*x47 - x19*x47 - x20*x47 + x21*x47 + x22*x47 - x23*x47 - x24*x47 + x29*x83 - x33*x39 + x45*x50 - x46*x49;
    res[6]=(*this)[10]*x26 + (*this)[11]*x40 + (*this)[14]*x121 - (*this)[15]*x56 + (*this)[2]*x100 - (*this)[2]*x124 + (*this)[3]*x110 - (*this)[3]*x114 + (*this)[3]*x123 + (*this)[4]*x115 + (*this)[5]*A[4]*x45 + (*this)[6]*x31 - (*this)[7]*x109 + (*this)[7]*x52 - (*this)[8]*x28 + (*this)[9]*x96 - A[1]*x116 + A[1]*x117 + A[2]*x0 - A[2]*x1 - A[2]*x10 - A[2]*x11 + A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 - A[2]*x8 + A[2]*x9 + A[4]*x131 + x113*x129 - x118*x41 - x120*x60 + x129*x36 - x129*x69 + x134*x45 - x17*x74 - x18*x74 + x19*x74 - x20*x74 + x21*x74 - x22*x74 + x23*x74 - x24*x74 - x25*x82 - x32*x88 + x33*x42 - x33*x77 + x38*x97 + x45*x61 + x46*x60;
    res[7]=-(*this)[11]*x100 + (*this)[12]*x137 + (*this)[14]*x139 - (*this)[15]*x90 + (*this)[2]*x40 + (*this)[4]*x110 + (*this)[4]*x122 - (*this)[5]*x143 + (*this)[5]*x56 + (*this)[6]*x104 - (*this)[6]*x52 + (*this)[8]*x95 - (*this)[9]*x26 + (*this)[9]*x57 + A[2]*x138 + A[3]*x0 + A[3]*x1 + A[3]*x10 - A[3]*x11 + A[3]*x12 + A[3]*x13 + A[3]*x14 - A[3]*x15 - A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 - A[3]*x6 - A[3]*x7 - A[3]*x8 - A[3]*x9 + A[4]*x130 - A[4]*x141 + A[5]*x136 + x106*x45 - x107*x49 + x119*x38 + x129*x39 - x132*x60 + x133*x42 + x134*x29 + x135*x17 - x135*x18 - x135*x19 - x135*x20 + x135*x21 - x135*x22 - x135*x23 + x135*x24 + x140*x41 - x142*x49 + x25*x50 + x32*x60 - x33*x65 - x33*x66 - x35*x97;
    res[8]=-(*this)[10]*x53 + (*this)[11]*x121 - (*this)[11]*x128 + (*this)[12]*x114 + (*this)[15]*x44 - (*this)[15]*x57 - (*this)[1]*x115 + (*this)[1]*x126 + (*this)[2]*x145 + (*this)[3]*x80 - (*this)[5]*x125 - (*this)[6]*x28 + (*this)[8]*A[5]*x45 + (*this)[8]*x31 - (*this)[8]*x78 + (*this)[9]*x103 - A[2]*x131 + A[3]*x130 + A[3]*x141 + A[4]*x0 - A[4]*x1 + A[4]*x10 - A[4]*x11 - A[4]*x12 - A[4]*x13 + A[4]*x14 + A[4]*x15 + A[4]*x2 + A[4]*x3 - A[4]*x4 - A[4]*x5 + A[4]*x6 - A[4]*x7 - A[4]*x8 + A[4]*x9 - A[5]*x116 - A[5]*x117 + A[6]*x127 - x120*x38 + x132*x35 - x133*x83 + x142*x64 + x144*x17 - x144*x18 + x144*x19 + x144*x20 + x144*x21 + x144*x22 - x144*x23 - x144*x24 - x146*x74 + x25*x66 - x29*x69 - x33*x50 + x38*x46 - x39*x45;
    res[9]=(*this)[10]*x52 + (*this)[11]*x137 + (*this)[11]*x147 - (*this)[12]*x100 - (*this)[14]*x115 - (*this)[1]*x40 - (*this)[2]*A[6]*x60 + (*this)[2]*x148 + (*this)[3]*A[3]*x38 + (*this)[3]*x145 + (*this)[4]*x121 + (*this)[5]*x28 - (*this)[5]*x84 - (*this)[6]*x96 - (*this)[7]*A[3]*x45 - (*this)[7]*x26 - (*this)[7]*x44 - (*this)[8]*x103 + (*this)[9]*x31 + A[1]*x131 + A[5]*x0 - A[5]*x1 + A[5]*x10 - A[5]*x11 - A[5]*x12 + A[5]*x13 - A[5]*x14 + A[5]*x15 + A[5]*x2 - A[5]*x3 + A[5]*x4 - A[5]*x5 - A[5]*x6 + A[5]*x7 + A[5]*x8 - A[5]*x9 + x107*x35 + x111*x17 + x111*x18 - x111*x19 + x111*x20 + x111*x21 - x111*x22 + x111*x23 - x111*x24 - x118*x60 + x120*x41 - x129*x82 - x129*x83 - x134*x33 + x146*x47 + x29*x65 + x32*x64 - x41*x46 + x42*x45 - x45*x77;
    res[10]=-(*this)[11]*A[1]*x35 - (*this)[11]*x148 + (*this)[12]*A[1]*x49 - (*this)[12]*x40 - (*this)[14]*x122 + (*this)[15]*(*this)[8]*x111 - (*this)[15]*x84 + (*this)[1]*x100 + (*this)[2]*x147 + (*this)[3]*x128 + (*this)[4]*x139 + (*this)[6]*x26 - (*this)[6]*x57 - (*this)[7]*x101 - (*this)[7]*x58 - (*this)[8]*x143 + (*this)[8]*x53 + (*this)[9]*x104 - (*this)[9]*x52 + A[2]*x136 - A[5]*x138 + A[6]*x0 + A[6]*x1 - A[6]*x10 - A[6]*x11 - A[6]*x12 + A[6]*x13 + A[6]*x14 - A[6]*x15 + A[6]*x2 - A[6]*x3 - A[6]*x4 + A[6]*x5 - A[6]*x6 - A[6]*x7 + A[6]*x8 + A[6]*x9 + x107*x38 + x119*x49 + x129*x50 + x132*x41 + x140*x60 - x17*x87 + x18*x87 + x19*x87 + x20*x87 + x21*x87 - x22*x87 - x23*x87 + x24*x87 - x25*x39 + x29*x42 - x29*x77 - x32*x41 - x45*x65 - x45*x66;
    res[11]=-(*this)[0]*x158 + (*this)[11]*x149 + (*this)[13]*x152 - (*this)[14]*x159 - (*this)[15]*x160 + (*this)[3]*x156 + (*this)[4]*x150 + (*this)[5]*x154 + (*this)[7]*x157 - (*this)[8]*x161 - (*this)[9]*x162 - A[0]*x92 + A[0]*x94 + x151*x64 + x153*x38 + x155*x88;
    res[12]=(*this)[0]*x154 + (*this)[12]*x149 - (*this)[13]*x150 + (*this)[14]*x156 - (*this)[15]*x161 - (*this)[3]*x159 - (*this)[4]*x152 - (*this)[5]*x158 - (*this)[7]*x162 - (*this)[8]*x160 + (*this)[9]*x157 + A[0]*x102 + x151*x60 - x155*x41 + x163*x48 + x163*x73;
    res[13]=(*this)[0]*x164 + (*this)[10]*x154 + (*this)[13]*x149 - (*this)[15]*x162 - (*this)[1]*x156 - (*this)[5]*x165 + (*this)[7]*x161 - (*this)[8]*x157 - (*this)[9]*x160 + A[0]*x30 - A[0]*x43 + A[0]*x59 - A[0]*x72 + x151*x49 - x153*x88 + x163*x75;
    res[14]=(*this)[0]*x157 - (*this)[11]*x159 - (*this)[12]*x156 - (*this)[15]*x165 - (*this)[1]*x150 + (*this)[2]*x152 + (*this)[5]*x162 - (*this)[7]*x158 + (*this)[8]*x164 - (*this)[9]*x154 + A[0]*x37 + A[0]*x54 - A[0]*x67 + A[0]*x76 - x151*x81 + x153*x60;
    res[15]=A[7]*x0 + A[7]*x1 + A[7]*x10 - A[7]*x11 + A[7]*x12 - A[7]*x13 - A[7]*x14 - A[7]*x15 - A[7]*x2 - A[7]*x3 - A[7]*x4 - A[7]*x5 + A[7]*x6 + A[7]*x7 + A[7]*x8 + A[7]*x9 - x163*x17 - x163*x18 - x163*x19 + x163*x20 - x163*x21 - x163*x22 - x163*x23 - x163*x24;
    return res;
};

//-----------------------------------
// (R130MV, scalar) binary operations
//-----------------------------------

template<typename T>
inline R130MV<T> operator*(const R130MV<T> &A, const T &B) {
    R130MV<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    res[4]=A[4]*B;
    res[5]=A[5]*B;
    res[6]=A[6]*B;
    res[7]=A[7]*B;
    res[8]=A[8]*B;
    res[9]=A[9]*B;
    res[10]=A[10]*B;
    res[11]=A[11]*B;
    res[12]=A[12]*B;
    res[13]=A[13]*B;
    res[14]=A[14]*B;
    res[15]=A[15]*B;
    return res;
};
template<typename T>
inline R130MV<T> operator*(const T &A, const R130MV<T> &B) {
    return B*A;
};
template<typename T>
inline R130MV<T> operator/(const R130MV<T> &A, const T &B) {
    return A * (1.0 / B);
};

//------------------------------------
// (Rotation, Boost) binary operations
//------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[1];
    res[6]=B[2];
    res[7]=B[3];
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const Boost<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[1];
    res[6]=-B[2];
    res[7]=-B[3];
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotation<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1] - A[2]*B[3] + A[3]*B[2];
    res[2]=A[0]*B[2] + A[1]*B[3] - A[3]*B[1];
    res[3]=A[0]*B[3] - A[1]*B[2] + A[2]*B[1];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    res[6]=A[3]*B[0];
    res[7]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotation<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[0]*B[1];
    res[2]=A[0]*B[2];
    res[3]=A[0]*B[3];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    res[6]=A[3]*B[0];
    res[7]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const Rotation<T> &A, const Boost<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline Boost<T> Rotation<T>::conjugate(const Boost<T> &A) const {
    Boost<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[3]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[3];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*(x0 + x1 + x2 + x3);
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x6 - A[2]*x9 + x7*x8;
    res[2]=(*this)[3]*x5 + A[1]*x6 + A[1]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x8;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x7 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//-------------------------------------
// (Rotation, R130B0) binary operations
//-------------------------------------

template<typename T>
inline Rotation<T> operator+(const Rotation<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const Rotation<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const Rotation<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline Rotation<T> operator/(const Rotation<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    return res;
};

template<typename T>
inline Rotation<T> operator|(const Rotation<T> &A, const R130B0<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotation<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B0<T> Rotation<T>::conjugate(const R130B0<T> &A) const {
    R130B0<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    return res;
};

//-------------------------------------
// (Rotation, R130B1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1] - A[2]*B[3] + A[3]*B[2];
    res[3]=A[0]*B[2] + A[1]*B[3] - A[3]*B[1];
    res[4]=A[0]*B[3] - A[1]*B[2] + A[2]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const Rotation<T> &A, const R130B1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130B1<T> Rotation<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[3]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[3];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*(x0 + x1 + x2 + x3);
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x6 - A[2]*x9 + x7*x8;
    res[2]=(*this)[3]*x5 + A[1]*x6 + A[1]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x8;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x7 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//----------------------------------------
// (Rotation, R130B1Sm1) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[3]=A[0]*B[1] + A[1]*B[2] - A[3]*B[0];
    res[4]=A[0]*B[2] - A[1]*B[1] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[0]*B[1];
    res[4]=A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sm1<T> operator^(const Rotation<T> &A, const R130B1Sm1<T> &B) {
    R130B1Sm1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B1Sm1<T> Rotation<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[2]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[2];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x5 + A[0]*x0 + A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x6 - A[1]*x9 + x7*x8;
    res[1]=(*this)[3]*x5 + A[0]*x6 + A[0]*x9 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 - x10*x8;
    res[2]=-(*this)[0]*A[0]*x4 + (*this)[1]*A[0]*x7 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x4 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3;
    return res;
};

//----------------------------------------
// (Rotation, R130B1Sp1) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
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
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
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
    res[12]=A[1]*B[0];
    res[13]=A[2]*B[0];
    res[14]=A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotation<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1Sp1<T> Rotation<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    return res;
};

//-------------------------------------
// (Rotation, R130B2) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=A[1] + B[3];
    res[9]=A[2] + B[4];
    res[10]=A[3] + B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B2<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=A[1] - B[3];
    res[9]=A[2] - B[4];
    res[10]=A[3] - B[5];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotation<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=-A[1]*B[3] - A[2]*B[4] - A[3]*B[5];
    res[1]=A[0]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[2]=A[0]*B[1] + A[1]*B[2] - A[3]*B[0];
    res[3]=A[0]*B[2] - A[1]*B[1] + A[2]*B[0];
    res[4]=A[0]*B[3] - A[2]*B[5] + A[3]*B[4];
    res[5]=A[0]*B[4] + A[1]*B[5] - A[3]*B[3];
    res[6]=A[0]*B[5] - A[1]*B[4] + A[2]*B[3];
    res[7]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotation<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=-A[1]*B[3] - A[2]*B[4] - A[3]*B[5];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=A[0]*B[4];
    res[6]=A[0]*B[5];
    res[7]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotation<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    res[3]=-A[2]*B[5] + A[3]*B[4];
    res[4]=A[1]*B[5] - A[3]*B[3];
    res[5]=-A[1]*B[4] + A[2]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> Rotation<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[2]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[2];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[1]*x7;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[0]*x4;
    res[0]=(*this)[0]*x5 + A[0]*x0 + A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x6 - A[1]*x9 + x7*x8;
    res[1]=(*this)[3]*x5 + A[0]*x6 + A[0]*x9 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 - x10*x8;
    res[2]=A[0]*x12 - A[0]*x14 + A[1]*x11 + A[1]*x13 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3;
    res[3]=A[3]*x0 + A[3]*x1 - A[3]*x2 - A[3]*x3 + A[4]*x6 - A[4]*x9 + A[5]*x12 + A[5]*x14;
    res[4]=A[3]*x6 + A[3]*x9 + A[4]*x0 - A[4]*x1 + A[4]*x2 - A[4]*x3 - A[5]*x11 + A[5]*x13;
    res[5]=A[3]*x12 - A[3]*x14 + A[4]*x11 + A[4]*x13 + A[5]*x0 - A[5]*x1 - A[5]*x2 + A[5]*x3;
    return res;
};

//----------------------------------------
// (Rotation, R130B2Sm1) binary operations
//----------------------------------------

template<typename T>
inline Rotation<T> operator+(const Rotation<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const Rotation<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const Rotation<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=A[0]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[2]=A[0]*B[1] + A[1]*B[2] - A[3]*B[0];
    res[3]=A[0]*B[2] - A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline Rotation<T> operator|(const Rotation<T> &A, const R130B2Sm1<T> &B) {
    Rotation<T> res;
    res[0]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const Rotation<T> &A, const R130B2Sm1<T> &B) {
    R130B2Sm1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> Rotation<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[2]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[2];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x5 + A[0]*x0 + A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x6 - A[1]*x9 + x7*x8;
    res[1]=(*this)[3]*x5 + A[0]*x6 + A[0]*x9 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 - x10*x8;
    res[2]=-(*this)[0]*A[0]*x4 + (*this)[1]*A[0]*x7 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x4 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3;
    return res;
};

//----------------------------------------
// (Rotation, R130B2Sp1) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=B[0];
    res[6]=B[1];
    res[7]=B[2];
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-B[0];
    res[6]=-B[1];
    res[7]=-B[2];
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[6]=A[0]*B[1] + A[1]*B[2] - A[3]*B[0];
    res[7]=A[0]*B[2] - A[1]*B[1] + A[2]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B2Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[0]*B[0];
    res[6]=A[0]*B[1];
    res[7]=A[0]*B[2];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2Sp1<T> operator^(const Rotation<T> &A, const R130B2Sp1<T> &B) {
    R130B2Sp1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2Sp1<T> Rotation<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2Sp1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[2]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[2];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x5 + A[0]*x0 + A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x6 - A[1]*x9 + x7*x8;
    res[1]=(*this)[3]*x5 + A[0]*x6 + A[0]*x9 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 - x10*x8;
    res[2]=-(*this)[0]*A[0]*x4 + (*this)[1]*A[0]*x7 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x4 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3;
    return res;
};

//-------------------------------------
// (Rotation, R130B3) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[1] - A[2]*B[3] + A[3]*B[2];
    res[13]=A[0]*B[2] + A[1]*B[3] - A[3]*B[1];
    res[14]=A[0]*B[3] - A[1]*B[2] + A[2]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const Rotation<T> &A, const R130B3<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline R130B3<T> Rotation<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[3]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[3];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*(x0 + x1 + x2 + x3);
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x6 - A[2]*x9 + x7*x8;
    res[2]=(*this)[3]*x5 + A[1]*x6 + A[1]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x8;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x7 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//----------------------------------------
// (Rotation, R130B3Sm1) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
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
    res[12]=A[0]*B[0] - A[2]*B[2] + A[3]*B[1];
    res[13]=A[0]*B[1] + A[1]*B[2] - A[3]*B[0];
    res[14]=A[0]*B[2] - A[1]*B[1] + A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
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
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sm1<T> operator^(const Rotation<T> &A, const R130B3Sm1<T> &B) {
    R130B3Sm1<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    return res;
};

template<typename T>
inline R130B3Sm1<T> Rotation<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3Sm1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[2]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[2];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x5 + A[0]*x0 + A[0]*x1 - A[0]*x2 - A[0]*x3 + A[1]*x6 - A[1]*x9 + x7*x8;
    res[1]=(*this)[3]*x5 + A[0]*x6 + A[0]*x9 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 - x10*x8;
    res[2]=-(*this)[0]*A[0]*x4 + (*this)[1]*A[0]*x7 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x4 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3;
    return res;
};

//----------------------------------------
// (Rotation, R130B3Sp1) binary operations
//----------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotation<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3Sp1<T> Rotation<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3Sp1<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    return res;
};

//-------------------------------------
// (Rotation, R130B4) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=A[1];
    res[9]=A[2];
    res[10]=A[3];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=-B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[0];
    res[6]=-A[2]*B[0];
    res[7]=-A[3]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*B[0];
    res[6]=-A[2]*B[0];
    res[7]=-A[3]*B[0];
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotation<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130B4<T> Rotation<T>::conjugate(const R130B4<T> &A) const {
    R130B4<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2));
    return res;
};

//-------------------------------------
// (Rotation, R130MV) binary operations
//-------------------------------------

template<typename T>
inline Rotation<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=(*this)[0];
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator+(const Rotation<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=B[5];
    res[6]=B[6];
    res[7]=B[7];
    res[8]=A[1] + B[8];
    res[9]=A[2] + B[9];
    res[10]=A[3] + B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotation<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=-B[5];
    res[6]=-B[6];
    res[7]=-B[7];
    res[8]=A[1] - B[8];
    res[9]=A[2] - B[9];
    res[10]=A[3] - B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=-B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotation<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[1]*B[8] - A[2]*B[9] - A[3]*B[10];
    res[1]=A[0]*B[1] - A[1]*B[12] - A[2]*B[13] - A[3]*B[14];
    res[2]=A[0]*B[2] + A[1]*B[11] - A[2]*B[4] + A[3]*B[3];
    res[3]=A[0]*B[3] + A[1]*B[4] + A[2]*B[11] - A[3]*B[2];
    res[4]=A[0]*B[4] - A[1]*B[3] + A[2]*B[2] + A[3]*B[11];
    res[5]=A[0]*B[5] - A[1]*B[15] - A[2]*B[7] + A[3]*B[6];
    res[6]=A[0]*B[6] + A[1]*B[7] - A[2]*B[15] - A[3]*B[5];
    res[7]=A[0]*B[7] - A[1]*B[6] + A[2]*B[5] - A[3]*B[15];
    res[8]=A[0]*B[8] + A[1]*B[0] - A[2]*B[10] + A[3]*B[9];
    res[9]=A[0]*B[9] + A[1]*B[10] + A[2]*B[0] - A[3]*B[8];
    res[10]=A[0]*B[10] - A[1]*B[9] + A[2]*B[8] + A[3]*B[0];
    res[11]=A[0]*B[11] - A[1]*B[2] - A[2]*B[3] - A[3]*B[4];
    res[12]=A[0]*B[12] + A[1]*B[1] - A[2]*B[14] + A[3]*B[13];
    res[13]=A[0]*B[13] + A[1]*B[14] + A[2]*B[1] - A[3]*B[12];
    res[14]=A[0]*B[14] - A[1]*B[13] + A[2]*B[12] + A[3]*B[1];
    res[15]=A[0]*B[15] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotation<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] - A[1]*B[8] - A[2]*B[9] - A[3]*B[10];
    res[1]=A[0]*B[1] - A[1]*B[12] - A[2]*B[13] - A[3]*B[14];
    res[2]=A[0]*B[2] + A[1]*B[11];
    res[3]=A[0]*B[3] + A[2]*B[11];
    res[4]=A[0]*B[4] + A[3]*B[11];
    res[5]=A[0]*B[5] - A[1]*B[15];
    res[6]=A[0]*B[6] - A[2]*B[15];
    res[7]=A[0]*B[7] - A[3]*B[15];
    res[8]=A[0]*B[8] + A[1]*B[0];
    res[9]=A[0]*B[9] + A[2]*B[0];
    res[10]=A[0]*B[10] + A[3]*B[0];
    res[11]=A[0]*B[11] - A[1]*B[2] - A[2]*B[3] - A[3]*B[4];
    res[12]=A[0]*B[12] + A[1]*B[1];
    res[13]=A[0]*B[13] + A[2]*B[1];
    res[14]=A[0]*B[14] + A[3]*B[1];
    res[15]=A[0]*B[15] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotation<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[2]*B[4] + A[3]*B[3];
    res[3]=A[1]*B[4] - A[3]*B[2];
    res[4]=-A[1]*B[3] + A[2]*B[2];
    res[5]=-A[2]*B[7] + A[3]*B[6];
    res[6]=A[1]*B[7] - A[3]*B[5];
    res[7]=-A[1]*B[6] + A[2]*B[5];
    res[8]=-A[2]*B[10] + A[3]*B[9];
    res[9]=A[1]*B[10] - A[3]*B[8];
    res[10]=-A[1]*B[9] + A[2]*B[8];
    res[11]=0;
    res[12]=-A[2]*B[14] + A[3]*B[13];
    res[13]=A[1]*B[14] - A[3]*B[12];
    res[14]=-A[1]*B[13] + A[2]*B[12];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> Rotation<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[4]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[4];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[1]*x7;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[0]*x4;
    res[0]=A[0]*(x0 + x1 + x2 + x3);
    res[1]=A[1]*(x0 + x1 + x2 + x3);
    res[2]=(*this)[0]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[3]*x6 - A[3]*x9 + x7*x8;
    res[3]=(*this)[3]*x5 + A[2]*x6 + A[2]*x9 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - x10*x8;
    res[4]=A[2]*x12 - A[2]*x14 + A[3]*x11 + A[3]*x13 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3;
    res[5]=A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[6]*x6 - A[6]*x9 + A[7]*x12 + A[7]*x14;
    res[6]=A[5]*x6 + A[5]*x9 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3 - A[7]*x11 + A[7]*x13;
    res[7]=A[5]*x12 - A[5]*x14 + A[6]*x11 + A[6]*x13 + A[7]*x0 - A[7]*x1 - A[7]*x2 + A[7]*x3;
    res[8]=A[10]*x12 + A[10]*x14 + A[8]*x0 + A[8]*x1 - A[8]*x2 - A[8]*x3 + A[9]*x6 - A[9]*x9;
    res[9]=-A[10]*x11 + A[10]*x13 + A[8]*x6 + A[8]*x9 + A[9]*x0 - A[9]*x1 + A[9]*x2 - A[9]*x3;
    res[10]=A[10]*x0 - A[10]*x1 - A[10]*x2 + A[10]*x3 + A[8]*x12 - A[8]*x14 + A[9]*x11 + A[9]*x13;
    res[11]=A[11]*(x0 + x1 + x2 + x3);
    res[12]=A[12]*x0 + A[12]*x1 - A[12]*x2 - A[12]*x3 + A[13]*x6 - A[13]*x9 + A[14]*x12 + A[14]*x14;
    res[13]=A[12]*x6 + A[12]*x9 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 - A[14]*x11 + A[14]*x13;
    res[14]=A[12]*x12 - A[12]*x14 + A[13]*x11 + A[13]*x13 + A[14]*x0 - A[14]*x1 - A[14]*x2 + A[14]*x3;
    res[15]=A[15]*(x0 + x1 + x2 + x3);
    return res;
};

//---------------------------------------
// (Rotation, Rotation) binary operations
//---------------------------------------

template<typename T>
inline Rotation<T> operator+(const Rotation<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator-(const Rotation<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    return res;
};

template<typename T>
inline Rotation<T> operator*(const Rotation<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0] - A[2]*B[3] + A[3]*B[2];
    res[2]=A[0]*B[2] + A[1]*B[3] + A[2]*B[0] - A[3]*B[1];
    res[3]=A[0]*B[3] - A[1]*B[2] + A[2]*B[1] + A[3]*B[0];
    return res;
};

template<typename T>
inline Rotation<T> operator|(const Rotation<T> &A, const Rotation<T> &B) {
    Rotation<T> res;
    res[0]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0];
    res[2]=A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[3] + A[3]*B[0];
    return res;
};

template<typename T>
inline R130B2Sm1<T> operator^(const Rotation<T> &A, const Rotation<T> &B) {
    R130B2Sm1<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    return res;
};

template<typename T>
inline Rotation<T> Rotation<T>::conjugate(const Rotation<T> &A) const {
    Rotation<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[3]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[3];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*(x0 + x1 + x2 + x3);
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x6 - A[2]*x9 + x7*x8;
    res[2]=(*this)[3]*x5 + A[1]*x6 + A[1]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x8;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x7 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    return res;
};

//------------------------------------
// (Rotation, Rotor) binary operations
//------------------------------------

template<typename T>
inline Rotation<T>::operator Rotor<T>() const {
    Rotor<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=(*this)[1];
    res[5]=(*this)[2];
    res[6]=(*this)[3];
    res[7]=0;
    return res;
};

template<typename T>
inline Rotor<T> operator+(const Rotation<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=A[1] + B[4];
    res[5]=A[2] + B[5];
    res[6]=A[3] + B[6];
    res[7]=B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotation<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=A[1] - B[4];
    res[5]=A[2] - B[5];
    res[6]=A[3] - B[6];
    res[7]=-B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotation<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] - A[1]*B[4] - A[2]*B[5] - A[3]*B[6];
    res[1]=A[0]*B[1] - A[1]*B[7] - A[2]*B[3] + A[3]*B[2];
    res[2]=A[0]*B[2] + A[1]*B[3] - A[2]*B[7] - A[3]*B[1];
    res[3]=A[0]*B[3] - A[1]*B[2] + A[2]*B[1] - A[3]*B[7];
    res[4]=A[0]*B[4] + A[1]*B[0] - A[2]*B[6] + A[3]*B[5];
    res[5]=A[0]*B[5] + A[1]*B[6] + A[2]*B[0] - A[3]*B[4];
    res[6]=A[0]*B[6] - A[1]*B[5] + A[2]*B[4] + A[3]*B[0];
    res[7]=A[0]*B[7] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotation<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] - A[1]*B[4] - A[2]*B[5] - A[3]*B[6];
    res[1]=A[0]*B[1] - A[1]*B[7];
    res[2]=A[0]*B[2] - A[2]*B[7];
    res[3]=A[0]*B[3] - A[3]*B[7];
    res[4]=A[0]*B[4] + A[1]*B[0];
    res[5]=A[0]*B[5] + A[2]*B[0];
    res[6]=A[0]*B[6] + A[3]*B[0];
    res[7]=A[0]*B[7] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotation<T> &A, const Rotor<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=-A[2]*B[6] + A[3]*B[5];
    res[4]=A[1]*B[6] - A[3]*B[4];
    res[5]=-A[1]*B[5] + A[2]*B[4];
    return res;
};

template<typename T>
inline Rotor<T> Rotation<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = 2.0*(*this)[2];
    T x5 = A[3]*x4;
    T x6 = (*this)[1]*x4;
    T x7 = 2.0*(*this)[3];
    T x8 = (*this)[1]*A[3];
    T x9 = (*this)[0]*x7;
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[1]*x7;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[0]*x4;
    res[0]=A[0]*(x0 + x1 + x2 + x3);
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x6 - A[2]*x9 + x7*x8;
    res[2]=(*this)[3]*x5 + A[1]*x6 + A[1]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x8;
    res[3]=A[1]*x12 - A[1]*x14 + A[2]*x11 + A[2]*x13 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    res[4]=A[4]*x0 + A[4]*x1 - A[4]*x2 - A[4]*x3 + A[5]*x6 - A[5]*x9 + A[6]*x12 + A[6]*x14;
    res[5]=A[4]*x6 + A[4]*x9 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[5]*x3 - A[6]*x11 + A[6]*x13;
    res[6]=A[4]*x12 - A[4]*x14 + A[5]*x11 + A[5]*x13 + A[6]*x0 - A[6]*x1 - A[6]*x2 + A[6]*x3;
    res[7]=A[7]*(x0 + x1 + x2 + x3);
    return res;
};

//-------------------------------------
// (Rotation, scalar) binary operations
//-------------------------------------

template<typename T>
inline Rotation<T> operator*(const Rotation<T> &A, const T &B) {
    Rotation<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    return res;
};
template<typename T>
inline Rotation<T> operator*(const T &A, const Rotation<T> &B) {
    return B*A;
};
template<typename T>
inline Rotation<T> operator/(const Rotation<T> &A, const T &B) {
    return A * (1.0 / B);
};

//---------------------------------
// (Rotor, Boost) binary operations
//---------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0] - A[5]*B[3] + A[6]*B[2];
    res[2]=A[0]*B[2] + A[2]*B[0] + A[4]*B[3] - A[6]*B[1];
    res[3]=A[0]*B[3] + A[3]*B[0] - A[4]*B[2] + A[5]*B[1];
    res[4]=A[2]*B[3] - A[3]*B[2] + A[4]*B[0] + A[7]*B[1];
    res[5]=-A[1]*B[3] + A[3]*B[1] + A[5]*B[0] + A[7]*B[2];
    res[6]=A[1]*B[2] - A[2]*B[1] + A[6]*B[0] + A[7]*B[3];
    res[7]=A[4]*B[1] + A[5]*B[2] + A[6]*B[3] + A[7]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const Boost<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[1]=A[0]*B[1] + A[1]*B[0];
    res[2]=A[0]*B[2] + A[2]*B[0];
    res[3]=A[0]*B[3] + A[3]*B[0];
    res[4]=A[4]*B[0] + A[7]*B[1];
    res[5]=A[5]*B[0] + A[7]*B[2];
    res[6]=A[6]*B[0] + A[7]*B[3];
    res[7]=A[4]*B[1] + A[5]*B[2] + A[6]*B[3] + A[7]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotor<T> &A, const Boost<T> &B) {
    R130B2<T> res;
    res[0]=-A[5]*B[3] + A[6]*B[2];
    res[1]=A[4]*B[3] - A[6]*B[1];
    res[2]=-A[4]*B[2] + A[5]*B[1];
    res[3]=A[2]*B[3] - A[3]*B[2];
    res[4]=-A[1]*B[3] + A[3]*B[1];
    res[5]=A[1]*B[2] - A[2]*B[1];
    return res;
};

template<typename T>
inline Rotor<T> Rotor<T>::conjugate(const Boost<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[6], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = std::pow((*this)[3], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*A[3];
    T x9 = (*this)[5]*x8;
    T x10 = (*this)[7]*x8;
    T x11 = 2.0*A[2];
    T x12 = (*this)[4]*(*this)[5];
    T x13 = (*this)[4]*x8;
    T x14 = (*this)[0]*(*this)[6];
    T x15 = (*this)[1]*(*this)[2];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*(*this)[7];
    T x18 = 2.0*A[1];
    T x19 = (*this)[0]*x11;
    T x20 = (*this)[1]*x11;
    T x21 = (*this)[6]*x18;
    T x22 = (*this)[6]*x11;
    T x23 = (*this)[0]*x18;
    T x24 = (*this)[1]*x18;
    T x25 = (*this)[2]*x11;
    T x26 = (*this)[2]*x18;
    T x27 = (*this)[0]*x8;
    T x28 = (*this)[6]*x8;
    T x29 = 2.0*A[0];
    res[0]=A[0]*(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7);
    res[1]=(*this)[0]*x9 - (*this)[1]*x16 + (*this)[2]*x10 + (*this)[6]*x13 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 + x11*x12 - x11*x14 - x11*x15 - x11*x17;
    res[2]=-(*this)[0]*x13 - (*this)[1]*x10 - (*this)[2]*x16 + (*this)[6]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7 + x12*x18 + x14*x18 - x15*x18 + x17*x18;
    res[3]=-(*this)[3]*x24 - (*this)[3]*x25 + (*this)[4]*x19 + (*this)[4]*x21 + (*this)[5]*x22 - (*this)[5]*x23 + (*this)[7]*x20 - (*this)[7]*x26 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7;
    res[4]=-(*this)[1]*x28 - (*this)[2]*x27 - (*this)[3]*x13 + (*this)[3]*x19 + (*this)[3]*x21 - (*this)[4]*x24 - (*this)[4]*x25 - (*this)[5]*x20 + (*this)[5]*x26 - (*this)[7]*x22 + (*this)[7]*x23 + (*this)[7]*x9;
    res[5]=(*this)[1]*x27 - (*this)[2]*x28 + (*this)[3]*x22 - (*this)[3]*x23 - (*this)[3]*x9 - (*this)[4]*x10 + (*this)[4]*x20 - (*this)[4]*x26 - (*this)[5]*x24 - (*this)[5]*x25 + (*this)[7]*x19 + (*this)[7]*x21;
    res[6]=(*this)[0]*x10 + (*this)[1]*x13 - (*this)[1]*x19 - (*this)[1]*x21 - (*this)[2]*x22 + (*this)[2]*x23 + (*this)[2]*x9 - (*this)[3]*(*this)[4]*x18 - (*this)[3]*(*this)[5]*x11 + (*this)[4]*(*this)[7]*x11 - (*this)[5]*(*this)[7]*x18 - (*this)[6]*x16;
    res[7]=x29*((*this)[0]*(*this)[7] - (*this)[1]*(*this)[4] - (*this)[2]*(*this)[5] - (*this)[3]*(*this)[6]);
    return res;
};

//----------------------------------
// (Rotor, R130B0) binary operations
//----------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const R130B0<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const R130B0<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const R130B0<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[4]*B[0];
    res[5]=A[5]*B[0];
    res[6]=A[6]*B[0];
    res[7]=A[7]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator/(const Rotor<T> &A, const R130B0<T> &B) {
    Rotor<T> res;
    T x0 = 1.0/B[0];
    res[0]=A[0]*x0;
    res[1]=A[1]*x0;
    res[2]=A[2]*x0;
    res[3]=A[3]*x0;
    res[4]=A[4]*x0;
    res[5]=A[5]*x0;
    res[6]=A[6]*x0;
    res[7]=A[7]*x0;
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const R130B0<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0];
    res[1]=A[1]*B[0];
    res[2]=A[2]*B[0];
    res[3]=A[3]*B[0];
    res[4]=A[4]*B[0];
    res[5]=A[5]*B[0];
    res[6]=A[6]*B[0];
    res[7]=A[7]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B0<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> Rotor<T>::conjugate(const R130B0<T> &A) const {
    R130MV<T> res;
    res[0]=A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) - std::pow((*this)[7], 2));
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
    res[15]=2.0*A[0]*((*this)[0]*(*this)[7] - (*this)[1]*(*this)[4] - (*this)[2]*(*this)[5] - (*this)[3]*(*this)[6]);
    return res;
};

//----------------------------------
// (Rotor, R130B1) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=B[1];
    res[3]=B[2];
    res[4]=B[3];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=-B[1];
    res[3]=-B[2];
    res[4]=-B[3];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[2]=A[0]*B[1] + A[1]*B[0] - A[5]*B[3] + A[6]*B[2];
    res[3]=A[0]*B[2] + A[2]*B[0] + A[4]*B[3] - A[6]*B[1];
    res[4]=A[0]*B[3] + A[3]*B[0] - A[4]*B[2] + A[5]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[4]*B[1] - A[5]*B[2] - A[6]*B[3] - A[7]*B[0];
    res[12]=A[2]*B[3] - A[3]*B[2] + A[4]*B[0] + A[7]*B[1];
    res[13]=-A[1]*B[3] + A[3]*B[1] + A[5]*B[0] + A[7]*B[2];
    res[14]=A[1]*B[2] - A[2]*B[1] + A[6]*B[0] + A[7]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[0]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[4]*B[1] - A[5]*B[2] - A[6]*B[3];
    res[12]=A[2]*B[3] - A[3]*B[2] + A[4]*B[0];
    res[13]=-A[1]*B[3] + A[3]*B[1] + A[5]*B[0];
    res[14]=A[1]*B[2] - A[2]*B[1] + A[6]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3];
    res[2]=A[1]*B[0] - A[5]*B[3] + A[6]*B[2];
    res[3]=A[2]*B[0] + A[4]*B[3] - A[6]*B[1];
    res[4]=A[3]*B[0] - A[4]*B[2] + A[5]*B[1];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[7]*B[0];
    res[12]=A[7]*B[1];
    res[13]=A[7]*B[2];
    res[14]=A[7]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> Rotor<T>::conjugate(const R130B1<T> &A) const {
    R130B1<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = std::pow((*this)[4], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = std::pow((*this)[6], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = A[2]*x8;
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[4]*A[3];
    T x12 = 2.0*(*this)[3];
    T x13 = (*this)[5]*x12;
    T x14 = (*this)[0]*x8;
    T x15 = (*this)[0]*A[2];
    T x16 = (*this)[0]*A[3];
    T x17 = A[3]*x8;
    T x18 = (*this)[6]*x10;
    T x19 = A[2]*x12;
    T x20 = 2.0*(*this)[7];
    T x21 = (*this)[4]*x20;
    T x22 = (*this)[5]*A[2];
    T x23 = (*this)[6]*A[3];
    T x24 = 2.0*(*this)[5];
    T x25 = 2.0*(*this)[4];
    T x26 = 2.0*(*this)[6];
    T x27 = A[3]*x10;
    T x28 = (*this)[0]*A[1];
    T x29 = A[1]*x8;
    T x30 = A[0]*x12;
    T x31 = A[0]*x10;
    T x32 = A[0]*x8;
    T x33 = (*this)[7]*A[1];
    T x34 = A[0]*x20;
    res[0]=-(*this)[4]*x19 - (*this)[5]*x17 + (*this)[6]*x9 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 + A[1]*x13 - A[1]*x14 - A[1]*x18 - A[1]*x21 + x10*x11 - x10*x15 - x12*x16 - x20*x22 - x20*x23;
    res[1]=(*this)[2]*x9 + (*this)[3]*x17 + (*this)[7]*x19 - (*this)[7]*x27 - A[0]*x13 - A[0]*x14 + A[0]*x18 - A[0]*x21 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 - A[1]*x5 - A[1]*x6 + A[1]*x7 + x11*x26 - x15*x26 + x16*x24 + x22*x25;
    res[2]=-2.0*(*this)[0]*x11 - (*this)[0]*x31 + (*this)[2]*x29 + (*this)[3]*x27 + (*this)[4]*A[1]*x24 + (*this)[4]*x30 - (*this)[5]*x34 - (*this)[6]*x32 + (*this)[7]*x17 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 - x12*x33 + x23*x24 + x26*x28;
    res[3]=-(*this)[0]*x30 + (*this)[3]*A[2]*x10 + (*this)[3]*x29 - (*this)[4]*x31 + (*this)[5]*x32 + (*this)[6]*A[1]*x25 - (*this)[6]*x34 - (*this)[7]*x9 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 - A[3]*x4 - A[3]*x5 + A[3]*x6 + A[3]*x7 + x10*x33 + x15*x25 + x22*x26 - x24*x28;
    return res;
};

//-------------------------------------
// (Rotor, R130B1Sm1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=B[0];
    res[3]=B[1];
    res[4]=B[2];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=-B[0];
    res[3]=-B[1];
    res[4]=-B[2];
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[2]=A[0]*B[0] - A[5]*B[2] + A[6]*B[1];
    res[3]=A[0]*B[1] + A[4]*B[2] - A[6]*B[0];
    res[4]=A[0]*B[2] - A[4]*B[1] + A[5]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[4]*B[0] - A[5]*B[1] - A[6]*B[2];
    res[12]=A[2]*B[2] - A[3]*B[1] + A[7]*B[0];
    res[13]=-A[1]*B[2] + A[3]*B[0] + A[7]*B[1];
    res[14]=A[1]*B[1] - A[2]*B[0] + A[7]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[0]*B[0];
    res[3]=A[0]*B[1];
    res[4]=A[0]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[4]*B[0] - A[5]*B[1] - A[6]*B[2];
    res[12]=A[2]*B[2] - A[3]*B[1];
    res[13]=-A[1]*B[2] + A[3]*B[0];
    res[14]=A[1]*B[1] - A[2]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B1Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[2]=-A[5]*B[2] + A[6]*B[1];
    res[3]=A[4]*B[2] - A[6]*B[0];
    res[4]=-A[4]*B[1] + A[5]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[7]*B[0];
    res[13]=A[7]*B[1];
    res[14]=A[7]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> Rotor<T>::conjugate(const R130B1Sm1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*(*this)[1];
    T x1 = A[1]*x0;
    T x2 = 2.0*(*this)[2];
    T x3 = (*this)[4]*A[2];
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[5]*A[0];
    T x6 = (*this)[0]*A[0];
    T x7 = (*this)[0]*A[1];
    T x8 = (*this)[0]*A[2];
    T x9 = A[2]*x0;
    T x10 = A[0]*x2;
    T x11 = A[1]*x4;
    T x12 = 2.0*(*this)[7];
    T x13 = (*this)[4]*A[0];
    T x14 = (*this)[5]*A[1];
    T x15 = (*this)[6]*A[2];
    T x16 = std::pow((*this)[0], 2);
    T x17 = std::pow((*this)[1], 2);
    T x18 = std::pow((*this)[4], 2);
    T x19 = std::pow((*this)[7], 2);
    T x20 = std::pow((*this)[2], 2);
    T x21 = std::pow((*this)[3], 2);
    T x22 = std::pow((*this)[5], 2);
    T x23 = std::pow((*this)[6], 2);
    T x24 = 2.0*(*this)[5];
    T x25 = 2.0*(*this)[4];
    T x26 = 2.0*(*this)[6];
    T x27 = A[2]*x2;
    T x28 = A[0]*x0;
    T x29 = 2.0*(*this)[0];
    res[0]=-(*this)[4]*x11 - (*this)[5]*x9 + (*this)[6]*x1 - (*this)[6]*x10 - x0*x6 - x12*x13 - x12*x14 - x12*x15 + x2*x3 - x2*x7 + x4*x5 - x4*x8;
    res[1]=(*this)[2]*x1 + (*this)[3]*x9 + (*this)[7]*x11 - (*this)[7]*x27 + A[0]*x16 + A[0]*x17 + A[0]*x18 + A[0]*x19 - A[0]*x20 - A[0]*x21 - A[0]*x22 - A[0]*x23 + x14*x25 + x24*x8 + x26*x3 - x26*x7;
    res[2]=(*this)[2]*x28 + (*this)[3]*x27 - (*this)[7]*A[0]*x4 + (*this)[7]*x9 + A[1]*x16 - A[1]*x17 - A[1]*x18 + A[1]*x19 + A[1]*x20 - A[1]*x21 + A[1]*x22 - A[1]*x23 + x15*x24 + x25*x5 + x26*x6 - x29*x3;
    res[3]=(*this)[3]*A[1]*x2 + (*this)[3]*x28 - (*this)[7]*x1 + (*this)[7]*x10 + A[2]*x16 - A[2]*x17 - A[2]*x18 + A[2]*x19 - A[2]*x20 + A[2]*x21 - A[2]*x22 + A[2]*x23 + x13*x26 + x14*x26 + x25*x7 - x29*x5;
    return res;
};

//-------------------------------------
// (Rotor, R130B1Sp1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=-B[0];
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[7]*B[0];
    res[12]=A[4]*B[0];
    res[13]=A[5]*B[0];
    res[14]=A[6]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[0]*B[0];
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
    res[12]=A[4]*B[0];
    res[13]=A[5]*B[0];
    res[14]=A[6]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B1Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[1]*B[0];
    res[3]=A[2]*B[0];
    res[4]=A[3]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[7]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130B1<T> Rotor<T>::conjugate(const R130B1Sp1<T> &A) const {
    R130B1<T> res;
    T x0 = 2.0*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) + std::pow((*this)[7], 2));
    res[1]=x0*(-(*this)[0]*(*this)[1] + (*this)[2]*(*this)[6] - (*this)[3]*(*this)[5] - (*this)[4]*(*this)[7]);
    res[2]=x0*(-(*this)[0]*(*this)[2] - (*this)[1]*(*this)[6] + (*this)[3]*(*this)[4] - (*this)[5]*(*this)[7]);
    res[3]=x0*(-(*this)[0]*(*this)[3] + (*this)[1]*(*this)[5] - (*this)[2]*(*this)[4] - (*this)[6]*(*this)[7]);
    return res;
};

//----------------------------------
// (Rotor, R130B2) binary operations
//----------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    res[4]=A[4] + B[3];
    res[5]=A[5] + B[4];
    res[6]=A[6] + B[5];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    res[4]=A[4] - B[3];
    res[5]=A[5] - B[4];
    res[6]=A[6] - B[5];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2] - A[4]*B[3] - A[5]*B[4] - A[6]*B[5];
    res[1]=A[0]*B[0] - A[2]*B[5] + A[3]*B[4] - A[5]*B[2] + A[6]*B[1] - A[7]*B[3];
    res[2]=A[0]*B[1] + A[1]*B[5] - A[3]*B[3] + A[4]*B[2] - A[6]*B[0] - A[7]*B[4];
    res[3]=A[0]*B[2] - A[1]*B[4] + A[2]*B[3] - A[4]*B[1] + A[5]*B[0] - A[7]*B[5];
    res[4]=A[0]*B[3] + A[2]*B[2] - A[3]*B[1] - A[5]*B[5] + A[6]*B[4] + A[7]*B[0];
    res[5]=A[0]*B[4] - A[1]*B[2] + A[3]*B[0] + A[4]*B[5] - A[6]*B[3] + A[7]*B[1];
    res[6]=A[0]*B[5] + A[1]*B[1] - A[2]*B[0] - A[4]*B[4] + A[5]*B[3] + A[7]*B[2];
    res[7]=A[1]*B[3] + A[2]*B[4] + A[3]*B[5] + A[4]*B[0] + A[5]*B[1] + A[6]*B[2];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const R130B2<T> &B) {
    Rotor<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2] - A[4]*B[3] - A[5]*B[4] - A[6]*B[5];
    res[1]=A[0]*B[0] - A[7]*B[3];
    res[2]=A[0]*B[1] - A[7]*B[4];
    res[3]=A[0]*B[2] - A[7]*B[5];
    res[4]=A[0]*B[3] + A[7]*B[0];
    res[5]=A[0]*B[4] + A[7]*B[1];
    res[6]=A[0]*B[5] + A[7]*B[2];
    res[7]=A[1]*B[3] + A[2]*B[4] + A[3]*B[5] + A[4]*B[0] + A[5]*B[1] + A[6]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotor<T> &A, const R130B2<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[5] + A[3]*B[4] - A[5]*B[2] + A[6]*B[1];
    res[1]=A[1]*B[5] - A[3]*B[3] + A[4]*B[2] - A[6]*B[0];
    res[2]=-A[1]*B[4] + A[2]*B[3] - A[4]*B[1] + A[5]*B[0];
    res[3]=A[2]*B[2] - A[3]*B[1] - A[5]*B[5] + A[6]*B[4];
    res[4]=-A[1]*B[2] + A[3]*B[0] + A[4]*B[5] - A[6]*B[3];
    res[5]=A[1]*B[1] - A[2]*B[0] - A[4]*B[4] + A[5]*B[3];
    return res;
};

template<typename T>
inline R130B2<T> Rotor<T>::conjugate(const R130B2<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[4], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = std::pow((*this)[6], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*(*this)[0];
    T x9 = (*this)[2]*A[5];
    T x10 = (*this)[5]*A[2];
    T x11 = 2.0*(*this)[1];
    T x12 = A[3]*x11;
    T x13 = (*this)[5]*A[4];
    T x14 = (*this)[6]*A[5];
    T x15 = 2.0*(*this)[2];
    T x16 = (*this)[4]*A[4];
    T x17 = (*this)[7]*A[2];
    T x18 = 2.0*(*this)[4];
    T x19 = (*this)[3]*A[5];
    T x20 = (*this)[5]*x18;
    T x21 = (*this)[6]*x18;
    T x22 = 2.0*(*this)[6];
    T x23 = (*this)[7]*A[4];
    T x24 = (*this)[3]*x8;
    T x25 = (*this)[6]*x8;
    T x26 = (*this)[7]*A[3];
    T x27 = (*this)[2]*x11;
    T x28 = (*this)[3]*A[2];
    T x29 = A[3]*x15;
    T x30 = (*this)[3]*x22;
    T x31 = 2.0*(*this)[7];
    T x32 = (*this)[3]*A[1];
    T x33 = (*this)[5]*A[5];
    T x34 = 2.0*(*this)[5];
    T x35 = (*this)[3]*A[0];
    T x36 = (*this)[7]*A[5];
    T x37 = (*this)[1]*x8;
    T x38 = (*this)[4]*x8;
    T x39 = A[1]*x11;
    T x40 = (*this)[6]*x15;
    T x41 = 2.0*(*this)[3];
    T x42 = A[1]*x22;
    T x43 = (*this)[2]*x8;
    T x44 = A[0]*x8;
    T x45 = (*this)[4]*x11;
    T x46 = A[0]*x15;
    T x47 = (*this)[6]*x11;
    T x48 = A[1]*x15;
    T x49 = (*this)[7]*A[1];
    T x50 = (*this)[5]*A[0];
    res[0]=(*this)[4]*x12 - (*this)[5]*x29 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7 + A[1]*x20 - A[1]*x25 - A[1]*x27 + A[2]*x21 - A[3]*x30 - A[4]*x24 + x10*x8 + x11*x13 + x11*x14 - x11*x28 + x15*x16 + x15*x17 + x18*x19 + x22*x23 - x26*x8 - x31*x32 - x31*x33 + x8*x9;
    res[1]=(*this)[4]*x29 + (*this)[5]*x12 + A[0]*x20 + A[0]*x25 - A[0]*x27 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 - A[1]*x6 - A[1]*x7 - A[2]*x38 + A[3]*x24 - A[4]*x30 - A[5]*x37 + x10*x22 - x11*x16 - x11*x17 + x13*x15 - x15*x28 + x18*x36 + x19*x34 - x22*x26 + x22*x9 - x23*x8 + x31*x35;
    res[2]=(*this)[3]*A[3]*x18 + (*this)[5]*x42 - (*this)[5]*x44 + (*this)[6]*x12 + (*this)[7]*x39 - (*this)[7]*x46 + A[0]*x21 + A[1]*x38 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7 - A[3]*x43 + A[4]*x37 + A[4]*x40 - A[5]*x45 - x11*x35 + x13*x41 + x14*x41 - x15*x32 - x16*x31 + x26*x34 - x34*x9 - x36*x8;
    res[3]=-(*this)[4]*x48 - (*this)[5]*x39 + (*this)[5]*x46 - (*this)[7]*x42 + (*this)[7]*x44 + A[0]*x30 - A[0]*x45 + A[1]*x24 - A[2]*x43 - A[2]*x47 + A[3]*x0 + A[3]*x1 + A[3]*x2 + A[3]*x3 - A[3]*x4 - A[3]*x5 - A[3]*x6 - A[3]*x7 - A[4]*x25 - A[4]*x27 + x10*x31 - x11*x19 + x13*x18 + x14*x18 - x18*x28 - x23*x41 + x31*x9 + x33*x8;
    res[4]=-(*this)[2]*x12 + (*this)[4]*x39 - (*this)[4]*x46 - (*this)[5]*x48 + (*this)[7]*A[0]*x22 - A[0]*x24 + A[1]*x30 + A[2]*x37 - A[2]*x40 + A[3]*x20 + A[3]*x25 + A[4]*x0 - A[4]*x1 + A[4]*x2 - A[4]*x3 + A[4]*x4 + A[4]*x5 - A[4]*x6 - A[4]*x7 - A[5]*x38 - x10*x41 - x11*x36 - x11*x50 + x14*x34 - x17*x18 + x26*x41 - x41*x9 + x49*x8;
    res[5]=-(*this)[3]*A[4]*x15 - (*this)[3]*x12 - (*this)[5]*A[3]*x8 + A[0]*x43 - A[0]*x47 - A[1]*x37 - A[1]*x40 + A[2]*x45 + A[3]*x21 + A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[5]*x4 - A[5]*x5 + A[5]*x6 - A[5]*x7 + x10*x15 + x11*x23 + x13*x22 - x15*x26 + x16*x8 + x17*x8 - x18*x35 + x18*x49 - x22*x28 - x31*x50 - x32*x34;
    return res;
};

//-------------------------------------
// (Rotor, R130B2Sm1) binary operations
//-------------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const R130B2Sm1<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4] + B[0];
    res[5]=A[5] + B[1];
    res[6]=A[6] + B[2];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const R130B2Sm1<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4] - B[0];
    res[5]=A[5] - B[1];
    res[6]=A[6] - B[2];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const R130B2Sm1<T> &B) {
    Rotor<T> res;
    res[0]=-A[4]*B[0] - A[5]*B[1] - A[6]*B[2];
    res[1]=-A[2]*B[2] + A[3]*B[1] - A[7]*B[0];
    res[2]=A[1]*B[2] - A[3]*B[0] - A[7]*B[1];
    res[3]=-A[1]*B[1] + A[2]*B[0] - A[7]*B[2];
    res[4]=A[0]*B[0] - A[5]*B[2] + A[6]*B[1];
    res[5]=A[0]*B[1] + A[4]*B[2] - A[6]*B[0];
    res[6]=A[0]*B[2] - A[4]*B[1] + A[5]*B[0];
    res[7]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const R130B2Sm1<T> &B) {
    Rotor<T> res;
    res[0]=-A[4]*B[0] - A[5]*B[1] - A[6]*B[2];
    res[1]=-A[7]*B[0];
    res[2]=-A[7]*B[1];
    res[3]=-A[7]*B[2];
    res[4]=A[0]*B[0];
    res[5]=A[0]*B[1];
    res[6]=A[0]*B[2];
    res[7]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotor<T> &A, const R130B2Sm1<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[2] + A[3]*B[1];
    res[1]=A[1]*B[2] - A[3]*B[0];
    res[2]=-A[1]*B[1] + A[2]*B[0];
    res[3]=-A[5]*B[2] + A[6]*B[1];
    res[4]=A[4]*B[2] - A[6]*B[0];
    res[5]=-A[4]*B[1] + A[5]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> Rotor<T>::conjugate(const R130B2Sm1<T> &A) const {
    R130B2<T> res;
    T x0 = 2.0*A[2];
    T x1 = (*this)[2]*x0;
    T x2 = 2.0*(*this)[1];
    T x3 = A[0]*x2;
    T x4 = (*this)[5]*A[1];
    T x5 = (*this)[1]*x0;
    T x6 = 2.0*A[1];
    T x7 = (*this)[2]*(*this)[4];
    T x8 = (*this)[3]*x0;
    T x9 = (*this)[7]*x6;
    T x10 = (*this)[0]*(*this)[3];
    T x11 = 2.0*A[0];
    T x12 = (*this)[0]*(*this)[7];
    T x13 = (*this)[5]*x11;
    T x14 = (*this)[3]*(*this)[6];
    T x15 = (*this)[7]*x0;
    T x16 = 2.0*x4;
    T x17 = A[1]*x2;
    T x18 = (*this)[6]*x11;
    T x19 = (*this)[6]*x6;
    T x20 = (*this)[3]*x11;
    T x21 = (*this)[2]*x11;
    T x22 = std::pow((*this)[0], 2);
    T x23 = std::pow((*this)[2], 2);
    T x24 = std::pow((*this)[3], 2);
    T x25 = std::pow((*this)[4], 2);
    T x26 = std::pow((*this)[1], 2);
    T x27 = std::pow((*this)[5], 2);
    T x28 = std::pow((*this)[6], 2);
    T x29 = std::pow((*this)[7], 2);
    T x30 = (*this)[5]*x0;
    T x31 = (*this)[4]*x0;
    res[0]=(*this)[0]*x1 - (*this)[2]*x13 + (*this)[4]*x3 + (*this)[4]*x8 - (*this)[5]*x15 + (*this)[6]*x5 + (*this)[6]*x9 - x10*x6 - x11*x12 - x11*x14 + x2*x4 + x6*x7;
    res[1]=-(*this)[0]*x5 - (*this)[0]*x9 + (*this)[2]*x16 + (*this)[4]*x15 - (*this)[4]*x17 + (*this)[5]*x3 + (*this)[5]*x8 + (*this)[6]*x1 - (*this)[7]*x18 + x10*x11 + x11*x7 - x14*x6;
    res[2]=(*this)[0]*x17 - (*this)[0]*x21 + (*this)[2]*x19 + (*this)[3]*x16 + (*this)[4]*x20 - (*this)[4]*x5 - (*this)[4]*x9 - (*this)[5]*x1 + (*this)[6]*x3 + (*this)[6]*x8 + (*this)[7]*x13 - x0*x12;
    res[3]=-(*this)[0]*x19 + (*this)[0]*x30 - (*this)[2]*x17 - (*this)[3]*x5 - (*this)[3]*x9 + (*this)[4]*x16 + (*this)[6]*x31 + (*this)[7]*x1 + A[0]*x22 + A[0]*x23 + A[0]*x24 + A[0]*x25 - A[0]*x26 - A[0]*x27 - A[0]*x28 - A[0]*x29;
    res[4]=(*this)[0]*x18 - (*this)[0]*x31 - (*this)[2]*x3 - (*this)[3]*x1 + (*this)[4]*x13 + (*this)[6]*x30 + (*this)[7]*x20 - (*this)[7]*x5 + A[1]*x22 - A[1]*x23 + A[1]*x24 - A[1]*x25 + A[1]*x26 + A[1]*x27 - A[1]*x28 - A[1]*x29;
    res[5]=(*this)[0]*(*this)[4]*x6 - (*this)[0]*x13 - (*this)[2]*(*this)[3]*x6 - (*this)[3]*x3 + (*this)[4]*x18 + (*this)[6]*x16 + (*this)[7]*x17 - (*this)[7]*x21 + A[2]*x22 + A[2]*x23 - A[2]*x24 - A[2]*x25 + A[2]*x26 - A[2]*x27 + A[2]*x28 - A[2]*x29;
    return res;
};

//-------------------------------------
// (Rotor, R130B2Sp1) binary operations
//-------------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const R130B2Sp1<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1] + B[0];
    res[2]=A[2] + B[1];
    res[3]=A[3] + B[2];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const R130B2Sp1<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1] - B[0];
    res[2]=A[2] - B[1];
    res[3]=A[3] - B[2];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const R130B2Sp1<T> &B) {
    Rotor<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0] - A[5]*B[2] + A[6]*B[1];
    res[2]=A[0]*B[1] + A[4]*B[2] - A[6]*B[0];
    res[3]=A[0]*B[2] - A[4]*B[1] + A[5]*B[0];
    res[4]=A[2]*B[2] - A[3]*B[1] + A[7]*B[0];
    res[5]=-A[1]*B[2] + A[3]*B[0] + A[7]*B[1];
    res[6]=A[1]*B[1] - A[2]*B[0] + A[7]*B[2];
    res[7]=A[4]*B[0] + A[5]*B[1] + A[6]*B[2];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const R130B2Sp1<T> &B) {
    Rotor<T> res;
    res[0]=A[1]*B[0] + A[2]*B[1] + A[3]*B[2];
    res[1]=A[0]*B[0];
    res[2]=A[0]*B[1];
    res[3]=A[0]*B[2];
    res[4]=A[7]*B[0];
    res[5]=A[7]*B[1];
    res[6]=A[7]*B[2];
    res[7]=A[4]*B[0] + A[5]*B[1] + A[6]*B[2];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotor<T> &A, const R130B2Sp1<T> &B) {
    R130B2<T> res;
    res[0]=-A[5]*B[2] + A[6]*B[1];
    res[1]=A[4]*B[2] - A[6]*B[0];
    res[2]=-A[4]*B[1] + A[5]*B[0];
    res[3]=A[2]*B[2] - A[3]*B[1];
    res[4]=-A[1]*B[2] + A[3]*B[0];
    res[5]=A[1]*B[1] - A[2]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> Rotor<T>::conjugate(const R130B2Sp1<T> &A) const {
    R130B2<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[2], 2);
    T x2 = std::pow((*this)[3], 2);
    T x3 = std::pow((*this)[4], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = std::pow((*this)[6], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*A[2];
    T x9 = (*this)[5]*x8;
    T x10 = (*this)[7]*x8;
    T x11 = 2.0*A[1];
    T x12 = (*this)[4]*(*this)[5];
    T x13 = (*this)[4]*x8;
    T x14 = (*this)[0]*(*this)[6];
    T x15 = (*this)[1]*(*this)[2];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*(*this)[7];
    T x18 = 2.0*A[0];
    T x19 = (*this)[0]*x11;
    T x20 = (*this)[1]*x11;
    T x21 = (*this)[6]*x18;
    T x22 = (*this)[6]*x11;
    T x23 = (*this)[0]*x18;
    T x24 = (*this)[1]*x18;
    T x25 = (*this)[2]*x11;
    T x26 = (*this)[2]*x18;
    T x27 = (*this)[0]*x8;
    T x28 = (*this)[6]*x8;
    res[0]=(*this)[0]*x9 - (*this)[1]*x16 + (*this)[2]*x10 + (*this)[6]*x13 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7 + x11*x12 - x11*x14 - x11*x15 - x11*x17;
    res[1]=-(*this)[0]*x13 - (*this)[1]*x10 - (*this)[2]*x16 + (*this)[6]*x9 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 - A[1]*x6 - A[1]*x7 + x12*x18 + x14*x18 - x15*x18 + x17*x18;
    res[2]=-(*this)[3]*x24 - (*this)[3]*x25 + (*this)[4]*x19 + (*this)[4]*x21 + (*this)[5]*x22 - (*this)[5]*x23 + (*this)[7]*x20 - (*this)[7]*x26 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7;
    res[3]=-(*this)[1]*x28 - (*this)[2]*x27 - (*this)[3]*x13 + (*this)[3]*x19 + (*this)[3]*x21 - (*this)[4]*x24 - (*this)[4]*x25 - (*this)[5]*x20 + (*this)[5]*x26 - (*this)[7]*x22 + (*this)[7]*x23 + (*this)[7]*x9;
    res[4]=(*this)[1]*x27 - (*this)[2]*x28 + (*this)[3]*x22 - (*this)[3]*x23 - (*this)[3]*x9 - (*this)[4]*x10 + (*this)[4]*x20 - (*this)[4]*x26 - (*this)[5]*x24 - (*this)[5]*x25 + (*this)[7]*x19 + (*this)[7]*x21;
    res[5]=(*this)[0]*x10 + (*this)[1]*x13 - (*this)[1]*x19 - (*this)[1]*x21 - (*this)[2]*x22 + (*this)[2]*x23 + (*this)[2]*x9 - (*this)[3]*(*this)[4]*x18 - (*this)[3]*(*this)[5]*x11 + (*this)[4]*(*this)[7]*x11 - (*this)[5]*(*this)[7]*x18 - (*this)[6]*x16;
    return res;
};

//----------------------------------
// (Rotor, R130B3) binary operations
//----------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=B[0];
    res[12]=B[1];
    res[13]=B[2];
    res[14]=B[3];
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=-B[0];
    res[12]=-B[1];
    res[13]=-B[2];
    res[14]=-B[3];
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[4]*B[1] - A[5]*B[2] - A[6]*B[3] + A[7]*B[0];
    res[2]=-A[2]*B[3] + A[3]*B[2] + A[4]*B[0] - A[7]*B[1];
    res[3]=A[1]*B[3] - A[3]*B[1] + A[5]*B[0] - A[7]*B[2];
    res[4]=-A[1]*B[2] + A[2]*B[1] + A[6]*B[0] - A[7]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=A[0]*B[1] - A[1]*B[0] - A[5]*B[3] + A[6]*B[2];
    res[13]=A[0]*B[2] - A[2]*B[0] + A[4]*B[3] - A[6]*B[1];
    res[14]=A[0]*B[3] - A[3]*B[0] - A[4]*B[2] + A[5]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[4]*B[1] - A[5]*B[2] - A[6]*B[3];
    res[2]=-A[2]*B[3] + A[3]*B[2] + A[4]*B[0];
    res[3]=A[1]*B[3] - A[3]*B[1] + A[5]*B[0];
    res[4]=-A[1]*B[2] + A[2]*B[1] + A[6]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=A[0]*B[1];
    res[13]=A[0]*B[2];
    res[14]=A[0]*B[3];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B3<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[7]*B[0];
    res[2]=-A[7]*B[1];
    res[3]=-A[7]*B[2];
    res[4]=-A[7]*B[3];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    res[12]=-A[1]*B[0] - A[5]*B[3] + A[6]*B[2];
    res[13]=-A[2]*B[0] + A[4]*B[3] - A[6]*B[1];
    res[14]=-A[3]*B[0] - A[4]*B[2] + A[5]*B[1];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> Rotor<T>::conjugate(const R130B3<T> &A) const {
    R130B3<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[1], 2);
    T x2 = std::pow((*this)[2], 2);
    T x3 = std::pow((*this)[3], 2);
    T x4 = std::pow((*this)[4], 2);
    T x5 = std::pow((*this)[5], 2);
    T x6 = std::pow((*this)[6], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*(*this)[0];
    T x9 = (*this)[1]*x8;
    T x10 = (*this)[2]*A[2];
    T x11 = A[3]*x8;
    T x12 = 2.0*(*this)[5];
    T x13 = (*this)[1]*A[3];
    T x14 = 2.0*A[1];
    T x15 = (*this)[2]*(*this)[6];
    T x16 = 2.0*A[2];
    T x17 = (*this)[3]*x16;
    T x18 = (*this)[4]*(*this)[7];
    T x19 = A[2]*x12;
    T x20 = 2.0*A[3];
    T x21 = (*this)[6]*x20;
    T x22 = (*this)[1]*(*this)[6];
    T x23 = (*this)[2]*x20;
    T x24 = (*this)[3]*x12;
    T x25 = 2.0*x10;
    T x26 = 2.0*x13;
    T x27 = 2.0*A[0];
    T x28 = (*this)[6]*x8;
    T x29 = A[0]*x8;
    T x30 = (*this)[1]*x14;
    T x31 = A[0]*x12;
    T x32 = (*this)[4]*x27;
    T x33 = (*this)[7]*x14;
    res[0]=(*this)[3]*x11 + (*this)[4]*x17 - (*this)[4]*x23 + (*this)[7]*x19 + (*this)[7]*x21 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[1]*x24 + A[1]*x9 + x10*x8 + x12*x13 + x14*x15 + x14*x18 - x16*x22;
    res[1]=(*this)[1]*x25 + (*this)[3]*x26 + (*this)[4]*x19 + (*this)[4]*x21 + (*this)[5]*x11 + (*this)[7]*x17 - (*this)[7]*x23 + A[0]*x24 + A[0]*x9 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 - A[1]*x5 - A[1]*x6 + A[1]*x7 - A[2]*x28 - x15*x27 + x18*x27;
    res[2]=(*this)[2]*x29 + (*this)[2]*x30 + (*this)[3]*x23 - (*this)[3]*x32 - (*this)[3]*x33 + (*this)[4]*A[1]*x12 - (*this)[4]*x11 + (*this)[6]*A[3]*x12 + (*this)[7]*x26 + (*this)[7]*x31 + A[1]*x28 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 + x22*x27;
    res[3]=-(*this)[1]*(*this)[7]*x16 - (*this)[1]*x31 + (*this)[2]*x32 + (*this)[2]*x33 + (*this)[3]*x25 + (*this)[3]*x29 + (*this)[3]*x30 + (*this)[4]*(*this)[6]*x14 + (*this)[4]*A[2]*x8 - (*this)[5]*A[1]*x8 + (*this)[6]*(*this)[7]*x27 + (*this)[6]*x19 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 - A[3]*x4 - A[3]*x5 + A[3]*x6 + A[3]*x7;
    return res;
};

//-------------------------------------
// (Rotor, R130B3Sm1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=B[0];
    res[13]=B[1];
    res[14]=B[2];
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=0;
    res[12]=-B[0];
    res[13]=-B[1];
    res[14]=-B[2];
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[4]*B[0] - A[5]*B[1] - A[6]*B[2];
    res[2]=-A[2]*B[2] + A[3]*B[1] - A[7]*B[0];
    res[3]=A[1]*B[2] - A[3]*B[0] - A[7]*B[1];
    res[4]=-A[1]*B[1] + A[2]*B[0] - A[7]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[12]=A[0]*B[0] - A[5]*B[2] + A[6]*B[1];
    res[13]=A[0]*B[1] + A[4]*B[2] - A[6]*B[0];
    res[14]=A[0]*B[2] - A[4]*B[1] + A[5]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=-A[4]*B[0] - A[5]*B[1] - A[6]*B[2];
    res[2]=-A[2]*B[2] + A[3]*B[1];
    res[3]=A[1]*B[2] - A[3]*B[0];
    res[4]=-A[1]*B[1] + A[2]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=0;
    res[12]=A[0]*B[0];
    res[13]=A[0]*B[1];
    res[14]=A[0]*B[2];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B3Sm1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=-A[7]*B[0];
    res[3]=-A[7]*B[1];
    res[4]=-A[7]*B[2];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=-A[1]*B[0] - A[2]*B[1] - A[3]*B[2];
    res[12]=-A[5]*B[2] + A[6]*B[1];
    res[13]=A[4]*B[2] - A[6]*B[0];
    res[14]=-A[4]*B[1] + A[5]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> Rotor<T>::conjugate(const R130B3Sm1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*(*this)[0];
    T x1 = A[0]*x0;
    T x2 = (*this)[2]*A[1];
    T x3 = A[2]*x0;
    T x4 = 2.0*(*this)[5];
    T x5 = (*this)[1]*A[2];
    T x6 = 2.0*A[0];
    T x7 = (*this)[2]*x6;
    T x8 = 2.0*A[1];
    T x9 = (*this)[3]*x8;
    T x10 = (*this)[7]*x6;
    T x11 = A[1]*x4;
    T x12 = 2.0*A[2];
    T x13 = (*this)[6]*x12;
    T x14 = (*this)[1]*x8;
    T x15 = (*this)[2]*x12;
    T x16 = A[0]*x4;
    T x17 = std::pow((*this)[0], 2);
    T x18 = std::pow((*this)[1], 2);
    T x19 = std::pow((*this)[4], 2);
    T x20 = std::pow((*this)[7], 2);
    T x21 = std::pow((*this)[2], 2);
    T x22 = std::pow((*this)[3], 2);
    T x23 = std::pow((*this)[5], 2);
    T x24 = std::pow((*this)[6], 2);
    T x25 = 2.0*x2;
    T x26 = 2.0*x5;
    T x27 = A[1]*x0;
    res[0]=(*this)[1]*x1 - (*this)[3]*x16 + (*this)[3]*x3 + (*this)[4]*x10 - (*this)[4]*x15 + (*this)[4]*x9 - (*this)[6]*x14 + (*this)[6]*x7 + (*this)[7]*x11 + (*this)[7]*x13 + x0*x2 + x4*x5;
    res[1]=(*this)[1]*x25 + (*this)[3]*x26 + (*this)[4]*x11 + (*this)[4]*x13 + (*this)[5]*x3 - (*this)[6]*x27 - (*this)[7]*x15 + (*this)[7]*x9 + A[0]*x17 + A[0]*x18 + A[0]*x19 + A[0]*x20 - A[0]*x21 - A[0]*x22 - A[0]*x23 - A[0]*x24;
    res[2]=(*this)[1]*x7 - (*this)[3]*x10 + (*this)[3]*x15 + (*this)[4]*x16 - (*this)[4]*x3 + (*this)[6]*A[2]*x4 + (*this)[6]*x1 + (*this)[7]*x26 + A[1]*x17 - A[1]*x18 - A[1]*x19 + A[1]*x20 + A[1]*x21 - A[1]*x22 + A[1]*x23 - A[1]*x24;
    res[3]=(*this)[1]*(*this)[3]*x6 + (*this)[3]*x25 + (*this)[4]*(*this)[6]*x6 + (*this)[4]*x27 - (*this)[5]*x1 + (*this)[6]*x11 - (*this)[7]*x14 + (*this)[7]*x7 + A[2]*x17 - A[2]*x18 - A[2]*x19 + A[2]*x20 - A[2]*x21 + A[2]*x22 - A[2]*x23 + A[2]*x24;
    return res;
};

//-------------------------------------
// (Rotor, R130B3Sp1) binary operations
//-------------------------------------

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=A[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1];
    res[6]=A[2];
    res[7]=A[3];
    res[8]=A[4];
    res[9]=A[5];
    res[10]=A[6];
    res[11]=-B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=A[7];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[7]*B[0];
    res[2]=A[4]*B[0];
    res[3]=A[5]*B[0];
    res[4]=A[6]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=-A[1]*B[0];
    res[13]=-A[2]*B[0];
    res[14]=-A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=0;
    res[2]=A[4]*B[0];
    res[3]=A[5]*B[0];
    res[4]=A[6]*B[0];
    res[5]=0;
    res[6]=0;
    res[7]=0;
    res[8]=0;
    res[9]=0;
    res[10]=0;
    res[11]=A[0]*B[0];
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B3Sp1<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[7]*B[0];
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
    res[12]=-A[1]*B[0];
    res[13]=-A[2]*B[0];
    res[14]=-A[3]*B[0];
    res[15]=0;
    return res;
};

template<typename T>
inline R130B3<T> Rotor<T>::conjugate(const R130B3Sp1<T> &A) const {
    R130B3<T> res;
    T x0 = 2.0*A[0];
    res[0]=A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) + std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) + std::pow((*this)[7], 2));
    res[1]=x0*((*this)[0]*(*this)[1] - (*this)[2]*(*this)[6] + (*this)[3]*(*this)[5] + (*this)[4]*(*this)[7]);
    res[2]=x0*((*this)[0]*(*this)[2] + (*this)[1]*(*this)[6] - (*this)[3]*(*this)[4] + (*this)[5]*(*this)[7]);
    res[3]=x0*((*this)[0]*(*this)[3] - (*this)[1]*(*this)[5] + (*this)[2]*(*this)[4] + (*this)[6]*(*this)[7]);
    return res;
};

//----------------------------------
// (Rotor, R130B4) binary operations
//----------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const R130B4<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7] + B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const R130B4<T> &B) {
    Rotor<T> res;
    res[0]=A[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4];
    res[5]=A[5];
    res[6]=A[6];
    res[7]=A[7] - B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const R130B4<T> &B) {
    Rotor<T> res;
    res[0]=-A[7]*B[0];
    res[1]=-A[4]*B[0];
    res[2]=-A[5]*B[0];
    res[3]=-A[6]*B[0];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    res[6]=A[3]*B[0];
    res[7]=A[0]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const R130B4<T> &B) {
    Rotor<T> res;
    res[0]=-A[7]*B[0];
    res[1]=-A[4]*B[0];
    res[2]=-A[5]*B[0];
    res[3]=-A[6]*B[0];
    res[4]=A[1]*B[0];
    res[5]=A[2]*B[0];
    res[6]=A[3]*B[0];
    res[7]=A[0]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130B4<T> &B) {
    R130MV<T> res;
    res[0]=0;
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
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> Rotor<T>::conjugate(const R130B4<T> &A) const {
    R130MV<T> res;
    res[0]=2.0*A[0]*(-(*this)[0]*(*this)[7] + (*this)[1]*(*this)[4] + (*this)[2]*(*this)[5] + (*this)[3]*(*this)[6]);
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
    res[15]=A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2) + std::pow((*this)[4], 2) + std::pow((*this)[5], 2) + std::pow((*this)[6], 2) - std::pow((*this)[7], 2));
    return res;
};

//----------------------------------
// (Rotor, R130MV) binary operations
//----------------------------------

template<typename T>
inline Rotor<T>::operator R130MV<T>() const {
    R130MV<T> res;
    res[0]=(*this)[0];
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[1];
    res[6]=(*this)[2];
    res[7]=(*this)[3];
    res[8]=(*this)[4];
    res[9]=(*this)[5];
    res[10]=(*this)[6];
    res[11]=0;
    res[12]=0;
    res[13]=0;
    res[14]=0;
    res[15]=(*this)[7];
    return res;
};

template<typename T>
inline R130MV<T> operator+(const Rotor<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] + B[0];
    res[1]=B[1];
    res[2]=B[2];
    res[3]=B[3];
    res[4]=B[4];
    res[5]=A[1] + B[5];
    res[6]=A[2] + B[6];
    res[7]=A[3] + B[7];
    res[8]=A[4] + B[8];
    res[9]=A[5] + B[9];
    res[10]=A[6] + B[10];
    res[11]=B[11];
    res[12]=B[12];
    res[13]=B[13];
    res[14]=B[14];
    res[15]=A[7] + B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator-(const Rotor<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0] - B[0];
    res[1]=-B[1];
    res[2]=-B[2];
    res[3]=-B[3];
    res[4]=-B[4];
    res[5]=A[1] - B[5];
    res[6]=A[2] - B[6];
    res[7]=A[3] - B[7];
    res[8]=A[4] - B[8];
    res[9]=A[5] - B[9];
    res[10]=A[6] - B[10];
    res[11]=-B[11];
    res[12]=-B[12];
    res[13]=-B[13];
    res[14]=-B[14];
    res[15]=A[7] - B[15];
    return res;
};

template<typename T>
inline R130MV<T> operator*(const Rotor<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7] - A[4]*B[8] - A[5]*B[9] - A[6]*B[10] - A[7]*B[15];
    res[1]=A[0]*B[1] + A[1]*B[2] + A[2]*B[3] + A[3]*B[4] - A[4]*B[12] - A[5]*B[13] - A[6]*B[14] + A[7]*B[11];
    res[2]=A[0]*B[2] + A[1]*B[1] - A[2]*B[14] + A[3]*B[13] + A[4]*B[11] - A[5]*B[4] + A[6]*B[3] - A[7]*B[12];
    res[3]=A[0]*B[3] + A[1]*B[14] + A[2]*B[1] - A[3]*B[12] + A[4]*B[4] + A[5]*B[11] - A[6]*B[2] - A[7]*B[13];
    res[4]=A[0]*B[4] - A[1]*B[13] + A[2]*B[12] + A[3]*B[1] - A[4]*B[3] + A[5]*B[2] + A[6]*B[11] - A[7]*B[14];
    res[5]=A[0]*B[5] + A[1]*B[0] - A[2]*B[10] + A[3]*B[9] - A[4]*B[15] - A[5]*B[7] + A[6]*B[6] - A[7]*B[8];
    res[6]=A[0]*B[6] + A[1]*B[10] + A[2]*B[0] - A[3]*B[8] + A[4]*B[7] - A[5]*B[15] - A[6]*B[5] - A[7]*B[9];
    res[7]=A[0]*B[7] - A[1]*B[9] + A[2]*B[8] + A[3]*B[0] - A[4]*B[6] + A[5]*B[5] - A[6]*B[15] - A[7]*B[10];
    res[8]=A[0]*B[8] + A[1]*B[15] + A[2]*B[7] - A[3]*B[6] + A[4]*B[0] - A[5]*B[10] + A[6]*B[9] + A[7]*B[5];
    res[9]=A[0]*B[9] - A[1]*B[7] + A[2]*B[15] + A[3]*B[5] + A[4]*B[10] + A[5]*B[0] - A[6]*B[8] + A[7]*B[6];
    res[10]=A[0]*B[10] + A[1]*B[6] - A[2]*B[5] + A[3]*B[15] - A[4]*B[9] + A[5]*B[8] + A[6]*B[0] + A[7]*B[7];
    res[11]=A[0]*B[11] - A[1]*B[12] - A[2]*B[13] - A[3]*B[14] - A[4]*B[2] - A[5]*B[3] - A[6]*B[4] - A[7]*B[1];
    res[12]=A[0]*B[12] - A[1]*B[11] + A[2]*B[4] - A[3]*B[3] + A[4]*B[1] - A[5]*B[14] + A[6]*B[13] + A[7]*B[2];
    res[13]=A[0]*B[13] - A[1]*B[4] - A[2]*B[11] + A[3]*B[2] + A[4]*B[14] + A[5]*B[1] - A[6]*B[12] + A[7]*B[3];
    res[14]=A[0]*B[14] + A[1]*B[3] - A[2]*B[2] - A[3]*B[11] - A[4]*B[13] + A[5]*B[12] + A[6]*B[1] + A[7]*B[4];
    res[15]=A[0]*B[15] + A[1]*B[8] + A[2]*B[9] + A[3]*B[10] + A[4]*B[5] + A[5]*B[6] + A[6]*B[7] + A[7]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator|(const Rotor<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=A[0]*B[0] + A[1]*B[5] + A[2]*B[6] + A[3]*B[7] - A[4]*B[8] - A[5]*B[9] - A[6]*B[10] - A[7]*B[15];
    res[1]=A[0]*B[1] - A[4]*B[12] - A[5]*B[13] - A[6]*B[14];
    res[2]=A[0]*B[2] - A[2]*B[14] + A[3]*B[13] + A[4]*B[11];
    res[3]=A[0]*B[3] + A[1]*B[14] - A[3]*B[12] + A[5]*B[11];
    res[4]=A[0]*B[4] - A[1]*B[13] + A[2]*B[12] + A[6]*B[11];
    res[5]=A[0]*B[5] + A[1]*B[0] - A[4]*B[15] - A[7]*B[8];
    res[6]=A[0]*B[6] + A[2]*B[0] - A[5]*B[15] - A[7]*B[9];
    res[7]=A[0]*B[7] + A[3]*B[0] - A[6]*B[15] - A[7]*B[10];
    res[8]=A[0]*B[8] + A[1]*B[15] + A[4]*B[0] + A[7]*B[5];
    res[9]=A[0]*B[9] + A[2]*B[15] + A[5]*B[0] + A[7]*B[6];
    res[10]=A[0]*B[10] + A[3]*B[15] + A[6]*B[0] + A[7]*B[7];
    res[11]=A[0]*B[11] - A[4]*B[2] - A[5]*B[3] - A[6]*B[4];
    res[12]=A[0]*B[12] + A[2]*B[4] - A[3]*B[3] + A[4]*B[1];
    res[13]=A[0]*B[13] - A[1]*B[4] + A[3]*B[2] + A[5]*B[1];
    res[14]=A[0]*B[14] + A[1]*B[3] - A[2]*B[2] + A[6]*B[1];
    res[15]=A[0]*B[15] + A[1]*B[8] + A[2]*B[9] + A[3]*B[10] + A[4]*B[5] + A[5]*B[6] + A[6]*B[7] + A[7]*B[0];
    return res;
};

template<typename T>
inline R130MV<T> operator^(const Rotor<T> &A, const R130MV<T> &B) {
    R130MV<T> res;
    res[0]=0;
    res[1]=A[1]*B[2] + A[2]*B[3] + A[3]*B[4] + A[7]*B[11];
    res[2]=A[1]*B[1] - A[5]*B[4] + A[6]*B[3] - A[7]*B[12];
    res[3]=A[2]*B[1] + A[4]*B[4] - A[6]*B[2] - A[7]*B[13];
    res[4]=A[3]*B[1] - A[4]*B[3] + A[5]*B[2] - A[7]*B[14];
    res[5]=-A[2]*B[10] + A[3]*B[9] - A[5]*B[7] + A[6]*B[6];
    res[6]=A[1]*B[10] - A[3]*B[8] + A[4]*B[7] - A[6]*B[5];
    res[7]=-A[1]*B[9] + A[2]*B[8] - A[4]*B[6] + A[5]*B[5];
    res[8]=A[2]*B[7] - A[3]*B[6] - A[5]*B[10] + A[6]*B[9];
    res[9]=-A[1]*B[7] + A[3]*B[5] + A[4]*B[10] - A[6]*B[8];
    res[10]=A[1]*B[6] - A[2]*B[5] - A[4]*B[9] + A[5]*B[8];
    res[11]=-A[1]*B[12] - A[2]*B[13] - A[3]*B[14] - A[7]*B[1];
    res[12]=-A[1]*B[11] - A[5]*B[14] + A[6]*B[13] + A[7]*B[2];
    res[13]=-A[2]*B[11] + A[4]*B[14] - A[6]*B[12] + A[7]*B[3];
    res[14]=-A[3]*B[11] - A[4]*B[13] + A[5]*B[12] + A[7]*B[4];
    res[15]=0;
    return res;
};

template<typename T>
inline R130MV<T> Rotor<T>::conjugate(const R130MV<T> &A) const {
    R130MV<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[6], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = std::pow((*this)[3], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*A[15];
    T x9 = (*this)[1]*(*this)[4];
    T x10 = (*this)[2]*(*this)[5];
    T x11 = (*this)[3]*(*this)[6];
    T x12 = (*this)[0]*(*this)[7];
    T x13 = 2.0*(*this)[1];
    T x14 = A[3]*x13;
    T x15 = 2.0*(*this)[2];
    T x16 = (*this)[4]*A[4];
    T x17 = 2.0*(*this)[3];
    T x18 = (*this)[5]*x17;
    T x19 = (*this)[0]*x13;
    T x20 = (*this)[0]*A[3];
    T x21 = (*this)[0]*A[4];
    T x22 = A[4]*x13;
    T x23 = (*this)[6]*x15;
    T x24 = A[3]*x17;
    T x25 = 2.0*(*this)[7];
    T x26 = (*this)[4]*x25;
    T x27 = (*this)[5]*A[3];
    T x28 = (*this)[6]*A[4];
    T x29 = 2.0*(*this)[5];
    T x30 = 2.0*(*this)[4];
    T x31 = 2.0*(*this)[6];
    T x32 = A[4]*x15;
    T x33 = (*this)[0]*A[2];
    T x34 = A[2]*x13;
    T x35 = A[1]*x17;
    T x36 = (*this)[4]*x29;
    T x37 = A[1]*x15;
    T x38 = A[1]*x13;
    T x39 = (*this)[7]*A[2];
    T x40 = A[1]*x25;
    T x41 = (*this)[3]*x15;
    T x42 = (*this)[6]*x30;
    T x43 = (*this)[0]*x15;
    T x44 = A[7]*x29;
    T x45 = 2.0*A[8];
    T x46 = (*this)[5]*x13;
    T x47 = (*this)[6]*x13;
    T x48 = (*this)[4]*x15;
    T x49 = (*this)[7]*A[7];
    T x50 = (*this)[4]*x17;
    T x51 = (*this)[6]*x25;
    T x52 = (*this)[0]*x17;
    T x53 = (*this)[0]*x31;
    T x54 = (*this)[2]*x13;
    T x55 = (*this)[3]*x13;
    T x56 = (*this)[7]*x17;
    T x57 = (*this)[5]*x25;
    T x58 = 2.0*A[9];
    T x59 = (*this)[0]*x30;
    T x60 = (*this)[7]*x13;
    T x61 = 2.0*A[10];
    T x62 = (*this)[6]*x29;
    T x63 = (*this)[0]*x29;
    T x64 = (*this)[7]*x15;
    T x65 = 2.0*A[5];
    T x66 = 2.0*A[6];
    T x67 = 2.0*A[7];
    T x68 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7 + x10*x8 + x11*x8 - x12*x8 + x8*x9;
    res[1]=-(*this)[4]*x24 - (*this)[5]*x22 + (*this)[6]*x14 + A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 + A[1]*x4 + A[1]*x5 + A[1]*x6 + A[1]*x7 + A[2]*x18 - A[2]*x19 - A[2]*x23 - A[2]*x26 + x15*x16 - x15*x20 - x17*x21 - x25*x27 - x25*x28;
    res[2]=(*this)[2]*x14 + (*this)[3]*x22 + (*this)[7]*x24 - (*this)[7]*x32 - A[1]*x18 - A[1]*x19 + A[1]*x23 - A[1]*x26 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + x16*x31 - x20*x31 + x21*x29 + x27*x30;
    res[3]=-2.0*(*this)[0]*x16 - (*this)[0]*x37 + (*this)[2]*x34 + (*this)[3]*x32 + (*this)[4]*x35 - (*this)[5]*x40 - (*this)[6]*x38 + (*this)[7]*x22 + A[2]*x36 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 - A[3]*x6 + A[3]*x7 - x17*x39 + x28*x29 + x31*x33;
    res[4]=-(*this)[0]*x35 + (*this)[3]*x34 - (*this)[4]*x37 + (*this)[5]*x38 - (*this)[6]*x40 - (*this)[7]*x14 + A[2]*x42 + A[3]*x41 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3 - A[4]*x4 - A[4]*x5 + A[4]*x6 + A[4]*x7 + x15*x39 + x20*x30 + x27*x31 - x29*x33;
    res[5]=(*this)[0]*x44 + A[10]*x43 + A[10]*x47 + A[10]*x50 - A[10]*x57 + A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 - A[5]*x4 + A[5]*x5 + A[5]*x6 - A[5]*x7 + A[6]*x36 - A[6]*x53 - A[6]*x54 - A[6]*x56 + A[7]*x42 - A[7]*x55 + A[9]*x46 + A[9]*x48 + A[9]*x51 - A[9]*x52 - x10*x45 - x11*x45 - x12*x45 + x15*x49 + x45*x9;
    res[6]=(*this)[6]*x44 + A[10]*x18 - A[10]*x19 + A[10]*x23 + A[10]*x26 + A[5]*x36 + A[5]*x53 - A[5]*x54 + A[5]*x56 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3 + A[6]*x4 - A[6]*x5 + A[6]*x6 - A[6]*x7 - A[7]*x41 - A[7]*x59 + A[8]*x46 + A[8]*x48 - A[8]*x51 + A[8]*x52 + x10*x58 - x11*x58 - x12*x58 - x13*x49 - x58*x9;
    res[7]=A[5]*x42 - A[5]*x55 - A[5]*x63 - A[5]*x64 - A[6]*x41 + A[6]*x59 + A[6]*x60 + A[6]*x62 + A[7]*x0 - A[7]*x1 - A[7]*x2 + A[7]*x3 + A[7]*x4 + A[7]*x5 - A[7]*x6 - A[7]*x7 - A[8]*x43 + A[8]*x47 + A[8]*x50 + A[8]*x57 + A[9]*x18 + A[9]*x19 + A[9]*x23 - A[9]*x26 - x10*x61 + x11*x61 - x12*x61 - x61*x9;
    res[8]=A[10]*x42 - A[10]*x55 + A[10]*x63 + A[10]*x64 - A[6]*x46 - A[6]*x48 - A[6]*x51 + A[6]*x52 - A[7]*x43 - A[7]*x47 - A[7]*x50 + A[7]*x57 + A[8]*x0 + A[8]*x1 - A[8]*x2 - A[8]*x3 - A[8]*x4 + A[8]*x5 + A[8]*x6 - A[8]*x7 + A[9]*x36 - A[9]*x53 - A[9]*x54 - A[9]*x56 + x10*x65 + x11*x65 + x12*x65 - x65*x9;
    res[9]=-A[10]*x41 - A[10]*x59 - A[10]*x60 + A[10]*x62 - A[5]*x46 - A[5]*x48 + A[5]*x51 - A[5]*x52 - A[7]*x18 + A[7]*x19 - A[7]*x23 - A[7]*x26 + A[8]*x36 + A[8]*x53 - A[8]*x54 + A[8]*x56 + A[9]*x0 - A[9]*x1 + A[9]*x2 - A[9]*x3 + A[9]*x4 - A[9]*x5 + A[9]*x6 - A[9]*x7 - x10*x66 + x11*x66 + x12*x66 + x66*x9;
    res[10]=A[10]*x0 - A[10]*x1 - A[10]*x2 + A[10]*x3 + A[10]*x4 + A[10]*x5 - A[10]*x6 - A[10]*x7 + A[5]*x43 - A[5]*x47 - A[5]*x50 - A[5]*x57 - A[6]*x18 - A[6]*x19 - A[6]*x23 + A[6]*x26 + A[8]*x42 - A[8]*x55 - A[8]*x63 - A[8]*x64 - A[9]*x41 + A[9]*x59 + A[9]*x60 + A[9]*x62 + x10*x67 - x11*x67 + x12*x67 + x67*x9;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[11]*x4 + A[11]*x5 + A[11]*x6 + A[11]*x7 - A[12]*x18 + A[12]*x19 + A[12]*x23 + A[12]*x26 + A[13]*x43 - A[13]*x47 + A[13]*x50 + A[13]*x57 + A[14]*x46 - A[14]*x48 + A[14]*x51 + A[14]*x52;
    res[12]=A[11]*x18 + A[11]*x19 - A[11]*x23 + A[11]*x26 + A[12]*x0 + A[12]*x1 - A[12]*x2 - A[12]*x3 + A[12]*x4 - A[12]*x5 - A[12]*x6 + A[12]*x7 + A[13]*x36 - A[13]*x53 + A[13]*x54 + A[13]*x56 + A[14]*x42 + A[14]*x55 + A[14]*x63 - A[14]*x64;
    res[13]=A[11]*x43 + A[11]*x47 - A[11]*x50 + A[11]*x57 + A[12]*x36 + A[12]*x53 + A[12]*x54 - A[12]*x56 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 - A[13]*x4 + A[13]*x5 - A[13]*x6 + A[13]*x7 + A[14]*x41 - A[14]*x59 + A[14]*x60 + A[14]*x62;
    res[14]=-A[11]*x46 + A[11]*x48 + A[11]*x51 + A[11]*x52 + A[12]*x42 + A[12]*x55 - A[12]*x63 + A[12]*x64 + A[13]*x41 + A[13]*x59 - A[13]*x60 + A[13]*x62 + A[14]*x0 - A[14]*x1 - A[14]*x2 + A[14]*x3 - A[14]*x4 - A[14]*x5 + A[14]*x6 + A[14]*x7;
    res[15]=A[15]*x0 + A[15]*x1 + A[15]*x2 + A[15]*x3 - A[15]*x4 - A[15]*x5 - A[15]*x6 - A[15]*x7 - x10*x68 - x11*x68 + x12*x68 - x68*x9;
    return res;
};

//------------------------------------
// (Rotor, Rotation) binary operations
//------------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4] + B[1];
    res[5]=A[5] + B[2];
    res[6]=A[6] + B[3];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1];
    res[2]=A[2];
    res[3]=A[3];
    res[4]=A[4] - B[1];
    res[5]=A[5] - B[2];
    res[6]=A[6] - B[3];
    res[7]=A[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] - A[4]*B[1] - A[5]*B[2] - A[6]*B[3];
    res[1]=A[1]*B[0] - A[2]*B[3] + A[3]*B[2] - A[7]*B[1];
    res[2]=A[1]*B[3] + A[2]*B[0] - A[3]*B[1] - A[7]*B[2];
    res[3]=-A[1]*B[2] + A[2]*B[1] + A[3]*B[0] - A[7]*B[3];
    res[4]=A[0]*B[1] + A[4]*B[0] - A[5]*B[3] + A[6]*B[2];
    res[5]=A[0]*B[2] + A[4]*B[3] + A[5]*B[0] - A[6]*B[1];
    res[6]=A[0]*B[3] - A[4]*B[2] + A[5]*B[1] + A[6]*B[0];
    res[7]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[7]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const Rotation<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] - A[4]*B[1] - A[5]*B[2] - A[6]*B[3];
    res[1]=A[1]*B[0] - A[7]*B[1];
    res[2]=A[2]*B[0] - A[7]*B[2];
    res[3]=A[3]*B[0] - A[7]*B[3];
    res[4]=A[0]*B[1] + A[4]*B[0];
    res[5]=A[0]*B[2] + A[5]*B[0];
    res[6]=A[0]*B[3] + A[6]*B[0];
    res[7]=A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[7]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotor<T> &A, const Rotation<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[3] + A[3]*B[2];
    res[1]=A[1]*B[3] - A[3]*B[1];
    res[2]=-A[1]*B[2] + A[2]*B[1];
    res[3]=-A[5]*B[3] + A[6]*B[2];
    res[4]=A[4]*B[3] - A[6]*B[1];
    res[5]=-A[4]*B[2] + A[5]*B[1];
    return res;
};

template<typename T>
inline Rotor<T> Rotor<T>::conjugate(const Rotation<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[6], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = std::pow((*this)[3], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*A[3];
    T x9 = (*this)[2]*x8;
    T x10 = 2.0*(*this)[1];
    T x11 = A[1]*x10;
    T x12 = (*this)[5]*A[2];
    T x13 = (*this)[1]*x8;
    T x14 = 2.0*A[2];
    T x15 = (*this)[2]*(*this)[4];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[7]*x14;
    T x18 = (*this)[0]*(*this)[3];
    T x19 = 2.0*A[1];
    T x20 = (*this)[0]*(*this)[7];
    T x21 = (*this)[5]*x19;
    T x22 = (*this)[3]*(*this)[6];
    T x23 = (*this)[7]*x8;
    T x24 = 2.0*x12;
    T x25 = A[2]*x10;
    T x26 = (*this)[6]*x19;
    T x27 = (*this)[6]*x14;
    T x28 = (*this)[3]*x19;
    T x29 = (*this)[2]*x19;
    T x30 = (*this)[5]*x8;
    T x31 = (*this)[4]*x8;
    T x32 = 2.0*A[0];
    res[0]=A[0]*(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7);
    res[1]=(*this)[0]*x9 - (*this)[2]*x21 + (*this)[4]*x11 + (*this)[4]*x16 - (*this)[5]*x23 + (*this)[6]*x13 + (*this)[6]*x17 + x10*x12 + x14*x15 - x14*x18 - x19*x20 - x19*x22;
    res[2]=-(*this)[0]*x13 - (*this)[0]*x17 + (*this)[2]*x24 + (*this)[4]*x23 - (*this)[4]*x25 + (*this)[5]*x11 + (*this)[5]*x16 + (*this)[6]*x9 - (*this)[7]*x26 - x14*x22 + x15*x19 + x18*x19;
    res[3]=(*this)[0]*x25 - (*this)[0]*x29 + (*this)[2]*x27 + (*this)[3]*x24 - (*this)[4]*x13 - (*this)[4]*x17 + (*this)[4]*x28 - (*this)[5]*x9 + (*this)[6]*x11 + (*this)[6]*x16 + (*this)[7]*x21 - x20*x8;
    res[4]=-(*this)[0]*x27 + (*this)[0]*x30 - (*this)[2]*x25 - (*this)[3]*x13 - (*this)[3]*x17 + (*this)[4]*x24 + (*this)[6]*x31 + (*this)[7]*x9 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7;
    res[5]=(*this)[0]*x26 - (*this)[0]*x31 - (*this)[2]*x11 - (*this)[3]*x9 + (*this)[4]*x21 + (*this)[6]*x30 - (*this)[7]*x13 + (*this)[7]*x28 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7;
    res[6]=(*this)[0]*(*this)[4]*x14 - (*this)[0]*x21 - (*this)[2]*(*this)[3]*x14 - (*this)[3]*x11 + (*this)[4]*x26 + (*this)[6]*x24 + (*this)[7]*x25 - (*this)[7]*x29 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7;
    res[7]=-(*this)[2]*(*this)[5]*x32 - (*this)[4]*A[0]*x10 + x20*x32 - x22*x32;
    return res;
};

//---------------------------------
// (Rotor, Rotor) binary operations
//---------------------------------

template<typename T>
inline Rotor<T> operator+(const Rotor<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] + B[0];
    res[1]=A[1] + B[1];
    res[2]=A[2] + B[2];
    res[3]=A[3] + B[3];
    res[4]=A[4] + B[4];
    res[5]=A[5] + B[5];
    res[6]=A[6] + B[6];
    res[7]=A[7] + B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator-(const Rotor<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0] - B[0];
    res[1]=A[1] - B[1];
    res[2]=A[2] - B[2];
    res[3]=A[3] - B[3];
    res[4]=A[4] - B[4];
    res[5]=A[5] - B[5];
    res[6]=A[6] - B[6];
    res[7]=A[7] - B[7];
    return res;
};

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] - A[4]*B[4] - A[5]*B[5] - A[6]*B[6] - A[7]*B[7];
    res[1]=A[0]*B[1] + A[1]*B[0] - A[2]*B[6] + A[3]*B[5] - A[4]*B[7] - A[5]*B[3] + A[6]*B[2] - A[7]*B[4];
    res[2]=A[0]*B[2] + A[1]*B[6] + A[2]*B[0] - A[3]*B[4] + A[4]*B[3] - A[5]*B[7] - A[6]*B[1] - A[7]*B[5];
    res[3]=A[0]*B[3] - A[1]*B[5] + A[2]*B[4] + A[3]*B[0] - A[4]*B[2] + A[5]*B[1] - A[6]*B[7] - A[7]*B[6];
    res[4]=A[0]*B[4] + A[1]*B[7] + A[2]*B[3] - A[3]*B[2] + A[4]*B[0] - A[5]*B[6] + A[6]*B[5] + A[7]*B[1];
    res[5]=A[0]*B[5] - A[1]*B[3] + A[2]*B[7] + A[3]*B[1] + A[4]*B[6] + A[5]*B[0] - A[6]*B[4] + A[7]*B[2];
    res[6]=A[0]*B[6] + A[1]*B[2] - A[2]*B[1] + A[3]*B[7] - A[4]*B[5] + A[5]*B[4] + A[6]*B[0] + A[7]*B[3];
    res[7]=A[0]*B[7] + A[1]*B[4] + A[2]*B[5] + A[3]*B[6] + A[4]*B[1] + A[5]*B[2] + A[6]*B[3] + A[7]*B[0];
    return res;
};

template<typename T>
inline Rotor<T> operator|(const Rotor<T> &A, const Rotor<T> &B) {
    Rotor<T> res;
    res[0]=A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] - A[4]*B[4] - A[5]*B[5] - A[6]*B[6] - A[7]*B[7];
    res[1]=A[0]*B[1] + A[1]*B[0] - A[4]*B[7] - A[7]*B[4];
    res[2]=A[0]*B[2] + A[2]*B[0] - A[5]*B[7] - A[7]*B[5];
    res[3]=A[0]*B[3] + A[3]*B[0] - A[6]*B[7] - A[7]*B[6];
    res[4]=A[0]*B[4] + A[1]*B[7] + A[4]*B[0] + A[7]*B[1];
    res[5]=A[0]*B[5] + A[2]*B[7] + A[5]*B[0] + A[7]*B[2];
    res[6]=A[0]*B[6] + A[3]*B[7] + A[6]*B[0] + A[7]*B[3];
    res[7]=A[0]*B[7] + A[1]*B[4] + A[2]*B[5] + A[3]*B[6] + A[4]*B[1] + A[5]*B[2] + A[6]*B[3] + A[7]*B[0];
    return res;
};

template<typename T>
inline R130B2<T> operator^(const Rotor<T> &A, const Rotor<T> &B) {
    R130B2<T> res;
    res[0]=-A[2]*B[6] + A[3]*B[5] - A[5]*B[3] + A[6]*B[2];
    res[1]=A[1]*B[6] - A[3]*B[4] + A[4]*B[3] - A[6]*B[1];
    res[2]=-A[1]*B[5] + A[2]*B[4] - A[4]*B[2] + A[5]*B[1];
    res[3]=A[2]*B[3] - A[3]*B[2] - A[5]*B[6] + A[6]*B[5];
    res[4]=-A[1]*B[3] + A[3]*B[1] + A[4]*B[6] - A[6]*B[4];
    res[5]=A[1]*B[2] - A[2]*B[1] - A[4]*B[5] + A[5]*B[4];
    return res;
};

template<typename T>
inline Rotor<T> Rotor<T>::conjugate(const Rotor<T> &A) const {
    Rotor<T> res;
    T x0 = std::pow((*this)[0], 2);
    T x1 = std::pow((*this)[4], 2);
    T x2 = std::pow((*this)[5], 2);
    T x3 = std::pow((*this)[6], 2);
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = std::pow((*this)[3], 2);
    T x7 = std::pow((*this)[7], 2);
    T x8 = 2.0*A[7];
    T x9 = (*this)[1]*(*this)[4];
    T x10 = (*this)[2]*(*this)[5];
    T x11 = (*this)[3]*(*this)[6];
    T x12 = (*this)[0]*(*this)[7];
    T x13 = 2.0*(*this)[0];
    T x14 = (*this)[2]*A[6];
    T x15 = (*this)[5]*A[3];
    T x16 = 2.0*A[4];
    T x17 = 2.0*(*this)[1];
    T x18 = (*this)[5]*A[5];
    T x19 = (*this)[6]*A[6];
    T x20 = 2.0*(*this)[2];
    T x21 = (*this)[4]*A[5];
    T x22 = (*this)[7]*A[3];
    T x23 = 2.0*(*this)[4];
    T x24 = (*this)[3]*A[6];
    T x25 = (*this)[5]*x23;
    T x26 = (*this)[6]*x23;
    T x27 = 2.0*(*this)[7];
    T x28 = (*this)[6]*A[5];
    T x29 = (*this)[3]*x13;
    T x30 = (*this)[6]*x13;
    T x31 = (*this)[2]*x17;
    T x32 = (*this)[3]*A[3];
    T x33 = (*this)[3]*x27;
    T x34 = (*this)[5]*A[6];
    T x35 = (*this)[1]*x16;
    T x36 = (*this)[4]*x16;
    T x37 = 2.0*A[5];
    T x38 = 2.0*(*this)[6];
    T x39 = 2.0*(*this)[5];
    T x40 = (*this)[7]*A[6];
    T x41 = (*this)[1]*x13;
    T x42 = (*this)[4]*x13;
    T x43 = (*this)[7]*x16;
    T x44 = A[2]*x17;
    T x45 = 2.0*(*this)[3];
    T x46 = 2.0*A[6];
    T x47 = (*this)[2]*x13;
    T x48 = (*this)[5]*A[1];
    T x49 = A[1]*x17;
    T x50 = A[2]*x20;
    T x51 = A[1]*x20;
    T x52 = 2.0*A[1];
    T x53 = (*this)[6]*A[3];
    T x54 = (*this)[6]*x27;
    T x55 = 2.0*A[2];
    T x56 = 2.0*A[3];
    T x57 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7 + x10*x8 + x11*x8 - x12*x8 + x8*x9;
    res[1]=A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 + A[2]*x25 - A[2]*x30 - A[2]*x31 - A[2]*x33 + A[3]*x26 - A[5]*x29 - x10*x16 - x11*x16 - x12*x16 + x13*x14 + x13*x15 + x16*x9 + x17*x18 + x17*x19 - x17*x32 + x20*x21 + x20*x22 + x23*x24 + x27*x28 - x27*x34;
    res[2]=(*this)[2]*x36 + (*this)[5]*x35 - (*this)[6]*x43 + A[1]*x25 + A[1]*x30 - A[1]*x31 + A[1]*x33 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7 - A[3]*x42 + A[4]*x29 - A[6]*x41 + x10*x37 - x11*x37 - x12*x37 + x14*x38 + x15*x38 - x17*x22 - x20*x32 + x23*x40 + x24*x39 - x37*x9;
    res[3]=(*this)[3]*x36 - (*this)[3]*x49 - (*this)[3]*x50 + (*this)[5]*A[2]*x38 + (*this)[5]*x43 + (*this)[6]*x35 + (*this)[7]*x44 - (*this)[7]*x51 + A[1]*x26 + A[2]*x42 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7 - A[4]*x47 + A[5]*x41 - x10*x46 + x11*x46 - x12*x46 - x13*x48 + x18*x45 + x20*x28 - x21*x27 - x46*x9;
    res[4]=-(*this)[4]*x50 - (*this)[5]*x44 + A[2]*x29 - A[2]*x54 - A[3]*x47 + A[4]*x0 + A[4]*x1 - A[4]*x2 - A[4]*x3 - A[4]*x4 + A[4]*x5 + A[4]*x6 - A[4]*x7 - A[5]*x31 - A[5]*x33 + x10*x52 + x11*x52 + x12*x52 - x13*x28 + x13*x34 + x14*x27 + x15*x27 - x17*x24 - x17*x53 + x18*x23 + x19*x23 - x23*x32 - x52*x9;
    res[5]=-(*this)[2]*x35 + (*this)[3]*x43 - (*this)[4]*x51 + (*this)[5]*x36 - A[1]*x29 + A[1]*x54 + A[3]*x41 + A[4]*x30 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[5]*x3 + A[5]*x4 - A[5]*x5 + A[5]*x6 - A[5]*x7 - A[6]*x42 - x10*x55 + x11*x55 + x12*x55 - x14*x45 - x15*x45 - x17*x40 - x17*x48 + x19*x39 - x20*x53 - x22*x23 + x55*x9;
    res[6]=-(*this)[2]*x43 - (*this)[3]*A[1]*x23 - (*this)[3]*A[2]*x39 - (*this)[3]*A[5]*x20 - (*this)[3]*x35 - (*this)[5]*A[4]*x13 + (*this)[6]*x36 - (*this)[6]*x49 - (*this)[6]*x50 + (*this)[7]*A[2]*x23 + (*this)[7]*A[5]*x17 + A[1]*x47 - A[2]*x41 + A[6]*x0 - A[6]*x1 - A[6]*x2 + A[6]*x3 + A[6]*x4 + A[6]*x5 - A[6]*x6 - A[6]*x7 + x10*x56 - x11*x56 + x12*x56 + x13*x21 + x18*x38 - x27*x48 + x56*x9;
    res[7]=A[7]*x0 + A[7]*x1 + A[7]*x2 + A[7]*x3 - A[7]*x4 - A[7]*x5 - A[7]*x6 - A[7]*x7 - x10*x57 - x11*x57 + x12*x57 - x57*x9;
    return res;
};

//----------------------------------
// (Rotor, scalar) binary operations
//----------------------------------

template<typename T>
inline Rotor<T> operator*(const Rotor<T> &A, const T &B) {
    Rotor<T> res;
    res[0]=A[0]*B;
    res[1]=A[1]*B;
    res[2]=A[2]*B;
    res[3]=A[3]*B;
    res[4]=A[4]*B;
    res[5]=A[5]*B;
    res[6]=A[6]*B;
    res[7]=A[7]*B;
    return res;
};
template<typename T>
inline Rotor<T> operator*(const T &A, const Rotor<T> &B) {
    return B*A;
};
template<typename T>
inline Rotor<T> operator/(const Rotor<T> &A, const T &B) {
    return A * (1.0 / B);
};


} // namespace stga3

#endif // LI_STGA3_BinaryOperators_H
