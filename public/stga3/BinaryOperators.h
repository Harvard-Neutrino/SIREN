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
    T x10 = (*this)[0]*A[3];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x5 - x4*x6;
    res[6]=-A[1]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - x6*x7;
    res[7]=A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 - x4*x8 - x7*x9;
    res[8]=-x10*x7 + x11*x9;
    res[9]=x10*x4 - x11*x8;
    res[10]=(*this)[0]*A[1]*x7 - (*this)[0]*A[2]*x4;
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
    T x9 = A[0]*x4;
    T x10 = A[1]*x8;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[1]*x5 - x4*x6 - x4*x7;
    res[1]=-A[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + x6*x8 + x7*x8;
    res[2]=(*this)[2]*x10 + 2.0*(*this)[2]*x7 - (*this)[2]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3;
    res[3]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 - (*this)[3]*x9 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
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
    T x6 = 2.0*(*this)[1];
    T x7 = std::pow((*this)[2], 2);
    T x8 = std::pow((*this)[3], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=-x0*x1 - x0*x2 - x0*x3;
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x7 - A[0]*x8 + x2*x6 + x3*x6;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x7 - A[1]*x8 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x7 + A[2]*x8 + x1*x10 + x10*x2;
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
    T x3 = 2.0*(*this)[0];
    T x4 = A[5]*x3;
    T x5 = (*this)[3]*x3;
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x6;
    T x8 = (*this)[3]*A[2];
    T x9 = std::pow((*this)[1], 2);
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[1]*x3;
    T x12 = (*this)[2]*x3;
    T x13 = (*this)[3]*x6;
    T x14 = (*this)[3]*x10;
    res[0]=(*this)[2]*x4 + A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x9 - A[1]*x7 - A[4]*x5 - x6*x8;
    res[1]=-(*this)[1]*x4 - A[0]*x7 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x9 + A[3]*x5 - x10*x8;
    res[2]=-A[0]*x13 - A[1]*x14 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x9 - A[3]*x12 + A[4]*x11;
    res[3]=A[1]*x5 - A[2]*x12 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x9 - A[4]*x7 - A[5]*x13;
    res[4]=-A[0]*x5 + A[2]*x11 - A[3]*x7 + A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x9 - A[5]*x14;
    res[5]=A[0]*x12 - A[1]*x11 - A[3]*x13 - A[4]*x14 + A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x9;
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
    T x8 = 2.0*(*this)[2];
    T x9 = 2.0*(*this)[1];
    T x10 = (*this)[3]*A[2];
    T x11 = std::pow((*this)[1], 2);
    res[0]=(*this)[2]*x1 - A[1]*x2;
    res[1]=-(*this)[1]*x1 + A[0]*x2;
    res[2]=x0*x3 - x0*x4;
    res[3]=-A[0]*x11 + A[0]*x5 + A[0]*x6 + A[0]*x7 - x10*x9 - x3*x8;
    res[4]=A[1]*x11 + A[1]*x5 - A[1]*x6 + A[1]*x7 - x10*x8 - x4*x9;
    res[5]=-(*this)[3]*A[0]*x9 - (*this)[3]*A[1]*x8 + A[2]*x11 + A[2]*x5 + A[2]*x6 - A[2]*x7;
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
    T x3 = 2.0*(*this)[1];
    T x4 = (*this)[2]*x3;
    T x5 = (*this)[3]*A[2];
    T x6 = std::pow((*this)[1], 2);
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[0];
    T x9 = (*this)[3]*A[1];
    T x10 = (*this)[0]*A[2];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x6 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x6 - x5*x7;
    res[2]=A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x6 - x3*x8 - x7*x9;
    res[3]=-x10*x7 + x11*x9;
    res[4]=x10*x3 - x11*x8;
    res[5]=(*this)[0]*A[0]*x7 - (*this)[0]*A[1]*x3;
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
    T x6 = 2.0*(*this)[1];
    T x7 = std::pow((*this)[2], 2);
    T x8 = std::pow((*this)[3], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=x0*x1 + x0*x2 + x0*x3;
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x7 - A[0]*x8 + x2*x6 + x3*x6;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x7 - A[1]*x8 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x7 + A[2]*x8 + x1*x10 + x10*x2;
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
    T x9 = A[1]*x4;
    T x10 = A[2]*x8;
    T x11 = 2.0*(*this)[2];
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[2]*x8;
    T x15 = (*this)[3]*A[7];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*x11;
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x5 - x4*x6 - x4*x7;
    res[2]=-A[1]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + x6*x8 + x7*x8;
    res[3]=(*this)[2]*x10 - (*this)[2]*x9 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 + x11*x7;
    res[4]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 - (*this)[3]*x9 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3;
    res[5]=A[10]*x12 + A[5]*x0 - A[5]*x1 + A[5]*x2 + A[5]*x3 - A[6]*x14 - A[9]*x13 - x15*x8;
    res[6]=-A[10]*x5 - A[5]*x14 + A[6]*x0 + A[6]*x1 - A[6]*x2 + A[6]*x3 + A[8]*x13 - x11*x15;
    res[7]=-A[5]*x16 - A[6]*x17 + A[7]*x0 + A[7]*x1 + A[7]*x2 - A[7]*x3 - A[8]*x12 + A[9]*x5;
    res[8]=-A[10]*x16 + A[6]*x13 - A[7]*x12 + A[8]*x0 - A[8]*x1 + A[8]*x2 + A[8]*x3 - A[9]*x14;
    res[9]=-A[10]*x17 - A[5]*x13 + A[7]*x5 - A[8]*x14 + A[9]*x0 + A[9]*x1 - A[9]*x2 + A[9]*x3;
    res[10]=A[10]*x0 + A[10]*x1 + A[10]*x2 - A[10]*x3 + A[5]*x12 - A[6]*x5 - A[8]*x16 - A[9]*x17;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[12]*x5 + A[13]*x12 + A[14]*x13;
    res[12]=A[11]*x5 + A[12]*x0 + A[12]*x1 - A[12]*x2 - A[12]*x3 + A[13]*x14 + A[14]*x16;
    res[13]=A[11]*x12 + A[12]*x14 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 + A[14]*x17;
    res[14]=A[11]*x13 + A[12]*x16 + A[13]*x17 + A[14]*x0 - A[14]*x1 - A[14]*x2 + A[14]*x3;
    res[15]=A[15]*x0 - A[15]*x1 - A[15]*x2 - A[15]*x3;
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
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[2]*x5 - A[2]*x6;
    res[6]=-(*this)[1]*x5 + A[1]*x6;
    res[7]=x4*x7 - x4*x8;
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
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=(*this)[2]*x5 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x8 - A[5]*x6 - x7*x9;
    res[2]=-(*this)[1]*x5 - A[1]*x8 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 + A[4]*x6 - x10*x9;
    res[3]=-A[1]*x13 - A[2]*x14 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 - A[4]*x12 + A[5]*x11;
    res[4]=A[2]*x6 - A[3]*x12 + A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3 - A[5]*x8 - A[6]*x13;
    res[5]=-A[1]*x6 + A[3]*x11 - A[4]*x8 + A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x3 - A[6]*x14;
    res[6]=A[1]*x12 - A[2]*x11 - A[4]*x13 - A[5]*x14 + A[6]*x0 + A[6]*x1 + A[6]*x2 - A[6]*x3;
    res[7]=A[7]*x0 - A[7]*x1 - A[7]*x2 - A[7]*x3;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x5 - x4*x6;
    res[6]=-A[1]*x5 + A[2]*x0 - A[2]*x1 + A[2]*x2 + A[2]*x3 - x6*x7;
    res[7]=A[3]*x0 + A[3]*x1 - A[3]*x2 + A[3]*x3 - x4*x8 - x7*x9;
    res[8]=x10*x7 - x11*x9;
    res[9]=-x10*x4 + x11*x8;
    res[10]=-(*this)[0]*A[1]*x7 + (*this)[0]*A[2]*x4;
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
    res[0]=-A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2));
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
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*x0;
    T x2 = (*this)[2]*A[2];
    T x3 = (*this)[3]*A[3];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = A[0]*x0;
    T x10 = A[1]*x8;
    res[0]=-A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7 + A[1]*x1 + x0*x2 + x0*x3;
    res[1]=-A[0]*x1 + A[1]*x4 + A[1]*x5 - A[1]*x6 - A[1]*x7 + x2*x8 + x3*x8;
    res[2]=(*this)[2]*x10 + 2.0*(*this)[2]*x3 - (*this)[2]*x9 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7;
    res[3]=(*this)[3]*x10 + 2.0*(*this)[3]*x2 - (*this)[3]*x9 + A[3]*x4 - A[3]*x5 - A[3]*x6 + A[3]*x7;
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
    T x6 = 2.0*(*this)[1];
    T x7 = std::pow((*this)[2], 2);
    T x8 = std::pow((*this)[3], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=x0*x1 + x0*x2 + x0*x3;
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x7 - A[0]*x8 + x2*x6 + x3*x6;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x7 - A[1]*x8 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x7 + A[2]*x8 + x1*x10 + x10*x2;
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
    T x3 = 2.0*(*this)[0];
    T x4 = A[5]*x3;
    T x5 = (*this)[3]*x3;
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x6;
    T x8 = (*this)[3]*A[2];
    T x9 = std::pow((*this)[1], 2);
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[1]*x3;
    T x12 = (*this)[2]*x3;
    T x13 = (*this)[3]*x6;
    T x14 = (*this)[3]*x10;
    res[0]=(*this)[2]*x4 + A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x9 - A[1]*x7 - A[4]*x5 - x6*x8;
    res[1]=-(*this)[1]*x4 - A[0]*x7 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x9 + A[3]*x5 - x10*x8;
    res[2]=-A[0]*x13 - A[1]*x14 + A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x9 - A[3]*x12 + A[4]*x11;
    res[3]=-A[1]*x5 + A[2]*x12 - A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x9 + A[4]*x7 + A[5]*x13;
    res[4]=A[0]*x5 - A[2]*x11 + A[3]*x7 - A[4]*x0 + A[4]*x1 - A[4]*x2 - A[4]*x9 + A[5]*x14;
    res[5]=-A[0]*x12 + A[1]*x11 + A[3]*x13 + A[4]*x14 - A[5]*x0 - A[5]*x1 + A[5]*x2 - A[5]*x9;
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
    T x6 = 2.0*(*this)[2];
    T x7 = 2.0*(*this)[1];
    T x8 = (*this)[3]*A[2];
    T x9 = std::pow((*this)[0], 2);
    T x10 = std::pow((*this)[2], 2);
    T x11 = std::pow((*this)[3], 2);
    res[0]=(*this)[2]*x1 - A[1]*x2;
    res[1]=-(*this)[1]*x1 + A[0]*x2;
    res[2]=x0*x3 - x0*x4;
    res[3]=-A[0]*x10 - A[0]*x11 + A[0]*x5 - A[0]*x9 + x3*x6 + x7*x8;
    res[4]=A[1]*x10 - A[1]*x11 - A[1]*x5 - A[1]*x9 + x4*x7 + x6*x8;
    res[5]=(*this)[3]*A[0]*x7 + (*this)[3]*A[1]*x6 - A[2]*x10 + A[2]*x11 - A[2]*x5 - A[2]*x9;
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
    T x3 = 2.0*(*this)[1];
    T x4 = (*this)[2]*x3;
    T x5 = (*this)[3]*A[2];
    T x6 = std::pow((*this)[1], 2);
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[0];
    T x9 = (*this)[3]*A[1];
    T x10 = (*this)[0]*A[2];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x6 - A[1]*x4 - x3*x5;
    res[1]=-A[0]*x4 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x6 - x5*x7;
    res[2]=A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x6 - x3*x8 - x7*x9;
    res[3]=x10*x7 - x11*x9;
    res[4]=-x10*x3 + x11*x8;
    res[5]=-(*this)[0]*A[0]*x7 + (*this)[0]*A[1]*x3;
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
    T x6 = 2.0*(*this)[1];
    T x7 = std::pow((*this)[0], 2);
    T x8 = std::pow((*this)[1], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=x0*x1 + x0*x2 + x0*x3;
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x7 - A[0]*x8 - x2*x6 - x3*x6;
    res[2]=-A[1]*x4 + A[1]*x5 - A[1]*x7 + A[1]*x8 - x1*x9 - x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x7 + A[2]*x8 - x1*x10 - x10*x2;
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
    T x9 = A[1]*x4;
    T x10 = A[2]*x8;
    T x11 = 2.0*(*this)[2];
    T x12 = (*this)[2]*x4;
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[2]*x8;
    T x15 = (*this)[3]*A[7];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*x11;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3;
    res[1]=-A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x5 + x4*x6 + x4*x7;
    res[2]=-A[1]*x5 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x3 + x6*x8 + x7*x8;
    res[3]=(*this)[2]*x10 - (*this)[2]*x9 - A[3]*x0 + A[3]*x1 - A[3]*x2 + A[3]*x3 + x11*x7;
    res[4]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 - (*this)[3]*x9 - A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3;
    res[5]=A[10]*x12 - A[5]*x0 + A[5]*x1 + A[5]*x2 + A[5]*x3 - A[6]*x14 - A[9]*x13 - x15*x8;
    res[6]=-A[10]*x5 - A[5]*x14 + A[6]*x0 - A[6]*x1 + A[6]*x2 + A[6]*x3 + A[8]*x13 - x11*x15;
    res[7]=-A[5]*x16 - A[6]*x17 + A[7]*x0 + A[7]*x1 - A[7]*x2 + A[7]*x3 - A[8]*x12 + A[9]*x5;
    res[8]=A[10]*x16 - A[6]*x13 + A[7]*x12 + A[8]*x0 - A[8]*x1 - A[8]*x2 - A[8]*x3 + A[9]*x14;
    res[9]=A[10]*x17 + A[5]*x13 - A[7]*x5 + A[8]*x14 - A[9]*x0 + A[9]*x1 - A[9]*x2 - A[9]*x3;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 - A[10]*x3 - A[5]*x12 + A[6]*x5 + A[8]*x16 + A[9]*x17;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[12]*x5 + A[13]*x12 + A[14]*x13;
    res[12]=-A[11]*x5 - A[12]*x0 + A[12]*x1 + A[12]*x2 - A[12]*x3 - A[13]*x14 - A[14]*x16;
    res[13]=-A[11]*x12 - A[12]*x14 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 - A[14]*x17;
    res[14]=-A[11]*x13 - A[12]*x16 - A[13]*x17 + A[14]*x0 + A[14]*x1 - A[14]*x2 - A[14]*x3;
    res[15]=-A[15]*x0 - A[15]*x1 - A[15]*x2 + A[15]*x3;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[2]*x5 - A[2]*x6;
    res[6]=-(*this)[1]*x5 + A[1]*x6;
    res[7]=x4*x7 - x4*x8;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3;
    res[1]=(*this)[2]*x5 - A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3 - A[2]*x8 - A[5]*x6 - x7*x9;
    res[2]=-(*this)[1]*x5 - A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 + A[2]*x3 + A[4]*x6 - x10*x9;
    res[3]=-A[1]*x13 - A[2]*x14 + A[3]*x0 + A[3]*x1 - A[3]*x2 + A[3]*x3 - A[4]*x12 + A[5]*x11;
    res[4]=-A[2]*x6 + A[3]*x12 + A[4]*x0 - A[4]*x1 - A[4]*x2 - A[4]*x3 + A[5]*x8 + A[6]*x13;
    res[5]=A[1]*x6 - A[3]*x11 + A[4]*x8 - A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[6]*x14;
    res[6]=-A[1]*x12 + A[2]*x11 + A[4]*x13 + A[5]*x14 - A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3;
    res[7]=-A[7]*x0 - A[7]*x1 - A[7]*x2 + A[7]*x3;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x2;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-A[0]*x7 - A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x5;
    res[3]=-A[3]*x0 - A[3]*x1 + A[3]*x5 + A[4]*x3 + A[5]*x7;
    res[4]=A[3]*x3 + A[4]*x0 - A[4]*x1 - A[4]*x5 + A[5]*x8;
    res[5]=A[3]*x7 + A[4]*x8 - A[5]*x0 + A[5]*x1 - A[5]*x5;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-(*this)[2]*A[0]*x2 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x5;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-(*this)[2]*A[0]*x2 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x5;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
    res[1]=-A[1]*x0 - A[1]*x1 - A[1]*x2;
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=-A[5]*x0 + A[5]*x1 + A[5]*x2 - A[6]*x4 - A[7]*x7;
    res[6]=-A[5]*x4 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[7]*x8;
    res[7]=-A[5]*x7 - A[6]*x8 + A[7]*x0 + A[7]*x1 - A[7]*x2;
    res[8]=A[10]*x7 + A[8]*x0 - A[8]*x1 - A[8]*x2 + A[9]*x4;
    res[9]=A[10]*x8 + A[8]*x4 - A[9]*x0 + A[9]*x1 - A[9]*x2;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 + A[8]*x7 + A[9]*x8;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2;
    res[12]=-A[12]*x0 + A[12]*x1 + A[12]*x2 - A[13]*x4 - A[14]*x7;
    res[13]=-A[12]*x4 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[14]*x8;
    res[14]=-A[12]*x7 - A[13]*x8 + A[14]*x0 + A[14]*x1 - A[14]*x2;
    res[15]=-A[15]*x0 - A[15]*x1 - A[15]*x2;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-A[1]*x7 - A[2]*x8 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    res[4]=A[4]*x0 - A[4]*x1 - A[4]*x2 + A[5]*x4 + A[6]*x7;
    res[5]=A[4]*x4 - A[5]*x0 + A[5]*x1 - A[5]*x2 + A[6]*x8;
    res[6]=A[4]*x7 + A[5]*x8 - A[6]*x0 - A[6]*x1 + A[6]*x2;
    res[7]=-A[7]*x0 - A[7]*x1 - A[7]*x2;
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
    T x6 = 2.0*(*this)[0];
    T x7 = (*this)[1]*A[2];
    T x8 = (*this)[2]*A[3];
    T x9 = 2.0*(*this)[3];
    T x10 = (*this)[4]*x9;
    T x11 = (*this)[5]*A[3];
    T x12 = A[1]*x6;
    T x13 = 2.0*(*this)[1];
    T x14 = 2.0*(*this)[4];
    T x15 = 2.0*(*this)[2];
    T x16 = (*this)[5]*A[1];
    T x17 = (*this)[5]*A[2];
    T x18 = (*this)[0]*x9;
    T x19 = (*this)[1]*A[1];
    T x20 = (*this)[1]*x14;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[0]*x4 - A[0]*x5;
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[2]*x10 + x11*x9 - x6*x7 - x6*x8;
    res[2]=-(*this)[1]*x12 + A[1]*x10 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 + x11*x14 - x13*x8;
    res[3]=-(*this)[2]*x12 - A[3]*x0 - A[3]*x1 + A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 + x14*x17 - x15*x7 + x16*x9;
    res[4]=-(*this)[4]*A[2]*x6 - A[1]*x18 - x11*x6 + x14*x19 + x15*x16 - x7*x9 - x8*x9;
    res[5]=-(*this)[4]*x12 + A[2]*x18 - x11*x13 - x14*x7 - x14*x8 + x15*x17 - x19*x9;
    res[6]=-(*this)[2]*A[1]*x9 - (*this)[2]*A[2]*x14 - (*this)[5]*x12 - 2.0*(*this)[5]*x7 + A[3]*x18 + A[3]*x20 - x11*x15;
    res[7]=-(*this)[5]*A[0]*x15 - A[0]*x18 - A[0]*x20;
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
    res[0]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) - std::pow((*this)[3], 2) - std::pow((*this)[4], 2) - std::pow((*this)[5], 2));
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
    T x7 = A[3]*x6;
    T x8 = A[2]*x6;
    T x9 = 2.0*(*this)[1];
    T x10 = (*this)[3]*A[3];
    T x11 = (*this)[5]*x9;
    T x12 = 2.0*(*this)[2];
    T x13 = (*this)[3]*A[2];
    T x14 = (*this)[4]*x12;
    T x15 = 2.0*(*this)[4];
    T x16 = 2.0*(*this)[5];
    T x17 = A[1]*x6;
    T x18 = A[0]*x6;
    T x19 = (*this)[2]*x9;
    T x20 = (*this)[3]*A[0];
    T x21 = (*this)[3]*A[1];
    T x22 = (*this)[5]*x15;
    res[0]=-(*this)[4]*x7 + (*this)[5]*x8 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 - A[1]*x11 + A[1]*x14 + x10*x9 - x12*x13;
    res[1]=(*this)[1]*x8 + (*this)[2]*x7 + A[0]*x11 - A[0]*x14 + A[1]*x0 - A[1]*x1 - A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + x10*x16 + x13*x15;
    res[2]=(*this)[1]*x17 - (*this)[5]*x18 - A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[3]*x19 + A[3]*x22 + x12*x20 + x15*x21;
    res[3]=(*this)[2]*x17 + (*this)[4]*x18 + A[2]*x19 + A[2]*x22 - A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 + x16*x21 - x20*x9;
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
    T x1 = A[2]*x0;
    T x2 = A[1]*x0;
    T x3 = 2.0*(*this)[1];
    T x4 = (*this)[3]*A[2];
    T x5 = (*this)[5]*A[0];
    T x6 = 2.0*(*this)[2];
    T x7 = (*this)[3]*A[1];
    T x8 = (*this)[4]*A[0];
    T x9 = std::pow((*this)[0], 2);
    T x10 = std::pow((*this)[3], 2);
    T x11 = 2.0*(*this)[4];
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[2], 2);
    T x14 = std::pow((*this)[4], 2);
    T x15 = std::pow((*this)[5], 2);
    T x16 = A[0]*x0;
    T x17 = (*this)[2]*x3;
    T x18 = 2.0*(*this)[3];
    T x19 = (*this)[5]*x11;
    res[0]=-(*this)[4]*x1 + (*this)[5]*x2 + x3*x4 - x3*x5 - x6*x7 + x6*x8;
    res[1]=(*this)[1]*x2 + (*this)[2]*x1 + 2.0*(*this)[5]*x4 + A[0]*x10 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x9 + x11*x7;
    res[2]=(*this)[1]*x16 - A[1]*x10 + A[1]*x12 - A[1]*x13 + A[1]*x14 - A[1]*x15 - A[1]*x9 + A[2]*x17 + A[2]*x19 + x18*x8;
    res[3]=(*this)[2]*x16 + A[1]*x17 + A[1]*x19 - A[2]*x10 - A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 - A[2]*x9 + x18*x5;
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
    res[2]=-x0*((*this)[0]*(*this)[5] - (*this)[2]*(*this)[3]);
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
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*x3;
    T x5 = (*this)[2]*A[2];
    T x6 = A[3]*x3;
    T x7 = (*this)[4]*A[4];
    T x8 = (*this)[5]*A[5];
    T x9 = 2.0*(*this)[3];
    T x10 = (*this)[1]*A[4];
    T x11 = 2.0*A[3];
    T x12 = (*this)[1]*(*this)[4];
    T x13 = (*this)[2]*A[5];
    T x14 = (*this)[2]*(*this)[5];
    T x15 = (*this)[4]*x9;
    T x16 = (*this)[5]*A[2];
    T x17 = std::pow((*this)[0], 2);
    T x18 = std::pow((*this)[4], 2);
    T x19 = std::pow((*this)[5], 2);
    T x20 = (*this)[3]*x3;
    T x21 = 2.0*(*this)[1];
    T x22 = A[3]*x9;
    T x23 = 2.0*(*this)[4];
    T x24 = 2.0*x14;
    T x25 = A[0]*x3;
    T x26 = (*this)[2]*A[1];
    T x27 = 2.0*x12;
    T x28 = 2.0*(*this)[5];
    T x29 = 2.0*(*this)[2];
    T x30 = A[0]*x9;
    T x31 = (*this)[5]*A[1];
    res[0]=(*this)[3]*x6 + A[0]*x0 + A[0]*x1 - A[0]*x17 - A[0]*x18 - A[0]*x19 + A[0]*x2 + A[1]*x15 - A[1]*x4 + x10*x9 - x11*x12 - x11*x14 + x13*x9 + x16*x9 - x3*x5 + x3*x7 + x3*x8;
    res[1]=(*this)[1]*x22 + (*this)[4]*x6 + A[0]*x15 - A[0]*x4 - A[1]*x0 + A[1]*x1 + A[1]*x17 + A[1]*x18 - A[1]*x19 - A[1]*x2 - A[4]*x20 - A[4]*x24 + x13*x23 + x16*x23 - x21*x5 + x21*x7 + x21*x8;
    res[2]=(*this)[2]*x22 - (*this)[2]*x25 + (*this)[5]*x30 + (*this)[5]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x17 - A[2]*x18 + A[2]*x19 - A[2]*x2 - A[5]*x20 - A[5]*x27 + x10*x28 - x21*x26 + x23*x31 + x29*x7 + x29*x8;
    res[3]=-(*this)[1]*A[1]*x9 - (*this)[4]*A[1]*x3 - A[0]*x20 + A[0]*x24 + A[0]*x27 + A[3]*x0 + A[3]*x1 - A[3]*x17 - A[3]*x18 - A[3]*x19 + A[3]*x2 - x10*x3 - x13*x3 - x16*x3 - x5*x9 + x7*x9 + x8*x9;
    res[4]=-(*this)[1]*x30 - (*this)[1]*x6 - (*this)[4]*x25 + A[1]*x20 + A[1]*x24 - A[1]*x27 + A[3]*x15 - A[4]*x0 + A[4]*x1 + A[4]*x17 + A[4]*x18 - A[4]*x19 - A[4]*x2 - x13*x21 - x16*x21 - x23*x5 + x23*x8;
    res[5]=-(*this)[2]*x30 - (*this)[2]*x6 + (*this)[5]*x22 - (*this)[5]*x25 + A[2]*x20 + A[2]*x27 + A[5]*x0 - A[5]*x1 + A[5]*x17 - A[5]*x18 + A[5]*x19 - A[5]*x2 - x10*x29 - x16*x29 - x21*x31 - x23*x26 + x28*x7;
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
    T x6 = 2.0*A[0];
    T x7 = (*this)[1]*(*this)[4];
    T x8 = (*this)[2]*A[2];
    T x9 = (*this)[2]*(*this)[5];
    T x10 = (*this)[3]*x0;
    T x11 = A[0]*x4;
    T x12 = 2.0*(*this)[1];
    T x13 = 2.0*(*this)[4];
    T x14 = 2.0*(*this)[5];
    T x15 = 2.0*(*this)[2];
    T x16 = std::pow((*this)[1], 2);
    T x17 = std::pow((*this)[2], 2);
    T x18 = std::pow((*this)[3], 2);
    T x19 = std::pow((*this)[0], 2);
    T x20 = std::pow((*this)[4], 2);
    T x21 = std::pow((*this)[5], 2);
    res[0]=(*this)[3]*x1 + x0*x2 + x0*x3 + x4*x5 + x4*x8 - x6*x7 - x6*x9;
    res[1]=(*this)[1]*x11 + (*this)[4]*x1 - A[1]*x10 - 2.0*A[1]*x9 + x12*x2 + x12*x3 + x13*x8;
    res[2]=(*this)[2]*x11 + (*this)[5]*x1 - A[2]*x10 - 2.0*A[2]*x7 + x14*x5 + x15*x2 + x15*x3;
    res[3]=A[0]*x16 + A[0]*x17 + A[0]*x18 - A[0]*x19 - A[0]*x20 - A[0]*x21 - x0*x5 - x0*x8 + x2*x4 + x3*x4;
    res[4]=-(*this)[1]*x1 + (*this)[4]*x11 - A[1]*x16 + A[1]*x17 - A[1]*x18 + A[1]*x19 + A[1]*x20 - A[1]*x21 - x12*x8 + x13*x3;
    res[5]=-(*this)[2]*x1 + (*this)[5]*x11 + A[2]*x16 - A[2]*x17 - A[2]*x18 + A[2]*x19 - A[2]*x20 + A[2]*x21 + x14*x2 - x15*x5;
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
    T x3 = 2.0*(*this)[0];
    T x4 = (*this)[1]*A[1];
    T x5 = (*this)[2]*A[2];
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[4]*x6;
    T x8 = (*this)[5]*A[2];
    T x9 = std::pow((*this)[0], 2);
    T x10 = std::pow((*this)[4], 2);
    T x11 = std::pow((*this)[5], 2);
    T x12 = A[0]*x3;
    T x13 = 2.0*(*this)[1];
    T x14 = 2.0*(*this)[4];
    T x15 = 2.0*(*this)[2];
    T x16 = (*this)[5]*A[0];
    T x17 = (*this)[5]*A[1];
    T x18 = (*this)[0]*x6;
    T x19 = (*this)[1]*A[0];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 + A[0]*x2 - A[0]*x9 + A[1]*x7 - x3*x4 - x3*x5 + x6*x8;
    res[1]=-(*this)[1]*x12 + A[0]*x7 - A[1]*x0 + A[1]*x1 + A[1]*x10 - A[1]*x11 - A[1]*x2 + A[1]*x9 - x13*x5 + x14*x8;
    res[2]=-(*this)[2]*x12 + A[2]*x0 - A[2]*x1 - A[2]*x10 + A[2]*x11 - A[2]*x2 + A[2]*x9 + x14*x17 - x15*x4 + x16*x6;
    res[3]=-(*this)[4]*A[1]*x3 - A[0]*x18 + x14*x19 + x15*x16 - x3*x8 - x4*x6 - x5*x6;
    res[4]=-(*this)[4]*x12 + A[1]*x18 - x13*x8 - x14*x4 - x14*x5 + x15*x17 - x19*x6;
    res[5]=(*this)[1]*A[2]*x14 - (*this)[2]*A[0]*x6 - (*this)[2]*A[1]*x14 - (*this)[5]*x12 - 2.0*(*this)[5]*x4 + A[2]*x18 - x15*x8;
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
    T x8 = A[2]*x6;
    T x9 = 2.0*(*this)[1];
    T x10 = (*this)[3]*A[3];
    T x11 = (*this)[5]*x9;
    T x12 = 2.0*(*this)[2];
    T x13 = (*this)[3]*A[2];
    T x14 = (*this)[4]*x12;
    T x15 = 2.0*(*this)[4];
    T x16 = 2.0*(*this)[5];
    T x17 = A[1]*x6;
    T x18 = A[0]*x6;
    T x19 = (*this)[2]*x9;
    T x20 = (*this)[3]*A[0];
    T x21 = (*this)[3]*A[1];
    T x22 = (*this)[5]*x15;
    res[0]=(*this)[4]*x7 - (*this)[5]*x8 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[1]*x11 - A[1]*x14 - x10*x9 + x12*x13;
    res[1]=(*this)[1]*x8 + (*this)[2]*x7 - A[0]*x11 + A[0]*x14 + A[1]*x0 - A[1]*x1 - A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + x10*x16 + x13*x15;
    res[2]=(*this)[1]*x17 + (*this)[5]*x18 - A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[3]*x19 + A[3]*x22 - x12*x20 + x15*x21;
    res[3]=(*this)[2]*x17 - (*this)[4]*x18 + A[2]*x19 + A[2]*x22 - A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 + x16*x21 + x20*x9;
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
    T x2 = A[1]*x0;
    T x3 = 2.0*(*this)[1];
    T x4 = (*this)[3]*A[2];
    T x5 = (*this)[5]*A[0];
    T x6 = 2.0*(*this)[2];
    T x7 = (*this)[3]*A[1];
    T x8 = (*this)[4]*A[0];
    T x9 = std::pow((*this)[0], 2);
    T x10 = std::pow((*this)[3], 2);
    T x11 = 2.0*(*this)[4];
    T x12 = std::pow((*this)[1], 2);
    T x13 = std::pow((*this)[2], 2);
    T x14 = std::pow((*this)[4], 2);
    T x15 = std::pow((*this)[5], 2);
    T x16 = A[0]*x0;
    T x17 = (*this)[2]*x3;
    T x18 = 2.0*(*this)[3];
    T x19 = (*this)[5]*x11;
    res[0]=(*this)[4]*x1 - (*this)[5]*x2 - x3*x4 + x3*x5 + x6*x7 - x6*x8;
    res[1]=(*this)[1]*x2 + (*this)[2]*x1 + 2.0*(*this)[5]*x4 + A[0]*x10 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x9 + x11*x7;
    res[2]=(*this)[1]*x16 - A[1]*x10 + A[1]*x12 - A[1]*x13 + A[1]*x14 - A[1]*x15 - A[1]*x9 + A[2]*x17 + A[2]*x19 + x18*x8;
    res[3]=(*this)[2]*x16 + A[1]*x17 + A[1]*x19 - A[2]*x10 - A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 - A[2]*x9 + x18*x5;
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
    res[1]=-x0*((*this)[1]*(*this)[5] - (*this)[2]*(*this)[4]);
    res[2]=x0*((*this)[0]*(*this)[5] - (*this)[2]*(*this)[3]);
    res[3]=-x0*((*this)[0]*(*this)[4] - (*this)[1]*(*this)[3]);
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
    res[15]=-A[0]*(std::pow((*this)[0], 2) + std::pow((*this)[1], 2) + std::pow((*this)[2], 2) - std::pow((*this)[3], 2) - std::pow((*this)[4], 2) - std::pow((*this)[5], 2));
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
    T x3 = 2.0*A[15];
    T x4 = (*this)[0]*(*this)[3];
    T x5 = (*this)[1]*(*this)[4];
    T x6 = (*this)[2]*(*this)[5];
    T x7 = std::pow((*this)[0], 2);
    T x8 = std::pow((*this)[1], 2);
    T x9 = std::pow((*this)[2], 2);
    T x10 = 2.0*(*this)[0];
    T x11 = A[4]*x10;
    T x12 = A[3]*x10;
    T x13 = 2.0*(*this)[1];
    T x14 = (*this)[3]*A[4];
    T x15 = (*this)[5]*x13;
    T x16 = 2.0*(*this)[2];
    T x17 = (*this)[3]*A[3];
    T x18 = (*this)[4]*x16;
    T x19 = 2.0*(*this)[4];
    T x20 = 2.0*(*this)[5];
    T x21 = A[2]*x10;
    T x22 = A[1]*x10;
    T x23 = (*this)[2]*x13;
    T x24 = (*this)[3]*A[1];
    T x25 = (*this)[3]*A[2];
    T x26 = (*this)[5]*x19;
    T x27 = (*this)[1]*x10;
    T x28 = (*this)[2]*x10;
    T x29 = 2.0*A[8];
    T x30 = (*this)[4]*x10;
    T x31 = (*this)[5]*x10;
    T x32 = (*this)[3]*x13;
    T x33 = (*this)[3]*x16;
    T x34 = (*this)[3]*x19;
    T x35 = (*this)[3]*x20;
    T x36 = 2.0*A[9];
    T x37 = 2.0*A[10];
    T x38 = 2.0*A[5];
    T x39 = 2.0*A[6];
    T x40 = 2.0*A[7];
    T x41 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x7 - A[0]*x8 - A[0]*x9 + x3*x4 + x3*x5 + x3*x6;
    res[1]=-(*this)[4]*x11 + (*this)[5]*x12 + A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x7 + A[1]*x8 + A[1]*x9 - A[2]*x15 + A[2]*x18 + x13*x14 - x16*x17;
    res[2]=(*this)[1]*x12 + (*this)[2]*x11 + A[1]*x15 - A[1]*x18 + A[2]*x0 - A[2]*x1 - A[2]*x2 + A[2]*x7 - A[2]*x8 - A[2]*x9 + x14*x20 + x17*x19;
    res[3]=(*this)[1]*x21 - (*this)[5]*x22 - A[3]*x0 + A[3]*x1 - A[3]*x2 - A[3]*x7 + A[3]*x8 - A[3]*x9 + A[4]*x23 + A[4]*x26 + x16*x24 + x19*x25;
    res[4]=(*this)[2]*x21 + (*this)[4]*x22 + A[3]*x23 + A[3]*x26 - A[4]*x0 - A[4]*x1 + A[4]*x2 - A[4]*x7 - A[4]*x8 + A[4]*x9 - x13*x24 + x20*x25;
    res[5]=A[10]*x31 + A[10]*x33 + A[5]*x0 - A[5]*x1 - A[5]*x2 - A[5]*x7 + A[5]*x8 + A[5]*x9 - A[6]*x27 + A[6]*x34 - A[7]*x28 + A[7]*x35 + A[9]*x30 + A[9]*x32 + x29*x4 - x29*x5 - x29*x6;
    res[6]=A[10]*x15 + A[10]*x18 - A[5]*x27 + A[5]*x34 - A[6]*x0 + A[6]*x1 - A[6]*x2 + A[6]*x7 - A[6]*x8 + A[6]*x9 - A[7]*x23 + A[7]*x26 + A[8]*x30 + A[8]*x32 - x36*x4 + x36*x5 - x36*x6;
    res[7]=-A[5]*x28 + A[5]*x35 - A[6]*x23 + A[6]*x26 - A[7]*x0 - A[7]*x1 + A[7]*x2 + A[7]*x7 + A[7]*x8 - A[7]*x9 + A[8]*x31 + A[8]*x33 + A[9]*x15 + A[9]*x18 - x37*x4 - x37*x5 + x37*x6;
    res[8]=-A[10]*x28 + A[10]*x35 - A[6]*x30 - A[6]*x32 - A[7]*x31 - A[7]*x33 + A[8]*x0 - A[8]*x1 - A[8]*x2 - A[8]*x7 + A[8]*x8 + A[8]*x9 - A[9]*x27 + A[9]*x34 - x38*x4 + x38*x5 + x38*x6;
    res[9]=-A[10]*x23 + A[10]*x26 - A[5]*x30 - A[5]*x32 - A[7]*x15 - A[7]*x18 - A[8]*x27 + A[8]*x34 - A[9]*x0 + A[9]*x1 - A[9]*x2 + A[9]*x7 - A[9]*x8 + A[9]*x9 + x39*x4 - x39*x5 + x39*x6;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 + A[10]*x7 + A[10]*x8 - A[10]*x9 - A[5]*x31 - A[5]*x33 - A[6]*x15 - A[6]*x18 - A[8]*x28 + A[8]*x35 - A[9]*x23 + A[9]*x26 + x4*x40 + x40*x5 - x40*x6;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x7 + A[11]*x8 + A[11]*x9 + A[12]*x15 - A[12]*x18 - A[13]*x31 + A[13]*x33 + A[14]*x30 - A[14]*x32;
    res[12]=-A[11]*x15 + A[11]*x18 + A[12]*x0 - A[12]*x1 - A[12]*x2 + A[12]*x7 - A[12]*x8 - A[12]*x9 + A[13]*x27 + A[13]*x34 + A[14]*x28 + A[14]*x35;
    res[13]=A[11]*x31 - A[11]*x33 + A[12]*x27 + A[12]*x34 - A[13]*x0 + A[13]*x1 - A[13]*x2 - A[13]*x7 + A[13]*x8 - A[13]*x9 + A[14]*x23 + A[14]*x26;
    res[14]=-A[11]*x30 + A[11]*x32 + A[12]*x28 + A[12]*x35 + A[13]*x23 + A[13]*x26 - A[14]*x0 - A[14]*x1 + A[14]*x2 - A[14]*x7 - A[14]*x8 + A[14]*x9;
    res[15]=A[15]*x0 + A[15]*x1 + A[15]*x2 - A[15]*x7 - A[15]*x8 - A[15]*x9 - x4*x41 - x41*x5 - x41*x6;
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
    T x12 = 2.0*A[1];
    T x13 = (*this)[1]*(*this)[4];
    T x14 = (*this)[2]*A[3];
    T x15 = (*this)[2]*(*this)[5];
    T x16 = (*this)[3]*x6;
    T x17 = A[1]*x10;
    T x18 = 2.0*(*this)[1];
    T x19 = 2.0*(*this)[4];
    T x20 = 2.0*x15;
    T x21 = 2.0*x13;
    T x22 = 2.0*(*this)[5];
    T x23 = 2.0*(*this)[2];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x3 - A[0]*x4 - A[0]*x5;
    res[1]=(*this)[3]*x7 + x10*x11 + x10*x14 - x12*x13 - x12*x15 + x6*x8 + x6*x9;
    res[2]=(*this)[1]*x17 + (*this)[4]*x7 - A[2]*x16 - A[2]*x20 + x14*x19 + x18*x8 + x18*x9;
    res[3]=(*this)[2]*x17 + (*this)[5]*x7 - A[3]*x16 - A[3]*x21 + x11*x22 + x23*x8 + x23*x9;
    res[4]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + x10*x8 + x10*x9 - x11*x6 - x14*x6;
    res[5]=-(*this)[1]*x7 + (*this)[4]*x17 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - x14*x18 + x19*x9;
    res[6]=-(*this)[2]*x7 + (*this)[5]*x17 - A[3]*x0 - A[3]*x1 + A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 - x11*x23 + x22*x8;
    res[7]=-A[0]*x16 - A[0]*x20 - A[0]*x21;
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
    T x3 = 2.0*A[7];
    T x4 = (*this)[0]*(*this)[3];
    T x5 = (*this)[1]*(*this)[4];
    T x6 = (*this)[2]*(*this)[5];
    T x7 = std::pow((*this)[0], 2);
    T x8 = std::pow((*this)[1], 2);
    T x9 = std::pow((*this)[2], 2);
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[2]*A[3];
    T x13 = 2.0*A[4];
    T x14 = (*this)[4]*A[5];
    T x15 = (*this)[5]*A[6];
    T x16 = 2.0*(*this)[3];
    T x17 = (*this)[1]*A[5];
    T x18 = (*this)[2]*A[6];
    T x19 = (*this)[4]*x16;
    T x20 = (*this)[5]*A[3];
    T x21 = 2.0*A[5];
    T x22 = (*this)[0]*x13;
    T x23 = 2.0*(*this)[1];
    T x24 = (*this)[3]*x13;
    T x25 = 2.0*(*this)[4];
    T x26 = A[1]*x10;
    T x27 = 2.0*A[6];
    T x28 = (*this)[2]*A[2];
    T x29 = 2.0*(*this)[5];
    T x30 = 2.0*(*this)[2];
    T x31 = A[1]*x16;
    T x32 = (*this)[5]*A[2];
    T x33 = 2.0*A[1];
    T x34 = 2.0*A[2];
    T x35 = 2.0*A[3];
    T x36 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 - A[0]*x7 - A[0]*x8 - A[0]*x9 + x3*x4 + x3*x5 + x3*x6;
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x7 + A[1]*x8 + A[1]*x9 - A[2]*x11 + A[2]*x19 - x10*x12 + x10*x14 + x10*x15 + x13*x4 - x13*x5 - x13*x6 + x16*x17 + x16*x18 + x16*x20;
    res[2]=(*this)[1]*x24 + (*this)[4]*x22 - A[1]*x11 + A[1]*x19 - A[2]*x0 + A[2]*x1 - A[2]*x2 + A[2]*x7 - A[2]*x8 + A[2]*x9 - x12*x23 + x15*x23 + x18*x25 + x20*x25 - x21*x4 + x21*x5 - x21*x6;
    res[3]=(*this)[2]*x24 - (*this)[2]*x26 + (*this)[5]*x22 + (*this)[5]*x31 - A[3]*x0 - A[3]*x1 + A[3]*x2 + A[3]*x7 + A[3]*x8 - A[3]*x9 + x14*x30 + x17*x29 - x23*x28 + x25*x32 - x27*x4 - x27*x5 + x27*x6;
    res[4]=-(*this)[1]*A[2]*x16 - (*this)[4]*A[2]*x10 + A[4]*x0 - A[4]*x1 - A[4]*x2 - A[4]*x7 + A[4]*x8 + A[4]*x9 - x10*x17 - x10*x18 - x10*x20 - x12*x16 + x14*x16 + x15*x16 - x33*x4 + x33*x5 + x33*x6;
    res[5]=-(*this)[1]*x22 - (*this)[1]*x31 + (*this)[4]*x24 - (*this)[4]*x26 - A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x7 - A[5]*x8 + A[5]*x9 - x12*x25 + x15*x25 - x18*x23 - x20*x23 + x34*x4 - x34*x5 + x34*x6;
    res[6]=-(*this)[2]*x22 - (*this)[2]*x31 + (*this)[5]*x24 - (*this)[5]*x26 - A[6]*x0 - A[6]*x1 + A[6]*x2 + A[6]*x7 + A[6]*x8 - A[6]*x9 + x14*x29 - x17*x30 - x23*x32 - x25*x28 + x35*x4 + x35*x5 - x35*x6;
    res[7]=A[7]*x0 + A[7]*x1 + A[7]*x2 - A[7]*x7 - A[7]*x8 - A[7]*x9 - x36*x4 - x36*x5 - x36*x6;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x1;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=A[0]*x7 + A[1]*x8 - A[2]*x0 - A[2]*x4 + A[2]*x5;
    res[3]=A[3]*x0 - A[3]*x4 - A[3]*x5 + A[4]*x2 + A[5]*x7;
    res[4]=A[3]*x2 - A[4]*x0 + A[4]*x4 - A[4]*x5 + A[5]*x8;
    res[5]=A[3]*x7 + A[4]*x8 - A[5]*x0 - A[5]*x4 + A[5]*x5;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
    res[1]=A[1]*x0 + A[1]*x1 + A[1]*x2;
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=A[5]*x0 - A[5]*x1 - A[5]*x2 + A[6]*x4 + A[7]*x7;
    res[6]=A[5]*x4 - A[6]*x0 + A[6]*x1 - A[6]*x2 + A[7]*x8;
    res[7]=A[5]*x7 + A[6]*x8 - A[7]*x0 - A[7]*x1 + A[7]*x2;
    res[8]=A[10]*x7 + A[8]*x0 - A[8]*x1 - A[8]*x2 + A[9]*x4;
    res[9]=A[10]*x8 + A[8]*x4 - A[9]*x0 + A[9]*x1 - A[9]*x2;
    res[10]=-A[10]*x0 - A[10]*x1 + A[10]*x2 + A[8]*x7 + A[9]*x8;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2;
    res[12]=A[12]*x0 - A[12]*x1 - A[12]*x2 + A[13]*x4 + A[14]*x7;
    res[13]=A[12]*x4 - A[13]*x0 + A[13]*x1 - A[13]*x2 + A[14]*x8;
    res[14]=A[12]*x7 + A[13]*x8 - A[14]*x0 - A[14]*x1 + A[14]*x2;
    res[15]=A[15]*x0 + A[15]*x1 + A[15]*x2;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=A[1]*x7 + A[2]*x8 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    res[4]=A[4]*x0 - A[4]*x1 - A[4]*x2 + A[5]*x4 + A[6]*x7;
    res[5]=A[4]*x4 - A[5]*x0 + A[5]*x1 - A[5]*x2 + A[6]*x8;
    res[6]=A[4]*x7 + A[5]*x8 - A[6]*x0 - A[6]*x1 + A[6]*x2;
    res[7]=A[7]*x0 + A[7]*x1 + A[7]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x2;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-A[0]*x7 - A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x5;
    res[3]=A[3]*x0 + A[3]*x1 - A[3]*x5 - A[4]*x3 - A[5]*x7;
    res[4]=-A[3]*x3 - A[4]*x0 + A[4]*x1 + A[4]*x5 - A[5]*x8;
    res[5]=-A[3]*x7 - A[4]*x8 + A[5]*x0 - A[5]*x1 + A[5]*x5;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-(*this)[2]*A[0]*x2 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x5;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-(*this)[2]*A[0]*x2 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x5;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
    res[1]=A[1]*x0 + A[1]*x1 + A[1]*x2;
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=-A[5]*x0 + A[5]*x1 + A[5]*x2 - A[6]*x4 - A[7]*x7;
    res[6]=-A[5]*x4 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[7]*x8;
    res[7]=-A[5]*x7 - A[6]*x8 + A[7]*x0 + A[7]*x1 - A[7]*x2;
    res[8]=-A[10]*x7 - A[8]*x0 + A[8]*x1 + A[8]*x2 - A[9]*x4;
    res[9]=-A[10]*x8 - A[8]*x4 + A[9]*x0 - A[9]*x1 + A[9]*x2;
    res[10]=A[10]*x0 + A[10]*x1 - A[10]*x2 - A[8]*x7 - A[9]*x8;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2;
    res[12]=A[12]*x0 - A[12]*x1 - A[12]*x2 + A[13]*x4 + A[14]*x7;
    res[13]=A[12]*x4 - A[13]*x0 + A[13]*x1 - A[13]*x2 + A[14]*x8;
    res[14]=A[12]*x7 + A[13]*x8 - A[14]*x0 - A[14]*x1 + A[14]*x2;
    res[15]=-A[15]*x0 - A[15]*x1 - A[15]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
    res[1]=-A[1]*x0 + A[1]*x1 + A[1]*x2 - A[2]*x4 - x3*x5;
    res[2]=-A[1]*x4 + A[2]*x0 - A[2]*x1 + A[2]*x2 - x5*x6;
    res[3]=-A[1]*x7 - A[2]*x8 + A[3]*x0 + A[3]*x1 - A[3]*x2;
    res[4]=-A[4]*x0 + A[4]*x1 + A[4]*x2 - A[5]*x4 - A[6]*x7;
    res[5]=-A[4]*x4 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[6]*x8;
    res[6]=-A[4]*x7 - A[5]*x8 + A[6]*x0 + A[6]*x1 - A[6]*x2;
    res[7]=-A[7]*x0 - A[7]*x1 - A[7]*x2;
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
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=-A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x5 + x4*x6;
    res[6]=A[1]*x5 - A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + x6*x7;
    res[7]=-A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + x4*x8 + x7*x9;
    res[8]=x10*x7 - x11*x9;
    res[9]=-x10*x4 + x11*x8;
    res[10]=-(*this)[0]*A[1]*x7 + (*this)[0]*A[2]*x4;
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
    T x0 = 2.0*(*this)[0];
    T x1 = (*this)[1]*x0;
    T x2 = (*this)[2]*A[2];
    T x3 = (*this)[3]*A[3];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[1], 2);
    T x6 = std::pow((*this)[2], 2);
    T x7 = std::pow((*this)[3], 2);
    T x8 = 2.0*(*this)[1];
    T x9 = A[0]*x0;
    T x10 = A[1]*x8;
    res[0]=-A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7 - A[1]*x1 - x0*x2 - x0*x3;
    res[1]=A[0]*x1 + A[1]*x4 + A[1]*x5 - A[1]*x6 - A[1]*x7 + x2*x8 + x3*x8;
    res[2]=(*this)[2]*x10 + 2.0*(*this)[2]*x3 + (*this)[2]*x9 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7;
    res[3]=(*this)[3]*x10 + 2.0*(*this)[3]*x2 + (*this)[3]*x9 + A[3]*x4 - A[3]*x5 - A[3]*x6 + A[3]*x7;
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
    T x6 = 2.0*(*this)[1];
    T x7 = std::pow((*this)[2], 2);
    T x8 = std::pow((*this)[3], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=-x0*x1 - x0*x2 - x0*x3;
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x7 - A[0]*x8 + x2*x6 + x3*x6;
    res[2]=A[1]*x4 - A[1]*x5 + A[1]*x7 - A[1]*x8 + x1*x9 + x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x7 + A[2]*x8 + x1*x10 + x10*x2;
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
    T x1 = 2.0*(*this)[2];
    T x2 = (*this)[0]*A[5];
    T x3 = 2.0*(*this)[3];
    T x4 = (*this)[0]*x3;
    T x5 = (*this)[1]*x1;
    T x6 = (*this)[1]*x3;
    T x7 = std::pow((*this)[0], 2);
    T x8 = std::pow((*this)[2], 2);
    T x9 = std::pow((*this)[3], 2);
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*x1;
    T x12 = (*this)[0]*x10;
    T x13 = (*this)[0]*x1;
    res[0]=A[0]*x0 - A[0]*x7 - A[0]*x8 - A[0]*x9 + A[1]*x5 + A[2]*x6 - A[4]*x4 + x1*x2;
    res[1]=A[0]*x5 - A[1]*x0 - A[1]*x7 + A[1]*x8 - A[1]*x9 + A[2]*x11 + A[3]*x4 - x10*x2;
    res[2]=A[0]*x6 + A[1]*x11 - A[2]*x0 - A[2]*x7 - A[2]*x8 + A[2]*x9 - A[3]*x13 + A[4]*x12;
    res[3]=-A[1]*x4 + A[2]*x13 - A[3]*x0 + A[3]*x7 + A[3]*x8 + A[3]*x9 - A[4]*x5 - A[5]*x6;
    res[4]=A[0]*x4 - A[2]*x12 - A[3]*x5 + A[4]*x0 + A[4]*x7 - A[4]*x8 + A[4]*x9 - A[5]*x11;
    res[5]=-A[0]*x13 + A[1]*x12 - A[3]*x6 - A[4]*x11 + A[5]*x0 + A[5]*x7 + A[5]*x8 - A[5]*x9;
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
    T x8 = 2.0*(*this)[2];
    T x9 = 2.0*(*this)[1];
    T x10 = (*this)[3]*A[2];
    T x11 = std::pow((*this)[1], 2);
    res[0]=(*this)[2]*x1 - A[1]*x2;
    res[1]=-(*this)[1]*x1 + A[0]*x2;
    res[2]=x0*x3 - x0*x4;
    res[3]=-A[0]*x11 + A[0]*x5 + A[0]*x6 + A[0]*x7 - x10*x9 - x3*x8;
    res[4]=A[1]*x11 + A[1]*x5 - A[1]*x6 + A[1]*x7 - x10*x8 - x4*x9;
    res[5]=-(*this)[3]*A[0]*x9 - (*this)[3]*A[1]*x8 + A[2]*x11 + A[2]*x5 + A[2]*x6 - A[2]*x7;
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
    T x1 = 2.0*(*this)[1];
    T x2 = (*this)[2]*x1;
    T x3 = (*this)[3]*A[2];
    T x4 = std::pow((*this)[0], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = std::pow((*this)[3], 2);
    T x7 = 2.0*(*this)[2];
    T x8 = (*this)[3]*A[0];
    T x9 = (*this)[3]*A[1];
    T x10 = (*this)[0]*A[2];
    T x11 = 2.0*(*this)[0];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 - A[0]*x6 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 - A[1]*x4 + A[1]*x5 - A[1]*x6 + x3*x7;
    res[2]=-A[2]*x0 - A[2]*x4 - A[2]*x5 + A[2]*x6 + x1*x8 + x7*x9;
    res[3]=x10*x7 - x11*x9;
    res[4]=-x1*x10 + x11*x8;
    res[5]=-(*this)[0]*A[0]*x7 + (*this)[0]*A[1]*x1;
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
    T x6 = 2.0*(*this)[1];
    T x7 = std::pow((*this)[0], 2);
    T x8 = std::pow((*this)[1], 2);
    T x9 = 2.0*(*this)[2];
    T x10 = 2.0*(*this)[3];
    res[0]=-x0*x1 - x0*x2 - x0*x3;
    res[1]=A[0]*x4 + A[0]*x5 - A[0]*x7 - A[0]*x8 - x2*x6 - x3*x6;
    res[2]=-A[1]*x4 + A[1]*x5 - A[1]*x7 + A[1]*x8 - x1*x9 - x3*x9;
    res[3]=A[2]*x4 - A[2]*x5 - A[2]*x7 + A[2]*x8 - x1*x10 - x10*x2;
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
    res[0]=-A[0]*(std::pow((*this)[0], 2) - std::pow((*this)[1], 2) - std::pow((*this)[2], 2) - std::pow((*this)[3], 2));
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
    T x13 = (*this)[3]*x4;
    T x14 = (*this)[2]*x8;
    T x15 = (*this)[3]*A[7];
    T x16 = (*this)[3]*x8;
    T x17 = (*this)[3]*x11;
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=-A[1]*x0 - A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x5 - x4*x6 - x4*x7;
    res[2]=A[1]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + x6*x8 + x7*x8;
    res[3]=(*this)[2]*x10 + (*this)[2]*x9 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 + x11*x7;
    res[4]=(*this)[3]*x10 + 2.0*(*this)[3]*x6 + (*this)[3]*x9 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3;
    res[5]=A[10]*x12 - A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[6]*x14 - A[9]*x13 + x15*x8;
    res[6]=-A[10]*x5 + A[5]*x14 - A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3 + A[8]*x13 + x11*x15;
    res[7]=A[5]*x16 + A[6]*x17 - A[7]*x0 - A[7]*x1 - A[7]*x2 + A[7]*x3 - A[8]*x12 + A[9]*x5;
    res[8]=-A[10]*x16 - A[6]*x13 + A[7]*x12 + A[8]*x0 - A[8]*x1 + A[8]*x2 + A[8]*x3 - A[9]*x14;
    res[9]=-A[10]*x17 + A[5]*x13 - A[7]*x5 - A[8]*x14 + A[9]*x0 + A[9]*x1 - A[9]*x2 + A[9]*x3;
    res[10]=A[10]*x0 + A[10]*x1 + A[10]*x2 - A[10]*x3 - A[5]*x12 + A[6]*x5 - A[8]*x16 - A[9]*x17;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 - A[12]*x5 - A[13]*x12 - A[14]*x13;
    res[12]=A[11]*x5 - A[12]*x0 - A[12]*x1 + A[12]*x2 + A[12]*x3 - A[13]*x14 - A[14]*x16;
    res[13]=A[11]*x12 - A[12]*x14 - A[13]*x0 + A[13]*x1 - A[13]*x2 + A[13]*x3 - A[14]*x17;
    res[14]=A[11]*x13 - A[12]*x16 - A[13]*x17 - A[14]*x0 + A[14]*x1 + A[14]*x2 - A[14]*x3;
    res[15]=-A[15]*x0 + A[15]*x1 + A[15]*x2 + A[15]*x3;
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
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=0;
    res[2]=0;
    res[3]=0;
    res[4]=0;
    res[5]=(*this)[2]*x5 - A[2]*x6;
    res[6]=-(*this)[1]*x5 + A[1]*x6;
    res[7]=x4*x7 - x4*x8;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*x6;
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[3]*x4;
    T x12 = (*this)[0]*x10;
    T x13 = (*this)[0]*x4;
    res[0]=A[0]*x0 - A[0]*x1 - A[0]*x2 - A[0]*x3;
    res[1]=-A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[2]*x8 + A[3]*x9 - A[5]*x7 + x4*x5;
    res[2]=A[1]*x8 - A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[3]*x11 + A[4]*x7 - x10*x5;
    res[3]=A[1]*x9 + A[2]*x11 - A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 - A[4]*x13 + A[5]*x12;
    res[4]=-A[2]*x7 + A[3]*x13 + A[4]*x0 - A[4]*x1 + A[4]*x2 + A[4]*x3 - A[5]*x8 - A[6]*x9;
    res[5]=A[1]*x7 - A[3]*x12 - A[4]*x8 + A[5]*x0 + A[5]*x1 - A[5]*x2 + A[5]*x3 - A[6]*x11;
    res[6]=-A[1]*x13 + A[2]*x12 - A[4]*x9 - A[5]*x11 + A[6]*x0 + A[6]*x1 + A[6]*x2 - A[6]*x3;
    res[7]=-A[7]*x0 + A[7]*x1 + A[7]*x2 + A[7]*x3;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    T x7 = (*this)[2]*x1;
    T x8 = (*this)[2]*x6;
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=A[0]*x7 + A[1]*x8 - A[2]*x0 - A[2]*x4 + A[2]*x5;
    res[3]=-A[3]*x0 + A[3]*x4 + A[3]*x5 - A[4]*x2 - A[5]*x7;
    res[4]=-A[3]*x2 + A[4]*x0 - A[4]*x4 + A[4]*x5 - A[5]*x8;
    res[5]=-A[3]*x7 - A[4]*x8 + A[5]*x0 + A[5]*x4 - A[5]*x5;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-(*this)[2]*A[0]*x2 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x5;
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
    T x1 = 2.0*(*this)[0];
    T x2 = (*this)[1]*x1;
    T x3 = (*this)[2]*A[2];
    T x4 = std::pow((*this)[1], 2);
    T x5 = std::pow((*this)[2], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 - A[0]*x4 - A[0]*x5 + A[1]*x2 + x1*x3;
    res[1]=A[0]*x2 - A[1]*x0 + A[1]*x4 - A[1]*x5 + x3*x6;
    res[2]=(*this)[2]*A[0]*x1 + (*this)[2]*A[1]*x6 - A[2]*x0 - A[2]*x4 + A[2]*x5;
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
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2;
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
    T x2 = 2.0*(*this)[0];
    T x3 = (*this)[1]*x2;
    T x4 = (*this)[2]*A[2];
    T x5 = std::pow((*this)[0], 2);
    T x6 = 2.0*(*this)[1];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x5 - A[1]*x3 - x2*x4;
    res[1]=-A[0]*x3 - A[1]*x0 + A[1]*x1 + A[1]*x5 - x4*x6;
    res[2]=-(*this)[2]*A[0]*x2 - (*this)[2]*A[1]*x6 + A[2]*x0 - A[2]*x1 + A[2]*x5;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
    res[1]=-A[1]*x0 - A[1]*x1 - A[1]*x2;
    res[2]=A[2]*x0 - A[2]*x1 - A[2]*x2 + A[3]*x4 + x3*x5;
    res[3]=A[2]*x4 - A[3]*x0 + A[3]*x1 - A[3]*x2 + x5*x6;
    res[4]=A[2]*x7 + A[3]*x8 - A[4]*x0 - A[4]*x1 + A[4]*x2;
    res[5]=A[5]*x0 - A[5]*x1 - A[5]*x2 + A[6]*x4 + A[7]*x7;
    res[6]=A[5]*x4 - A[6]*x0 + A[6]*x1 - A[6]*x2 + A[7]*x8;
    res[7]=A[5]*x7 + A[6]*x8 - A[7]*x0 - A[7]*x1 + A[7]*x2;
    res[8]=-A[10]*x7 - A[8]*x0 + A[8]*x1 + A[8]*x2 - A[9]*x4;
    res[9]=-A[10]*x8 - A[8]*x4 + A[9]*x0 - A[9]*x1 + A[9]*x2;
    res[10]=A[10]*x0 + A[10]*x1 - A[10]*x2 - A[8]*x7 - A[9]*x8;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2;
    res[12]=-A[12]*x0 + A[12]*x1 + A[12]*x2 - A[13]*x4 - A[14]*x7;
    res[13]=-A[12]*x4 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[14]*x8;
    res[14]=-A[12]*x7 - A[13]*x8 + A[14]*x0 + A[14]*x1 - A[14]*x2;
    res[15]=A[15]*x0 + A[15]*x1 + A[15]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
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
    res[0]=-A[0]*x0 - A[0]*x1 - A[0]*x2;
    res[1]=A[1]*x0 - A[1]*x1 - A[1]*x2 + A[2]*x4 + x3*x5;
    res[2]=A[1]*x4 - A[2]*x0 + A[2]*x1 - A[2]*x2 + x5*x6;
    res[3]=A[1]*x7 + A[2]*x8 - A[3]*x0 - A[3]*x1 + A[3]*x2;
    res[4]=-A[4]*x0 + A[4]*x1 + A[4]*x2 - A[5]*x4 - A[6]*x7;
    res[5]=-A[4]*x4 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[6]*x8;
    res[6]=-A[4]*x7 - A[5]*x8 + A[6]*x0 + A[6]*x1 - A[6]*x2;
    res[7]=A[7]*x0 + A[7]*x1 + A[7]*x2;
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
    T x20 = 2.0*A[3];
    T x21 = (*this)[11]*x20;
    T x22 = 2.0*A[2];
    T x23 = (*this)[10]*x22;
    T x24 = 2.0*A[1];
    T x25 = (*this)[3]*x24;
    T x26 = (*this)[11]*x24;
    T x27 = (*this)[11]*x22;
    T x28 = (*this)[12]*x24;
    T x29 = (*this)[12]*x20;
    T x30 = (*this)[12]*x22;
    T x31 = (*this)[15]*x22;
    T x32 = (*this)[5]*x20;
    T x33 = (*this)[7]*x24;
    T x34 = (*this)[14]*x20;
    T x35 = (*this)[5]*x22;
    T x36 = (*this)[14]*x24;
    T x37 = (*this)[5]*x24;
    T x38 = (*this)[6]*x22;
    T x39 = (*this)[7]*x20;
    T x40 = (*this)[9]*x20;
    T x41 = (*this)[3]*x20;
    T x42 = (*this)[8]*x22;
    T x43 = (*this)[9]*x24;
    T x44 = (*this)[10]*x20;
    T x45 = (*this)[10]*x24;
    T x46 = (*this)[13]*x24;
    T x47 = (*this)[9]*x22;
    T x48 = (*this)[15]*x20;
    T x49 = (*this)[15]*x24;
    T x50 = (*this)[1]*x20;
    T x51 = (*this)[2]*x24;
    T x52 = (*this)[4]*x20;
    T x53 = (*this)[7]*x22;
    T x54 = (*this)[14]*x22;
    T x55 = (*this)[3]*x22;
    T x56 = (*this)[13]*x22;
    T x57 = (*this)[1]*(*this)[4];
    T x58 = A[0]*x16;
    T x59 = 2.0*A[0];
    T x60 = (*this)[10]*x59;
    T x61 = (*this)[5]*x59;
    T x62 = (*this)[6]*x59;
    T x63 = (*this)[14]*x59;
    T x64 = (*this)[1]*x59;
    T x65 = (*this)[2]*x59;
    T x66 = (*this)[3]*x59;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - A[0]*x9;
    res[1]=-(*this)[10]*x21 + (*this)[10]*x25 - (*this)[13]*x31 - (*this)[13]*x32 + (*this)[13]*x33 + (*this)[14]*x35 - (*this)[15]*x28 - (*this)[15]*x34 - (*this)[1]*x37 - (*this)[1]*x38 - (*this)[1]*x39 + (*this)[2]*x17 - (*this)[2]*x23 + (*this)[2]*x40 + (*this)[3]*x18 + (*this)[4]*x19 + (*this)[4]*x42 - (*this)[4]*x43 + (*this)[6]*x29 - (*this)[6]*x36 - (*this)[7]*x30 - (*this)[8]*x26 - (*this)[8]*x41 - (*this)[9]*x27;
    res[2]=(*this)[10]*x29 - (*this)[10]*x36 + (*this)[13]*x19 + (*this)[13]*x42 - (*this)[13]*x43 - (*this)[14]*x18 + (*this)[15]*x26 + (*this)[15]*x41 + (*this)[1]*x17 - (*this)[1]*x23 + (*this)[1]*x40 - (*this)[2]*x37 - (*this)[2]*x38 - (*this)[2]*x39 - (*this)[3]*x35 - (*this)[4]*x31 - (*this)[4]*x32 + (*this)[4]*x33 - (*this)[6]*x21 + (*this)[6]*x25 + (*this)[7]*x27 + (*this)[8]*x28 + (*this)[8]*x34 + (*this)[9]*x30;
    res[3]=-(*this)[11]*x33 - (*this)[12]*x19 - (*this)[12]*x42 + (*this)[13]*x44 + (*this)[13]*x47 + (*this)[14]*x17 - (*this)[14]*x23 + (*this)[14]*x40 + (*this)[15]*x27 + (*this)[1]*x18 + (*this)[1]*x45 + (*this)[2]*x35 - (*this)[2]*x48 - (*this)[3]*x38 - (*this)[3]*x39 + (*this)[4]*x49 + (*this)[4]*x53 + (*this)[5]*x21 - (*this)[5]*x25 - (*this)[6]*x51 - (*this)[6]*x52 + (*this)[8]*x46 - (*this)[8]*x50 + (*this)[9]*x28;
    res[4]=(*this)[10]*x28 + (*this)[10]*x34 - (*this)[11]*x35 + (*this)[12]*x18 - (*this)[13]*x17 + (*this)[13]*x23 - (*this)[13]*x40 + (*this)[14]*x47 + (*this)[15]*x21 - (*this)[15]*x25 + (*this)[1]*x19 + (*this)[1]*x42 - (*this)[1]*x43 + (*this)[2]*x31 + (*this)[2]*x32 - (*this)[2]*x33 - (*this)[3]*x53 - (*this)[4]*x37 - (*this)[4]*x38 - (*this)[4]*x39 + (*this)[6]*x26 + (*this)[6]*x41 - (*this)[8]*x29 + (*this)[8]*x36;
    res[5]=-(*this)[10]*x18 + (*this)[13]*x30 + (*this)[13]*x50 + (*this)[14]*x29 - (*this)[1]*x54 - (*this)[2]*x52 - (*this)[2]*x55 - (*this)[3]*x21 + (*this)[4]*x27 - (*this)[6]*x35 + (*this)[6]*x48 - (*this)[7]*x31 - (*this)[7]*x32 + (*this)[8]*x44 + (*this)[9]*x19 + (*this)[9]*x42 + A[1]*x0 - A[1]*x1 - A[1]*x10 - A[1]*x11 + A[1]*x12 - A[1]*x13 + A[1]*x14 + A[1]*x15 - A[1]*x2 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 + A[1]*x8 - A[1]*x9;
    res[6]=(*this)[10]*x17 + (*this)[10]*x40 + (*this)[13]*x28 + (*this)[13]*x34 - (*this)[15]*x32 + (*this)[15]*x33 - (*this)[1]*x29 + (*this)[1]*x36 + (*this)[2]*x21 - (*this)[2]*x25 - (*this)[4]*x26 - (*this)[4]*x41 - (*this)[6]*x37 - (*this)[6]*x39 - (*this)[8]*x19 + (*this)[8]*x43 + A[2]*x0 - A[2]*x1 - A[2]*x10 - A[2]*x11 + A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 - A[2]*x2 + A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 - A[2]*x8 + A[2]*x9;
    res[7]=(*this)[11]*x25 + (*this)[13]*x54 + (*this)[14]*x28 + (*this)[15]*x35 + (*this)[1]*x30 - (*this)[1]*x46 - (*this)[2]*x27 - (*this)[4]*x51 - (*this)[4]*x55 - (*this)[5]*x33 - (*this)[6]*x49 - (*this)[7]*x38 + (*this)[8]*x18 + (*this)[8]*x45 - (*this)[9]*x17 + (*this)[9]*x23 + A[3]*x0 + A[3]*x1 + A[3]*x10 - A[3]*x11 + A[3]*x12 + A[3]*x13 + A[3]*x14 - A[3]*x15 - A[3]*x2 + A[3]*x3 + A[3]*x4 - A[3]*x5 - A[3]*x6 - A[3]*x7 - A[3]*x8 - A[3]*x9;
    res[8]=-(*this)[10]*x32 + (*this)[10]*x33 + (*this)[13]*x21 - (*this)[13]*x25 - (*this)[14]*x27 + (*this)[15]*x17 - (*this)[15]*x23 + (*this)[15]*x40 + (*this)[1]*x26 + (*this)[1]*x41 + (*this)[2]*x28 + (*this)[2]*x34 + (*this)[2]*x56 + (*this)[3]*x30 + (*this)[4]*x29 - (*this)[4]*x36 - (*this)[6]*x19 - (*this)[6]*x42 + (*this)[6]*x43 + (*this)[7]*x18 - (*this)[8]*x37 - (*this)[8]*x39 - (*this)[9]*x35 - x22*x57;
    res[9]=-(*this)[11]*x29 + (*this)[12]*x25 + (*this)[13]*x52 + (*this)[13]*x55 + (*this)[14]*x26 + (*this)[15]*x18 + (*this)[15]*x45 + (*this)[1]*x27 - (*this)[2]*x30 + (*this)[2]*x46 - (*this)[2]*x50 + (*this)[3]*x34 - (*this)[4]*x54 + (*this)[5]*x19 - (*this)[6]*(*this)[8]*x24 - (*this)[6]*x44 - (*this)[7]*x17 + (*this)[7]*x23 - (*this)[7]*x40 + (*this)[8]*x35 - (*this)[8]*x48 - (*this)[9]*x37 - (*this)[9]*x38 + x24*x57;
    res[10]=-(*this)[10]*x37 - (*this)[10]*x39 + (*this)[12]*x27 - (*this)[13]*x26 - (*this)[13]*x41 + (*this)[15]*x19 + (*this)[15]*x42 - (*this)[15]*x43 + (*this)[1]*(*this)[2]*x22 + (*this)[1]*x21 - (*this)[1]*x25 - (*this)[2]*x29 + (*this)[2]*x36 + (*this)[3]*x54 + (*this)[4]*x28 + (*this)[4]*x34 + (*this)[4]*x56 - (*this)[5]*x18 + (*this)[6]*x17 - (*this)[6]*x23 + (*this)[6]*x40 - (*this)[7]*x47 + (*this)[8]*x32 - (*this)[8]*x33;
    res[11]=(*this)[11]*x58 + (*this)[12]*x61 + (*this)[13]*x62 - (*this)[15]*x64 + (*this)[4]*x60 + (*this)[7]*x63 + (*this)[8]*x65 + (*this)[9]*x66;
    res[12]=(*this)[11]*x61 + (*this)[12]*x58 - (*this)[13]*x60 + (*this)[15]*x65 - (*this)[4]*x62 + (*this)[7]*x66 - (*this)[8]*x64 + (*this)[9]*x63;
    res[13]=(*this)[11]*x62 + (*this)[12]*x60 + (*this)[13]*x58 + (*this)[15]*x66 + (*this)[4]*x61 - (*this)[7]*x65 - (*this)[8]*x63 - (*this)[9]*x64;
    res[14]=(*this)[11]*(*this)[7]*x59 - (*this)[12]*(*this)[9]*x59 + (*this)[13]*(*this)[8]*x59 + (*this)[14]*x58 + (*this)[15]*(*this)[4]*x59 - (*this)[1]*x60 + (*this)[2]*x62 - (*this)[3]*x61;
    res[15]=-(*this)[11]*x64 - (*this)[12]*x65 - (*this)[13]*x66 + (*this)[15]*x58 - (*this)[4]*x63 - (*this)[7]*x60 - (*this)[8]*x61 - (*this)[9]*x62;
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
    T x8 = 2.0*A[1];
    T x9 = (*this)[0]*(*this)[5];
    T x10 = 2.0*A[2];
    T x11 = (*this)[0]*x10;
    T x12 = 2.0*A[3];
    T x13 = (*this)[0]*x12;
    T x14 = (*this)[15]*x12;
    T x15 = (*this)[5]*x10;
    T x16 = (*this)[10]*(*this)[6];
    T x17 = (*this)[11]*(*this)[12];
    T x18 = (*this)[11]*x10;
    T x19 = (*this)[11]*x12;
    T x20 = (*this)[12]*x12;
    T x21 = (*this)[12]*x10;
    T x22 = 2.0*(*this)[2];
    T x23 = (*this)[13]*A[3];
    T x24 = (*this)[13]*(*this)[4];
    T x25 = (*this)[14]*x10;
    T x26 = (*this)[14]*(*this)[3];
    T x27 = (*this)[15]*(*this)[8];
    T x28 = (*this)[15]*x10;
    T x29 = (*this)[1]*x22;
    T x30 = (*this)[3]*x10;
    T x31 = (*this)[4]*x12;
    T x32 = (*this)[5]*x12;
    T x33 = (*this)[8]*x12;
    T x34 = (*this)[8]*x10;
    T x35 = (*this)[7]*(*this)[9];
    T x36 = std::pow((*this)[11], 2);
    T x37 = std::pow((*this)[12], 2);
    T x38 = std::pow((*this)[13], 2);
    T x39 = std::pow((*this)[14], 2);
    T x40 = std::pow((*this)[1], 2);
    T x41 = std::pow((*this)[2], 2);
    T x42 = std::pow((*this)[3], 2);
    T x43 = std::pow((*this)[4], 2);
    T x44 = 2.0*A[0];
    T x45 = 2.0*x23;
    T x46 = A[3]*x22;
    T x47 = (*this)[10]*x8;
    T x48 = (*this)[6]*x44;
    T x49 = (*this)[5]*x44;
    T x50 = (*this)[10]*(*this)[9];
    T x51 = (*this)[11]*x44;
    T x52 = (*this)[11]*x8;
    T x53 = (*this)[12]*x8;
    T x54 = (*this)[12]*x44;
    T x55 = (*this)[1]*x8;
    T x56 = A[0]*x22;
    T x57 = (*this)[7]*x8;
    T x58 = (*this)[15]*x44;
    T x59 = (*this)[1]*x44;
    T x60 = A[1]*x22;
    T x61 = (*this)[6]*x8;
    T x62 = (*this)[6]*(*this)[7];
    T x63 = (*this)[7]*x44;
    T x64 = (*this)[9]*x8;
    T x65 = (*this)[10]*x10;
    T x66 = (*this)[3]*x44;
    T x67 = (*this)[9]*x12;
    T x68 = (*this)[6]*x10;
    T x69 = (*this)[4]*x44;
    T x70 = (*this)[14]*x8;
    T x71 = (*this)[13]*x8;
    T x72 = (*this)[9]*x10;
    T x73 = (*this)[4]*x8;
    T x74 = (*this)[3]*x8;
    T x75 = (*this)[3]*x12;
    T x76 = (*this)[7]*x10;
    T x77 = (*this)[14]*x12;
    T x78 = (*this)[13]*x44;
    T x79 = (*this)[14]*x44;
    res[0]=0;
    res[1]=-(*this)[10]*x14 + (*this)[10]*x15 - (*this)[13]*x18 - (*this)[14]*x19 + (*this)[1]*x30 + (*this)[1]*x31 - (*this)[2]*x25 - (*this)[3]*x20 + (*this)[4]*x21 - (*this)[6]*x11 + (*this)[6]*x33 - (*this)[7]*x13 - (*this)[7]*x34 - (*this)[9]*x28 - (*this)[9]*x32 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x36 - A[0]*x37 - A[0]*x38 - A[0]*x39 + A[0]*x4 - A[0]*x40 - A[0]*x41 - A[0]*x42 - A[0]*x43 + A[0]*x5 + A[0]*x6 + A[0]*x7 + A[1]*x29 - x16*x8 - x17*x8 + x22*x23 - x24*x8 + x26*x8 - x27*x8 + x35*x8 - x8*x9;
    res[2]=-(*this)[10]*x11 + (*this)[10]*x33 + (*this)[13]*x21 + (*this)[14]*x20 - (*this)[1]*x25 + (*this)[1]*x45 + (*this)[2]*x30 + (*this)[3]*x19 - (*this)[4]*x18 + (*this)[4]*x46 - (*this)[6]*x14 + (*this)[6]*x15 + (*this)[7]*x28 + (*this)[7]*x32 + (*this)[9]*x13 + (*this)[9]*x34 - A[0]*x29 + A[1]*x0 - A[1]*x1 + A[1]*x2 + A[1]*x3 + A[1]*x36 + A[1]*x37 - A[1]*x38 - A[1]*x39 - A[1]*x4 + A[1]*x40 + A[1]*x41 - A[1]*x42 - A[1]*x43 - A[1]*x5 + A[1]*x6 - A[1]*x7 + x16*x44 + x17*x44 - x24*x44 + x26*x44 - x27*x44 - x35*x44 - x44*x9;
    res[3]=-(*this)[0]*x33 + (*this)[0]*x47 - (*this)[0]*x48 - (*this)[10]*x49 - (*this)[11]*x46 + (*this)[13]*x51 + (*this)[13]*x53 + (*this)[14]*x45 + (*this)[14]*x55 - (*this)[14]*x56 - (*this)[15]*x57 - (*this)[1]*x20 + (*this)[3]*x31 - (*this)[3]*x59 + (*this)[3]*x60 + (*this)[4]*x52 + (*this)[4]*x54 + (*this)[5]*x14 + (*this)[5]*x61 + (*this)[8]*x63 + (*this)[8]*x64 - (*this)[9]*x58 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[2]*x36 - A[2]*x37 + A[2]*x38 - A[2]*x39 + A[2]*x4 + A[2]*x40 - A[2]*x41 + A[2]*x42 - A[2]*x43 - A[2]*x5 - A[2]*x6 + A[2]*x7 + x12*x50 + x12*x62;
    res[4]=-(*this)[0]*x63 - (*this)[0]*x64 - (*this)[10]*x58 + (*this)[13]*x25 - (*this)[13]*x55 + (*this)[13]*x56 + (*this)[14]*x51 + (*this)[14]*x53 - (*this)[15]*x15 + (*this)[15]*x61 + (*this)[1]*x21 + (*this)[2]*x18 - (*this)[3]*x52 - (*this)[3]*x54 + (*this)[4]*x30 - (*this)[4]*x59 + (*this)[4]*x60 + (*this)[5]*x57 + (*this)[8]*x11 + (*this)[8]*x47 - (*this)[8]*x48 + (*this)[9]*x49 + A[3]*x0 + A[3]*x1 + A[3]*x2 - A[3]*x3 + A[3]*x36 - A[3]*x37 - A[3]*x38 + A[3]*x39 - A[3]*x4 + A[3]*x40 - A[3]*x41 - A[3]*x42 + A[3]*x43 + A[3]*x5 - A[3]*x6 - A[3]*x7 + x10*x50 + x10*x62;
    res[5]=(*this)[0]*x45 + (*this)[0]*x55 - (*this)[0]*x56 + (*this)[10]*x20 + (*this)[10]*x66 + (*this)[13]*x34 - (*this)[13]*x63 - (*this)[13]*x64 - (*this)[14]*x11 + (*this)[14]*x33 - (*this)[14]*x47 + (*this)[14]*x48 - (*this)[15]*x52 - (*this)[15]*x54 - (*this)[1]*x49 - (*this)[1]*x65 + (*this)[1]*x67 + (*this)[2]*x68 - (*this)[3]*x14 + (*this)[3]*x15 - (*this)[3]*x61 + (*this)[4]*x28 - (*this)[4]*x57 + (*this)[5]*x31 + (*this)[5]*x60 + (*this)[6]*x19 - (*this)[7]*x18 + (*this)[7]*x46 + (*this)[8]*x51 + (*this)[8]*x53 + (*this)[9]*x21 - (*this)[9]*x69;
    res[6]=-(*this)[0]*x66 + (*this)[0]*x70 - (*this)[10]*x25 + (*this)[10]*x45 - (*this)[10]*x56 - (*this)[12]*x13 - (*this)[13]*x58 + (*this)[13]*x72 - (*this)[14]*x49 + (*this)[14]*x67 - (*this)[15]*x18 + (*this)[15]*x46 - (*this)[15]*x73 + (*this)[1]*x11 - (*this)[1]*x33 + (*this)[1]*x47 - (*this)[1]*x48 - (*this)[2]*x15 - (*this)[4]*x76 - (*this)[5]*x19 + (*this)[5]*x74 + (*this)[6]*x30 + (*this)[6]*x31 + (*this)[6]*x60 + (*this)[7]*x52 + (*this)[7]*x54 + (*this)[7]*x75 - (*this)[8]*x21 + (*this)[8]*x69 + (*this)[8]*x71 + (*this)[9]*x51 + (*this)[9]*x53;
    res[7]=(*this)[0]*x21 - (*this)[0]*x69 - (*this)[0]*x71 + (*this)[10]*x51 + (*this)[10]*x77 - (*this)[11]*x14 + (*this)[11]*x15 - (*this)[12]*x33 + (*this)[12]*x47 + (*this)[13]*x49 + (*this)[13]*x65 - (*this)[14]*x58 + (*this)[15]*x74 + (*this)[1]*x13 + (*this)[1]*x34 - (*this)[1]*x63 - (*this)[2]*x28 + (*this)[4]*x68 - (*this)[5]*x46 + (*this)[5]*x73 - (*this)[6]*x52 - (*this)[6]*x54 - (*this)[6]*x75 + (*this)[7]*x30 + (*this)[7]*x31 + (*this)[7]*x60 - (*this)[8]*x66 + (*this)[8]*x70 + (*this)[9]*x25 - (*this)[9]*x45 - (*this)[9]*x55 + (*this)[9]*x56;
    res[8]=(*this)[0]*x52 + (*this)[0]*x54 - (*this)[10]*x18 + (*this)[10]*x46 - (*this)[10]*x78 - (*this)[13]*x15 + (*this)[13]*x61 - (*this)[14]*x32 + (*this)[14]*x57 - (*this)[15]*x25 + (*this)[15]*x45 + (*this)[15]*x55 - (*this)[15]*x56 - (*this)[1]*(*this)[6]*x12 + (*this)[1]*x76 + (*this)[2]*x72 + (*this)[3]*x13 - (*this)[3]*x63 - (*this)[3]*x64 - (*this)[4]*x11 - (*this)[4]*x47 + (*this)[4]*x48 - (*this)[5]*x51 - (*this)[5]*x53 - (*this)[6]*x21 - (*this)[7]*x20 + (*this)[8]*x30 + (*this)[8]*x31 - (*this)[8]*x59 + (*this)[8]*x60 + (*this)[9]*x19 + (*this)[9]*x79;
    res[9]=-(*this)[0]*x46 + (*this)[0]*x73 + (*this)[0]*x78 + (*this)[10]*x54 + (*this)[10]*x75 + (*this)[11]*x11 - (*this)[11]*x33 + (*this)[11]*x47 - (*this)[12]*x14 + (*this)[12]*x15 - (*this)[13]*x68 + (*this)[15]*x70 + (*this)[1]*x28 + (*this)[1]*x32 - (*this)[2]*x34 - (*this)[3]*x58 - (*this)[4]*x49 - (*this)[4]*x65 - (*this)[5]*x71 - (*this)[6]*x51 - (*this)[6]*x53 - (*this)[6]*x77 + (*this)[7]*x25 - (*this)[7]*x45 - (*this)[7]*x55 + (*this)[7]*x56 + (*this)[8]*x74 - (*this)[8]*x79 + (*this)[9]*x30 + (*this)[9]*x31 - (*this)[9]*x59 + (*this)[9]*x60;
    res[10]=-(*this)[0]*x74 + (*this)[0]*x79 + (*this)[10]*x30 + (*this)[10]*x31 - (*this)[10]*x59 + (*this)[10]*x60 + (*this)[11]*x13 - (*this)[13]*x76 + (*this)[15]*x21 - (*this)[15]*x71 + (*this)[1]*x14 - (*this)[1]*x15 + (*this)[2]*x11 + (*this)[3]*x49 - (*this)[3]*x67 - (*this)[4]*x58 + (*this)[4]*x72 + (*this)[5]*x20 - (*this)[5]*x70 - (*this)[6]*x25 + (*this)[6]*x45 + (*this)[6]*x55 - (*this)[6]*x56 - (*this)[7]*x51 - (*this)[7]*x53 - (*this)[7]*x77 + (*this)[8]*x18 - (*this)[8]*x46 + (*this)[8]*x73 + (*this)[8]*x78 - (*this)[9]*x52 - (*this)[9]*x54;
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
    T x0 = 2.0*A[0];
    T x1 = (*this)[0]*x0;
    T x2 = 2.0*A[1];
    T x3 = (*this)[0]*x2;
    T x4 = 2.0*A[2];
    T x5 = (*this)[0]*x4;
    T x6 = (*this)[15]*x4;
    T x7 = (*this)[5]*x2;
    T x8 = (*this)[6]*x0;
    T x9 = (*this)[11]*x0;
    T x10 = (*this)[11]*x2;
    T x11 = (*this)[11]*x4;
    T x12 = (*this)[12]*x4;
    T x13 = (*this)[12]*x2;
    T x14 = 2.0*(*this)[2];
    T x15 = (*this)[13]*A[2];
    T x16 = (*this)[13]*x0;
    T x17 = (*this)[14]*x2;
    T x18 = (*this)[14]*x0;
    T x19 = (*this)[15]*x0;
    T x20 = (*this)[15]*x2;
    T x21 = A[0]*x14;
    T x22 = (*this)[3]*x2;
    T x23 = (*this)[4]*x4;
    T x24 = (*this)[5]*x4;
    T x25 = (*this)[8]*x4;
    T x26 = (*this)[8]*x2;
    T x27 = (*this)[9]*x0;
    T x28 = std::pow((*this)[0], 2);
    T x29 = std::pow((*this)[11], 2);
    T x30 = std::pow((*this)[12], 2);
    T x31 = std::pow((*this)[15], 2);
    T x32 = std::pow((*this)[1], 2);
    T x33 = std::pow((*this)[2], 2);
    T x34 = std::pow((*this)[5], 2);
    T x35 = std::pow((*this)[8], 2);
    T x36 = 2.0*x15;
    T x37 = A[2]*x14;
    T x38 = std::pow((*this)[10], 2);
    T x39 = std::pow((*this)[13], 2);
    T x40 = std::pow((*this)[14], 2);
    T x41 = std::pow((*this)[3], 2);
    T x42 = std::pow((*this)[4], 2);
    T x43 = std::pow((*this)[6], 2);
    T x44 = std::pow((*this)[7], 2);
    T x45 = std::pow((*this)[9], 2);
    T x46 = (*this)[10]*(*this)[9];
    T x47 = (*this)[6]*(*this)[7];
    T x48 = (*this)[8]*x0;
    T x49 = (*this)[7]*x0;
    T x50 = (*this)[10]*(*this)[1];
    T x51 = (*this)[9]*x4;
    T x52 = (*this)[6]*x2;
    T x53 = (*this)[13]*x2;
    T x54 = (*this)[5]*x0;
    T x55 = (*this)[3]*x4;
    T x56 = (*this)[7]*x2;
    T x57 = (*this)[10]*x0;
    T x58 = (*this)[14]*x4;
    T x59 = (*this)[9]*x2;
    res[0]=0;
    res[1]=-(*this)[10]*x6 + (*this)[10]*x7 - (*this)[10]*x8 - (*this)[12]*x9 - (*this)[13]*x10 - (*this)[14]*x11 + (*this)[1]*x21 + (*this)[1]*x22 + (*this)[1]*x23 - (*this)[2]*x17 - (*this)[3]*x12 + (*this)[3]*x18 + (*this)[4]*x13 - (*this)[4]*x16 - (*this)[5]*x1 + (*this)[6]*x25 - (*this)[6]*x3 - (*this)[7]*x26 + (*this)[7]*x27 - (*this)[7]*x5 - (*this)[8]*x19 - (*this)[9]*x20 - (*this)[9]*x24 + x14*x15;
    res[2]=(*this)[10]*x25 - (*this)[10]*x3 + (*this)[13]*x13 + (*this)[14]*x12 - (*this)[1]*x17 + (*this)[1]*x36 + (*this)[2]*x22 + (*this)[3]*x11 - (*this)[4]*x10 + (*this)[4]*x37 - (*this)[6]*x6 + (*this)[6]*x7 + (*this)[7]*x20 + (*this)[7]*x24 + (*this)[9]*x26 + (*this)[9]*x5 + A[0]*x28 + A[0]*x29 + A[0]*x30 + A[0]*x31 + A[0]*x32 + A[0]*x33 + A[0]*x34 + A[0]*x35 - A[0]*x38 - A[0]*x39 - A[0]*x40 - A[0]*x41 - A[0]*x42 - A[0]*x43 - A[0]*x44 - A[0]*x45;
    res[3]=-(*this)[0]*x25 + (*this)[10]*x1 - (*this)[11]*x37 + (*this)[12]*x16 + (*this)[14]*x36 - (*this)[1]*x12 + (*this)[1]*x18 + (*this)[3]*x21 + (*this)[3]*x23 + (*this)[4]*x9 + (*this)[5]*x6 + (*this)[5]*x8 - (*this)[7]*x19 + (*this)[8]*x27 + A[1]*x28 + A[1]*x29 - A[1]*x30 + A[1]*x31 + A[1]*x32 - A[1]*x33 - A[1]*x34 - A[1]*x35 - A[1]*x38 + A[1]*x39 - A[1]*x40 + A[1]*x41 - A[1]*x42 + A[1]*x43 - A[1]*x44 + A[1]*x45 + x4*x46 + x4*x47;
    res[4]=-(*this)[0]*x27 + (*this)[10]*x48 + (*this)[12]*x18 + (*this)[13]*x17 - (*this)[15]*x7 + (*this)[15]*x8 + (*this)[1]*x13 - (*this)[1]*x16 + (*this)[2]*x10 - (*this)[3]*x9 + (*this)[4]*x21 + (*this)[4]*x22 + (*this)[5]*x49 + (*this)[8]*x3 + A[2]*x28 + A[2]*x29 - A[2]*x30 + A[2]*x31 + A[2]*x32 - A[2]*x33 - A[2]*x34 - A[2]*x35 + A[2]*x38 - A[2]*x39 + A[2]*x40 - A[2]*x41 + A[2]*x42 - A[2]*x43 + A[2]*x44 - A[2]*x45 + x2*x46 + x2*x47;
    res[5]=(*this)[0]*x36 + (*this)[10]*x12 - (*this)[10]*x18 + (*this)[12]*x48 + (*this)[13]*x26 - (*this)[13]*x27 + (*this)[14]*x25 - (*this)[14]*x3 - (*this)[15]*x9 + (*this)[1]*x1 + (*this)[1]*x51 + (*this)[2]*x52 - (*this)[3]*x6 + (*this)[3]*x7 - (*this)[3]*x8 + (*this)[4]*x20 - (*this)[4]*x49 + (*this)[5]*x21 + (*this)[5]*x23 + (*this)[6]*x11 - (*this)[7]*x10 + (*this)[7]*x37 + (*this)[9]*x13 - x2*x50;
    res[6]=(*this)[0]*x18 - (*this)[10]*x17 + (*this)[10]*x36 + (*this)[12]*x27 - (*this)[12]*x5 + (*this)[14]*x51 - (*this)[15]*x10 + (*this)[15]*x37 - (*this)[1]*x25 + (*this)[1]*x3 - (*this)[2]*x7 + (*this)[3]*x54 - (*this)[4]*x19 - (*this)[4]*x56 - (*this)[5]*x11 + (*this)[6]*x21 + (*this)[6]*x22 + (*this)[6]*x23 + (*this)[7]*x55 + (*this)[7]*x9 - (*this)[8]*x13 + (*this)[8]*x16 + (*this)[9]*x53 + x0*x50;
    res[7]=(*this)[0]*x13 + (*this)[10]*x53 + (*this)[10]*x58 - (*this)[11]*x6 + (*this)[11]*x7 - (*this)[11]*x8 - (*this)[12]*x25 + (*this)[12]*x57 - (*this)[13]*x1 + (*this)[1]*x26 - (*this)[1]*x27 + (*this)[1]*x5 - (*this)[2]*x20 + (*this)[3]*x19 + (*this)[4]*x52 + (*this)[4]*x54 - (*this)[5]*x37 - (*this)[6]*x55 + (*this)[7]*x21 + (*this)[7]*x22 + (*this)[7]*x23 + (*this)[8]*x18 + (*this)[9]*x17 - (*this)[9]*x36;
    res[8]=-(*this)[10]*x10 + (*this)[10]*x37 + (*this)[11]*x1 - (*this)[12]*x54 - (*this)[13]*x7 + (*this)[13]*x8 - (*this)[14]*x24 - (*this)[15]*x17 + (*this)[15]*x36 - (*this)[1]*(*this)[6]*x4 + (*this)[1]*x19 + (*this)[1]*x56 + (*this)[2]*x59 - (*this)[3]*x27 + (*this)[3]*x5 - (*this)[4]*x3 - (*this)[4]*x57 - (*this)[6]*x13 - (*this)[7]*x12 + (*this)[7]*x18 + (*this)[8]*x21 + (*this)[8]*x22 + (*this)[8]*x23 + (*this)[9]*x11;
    res[9]=-(*this)[0]*x37 - (*this)[10]*(*this)[4]*x2 + (*this)[10]*x55 + (*this)[10]*x9 - (*this)[11]*x25 + (*this)[11]*x3 - (*this)[12]*x6 + (*this)[12]*x7 - (*this)[12]*x8 - (*this)[13]*x52 + (*this)[15]*x18 + (*this)[1]*x20 + (*this)[1]*x24 - (*this)[1]*x49 - (*this)[2]*x26 + (*this)[3]*x48 + (*this)[4]*x1 - (*this)[5]*x16 - (*this)[6]*x58 + (*this)[7]*x17 - (*this)[7]*x36 + (*this)[9]*x21 + (*this)[9]*x22 + (*this)[9]*x23;
    res[10]=(*this)[10]*x21 + (*this)[10]*x22 + (*this)[10]*x23 - (*this)[11]*x27 + (*this)[11]*x5 - (*this)[12]*x49 + (*this)[15]*x13 - (*this)[15]*x16 + (*this)[1]*x6 - (*this)[1]*x7 + (*this)[1]*x8 + (*this)[2]*x3 - (*this)[3]*x1 - (*this)[3]*x51 + (*this)[4]*x48 + (*this)[4]*x59 + (*this)[5]*x12 - (*this)[5]*x18 - (*this)[6]*x17 + (*this)[6]*x36 - (*this)[7]*x53 - (*this)[7]*x58 + (*this)[8]*x10 - (*this)[8]*x37;
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
    res[2]=-x0*((*this)[0]*(*this)[5] - (*this)[10]*(*this)[6] - (*this)[11]*(*this)[12] + (*this)[13]*(*this)[4] - (*this)[14]*(*this)[3] + (*this)[15]*(*this)[8] + (*this)[1]*(*this)[2] + (*this)[7]*(*this)[9]);
    res[3]=-x0*((*this)[0]*(*this)[6] + (*this)[10]*(*this)[5] - (*this)[11]*(*this)[13] - (*this)[12]*(*this)[4] + (*this)[14]*(*this)[2] + (*this)[15]*(*this)[9] + (*this)[1]*(*this)[3] - (*this)[7]*(*this)[8]);
    res[4]=-x0*((*this)[0]*(*this)[7] + (*this)[10]*(*this)[15] - (*this)[11]*(*this)[14] + (*this)[12]*(*this)[3] - (*this)[13]*(*this)[2] + (*this)[1]*(*this)[4] - (*this)[5]*(*this)[9] + (*this)[6]*(*this)[8]);
    res[5]=-x0*((*this)[0]*(*this)[2] - (*this)[10]*(*this)[3] - (*this)[11]*(*this)[8] + (*this)[12]*(*this)[15] + (*this)[13]*(*this)[7] - (*this)[14]*(*this)[6] + (*this)[1]*(*this)[5] + (*this)[4]*(*this)[9]);
    res[6]=-x0*((*this)[0]*(*this)[3] + (*this)[10]*(*this)[2] - (*this)[11]*(*this)[9] - (*this)[12]*(*this)[7] + (*this)[13]*(*this)[15] + (*this)[14]*(*this)[5] + (*this)[1]*(*this)[6] - (*this)[4]*(*this)[8]);
    res[7]=-x0*((*this)[0]*(*this)[4] - (*this)[10]*(*this)[11] + (*this)[12]*(*this)[6] - (*this)[13]*(*this)[5] + (*this)[14]*(*this)[15] + (*this)[1]*(*this)[7] - (*this)[2]*(*this)[9] + (*this)[3]*(*this)[8]);
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
    T x1 = A[3]*x0;
    T x2 = A[4]*x0;
    T x3 = (*this)[14]*x0;
    T x4 = A[0]*x0;
    T x5 = (*this)[3]*x0;
    T x6 = A[2]*x0;
    T x7 = 2.0*(*this)[10];
    T x8 = (*this)[11]*x7;
    T x9 = (*this)[12]*x7;
    T x10 = A[3]*x7;
    T x11 = A[5]*x7;
    T x12 = A[1]*x7;
    T x13 = (*this)[3]*A[0];
    T x14 = 2.0*A[3];
    T x15 = (*this)[5]*x14;
    T x16 = 2.0*(*this)[6];
    T x17 = (*this)[11]*x16;
    T x18 = 2.0*(*this)[11];
    T x19 = (*this)[7]*A[5];
    T x20 = A[0]*x18;
    T x21 = 2.0*(*this)[9];
    T x22 = (*this)[11]*x21;
    T x23 = 2.0*(*this)[12];
    T x24 = (*this)[15]*A[0];
    T x25 = (*this)[12]*x16;
    T x26 = (*this)[7]*A[1];
    T x27 = (*this)[12]*x21;
    T x28 = 2.0*(*this)[13];
    T x29 = (*this)[15]*x28;
    T x30 = (*this)[5]*A[2];
    T x31 = (*this)[7]*A[0];
    T x32 = (*this)[8]*x28;
    T x33 = 2.0*(*this)[14];
    T x34 = (*this)[15]*x33;
    T x35 = (*this)[5]*A[1];
    T x36 = (*this)[14]*A[0];
    T x37 = (*this)[8]*x33;
    T x38 = A[3]*x33;
    T x39 = (*this)[15]*(*this)[2];
    T x40 = 2.0*A[4];
    T x41 = (*this)[15]*(*this)[3];
    T x42 = 2.0*(*this)[4];
    T x43 = (*this)[15]*x42;
    T x44 = 2.0*(*this)[1];
    T x45 = (*this)[5]*A[0];
    T x46 = A[1]*x16;
    T x47 = (*this)[7]*x44;
    T x48 = (*this)[8]*A[3];
    T x49 = (*this)[9]*x44;
    T x50 = A[5]*x16;
    T x51 = (*this)[2]*x40;
    T x52 = (*this)[2]*A[2];
    T x53 = 2.0*(*this)[3];
    T x54 = (*this)[5]*A[5];
    T x55 = (*this)[3]*(*this)[7];
    T x56 = (*this)[8]*A[2];
    T x57 = (*this)[5]*A[4];
    T x58 = A[3]*x16;
    T x59 = (*this)[8]*x42;
    T x60 = (*this)[4]*x21;
    T x61 = (*this)[8]*x23;
    T x62 = (*this)[9]*x28;
    T x63 = (*this)[15]*x44;
    T x64 = 2.0*A[2];
    T x65 = 2.0*(*this)[2];
    T x66 = 2.0*x52;
    T x67 = (*this)[2]*x21;
    T x68 = (*this)[8]*x40;
    T x69 = (*this)[3]*x21;
    T x70 = A[1]*x0;
    T x71 = A[5]*x0;
    T x72 = A[2]*x7;
    T x73 = (*this)[1]*A[0];
    T x74 = A[4]*x7;
    T x75 = (*this)[15]*x18;
    T x76 = (*this)[8]*x18;
    T x77 = (*this)[15]*x23;
    T x78 = A[3]*x28;
    T x79 = A[4]*x16;
    T x80 = (*this)[7]*A[4];
    T x81 = (*this)[9]*x33;
    T x82 = (*this)[2]*A[0];
    T x83 = 2.0*x13;
    T x84 = A[2]*x16;
    T x85 = (*this)[12]*x14;
    T x86 = 2.0*A[1];
    T x87 = (*this)[8]*A[1];
    T x88 = (*this)[8]*A[5];
    T x89 = A[2]*x42;
    T x90 = std::pow((*this)[0], 2);
    T x91 = std::pow((*this)[12], 2);
    T x92 = std::pow((*this)[1], 2);
    T x93 = std::pow((*this)[3], 2);
    T x94 = std::pow((*this)[4], 2);
    T x95 = std::pow((*this)[6], 2);
    T x96 = std::pow((*this)[7], 2);
    T x97 = std::pow((*this)[8], 2);
    T x98 = A[5]*x28;
    T x99 = A[4]*x33;
    T x100 = (*this)[11]*x44;
    T x101 = (*this)[3]*x18;
    T x102 = (*this)[11]*x42;
    T x103 = (*this)[12]*x28;
    T x104 = (*this)[12]*A[2];
    T x105 = (*this)[3]*A[4];
    T x106 = (*this)[12]*A[5];
    T x107 = A[2]*x28;
    T x108 = (*this)[2]*A[4];
    T x109 = A[1]*x33;
    T x110 = A[5]*x33;
    T x111 = 2.0*(*this)[15];
    T x112 = (*this)[15]*x21;
    T x113 = A[5]*x44;
    T x114 = (*this)[4]*x44;
    T x115 = (*this)[2]*A[1];
    T x116 = 2.0*(*this)[7];
    T x117 = 2.0*x48;
    T x118 = 2.0*(*this)[8];
    T x119 = std::pow((*this)[10], 2);
    T x120 = std::pow((*this)[11], 2);
    T x121 = std::pow((*this)[13], 2);
    T x122 = std::pow((*this)[14], 2);
    T x123 = std::pow((*this)[15], 2);
    T x124 = std::pow((*this)[2], 2);
    T x125 = std::pow((*this)[5], 2);
    T x126 = std::pow((*this)[9], 2);
    T x127 = (*this)[8]*A[0];
    T x128 = (*this)[12]*x18;
    T x129 = A[0]*x33;
    T x130 = (*this)[12]*x44;
    T x131 = (*this)[2]*A[5];
    T x132 = A[1]*x28;
    T x133 = A[4]*x28;
    T x134 = (*this)[3]*x44;
    T x135 = (*this)[3]*A[1];
    T x136 = 2.0*(*this)[5];
    T x137 = (*this)[15]*(*this)[7];
    T x138 = (*this)[4]*x14;
    T x139 = A[2]*x33;
    T x140 = (*this)[3]*x14;
    res[0]=0;
    res[1]=-(*this)[11]*x15 - (*this)[12]*x1 - (*this)[13]*x10 - (*this)[13]*x2 + (*this)[1]*x11 - (*this)[1]*x46 - (*this)[2]*x12 + (*this)[2]*x4 + (*this)[2]*x50 - (*this)[4]*x58 + (*this)[4]*x6 - (*this)[7]*x51 - (*this)[8]*x20 + (*this)[9]*x38 - A[0]*x60 - A[1]*x22 - A[1]*x29 + A[1]*x5 + A[1]*x59 + A[2]*x25 - A[2]*x34 - A[2]*x47 - A[2]*x8 - A[4]*x17 - A[4]*x37 + A[4]*x49 + A[4]*x9 - A[5]*x27 - A[5]*x3 + A[5]*x32 - A[5]*x43 + x13*x7 - x14*x39 + x14*x55 - x16*x36 - x18*x19 + x21*x52 - x23*x24 - x23*x26 - x28*x30 + x28*x31 + x33*x35 - x40*x41 + x42*x57 - x44*x45 + x44*x48 - x53*x54 - x53*x56;
    res[2]=(*this)[11]*x1 + (*this)[12]*x15 - (*this)[13]*x58 + (*this)[13]*x6 + (*this)[15]*x20 - (*this)[1]*x12 + (*this)[1]*x4 + (*this)[1]*x50 + (*this)[2]*x11 - (*this)[2]*x46 + (*this)[3]*x68 - (*this)[4]*x10 - (*this)[4]*x2 - (*this)[7]*x38 - (*this)[7]*x66 + A[0]*x61 - A[0]*x62 + A[1]*x27 - A[1]*x3 + A[1]*x32 - A[1]*x43 - A[2]*x17 + A[2]*x37 + A[2]*x49 + A[2]*x9 - A[3]*x63 - A[3]*x69 + A[4]*x25 + A[4]*x34 - A[4]*x47 + A[4]*x67 - A[4]*x8 + A[5]*x22 - A[5]*x29 + A[5]*x5 + A[5]*x59 + x13*x16 + x18*x26 + x19*x23 + x28*x57 - x30*x42 + x31*x42 + x33*x54 - x35*x53 - x36*x7 + x41*x64 - x45*x65 + x48*x65;
    res[3]=(*this)[11]*x2 - (*this)[12]*x6 + (*this)[13]*x72 + (*this)[13]*x79 - (*this)[14]*x12 + (*this)[14]*x4 + (*this)[14]*x50 - (*this)[15]*x38 - (*this)[15]*x66 + (*this)[1]*x70 - (*this)[2]*x71 + (*this)[3]*x11 - (*this)[3]*x46 + (*this)[4]*x1 - (*this)[4]*x74 - (*this)[4]*x84 + (*this)[5]*x78 - (*this)[5]*x83 - (*this)[8]*x51 + A[0]*x27 + A[0]*x32 - A[1]*x61 + A[1]*x62 + A[1]*x75 + A[2]*x81 + A[3]*x25 + A[3]*x47 + A[3]*x67 + A[3]*x8 - A[4]*x63 + A[4]*x69 + A[5]*x60 - A[5]*x76 + A[5]*x77 - x16*x82 + x18*x30 - x18*x31 + x19*x28 - x23*x57 + x24*x42 + x26*x42 - x33*x80 + x35*x65 - x44*x54 - x44*x56 + x48*x53 - x55*x64 + x7*x73;
    res[4]=(*this)[11]*x71 + (*this)[12]*x70 + (*this)[13]*x12 - (*this)[13]*x4 - (*this)[13]*x50 + (*this)[14]*x72 + (*this)[14]*x79 - (*this)[15]*x83 - (*this)[1]*x58 + (*this)[1]*x6 + (*this)[2]*x10 + (*this)[2]*x2 + (*this)[3]*x74 + (*this)[3]*x84 + (*this)[4]*x11 - (*this)[4]*x46 + (*this)[5]*x38 + (*this)[5]*x66 + (*this)[7]*x85 - (*this)[7]*x89 + A[0]*x17 + A[0]*x37 - A[0]*x49 + A[0]*x9 + A[1]*x81 - A[2]*x62 + A[2]*x75 - A[3]*x22 + A[3]*x29 - A[3]*x5 + A[4]*x60 + A[4]*x76 - A[4]*x77 - A[5]*x63 - A[5]*x69 - x18*x35 + x19*x33 - x23*x54 - x23*x56 + x28*x80 - x31*x65 + x39*x86 - x42*x45 + x42*x48 + x44*x57 + x44*x87 - x55*x86 - x65*x88;
    res[5]=-(*this)[10]*x70 + (*this)[11]*x98 - (*this)[11]*x99 - (*this)[15]*x1 + (*this)[15]*x74 + (*this)[15]*x84 + (*this)[1]*x107 - (*this)[1]*x109 + (*this)[2]*x110 + (*this)[2]*x85 + (*this)[3]*x113 - (*this)[3]*x78 - (*this)[4]*x38 + (*this)[5]*x11 + (*this)[5]*x117 + (*this)[6]*x71 - (*this)[7]*x10 - (*this)[7]*x2 + (*this)[8]*x79 - (*this)[9]*x58 + (*this)[9]*x6 - A[0]*x119 - A[0]*x120 - A[0]*x121 - A[0]*x122 - A[0]*x123 - A[0]*x124 - A[0]*x125 - A[0]*x126 + A[0]*x90 + A[0]*x91 + A[0]*x92 + A[0]*x93 + A[0]*x94 + A[0]*x95 + A[0]*x96 + A[0]*x97 + A[1]*x102 + A[1]*x103 - A[2]*x101 + A[3]*x100 - A[4]*x114 - A[5]*x112 + x104*x33 + x105*x23 + x106*x42 + x108*x28 - x111*x26 - x115*x53 - x116*x30 + x118*x19 - x16*x35 + x21*x57 + x21*x87 - x42*x52 + x56*x7;
    res[6]=(*this)[10]*x4 + (*this)[11]*x38 + (*this)[14]*x107 - (*this)[15]*x10 - (*this)[15]*x2 - (*this)[2]*x113 + (*this)[2]*x78 + (*this)[3]*x110 + (*this)[3]*x85 - (*this)[3]*x89 + (*this)[4]*x98 - (*this)[4]*x99 + (*this)[5]*A[3]*x21 + (*this)[6]*x11 + (*this)[7]*x1 - (*this)[7]*x74 - (*this)[7]*x84 - (*this)[8]*x6 + (*this)[9]*x72 + (*this)[9]*x79 - A[0]*x102 + A[0]*x103 - A[1]*x119 - A[1]*x120 + A[1]*x121 - A[1]*x122 - A[1]*x123 + A[1]*x124 + A[1]*x125 + A[1]*x126 + A[1]*x90 - A[1]*x91 + A[1]*x92 - A[1]*x93 + A[1]*x94 - A[1]*x95 + A[1]*x96 - A[1]*x97 + A[3]*x114 + A[4]*x100 - x0*x54 - x104*x44 + x105*x28 - x106*x18 - x108*x23 - x111*x30 + x111*x31 + x111*x88 - x118*x57 + x127*x21 - x13*x65 - x16*x45 + x16*x48 + x18*x52 + x19*x21 + x33*x73;
    res[7]=-(*this)[11]*x78 + (*this)[12]*x129 + (*this)[14]*x132 - (*this)[15]*x68 - (*this)[15]*x71 + (*this)[2]*x38 - (*this)[3]*x98 + (*this)[4]*x110 + (*this)[4]*x133 + (*this)[4]*x85 + (*this)[5]*x10 - (*this)[6]*x1 + (*this)[6]*x74 + (*this)[7]*x11 + (*this)[7]*x117 + (*this)[8]*x70 + (*this)[9]*x12 - (*this)[9]*x4 - (*this)[9]*x50 + A[1]*x130 + A[2]*x119 - A[2]*x120 - A[2]*x121 + A[2]*x122 - A[2]*x123 + A[2]*x124 + A[2]*x125 - A[2]*x126 + A[2]*x90 - A[2]*x91 + A[2]*x92 + A[2]*x93 - A[2]*x94 + A[2]*x95 - A[2]*x96 - A[2]*x97 + A[3]*x112 - A[3]*x134 + A[4]*x128 + A[5]*x100 + x0*x57 + x105*x33 + x108*x44 + x111*x35 - x115*x18 - x118*x54 + x127*x7 + x13*x18 - x131*x23 - x135*x42 - x136*x31 - x16*x24 - x16*x26 + x21*x80 - x28*x73 - x42*x82;
    res[8]=-(*this)[10]*x2 + (*this)[11]*x107 - (*this)[11]*x109 + (*this)[12]*x89 - (*this)[15]*x12 + (*this)[15]*x4 + (*this)[15]*x50 - (*this)[1]*x98 + (*this)[1]*x99 + (*this)[3]*x51 - (*this)[4]*x129 - (*this)[6]*x6 + (*this)[8]*A[4]*x21 + (*this)[8]*x11 - (*this)[8]*x46 + (*this)[9]*A[0]*x16 + (*this)[9]*x71 + A[0]*x100 - A[1]*x114 + A[2]*x112 + A[2]*x134 - A[3]*x119 + A[3]*x120 + A[3]*x121 + A[3]*x122 - A[3]*x123 + A[3]*x124 - A[3]*x125 - A[3]*x126 + A[3]*x90 - A[3]*x91 - A[3]*x92 - A[3]*x93 - A[3]*x94 + A[3]*x95 + A[3]*x96 + A[3]*x97 - A[4]*x102 - A[4]*x103 + A[5]*x101 + x0*x26 - x106*x33 + x115*x28 - x116*x56 - x118*x45 - x13*x28 + x131*x42 + x135*x23 - x136*x19 - x137*x40 - x16*x57 - x21*x35 + x23*x82 - x30*x7 + x31*x7 + x33*x52;
    res[9]=(*this)[10]*x1 + (*this)[11]*x129 + (*this)[11]*x138 - (*this)[12]*x78 - (*this)[14]*x98 + (*this)[15]*x70 - (*this)[1]*x38 + (*this)[2]*x140 + (*this)[3]*A[5]*x42 + (*this)[3]*x132 + (*this)[3]*x139 + (*this)[4]*x107 - (*this)[4]*x109 - (*this)[5]*x58 + (*this)[5]*x6 - (*this)[6]*x72 - (*this)[7]*A[2]*x21 + (*this)[7]*x12 - (*this)[7]*x4 - (*this)[7]*x50 - (*this)[8]*x71 + (*this)[9]*x11 - (*this)[9]*x46 + A[0]*x114 + A[1]*x100 - A[4]*x119 + A[4]*x120 - A[4]*x121 + A[4]*x122 - A[4]*x123 - A[4]*x124 + A[4]*x125 + A[4]*x126 + A[4]*x90 + A[4]*x91 - A[4]*x92 + A[4]*x93 - A[4]*x94 - A[4]*x95 + A[4]*x96 - A[4]*x97 - x104*x18 + x106*x44 - x111*x54 - x111*x56 - x115*x23 + x118*x35 - x127*x16 + x13*x23 - x131*x18 + x137*x14 - x21*x45 + x21*x48 + x24*x7 + x28*x82 - x44*x52;
    res[10]=-(*this)[11]*A[0]*x28 - (*this)[11]*x140 + (*this)[12]*A[0]*x42 - (*this)[12]*x38 - (*this)[14]*x133 + (*this)[15]*(*this)[8]*x86 - (*this)[15]*x58 + (*this)[15]*x6 + (*this)[1]*x78 + (*this)[2]*x138 - (*this)[3]*x107 + (*this)[3]*x109 + (*this)[4]*x132 + (*this)[4]*x139 - (*this)[6]*x12 + (*this)[6]*x4 - (*this)[7]*x15 - (*this)[7]*x72 - (*this)[7]*x79 + (*this)[8]*x2 - (*this)[9]*x1 + (*this)[9]*x74 + (*this)[9]*x84 + A[1]*x128 + A[2]*x100 - A[4]*x130 + A[5]*x119 + A[5]*x120 + A[5]*x121 - A[5]*x122 - A[5]*x123 - A[5]*x124 + A[5]*x125 - A[5]*x126 + A[5]*x90 + A[5]*x91 - A[5]*x92 - A[5]*x93 + A[5]*x94 + A[5]*x95 - A[5]*x96 - A[5]*x97 - x0*x35 + x105*x42 + x108*x18 + x111*x57 + x115*x44 + x118*x30 - x118*x31 - x13*x44 - x21*x24 - x21*x26 - x23*x52 + x33*x82 - x45*x7 + x48*x7;
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
    T x0 = 2.0*A[0];
    T x1 = (*this)[0]*x0;
    T x2 = 2.0*A[1];
    T x3 = (*this)[0]*x2;
    T x4 = 2.0*A[2];
    T x5 = (*this)[0]*x4;
    T x6 = 2.0*(*this)[10];
    T x7 = (*this)[12]*A[1];
    T x8 = A[0]*x6;
    T x9 = A[2]*x6;
    T x10 = (*this)[5]*x0;
    T x11 = (*this)[6]*x2;
    T x12 = (*this)[7]*x4;
    T x13 = 2.0*(*this)[9];
    T x14 = A[2]*x13;
    T x15 = 2.0*(*this)[8];
    T x16 = A[2]*x15;
    T x17 = A[1]*x15;
    T x18 = (*this)[14]*A[0];
    T x19 = (*this)[15]*x0;
    T x20 = (*this)[15]*x2;
    T x21 = (*this)[15]*x4;
    T x22 = A[0]*x15;
    T x23 = A[1]*x13;
    T x24 = (*this)[6]*x4;
    T x25 = (*this)[7]*x2;
    T x26 = (*this)[5]*x4;
    T x27 = (*this)[7]*x0;
    T x28 = (*this)[5]*x2;
    T x29 = (*this)[6]*x0;
    T x30 = A[1]*x6;
    T x31 = 2.0*x7;
    T x32 = 2.0*x18;
    T x33 = A[0]*x13;
    T x34 = (*this)[13]*x4;
    T x35 = (*this)[14]*x2;
    T x36 = (*this)[11]*(*this)[1];
    T x37 = (*this)[12]*x0;
    T x38 = (*this)[12]*x4;
    T x39 = (*this)[13]*(*this)[2];
    T x40 = (*this)[13]*(*this)[3];
    T x41 = (*this)[14]*x4;
    T x42 = (*this)[1]*x4;
    T x43 = (*this)[1]*(*this)[4];
    T x44 = (*this)[11]*x0;
    T x45 = (*this)[4]*x2;
    T x46 = (*this)[2]*x2;
    T x47 = (*this)[3]*x0;
    T x48 = std::pow((*this)[0], 2);
    T x49 = std::pow((*this)[11], 2);
    T x50 = std::pow((*this)[13], 2);
    T x51 = std::pow((*this)[14], 2);
    T x52 = std::pow((*this)[2], 2);
    T x53 = std::pow((*this)[6], 2);
    T x54 = std::pow((*this)[7], 2);
    T x55 = std::pow((*this)[8], 2);
    T x56 = (*this)[3]*x4;
    T x57 = (*this)[2]*x4;
    T x58 = std::pow((*this)[10], 2);
    T x59 = std::pow((*this)[12], 2);
    T x60 = std::pow((*this)[15], 2);
    T x61 = std::pow((*this)[1], 2);
    T x62 = std::pow((*this)[3], 2);
    T x63 = std::pow((*this)[4], 2);
    T x64 = std::pow((*this)[5], 2);
    T x65 = std::pow((*this)[9], 2);
    res[0]=0;
    res[1]=-(*this)[11]*x10 - (*this)[11]*x11 - (*this)[11]*x12 - (*this)[12]*x1 - (*this)[12]*x14 + (*this)[13]*x16 - (*this)[13]*x3 - (*this)[13]*x8 - (*this)[14]*x17 - (*this)[14]*x5 + (*this)[1]*x22 + (*this)[1]*x23 + (*this)[1]*x9 - (*this)[2]*x19 + (*this)[2]*x24 - (*this)[2]*x25 - (*this)[3]*x20 - (*this)[3]*x26 + (*this)[3]*x27 - (*this)[4]*x21 + (*this)[4]*x28 - (*this)[4]*x29 + x13*x18 + x6*x7;
    res[2]=(*this)[11]*x1 + (*this)[11]*x14 - (*this)[11]*x30 + (*this)[12]*x10 + (*this)[12]*x12 - (*this)[13]*x21 + (*this)[13]*x28 - (*this)[13]*x29 + (*this)[14]*x20 + (*this)[14]*x26 - (*this)[1]*x19 + (*this)[1]*x24 - (*this)[1]*x25 + (*this)[2]*x22 + (*this)[2]*x23 + (*this)[2]*x9 + (*this)[3]*x17 - (*this)[3]*x33 + (*this)[3]*x5 + (*this)[4]*x16 - (*this)[4]*x3 - (*this)[4]*x8 + (*this)[6]*x31 - (*this)[7]*x32;
    res[3]=-(*this)[11]*x16 + (*this)[11]*x3 + (*this)[11]*x8 + (*this)[12]*x21 + (*this)[12]*x29 + (*this)[13]*x10 + (*this)[13]*x11 + (*this)[13]*x12 + (*this)[14]*x24 - (*this)[14]*x25 - (*this)[15]*x32 - (*this)[1]*x20 - (*this)[1]*x26 + (*this)[1]*x27 - (*this)[2]*x17 + (*this)[2]*x33 - (*this)[2]*x5 + (*this)[3]*x22 + (*this)[3]*x23 + (*this)[3]*x9 + (*this)[4]*x1 + (*this)[4]*x14 - (*this)[4]*x30 - (*this)[5]*x31;
    res[4]=(*this)[11]*x17 - (*this)[11]*x33 + (*this)[11]*x5 - (*this)[12]*x26 + (*this)[12]*x27 + (*this)[13]*x19 - (*this)[13]*x24 + (*this)[13]*x25 + (*this)[14]*x11 + (*this)[14]*x12 - (*this)[15]*x31 - (*this)[1]*x21 + (*this)[1]*x28 - (*this)[1]*x29 - (*this)[2]*x16 + (*this)[2]*x3 + (*this)[2]*x8 - (*this)[3]*x1 - (*this)[3]*x14 + (*this)[3]*x30 + (*this)[4]*x22 + (*this)[4]*x23 + (*this)[4]*x9 + (*this)[5]*x32;
    res[5]=(*this)[0]*x24 + (*this)[11]*x34 - (*this)[11]*x35 - (*this)[15]*x1 - (*this)[15]*x14 + (*this)[15]*x30 + (*this)[2]*x37 + (*this)[2]*x41 + (*this)[3]*x31 + (*this)[3]*x42 - (*this)[4]*x32 + (*this)[4]*x38 + (*this)[5]*x22 + (*this)[5]*x23 + (*this)[5]*x9 + (*this)[6]*x17 - (*this)[6]*x33 + (*this)[7]*x16 - (*this)[7]*x3 - (*this)[7]*x8 + x0*x36 - x0*x40 + x2*x39 - x2*x43;
    res[6]=(*this)[0]*x27 + (*this)[11]*x32 - (*this)[11]*x38 + (*this)[15]*x16 - (*this)[15]*x3 - (*this)[15]*x8 - (*this)[2]*x31 - (*this)[2]*x42 + (*this)[3]*x37 + (*this)[3]*x41 + (*this)[4]*x34 - (*this)[4]*x35 - (*this)[5]*x17 + (*this)[5]*x33 - (*this)[5]*x5 + (*this)[6]*x22 + (*this)[6]*x23 + (*this)[6]*x9 + (*this)[7]*x14 - (*this)[7]*x30 + x0*x39 + x0*x43 + x2*x36 + x2*x40;
    res[7]=(*this)[0]*x28 + (*this)[11]*x31 - (*this)[13]*x44 + (*this)[13]*x45 - (*this)[15]*x17 + (*this)[15]*x33 - (*this)[15]*x5 + (*this)[1]*x46 - (*this)[1]*x47 + (*this)[2]*x32 - (*this)[2]*x38 - (*this)[3]*x34 + (*this)[3]*x35 + (*this)[4]*x37 + (*this)[4]*x41 - (*this)[5]*x16 + (*this)[5]*x8 - (*this)[6]*x1 - (*this)[6]*x14 + (*this)[6]*x30 + (*this)[7]*x22 + (*this)[7]*x23 + (*this)[7]*x9 + x36*x4;
    res[8]=(*this)[0]*x14 - (*this)[0]*x30 - (*this)[11]*x45 + (*this)[11]*x56 - (*this)[13]*x31 - (*this)[14]*x38 + (*this)[15]*x24 - (*this)[1]*x34 + (*this)[1]*x35 + (*this)[3]*x46 + (*this)[4]*x57 - (*this)[5]*x12 - (*this)[6]*x28 - (*this)[7]*x20 + (*this)[8]*x9 + (*this)[9]*x17 + A[0]*x48 + A[0]*x49 + A[0]*x50 + A[0]*x51 + A[0]*x52 + A[0]*x53 + A[0]*x54 + A[0]*x55 - A[0]*x58 - A[0]*x59 - A[0]*x60 - A[0]*x61 - A[0]*x62 - A[0]*x63 - A[0]*x64 - A[0]*x65;
    res[9]=-(*this)[0]*x16 + (*this)[0]*x8 - (*this)[11]*x57 - (*this)[13]*x37 - (*this)[14]*x34 + (*this)[15]*x27 - (*this)[1]*x32 + (*this)[1]*x38 + (*this)[2]*x47 + (*this)[4]*x44 + (*this)[4]*x56 - (*this)[5]*x21 - (*this)[6]*x10 - (*this)[7]*x24 + (*this)[9]*x22 + (*this)[9]*x9 + A[1]*x48 + A[1]*x49 - A[1]*x50 + A[1]*x51 - A[1]*x52 - A[1]*x53 + A[1]*x54 - A[1]*x55 - A[1]*x58 + A[1]*x59 - A[1]*x60 - A[1]*x61 + A[1]*x62 - A[1]*x63 + A[1]*x64 + A[1]*x65;
    res[10]=(*this)[0]*x17 - (*this)[0]*x33 + (*this)[11]*x46 - (*this)[12]*x32 + (*this)[13]*(*this)[1]*x0 - (*this)[13]*x35 + (*this)[15]*x28 - (*this)[1]*x31 + (*this)[2]*(*this)[4]*x0 - (*this)[3]*x44 + (*this)[3]*x45 - (*this)[5]*x27 - (*this)[6]*x19 - (*this)[7]*x11 + (*this)[8]*x8 + (*this)[9]*x30 + A[2]*x48 + A[2]*x49 + A[2]*x50 - A[2]*x51 - A[2]*x52 + A[2]*x53 - A[2]*x54 - A[2]*x55 + A[2]*x58 + A[2]*x59 - A[2]*x60 - A[2]*x61 - A[2]*x62 + A[2]*x63 + A[2]*x64 - A[2]*x65;
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
    T x4 = 2.0*A[2];
    T x5 = (*this)[11]*x4;
    T x6 = 2.0*A[1];
    T x7 = (*this)[10]*x6;
    T x8 = 2.0*A[0];
    T x9 = (*this)[3]*x8;
    T x10 = (*this)[11]*x8;
    T x11 = (*this)[11]*x6;
    T x12 = (*this)[12]*x8;
    T x13 = (*this)[12]*x4;
    T x14 = (*this)[12]*x6;
    T x15 = (*this)[15]*x6;
    T x16 = (*this)[5]*x4;
    T x17 = (*this)[7]*x8;
    T x18 = (*this)[14]*x4;
    T x19 = (*this)[5]*x6;
    T x20 = (*this)[14]*x8;
    T x21 = (*this)[5]*x8;
    T x22 = (*this)[6]*x6;
    T x23 = (*this)[7]*x4;
    T x24 = (*this)[9]*x4;
    T x25 = (*this)[3]*x4;
    T x26 = (*this)[8]*x6;
    T x27 = (*this)[9]*x8;
    T x28 = (*this)[10]*x4;
    T x29 = (*this)[10]*x8;
    T x30 = (*this)[13]*x8;
    T x31 = (*this)[9]*x6;
    T x32 = (*this)[15]*x4;
    T x33 = (*this)[15]*x8;
    T x34 = (*this)[1]*x4;
    T x35 = (*this)[2]*x8;
    T x36 = (*this)[4]*x4;
    T x37 = (*this)[7]*x6;
    T x38 = std::pow((*this)[0], 2);
    T x39 = std::pow((*this)[12], 2);
    T x40 = std::pow((*this)[1], 2);
    T x41 = std::pow((*this)[3], 2);
    T x42 = std::pow((*this)[4], 2);
    T x43 = std::pow((*this)[6], 2);
    T x44 = std::pow((*this)[7], 2);
    T x45 = std::pow((*this)[8], 2);
    T x46 = (*this)[14]*x6;
    T x47 = (*this)[3]*x6;
    T x48 = std::pow((*this)[10], 2);
    T x49 = std::pow((*this)[11], 2);
    T x50 = std::pow((*this)[13], 2);
    T x51 = std::pow((*this)[14], 2);
    T x52 = std::pow((*this)[15], 2);
    T x53 = std::pow((*this)[2], 2);
    T x54 = std::pow((*this)[5], 2);
    T x55 = std::pow((*this)[9], 2);
    T x56 = (*this)[13]*x6;
    T x57 = (*this)[1]*(*this)[4];
    res[0]=0;
    res[1]=-(*this)[10]*x5 + (*this)[10]*x9 - (*this)[13]*x15 - (*this)[13]*x16 + (*this)[13]*x17 + (*this)[14]*x19 - (*this)[15]*x12 - (*this)[15]*x18 - (*this)[1]*x21 - (*this)[1]*x22 - (*this)[1]*x23 + (*this)[2]*x1 + (*this)[2]*x24 - (*this)[2]*x7 + (*this)[3]*x2 + (*this)[4]*x26 - (*this)[4]*x27 + (*this)[4]*x3 + (*this)[6]*x13 - (*this)[6]*x20 - (*this)[7]*x14 - (*this)[8]*x10 - (*this)[8]*x25 - (*this)[9]*x11;
    res[2]=(*this)[10]*x13 - (*this)[10]*x20 + (*this)[13]*x26 - (*this)[13]*x27 + (*this)[13]*x3 - (*this)[14]*x2 + (*this)[15]*x10 + (*this)[15]*x25 + (*this)[1]*x1 + (*this)[1]*x24 - (*this)[1]*x7 - (*this)[2]*x21 - (*this)[2]*x22 - (*this)[2]*x23 - (*this)[3]*x19 - (*this)[4]*x15 - (*this)[4]*x16 + (*this)[4]*x17 - (*this)[6]*x5 + (*this)[6]*x9 + (*this)[7]*x11 + (*this)[8]*x12 + (*this)[8]*x18 + (*this)[9]*x14;
    res[3]=-(*this)[11]*x17 - (*this)[12]*x26 - (*this)[12]*x3 + (*this)[13]*x28 + (*this)[13]*x31 + (*this)[14]*x1 + (*this)[14]*x24 - (*this)[14]*x7 + (*this)[15]*x11 + (*this)[1]*x2 + (*this)[1]*x29 + (*this)[2]*x19 - (*this)[2]*x32 - (*this)[3]*x22 - (*this)[3]*x23 + (*this)[4]*x33 + (*this)[4]*x37 + (*this)[5]*x5 - (*this)[5]*x9 - (*this)[6]*x35 - (*this)[6]*x36 + (*this)[8]*x30 - (*this)[8]*x34 + (*this)[9]*x12;
    res[4]=(*this)[10]*x12 + (*this)[10]*x18 - (*this)[11]*x19 + (*this)[12]*x2 - (*this)[13]*x1 - (*this)[13]*x24 + (*this)[13]*x7 + (*this)[14]*x31 + (*this)[15]*x5 - (*this)[15]*x9 + (*this)[1]*x26 - (*this)[1]*x27 + (*this)[1]*x3 + (*this)[2]*x15 + (*this)[2]*x16 - (*this)[2]*x17 - (*this)[3]*x37 - (*this)[4]*x21 - (*this)[4]*x22 - (*this)[4]*x23 + (*this)[6]*x10 + (*this)[6]*x25 - (*this)[8]*x13 + (*this)[8]*x20;
    res[5]=-(*this)[10]*x2 + (*this)[13]*x14 + (*this)[13]*x34 + (*this)[14]*x13 - (*this)[1]*x46 - (*this)[2]*x36 - (*this)[2]*x47 - (*this)[3]*x5 + (*this)[4]*x11 - (*this)[6]*x19 + (*this)[6]*x32 - (*this)[7]*x15 - (*this)[7]*x16 + (*this)[8]*x28 + (*this)[9]*x26 + (*this)[9]*x3 + A[0]*x38 + A[0]*x39 + A[0]*x40 + A[0]*x41 + A[0]*x42 + A[0]*x43 + A[0]*x44 + A[0]*x45 - A[0]*x48 - A[0]*x49 - A[0]*x50 - A[0]*x51 - A[0]*x52 - A[0]*x53 - A[0]*x54 - A[0]*x55;
    res[6]=(*this)[10]*x1 + (*this)[10]*x24 + (*this)[13]*x12 + (*this)[13]*x18 - (*this)[15]*x16 + (*this)[15]*x17 - (*this)[1]*x13 + (*this)[1]*x20 + (*this)[2]*x5 - (*this)[2]*x9 - (*this)[4]*x10 - (*this)[4]*x25 - (*this)[6]*x21 - (*this)[6]*x23 + (*this)[8]*x27 - (*this)[8]*x3 + A[1]*x38 - A[1]*x39 + A[1]*x40 - A[1]*x41 + A[1]*x42 - A[1]*x43 + A[1]*x44 - A[1]*x45 - A[1]*x48 - A[1]*x49 + A[1]*x50 - A[1]*x51 - A[1]*x52 + A[1]*x53 + A[1]*x54 + A[1]*x55;
    res[7]=(*this)[11]*x9 + (*this)[13]*x46 + (*this)[14]*x12 + (*this)[15]*x19 + (*this)[1]*x14 - (*this)[1]*x30 - (*this)[2]*x11 - (*this)[4]*x35 - (*this)[4]*x47 - (*this)[5]*x17 - (*this)[6]*x33 - (*this)[7]*x22 + (*this)[8]*x2 + (*this)[8]*x29 - (*this)[9]*x1 + (*this)[9]*x7 + A[2]*x38 - A[2]*x39 + A[2]*x40 + A[2]*x41 - A[2]*x42 + A[2]*x43 - A[2]*x44 - A[2]*x45 + A[2]*x48 - A[2]*x49 - A[2]*x50 + A[2]*x51 - A[2]*x52 + A[2]*x53 + A[2]*x54 - A[2]*x55;
    res[8]=-(*this)[10]*x16 + (*this)[10]*x17 + (*this)[13]*x5 - (*this)[13]*x9 - (*this)[14]*x11 + (*this)[15]*x1 + (*this)[15]*x24 - (*this)[15]*x7 + (*this)[1]*x10 + (*this)[1]*x25 + (*this)[2]*x12 + (*this)[2]*x18 + (*this)[2]*x56 + (*this)[3]*x14 + (*this)[4]*x13 - (*this)[4]*x20 - (*this)[6]*x26 + (*this)[6]*x27 - (*this)[6]*x3 + (*this)[7]*x2 - (*this)[8]*x21 - (*this)[8]*x23 - (*this)[9]*x19 - x57*x6;
    res[9]=-(*this)[11]*x13 + (*this)[12]*x9 + (*this)[13]*x36 + (*this)[13]*x47 + (*this)[14]*x10 + (*this)[15]*x2 + (*this)[15]*x29 + (*this)[1]*x11 - (*this)[2]*x14 + (*this)[2]*x30 - (*this)[2]*x34 + (*this)[3]*x18 - (*this)[4]*x46 + (*this)[5]*x3 - (*this)[6]*(*this)[8]*x8 - (*this)[6]*x28 - (*this)[7]*x1 - (*this)[7]*x24 + (*this)[7]*x7 + (*this)[8]*x19 - (*this)[8]*x32 - (*this)[9]*x21 - (*this)[9]*x22 + x57*x8;
    res[10]=-(*this)[10]*x21 - (*this)[10]*x23 + (*this)[12]*x11 - (*this)[13]*x10 - (*this)[13]*x25 + (*this)[15]*x26 - (*this)[15]*x27 + (*this)[15]*x3 + (*this)[1]*(*this)[2]*x6 + (*this)[1]*x5 - (*this)[1]*x9 - (*this)[2]*x13 + (*this)[2]*x20 + (*this)[3]*x46 + (*this)[4]*x12 + (*this)[4]*x18 + (*this)[4]*x56 - (*this)[5]*x2 + (*this)[6]*x1 + (*this)[6]*x24 - (*this)[6]*x7 - (*this)[7]*x31 + (*this)[8]*x16 - (*this)[8]*x17;
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
    T x2 = 2.0*A[1];
    T x3 = (*this)[12]*x2;
    T x4 = 2.0*A[2];
    T x5 = (*this)[0]*x4;
    T x6 = 2.0*A[3];
    T x7 = (*this)[0]*x6;
    T x8 = 2.0*(*this)[10];
    T x9 = (*this)[12]*A[2];
    T x10 = A[1]*x8;
    T x11 = A[3]*x8;
    T x12 = (*this)[10]*x0;
    T x13 = 2.0*(*this)[11];
    T x14 = (*this)[5]*A[1];
    T x15 = (*this)[6]*A[2];
    T x16 = (*this)[7]*A[3];
    T x17 = (*this)[12]*x0;
    T x18 = (*this)[12]*x6;
    T x19 = (*this)[13]*x0;
    T x20 = (*this)[13]*x6;
    T x21 = (*this)[7]*x0;
    T x22 = (*this)[14]*x4;
    T x23 = (*this)[14]*x2;
    T x24 = (*this)[15]*x0;
    T x25 = (*this)[15]*x2;
    T x26 = (*this)[15]*x4;
    T x27 = (*this)[4]*x6;
    T x28 = (*this)[1]*x2;
    T x29 = (*this)[1]*x4;
    T x30 = (*this)[6]*x6;
    T x31 = (*this)[7]*x4;
    T x32 = (*this)[2]*x0;
    T x33 = (*this)[5]*x6;
    T x34 = (*this)[7]*x2;
    T x35 = (*this)[3]*x0;
    T x36 = (*this)[4]*x4;
    T x37 = (*this)[4]*x2;
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
    T x55 = A[2]*x8;
    T x56 = A[1]*x13;
    T x57 = A[2]*x13;
    T x58 = A[3]*x13;
    T x59 = 2.0*x9;
    T x60 = (*this)[2]*(*this)[3];
    T x61 = 2.0*(*this)[5];
    T x62 = (*this)[8]*(*this)[9];
    T x63 = 2.0*(*this)[6];
    T x64 = (*this)[0]*x2;
    T x65 = 2.0*(*this)[7];
    T x66 = 2.0*(*this)[1];
    res[0]=-(*this)[0]*x3 + (*this)[11]*x1 - (*this)[13]*x10 - (*this)[13]*x5 - (*this)[14]*x21 - (*this)[14]*x7 + (*this)[15]*x27 + (*this)[1]*x11 + (*this)[1]*x24 + (*this)[2]*x25 - (*this)[2]*x30 + (*this)[2]*x31 + (*this)[3]*x26 + (*this)[3]*x33 - (*this)[3]*x34 + (*this)[4]*x12 - (*this)[5]*x17 - (*this)[5]*x36 - (*this)[6]*x19 + (*this)[6]*x37 + (*this)[8]*x20 - (*this)[8]*x22 + (*this)[8]*x28 + (*this)[8]*x32 - (*this)[9]*x18 + (*this)[9]*x23 + (*this)[9]*x29 + (*this)[9]*x35 + x13*x14 + x13*x15 + x13*x16 + x8*x9;
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
    res[11]=-(*this)[12]*x56 - (*this)[13]*x37 - (*this)[13]*x57 - (*this)[14]*x58 + (*this)[15]*x11 + (*this)[1]*x27 + (*this)[2]*x20 - (*this)[2]*x22 + (*this)[2]*x28 - (*this)[3]*x18 + (*this)[3]*x23 + (*this)[3]*x29 + (*this)[4]*x59 - (*this)[5]*x55 + (*this)[6]*x10 + (*this)[8]*x25 - (*this)[8]*x30 + (*this)[8]*x31 + (*this)[9]*x26 + (*this)[9]*x33 - (*this)[9]*x34 + A[0]*x38 + A[0]*x39 + A[0]*x40 + A[0]*x41 + A[0]*x42 + A[0]*x43 + A[0]*x44 + A[0]*x45 + A[0]*x46 + A[0]*x47 + A[0]*x48 + A[0]*x49 + A[0]*x50 + A[0]*x51 + A[0]*x52 + A[0]*x53 + x14*x54 + x15*x54 + x16*x54;
    res[12]=-(*this)[0]*x55 + (*this)[11]*x17 - (*this)[13]*x59 - (*this)[14]*x18 + (*this)[14]*x29 + (*this)[14]*x35 - (*this)[15]*x30 - (*this)[1]*x20 - (*this)[1]*x32 - (*this)[2]*x27 - (*this)[3]*x58 - (*this)[4]*x19 + (*this)[4]*x57 + (*this)[5]*x1 - (*this)[6]*x12 + (*this)[7]*x26 + (*this)[8]*x11 + (*this)[8]*x24 + (*this)[9]*x21 + (*this)[9]*x7 + A[1]*x38 - A[1]*x39 - A[1]*x40 - A[1]*x41 + A[1]*x42 + A[1]*x43 + A[1]*x44 - A[1]*x45 - A[1]*x46 + A[1]*x47 + A[1]*x48 + A[1]*x49 - A[1]*x50 - A[1]*x51 + A[1]*x52 - A[1]*x53 + x15*x61 + x16*x61 - x4*x60 + x4*x62;
    res[13]=(*this)[0]*x10 + (*this)[11]*x19 - (*this)[13]*x3 - (*this)[14]*x20 - (*this)[14]*x32 + (*this)[15]*x33 + (*this)[1]*x18 - (*this)[1]*x23 - (*this)[1]*x35 + (*this)[2]*x58 - (*this)[3]*x27 + (*this)[4]*x17 - (*this)[4]*x56 + (*this)[5]*x12 + (*this)[6]*x1 - (*this)[7]*x25 - (*this)[8]*x21 - (*this)[8]*x7 + (*this)[9]*x11 + (*this)[9]*x24 + A[2]*x38 - A[2]*x39 - A[2]*x40 + A[2]*x41 - A[2]*x42 + A[2]*x43 + A[2]*x44 - A[2]*x45 + A[2]*x46 - A[2]*x47 + A[2]*x48 - A[2]*x49 + A[2]*x50 - A[2]*x51 - A[2]*x52 + A[2]*x53 + x14*x63 + x16*x63 - x2*x60 + x2*x62;
    res[14]=(*this)[11]*(*this)[14]*x0 - (*this)[12]*x23 - (*this)[12]*x35 - (*this)[13]*x22 + (*this)[13]*x28 + (*this)[13]*x32 + (*this)[15]*x12 - (*this)[1]*(*this)[4]*x0 - (*this)[1]*x59 - (*this)[2]*x37 - (*this)[2]*x57 - (*this)[3]*x36 + (*this)[3]*x56 - (*this)[5]*(*this)[9]*x0 - (*this)[5]*x26 + (*this)[6]*(*this)[8]*x0 + (*this)[6]*x25 + (*this)[7]*x1 + (*this)[8]*x10 + (*this)[8]*x5 + (*this)[9]*x55 - (*this)[9]*x64 + A[3]*x38 + A[3]*x39 - A[3]*x40 + A[3]*x41 + A[3]*x42 - A[3]*x43 + A[3]*x44 - A[3]*x45 + A[3]*x46 + A[3]*x47 - A[3]*x48 - A[3]*x49 - A[3]*x50 + A[3]*x51 - A[3]*x52 - A[3]*x53 + x14*x65 + x15*x65;
    res[15]=-(*this)[0]*x27 + (*this)[11]*x11 + (*this)[11]*x24 - (*this)[12]*x25 - (*this)[13]*x26 + (*this)[13]*x34 - (*this)[14]*(*this)[15]*x6 - (*this)[14]*x12 - (*this)[1]*x1 - (*this)[2]*(*this)[9]*x6 + (*this)[2]*x55 - (*this)[2]*x64 + (*this)[3]*(*this)[8]*x6 - (*this)[3]*x10 - (*this)[3]*x5 - (*this)[4]*x21 - (*this)[5]*x20 + (*this)[5]*x22 - (*this)[5]*x32 + (*this)[6]*x18 - (*this)[6]*x23 - (*this)[6]*x35 - (*this)[7]*x59 - (*this)[8]*x17 - (*this)[8]*x36 + (*this)[8]*x56 - (*this)[9]*x19 + (*this)[9]*x37 + (*this)[9]*x57 - x14*x66 - x15*x66 - x16*x66;
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
    T x0 = 2.0*A[0];
    T x1 = (*this)[12]*x0;
    T x2 = 2.0*A[1];
    T x3 = (*this)[0]*x2;
    T x4 = 2.0*A[2];
    T x5 = (*this)[0]*x4;
    T x6 = 2.0*(*this)[10];
    T x7 = (*this)[12]*A[1];
    T x8 = A[0]*x6;
    T x9 = A[2]*x6;
    T x10 = 2.0*(*this)[11];
    T x11 = (*this)[5]*A[0];
    T x12 = (*this)[6]*A[1];
    T x13 = (*this)[7]*A[2];
    T x14 = (*this)[12]*x4;
    T x15 = (*this)[13]*x4;
    T x16 = (*this)[14]*x2;
    T x17 = (*this)[14]*x0;
    T x18 = (*this)[15]*x0;
    T x19 = (*this)[15]*x2;
    T x20 = (*this)[4]*x4;
    T x21 = (*this)[1]*x0;
    T x22 = (*this)[1]*x2;
    T x23 = (*this)[6]*x4;
    T x24 = (*this)[7]*x2;
    T x25 = (*this)[5]*x4;
    T x26 = (*this)[7]*x0;
    T x27 = (*this)[4]*x2;
    T x28 = (*this)[4]*x0;
    T x29 = 2.0*(*this)[0];
    T x30 = A[1]*x6;
    T x31 = A[0]*x10;
    T x32 = A[1]*x10;
    T x33 = A[2]*x10;
    T x34 = 2.0*x7;
    T x35 = std::pow((*this)[0], 2);
    T x36 = std::pow((*this)[13], 2);
    T x37 = std::pow((*this)[14], 2);
    T x38 = std::pow((*this)[15], 2);
    T x39 = std::pow((*this)[3], 2);
    T x40 = std::pow((*this)[4], 2);
    T x41 = std::pow((*this)[5], 2);
    T x42 = std::pow((*this)[8], 2);
    T x43 = (*this)[2]*(*this)[3];
    T x44 = 2.0*(*this)[5];
    T x45 = (*this)[8]*(*this)[9];
    T x46 = std::pow((*this)[10], 2);
    T x47 = std::pow((*this)[11], 2);
    T x48 = std::pow((*this)[12], 2);
    T x49 = std::pow((*this)[1], 2);
    T x50 = std::pow((*this)[2], 2);
    T x51 = std::pow((*this)[6], 2);
    T x52 = std::pow((*this)[7], 2);
    T x53 = std::pow((*this)[9], 2);
    T x54 = 2.0*(*this)[6];
    T x55 = (*this)[0]*x0;
    T x56 = 2.0*(*this)[7];
    T x57 = 2.0*(*this)[1];
    res[0]=-(*this)[0]*x1 - (*this)[13]*x3 - (*this)[13]*x8 - (*this)[14]*x5 + (*this)[15]*x20 + (*this)[1]*x9 + (*this)[2]*x18 - (*this)[2]*x23 + (*this)[2]*x24 + (*this)[3]*x19 + (*this)[3]*x25 - (*this)[3]*x26 - (*this)[5]*x27 + (*this)[6]*x28 + (*this)[8]*x15 - (*this)[8]*x16 + (*this)[8]*x21 - (*this)[9]*x14 + (*this)[9]*x17 + (*this)[9]*x22 + x10*x11 + x10*x12 + x10*x13 + x6*x7;
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
    res[11]=-(*this)[12]*x31 - (*this)[13]*x28 - (*this)[13]*x32 - (*this)[14]*x33 + (*this)[15]*x9 + (*this)[1]*x20 + (*this)[2]*x15 - (*this)[2]*x16 + (*this)[2]*x21 - (*this)[3]*x14 + (*this)[3]*x17 + (*this)[3]*x22 + (*this)[4]*x34 - (*this)[5]*x30 + (*this)[6]*x8 + (*this)[8]*x18 - (*this)[8]*x23 + (*this)[8]*x24 + (*this)[9]*x19 + (*this)[9]*x25 - (*this)[9]*x26 + x11*x29 + x12*x29 + x13*x29;
    res[12]=-(*this)[0]*x30 - (*this)[13]*x34 - (*this)[14]*x14 + (*this)[14]*x22 - (*this)[15]*x23 - (*this)[1]*x15 - (*this)[2]*x20 - (*this)[3]*x33 + (*this)[4]*x32 + (*this)[7]*x19 + (*this)[8]*x9 + (*this)[9]*x5 + A[0]*x35 + A[0]*x36 + A[0]*x37 + A[0]*x38 + A[0]*x39 + A[0]*x40 + A[0]*x41 + A[0]*x42 - A[0]*x46 - A[0]*x47 - A[0]*x48 - A[0]*x49 - A[0]*x50 - A[0]*x51 - A[0]*x52 - A[0]*x53 + x12*x44 + x13*x44 - x2*x43 + x2*x45;
    res[13]=(*this)[0]*x8 - (*this)[13]*x1 - (*this)[14]*x15 + (*this)[15]*x25 + (*this)[1]*x14 - (*this)[1]*x17 + (*this)[2]*x33 - (*this)[3]*x20 - (*this)[4]*x31 - (*this)[7]*x18 - (*this)[8]*x5 + (*this)[9]*x9 + A[1]*x35 - A[1]*x36 + A[1]*x37 + A[1]*x38 - A[1]*x39 + A[1]*x40 - A[1]*x41 - A[1]*x42 - A[1]*x46 - A[1]*x47 + A[1]*x48 - A[1]*x49 + A[1]*x50 + A[1]*x51 - A[1]*x52 + A[1]*x53 - x0*x43 + x0*x45 + x11*x54 + x13*x54;
    res[14]=-(*this)[12]*x17 - (*this)[13]*x16 + (*this)[13]*x21 - (*this)[1]*x34 - (*this)[2]*x28 - (*this)[2]*x32 - (*this)[3]*x27 + (*this)[3]*x31 - (*this)[5]*x19 + (*this)[6]*x18 + (*this)[8]*x3 + (*this)[8]*x8 + (*this)[9]*x30 - (*this)[9]*x55 + A[2]*x35 + A[2]*x36 - A[2]*x37 + A[2]*x38 + A[2]*x39 - A[2]*x40 - A[2]*x41 - A[2]*x42 + A[2]*x46 - A[2]*x47 + A[2]*x48 - A[2]*x49 + A[2]*x50 - A[2]*x51 + A[2]*x52 - A[2]*x53 + x11*x56 + x12*x56;
    res[15]=-(*this)[0]*x20 + (*this)[11]*x9 - (*this)[12]*x18 - (*this)[13]*x19 + (*this)[13]*x26 - (*this)[14]*(*this)[15]*x4 - (*this)[2]*(*this)[9]*x4 + (*this)[2]*x30 - (*this)[2]*x55 + (*this)[3]*(*this)[8]*x4 - (*this)[3]*x3 - (*this)[3]*x8 - (*this)[5]*x15 + (*this)[5]*x16 + (*this)[6]*x14 - (*this)[6]*x17 - (*this)[7]*x34 - (*this)[8]*x27 + (*this)[8]*x31 + (*this)[9]*x28 + (*this)[9]*x32 - x11*x57 - x12*x57 - x13*x57;
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
    res[15]=-x0*((*this)[0]*(*this)[1] + (*this)[10]*(*this)[14] - (*this)[11]*(*this)[15] + (*this)[12]*(*this)[8] + (*this)[13]*(*this)[9] + (*this)[2]*(*this)[5] + (*this)[3]*(*this)[6] + (*this)[4]*(*this)[7]);
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
    res[0]=-x0*((*this)[0]*(*this)[15] - (*this)[10]*(*this)[7] + (*this)[11]*(*this)[1] + (*this)[12]*(*this)[2] + (*this)[13]*(*this)[3] + (*this)[14]*(*this)[4] - (*this)[5]*(*this)[8] - (*this)[6]*(*this)[9]);
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
    res[11]=-x0*((*this)[0]*(*this)[1] + (*this)[10]*(*this)[14] + (*this)[11]*(*this)[15] + (*this)[12]*(*this)[8] + (*this)[13]*(*this)[9] - (*this)[2]*(*this)[5] - (*this)[3]*(*this)[6] - (*this)[4]*(*this)[7]);
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
    T x8 = 2.0*A[11];
    T x9 = (*this)[0]*(*this)[11];
    T x10 = 2.0*A[12];
    T x11 = (*this)[0]*(*this)[12];
    T x12 = 2.0*A[13];
    T x13 = (*this)[0]*(*this)[13];
    T x14 = 2.0*A[14];
    T x15 = (*this)[0]*(*this)[14];
    T x16 = 2.0*A[15];
    T x17 = (*this)[0]*(*this)[15];
    T x18 = 2.0*(*this)[10];
    T x19 = (*this)[12]*x18;
    T x20 = (*this)[13]*x18;
    T x21 = (*this)[1]*x18;
    T x22 = (*this)[10]*x8;
    T x23 = (*this)[7]*x18;
    T x24 = 2.0*(*this)[11];
    T x25 = (*this)[1]*x24;
    T x26 = (*this)[5]*x24;
    T x27 = (*this)[6]*x24;
    T x28 = (*this)[7]*x24;
    T x29 = (*this)[12]*(*this)[2];
    T x30 = (*this)[12]*(*this)[5];
    T x31 = (*this)[12]*(*this)[9];
    T x32 = (*this)[13]*(*this)[3];
    T x33 = (*this)[13]*(*this)[6];
    T x34 = (*this)[13]*(*this)[8];
    T x35 = (*this)[14]*(*this)[4];
    T x36 = (*this)[14]*(*this)[7];
    T x37 = (*this)[14]*(*this)[8];
    T x38 = (*this)[14]*(*this)[9];
    T x39 = (*this)[15]*(*this)[1];
    T x40 = (*this)[15]*(*this)[2];
    T x41 = (*this)[15]*(*this)[3];
    T x42 = (*this)[15]*(*this)[4];
    T x43 = (*this)[1]*(*this)[8];
    T x44 = (*this)[1]*(*this)[9];
    T x45 = (*this)[2]*(*this)[6];
    T x46 = (*this)[2]*(*this)[7];
    T x47 = (*this)[2]*(*this)[8];
    T x48 = (*this)[3]*(*this)[5];
    T x49 = (*this)[3]*(*this)[7];
    T x50 = (*this)[3]*(*this)[9];
    T x51 = (*this)[4]*(*this)[5];
    T x52 = (*this)[4]*(*this)[6];
    T x53 = (*this)[5]*(*this)[8];
    T x54 = (*this)[6]*(*this)[9];
    T x55 = std::pow((*this)[12], 2);
    T x56 = std::pow((*this)[13], 2);
    T x57 = std::pow((*this)[14], 2);
    T x58 = std::pow((*this)[15], 2);
    T x59 = std::pow((*this)[1], 2);
    T x60 = std::pow((*this)[5], 2);
    T x61 = std::pow((*this)[6], 2);
    T x62 = std::pow((*this)[7], 2);
    T x63 = 2.0*A[8];
    T x64 = 2.0*A[9];
    T x65 = 2.0*A[10];
    T x66 = 2.0*(*this)[0];
    T x67 = A[5]*x66;
    T x68 = (*this)[3]*x66;
    T x69 = (*this)[4]*x66;
    T x70 = (*this)[5]*x66;
    T x71 = (*this)[6]*x66;
    T x72 = A[4]*x66;
    T x73 = (*this)[11]*x18;
    T x74 = A[4]*x18;
    T x75 = (*this)[2]*x18;
    T x76 = (*this)[3]*A[5];
    T x77 = (*this)[5]*A[3];
    T x78 = (*this)[6]*x18;
    T x79 = (*this)[12]*x24;
    T x80 = (*this)[13]*A[3];
    T x81 = A[4]*x24;
    T x82 = A[5]*x24;
    T x83 = (*this)[9]*x24;
    T x84 = 2.0*(*this)[12];
    T x85 = A[5]*x84;
    T x86 = A[4]*x84;
    T x87 = (*this)[4]*A[3];
    T x88 = (*this)[6]*x84;
    T x89 = (*this)[7]*x84;
    T x90 = 2.0*(*this)[13];
    T x91 = (*this)[15]*x90;
    T x92 = (*this)[2]*A[4];
    T x93 = (*this)[4]*x90;
    T x94 = (*this)[5]*x90;
    T x95 = (*this)[7]*A[5];
    T x96 = 2.0*(*this)[14];
    T x97 = (*this)[15]*x96;
    T x98 = A[3]*x96;
    T x99 = (*this)[3]*A[2];
    T x100 = (*this)[5]*A[6];
    T x101 = (*this)[6]*x96;
    T x102 = 2.0*(*this)[8];
    T x103 = (*this)[15]*x102;
    T x104 = 2.0*(*this)[9];
    T x105 = (*this)[15]*A[3];
    T x106 = 2.0*(*this)[1];
    T x107 = (*this)[2]*A[2];
    T x108 = (*this)[3]*A[3];
    T x109 = (*this)[4]*A[4];
    T x110 = (*this)[5]*A[5];
    T x111 = (*this)[6]*A[6];
    T x112 = (*this)[7]*x106;
    T x113 = (*this)[2]*A[7];
    T x114 = (*this)[3]*A[7];
    T x115 = (*this)[4]*A[6];
    T x116 = (*this)[4]*x104;
    T x117 = (*this)[5]*A[4];
    T x118 = (*this)[6]*A[4];
    T x119 = (*this)[7]*x102;
    T x120 = (*this)[7]*x104;
    T x121 = (*this)[0]*x18;
    T x122 = 2.0*A[7];
    T x123 = 2.0*A[6];
    T x124 = (*this)[14]*x18;
    T x125 = (*this)[4]*x18;
    T x126 = A[4]*x90;
    T x127 = (*this)[9]*x90;
    T x128 = A[1]*x96;
    T x129 = (*this)[14]*x65;
    T x130 = 2.0*(*this)[15];
    T x131 = 2.0*(*this)[7];
    T x132 = A[1]*x106;
    T x133 = (*this)[1]*x65;
    T x134 = 2.0*(*this)[2];
    T x135 = 2.0*(*this)[4];
    T x136 = (*this)[2]*(*this)[9];
    T x137 = 2.0*(*this)[6];
    T x138 = (*this)[3]*(*this)[8];
    T x139 = (*this)[4]*x65;
    T x140 = (*this)[8]*x104;
    T x141 = 2.0*A[5];
    T x142 = (*this)[1]*x66;
    T x143 = (*this)[2]*x66;
    T x144 = (*this)[3]*x18;
    T x145 = (*this)[5]*A[1];
    T x146 = A[1]*x24;
    T x147 = (*this)[15]*x24;
    T x148 = (*this)[4]*A[2];
    T x149 = (*this)[8]*x24;
    T x150 = A[2]*x84;
    T x151 = (*this)[15]*x84;
    T x152 = A[1]*x84;
    T x153 = (*this)[8]*x84;
    T x154 = (*this)[7]*x90;
    T x155 = (*this)[1]*A[2];
    T x156 = A[2]*x131;
    T x157 = (*this)[15]*A[1];
    T x158 = (*this)[1]*x63;
    T x159 = 2.0*(*this)[3];
    T x160 = A[2]*x137;
    T x161 = (*this)[7]*x66;
    T x162 = (*this)[8]*x66;
    T x163 = (*this)[9]*x66;
    T x164 = (*this)[8]*x18;
    T x165 = (*this)[9]*x18;
    T x166 = (*this)[2]*x24;
    T x167 = (*this)[1]*x84;
    T x168 = (*this)[14]*x90;
    T x169 = (*this)[2]*x90;
    T x170 = (*this)[5]*A[8];
    T x171 = A[9]*x106;
    T x172 = 2.0*(*this)[5];
    T x173 = (*this)[7]*x122;
    T x174 = (*this)[8]*x63;
    T x175 = (*this)[9]*x64;
    T x176 = A[1]*x102;
    T x177 = 2.0*A[4];
    T x178 = 2.0*A[3];
    T x179 = (*this)[15]*x18;
    T x180 = (*this)[5]*x18;
    T x181 = (*this)[13]*x24;
    T x182 = (*this)[14]*x24;
    T x183 = A[6]*x84;
    T x184 = A[7]*x84;
    T x185 = (*this)[3]*x84;
    T x186 = (*this)[4]*x84;
    T x187 = (*this)[1]*x90;
    T x188 = (*this)[1]*x96;
    T x189 = (*this)[15]*x122;
    T x190 = (*this)[15]*(*this)[7];
    T x191 = (*this)[9]*x65;
    T x192 = (*this)[2]*(*this)[3];
    T x193 = (*this)[8]*x64;
    T x194 = (*this)[8]*x65;
    T x195 = 2.0*A[2];
    T x196 = (*this)[9]*x63;
    T x197 = (*this)[2]*x96;
    T x198 = (*this)[3]*x96;
    T x199 = (*this)[15]*A[5];
    T x200 = (*this)[4]*A[5];
    T x201 = (*this)[2]*x104;
    T x202 = 2.0*A[1];
    T x203 = (*this)[3]*x24;
    T x204 = (*this)[4]*x24;
    T x205 = (*this)[13]*x84;
    T x206 = (*this)[14]*x84;
    T x207 = (*this)[15]*x65;
    T x208 = (*this)[5]*(*this)[6];
    T x209 = (*this)[7]*x65;
    T x210 = (*this)[15]*(*this)[5];
    T x211 = (*this)[15]*(*this)[6];
    T x212 = (*this)[2]*(*this)[4];
    T x213 = (*this)[3]*(*this)[4];
    T x214 = (*this)[5]*(*this)[7];
    T x215 = (*this)[6]*(*this)[7];
    T x216 = 2.0*A[0];
    T x217 = (*this)[1]*x16;
    T x218 = (*this)[0]*(*this)[5];
    T x219 = (*this)[0]*(*this)[6];
    T x220 = (*this)[0]*x14;
    T x221 = (*this)[12]*x14;
    T x222 = (*this)[12]*x12;
    T x223 = (*this)[12]*x16;
    T x224 = (*this)[13]*x14;
    T x225 = (*this)[13]*(*this)[4];
    T x226 = (*this)[13]*x16;
    T x227 = (*this)[14]*(*this)[2];
    T x228 = (*this)[14]*(*this)[3];
    T x229 = (*this)[15]*(*this)[8];
    T x230 = (*this)[9]*x12;
    T x231 = (*this)[1]*(*this)[2];
    T x232 = (*this)[1]*x12;
    T x233 = (*this)[1]*x14;
    T x234 = (*this)[2]*x16;
    T x235 = (*this)[6]*x16;
    T x236 = (*this)[4]*x16;
    T x237 = (*this)[5]*(*this)[9];
    T x238 = (*this)[6]*(*this)[8];
    T x239 = (*this)[7]*(*this)[8];
    T x240 = (*this)[7]*(*this)[9];
    T x241 = (*this)[11]*x8;
    T x242 = (*this)[0]*(*this)[3];
    T x243 = (*this)[13]*x10;
    T x244 = (*this)[12]*x8;
    T x245 = (*this)[14]*x10;
    T x246 = (*this)[14]*x16;
    T x247 = (*this)[9]*x8;
    T x248 = (*this)[1]*x8;
    T x249 = (*this)[9]*x10;
    T x250 = (*this)[7]*x8;
    T x251 = (*this)[8]*x12;
    T x252 = (*this)[14]*x12;
    T x253 = (*this)[2]*x8;
    res[0]=(*this)[4]*x22 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 - A[0]*x55 - A[0]*x56 - A[0]*x57 - A[0]*x58 - A[0]*x59 + A[0]*x6 - A[0]*x60 - A[0]*x61 - A[0]*x62 + A[0]*x7 - A[12]*x20 + A[12]*x26 + A[13]*x19 + A[13]*x27 + A[14]*x21 + A[14]*x28 + A[15]*x23 - A[15]*x25 - x10*x11 + x10*x38 + x10*x40 + x10*x43 - x10*x49 + x10*x52 - x12*x13 - x12*x37 + x12*x41 + x12*x44 + x12*x46 - x12*x51 - x14*x15 - x14*x31 + x14*x34 + x14*x42 - x14*x45 + x14*x48 - x16*x17 - x16*x29 - x16*x32 - x16*x35 + x16*x53 + x16*x54 - x30*x8 - x33*x8 - x36*x8 + x39*x8 + x47*x8 + x50*x8 + x8*x9;
    res[1]=-(*this)[14]*x81 - (*this)[15]*x74 - (*this)[15]*x85 + (*this)[2]*x67 - (*this)[2]*x98 - (*this)[3]*x86 - (*this)[7]*x72 - (*this)[8]*x82 + A[10]*x21 - A[10]*x28 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[1]*x4 - A[1]*x5 - A[1]*x55 - A[1]*x56 - A[1]*x57 + A[1]*x58 - A[1]*x59 + A[1]*x6 + A[1]*x60 + A[1]*x61 + A[1]*x62 + A[1]*x7 - A[2]*x103 + A[2]*x120 - A[2]*x70 - A[2]*x78 - A[2]*x79 - A[2]*x93 - A[3]*x119 - A[3]*x71 - A[5]*x101 - A[5]*x116 + A[6]*x68 - A[6]*x75 - A[6]*x83 - A[6]*x89 - A[6]*x91 - A[7]*x112 + A[7]*x69 - A[7]*x73 + A[7]*x88 - A[7]*x94 - A[7]*x97 - A[8]*x20 - A[8]*x26 + A[9]*x19 - A[9]*x27 + x100*x96 - x102*x114 + x102*x115 + x102*x118 - x104*x105 + x104*x113 - x104*x117 + x106*x107 + x106*x108 + x106*x109 - x106*x110 - x106*x111 - x11*x63 - x13*x64 - x15*x65 + x18*x76 + x18*x77 - x24*x80 - x31*x65 + x34*x65 - x37*x64 + x38*x63 - x40*x63 - x41*x64 - x42*x65 + x43*x63 + x44*x64 + x45*x65 - x46*x64 - x48*x65 + x49*x63 + x51*x64 - x52*x63 + x84*x87 + x90*x92 + x90*x95 + x96*x99;
    res[2]=(*this)[14]*x86 + (*this)[15]*x82 + (*this)[1]*x126 + (*this)[1]*x67 - (*this)[1]*x98 - (*this)[2]*x132 + (*this)[3]*x128 + (*this)[3]*x81 + (*this)[5]*x129 + (*this)[6]*x133 + (*this)[8]*x139 + (*this)[8]*x74 + (*this)[8]*x85 + (*this)[9]*x72 + A[10]*x68 + A[10]*x75 + A[10]*x83 + A[10]*x89 - A[10]*x91 - A[1]*x103 - A[1]*x120 - A[1]*x70 + A[1]*x78 + A[1]*x79 - A[1]*x93 + A[2]*x0 - A[2]*x1 + A[2]*x2 + A[2]*x3 - A[2]*x4 - A[2]*x5 + A[2]*x55 - A[2]*x56 - A[2]*x57 + A[2]*x58 + A[2]*x59 + A[2]*x6 + A[2]*x60 - A[2]*x61 - A[2]*x62 - A[2]*x7 - A[3]*x121 + A[3]*x140 - A[5]*x124 - A[5]*x127 - A[6]*x21 + A[6]*x28 + A[7]*x19 - A[7]*x27 - A[8]*x125 - A[9]*x112 - A[9]*x69 - A[9]*x73 + A[9]*x88 + A[9]*x94 + A[9]*x97 + x105*x131 + x108*x134 - x110*x134 + x117*x131 - x118*x130 + x122*x13 + x122*x37 + x122*x41 + x122*x44 - x122*x46 - x122*x51 - x123*x15 + x123*x31 + x123*x34 - x123*x42 - x123*x45 - x123*x48 + x135*x92 + x135*x95 + x136*x64 + x137*x76 + x137*x77 + x138*x64 - x24*x87 + x30*x63 - x33*x63 - x36*x63 - x39*x63 + x47*x63 - x50*x63 + x63*x9 + x80*x84;
    res[3]=(*this)[13]*x146 + (*this)[13]*x150 + (*this)[14]*x126 - (*this)[15]*x156 - (*this)[1]*x86 - (*this)[2]*x128 - (*this)[3]*x132 + (*this)[4]*x152 - (*this)[5]*x133 + (*this)[5]*x160 + (*this)[6]*x129 + (*this)[7]*x158 - (*this)[8]*x72 + (*this)[9]*x139 + (*this)[9]*x74 - A[10]*x143 + A[10]*x144 - A[10]*x149 + A[10]*x151 + A[10]*x154 + A[1]*x119 - A[1]*x71 + A[2]*x121 + A[2]*x140 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 + A[3]*x4 - A[3]*x5 - A[3]*x55 + A[3]*x56 - A[3]*x57 + A[3]*x58 + A[3]*x59 - A[3]*x6 - A[3]*x60 + A[3]*x61 - A[3]*x62 + A[3]*x7 + A[5]*x21 - A[5]*x28 - A[6]*x124 + A[6]*x127 + A[6]*x142 + A[6]*x147 - A[6]*x153 + A[7]*x20 + A[7]*x26 + A[8]*x69 + A[8]*x73 + A[8]*x88 + A[8]*x94 - A[8]*x97 - A[9]*x125 + x100*x134 - x104*x157 + x109*x159 - x11*x122 - x111*x159 + x115*x131 + x117*x130 + x118*x131 + x122*x38 - x122*x40 - x122*x43 - x122*x49 - x122*x52 + x134*x99 + x136*x63 + x138*x63 + x141*x15 + x141*x31 + x141*x34 + x141*x42 - x141*x45 - x141*x48 - x145*x18 + x148*x24 + x155*x96 - x24*x92 - x30*x64 + x33*x64 - x36*x64 - x39*x64 - x47*x64 + x50*x64 + x64*x9;
    res[4]=(*this)[14]*x146 + (*this)[14]*x150 + (*this)[15]*x160 - (*this)[3]*x152 - (*this)[4]*x132 - (*this)[4]*x173 + (*this)[4]*x174 + (*this)[4]*x175 + (*this)[5]*x156 + (*this)[5]*x171 + (*this)[6]*A[3]*x131 - (*this)[6]*x158 - (*this)[6]*x176 + A[10]*x125 - A[1]*x161 + A[1]*x169 - A[2]*x163 + A[2]*x164 + A[3]*x162 + A[3]*x165 + A[3]*x166 + A[3]*x167 + A[3]*x168 + A[4]*x0 + A[4]*x1 + A[4]*x2 - A[4]*x3 - A[4]*x4 + A[4]*x5 - A[4]*x55 - A[4]*x56 + A[4]*x57 + A[4]*x58 + A[4]*x59 - A[4]*x6 - A[4]*x60 - A[4]*x61 + A[4]*x62 - A[4]*x7 + A[5]*x19 + A[5]*x27 + A[6]*x20 - A[6]*x26 + A[7]*x124 - A[7]*x127 + A[7]*x142 + A[7]*x147 - A[7]*x153 - A[8]*x68 + A[8]*x75 - A[8]*x83 + A[8]*x89 + A[8]*x91 + A[9]*x101 + A[9]*x143 + A[9]*x144 + A[9]*x149 - A[9]*x151 + A[9]*x154 + x104*x145 + x107*x135 + x11*x123 + x113*x172 + x114*x137 + x123*x38 + x123*x40 + x123*x43 - x123*x49 - x123*x52 - x13*x141 - x130*x77 + x141*x37 - x141*x41 - x141*x44 - x141*x46 - x141*x51 - x155*x90 - x157*x18 + x159*x87 + x170*x96 - x24*x99 - x30*x65 - x33*x65 + x36*x65 - x39*x65 - x47*x65 - x50*x65 + x65*x9;
    res[5]=(*this)[13]*x183 + (*this)[14]*x184 - (*this)[15]*x191 + (*this)[2]*x129 + (*this)[3]*x133 - (*this)[4]*x171 - (*this)[5]*x132 - (*this)[5]*x173 + (*this)[5]*x175 + (*this)[6]*x189 + (*this)[6]*x193 + (*this)[7]*x194 + (*this)[8]*x146 + (*this)[8]*x150 + A[10]*x180 + A[10]*x181 + A[10]*x186 + A[10]*x71 + A[1]*x101 - A[1]*x116 - A[1]*x143 + A[1]*x144 - A[1]*x151 - A[1]*x154 - A[2]*x124 - A[2]*x127 + A[2]*x142 - A[2]*x147 - A[3]*x21 - A[3]*x28 + A[4]*x19 + A[4]*x27 + A[5]*x0 - A[5]*x1 - A[5]*x2 - A[5]*x3 + A[5]*x4 + A[5]*x5 + A[5]*x55 - A[5]*x56 - A[5]*x57 - A[5]*x58 + A[5]*x59 + A[5]*x6 - A[5]*x60 + A[5]*x61 + A[5]*x62 - A[5]*x7 - A[6]*x121 + A[6]*x140 - A[6]*x188 + A[7]*x163 + A[7]*x164 + A[7]*x187 - A[8]*x23 + A[8]*x25 - A[9]*x161 + A[9]*x169 + A[9]*x179 - A[9]*x182 + A[9]*x185 - x100*x137 + x107*x172 - x113*x135 - x114*x24 + x115*x24 - x123*x190 - x123*x192 + x13*x177 - x131*x148 - x137*x99 - x15*x178 - x17*x63 + x177*x37 - x177*x41 + x177*x44 + x177*x46 + x177*x51 + x178*x31 + x178*x34 + x178*x42 + x178*x45 + x178*x48 + x29*x63 - x32*x63 - x35*x63 + x53*x63 - x54*x63;
    res[6]=(*this)[13]*x85 + (*this)[15]*x194 - (*this)[2]*x133 + (*this)[3]*x129 + (*this)[4]*x158 + (*this)[4]*x176 - (*this)[4]*x82 - (*this)[5]*x128 - (*this)[5]*x189 + (*this)[5]*x196 - (*this)[6]*x132 - (*this)[6]*x173 + (*this)[6]*x174 + (*this)[7]*x191 - A[10]*x70 + A[10]*x78 - A[10]*x79 + A[10]*x93 - A[1]*x68 - A[1]*x75 + A[1]*x83 + A[1]*x89 - A[1]*x91 + A[2]*x21 + A[2]*x28 - A[3]*x124 + A[3]*x127 + A[3]*x142 - A[3]*x153 + A[4]*x20 - A[4]*x26 + A[5]*x121 + A[5]*x140 + A[5]*x188 + A[6]*x0 - A[6]*x1 - A[6]*x2 + A[6]*x3 - A[6]*x4 + A[6]*x5 - A[6]*x55 + A[6]*x56 - A[6]*x57 - A[6]*x58 + A[6]*x59 - A[6]*x6 + A[6]*x60 - A[6]*x61 + A[6]*x62 + A[6]*x7 - A[7]*x162 + A[7]*x165 - A[7]*x167 + A[7]*x168 + A[8]*x161 + A[8]*x169 - A[8]*x179 + A[8]*x182 + A[8]*x185 - A[9]*x23 + A[9]*x25 - x105*x24 + x108*x137 - x11*x177 - x110*x137 + x113*x24 - x114*x135 + x130*x95 - x131*x87 - x134*x76 - x134*x77 + x15*x195 - x17*x64 + x177*x38 + x177*x40 - x177*x43 + x177*x49 + x177*x52 + x195*x31 + x195*x34 - x195*x42 + x195*x45 + x195*x48 - x29*x64 + x32*x64 - x35*x64 - x53*x64 + x54*x64;
    res[7]=(*this)[14]*x74 + (*this)[14]*x85 - (*this)[15]*x193 + (*this)[15]*x196 - (*this)[15]*x81 + (*this)[1]*x72 + (*this)[2]*x171 - (*this)[3]*x158 - (*this)[3]*x176 + (*this)[7]*x174 + (*this)[7]*x175 - (*this)[8]*x86 - (*this)[9]*x126 - (*this)[9]*x67 + A[10]*x23 + A[10]*x25 - A[1]*x112 + A[1]*x201 - A[1]*x69 + A[1]*x73 - A[1]*x88 + A[1]*x94 - A[1]*x97 + A[2]*x19 - A[2]*x27 + A[3]*x20 + A[3]*x26 + A[5]*x164 - A[5]*x187 + A[6]*x162 + A[6]*x165 - A[6]*x166 + A[6]*x167 + A[6]*x168 + A[7]*x0 + A[7]*x1 - A[7]*x2 + A[7]*x3 + A[7]*x4 - A[7]*x5 - A[7]*x55 - A[7]*x56 + A[7]*x57 - A[7]*x58 + A[7]*x59 - A[7]*x6 + A[7]*x60 + A[7]*x61 - A[7]*x62 - A[7]*x7 - A[8]*x181 + A[8]*x186 + A[8]*x197 - A[8]*x71 + A[9]*x198 + A[9]*x70 + A[9]*x78 + A[9]*x79 + A[9]*x93 + x100*x130 + x109*x131 + x11*x178 - x111*x131 - x115*x159 - x118*x159 - x13*x195 - x134*x200 - x137*x199 - x17*x65 + x170*x18 - x172*x92 - x172*x95 + x178*x38 - x178*x40 + x178*x43 + x178*x49 + x178*x52 + x195*x37 + x195*x41 - x195*x44 + x195*x46 + x195*x51 + x24*x76 - x29*x65 - x32*x65 + x35*x65 - x53*x65 - x54*x65;
    res[8]=(*this)[15]*A[7]*x104 + (*this)[2]*x139 + (*this)[3]*x183 + (*this)[4]*x184 - (*this)[5]*x209 + (*this)[6]*x207 - (*this)[7]*x86 + (*this)[8]*x175 + (*this)[9]*x81 + A[10]*x163 + A[10]*x164 - A[10]*x187 + A[10]*x203 - A[10]*x206 - A[1]*x20 - A[1]*x26 - A[2]*x125 + A[3]*x112 + A[3]*x201 - A[3]*x69 - A[3]*x73 - A[3]*x88 - A[3]*x97 + A[4]*x68 + A[4]*x91 + A[5]*x23 + A[5]*x25 + A[6]*x161 + A[6]*x169 - A[6]*x179 - A[6]*x182 - A[7]*x119 - A[7]*x180 + A[7]*x181 - A[7]*x71 + A[8]*x0 - A[8]*x1 + A[8]*x2 + A[8]*x3 - A[8]*x4 - A[8]*x5 - A[8]*x55 + A[8]*x56 + A[8]*x57 - A[8]*x58 - A[8]*x59 + A[8]*x6 - A[8]*x60 + A[8]*x61 + A[8]*x62 - A[8]*x7 - A[9]*x121 + A[9]*x188 - A[9]*x204 - A[9]*x205 - x100*x104 + x102*x108 + x102*x109 - x102*x111 + x106*x114 - x106*x115 - x106*x118 + x11*x202 + x113*x96 - x117*x96 + x141*x17 + x141*x29 - x141*x32 - x141*x35 - x141*x53 + x141*x54 + x18*x92 - x190*x64 + x192*x64 - x195*x30 + x195*x33 + x195*x36 + x195*x39 + x195*x47 - x195*x50 + x195*x9 + x202*x38 - x202*x40 - x202*x43 - x202*x49 + x202*x52 - x208*x64 - x77*x90;
    res[9]=(*this)[14]*x82 - (*this)[15]*x86 + (*this)[3]*x139 + (*this)[3]*x74 - (*this)[5]*x207 - (*this)[6]*A[5]*x102 - (*this)[6]*x209 - (*this)[7]*x126 - (*this)[7]*x67 - (*this)[8]*x81 + (*this)[9]*x174 - A[10]*x162 + A[10]*x165 - A[10]*x166 + A[10]*x167 - A[10]*x168 + A[1]*x19 - A[1]*x27 - A[2]*x112 + A[2]*x69 + A[2]*x73 - A[2]*x88 - A[2]*x94 + A[2]*x97 + A[5]*x169 + A[5]*x179 + A[6]*x23 + A[6]*x25 - A[7]*x103 - A[7]*x120 + A[7]*x70 - A[7]*x78 - A[7]*x79 + A[7]*x93 + A[8]*x121 - A[8]*x188 + A[8]*x204 - A[8]*x205 + A[9]*x0 - A[9]*x1 + A[9]*x2 - A[9]*x3 + A[9]*x4 - A[9]*x5 + A[9]*x55 - A[9]*x56 + A[9]*x57 - A[9]*x58 - A[9]*x59 - A[9]*x6 + A[9]*x60 - A[9]*x61 + A[9]*x62 + A[9]*x7 + x102*x99 + x104*x107 + x104*x109 - x104*x110 - x106*x113 + x106*x117 + x106*x200 + x114*x96 - x118*x96 + x123*x17 - x123*x29 + x123*x32 - x123*x35 + x123*x53 - x123*x54 + x13*x202 + x178*x30 - x178*x33 + x178*x36 + x178*x39 - x178*x47 + x178*x50 + x178*x9 - x18*x87 + x190*x63 + x192*x63 - x202*x37 - x202*x41 - x202*x44 + x202*x46 - x202*x51 - x208*x63 - x66*x92 + x76*x84;
    res[10]=-(*this)[13]*x82 + (*this)[2]*A[6]*x106 + (*this)[4]*x85 - (*this)[5]*A[2]*x96 + (*this)[6]*A[2]*x106 + (*this)[6]*x67 - (*this)[6]*x98 + A[10]*x0 + A[10]*x1 + A[10]*x2 - A[10]*x3 - A[10]*x4 + A[10]*x5 + A[10]*x55 + A[10]*x56 - A[10]*x57 - A[10]*x58 - A[10]*x59 - A[10]*x6 + A[10]*x60 + A[10]*x61 - A[10]*x62 - A[10]*x7 - A[1]*x21 - A[1]*x28 - A[2]*x68 - A[2]*x83 - A[2]*x89 - A[2]*x91 + A[3]*x143 + A[3]*x149 - A[3]*x154 + A[5]*x197 + A[6]*x103 - A[6]*x120 + A[6]*x198 - A[6]*x78 + A[6]*x79 - A[7]*x23 + A[7]*x25 - A[8]*x163 + A[8]*x164 + A[8]*x187 - A[8]*x203 - A[8]*x206 + A[9]*x162 + A[9]*x165 + A[9]*x166 - A[9]*x167 - A[9]*x168 - x100*x66 + x102*x148 - x102*x95 - x104*x199 + x104*x87 + x105*x84 - x106*x76 - x106*x77 + x107*x18 + x108*x18 + x109*x18 - x110*x18 + x115*x90 + x122*x17 - x122*x29 - x122*x32 + x122*x35 + x122*x53 + x122*x54 + x15*x202 + x177*x30 + x177*x33 - x177*x36 + x177*x39 - x177*x47 - x177*x50 + x177*x9 - x202*x31 + x202*x34 - x202*x42 - x202*x45 + x202*x48 + x210*x64 - x211*x63 + x212*x63 + x213*x64 - x214*x63 - x215*x64;
    res[11]=-(*this)[0]*x217 + (*this)[15]*x230 + (*this)[2]*x224 - (*this)[3]*x221 + (*this)[3]*x232 + (*this)[3]*x235 + (*this)[4]*x222 + (*this)[4]*x233 + (*this)[5]*x234 + (*this)[7]*x220 + (*this)[7]*x236 - (*this)[8]*x223 - (*this)[9]*x226 + A[0]*x125 + A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3 + A[11]*x4 + A[11]*x5 + A[11]*x55 + A[11]*x56 + A[11]*x57 + A[11]*x58 + A[11]*x59 + A[11]*x6 + A[11]*x60 + A[11]*x61 + A[11]*x62 + A[11]*x7 + A[12]*x78 - A[12]*x79 - A[13]*x180 - A[13]*x181 + A[14]*x179 - A[14]*x182 - A[15]*x124 - A[15]*x147 + x10*x218 - x10*x225 + x10*x228 + x10*x229 + x10*x231 - x10*x240 + x12*x219 - x12*x227 + x12*x239 + x14*x237 - x14*x238 + x216*x30 + x216*x33 + x216*x36 - x216*x39 + x216*x47 + x216*x50 + x216*x9;
    res[12]=(*this)[0]*x234 + (*this)[12]*x241 - (*this)[13]*x222 - (*this)[14]*x221 + (*this)[14]*x232 + (*this)[14]*x235 - (*this)[15]*x223 - (*this)[1]*x224 - (*this)[5]*x217 - (*this)[6]*x22 - (*this)[7]*x226 + (*this)[8]*x230 + (*this)[9]*x220 + (*this)[9]*x236 - A[0]*x20 + A[0]*x26 + A[12]*x0 - A[12]*x1 - A[12]*x2 - A[12]*x3 + A[12]*x4 + A[12]*x5 - A[12]*x55 + A[12]*x56 + A[12]*x57 + A[12]*x58 - A[12]*x59 + A[12]*x6 + A[12]*x60 - A[12]*x61 - A[12]*x62 - A[12]*x7 - A[13]*x121 + A[13]*x204 + A[14]*x164 - A[14]*x203 - A[15]*x144 - A[15]*x149 + x11*x216 + x12*x190 - x12*x192 + x12*x208 - x14*x211 - x14*x212 + x14*x214 + x216*x38 + x216*x40 - x216*x43 + x216*x49 - x216*x52 + x218*x8 - x225*x8 + x228*x8 + x229*x8 - x231*x8 + x240*x8;
    res[13]=(*this)[12]*x233 - (*this)[12]*x243 + (*this)[13]*x241 - (*this)[14]*x224 - (*this)[15]*x226 + (*this)[15]*x247 - (*this)[1]*x235 - (*this)[1]*x245 - (*this)[3]*x248 + (*this)[4]*x244 + (*this)[5]*x22 - (*this)[5]*x246 + (*this)[7]*x223 - (*this)[8]*x220 - (*this)[8]*x236 + (*this)[8]*x249 + A[0]*x19 + A[0]*x27 + A[12]*x121 - A[12]*x204 + A[13]*x0 - A[13]*x1 - A[13]*x2 + A[13]*x3 - A[13]*x4 + A[13]*x5 + A[13]*x55 - A[13]*x56 + A[13]*x57 + A[13]*x58 - A[13]*x59 - A[13]*x6 - A[13]*x60 + A[13]*x61 - A[13]*x62 + A[13]*x7 + A[14]*x165 + A[14]*x166 + A[15]*x75 - A[15]*x83 - x10*x190 - x10*x192 + x10*x208 + x13*x216 + x14*x210 - x14*x213 + x14*x215 + x16*x242 - x216*x37 + x216*x41 - x216*x44 - x216*x46 + x216*x51 + x219*x8 - x227*x8 - x239*x8;
    res[14]=(*this)[0]*x236 - (*this)[0]*x249 + (*this)[0]*x250 + (*this)[0]*x251 - (*this)[12]*x235 - (*this)[12]*x245 - (*this)[13]*x252 + (*this)[13]*x253 + (*this)[14]*x241 + (*this)[15]*x22 - (*this)[15]*x246 - (*this)[1]*x222 + (*this)[1]*x243 - (*this)[3]*x244 - (*this)[4]*x248 + (*this)[5]*x226 - (*this)[7]*x217 - A[0]*x21 + A[0]*x28 + A[12]*x164 + A[12]*x203 + A[13]*x165 - A[13]*x166 + A[14]*x0 + A[14]*x1 - A[14]*x2 + A[14]*x3 + A[14]*x4 - A[14]*x5 + A[14]*x55 + A[14]*x56 - A[14]*x57 + A[14]*x58 - A[14]*x59 - A[14]*x6 - A[14]*x60 - A[14]*x61 + A[14]*x62 - A[14]*x7 - A[15]*x73 + x10*x211 - x10*x212 + x10*x214 - x12*x210 - x12*x213 + x12*x215 - x136*x16 + x138*x16 + x15*x216 - x216*x31 + x216*x34 + x216*x42 + x216*x45 - x216*x48 - x237*x8 + x238*x8;
    res[15]=-(*this)[0]*(*this)[2]*x10 - (*this)[0]*x248 - (*this)[12]*(*this)[15]*x10 - (*this)[13]*(*this)[15]*x12 - (*this)[13]*x247 - (*this)[14]*(*this)[15]*x14 - (*this)[14]*x22 + (*this)[15]*x241 - (*this)[1]*(*this)[5]*x10 - (*this)[3]*(*this)[6]*x8 - (*this)[4]*x220 + (*this)[4]*x249 - (*this)[4]*x250 - (*this)[4]*x251 - (*this)[5]*x224 + (*this)[5]*x252 - (*this)[5]*x253 + (*this)[6]*x221 - (*this)[6]*x232 - (*this)[6]*x245 - (*this)[7]*x222 - (*this)[7]*x233 + (*this)[7]*x243 - (*this)[8]*x244 - A[0]*x23 - A[0]*x25 - A[12]*x144 + A[12]*x149 + A[13]*x75 + A[13]*x83 + A[14]*x73 + A[15]*x0 + A[15]*x1 - A[15]*x2 - A[15]*x3 - A[15]*x4 - A[15]*x5 + A[15]*x55 + A[15]*x56 + A[15]*x57 - A[15]*x58 + A[15]*x59 + A[15]*x6 - A[15]*x60 - A[15]*x61 - A[15]*x62 + A[15]*x7 - x12*x242 - x136*x14 + x138*x14 + x17*x216 - x216*x29 - x216*x32 - x216*x35 - x216*x53 - x216*x54;
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
    T x16 = 2.0*A[1];
    T x17 = (*this)[0]*x16;
    T x18 = 2.0*A[2];
    T x19 = (*this)[0]*x18;
    T x20 = 2.0*A[3];
    T x21 = (*this)[0]*x20;
    T x22 = 2.0*(*this)[10];
    T x23 = (*this)[12]*A[2];
    T x24 = A[1]*x22;
    T x25 = A[3]*x22;
    T x26 = (*this)[5]*x16;
    T x27 = (*this)[6]*x18;
    T x28 = (*this)[7]*x20;
    T x29 = 2.0*(*this)[9];
    T x30 = A[3]*x29;
    T x31 = 2.0*(*this)[8];
    T x32 = A[3]*x31;
    T x33 = A[2]*x31;
    T x34 = (*this)[14]*A[1];
    T x35 = (*this)[15]*x16;
    T x36 = (*this)[15]*x18;
    T x37 = (*this)[15]*x20;
    T x38 = A[1]*x31;
    T x39 = A[2]*x29;
    T x40 = (*this)[6]*x20;
    T x41 = (*this)[7]*x18;
    T x42 = (*this)[5]*x20;
    T x43 = (*this)[7]*x16;
    T x44 = (*this)[5]*x18;
    T x45 = (*this)[6]*x16;
    T x46 = A[2]*x22;
    T x47 = 2.0*x23;
    T x48 = 2.0*x34;
    T x49 = A[1]*x29;
    T x50 = (*this)[13]*x20;
    T x51 = (*this)[14]*x18;
    T x52 = (*this)[11]*(*this)[1];
    T x53 = (*this)[12]*x16;
    T x54 = (*this)[12]*x20;
    T x55 = (*this)[13]*(*this)[2];
    T x56 = (*this)[13]*(*this)[3];
    T x57 = (*this)[14]*x20;
    T x58 = (*this)[1]*x20;
    T x59 = (*this)[1]*(*this)[4];
    T x60 = (*this)[11]*x16;
    T x61 = (*this)[4]*x18;
    T x62 = (*this)[2]*x18;
    T x63 = (*this)[3]*x16;
    T x64 = (*this)[3]*x20;
    T x65 = (*this)[2]*x20;
    T x66 = 2.0*A[0];
    T x67 = (*this)[0]*x66;
    T x68 = A[0]*x22;
    T x69 = (*this)[5]*x66;
    T x70 = (*this)[6]*x66;
    T x71 = (*this)[7]*x66;
    T x72 = (*this)[15]*x66;
    T x73 = A[0]*x31;
    T x74 = A[0]*x29;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 - A[0]*x13 - A[0]*x14 - A[0]*x15 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[0]*x8 - A[0]*x9;
    res[1]=-(*this)[11]*x26 - (*this)[11]*x27 - (*this)[11]*x28 - (*this)[12]*x17 - (*this)[12]*x30 - (*this)[13]*x19 - (*this)[13]*x24 + (*this)[13]*x32 - (*this)[14]*x21 - (*this)[14]*x33 + (*this)[1]*x25 + (*this)[1]*x38 + (*this)[1]*x39 - (*this)[2]*x35 + (*this)[2]*x40 - (*this)[2]*x41 - (*this)[3]*x36 - (*this)[3]*x42 + (*this)[3]*x43 - (*this)[4]*x37 + (*this)[4]*x44 - (*this)[4]*x45 + x22*x23 + x29*x34;
    res[2]=(*this)[11]*x17 + (*this)[11]*x30 - (*this)[11]*x46 + (*this)[12]*x26 + (*this)[12]*x28 - (*this)[13]*x37 + (*this)[13]*x44 - (*this)[13]*x45 + (*this)[14]*x36 + (*this)[14]*x42 - (*this)[1]*x35 + (*this)[1]*x40 - (*this)[1]*x41 + (*this)[2]*x25 + (*this)[2]*x38 + (*this)[2]*x39 + (*this)[3]*x21 + (*this)[3]*x33 - (*this)[3]*x49 - (*this)[4]*x19 - (*this)[4]*x24 + (*this)[4]*x32 + (*this)[6]*x47 - (*this)[7]*x48;
    res[3]=(*this)[11]*x19 + (*this)[11]*x24 - (*this)[11]*x32 + (*this)[12]*x37 + (*this)[12]*x45 + (*this)[13]*x26 + (*this)[13]*x27 + (*this)[13]*x28 + (*this)[14]*x40 - (*this)[14]*x41 - (*this)[15]*x48 - (*this)[1]*x36 - (*this)[1]*x42 + (*this)[1]*x43 - (*this)[2]*x21 - (*this)[2]*x33 + (*this)[2]*x49 + (*this)[3]*x25 + (*this)[3]*x38 + (*this)[3]*x39 + (*this)[4]*x17 + (*this)[4]*x30 - (*this)[4]*x46 - (*this)[5]*x47;
    res[4]=(*this)[11]*x21 + (*this)[11]*x33 - (*this)[11]*x49 - (*this)[12]*x42 + (*this)[12]*x43 + (*this)[13]*x35 - (*this)[13]*x40 + (*this)[13]*x41 + (*this)[14]*x27 + (*this)[14]*x28 - (*this)[15]*x47 - (*this)[1]*x37 + (*this)[1]*x44 - (*this)[1]*x45 + (*this)[2]*x19 + (*this)[2]*x24 - (*this)[2]*x32 - (*this)[3]*x17 - (*this)[3]*x30 + (*this)[3]*x46 + (*this)[4]*x25 + (*this)[4]*x38 + (*this)[4]*x39 + (*this)[5]*x48;
    res[5]=(*this)[0]*x40 + (*this)[11]*x50 - (*this)[11]*x51 - (*this)[15]*x17 - (*this)[15]*x30 + (*this)[15]*x46 + (*this)[2]*x53 + (*this)[2]*x57 + (*this)[3]*x47 + (*this)[3]*x58 - (*this)[4]*x48 + (*this)[4]*x54 + (*this)[5]*x25 + (*this)[5]*x38 + (*this)[5]*x39 + (*this)[6]*x33 - (*this)[6]*x49 - (*this)[7]*x19 - (*this)[7]*x24 + (*this)[7]*x32 + x16*x52 - x16*x56 + x18*x55 - x18*x59;
    res[6]=(*this)[0]*x43 + (*this)[11]*x48 - (*this)[11]*x54 - (*this)[15]*x19 - (*this)[15]*x24 + (*this)[15]*x32 - (*this)[2]*x47 - (*this)[2]*x58 + (*this)[3]*x53 + (*this)[3]*x57 + (*this)[4]*x50 - (*this)[4]*x51 - (*this)[5]*x21 - (*this)[5]*x33 + (*this)[5]*x49 + (*this)[6]*x25 + (*this)[6]*x38 + (*this)[6]*x39 + (*this)[7]*x30 - (*this)[7]*x46 + x16*x55 + x16*x59 + x18*x52 + x18*x56;
    res[7]=(*this)[0]*x44 + (*this)[11]*x47 - (*this)[13]*x60 + (*this)[13]*x61 - (*this)[15]*x21 - (*this)[15]*x33 + (*this)[15]*x49 + (*this)[1]*x62 - (*this)[1]*x63 + (*this)[2]*x48 - (*this)[2]*x54 - (*this)[3]*x50 + (*this)[3]*x51 + (*this)[4]*x53 + (*this)[4]*x57 + (*this)[5]*x24 - (*this)[5]*x32 - (*this)[6]*x17 - (*this)[6]*x30 + (*this)[6]*x46 + (*this)[7]*x25 + (*this)[7]*x38 + (*this)[7]*x39 + x20*x52;
    res[8]=(*this)[0]*x30 - (*this)[0]*x46 - (*this)[11]*x61 + (*this)[11]*x64 - (*this)[13]*x47 - (*this)[14]*x54 + (*this)[15]*x40 - (*this)[1]*x50 + (*this)[1]*x51 + (*this)[3]*x62 + (*this)[4]*x65 - (*this)[5]*x28 - (*this)[6]*x44 - (*this)[7]*x36 + (*this)[8]*x25 + (*this)[9]*x33 + A[1]*x0 - A[1]*x1 + A[1]*x10 - A[1]*x11 - A[1]*x12 - A[1]*x13 + A[1]*x14 + A[1]*x15 + A[1]*x2 + A[1]*x3 - A[1]*x4 - A[1]*x5 + A[1]*x6 - A[1]*x7 - A[1]*x8 + A[1]*x9;
    res[9]=(*this)[0]*x24 - (*this)[0]*x32 - (*this)[11]*x65 - (*this)[13]*x53 - (*this)[14]*x50 + (*this)[15]*x43 - (*this)[1]*x48 + (*this)[1]*x54 + (*this)[2]*x63 + (*this)[4]*x60 + (*this)[4]*x64 - (*this)[5]*x37 - (*this)[6]*x26 - (*this)[7]*x40 + (*this)[9]*x25 + (*this)[9]*x38 + A[2]*x0 - A[2]*x1 + A[2]*x10 - A[2]*x11 - A[2]*x12 + A[2]*x13 - A[2]*x14 + A[2]*x15 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 - A[2]*x6 + A[2]*x7 + A[2]*x8 - A[2]*x9;
    res[10]=(*this)[0]*x33 - (*this)[0]*x49 + (*this)[11]*x62 - (*this)[12]*x48 + (*this)[13]*(*this)[1]*x16 - (*this)[13]*x51 + (*this)[15]*x44 - (*this)[1]*x47 + (*this)[2]*(*this)[4]*x16 - (*this)[3]*x60 + (*this)[3]*x61 - (*this)[5]*x43 - (*this)[6]*x35 - (*this)[7]*x27 + (*this)[8]*x24 + (*this)[9]*x46 + A[3]*x0 + A[3]*x1 - A[3]*x10 - A[3]*x11 - A[3]*x12 + A[3]*x13 + A[3]*x14 - A[3]*x15 + A[3]*x2 - A[3]*x3 - A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7 + A[3]*x8 + A[3]*x9;
    res[11]=(*this)[11]*x67 + (*this)[12]*x69 + (*this)[13]*x70 + (*this)[14]*x71 - (*this)[1]*x72 + (*this)[2]*x73 + (*this)[3]*x74 + (*this)[4]*x68;
    res[12]=(*this)[11]*x69 + (*this)[12]*x67 - (*this)[13]*x68 + (*this)[14]*x74 - (*this)[1]*x73 + (*this)[2]*x72 + (*this)[3]*x71 - (*this)[4]*x70;
    res[13]=(*this)[11]*x70 + (*this)[12]*x68 + (*this)[13]*x67 - (*this)[14]*x73 - (*this)[1]*x74 - (*this)[2]*x71 + (*this)[3]*x72 + (*this)[4]*x69;
    res[14]=(*this)[11]*x71 - (*this)[12]*x74 + (*this)[13]*x73 + (*this)[14]*x67 - (*this)[1]*x68 + (*this)[2]*x70 - (*this)[3]*x69 + (*this)[4]*x72;
    res[15]=-(*this)[12]*(*this)[2]*x66 - (*this)[14]*(*this)[4]*x66 + (*this)[15]*x67 - (*this)[5]*x73 - (*this)[6]*x74 - (*this)[7]*x68 - x52*x66 - x56*x66;
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
    T x8 = 2.0*A[7];
    T x9 = (*this)[0]*(*this)[15];
    T x10 = (*this)[10]*(*this)[7];
    T x11 = (*this)[11]*(*this)[1];
    T x12 = (*this)[12]*(*this)[2];
    T x13 = (*this)[13]*(*this)[3];
    T x14 = (*this)[14]*(*this)[4];
    T x15 = (*this)[5]*(*this)[8];
    T x16 = (*this)[6]*(*this)[9];
    T x17 = std::pow((*this)[12], 2);
    T x18 = std::pow((*this)[13], 2);
    T x19 = std::pow((*this)[14], 2);
    T x20 = std::pow((*this)[15], 2);
    T x21 = std::pow((*this)[1], 2);
    T x22 = std::pow((*this)[5], 2);
    T x23 = std::pow((*this)[6], 2);
    T x24 = std::pow((*this)[7], 2);
    T x25 = 2.0*(*this)[0];
    T x26 = A[4]*x25;
    T x27 = A[5]*x25;
    T x28 = (*this)[14]*x25;
    T x29 = A[1]*x25;
    T x30 = (*this)[3]*x25;
    T x31 = A[3]*x25;
    T x32 = 2.0*(*this)[10];
    T x33 = (*this)[11]*x32;
    T x34 = (*this)[12]*x32;
    T x35 = A[4]*x32;
    T x36 = A[6]*x32;
    T x37 = A[2]*x32;
    T x38 = (*this)[3]*A[1];
    T x39 = 2.0*A[4];
    T x40 = (*this)[5]*x39;
    T x41 = 2.0*(*this)[6];
    T x42 = (*this)[11]*x41;
    T x43 = 2.0*(*this)[11];
    T x44 = (*this)[7]*A[6];
    T x45 = A[1]*x43;
    T x46 = 2.0*(*this)[9];
    T x47 = (*this)[11]*x46;
    T x48 = 2.0*(*this)[12];
    T x49 = (*this)[15]*A[1];
    T x50 = (*this)[12]*x41;
    T x51 = (*this)[7]*A[2];
    T x52 = (*this)[12]*x46;
    T x53 = 2.0*(*this)[13];
    T x54 = (*this)[15]*x53;
    T x55 = (*this)[5]*A[3];
    T x56 = (*this)[7]*A[1];
    T x57 = (*this)[8]*x53;
    T x58 = 2.0*(*this)[14];
    T x59 = (*this)[15]*x58;
    T x60 = (*this)[5]*A[2];
    T x61 = (*this)[14]*A[1];
    T x62 = (*this)[8]*x58;
    T x63 = A[4]*x58;
    T x64 = (*this)[15]*(*this)[2];
    T x65 = 2.0*A[5];
    T x66 = (*this)[15]*(*this)[3];
    T x67 = 2.0*(*this)[4];
    T x68 = (*this)[15]*x67;
    T x69 = 2.0*(*this)[1];
    T x70 = (*this)[5]*A[1];
    T x71 = A[2]*x41;
    T x72 = (*this)[7]*x69;
    T x73 = (*this)[8]*A[4];
    T x74 = (*this)[9]*x69;
    T x75 = A[6]*x41;
    T x76 = (*this)[2]*x65;
    T x77 = (*this)[2]*A[3];
    T x78 = 2.0*(*this)[3];
    T x79 = (*this)[5]*A[6];
    T x80 = (*this)[3]*(*this)[7];
    T x81 = (*this)[8]*A[3];
    T x82 = (*this)[5]*A[5];
    T x83 = A[4]*x41;
    T x84 = (*this)[8]*x67;
    T x85 = (*this)[4]*x46;
    T x86 = (*this)[8]*x48;
    T x87 = (*this)[9]*x53;
    T x88 = (*this)[15]*x69;
    T x89 = 2.0*A[3];
    T x90 = 2.0*(*this)[2];
    T x91 = 2.0*x77;
    T x92 = (*this)[2]*x46;
    T x93 = (*this)[8]*x65;
    T x94 = (*this)[3]*x46;
    T x95 = A[2]*x25;
    T x96 = A[6]*x25;
    T x97 = A[3]*x32;
    T x98 = (*this)[1]*A[1];
    T x99 = A[5]*x32;
    T x100 = (*this)[15]*x43;
    T x101 = (*this)[8]*x43;
    T x102 = (*this)[15]*x48;
    T x103 = A[4]*x53;
    T x104 = A[5]*x41;
    T x105 = (*this)[7]*A[5];
    T x106 = (*this)[9]*x58;
    T x107 = (*this)[2]*A[1];
    T x108 = 2.0*x38;
    T x109 = A[3]*x41;
    T x110 = (*this)[12]*x39;
    T x111 = 2.0*A[2];
    T x112 = (*this)[8]*A[2];
    T x113 = (*this)[8]*A[6];
    T x114 = A[3]*x67;
    T x115 = A[6]*x53;
    T x116 = A[5]*x58;
    T x117 = (*this)[3]*x43;
    T x118 = (*this)[11]*x67;
    T x119 = (*this)[12]*x53;
    T x120 = (*this)[12]*A[3];
    T x121 = (*this)[3]*A[5];
    T x122 = (*this)[12]*A[6];
    T x123 = A[3]*x53;
    T x124 = A[5]*x53;
    T x125 = A[2]*x58;
    T x126 = A[6]*x58;
    T x127 = 2.0*(*this)[15];
    T x128 = (*this)[15]*x46;
    T x129 = A[6]*x69;
    T x130 = (*this)[4]*x69;
    T x131 = (*this)[2]*A[2];
    T x132 = 2.0*(*this)[7];
    T x133 = 2.0*x44;
    T x134 = (*this)[8]*A[1];
    T x135 = 2.0*A[6];
    T x136 = (*this)[12]*x43;
    T x137 = A[1]*x58;
    T x138 = (*this)[12]*x69;
    T x139 = A[2]*x53;
    T x140 = (*this)[2]*A[5];
    T x141 = (*this)[3]*x69;
    T x142 = (*this)[3]*A[2];
    T x143 = 2.0*x56;
    T x144 = 2.0*A[1];
    T x145 = (*this)[15]*(*this)[7];
    T x146 = A[6]*x67;
    T x147 = (*this)[4]*x39;
    T x148 = (*this)[3]*x39;
    T x149 = A[0]*x25;
    T x150 = (*this)[1]*x8;
    T x151 = (*this)[10]*x8;
    T x152 = A[0]*x32;
    T x153 = (*this)[11]*x8;
    T x154 = (*this)[5]*A[0];
    T x155 = (*this)[12]*x8;
    T x156 = A[0]*x41;
    T x157 = (*this)[13]*x8;
    T x158 = (*this)[7]*A[0];
    T x159 = (*this)[2]*x8;
    T x160 = (*this)[8]*A[0];
    T x161 = (*this)[6]*x8;
    T x162 = (*this)[4]*x8;
    T x163 = 2.0*A[0];
    T x164 = (*this)[3]*x8;
    T x165 = (*this)[14]*x8;
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x17 - A[0]*x18 - A[0]*x19 + A[0]*x2 - A[0]*x20 - A[0]*x21 - A[0]*x22 - A[0]*x23 - A[0]*x24 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 + x10*x8 - x11*x8 - x12*x8 - x13*x8 - x14*x8 + x15*x8 + x16*x8 - x8*x9;
    res[1]=-(*this)[11]*x40 - (*this)[12]*x26 - (*this)[13]*x27 - (*this)[13]*x35 + (*this)[1]*x36 - (*this)[1]*x71 + (*this)[2]*x29 - (*this)[2]*x37 + (*this)[2]*x75 + (*this)[4]*x31 - (*this)[4]*x83 - (*this)[7]*x76 - (*this)[8]*x45 + (*this)[9]*x63 - A[1]*x85 + A[2]*x30 - A[2]*x47 - A[2]*x54 + A[2]*x84 - A[3]*x33 + A[3]*x50 - A[3]*x59 - A[3]*x72 + A[5]*x34 - A[5]*x42 - A[5]*x62 + A[5]*x74 - A[6]*x28 - A[6]*x52 + A[6]*x57 - A[6]*x68 + x32*x38 - x39*x64 + x39*x80 - x41*x61 - x43*x44 + x46*x77 - x48*x49 - x48*x51 - x53*x55 + x53*x56 + x58*x60 - x65*x66 + x67*x82 - x69*x70 + x69*x73 - x78*x79 - x78*x81;
    res[2]=(*this)[11]*x26 + (*this)[12]*x40 + (*this)[13]*x31 - (*this)[13]*x83 + (*this)[15]*x45 + (*this)[1]*x29 - (*this)[1]*x37 + (*this)[1]*x75 + (*this)[2]*x36 - (*this)[2]*x71 + (*this)[3]*x93 - (*this)[4]*x27 - (*this)[4]*x35 - (*this)[7]*x63 - (*this)[7]*x91 + A[1]*x86 - A[1]*x87 - A[2]*x28 + A[2]*x52 + A[2]*x57 - A[2]*x68 + A[3]*x34 - A[3]*x42 + A[3]*x62 + A[3]*x74 - A[4]*x88 - A[4]*x94 - A[5]*x33 + A[5]*x50 + A[5]*x59 - A[5]*x72 + A[5]*x92 + A[6]*x30 + A[6]*x47 - A[6]*x54 + A[6]*x84 - x32*x61 + x38*x41 + x43*x51 + x44*x48 + x53*x82 - x55*x67 + x56*x67 + x58*x79 - x60*x78 + x66*x89 - x70*x90 + x73*x90;
    res[3]=(*this)[11]*x27 - (*this)[12]*x31 + (*this)[13]*x104 + (*this)[13]*x97 + (*this)[14]*x29 - (*this)[14]*x37 + (*this)[14]*x75 - (*this)[15]*x63 - (*this)[15]*x91 + (*this)[1]*x95 - (*this)[2]*x96 + (*this)[3]*x36 - (*this)[3]*x71 - (*this)[4]*x109 + (*this)[4]*x26 - (*this)[4]*x99 + (*this)[5]*x103 - (*this)[5]*x108 - (*this)[8]*x76 + A[1]*x52 + A[1]*x57 + A[2]*x100 - A[2]*x86 + A[2]*x87 + A[3]*x106 + A[4]*x33 + A[4]*x50 + A[4]*x72 + A[4]*x92 - A[5]*x88 + A[5]*x94 - A[6]*x101 + A[6]*x102 + A[6]*x85 - x105*x58 - x107*x41 + x32*x98 + x43*x55 - x43*x56 + x44*x53 - x48*x82 + x49*x67 + x51*x67 + x60*x90 - x69*x79 - x69*x81 + x73*x78 - x80*x89;
    res[4]=(*this)[11]*x96 + (*this)[12]*x95 - (*this)[13]*x29 + (*this)[13]*x37 - (*this)[13]*x75 + (*this)[14]*x104 + (*this)[14]*x97 - (*this)[15]*x108 + (*this)[1]*x31 - (*this)[1]*x83 + (*this)[2]*x27 + (*this)[2]*x35 + (*this)[3]*x109 + (*this)[3]*x99 + (*this)[4]*x36 - (*this)[4]*x71 + (*this)[5]*x63 + (*this)[5]*x91 + (*this)[7]*x110 - (*this)[7]*x114 + A[1]*x34 + A[1]*x42 + A[1]*x62 - A[1]*x74 + A[2]*x106 + A[3]*x100 - A[3]*x87 - A[4]*x30 - A[4]*x47 + A[4]*x54 + A[5]*x101 - A[5]*x102 + A[5]*x85 - A[6]*x88 - A[6]*x94 + x105*x53 + x111*x64 - x111*x80 + x112*x69 - x113*x90 - x43*x60 + x44*x58 - x48*x79 - x48*x81 - x56*x90 - x67*x70 + x67*x73 + x69*x82;
    res[5]=-(*this)[10]*x95 + (*this)[11]*x115 - (*this)[11]*x116 + (*this)[15]*x109 + (*this)[15]*x99 + (*this)[1]*x123 - (*this)[1]*x125 + (*this)[2]*x124 + (*this)[2]*x126 + (*this)[3]*x129 + (*this)[5]*x36 + (*this)[6]*x96 - (*this)[7]*x27 + (*this)[8]*x104 + (*this)[8]*x133 + (*this)[9]*x31 + A[1]*x0 - A[1]*x1 + A[1]*x17 - A[1]*x18 - A[1]*x19 - A[1]*x2 - A[1]*x20 + A[1]*x21 - A[1]*x22 + A[1]*x23 + A[1]*x24 - A[1]*x3 + A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 + A[2]*x118 + A[2]*x119 - A[3]*x117 - A[5]*x130 - A[6]*x128 - x10*x39 + x11*x39 + x112*x46 + x12*x39 + x120*x58 + x121*x48 + x122*x67 - x127*x51 - x13*x39 - x131*x78 - x132*x55 - x14*x39 + x15*x39 - x16*x39 + x32*x81 - x39*x9 - x41*x60 + x46*x82 - x67*x77;
    res[6]=(*this)[10]*x29 + (*this)[11]*x63 + (*this)[14]*x123 - (*this)[15]*x35 + (*this)[2]*x103 - (*this)[2]*x129 + (*this)[3]*x110 - (*this)[3]*x114 + (*this)[3]*x126 + (*this)[4]*x115 + (*this)[5]*A[4]*x46 + (*this)[6]*x36 - (*this)[7]*x109 + (*this)[7]*x26 - (*this)[8]*x31 + (*this)[9]*x97 - A[1]*x118 + A[1]*x119 + A[2]*x0 - A[2]*x1 - A[2]*x17 + A[2]*x18 - A[2]*x19 - A[2]*x2 - A[2]*x20 + A[2]*x21 + A[2]*x22 - A[2]*x23 + A[2]*x24 + A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 + A[4]*x130 - x10*x65 + x11*x65 + x113*x127 - x12*x65 - x120*x69 - x122*x43 - x127*x55 + x127*x56 + x13*x65 + x134*x46 - x14*x65 - x15*x65 + x16*x65 - x25*x79 - x38*x90 - x41*x70 + x41*x73 + x43*x77 + x44*x46 + x58*x98 - x65*x9;
    res[7]=-(*this)[11]*x103 + (*this)[12]*x137 + (*this)[14]*x139 - (*this)[15]*x93 + (*this)[2]*x63 + (*this)[4]*x110 + (*this)[4]*x124 - (*this)[5]*x143 + (*this)[5]*x35 - (*this)[6]*x26 + (*this)[6]*x99 + (*this)[8]*x95 - (*this)[9]*x29 + (*this)[9]*x37 + A[2]*x138 + A[3]*x0 + A[3]*x1 - A[3]*x17 - A[3]*x18 + A[3]*x19 - A[3]*x2 - A[3]*x20 + A[3]*x21 + A[3]*x22 + A[3]*x23 - A[3]*x24 + A[3]*x3 + A[3]*x4 - A[3]*x5 - A[3]*x6 - A[3]*x7 + A[4]*x128 - A[4]*x141 + A[5]*x136 + x10*x135 + x105*x46 - x107*x67 + x11*x135 - x12*x135 + x121*x58 + x127*x60 - x13*x135 - x131*x43 + x132*x73 + x134*x32 + x135*x14 - x135*x15 - x135*x16 - x135*x9 + x140*x69 - x142*x67 + x25*x82 + x38*x43 - x41*x49 - x41*x51 - x53*x98;
    res[8]=-(*this)[10]*x27 + (*this)[11]*x123 - (*this)[11]*x125 + (*this)[12]*x114 - (*this)[15]*x37 + (*this)[15]*x75 - (*this)[1]*x115 + (*this)[1]*x116 + (*this)[2]*x146 + (*this)[3]*x76 - (*this)[5]*x133 - (*this)[6]*x31 + (*this)[8]*A[5]*x46 + (*this)[8]*x36 - (*this)[8]*x71 + (*this)[9]*x96 - A[2]*x130 + A[3]*x128 + A[3]*x141 + A[4]*x0 - A[4]*x1 - A[4]*x17 + A[4]*x18 + A[4]*x19 + A[4]*x2 - A[4]*x20 - A[4]*x21 - A[4]*x22 + A[4]*x23 + A[4]*x24 + A[4]*x3 - A[4]*x4 - A[4]*x5 + A[4]*x6 - A[4]*x7 - A[5]*x118 - A[5]*x119 + A[6]*x117 + x10*x144 + x11*x144 + x12*x144 - x122*x58 - x13*x144 + x131*x53 - x132*x81 - x14*x144 + x142*x48 - x144*x15 + x144*x16 + x144*x9 - x145*x65 + x25*x51 - x32*x55 - x41*x82 - x46*x60 + x58*x77;
    res[9]=(*this)[10]*x26 + (*this)[11]*x137 + (*this)[11]*x147 - (*this)[12]*x103 - (*this)[14]*x115 - (*this)[1]*x63 - (*this)[2]*A[6]*x43 + (*this)[2]*x148 + (*this)[3]*A[3]*x58 + (*this)[3]*x146 + (*this)[4]*x123 + (*this)[5]*x31 - (*this)[5]*x83 - (*this)[6]*x97 - (*this)[7]*A[3]*x46 - (*this)[7]*x29 - (*this)[7]*x75 - (*this)[8]*x96 + (*this)[9]*x36 + A[1]*x130 + A[5]*x0 - A[5]*x1 + A[5]*x17 - A[5]*x18 + A[5]*x19 + A[5]*x2 - A[5]*x20 - A[5]*x21 + A[5]*x22 - A[5]*x23 + A[5]*x24 - A[5]*x3 + A[5]*x4 - A[5]*x5 - A[5]*x6 + A[5]*x7 + x10*x111 + x107*x53 + x11*x111 - x111*x12 + x111*x13 - x111*x14 + x111*x15 - x111*x16 + x111*x9 - x120*x43 + x122*x69 - x127*x79 - x127*x81 - x134*x41 + x145*x39 + x32*x49 + x38*x48 - x46*x70 + x46*x73 - x69*x77;
    res[10]=-(*this)[11]*A[1]*x53 - (*this)[11]*x148 + (*this)[12]*A[1]*x67 - (*this)[12]*x63 - (*this)[14]*x124 + (*this)[15]*(*this)[8]*x111 - (*this)[15]*x83 + (*this)[1]*x103 + (*this)[2]*x147 + (*this)[3]*x125 + (*this)[4]*x139 + (*this)[6]*x29 - (*this)[6]*x37 - (*this)[7]*x104 - (*this)[7]*x40 - (*this)[8]*x143 + (*this)[8]*x27 - (*this)[9]*x26 + (*this)[9]*x99 + A[2]*x136 - A[5]*x138 + A[6]*x0 + A[6]*x1 + A[6]*x17 + A[6]*x18 - A[6]*x19 + A[6]*x2 - A[6]*x20 - A[6]*x21 + A[6]*x22 + A[6]*x23 - A[6]*x24 - A[6]*x3 - A[6]*x4 + A[6]*x5 - A[6]*x6 - A[6]*x7 - x10*x89 + x107*x58 + x11*x89 - x12*x89 + x121*x67 + x127*x82 - x13*x89 + x131*x69 + x14*x89 + x140*x43 + x15*x89 + x16*x89 - x25*x60 - x32*x70 + x32*x73 - x38*x69 - x46*x49 - x46*x51 + x89*x9;
    res[11]=-(*this)[0]*x150 + (*this)[11]*x149 + (*this)[13]*x156 - (*this)[14]*x151 - (*this)[15]*x153 + (*this)[3]*x161 + (*this)[4]*x152 + (*this)[5]*x159 + (*this)[7]*x162 - (*this)[8]*x155 - (*this)[9]*x157 - A[0]*x88 + A[0]*x94 + x154*x48 + x158*x58 + x160*x90;
    res[12]=(*this)[0]*x159 + (*this)[12]*x149 - (*this)[13]*x152 + (*this)[14]*x161 - (*this)[15]*x155 - (*this)[3]*x151 - (*this)[4]*x156 - (*this)[5]*x150 - (*this)[7]*x157 - (*this)[8]*x153 + (*this)[9]*x162 + A[0]*x106 + x154*x43 - x160*x69 + x163*x64 + x163*x80;
    res[13]=(*this)[0]*x164 + (*this)[10]*x159 + (*this)[13]*x149 - (*this)[15]*x157 - (*this)[1]*x161 - (*this)[5]*x165 + (*this)[7]*x155 - (*this)[8]*x162 - (*this)[9]*x153 + A[0]*x34 + A[0]*x42 - A[0]*x62 - A[0]*x74 + x154*x67 - x158*x90 + x163*x66;
    res[14]=(*this)[0]*x162 - (*this)[11]*x151 - (*this)[12]*x161 - (*this)[15]*x165 - (*this)[1]*x152 + (*this)[2]*x156 + (*this)[5]*x157 - (*this)[7]*x150 + (*this)[8]*x164 - (*this)[9]*x159 + A[0]*x28 - A[0]*x52 + A[0]*x57 + A[0]*x68 - x154*x78 + x158*x43;
    res[15]=A[7]*x0 + A[7]*x1 + A[7]*x17 + A[7]*x18 + A[7]*x19 - A[7]*x2 - A[7]*x20 + A[7]*x21 - A[7]*x22 - A[7]*x23 - A[7]*x24 - A[7]*x3 - A[7]*x4 - A[7]*x5 + A[7]*x6 + A[7]*x7 - x10*x163 - x11*x163 - x12*x163 - x13*x163 - x14*x163 - x15*x163 - x16*x163 + x163*x9;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*A[3];
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3;
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x7 + A[2]*x8 + x6*x9;
    res[2]=(*this)[3]*x5 + A[1]*x7 + A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x9;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x6 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*A[3];
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3;
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x7 + A[2]*x8 + x6*x9;
    res[2]=(*this)[3]*x5 + A[1]*x7 + A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x9;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x6 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
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
    T x2 = 2.0*(*this)[2];
    T x3 = A[2]*x2;
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[0]*x4;
    T x6 = (*this)[1]*x2;
    T x7 = (*this)[1]*A[2];
    T x8 = std::pow((*this)[2], 2);
    T x9 = std::pow((*this)[3], 2);
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x3 + A[0]*x0 + A[0]*x1 - A[0]*x8 - A[0]*x9 - A[1]*x5 + A[1]*x6 + x4*x7;
    res[1]=(*this)[3]*x3 + A[0]*x5 + A[0]*x6 + A[1]*x0 - A[1]*x1 + A[1]*x8 - A[1]*x9 - x10*x7;
    res[2]=-(*this)[0]*A[0]*x2 + (*this)[1]*A[0]*x4 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x2 + A[2]*x0 - A[2]*x1 - A[2]*x8 + A[2]*x9;
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
    T x2 = 2.0*(*this)[2];
    T x3 = A[2]*x2;
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[0]*x4;
    T x6 = (*this)[1]*x2;
    T x7 = (*this)[1]*A[2];
    T x8 = std::pow((*this)[2], 2);
    T x9 = std::pow((*this)[3], 2);
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[0]*x2;
    T x13 = (*this)[1]*x4;
    T x14 = (*this)[3]*x2;
    res[0]=(*this)[0]*x3 + A[0]*x0 + A[0]*x1 - A[0]*x8 - A[0]*x9 - A[1]*x5 + A[1]*x6 + x4*x7;
    res[1]=(*this)[3]*x3 + A[0]*x5 + A[0]*x6 + A[1]*x0 - A[1]*x1 + A[1]*x8 - A[1]*x9 - x10*x7;
    res[2]=-A[0]*x12 + A[0]*x13 + A[1]*x11 + A[1]*x14 + A[2]*x0 - A[2]*x1 - A[2]*x8 + A[2]*x9;
    res[3]=A[3]*x0 + A[3]*x1 - A[3]*x8 - A[3]*x9 - A[4]*x5 + A[4]*x6 + A[5]*x12 + A[5]*x13;
    res[4]=A[3]*x5 + A[3]*x6 + A[4]*x0 - A[4]*x1 + A[4]*x8 - A[4]*x9 - A[5]*x11 + A[5]*x14;
    res[5]=-A[3]*x12 + A[3]*x13 + A[4]*x11 + A[4]*x14 + A[5]*x0 - A[5]*x1 - A[5]*x8 + A[5]*x9;
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
    T x2 = 2.0*(*this)[2];
    T x3 = A[2]*x2;
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[0]*x4;
    T x6 = (*this)[1]*x2;
    T x7 = (*this)[1]*A[2];
    T x8 = std::pow((*this)[2], 2);
    T x9 = std::pow((*this)[3], 2);
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x3 + A[0]*x0 + A[0]*x1 - A[0]*x8 - A[0]*x9 - A[1]*x5 + A[1]*x6 + x4*x7;
    res[1]=(*this)[3]*x3 + A[0]*x5 + A[0]*x6 + A[1]*x0 - A[1]*x1 + A[1]*x8 - A[1]*x9 - x10*x7;
    res[2]=-(*this)[0]*A[0]*x2 + (*this)[1]*A[0]*x4 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x2 + A[2]*x0 - A[2]*x1 - A[2]*x8 + A[2]*x9;
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
    T x2 = 2.0*(*this)[2];
    T x3 = A[2]*x2;
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[0]*x4;
    T x6 = (*this)[1]*x2;
    T x7 = (*this)[1]*A[2];
    T x8 = std::pow((*this)[2], 2);
    T x9 = std::pow((*this)[3], 2);
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x3 + A[0]*x0 + A[0]*x1 - A[0]*x8 - A[0]*x9 - A[1]*x5 + A[1]*x6 + x4*x7;
    res[1]=(*this)[3]*x3 + A[0]*x5 + A[0]*x6 + A[1]*x0 - A[1]*x1 + A[1]*x8 - A[1]*x9 - x10*x7;
    res[2]=-(*this)[0]*A[0]*x2 + (*this)[1]*A[0]*x4 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x2 + A[2]*x0 - A[2]*x1 - A[2]*x8 + A[2]*x9;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*A[3];
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3;
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x7 + A[2]*x8 + x6*x9;
    res[2]=(*this)[3]*x5 + A[1]*x7 + A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x9;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x6 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
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
    T x2 = 2.0*(*this)[2];
    T x3 = A[2]*x2;
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[0]*x4;
    T x6 = (*this)[1]*x2;
    T x7 = (*this)[1]*A[2];
    T x8 = std::pow((*this)[2], 2);
    T x9 = std::pow((*this)[3], 2);
    T x10 = 2.0*(*this)[0];
    res[0]=(*this)[0]*x3 + A[0]*x0 + A[0]*x1 - A[0]*x8 - A[0]*x9 - A[1]*x5 + A[1]*x6 + x4*x7;
    res[1]=(*this)[3]*x3 + A[0]*x5 + A[0]*x6 + A[1]*x0 - A[1]*x1 + A[1]*x8 - A[1]*x9 - x10*x7;
    res[2]=-(*this)[0]*A[0]*x2 + (*this)[1]*A[0]*x4 + (*this)[1]*A[1]*x10 + (*this)[3]*A[1]*x2 + A[2]*x0 - A[2]*x1 - A[2]*x8 + A[2]*x9;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*A[4];
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[0]*x4;
    T x13 = (*this)[1]*x6;
    T x14 = (*this)[3]*x4;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3;
    res[1]=A[1]*x0 + A[1]*x1 + A[1]*x2 + A[1]*x3;
    res[2]=(*this)[0]*x5 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 - A[3]*x7 + A[3]*x8 + x6*x9;
    res[3]=(*this)[3]*x5 + A[2]*x7 + A[2]*x8 + A[3]*x0 - A[3]*x1 + A[3]*x2 - A[3]*x3 - x10*x9;
    res[4]=-A[2]*x12 + A[2]*x13 + A[3]*x11 + A[3]*x14 + A[4]*x0 - A[4]*x1 - A[4]*x2 + A[4]*x3;
    res[5]=A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 - A[6]*x7 + A[6]*x8 + A[7]*x12 + A[7]*x13;
    res[6]=A[5]*x7 + A[5]*x8 + A[6]*x0 - A[6]*x1 + A[6]*x2 - A[6]*x3 - A[7]*x11 + A[7]*x14;
    res[7]=-A[5]*x12 + A[5]*x13 + A[6]*x11 + A[6]*x14 + A[7]*x0 - A[7]*x1 - A[7]*x2 + A[7]*x3;
    res[8]=A[10]*x12 + A[10]*x13 + A[8]*x0 + A[8]*x1 - A[8]*x2 - A[8]*x3 - A[9]*x7 + A[9]*x8;
    res[9]=-A[10]*x11 + A[10]*x14 + A[8]*x7 + A[8]*x8 + A[9]*x0 - A[9]*x1 + A[9]*x2 - A[9]*x3;
    res[10]=A[10]*x0 - A[10]*x1 - A[10]*x2 + A[10]*x3 - A[8]*x12 + A[8]*x13 + A[9]*x11 + A[9]*x14;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x2 + A[11]*x3;
    res[12]=A[12]*x0 + A[12]*x1 - A[12]*x2 - A[12]*x3 - A[13]*x7 + A[13]*x8 + A[14]*x12 + A[14]*x13;
    res[13]=A[12]*x7 + A[12]*x8 + A[13]*x0 - A[13]*x1 + A[13]*x2 - A[13]*x3 - A[14]*x11 + A[14]*x14;
    res[14]=-A[12]*x12 + A[12]*x13 + A[13]*x11 + A[13]*x14 + A[14]*x0 - A[14]*x1 - A[14]*x2 + A[14]*x3;
    res[15]=A[15]*x0 + A[15]*x1 + A[15]*x2 + A[15]*x3;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*A[3];
    T x10 = 2.0*(*this)[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3;
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x7 + A[2]*x8 + x6*x9;
    res[2]=(*this)[3]*x5 + A[1]*x7 + A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x9;
    res[3]=-(*this)[0]*A[1]*x4 + (*this)[1]*A[1]*x6 + (*this)[1]*A[2]*x10 + (*this)[3]*A[2]*x4 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
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
    T x6 = 2.0*(*this)[3];
    T x7 = (*this)[0]*x6;
    T x8 = (*this)[1]*x4;
    T x9 = (*this)[1]*A[3];
    T x10 = 2.0*(*this)[0];
    T x11 = (*this)[1]*x10;
    T x12 = (*this)[0]*x4;
    T x13 = (*this)[1]*x6;
    T x14 = (*this)[3]*x4;
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3;
    res[1]=(*this)[0]*x5 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[2]*x7 + A[2]*x8 + x6*x9;
    res[2]=(*this)[3]*x5 + A[1]*x7 + A[1]*x8 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - x10*x9;
    res[3]=-A[1]*x12 + A[1]*x13 + A[2]*x11 + A[2]*x14 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3;
    res[4]=A[4]*x0 + A[4]*x1 - A[4]*x2 - A[4]*x3 - A[5]*x7 + A[5]*x8 + A[6]*x12 + A[6]*x13;
    res[5]=A[4]*x7 + A[4]*x8 + A[5]*x0 - A[5]*x1 + A[5]*x2 - A[5]*x3 - A[6]*x11 + A[6]*x14;
    res[6]=-A[4]*x12 + A[4]*x13 + A[5]*x11 + A[5]*x14 + A[6]*x0 - A[6]*x1 - A[6]*x2 + A[6]*x3;
    res[7]=A[7]*x0 + A[7]*x1 + A[7]*x2 + A[7]*x3;
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
    T x10 = 2.0*A[2];
    T x11 = (*this)[0]*(*this)[6];
    T x12 = (*this)[1]*(*this)[2];
    T x13 = (*this)[3]*x8;
    T x14 = (*this)[7]*x8;
    T x15 = (*this)[3]*(*this)[7];
    T x16 = (*this)[4]*(*this)[5];
    T x17 = (*this)[4]*x8;
    T x18 = 2.0*A[1];
    T x19 = (*this)[0]*x10;
    T x20 = (*this)[0]*x18;
    T x21 = (*this)[1]*x18;
    T x22 = (*this)[1]*x10;
    T x23 = (*this)[2]*x10;
    T x24 = (*this)[2]*x18;
    T x25 = (*this)[6]*x18;
    T x26 = (*this)[6]*x10;
    T x27 = (*this)[0]*x8;
    T x28 = (*this)[6]*x8;
    T x29 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7;
    res[1]=(*this)[0]*x9 - (*this)[1]*x13 + (*this)[2]*x14 + (*this)[6]*x17 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7 - x10*x11 - x10*x12 - x10*x15 + x10*x16;
    res[2]=-(*this)[0]*x17 - (*this)[1]*x14 - (*this)[2]*x13 + (*this)[6]*x9 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7 + x11*x18 - x12*x18 + x15*x18 + x16*x18;
    res[3]=-(*this)[3]*x21 - (*this)[3]*x23 + (*this)[4]*x19 + (*this)[4]*x25 - (*this)[5]*x20 + (*this)[5]*x26 + (*this)[7]*x22 - (*this)[7]*x24 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7;
    res[4]=-(*this)[1]*x28 - (*this)[2]*x27 - (*this)[3]*x17 + (*this)[3]*x19 + (*this)[3]*x25 - (*this)[4]*x21 - (*this)[4]*x23 - (*this)[5]*x22 + (*this)[5]*x24 + (*this)[7]*x20 - (*this)[7]*x26 + (*this)[7]*x9;
    res[5]=(*this)[1]*x27 - (*this)[2]*x28 - (*this)[3]*x20 + (*this)[3]*x26 - (*this)[3]*x9 - (*this)[4]*x14 + (*this)[4]*x22 - (*this)[4]*x24 - (*this)[5]*x21 - (*this)[5]*x23 + (*this)[7]*x19 + (*this)[7]*x25;
    res[6]=(*this)[0]*x14 + (*this)[1]*x17 - (*this)[1]*x19 - (*this)[1]*x25 + (*this)[2]*x20 - (*this)[2]*x26 + (*this)[2]*x9 - (*this)[3]*(*this)[4]*x18 - (*this)[3]*(*this)[5]*x10 + (*this)[4]*(*this)[7]*x10 - (*this)[5]*(*this)[7]*x18 - (*this)[6]*x13;
    res[7]=(*this)[0]*(*this)[7]*x29 - (*this)[1]*(*this)[4]*x29 - (*this)[2]*(*this)[5]*x29 - (*this)[3]*(*this)[6]*x29;
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
    T x9 = (*this)[0]*x8;
    T x10 = 2.0*(*this)[2];
    T x11 = (*this)[0]*A[2];
    T x12 = 2.0*(*this)[3];
    T x13 = (*this)[0]*A[3];
    T x14 = A[3]*x8;
    T x15 = A[2]*x8;
    T x16 = (*this)[4]*A[3];
    T x17 = (*this)[6]*x10;
    T x18 = A[2]*x12;
    T x19 = (*this)[5]*x12;
    T x20 = 2.0*(*this)[7];
    T x21 = (*this)[4]*x20;
    T x22 = (*this)[5]*A[2];
    T x23 = (*this)[6]*A[3];
    T x24 = 2.0*(*this)[5];
    T x25 = 2.0*(*this)[6];
    T x26 = A[3]*x10;
    T x27 = 2.0*(*this)[4];
    T x28 = A[0]*x10;
    T x29 = (*this)[0]*A[1];
    T x30 = A[1]*x8;
    T x31 = A[0]*x8;
    T x32 = A[0]*x12;
    T x33 = (*this)[7]*A[1];
    T x34 = A[0]*x20;
    res[0]=-(*this)[4]*x18 - (*this)[5]*x14 + (*this)[6]*x15 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[1]*x17 + A[1]*x19 - A[1]*x21 - A[1]*x9 - x10*x11 + x10*x16 - x12*x13 - x20*x22 - x20*x23;
    res[1]=(*this)[2]*x15 + (*this)[3]*x14 + (*this)[7]*x18 - (*this)[7]*x26 + A[0]*x17 - A[0]*x19 - A[0]*x21 - A[0]*x9 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 - A[1]*x5 - A[1]*x6 + A[1]*x7 - x11*x25 + x13*x24 + x16*x25 + x22*x27;
    res[2]=-2.0*(*this)[0]*x16 - (*this)[0]*x28 + (*this)[2]*x30 + (*this)[3]*x26 + (*this)[4]*A[1]*x24 + (*this)[4]*x32 - (*this)[5]*x34 - (*this)[6]*x31 + (*this)[7]*x14 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 - x12*x33 + x23*x24 + x25*x29;
    res[3]=-(*this)[0]*x32 + (*this)[3]*A[2]*x10 + (*this)[3]*x30 - (*this)[4]*x28 + (*this)[5]*x31 + (*this)[6]*A[1]*x27 - (*this)[6]*x34 - (*this)[7]*x15 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 - A[3]*x4 - A[3]*x5 + A[3]*x6 + A[3]*x7 + x10*x33 + x11*x27 + x22*x25 - x24*x29;
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
    T x1 = (*this)[0]*A[0];
    T x2 = 2.0*(*this)[2];
    T x3 = (*this)[0]*A[1];
    T x4 = 2.0*(*this)[3];
    T x5 = (*this)[0]*A[2];
    T x6 = A[2]*x0;
    T x7 = A[1]*x0;
    T x8 = (*this)[4]*A[2];
    T x9 = A[0]*x2;
    T x10 = A[1]*x4;
    T x11 = (*this)[5]*A[0];
    T x12 = 2.0*(*this)[7];
    T x13 = (*this)[4]*A[0];
    T x14 = (*this)[5]*A[1];
    T x15 = (*this)[6]*A[2];
    T x16 = std::pow((*this)[0], 2);
    T x17 = std::pow((*this)[1], 2);
    T x18 = std::pow((*this)[4], 2);
    T x19 = std::pow((*this)[7], 2);
    T x20 = 2.0*(*this)[5];
    T x21 = 2.0*(*this)[6];
    T x22 = A[2]*x2;
    T x23 = 2.0*(*this)[4];
    T x24 = std::pow((*this)[2], 2);
    T x25 = std::pow((*this)[3], 2);
    T x26 = std::pow((*this)[5], 2);
    T x27 = std::pow((*this)[6], 2);
    T x28 = 2.0*(*this)[0];
    T x29 = A[0]*x0;
    res[0]=-(*this)[4]*x10 - (*this)[5]*x6 + (*this)[6]*x7 - (*this)[6]*x9 - x0*x1 + x11*x4 - x12*x13 - x12*x14 - x12*x15 - x2*x3 + x2*x8 - x4*x5;
    res[1]=(*this)[2]*x7 + (*this)[3]*x6 + (*this)[7]*x10 - (*this)[7]*x22 + A[0]*x16 + A[0]*x17 + A[0]*x18 + A[0]*x19 - A[0]*x24 - A[0]*x25 - A[0]*x26 - A[0]*x27 + x14*x23 + x20*x5 - x21*x3 + x21*x8;
    res[2]=(*this)[2]*x29 + (*this)[3]*x22 - (*this)[7]*A[0]*x4 + (*this)[7]*x6 + A[1]*x16 - A[1]*x17 - A[1]*x18 + A[1]*x19 + A[1]*x24 - A[1]*x25 + A[1]*x26 - A[1]*x27 + x1*x21 + x11*x23 + x15*x20 - x28*x8;
    res[3]=(*this)[3]*A[1]*x2 + (*this)[3]*x29 - (*this)[7]*x7 + (*this)[7]*x9 + A[2]*x16 - A[2]*x17 - A[2]*x18 + A[2]*x19 - A[2]*x24 + A[2]*x25 - A[2]*x26 + A[2]*x27 - x11*x28 + x13*x21 + x14*x21 + x23*x3;
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
    res[1]=-x0*((*this)[0]*(*this)[1] - (*this)[2]*(*this)[6] + (*this)[3]*(*this)[5] + (*this)[4]*(*this)[7]);
    res[2]=-x0*((*this)[0]*(*this)[2] + (*this)[1]*(*this)[6] - (*this)[3]*(*this)[4] + (*this)[5]*(*this)[7]);
    res[3]=-x0*((*this)[0]*(*this)[3] - (*this)[1]*(*this)[5] + (*this)[2]*(*this)[4] + (*this)[6]*(*this)[7]);
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
    T x4 = 2.0*(*this)[0];
    T x5 = (*this)[2]*A[5];
    T x6 = (*this)[3]*x4;
    T x7 = (*this)[5]*A[2];
    T x8 = (*this)[6]*x4;
    T x9 = (*this)[7]*A[3];
    T x10 = 2.0*(*this)[1];
    T x11 = (*this)[2]*x10;
    T x12 = (*this)[3]*A[2];
    T x13 = A[3]*x10;
    T x14 = (*this)[5]*A[4];
    T x15 = (*this)[6]*A[5];
    T x16 = 2.0*(*this)[2];
    T x17 = (*this)[4]*A[4];
    T x18 = A[3]*x16;
    T x19 = (*this)[7]*A[2];
    T x20 = 2.0*(*this)[4];
    T x21 = (*this)[3]*A[5];
    T x22 = 2.0*(*this)[6];
    T x23 = (*this)[3]*x22;
    T x24 = 2.0*(*this)[7];
    T x25 = (*this)[3]*A[1];
    T x26 = (*this)[5]*x20;
    T x27 = (*this)[6]*x20;
    T x28 = (*this)[5]*A[5];
    T x29 = (*this)[7]*A[4];
    T x30 = std::pow((*this)[1], 2);
    T x31 = std::pow((*this)[5], 2);
    T x32 = std::pow((*this)[6], 2);
    T x33 = std::pow((*this)[7], 2);
    T x34 = (*this)[1]*x4;
    T x35 = (*this)[4]*x4;
    T x36 = 2.0*(*this)[5];
    T x37 = (*this)[3]*A[0];
    T x38 = (*this)[7]*A[5];
    T x39 = (*this)[2]*x4;
    T x40 = A[0]*x4;
    T x41 = (*this)[4]*x10;
    T x42 = A[1]*x10;
    T x43 = (*this)[6]*x16;
    T x44 = A[0]*x16;
    T x45 = 2.0*(*this)[3];
    T x46 = A[1]*x22;
    T x47 = (*this)[6]*x10;
    T x48 = A[1]*x16;
    T x49 = (*this)[7]*A[1];
    T x50 = (*this)[5]*A[0];
    res[0]=(*this)[4]*x13 - (*this)[5]*x18 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x30 - A[0]*x31 - A[0]*x32 - A[0]*x33 - A[1]*x11 + A[1]*x26 - A[1]*x8 + A[2]*x27 - A[3]*x23 - A[4]*x6 - x10*x12 + x10*x14 + x10*x15 + x16*x17 + x16*x19 + x20*x21 + x22*x29 - x24*x25 - x24*x28 + x4*x5 + x4*x7 - x4*x9;
    res[1]=(*this)[4]*x18 + (*this)[5]*x13 - A[0]*x11 + A[0]*x26 + A[0]*x8 + A[1]*x0 - A[1]*x1 + A[1]*x2 - A[1]*x3 + A[1]*x30 + A[1]*x31 - A[1]*x32 - A[1]*x33 - A[2]*x35 + A[3]*x6 - A[4]*x23 - A[5]*x34 - x10*x17 - x10*x19 - x12*x16 + x14*x16 + x20*x38 + x21*x36 + x22*x5 + x22*x7 - x22*x9 + x24*x37 - x29*x4;
    res[2]=(*this)[3]*A[3]*x20 - (*this)[5]*x40 + (*this)[5]*x46 + (*this)[6]*x13 + (*this)[7]*x42 - (*this)[7]*x44 + A[0]*x27 + A[1]*x35 + A[2]*x0 + A[2]*x1 - A[2]*x2 - A[2]*x3 + A[2]*x30 - A[2]*x31 + A[2]*x32 - A[2]*x33 - A[3]*x39 + A[4]*x34 + A[4]*x43 - A[5]*x41 - x10*x37 + x14*x45 + x15*x45 - x16*x25 - x17*x24 - x36*x5 + x36*x9 - x38*x4;
    res[3]=-(*this)[4]*x48 - (*this)[5]*x42 + (*this)[5]*x44 + (*this)[7]*x40 - (*this)[7]*x46 + A[0]*x23 - A[0]*x41 + A[1]*x6 - A[2]*x39 - A[2]*x47 + A[3]*x0 + A[3]*x1 + A[3]*x2 + A[3]*x3 - A[3]*x30 - A[3]*x31 - A[3]*x32 - A[3]*x33 - A[4]*x11 - A[4]*x8 - x10*x21 - x12*x20 + x14*x20 + x15*x20 + x24*x5 + x24*x7 + x28*x4 - x29*x45;
    res[4]=-(*this)[2]*x13 + (*this)[4]*x42 - (*this)[4]*x44 - (*this)[5]*x48 + (*this)[7]*A[0]*x22 - A[0]*x6 + A[1]*x23 + A[2]*x34 - A[2]*x43 + A[3]*x26 + A[3]*x8 + A[4]*x0 - A[4]*x1 + A[4]*x2 - A[4]*x3 + A[4]*x30 + A[4]*x31 - A[4]*x32 - A[4]*x33 - A[5]*x35 - x10*x38 - x10*x50 + x15*x36 - x19*x20 + x4*x49 - x45*x5 - x45*x7 + x45*x9;
    res[5]=-(*this)[3]*A[4]*x16 - (*this)[3]*x13 - (*this)[5]*A[3]*x4 + A[0]*x39 - A[0]*x47 - A[1]*x34 - A[1]*x43 + A[2]*x41 + A[3]*x27 + A[5]*x0 + A[5]*x1 - A[5]*x2 - A[5]*x3 + A[5]*x30 - A[5]*x31 + A[5]*x32 - A[5]*x33 + x10*x29 - x12*x22 + x14*x22 + x16*x7 - x16*x9 + x17*x4 + x19*x4 - x20*x37 + x20*x49 - x24*x50 - x25*x36;
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
    T x2 = 2.0*A[1];
    T x3 = (*this)[0]*(*this)[3];
    T x4 = 2.0*A[0];
    T x5 = (*this)[0]*(*this)[7];
    T x6 = 2.0*(*this)[1];
    T x7 = A[0]*x6;
    T x8 = (*this)[5]*A[1];
    T x9 = (*this)[1]*x0;
    T x10 = (*this)[2]*(*this)[4];
    T x11 = (*this)[5]*x4;
    T x12 = (*this)[3]*x0;
    T x13 = (*this)[3]*(*this)[6];
    T x14 = (*this)[7]*x0;
    T x15 = (*this)[7]*x2;
    T x16 = A[1]*x6;
    T x17 = 2.0*x8;
    T x18 = (*this)[6]*x4;
    T x19 = (*this)[2]*x4;
    T x20 = (*this)[6]*x2;
    T x21 = (*this)[3]*x4;
    T x22 = std::pow((*this)[0], 2);
    T x23 = std::pow((*this)[2], 2);
    T x24 = std::pow((*this)[3], 2);
    T x25 = std::pow((*this)[4], 2);
    T x26 = (*this)[5]*x0;
    T x27 = (*this)[4]*x0;
    T x28 = std::pow((*this)[1], 2);
    T x29 = std::pow((*this)[5], 2);
    T x30 = std::pow((*this)[6], 2);
    T x31 = std::pow((*this)[7], 2);
    res[0]=(*this)[0]*x1 - (*this)[2]*x11 + (*this)[4]*x12 + (*this)[4]*x7 - (*this)[5]*x14 + (*this)[6]*x15 + (*this)[6]*x9 + x10*x2 - x13*x4 - x2*x3 - x4*x5 + x6*x8;
    res[1]=-(*this)[0]*x15 - (*this)[0]*x9 + (*this)[2]*x17 + (*this)[4]*x14 - (*this)[4]*x16 + (*this)[5]*x12 + (*this)[5]*x7 + (*this)[6]*x1 - (*this)[7]*x18 + x10*x4 - x13*x2 + x3*x4;
    res[2]=(*this)[0]*x16 - (*this)[0]*x19 + (*this)[2]*x20 + (*this)[3]*x17 - (*this)[4]*x15 + (*this)[4]*x21 - (*this)[4]*x9 - (*this)[5]*x1 + (*this)[6]*x12 + (*this)[6]*x7 + (*this)[7]*x11 - x0*x5;
    res[3]=-(*this)[0]*x20 + (*this)[0]*x26 - (*this)[2]*x16 - (*this)[3]*x15 - (*this)[3]*x9 + (*this)[4]*x17 + (*this)[6]*x27 + (*this)[7]*x1 + A[0]*x22 + A[0]*x23 + A[0]*x24 + A[0]*x25 - A[0]*x28 - A[0]*x29 - A[0]*x30 - A[0]*x31;
    res[4]=(*this)[0]*x18 - (*this)[0]*x27 - (*this)[2]*x7 - (*this)[3]*x1 + (*this)[4]*x11 + (*this)[6]*x26 + (*this)[7]*x21 - (*this)[7]*x9 + A[1]*x22 - A[1]*x23 + A[1]*x24 - A[1]*x25 + A[1]*x28 + A[1]*x29 - A[1]*x30 - A[1]*x31;
    res[5]=(*this)[0]*(*this)[4]*x2 - (*this)[0]*x11 - (*this)[2]*(*this)[3]*x2 - (*this)[3]*x7 + (*this)[4]*x18 + (*this)[6]*x17 + (*this)[7]*x16 - (*this)[7]*x19 + A[2]*x22 + A[2]*x23 - A[2]*x24 - A[2]*x25 + A[2]*x28 - A[2]*x29 + A[2]*x30 - A[2]*x31;
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
    T x4 = 2.0*A[2];
    T x5 = (*this)[5]*x4;
    T x6 = 2.0*A[1];
    T x7 = (*this)[0]*(*this)[6];
    T x8 = (*this)[1]*(*this)[2];
    T x9 = (*this)[3]*x4;
    T x10 = (*this)[7]*x4;
    T x11 = (*this)[3]*(*this)[7];
    T x12 = (*this)[4]*(*this)[5];
    T x13 = (*this)[4]*x4;
    T x14 = std::pow((*this)[1], 2);
    T x15 = std::pow((*this)[5], 2);
    T x16 = std::pow((*this)[6], 2);
    T x17 = std::pow((*this)[7], 2);
    T x18 = 2.0*A[0];
    T x19 = (*this)[0]*x6;
    T x20 = (*this)[0]*x18;
    T x21 = (*this)[1]*x18;
    T x22 = (*this)[1]*x6;
    T x23 = (*this)[2]*x6;
    T x24 = (*this)[2]*x18;
    T x25 = (*this)[6]*x18;
    T x26 = (*this)[6]*x6;
    T x27 = (*this)[0]*x4;
    T x28 = (*this)[6]*x4;
    res[0]=(*this)[0]*x5 - (*this)[1]*x9 + (*this)[2]*x10 + (*this)[6]*x13 + A[0]*x0 + A[0]*x1 - A[0]*x14 - A[0]*x15 - A[0]*x16 - A[0]*x17 + A[0]*x2 + A[0]*x3 - x11*x6 + x12*x6 - x6*x7 - x6*x8;
    res[1]=-(*this)[0]*x13 - (*this)[1]*x10 - (*this)[2]*x9 + (*this)[6]*x5 + A[1]*x0 - A[1]*x1 + A[1]*x14 + A[1]*x15 - A[1]*x16 - A[1]*x17 + A[1]*x2 - A[1]*x3 + x11*x18 + x12*x18 + x18*x7 - x18*x8;
    res[2]=-(*this)[3]*x21 - (*this)[3]*x23 + (*this)[4]*x19 + (*this)[4]*x25 - (*this)[5]*x20 + (*this)[5]*x26 + (*this)[7]*x22 - (*this)[7]*x24 + A[2]*x0 + A[2]*x1 + A[2]*x14 - A[2]*x15 + A[2]*x16 - A[2]*x17 - A[2]*x2 - A[2]*x3;
    res[3]=-(*this)[1]*x28 - (*this)[2]*x27 - (*this)[3]*x13 + (*this)[3]*x19 + (*this)[3]*x25 - (*this)[4]*x21 - (*this)[4]*x23 - (*this)[5]*x22 + (*this)[5]*x24 + (*this)[7]*x20 - (*this)[7]*x26 + (*this)[7]*x5;
    res[4]=(*this)[1]*x27 - (*this)[2]*x28 - (*this)[3]*x20 + (*this)[3]*x26 - (*this)[3]*x5 - (*this)[4]*x10 + (*this)[4]*x22 - (*this)[4]*x24 - (*this)[5]*x21 - (*this)[5]*x23 + (*this)[7]*x19 + (*this)[7]*x25;
    res[5]=(*this)[0]*x10 + (*this)[1]*x13 - (*this)[1]*x19 - (*this)[1]*x25 + (*this)[2]*x20 - (*this)[2]*x26 + (*this)[2]*x5 - (*this)[3]*(*this)[4]*x18 - (*this)[3]*(*this)[5]*x6 + (*this)[4]*(*this)[7]*x6 - (*this)[5]*(*this)[7]*x18 - (*this)[6]*x9;
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
    T x14 = 2.0*A[2];
    T x15 = (*this)[1]*(*this)[6];
    T x16 = 2.0*A[3];
    T x17 = (*this)[2]*x16;
    T x18 = 2.0*A[1];
    T x19 = (*this)[2]*(*this)[6];
    T x20 = (*this)[3]*x14;
    T x21 = (*this)[3]*x12;
    T x22 = (*this)[4]*(*this)[7];
    T x23 = A[2]*x12;
    T x24 = (*this)[6]*x16;
    T x25 = (*this)[6]*x8;
    T x26 = 2.0*x10;
    T x27 = 2.0*x13;
    T x28 = 2.0*A[0];
    T x29 = A[0]*x8;
    T x30 = (*this)[1]*x18;
    T x31 = (*this)[4]*x28;
    T x32 = (*this)[7]*x18;
    T x33 = A[0]*x12;
    res[0]=(*this)[3]*x11 - (*this)[4]*x17 + (*this)[4]*x20 + (*this)[7]*x23 + (*this)[7]*x24 + A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 + A[0]*x4 + A[0]*x5 + A[0]*x6 + A[0]*x7 - A[1]*x21 + A[1]*x9 + x10*x8 + x12*x13 - x14*x15 + x18*x19 + x18*x22;
    res[1]=(*this)[1]*x26 + (*this)[3]*x27 + (*this)[4]*x23 + (*this)[4]*x24 + (*this)[5]*x11 - (*this)[7]*x17 + (*this)[7]*x20 + A[0]*x21 + A[0]*x9 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 + A[1]*x4 - A[1]*x5 - A[1]*x6 + A[1]*x7 - A[2]*x25 - x19*x28 + x22*x28;
    res[2]=(*this)[2]*x29 + (*this)[2]*x30 + (*this)[3]*x17 - (*this)[3]*x31 - (*this)[3]*x32 + (*this)[4]*A[1]*x12 - (*this)[4]*x11 + (*this)[6]*A[3]*x12 + (*this)[7]*x27 + (*this)[7]*x33 + A[1]*x25 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 - A[2]*x4 + A[2]*x5 - A[2]*x6 + A[2]*x7 + x15*x28;
    res[3]=-(*this)[1]*(*this)[7]*x14 - (*this)[1]*x33 + (*this)[2]*x31 + (*this)[2]*x32 + (*this)[3]*x26 + (*this)[3]*x29 + (*this)[3]*x30 + (*this)[4]*(*this)[6]*x18 + (*this)[4]*A[2]*x8 - (*this)[5]*A[1]*x8 + (*this)[6]*(*this)[7]*x28 + (*this)[6]*x23 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 - A[3]*x4 - A[3]*x5 + A[3]*x6 + A[3]*x7;
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
    T x6 = 2.0*A[1];
    T x7 = (*this)[1]*x6;
    T x8 = 2.0*A[2];
    T x9 = (*this)[2]*x8;
    T x10 = 2.0*A[0];
    T x11 = (*this)[2]*x10;
    T x12 = (*this)[3]*x6;
    T x13 = A[0]*x4;
    T x14 = (*this)[7]*x10;
    T x15 = A[1]*x4;
    T x16 = (*this)[6]*x8;
    T x17 = std::pow((*this)[0], 2);
    T x18 = std::pow((*this)[1], 2);
    T x19 = std::pow((*this)[4], 2);
    T x20 = std::pow((*this)[7], 2);
    T x21 = A[1]*x0;
    T x22 = 2.0*x2;
    T x23 = 2.0*x5;
    T x24 = std::pow((*this)[2], 2);
    T x25 = std::pow((*this)[3], 2);
    T x26 = std::pow((*this)[5], 2);
    T x27 = std::pow((*this)[6], 2);
    res[0]=(*this)[1]*x1 - (*this)[3]*x13 + (*this)[3]*x3 + (*this)[4]*x12 + (*this)[4]*x14 - (*this)[4]*x9 + (*this)[6]*x11 - (*this)[6]*x7 + (*this)[7]*x15 + (*this)[7]*x16 + x0*x2 + x4*x5;
    res[1]=(*this)[1]*x22 + (*this)[3]*x23 + (*this)[4]*x15 + (*this)[4]*x16 + (*this)[5]*x3 - (*this)[6]*x21 + (*this)[7]*x12 - (*this)[7]*x9 + A[0]*x17 + A[0]*x18 + A[0]*x19 + A[0]*x20 - A[0]*x24 - A[0]*x25 - A[0]*x26 - A[0]*x27;
    res[2]=(*this)[1]*x11 - (*this)[3]*x14 + (*this)[3]*x9 + (*this)[4]*x13 - (*this)[4]*x3 + (*this)[6]*A[2]*x4 + (*this)[6]*x1 + (*this)[7]*x23 + A[1]*x17 - A[1]*x18 - A[1]*x19 + A[1]*x20 + A[1]*x24 - A[1]*x25 + A[1]*x26 - A[1]*x27;
    res[3]=(*this)[1]*(*this)[3]*x10 + (*this)[3]*x22 + (*this)[4]*(*this)[6]*x10 + (*this)[4]*x21 - (*this)[5]*x1 + (*this)[6]*x15 + (*this)[7]*x11 - (*this)[7]*x7 + A[2]*x17 - A[2]*x18 - A[2]*x19 + A[2]*x20 - A[2]*x24 + A[2]*x25 - A[2]*x26 + A[2]*x27;
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
    res[0]=-2.0*A[0]*((*this)[0]*(*this)[7] - (*this)[1]*(*this)[4] - (*this)[2]*(*this)[5] - (*this)[3]*(*this)[6]);
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
    T x4 = 2.0*A[15];
    T x5 = (*this)[0]*(*this)[7];
    T x6 = (*this)[1]*(*this)[4];
    T x7 = (*this)[2]*(*this)[5];
    T x8 = (*this)[3]*(*this)[6];
    T x9 = std::pow((*this)[1], 2);
    T x10 = std::pow((*this)[2], 2);
    T x11 = std::pow((*this)[3], 2);
    T x12 = std::pow((*this)[7], 2);
    T x13 = 2.0*(*this)[1];
    T x14 = (*this)[0]*x13;
    T x15 = 2.0*(*this)[2];
    T x16 = (*this)[0]*A[3];
    T x17 = 2.0*(*this)[3];
    T x18 = (*this)[0]*A[4];
    T x19 = A[4]*x13;
    T x20 = A[3]*x13;
    T x21 = (*this)[4]*A[4];
    T x22 = (*this)[6]*x15;
    T x23 = A[3]*x17;
    T x24 = (*this)[5]*x17;
    T x25 = 2.0*(*this)[7];
    T x26 = (*this)[4]*x25;
    T x27 = (*this)[5]*A[3];
    T x28 = (*this)[6]*A[4];
    T x29 = 2.0*(*this)[5];
    T x30 = 2.0*(*this)[6];
    T x31 = A[4]*x15;
    T x32 = 2.0*(*this)[4];
    T x33 = A[1]*x15;
    T x34 = (*this)[0]*A[2];
    T x35 = A[2]*x13;
    T x36 = A[1]*x13;
    T x37 = A[1]*x17;
    T x38 = (*this)[7]*A[2];
    T x39 = (*this)[4]*x29;
    T x40 = A[1]*x25;
    T x41 = (*this)[3]*x15;
    T x42 = (*this)[6]*x32;
    T x43 = (*this)[0]*x15;
    T x44 = (*this)[0]*x17;
    T x45 = A[7]*x29;
    T x46 = (*this)[0]*x30;
    T x47 = 2.0*A[8];
    T x48 = (*this)[2]*x13;
    T x49 = (*this)[3]*x13;
    T x50 = (*this)[5]*x13;
    T x51 = (*this)[6]*x13;
    T x52 = (*this)[4]*x15;
    T x53 = (*this)[7]*A[7];
    T x54 = (*this)[4]*x17;
    T x55 = (*this)[7]*x17;
    T x56 = (*this)[5]*x25;
    T x57 = (*this)[6]*x25;
    T x58 = (*this)[0]*x32;
    T x59 = 2.0*A[9];
    T x60 = (*this)[0]*x29;
    T x61 = 2.0*A[10];
    T x62 = (*this)[7]*x13;
    T x63 = (*this)[7]*x15;
    T x64 = (*this)[6]*x29;
    T x65 = 2.0*A[5];
    T x66 = 2.0*A[6];
    T x67 = 2.0*A[7];
    T x68 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 + A[0]*x2 + A[0]*x3 - A[0]*x9 - x4*x5 + x4*x6 + x4*x7 + x4*x8;
    res[1]=-(*this)[4]*x23 - (*this)[5]*x19 + (*this)[6]*x20 + A[1]*x0 + A[1]*x1 + A[1]*x10 + A[1]*x11 + A[1]*x12 + A[1]*x2 + A[1]*x3 + A[1]*x9 - A[2]*x14 - A[2]*x22 + A[2]*x24 - A[2]*x26 - x15*x16 + x15*x21 - x17*x18 - x25*x27 - x25*x28;
    res[2]=(*this)[2]*x20 + (*this)[3]*x19 + (*this)[7]*x23 - (*this)[7]*x31 - A[1]*x14 + A[1]*x22 - A[1]*x24 - A[1]*x26 + A[2]*x0 + A[2]*x1 - A[2]*x10 - A[2]*x11 + A[2]*x12 - A[2]*x2 - A[2]*x3 + A[2]*x9 - x16*x30 + x18*x29 + x21*x30 + x27*x32;
    res[3]=-2.0*(*this)[0]*x21 - (*this)[0]*x33 + (*this)[2]*x35 + (*this)[3]*x31 + (*this)[4]*x37 - (*this)[5]*x40 - (*this)[6]*x36 + (*this)[7]*x19 + A[2]*x39 + A[3]*x0 - A[3]*x1 + A[3]*x10 - A[3]*x11 + A[3]*x12 + A[3]*x2 - A[3]*x3 - A[3]*x9 - x17*x38 + x28*x29 + x30*x34;
    res[4]=-(*this)[0]*x37 + (*this)[3]*x35 - (*this)[4]*x33 + (*this)[5]*x36 - (*this)[6]*x40 - (*this)[7]*x20 + A[2]*x42 + A[3]*x41 + A[4]*x0 - A[4]*x1 - A[4]*x10 + A[4]*x11 + A[4]*x12 - A[4]*x2 + A[4]*x3 - A[4]*x9 + x15*x38 + x16*x32 + x27*x30 - x29*x34;
    res[5]=(*this)[0]*x45 + A[10]*x43 + A[10]*x51 + A[10]*x54 - A[10]*x56 + A[5]*x0 + A[5]*x1 + A[5]*x10 + A[5]*x11 - A[5]*x12 - A[5]*x2 - A[5]*x3 - A[5]*x9 + A[6]*x39 - A[6]*x46 - A[6]*x48 - A[6]*x55 + A[7]*x42 - A[7]*x49 - A[9]*x44 + A[9]*x50 + A[9]*x52 + A[9]*x57 + x15*x53 - x47*x5 + x47*x6 - x47*x7 - x47*x8;
    res[6]=(*this)[6]*x45 - A[10]*x14 + A[10]*x22 + A[10]*x24 + A[10]*x26 + A[5]*x39 + A[5]*x46 - A[5]*x48 + A[5]*x55 + A[6]*x0 - A[6]*x1 - A[6]*x10 + A[6]*x11 - A[6]*x12 + A[6]*x2 - A[6]*x3 + A[6]*x9 - A[7]*x41 - A[7]*x58 + A[8]*x44 + A[8]*x50 + A[8]*x52 - A[8]*x57 - x13*x53 - x5*x59 - x59*x6 + x59*x7 - x59*x8;
    res[7]=A[5]*x42 - A[5]*x49 - A[5]*x60 - A[5]*x63 - A[6]*x41 + A[6]*x58 + A[6]*x62 + A[6]*x64 + A[7]*x0 - A[7]*x1 + A[7]*x10 - A[7]*x11 - A[7]*x12 - A[7]*x2 + A[7]*x3 + A[7]*x9 - A[8]*x43 + A[8]*x51 + A[8]*x54 + A[8]*x56 + A[9]*x14 + A[9]*x22 + A[9]*x24 - A[9]*x26 - x5*x61 - x6*x61 - x61*x7 + x61*x8;
    res[8]=A[10]*x42 - A[10]*x49 + A[10]*x60 + A[10]*x63 + A[6]*x44 - A[6]*x50 - A[6]*x52 - A[6]*x57 - A[7]*x43 - A[7]*x51 - A[7]*x54 + A[7]*x56 + A[8]*x0 + A[8]*x1 + A[8]*x10 + A[8]*x11 - A[8]*x12 - A[8]*x2 - A[8]*x3 - A[8]*x9 + A[9]*x39 - A[9]*x46 - A[9]*x48 - A[9]*x55 + x5*x65 - x6*x65 + x65*x7 + x65*x8;
    res[9]=-A[10]*x41 - A[10]*x58 - A[10]*x62 + A[10]*x64 - A[5]*x44 - A[5]*x50 - A[5]*x52 + A[5]*x57 + A[7]*x14 - A[7]*x22 - A[7]*x24 - A[7]*x26 + A[8]*x39 + A[8]*x46 - A[8]*x48 + A[8]*x55 + A[9]*x0 - A[9]*x1 - A[9]*x10 + A[9]*x11 - A[9]*x12 + A[9]*x2 - A[9]*x3 + A[9]*x9 + x5*x66 + x6*x66 - x66*x7 + x66*x8;
    res[10]=A[10]*x0 - A[10]*x1 + A[10]*x10 - A[10]*x11 - A[10]*x12 - A[10]*x2 + A[10]*x3 + A[10]*x9 + A[5]*x43 - A[5]*x51 - A[5]*x54 - A[5]*x56 - A[6]*x14 - A[6]*x22 - A[6]*x24 + A[6]*x26 + A[8]*x42 - A[8]*x49 - A[8]*x60 - A[8]*x63 - A[9]*x41 + A[9]*x58 + A[9]*x62 + A[9]*x64 + x5*x67 + x6*x67 + x67*x7 - x67*x8;
    res[11]=A[11]*x0 + A[11]*x1 + A[11]*x10 + A[11]*x11 + A[11]*x12 + A[11]*x2 + A[11]*x3 + A[11]*x9 + A[12]*x14 + A[12]*x22 - A[12]*x24 + A[12]*x26 + A[13]*x43 - A[13]*x51 + A[13]*x54 + A[13]*x56 + A[14]*x44 + A[14]*x50 - A[14]*x52 + A[14]*x57;
    res[12]=A[11]*x14 - A[11]*x22 + A[11]*x24 + A[11]*x26 + A[12]*x0 + A[12]*x1 - A[12]*x10 - A[12]*x11 + A[12]*x12 - A[12]*x2 - A[12]*x3 + A[12]*x9 + A[13]*x39 - A[13]*x46 + A[13]*x48 + A[13]*x55 + A[14]*x42 + A[14]*x49 + A[14]*x60 - A[14]*x63;
    res[13]=A[11]*x43 + A[11]*x51 - A[11]*x54 + A[11]*x56 + A[12]*x39 + A[12]*x46 + A[12]*x48 - A[12]*x55 + A[13]*x0 - A[13]*x1 + A[13]*x10 - A[13]*x11 + A[13]*x12 + A[13]*x2 - A[13]*x3 - A[13]*x9 + A[14]*x41 - A[14]*x58 + A[14]*x62 + A[14]*x64;
    res[14]=A[11]*x44 - A[11]*x50 + A[11]*x52 + A[11]*x57 + A[12]*x42 + A[12]*x49 - A[12]*x60 + A[12]*x63 + A[13]*x41 + A[13]*x58 - A[13]*x62 + A[13]*x64 + A[14]*x0 - A[14]*x1 - A[14]*x10 + A[14]*x11 + A[14]*x12 - A[14]*x2 + A[14]*x3 - A[14]*x9;
    res[15]=A[15]*x0 + A[15]*x1 - A[15]*x10 - A[15]*x11 - A[15]*x12 + A[15]*x2 + A[15]*x3 - A[15]*x9 + x5*x68 - x6*x68 - x68*x7 - x68*x8;
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
    T x10 = 2.0*A[2];
    T x11 = (*this)[0]*(*this)[3];
    T x12 = 2.0*A[1];
    T x13 = (*this)[0]*(*this)[7];
    T x14 = 2.0*(*this)[1];
    T x15 = A[1]*x14;
    T x16 = (*this)[5]*A[2];
    T x17 = (*this)[1]*x8;
    T x18 = (*this)[2]*(*this)[4];
    T x19 = (*this)[5]*x12;
    T x20 = (*this)[3]*x8;
    T x21 = (*this)[3]*(*this)[6];
    T x22 = (*this)[7]*x8;
    T x23 = (*this)[7]*x10;
    T x24 = A[2]*x14;
    T x25 = 2.0*x16;
    T x26 = (*this)[6]*x12;
    T x27 = (*this)[2]*x12;
    T x28 = (*this)[6]*x10;
    T x29 = (*this)[3]*x12;
    T x30 = (*this)[5]*x8;
    T x31 = (*this)[4]*x8;
    T x32 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 + A[0]*x2 + A[0]*x3 - A[0]*x4 - A[0]*x5 - A[0]*x6 - A[0]*x7;
    res[1]=(*this)[0]*x9 - (*this)[2]*x19 + (*this)[4]*x15 + (*this)[4]*x20 - (*this)[5]*x22 + (*this)[6]*x17 + (*this)[6]*x23 - x10*x11 + x10*x18 - x12*x13 - x12*x21 + x14*x16;
    res[2]=-(*this)[0]*x17 - (*this)[0]*x23 + (*this)[2]*x25 + (*this)[4]*x22 - (*this)[4]*x24 + (*this)[5]*x15 + (*this)[5]*x20 + (*this)[6]*x9 - (*this)[7]*x26 - x10*x21 + x11*x12 + x12*x18;
    res[3]=(*this)[0]*x24 - (*this)[0]*x27 + (*this)[2]*x28 + (*this)[3]*x25 - (*this)[4]*x17 - (*this)[4]*x23 + (*this)[4]*x29 - (*this)[5]*x9 + (*this)[6]*x15 + (*this)[6]*x20 + (*this)[7]*x19 - x13*x8;
    res[4]=-(*this)[0]*x28 + (*this)[0]*x30 - (*this)[2]*x24 - (*this)[3]*x17 - (*this)[3]*x23 + (*this)[4]*x25 + (*this)[6]*x31 + (*this)[7]*x9 + A[1]*x0 + A[1]*x1 - A[1]*x2 - A[1]*x3 - A[1]*x4 + A[1]*x5 + A[1]*x6 - A[1]*x7;
    res[5]=(*this)[0]*x26 - (*this)[0]*x31 - (*this)[2]*x15 - (*this)[3]*x9 + (*this)[4]*x19 + (*this)[6]*x30 - (*this)[7]*x17 + (*this)[7]*x29 + A[2]*x0 - A[2]*x1 + A[2]*x2 - A[2]*x3 + A[2]*x4 - A[2]*x5 + A[2]*x6 - A[2]*x7;
    res[6]=(*this)[0]*(*this)[4]*x10 - (*this)[0]*x19 - (*this)[2]*(*this)[3]*x10 - (*this)[3]*x15 + (*this)[4]*x26 + (*this)[6]*x25 + (*this)[7]*x24 - (*this)[7]*x27 + A[3]*x0 - A[3]*x1 - A[3]*x2 + A[3]*x3 + A[3]*x4 + A[3]*x5 - A[3]*x6 - A[3]*x7;
    res[7]=-(*this)[2]*(*this)[5]*x32 - (*this)[4]*A[0]*x14 + x13*x32 - x21*x32;
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
    T x4 = 2.0*A[7];
    T x5 = (*this)[0]*(*this)[7];
    T x6 = (*this)[1]*(*this)[4];
    T x7 = (*this)[2]*(*this)[5];
    T x8 = (*this)[3]*(*this)[6];
    T x9 = std::pow((*this)[1], 2);
    T x10 = std::pow((*this)[2], 2);
    T x11 = std::pow((*this)[3], 2);
    T x12 = std::pow((*this)[7], 2);
    T x13 = 2.0*(*this)[0];
    T x14 = (*this)[2]*A[6];
    T x15 = (*this)[3]*x13;
    T x16 = (*this)[5]*A[3];
    T x17 = (*this)[6]*x13;
    T x18 = 2.0*A[4];
    T x19 = 2.0*(*this)[1];
    T x20 = (*this)[2]*x19;
    T x21 = (*this)[3]*A[3];
    T x22 = (*this)[5]*A[5];
    T x23 = (*this)[6]*A[6];
    T x24 = 2.0*(*this)[2];
    T x25 = (*this)[4]*A[5];
    T x26 = (*this)[7]*A[3];
    T x27 = 2.0*(*this)[4];
    T x28 = (*this)[3]*A[6];
    T x29 = 2.0*(*this)[7];
    T x30 = (*this)[3]*x29;
    T x31 = (*this)[5]*x27;
    T x32 = (*this)[6]*x27;
    T x33 = (*this)[5]*A[6];
    T x34 = (*this)[6]*A[5];
    T x35 = (*this)[1]*x13;
    T x36 = (*this)[4]*x13;
    T x37 = 2.0*A[5];
    T x38 = (*this)[1]*x18;
    T x39 = (*this)[4]*x18;
    T x40 = 2.0*(*this)[6];
    T x41 = 2.0*(*this)[5];
    T x42 = (*this)[7]*A[6];
    T x43 = (*this)[7]*x18;
    T x44 = (*this)[2]*x13;
    T x45 = (*this)[5]*A[1];
    T x46 = 2.0*A[6];
    T x47 = A[1]*x19;
    T x48 = A[2]*x19;
    T x49 = A[2]*x24;
    T x50 = A[1]*x24;
    T x51 = 2.0*(*this)[3];
    T x52 = 2.0*A[1];
    T x53 = (*this)[6]*A[3];
    T x54 = (*this)[6]*x29;
    T x55 = 2.0*A[2];
    T x56 = 2.0*A[3];
    T x57 = 2.0*A[0];
    res[0]=A[0]*x0 + A[0]*x1 - A[0]*x10 - A[0]*x11 - A[0]*x12 + A[0]*x2 + A[0]*x3 - A[0]*x9 - x4*x5 + x4*x6 + x4*x7 + x4*x8;
    res[1]=A[1]*x0 + A[1]*x1 + A[1]*x10 + A[1]*x11 - A[1]*x12 - A[1]*x2 - A[1]*x3 - A[1]*x9 - A[2]*x17 - A[2]*x20 - A[2]*x30 + A[2]*x31 + A[3]*x32 - A[5]*x15 + x13*x14 + x13*x16 - x18*x5 + x18*x6 - x18*x7 - x18*x8 - x19*x21 + x19*x22 + x19*x23 + x24*x25 + x24*x26 + x27*x28 - x29*x33 + x29*x34;
    res[2]=(*this)[2]*x39 + (*this)[5]*x38 - (*this)[6]*x43 + A[1]*x17 - A[1]*x20 + A[1]*x30 + A[1]*x31 + A[2]*x0 - A[2]*x1 - A[2]*x10 + A[2]*x11 - A[2]*x12 + A[2]*x2 - A[2]*x3 + A[2]*x9 - A[3]*x36 + A[4]*x15 - A[6]*x35 + x14*x40 + x16*x40 - x19*x26 - x21*x24 + x27*x42 + x28*x41 - x37*x5 - x37*x6 + x37*x7 - x37*x8;
    res[3]=(*this)[3]*x39 - (*this)[3]*x47 - (*this)[3]*x49 + (*this)[5]*A[2]*x40 + (*this)[5]*x43 + (*this)[6]*x38 + (*this)[7]*x48 - (*this)[7]*x50 + A[1]*x32 + A[2]*x36 + A[3]*x0 - A[3]*x1 + A[3]*x10 - A[3]*x11 - A[3]*x12 - A[3]*x2 + A[3]*x3 + A[3]*x9 - A[4]*x44 + A[5]*x35 - x13*x45 + x22*x51 + x24*x34 - x25*x29 - x46*x5 - x46*x6 - x46*x7 + x46*x8;
    res[4]=-(*this)[4]*x49 - (*this)[5]*x48 + A[2]*x15 - A[2]*x54 - A[3]*x44 + A[4]*x0 + A[4]*x1 + A[4]*x10 + A[4]*x11 - A[4]*x12 - A[4]*x2 - A[4]*x3 - A[4]*x9 - A[5]*x20 - A[5]*x30 + x13*x33 - x13*x34 + x14*x29 + x16*x29 - x19*x28 - x19*x53 - x21*x27 + x22*x27 + x23*x27 + x5*x52 - x52*x6 + x52*x7 + x52*x8;
    res[5]=-(*this)[2]*x38 + (*this)[3]*x43 - (*this)[4]*x50 + (*this)[5]*x39 - A[1]*x15 + A[1]*x54 + A[3]*x35 + A[4]*x17 + A[5]*x0 - A[5]*x1 - A[5]*x10 + A[5]*x11 - A[5]*x12 + A[5]*x2 - A[5]*x3 + A[5]*x9 - A[6]*x36 - x14*x51 - x16*x51 - x19*x42 - x19*x45 + x23*x41 - x24*x53 - x26*x27 + x5*x55 + x55*x6 - x55*x7 + x55*x8;
    res[6]=-(*this)[2]*x43 - (*this)[3]*A[1]*x27 - (*this)[3]*A[2]*x41 - (*this)[3]*A[5]*x24 - (*this)[3]*x38 - (*this)[5]*A[4]*x13 + (*this)[6]*x39 - (*this)[6]*x47 - (*this)[6]*x49 + (*this)[7]*A[2]*x27 + (*this)[7]*A[5]*x19 + A[1]*x44 - A[2]*x35 + A[6]*x0 - A[6]*x1 + A[6]*x10 - A[6]*x11 - A[6]*x12 - A[6]*x2 + A[6]*x3 + A[6]*x9 + x13*x25 + x22*x40 - x29*x45 + x5*x56 + x56*x6 + x56*x7 - x56*x8;
    res[7]=A[7]*x0 + A[7]*x1 - A[7]*x10 - A[7]*x11 - A[7]*x12 + A[7]*x2 + A[7]*x3 - A[7]*x9 + x5*x57 - x57*x6 - x57*x7 - x57*x8;
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
