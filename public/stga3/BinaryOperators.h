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

//----------------------------------
// (Boost, R130MV) binary operations
//----------------------------------

template<typename T>
inline Boost::operator R130MV<T>() const {
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

//---------------------------------
// (Boost, Rotor) binary operations
//---------------------------------

template<typename T>
inline Boost::operator Rotor<T>() const {
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

//----------------------------------
// (R130B0, Boost) binary operations
//----------------------------------

template<typename T>
inline R130B0::operator Boost<T>() const {
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

//-----------------------------------
// (R130B0, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B0::operator R130MV<T>() const {
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

//-------------------------------------
// (R130B0, Rotation) binary operations
//-------------------------------------

template<typename T>
inline R130B0::operator Rotation<T>() const {
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

//----------------------------------
// (R130B0, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130B0::operator Rotor<T>() const {
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

//-----------------------------------
// (R130B1, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B1::operator R130MV<T>() const {
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

//--------------------------------------
// (R130B1Sm1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sm1::operator R130B1<T>() const {
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

//--------------------------------------
// (R130B1Sm1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sm1::operator R130MV<T>() const {
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

//--------------------------------------
// (R130B1Sp1, R130B1) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sp1::operator R130B1<T>() const {
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

//--------------------------------------
// (R130B1Sp1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B1Sp1::operator R130MV<T>() const {
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

//-----------------------------------
// (R130B2, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B2::operator R130MV<T>() const {
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

//----------------------------------
// (R130B2, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130B2::operator Rotor<T>() const {
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

//--------------------------------------
// (R130B2Sm1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sm1::operator R130B2<T>() const {
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

//--------------------------------------
// (R130B2Sm1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sm1::operator R130MV<T>() const {
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

//----------------------------------------
// (R130B2Sm1, Rotation) binary operations
//----------------------------------------

template<typename T>
inline R130B2Sm1::operator Rotation<T>() const {
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

//-------------------------------------
// (R130B2Sm1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130B2Sm1::operator Rotor<T>() const {
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

//-------------------------------------
// (R130B2Sp1, Boost) binary operations
//-------------------------------------

template<typename T>
inline R130B2Sp1::operator Boost<T>() const {
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

//--------------------------------------
// (R130B2Sp1, R130B2) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sp1::operator R130B2<T>() const {
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

//--------------------------------------
// (R130B2Sp1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B2Sp1::operator R130MV<T>() const {
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

//-------------------------------------
// (R130B2Sp1, Rotor) binary operations
//-------------------------------------

template<typename T>
inline R130B2Sp1::operator Rotor<T>() const {
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

//-----------------------------------
// (R130B3, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B3::operator R130MV<T>() const {
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

//--------------------------------------
// (R130B3Sm1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sm1::operator R130B3<T>() const {
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

//--------------------------------------
// (R130B3Sm1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sm1::operator R130MV<T>() const {
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

//--------------------------------------
// (R130B3Sp1, R130B3) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sp1::operator R130B3<T>() const {
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

//--------------------------------------
// (R130B3Sp1, R130MV) binary operations
//--------------------------------------

template<typename T>
inline R130B3Sp1::operator R130MV<T>() const {
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

//-----------------------------------
// (R130B4, R130MV) binary operations
//-----------------------------------

template<typename T>
inline R130B4::operator R130MV<T>() const {
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

//----------------------------------
// (R130B4, Rotor) binary operations
//----------------------------------

template<typename T>
inline R130B4::operator Rotor<T>() const {
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

//-------------------------------------
// (Rotation, R130MV) binary operations
//-------------------------------------

template<typename T>
inline Rotation::operator R130MV<T>() const {
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

//------------------------------------
// (Rotation, Rotor) binary operations
//------------------------------------

template<typename T>
inline Rotation::operator Rotor<T>() const {
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

//----------------------------------
// (Rotor, R130MV) binary operations
//----------------------------------

template<typename T>
inline Rotor::operator R130MV<T>() const {
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


} // namespace stga3

#endif // LI_STGA3_BinaryOperators_H
