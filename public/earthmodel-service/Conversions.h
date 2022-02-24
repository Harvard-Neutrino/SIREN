#pragma once
#ifndef LI_Conversions_H
#define LI_Conversions_H

#include <limits>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Matrix3D.h"
#include "earthmodel-service/Quaternion.h"
#include "earthmodel-service/EulerAngles.h"

namespace earthmodel {

class Vector3D;
class Matrix3D;
class Quaternion;
class EulerAngles;

inline
Quaternion QuaternionFromEulerAngles(EulerAngles const & euler)
{
    double h_alpha = euler.GetAlpha() * 0.5;
    double h_beta = euler.GetBeta() * 0.5;
    double h_gamma = euler.GetGamma() * 0.5;
    EulerOrder order = euler.GetOrder();
    EulerAxis i = GetEulerAxisI(order);
    EulerAxis j = GetEulerAxisJ(order);
    EulerAxis k = GetEulerAxisK(order);
    EulerAxis h = GetEulerAxisH(order);
    EulerParity n = GetEulerParity(order);
    EulerRepetition s = GetEulerRepetition(order);
    EulerFrame f = GetEulerFrame(order);

    if(f == EulerFrame::Rotating) {
        std::swap(h_alpha, h_gamma);
    }

    if(n == EulerParity::Odd) {
        h_beta = -h_beta;
    }

    double ca = cos(h_alpha);
    double cb = cos(h_beta);
    double cg = cos(h_gamma);
    double sa = sin(h_alpha);
    double sb = sin(h_beta);
    double sg = sin(h_gamma);

    double cc = ca*cg;
    double cs = ca*sg;
    double sc = sa*cg;
    double ss = sa*sg;

    double res[4];

    if(s == EulerRepetition::Yes) {
        res[i] = cb * (cs + sc);
        res[j] = sb * (cc + ss);
        res[k] = sb * (cs - sc);
        res[3] = cb * (cc - ss);
    } else {
        res[i] = cb * sc - sb * cs;
        res[j] = cb * ss + sb * cc;
        res[k] = cb * cs - sb * sc;
        res[3] = cb * cc + sb * ss;
    }

    if(n == EulerParity::Odd) {
        res[j] = -res[j];
    }

    return Quaternion(
        res[(unsigned int)EulerAxis::X],
        res[(unsigned int)EulerAxis::Y],
        res[(unsigned int)EulerAxis::Z],
        res[3]
    );
}

inline
EulerAngles EulerAnglesFromMatrix3D(Matrix3D const & matrix, EulerOrder const & order)
{
    EulerAxis i = GetEulerAxisI(order);
    EulerAxis j = GetEulerAxisJ(order);
    EulerAxis k = GetEulerAxisK(order);
    EulerAxis h = GetEulerAxisH(order);
    EulerParity n = GetEulerParity(order);
    EulerRepetition s = GetEulerRepetition(order);
    EulerFrame f = GetEulerFrame(order);

    double alpha, beta, gamma;

    if(s == EulerRepetition::Yes) {
        double sy = sqrt(matrix[{i,j}] * matrix[{i,j}] + matrix[{i,k}] * matrix[{i,k}]);
        if(sy > 16 * std::numeric_limits<double>::epsilon()) {
            alpha = atan2(matrix[{i,j}], matrix[{i,k}]);
            beta = atan2(sy, matrix[{i,i}]);
            gamma = atan2(matrix[{j,i}], -matrix[{k,i}]);
        } else {
            alpha = atan2(-matrix[{j,k}], matrix[{j,j}]);
            beta = atan2(sy, matrix[{i,i}]);
            gamma = 0;
        }
    } else {
        double cy = sqrt(matrix[{i,i}] * matrix[{i,i}] + matrix[{j,i}] * matrix[{j,i}]);
        if(cy > 16 * std::numeric_limits<double>::epsilon()) {
            alpha = atan2(matrix[{k,j}], matrix[{k,k}]);
            beta = atan2(-matrix[{k,i}], cy);
            gamma = atan2(matrix[{j,i}], matrix[{i,i}]);
        } else {
            alpha = atan2(-matrix[{j,k}], matrix[{j,j}]);
            beta = atan2(-matrix[{k,i}], cy);
            gamma = 0;
        }
    }

    if(n == EulerParity::Odd) {
        alpha = -alpha;
        beta = -beta;
        gamma = -gamma;
    }

    if(f == EulerFrame::Rotating) {
        std::swap(alpha, gamma);
    }

    return EulerAngles(order, alpha, beta, gamma);
}

inline
Matrix3D Matrix3DFromEulerAngles(EulerAngles const & euler) {
    EulerOrder order = euler.GetOrder();
    EulerAxis i = GetEulerAxisI(order);
    EulerAxis j = GetEulerAxisJ(order);
    EulerAxis k = GetEulerAxisK(order);
    EulerAxis h = GetEulerAxisH(order);
    EulerParity n = GetEulerParity(order);
    EulerRepetition s = GetEulerRepetition(order);
    EulerFrame f = GetEulerFrame(order);

    double alpha = euler.GetAlpha();
    double beta = euler.GetBeta();
    double gamma = euler.GetGamma();

    if(f == EulerFrame::Rotating) {
        std::swap(alpha, beta);
    }

    if(n == EulerParity::Odd) {
        alpha = -alpha;
        beta = -beta;
        gamma = -gamma;
    }

    double ca = cos(alpha);
    double cb = cos(beta);
    double cg = cos(gamma);
    double sa = sin(alpha);
    double sb = sin(beta);
    double sg = sin(gamma);

    double cc = ca*cg;
    double cs = ca*sg;
    double sc = sa*cg;
    double ss = sa*sg;

    Matrix3D matrix;

    if(s == EulerRepetition::Yes) {
        matrix[{i,i}] = cb;
        matrix[{i,j}] = sb * sa;
        matrix[{i,k}] = sb * ca;
        matrix[{j,i}] = sb * sg;
        matrix[{j,j}] = -cb * ss + cc;
        matrix[{j,k}] = -cb * cs - sc;
        matrix[{k,i}] = -sb * cg;
        matrix[{k,j}] = cb * sc + cs;
        matrix[{k,k}] = cb * cc - ss;
    } else {
        matrix[{i,i}] = cb * cg;
        matrix[{i,j}] = sb * sc - cs;
        matrix[{i,k}] = sb * cc + ss;
        matrix[{j,i}] = cb * sg;
        matrix[{j,j}] = sb * ss + cc;
        matrix[{j,k}] = sb * cs - sc;
        matrix[{k,i}] = -sb;
        matrix[{k,j}] = cb * sa;
        matrix[{k,k}] = cb * ca;
    }

    return matrix;
}

inline
Matrix3D Matrix3DFromQuaternion(Quaternion const & quaternion) {
    Matrix3D matrix;
    double Nq = quaternion.DotProduct(quaternion);
    double s = (Nq > 0) ? (2.0 / Nq) : 0;
    double x =  quaternion.GetX();
    double y =  quaternion.GetY();
    double z =  quaternion.GetZ();
    double w =  quaternion.GetW();
    double xs = x * s,   ys = y * s,   zs = z * s;
    double wx = w * xs,  wy = w * ys,  wz = w * zs;
    double xx = x * xs,  xy = x * ys,  xz = x * zs;
    double yy = y * ys,  yz = y * zs,  zz = z * zs;

    matrix.SetXX(1 - (yy + zz));
    matrix.SetXY(xy - wz);
    matrix.SetXZ(xz + wy);
    matrix.SetYX(xy + wz);
    matrix.SetYY(1 - (xx + zz));
    matrix.SetYZ(yz - wx);
    matrix.SetZX(xz - wy);
    matrix.SetZY(yz + wx);
    matrix.SetZZ(1 - (xx + yy));
    return matrix;
}

inline
EulerAngles EulerAnglesFromQuaternion(Quaternion const & quaternion, EulerOrder const & order) {
    return EulerAnglesFromMatrix3D(quaternion.GetMatrix(), order);
}

} // namespace earthmodel

#endif // LI_Conversions_H

