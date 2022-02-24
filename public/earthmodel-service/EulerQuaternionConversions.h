#pragma once
#ifndef LI_EulerQuaternionConversions_H
#define LI_EulerQuaternionConversions_H

#include <limits>

#include "earthmodel-service/Quaternion.h"
#include "earthmodel-service/EulerAngles.h"

namespace earthmodel {

class Quaternion;
class EulerAngles;

inline
Quaternion QFromZXZr(double alpha, double beta, double gamma)
{
    alpha *= 0.5;
    beta *= 0.5;
    gamma *= 0.5;

    double cb = cos(beta);
    double sb = sin(beta);

    double res[4];
    return Quaternion(
            sb * cos(alpha - gamma),
            sb * sin(alpha - gamma),
            cb * sin(alpha + gamma),
            cb * cos(alpha + gamma)
        );
}

inline
EulerAngles ZXZrFromQ(Quaternion const & quaternion) {
    double Nq = quaternion.DotProduct(quaternion);
    double s = (Nq > 0) ? (2.0 / Nq) : 0;
    double x =  quaternion.GetX();
    double y =  quaternion.GetY();
    double z =  quaternion.GetZ();
    double w =  quaternion.GetW();
    double ww = w * w * s;
    double xs = x * s,   ys = y * s,   zs = z * s;
    double wx = w * xs,  wy = w * ys,  wz = w * zs;
    double xx = x * xs,  xy = x * ys,  xz = x * zs;
    double yy = y * ys,  yz = y * zs,  zz = z * zs;
    double check = sqrt((xx + yy) * (ww + zz));

    double alpha, beta, gamma;

    if(check > 16 * std::numeric_limits<double>::epsilon()) {
        alpha = atan2(wy + xz, wx - yz);
        beta = atan2(check, 1 - (xx + yy));
        gamma = atan2(xz - wy, wx + yz);
    } else {
        alpha = 0;
        beta = atan2(check, 1 - (xx + yy));
        gamma = atan2(wz - xy, 1 - (yy + zz));
    }
    return EulerAngles(EulerOrder::ZXZr, alpha, beta, gamma);
}

inline
Quaternion QFromXYZs(double alpha, double beta, double gamma)
{
    alpha *= 0.5;
    beta *= 0.5;
    gamma *= 0.5;

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

    double res[4];
    return Quaternion(
            cb * sc - sb * cs,
            sb * cc - cb * ss,
            cb * cs - sb * sc,
            cb * cc + sb * ss
        );
}

inline
EulerAngles XYZsFromQ(Quaternion const & quaternion) {
    double Nq = quaternion.DotProduct(quaternion);
    double s = (Nq > 0) ? (2.0 / Nq) : 0;
    double x =  quaternion.GetX();
    double y =  quaternion.GetY();
    double z =  quaternion.GetZ();
    double w =  quaternion.GetW();
    double ww = w * w * s;
    double xs = x * s,   ys = y * s,   zs = z * s;
    double wx = w * xs,  wy = w * ys,  wz = w * zs;
    double xx = x * xs,  xy = x * ys,  xz = x * zs;
    double yy = y * ys,  yz = y * zs,  zz = z * zs;
    double check = sqrt(1 - (wy - xz) * (wy - xz));

    double alpha, beta, gamma;

    if(check > 16 * std::numeric_limits<double>::epsilon()) {
        alpha = atan2(wx + yz, 1 - (xx + yy));
        beta = atan2(wy - xz, check);
        gamma = atan2(xy + wz, 1 - (yy + zz));
    } else {
        alpha = atan2(wx - yz, 1 - (xx + zz));
        beta = atan2(wy - xz, check);
        gamma = 0;
    }
    return EulerAngles(EulerOrder::XYZs, alpha, beta, gamma);
}

} // namespace earthmodel

#endif // LI_EulerQuaternionConversions_H

