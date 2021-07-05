#ifndef LI_Conversions_H
#define LI_Conversions_H

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
        res[i] = cb * (cs + ss);
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

} // namespace earthmodel

#endif // LI_Conversions_H

