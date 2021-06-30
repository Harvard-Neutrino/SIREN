#ifndef LI_Quaternion_H
#define LI_Quaternion_H

#include <sstream>
#include "earthmodel-service/Vector3D.h"

namespace earthmodel {

class Vector3D;

class Quaternion
{
public:
    // constructors
    Quaternion();
    Quaternion(const double x, const double y, const double z, const double w);
    Quaternion(const Quaternion& quaternion);
    Quaternion(Quaternion&& other);
    ~Quaternion();

    void SetEulerAngles(double alpha, double beta, double gamma);
    void SetPosition(Vector3D const & vec);

    void GetMatrix(Matrix3D &) const;
    Matrix3D GetMatrix() const;
    void SetMatrix();

    Quaternion & invert();
    Quaternion & normalize();

    double DotProduct(Quaternion const &);

    static Quaternion lerp(Quaternion const &, Quaternion const &);
    static Quaternion slerp(Quaternion const &, Quaternion const &);

    void SetAxisAngle(Vector3D const &, double);
    void GetAxisAngle(Vector3D &, double &);
    std::tuple<Vector3D, double> GetAxisAngle();

    void GetAnglesEulerZXZ(double & alpha, double & beta, double & gamma);


    //-------------------------------------//
    // operator functions and swap
    Quaternion& operator=(Quaternion const & quaternion);
    Quaternion& operator=(Quaternion const && quaternion);
    Quaternion& operator=(Quaternion && quaternion);
    bool operator==(const Quaternion& quaternion) const;
    bool operator!=(const Quaternion& quaternion) const;
    void swap(Quaternion& quaternion);
    friend std::ostream& operator<<(std::ostream& os, Quaternion const& quaternion);

    Quaternion operator*(Quaternion & const) const;
    Quaternion & operator*=(Quaternion & const);
    Quaternion operator*(double) const;
    Quaternion & operator*=(double);
private:
    double x_;
    double y_;
    double z_;
    double w_;
};

} // namespace earthmodel

#endif // LI_Quaternion_H

