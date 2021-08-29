#ifndef LI_Quaternion_H
#define LI_Quaternion_H

#include <sstream>
#include "earthmodel-service/Matrix3D.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EulerAngles.h"

namespace earthmodel {

class Matrix3D;
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

    void SetPosition(Vector3D const & vec);

    void GetMatrix(Matrix3D &) const;
    Matrix3D GetMatrix() const;
    void SetMatrix(Matrix3D const &);

    Quaternion & invert();
    Quaternion & conjugate();
    Quaternion & normalize();
    Quaternion inverted() const;
    Quaternion conjugated() const;
    Quaternion normalized() const;
    double magnitude() const;
    double magnitudesq() const;

    double DotProduct(Quaternion const &) const;

    static Quaternion lerp(Quaternion const &, Quaternion const &, double);
    static Quaternion slerp(Quaternion const &, Quaternion const &, double);

    void SetAxisAngle(Vector3D const &, double);
    void GetAxisAngle(Vector3D &, double &) const;
    std::tuple<Vector3D, double> GetAxisAngle() const;

    void GetEulerAngles(EulerAngles & euler, EulerOrder order) const;
    void SetEulerAngles(EulerAngles const & euler);
    void GetEulerAnglesZXZr(double & alpha, double & beta, double & gamma) const;
    void SetEulerAnglesZXZr(double alpha, double beta, double gamma);
    void GetEulerAnglesXYZs(double & alpha, double & beta, double & gamma) const;
    void SetEulerAnglesXYZs(double alpha, double beta, double gamma);

    //-------------------------------------//
    // operator functions and swap
    Quaternion& operator=(Quaternion const & quaternion);
    Quaternion& operator=(Quaternion const && quaternion);
    Quaternion& operator=(Quaternion && quaternion);
    bool operator==(const Quaternion& quaternion) const;
    bool operator!=(const Quaternion& quaternion) const;
    void swap(Quaternion& quaternion);
    friend std::ostream& operator<<(std::ostream& os, Quaternion const& quaternion);

    Quaternion operator*(Quaternion const &) const;
    Quaternion & operator*=(Quaternion const &);
    Quaternion operator*(double) const;
    Quaternion & operator*=(double);
    Quaternion operator+(Quaternion const &) const;
    Quaternion & operator+=(Quaternion const &);
    Quaternion operator+(double) const;
    Quaternion & operator+=(double);

    Quaternion operator~() const;
    Quaternion operator!() const;

    Quaternion rotate(Quaternion const & p, bool inv) const;
    Vector3D rotate(Vector3D const & p, bool inv) const;

    double GetX() const {return x_;}
    double GetY() const {return y_;}
    double GetZ() const {return z_;}
    double GetW() const {return w_;}
    void SetX(double x) {x_ = x;}
    void SetY(double y) {y_ = y;}
    void SetZ(double z) {z_ = z;}
    void SetW(double w) {w_ = w;}
private:
    double x_;
    double y_;
    double z_;
    double w_;
};

} // namespace earthmodel

#endif // LI_Quaternion_H

