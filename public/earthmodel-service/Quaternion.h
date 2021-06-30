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
    Matrix3D GetMatrix() const;
    void SetMatrix();


    //-------------------------------------//
    // operator functions and swap
    Quaternion& operator=(Quaternion const & quaternion);
    Quaternion& operator=(Quaternion const && quaternion);
    Quaternion& operator=(Quaternion && quaternion);
    bool operator==(const Quaternion& quaternion) const;
    bool operator!=(const Quaternion& quaternion) const;
    void swap(Quaternion& quaternion);
    friend std::ostream& operator<<(std::ostream& os, Quaternion const& quaternion);

    friend Quaternion operator*(Quaternion & const)
private:
    double x_;
    double y_;
    double z_;
    double w_;
};

} // namespace earthmodel

#endif // LI_Quaternion_H

