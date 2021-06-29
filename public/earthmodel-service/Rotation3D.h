#ifndef LI_Rotation3D_H
#define LI_Rotation3D_H

#include <sstream>
#include "earthmodel-service/Matrix3D.h"
#include "earthmodel-service/Vector3D.h"

namespace earthmodel {

class Rotation3D
{
public:
    // constructors
    Rotation3D();
    Rotation3D(
            const double xx,
            const double xy,
            const double xz,
            const double yx,
            const double yy,
            const double yz,
            const double zx,
            const double zy,
            const double zz
            );
    Rotation3D(const Rotation3D& matrix_3d);
    Rotation3D(Rotation3D&& other);
    ~Rotation3D();

    //-------------------------------------//
    // operator functions and swap
    Rotation3D& operator=(Rotation3D const & matrix_3d);
    Rotation3D& operator=(Rotation3D const && matrix_3d);
    Rotation3D& operator=(Rotation3D && matrix_3d);
    bool operator==(const Rotation3D& matrix_3d) const;
    bool operator!=(const Rotation3D& matrix_3d) const;
    void swap(Rotation3D& matrix_3d);
    friend std::ostream& operator<<(std::ostream& os, Rotation3D const& matrix_3d);

    //-------------------------------------//
    // basic arithmetic
    friend Rotation3D operator+(const Rotation3D& mat1, const Rotation3D& mat2);
    friend Rotation3D& operator+=(Rotation3D& mat1, const Rotation3D& mat2);
    friend Rotation3D operator-(const Rotation3D& mat1, const Rotation3D& mat2);
    friend Rotation3D operator*(const double factor1, const Rotation3D& mat1);
    friend Rotation3D operator*(const Rotation3D& mat1, const double factor1);
    friend Vector3D operator*(const Rotation3D& mat1, const Vector3D& vec1);
    friend Rotation3D operator/(const Rotation3D& mat1, const double factor1);
    friend Rotation3D& operator*=(Rotation3D& mat1, const double factor1);
    friend Rotation3D& operator/=(Rotation3D& mat1, const double factor1);
    friend Rotation3D operator*(const Rotation3D& mat1, const Rotation3D& mat2);
    friend Rotation3D scalar_product(const Rotation3D& mat1, const Rotation3D& mat2);
    friend Rotation3D matrix_product(const Rotation3D& mat1, const Rotation3D& mat2);
    Rotation3D operator-() const;

    //-------------------------------------//
    // getter
    //----------------------------------------------//
private:
    Matrix3D rotation_matrix_;
    Vector3D euler_angles_;
    Vector3D yaw_pitch_roll_;
};

} // namespace earthmodel

#endif // LI_Rotation3D_H

