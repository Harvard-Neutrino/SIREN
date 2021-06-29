#ifndef LI_Matrix3D_H
#define LI_Matrix3D_H

#include <sstream>
#include "earthmodel-service/Vector3D.h"

namespace earthmodel {

class Vector3D;

class Matrix3D
{
public:
    // constructors
    Matrix3D();
    Matrix3D(
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
    Matrix3D(const Matrix3D& matrix_3d);
    Matrix3D(Matrix3D&& other);
    ~Matrix3D();

    //-------------------------------------//
    // operator functions and swap
    Matrix3D& operator=(Matrix3D const & matrix_3d);
    Matrix3D& operator=(Matrix3D const && matrix_3d);
    Matrix3D& operator=(Matrix3D && matrix_3d);
    bool operator==(const Matrix3D& matrix_3d) const;
    bool operator!=(const Matrix3D& matrix_3d) const;
    void swap(Matrix3D& matrix_3d);
    friend std::ostream& operator<<(std::ostream& os, Matrix3D const& matrix_3d);

    //-------------------------------------//
    // basic arithmetic
    friend Matrix3D operator+(const Matrix3D& mat1, const Matrix3D& mat2);
    friend Matrix3D& operator+=(Matrix3D& mat1, const Matrix3D& mat2);
    friend Matrix3D operator-(const Matrix3D& mat1, const Matrix3D& mat2);
    friend Matrix3D operator*(const double factor1, const Matrix3D& mat1);
    friend Matrix3D operator*(const Matrix3D& mat1, const double factor1);
    friend Vector3D operator*(const Matrix3D& mat1, const Vector3D& vec1);
    friend Matrix3D operator/(const Matrix3D& mat1, const double factor1);
    friend Matrix3D& operator*=(Matrix3D& mat1, const double factor1);
    friend Matrix3D& operator/=(Matrix3D& mat1, const double factor1);
    friend Matrix3D operator*(const Matrix3D& mat1, const Matrix3D& mat2);
    friend Matrix3D scalar_product(const Matrix3D& mat1, const Matrix3D& mat2);
    friend Matrix3D matrix_product(const Matrix3D& mat1, const Matrix3D& mat2);
    Matrix3D operator-() const;

    //-------------------------------------//
    // getter
    //----------------------------------------------//
private:
    double xx_ = 0.0;
    double xy_ = 0.0;
    double xz_ = 0.0;
    double yx_ = 0.0;
    double yy_ = 0.0;
    double yz_ = 0.0;
    double zx_ = 0.0;
    double zy_ = 0.0;
    double zz_ = 0.0;
};

} // namespace earthmodel

#endif // LI_Matrix3D_H

