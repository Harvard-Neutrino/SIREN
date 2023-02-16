#pragma once
#ifndef LI_Matrix3D_H
#define LI_Matrix3D_H

#include <sstream>
#include <initializer_list>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

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

    double const & operator[](std::initializer_list<unsigned int>) const;
    double & operator[](std::initializer_list<unsigned int>);

    //-------------------------------------//
    // getter
    //----------------------------------------------//
    double GetXX() const {return xx_;};
    double GetXY() const {return xy_;};
    double GetXZ() const {return xz_;};
    double GetYX() const {return yx_;};
    double GetYY() const {return yy_;};
    double GetYZ() const {return yz_;};
    double GetZX() const {return zx_;};
    double GetZY() const {return zy_;};
    double GetZZ() const {return zz_;};

    void SetXX(double xx) {xx_ = xx;};
    void SetXY(double xy) {xy_ = xy;};
    void SetXZ(double xz) {xz_ = xz;};
    void SetYX(double yx) {yx_ = yx;};
    void SetYY(double yy) {yy_ = yy;};
    void SetYZ(double yz) {yz_ = yz;};
    void SetZX(double zx) {zx_ = zx;};
    void SetZY(double zy) {zy_ = zy;};
    void SetZZ(double zz) {zz_ = zz;};

    //-------------------------------------//
    // serialization
    //----------------------------------------------//
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("XX", xx_));
            archive(cereal::make_nvp("XY", xy_));
            archive(cereal::make_nvp("XZ", xz_));
            archive(cereal::make_nvp("YX", yx_));
            archive(cereal::make_nvp("YY", yy_));
            archive(cereal::make_nvp("YZ", yz_));
            archive(cereal::make_nvp("ZX", zx_));
            archive(cereal::make_nvp("ZY", zy_));
            archive(cereal::make_nvp("ZZ", zz_));
        } else {
            throw std::runtime_error("Matrix3D only supports version <= 0!");
        }
    }

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

CEREAL_CLASS_VERSION(earthmodel::Matrix3D, 0);

#endif // LI_Matrix3D_H

