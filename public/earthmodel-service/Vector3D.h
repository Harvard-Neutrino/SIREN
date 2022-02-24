
/******************************************************************************
 *                                                                            *
 * This file is adapted from a component of the simulation tool PROPOSAL.     *
 *                                                                            *
 * The PROPOSAL simulation tool is distributed under the terms of a modified  *
 * GNU Lesser General Public License version 3 (LGPLv3).                       *
 *                                                                            *
 * Modifcations to the LGPLv3 License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once
#ifndef LI_Vector3D_H
#define LI_Vector3D_H

#include <sstream>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include "earthmodel-service/Matrix3D.h"

#include <rk/geom3.hh>

namespace earthmodel {

class Matrix3D;

class Vector3D
{
public:
    // constructors
    Vector3D();
    Vector3D(const double x, const double y, const double z);
    Vector3D(const Vector3D& vector_3d);
    Vector3D(Vector3D&& other);
    Vector3D(std::array<double, 3> const & vec);
    Vector3D(geom3::UnitVector3 const & vec);
    Vector3D(geom3::Vector3 const & vec);
    Vector3D(geom3::Point3 const & vec);
    //Vector3D(const nlohmann::json&);
    ~Vector3D();

    //-------------------------------------//
    // operator functions and swap
    Vector3D& operator=(Vector3D const & vector_3d);
    Vector3D& operator=(Vector3D const && vector_3d);
    Vector3D& operator=(Vector3D && vector_3d);

    operator std::array<double, 3>() const;
    operator geom3::UnitVector3() const;
    operator geom3::Vector3() const;
    operator geom3::Point3() const;

    bool operator==(const Vector3D& vector_3d) const;
    bool operator!=(const Vector3D& vector_3d) const;
    bool operator<(const Vector3D& vector_3d) const;

    void swap(Vector3D& vector_3d);
    friend std::ostream& operator<<(std::ostream& os, Vector3D const& vector_3d);

    //-------------------------------------//
    // basic arithmetic
    Vector3D operator-() const;

    friend Vector3D operator+(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D& operator+=(Vector3D& vec1, const Vector3D& vec2);

    friend Vector3D operator-(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D& operator-=(Vector3D& vec1, const Vector3D& vec2);

    friend double operator*(const Vector3D& vec1, const Vector3D& vec2);

    friend Vector3D operator*(const Vector3D& vec1, const double factor1);
    friend Vector3D& operator*=(Vector3D& vec1, const double factor1);
    friend Vector3D operator*(const Vector3D& vec1, const Matrix3D& mat1);
    friend Vector3D& operator*=(Vector3D& vec1, const Matrix3D& mat1);

    friend Vector3D operator*(const double factor1, const Vector3D& vec1);
    friend Vector3D operator*(const Matrix3D& mat1, const Vector3D& vec1);

    friend Vector3D operator/(const Vector3D& vec1, const double factor1);
    friend Vector3D& operator/=(Vector3D& vec1, const double factor1);

    friend double scalar_product(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D vector_product(const Vector3D& vec1, const Vector3D& vec2);
    friend Vector3D cross_product(const Vector3D& vec1, const Vector3D& vec2);

    double magnitude() const;
    void normalize();
    Vector3D normalized() const;
    void deflect(const double, const double);
    void invert();
    Vector3D inverted() const;

    struct CartesianCoordinates {
        CartesianCoordinates() {};
        CartesianCoordinates(double, double, double);
        CartesianCoordinates(const CartesianCoordinates&);

        double x_, y_, z_;

        template<typename Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(cereal::make_nvp("X", x_));
                archive(cereal::make_nvp("Y", y_));
                archive(cereal::make_nvp("Z", z_));
            } else {
                throw std::runtime_error("CartesianCoordinates only supports version <= 0!");
            }
        }
    };

    struct SphericalCoordinates {
        SphericalCoordinates() {};
        SphericalCoordinates(double, double, double);
        SphericalCoordinates(const SphericalCoordinates&);

        double radius_, azimuth_, zenith_;

        template<typename Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(cereal::make_nvp("Radius", radius_));
                archive(cereal::make_nvp("Azimuth", azimuth_));
                archive(cereal::make_nvp("Zenith", zenith_));
            } else {
                throw std::runtime_error("SphericalCoordinates only supports version <= 0!");
            }
        }
    };

    //-------------------------------------//
    // conversions to spherical coordinates
    void CalculateCartesianFromSpherical();
    void CalculateSphericalCoordinates();

    //-------------------------------------//
    // setter
    void SetCartesianCoordinates(const double x, const double y, const double z)
    {
        cartesian_.x_ = x;
        cartesian_.y_ = y;
        cartesian_.z_ = z;
    }
    void SetSphericalCoordinates(const double radius, const double azimuth, const double zenith)
    {
        spherical_.radius_  = radius;
        spherical_.azimuth_ = azimuth;
        spherical_.zenith_  = zenith;
    }

    //-------------------------------------//
    // getter
    double GetX() const { return cartesian_.x_; }
    double GetY() const { return cartesian_.y_; }
    double GetZ() const { return cartesian_.z_; }
    double GetRadius() const { return spherical_.radius_; }
    double GetPhi() const { return spherical_.azimuth_; }
    double GetTheta() const { return spherical_.zenith_; }
    Vector3D::CartesianCoordinates GetCartesianCoordinates() const { return cartesian_; }
    Vector3D::SphericalCoordinates GetSphericalCoordinates() const { return spherical_; }

    //-------------------------------------//
    // serialization
    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("CartesianCoordinates", cartesian_));
            archive(cereal::make_nvp("SphericalCoordinates", spherical_));
        } else {
            throw std::runtime_error("Vector3D only supports version <= 0!");
        }
    }

    //----------------------------------------------//
private:
    CartesianCoordinates cartesian_;
    SphericalCoordinates spherical_;
};

double scalar_product(const Vector3D& vec1, const Vector3D& vec2);
Vector3D vector_product(const Vector3D& vec1, const Vector3D& vec2);
Vector3D cross_product(const Vector3D& vec1, const Vector3D& vec2);

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::Vector3D, 0);

#endif // LI_Vector3D_H

