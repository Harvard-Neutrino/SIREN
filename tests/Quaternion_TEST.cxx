
#include <cmath>
#include <random>
#include <iostream>

#include <gtest/gtest.h>

#include "earthmodel-service/Matrix3D.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Quaternion.h"

using namespace earthmodel;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Constructor, Default)
{
    ASSERT_NO_THROW(Quaternion());
}

TEST(Constructor, Coordinates)
{
    ASSERT_NO_THROW(Quaternion(0.0, 0.0, 0.0, 0.0));
    ASSERT_NO_THROW(Quaternion(1.0, 1.0, 1.0, 1.0));
    ASSERT_NO_THROW(Quaternion(-1.0, -1.0, -1.0, -1.0));
}

TEST(Constructor, ConstRef)
{
    Quaternion const A;
    ASSERT_NO_THROW(Quaternion(A));
}

TEST(Constructor, RRef)
{
    Quaternion const A;
    ASSERT_NO_THROW(Quaternion(Quaterion()));
}

TEST(Constructor, DefaultCheckWXYZ)
{
    Quaternion A;
    EXPECT_DOUBLE_EQ(1.0, A.GetW());
    EXPECT_DOUBLE_EQ(0.0, A.GetX());
    EXPECT_DOUBLE_EQ(0.0, A.GetY());
    EXPECT_DOUBLE_EQ(0.0, A.GetZ());
}

TEST(Constructor, CoordinatesCheckWXYZ)
{
    double w = 1.0;
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    Quaternion A(x, y, z, w);
    EXPECT_DOUBLE_EQ(w, A.GetW());
    EXPECT_DOUBLE_EQ(x, A.GetX());
    EXPECT_DOUBLE_EQ(y, A.GetY());
    EXPECT_DOUBLE_EQ(z, A.GetZ());
}

TEST(Constructor, ConstRefCheckWXYZ)
{
    double w = 1.0;
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    Quaternion const A(x, y, z, w);
    Quaternion B(A);
    EXPECT_DOUBLE_EQ(w, B.GetW());
    EXPECT_DOUBLE_EQ(x, B.GetX());
    EXPECT_DOUBLE_EQ(y, B.GetY());
    EXPECT_DOUBLE_EQ(z, B.GetZ());
}

TEST(Constructor, RRefCheckWXYZ)
{
    double w = 1.0;
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    Quaternion A(Quaternion(x, y, z, w));
    EXPECT_DOUBLE_EQ(w, A.GetW());
    EXPECT_DOUBLE_EQ(x, A.GetX());
    EXPECT_DOUBLE_EQ(y, A.GetY());
    EXPECT_DOUBLE_EQ(z, A.GetZ());
}

TEST(SetGet, CheckWXYZ)
{
    double w = 1.0;
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    Quaternion A;
    A.SetX(x);
    A.SetY(y);
    A.SetZ(z);
    A.SetW(w);
    EXPECT_DOUBLE_EQ(w, A.GetW());
    EXPECT_DOUBLE_EQ(x, A.GetX());
    EXPECT_DOUBLE_EQ(y, A.GetY());
    EXPECT_DOUBLE_EQ(z, A.GetZ());
}

TEST(SetGet, SetPosition)
{
    double w = 0.0;
    double x = 2.0;
    double y = 3.0;
    double z = 4.0;
    Vector3D vec(x, y, z);
    Quaternion A;
    A.SetPosition(vec);
    EXPECT_DOUBLE_EQ(w, A.GetW());
    EXPECT_DOUBLE_EQ(x, A.GetX());
    EXPECT_DOUBLE_EQ(y, A.GetY());
    EXPECT_DOUBLE_EQ(z, A.GetZ());
}

TEST(Matrix, IdentityRoundtrip)
{
    Quaternion A;
    Matrix3D B;
    B = A.GetMatrix();
    EXPECT_DOUBLE_EQ(1.0, B.GetXX());
    EXPECT_DOUBLE_EQ(0.0, B.GetXY());
    EXPECT_DOUBLE_EQ(0.0, B.GetXZ());
    EXPECT_DOUBLE_EQ(0.0, B.GetYX());
    EXPECT_DOUBLE_EQ(1.0, B.GetYY());
    EXPECT_DOUBLE_EQ(0.0, B.GetYZ());
    EXPECT_DOUBLE_EQ(0.0, B.GetZX());
    EXPECT_DOUBLE_EQ(0.0, B.GetZY());
    EXPECT_DOUBLE_EQ(1.0, B.GetZZ());

    A.SetMatrix(B);

    EXPECT_DOUBLE_EQ(1.0, A.GetW());
    EXPECT_DOUBLE_EQ(0.0, A.GetX());
    EXPECT_DOUBLE_EQ(0.0, A.GetY());
    EXPECT_DOUBLE_EQ(0.0, A.GetZ());
}

TEST(Quaternion, Magnitude)
{
    unsigned int n_rand = 100;
    for(unsigned int i=0; i<n_rand; ++i) {
        double w = RandomDouble() * 20 - 10;
        double x = RandomDouble() * 20 - 10;
        double y = RandomDouble() * 20 - 10;
        double z = RandomDouble() * 20 - 10;
        Quaternion A(x, y, z, w);
        double norm = A.magnitude();
        EXPECT_DOUBLE_EQ(std::sqrt(w*w + x*x + y*y + z*z), norm);

        Quaternion B = A.normalized();
        EXPECT_DOUBLE_EQ(1.0, B.magnitude());

        A.normalize();
        EXPECT_DOUBLE_EQ(1.0, A.magnitude());
    }
}

TEST(Quaternion, Invert)
{
    unsigned int n_rand = 100;
    for(unsigned int i=0; i<n_rand; ++i) {
        double w = RandomDouble() * 20 - 10;
        double x = RandomDouble() * 20 - 10;
        double y = RandomDouble() * 20 - 10;
        double z = RandomDouble() * 20 - 10;
        double norm2 = w*w + x*x + y*y + z*z;
        Quaternion A(x, y, z, w);

        Quaternion B = A.inverted();
        EXPECT_DOUBLE_EQ(w/norm2, B.GetW());
        EXPECT_DOUBLE_EQ(-x/norm2, B.GetX());
        EXPECT_DOUBLE_EQ(-y/norm2, B.GetY());
        EXPECT_DOUBLE_EQ(-z/norm2, B.GetZ());

        Quaternion C = A*B;
        EXPECT_NEAR(1.0, C.GetW(), 1e-12);
        EXPECT_NEAR(0.0, C.GetX(), 1e-12);
        EXPECT_NEAR(0.0, C.GetY(), 1e-12);
        EXPECT_NEAR(0.0, C.GetZ(), 1e-12);

        A.invert();
        EXPECT_DOUBLE_EQ(w/norm2, A.GetW());
        EXPECT_DOUBLE_EQ(-x/norm2, A.GetX());
        EXPECT_DOUBLE_EQ(-y/norm2, A.GetY());
        EXPECT_DOUBLE_EQ(-z/norm2, A.GetZ());
    }
}

TEST(Quaternion, Conjugate)
{
    unsigned int n_rand = 100;
    for(unsigned int i=0; i<n_rand; ++i) {
        double w = RandomDouble() * 20 - 10;
        double x = RandomDouble() * 20 - 10;
        double y = RandomDouble() * 20 - 10;
        double z = RandomDouble() * 20 - 10;
        double norm2 = w*w + x*x + y*y + z*z;
        Quaternion A(x, y, z, w);

        Quaternion B = A.conjugated();
        EXPECT_DOUBLE_EQ(w, B.GetW());
        EXPECT_DOUBLE_EQ(-x, B.GetX());
        EXPECT_DOUBLE_EQ(-y, B.GetY());
        EXPECT_DOUBLE_EQ(-z, B.GetZ());

        Quaternion C = A*B;
        EXPECT_NEAR(norm2, C.GetW(), 1e-12);
        EXPECT_NEAR(0.0, C.GetX(), 1e-12);
        EXPECT_NEAR(0.0, C.GetY(), 1e-12);
        EXPECT_NEAR(0.0, C.GetZ(), 1e-12);

        A.conjugate();
        EXPECT_DOUBLE_EQ(w, A.GetW());
        EXPECT_DOUBLE_EQ(-x, A.GetX());
        EXPECT_DOUBLE_EQ(-y, A.GetY());
        EXPECT_DOUBLE_EQ(-z, A.GetZ());
    }
}

TEST(Quaternion, Normalize)
{
    unsigned int n_rand = 100;
    for(unsigned int i=0; i<n_rand; ++i) {
        double w = RandomDouble() * 20 - 10;
        double x = RandomDouble() * 20 - 10;
        double y = RandomDouble() * 20 - 10;
        double z = RandomDouble() * 20 - 10;
        double norm2 = w*w + x*x + y*y + z*z;
        double norm = std::sqrt(norm2);
        Quaternion A(x, y, z, w);

        Quaternion B = A.normalized();
        EXPECT_DOUBLE_EQ(w/norm, B.GetW());
        EXPECT_DOUBLE_EQ(x/norm, B.GetX());
        EXPECT_DOUBLE_EQ(y/norm, B.GetY());
        EXPECT_DOUBLE_EQ(z/norm, B.GetZ());

        Quaternion C = A*B.inverted();
        EXPECT_NEAR(norm, C.GetW(), 1e-12);
        EXPECT_NEAR(0.0, C.GetX(), 1e-12);
        EXPECT_NEAR(0.0, C.GetY(), 1e-12);
        EXPECT_NEAR(0.0, C.GetZ(), 1e-12);

        A.normalize();
        EXPECT_DOUBLE_EQ(w/norm, A.GetW());
        EXPECT_DOUBLE_EQ(x/norm, A.GetX());
        EXPECT_DOUBLE_EQ(y/norm, A.GetY());
        EXPECT_DOUBLE_EQ(z/norm, A.GetZ());
    }
}

TEST(Quaternion, DotProduct)
{
    unsigned int n_rand = 100;
    for(unsigned int i=0; i<n_rand; ++i) {
        double wa = RandomDouble() * 20 - 10;
        double xa = RandomDouble() * 20 - 10;
        double ya = RandomDouble() * 20 - 10;
        double za = RandomDouble() * 20 - 10;
        Quaternion A(xa, ya, za, wa);

        double wb = RandomDouble() * 20 - 10;
        double xb = RandomDouble() * 20 - 10;
        double yb = RandomDouble() * 20 - 10;
        double zb = RandomDouble() * 20 - 10;
        Quaternion B(xb, yb, zb, wb);

        double C = A.DotProduct(B);
        double D = B.DotProduct(A);

        EXPECT_DOUBLE_EQ(C, D);

        double dot = xa*xb + ya*yb + za*zb + wa*wb;

        EXPECT_DOUBLE_EQ(dot, C);
        EXPECT_DOUBLE_EQ(dot, D);
    }
}

TEST(Quaternion, AxisAngle)
{
    unsigned int n_rand = 100;
    for(unsigned int i=0; i<n_rand; ++i) {
        double nz = RandomDouble() * 2 - 1;
        double phi = RandomDouble() * 2 * M_PI;
        double nrho = std::sqrt(1.0 - nz*nz);
        double nx = nrho * std::cos(phi);
        double ny = nrho * std::sin(phi);
        Vector3D axis(nx, ny, nz);

        double theta = RandomDouble() * 2 * M_PI;

        Quaternion A;
        A.SetAxisAngle(axis, theta);

        Vector3D vec;
        vec.SetCartesianCoordinates(A.GetX(), A.GetY(), A.GetZ());

        double sin_angle = vec.magnitude();
        double cos_angle = A.GetW();

        EXPECT_DOUBLE_EQ(std::sin(theta / 2.0), sin_angle);
        EXPECT_DOUBLE_EQ(std::cos(theta / 2.0), cos_angle);

        Vector3D normed_vec = vec.normalized();

        EXPECT_DOUBLE_EQ(axis.magnitude(), normed_vec.magnitude());

        EXPECT_DOUBLE_EQ(axis.GetX(), normed_vec.GetX());
        EXPECT_DOUBLE_EQ(axis.GetY(), normed_vec.GetY());
        EXPECT_DOUBLE_EQ(axis.GetZ(), normed_vec.GetZ());

        Vector3D new_axis;
        double new_angle;
        A.GetAxisAngle(new_axis, new_angle);

        EXPECT_DOUBLE_EQ(axis.GetX(), new_axis.GetX());
        EXPECT_DOUBLE_EQ(axis.GetY(), new_axis.GetY());
        EXPECT_DOUBLE_EQ(axis.GetZ(), new_axis.GetZ());
        EXPECT_DOUBLE_EQ(theta, new_angle);

        std::tuple<Vector3D, double> tup = A.GetAxisAngle();

        new_axis = std::get<0>(tup);
        new_angle = std::get<1>(tup);

        EXPECT_DOUBLE_EQ(axis.GetX(), new_axis.GetX());
        EXPECT_DOUBLE_EQ(axis.GetY(), new_axis.GetY());
        EXPECT_DOUBLE_EQ(axis.GetZ(), new_axis.GetZ());
        EXPECT_DOUBLE_EQ(theta, new_angle);
    }
}

TEST(Quaternion, Parity)
{
    Quaternion x(1, 0, 0, 0);
    Quaternion y(0, 1, 0, 0);
    Quaternion z(0, 0, 1, 0);
    Quaternion w(0, 0, 0, 1);

    Quaternion res;

    res = x * y;
    EXPECT_EQ(1.0, res.GetZ());

    res = y * x;
    EXPECT_EQ(-1.0, res.GetZ());

    res = y * z;
    EXPECT_EQ(1.0, res.GetX());

    res = z * y;
    EXPECT_EQ(-1.0, res.GetX());

    res = z * x;
    EXPECT_EQ(1.0, res.GetY());

    res = x * z;
    EXPECT_EQ(-1.0, res.GetY());

    res = w * x;
    EXPECT_EQ(1.0, res.GetX());

    res = x * w;
    EXPECT_EQ(1.0, res.GetX());

    res = w * y;
    EXPECT_EQ(1.0, res.GetY());

    res = y * w;
    EXPECT_EQ(1.0, res.GetY());

    res = w * z;
    EXPECT_EQ(1.0, res.GetZ());

    res = z * w;
    EXPECT_EQ(1.0, res.GetZ());

    res = w * w;
    EXPECT_EQ(1.0, res.GetW());

    res = x * x;
    EXPECT_EQ(-1.0, res.GetW());

    res = y * y;
    EXPECT_EQ(-1.0, res.GetW());

    res = z * z;
    EXPECT_EQ(-1.0, res.GetW());
}

TEST(Quaternion, SimpleRotation)
{
    Vector3D res;
    Vector3D vec;
    Vector3D axis;
    double angle;
    Quaternion rotor;
    Quaternion q_res;

    // pos x
    vec = Vector3D(1, 0, 0);

    // null rotation
    angle = 0;
    axis = Vector3D(0, 0, 1);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_DOUBLE_EQ(1.0, res.GetX());
    EXPECT_DOUBLE_EQ(0.0, res.GetY());
    EXPECT_DOUBLE_EQ(0.0, res.GetZ());
    res = rotor.rotate(vec, true);
    EXPECT_DOUBLE_EQ(1.0, res.GetX());
    EXPECT_DOUBLE_EQ(0.0, res.GetY());
    EXPECT_DOUBLE_EQ(0.0, res.GetZ());
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_DOUBLE_EQ(1.0, q_res.GetX());
    EXPECT_DOUBLE_EQ(0.0, q_res.GetY());
    EXPECT_DOUBLE_EQ(0.0, q_res.GetZ());

    angle = M_PI * 0.5;

    // about z
    axis = Vector3D(0, 0, 1);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(1.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(-1.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(1.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetZ(), 1e-12);

    // about y
    axis = Vector3D(0, 1, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(-1.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(-1.0, q_res.GetZ(), 1e-12);

    // about x
    axis = Vector3D(1, 0, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(1.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(1.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(1.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetZ(), 1e-12);

    // pos y
    vec = Vector3D(0, 1, 0);

    // null rotation
    angle = 0;
    axis = Vector3D(1, 0, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_DOUBLE_EQ(0.0, res.GetX());
    EXPECT_DOUBLE_EQ(1.0, res.GetY());
    EXPECT_DOUBLE_EQ(0.0, res.GetZ());
    res = rotor.rotate(vec, true);
    EXPECT_DOUBLE_EQ(0.0, res.GetX());
    EXPECT_DOUBLE_EQ(1.0, res.GetY());
    EXPECT_DOUBLE_EQ(0.0, res.GetZ());
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_DOUBLE_EQ(0.0, q_res.GetX());
    EXPECT_DOUBLE_EQ(1.0, q_res.GetY());
    EXPECT_DOUBLE_EQ(0.0, q_res.GetZ());

    angle = M_PI * 0.5;

    // about z
    axis = Vector3D(0, 0, 1);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(-1.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(1.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(-1.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetZ(), 1e-12);

    // about y
    axis = Vector3D(0, 1, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(1.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(1.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(1.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetZ(), 1e-12);

    // about x
    axis = Vector3D(1, 0, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(-1.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, q_res.GetZ(), 1e-12);

    // pos z
    vec = Vector3D(0, 0, 1);

    // null rotation
    angle = 0;
    axis = Vector3D(0, 1, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_DOUBLE_EQ(0.0, res.GetX());
    EXPECT_DOUBLE_EQ(0.0, res.GetY());
    EXPECT_DOUBLE_EQ(1.0, res.GetZ());
    res = rotor.rotate(vec, true);
    EXPECT_DOUBLE_EQ(0.0, res.GetX());
    EXPECT_DOUBLE_EQ(0.0, res.GetY());
    EXPECT_DOUBLE_EQ(1.0, res.GetZ());
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_DOUBLE_EQ(0.0, q_res.GetX());
    EXPECT_DOUBLE_EQ(0.0, q_res.GetY());
    EXPECT_DOUBLE_EQ(1.0, q_res.GetZ());

    angle = M_PI * 0.5;

    // about z
    axis = Vector3D(0, 0, 1);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, q_res.GetZ(), 1e-12);

    // about y
    axis = Vector3D(0, 1, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(1.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(-1.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(1.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetZ(), 1e-12);

    // about x
    axis = Vector3D(1, 0, 0);
    rotor.SetAxisAngle(axis, angle);
    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(-1.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    res = rotor.rotate(vec, true);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(1.0, res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, res.GetZ(), 1e-12);
    q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(-1.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetZ(), 1e-12);
}

TEST(Quaternion, RotationComposition)
{
    Vector3D res;
    Vector3D vec(0, 1, 0);
    double x_angle = M_PI * 0.5;
    double z_angle = M_PI * 0.5;
    Vector3D x_axis(1, 0, 0);
    Vector3D z_axis(0, 0, 1);
    Quaternion x_rotor;
    x_rotor.SetAxisAngle(x_axis, x_angle);
    Quaternion z_rotor;
    z_rotor.SetAxisAngle(z_axis, z_angle);

    Quaternion rotor = z_rotor * x_rotor;

    Quaternion q_res = rotor * Quaternion(vec.GetX(), vec.GetY(), vec.GetZ(), 0) * (rotor.inverted());
    EXPECT_NEAR(0.0, q_res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, q_res.GetZ(), 1e-12);
    EXPECT_NEAR(0.0, q_res.GetW(), 1e-12);

    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);

    res = x_rotor.rotate(vec, false);
    res = z_rotor.rotate(res, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);
}

TEST(Quaternion, EulerAngleConversion)
{
    Vector3D res;
    Vector3D vec(0, 1, 0);
    double x_angle = M_PI * 0.5;
    double z_angle = M_PI * 0.5;
    Vector3D x_axis(1, 0, 0);
    Vector3D z_axis(0, 0, 1);
    Quaternion x_rotor;
    x_rotor.SetAxisAngle(x_axis, x_angle);
    Quaternion z_rotor;
    z_rotor.SetAxisAngle(z_axis, z_angle);

    Quaternion rotor = z_rotor * x_rotor;

    Quaternion q;
    q.SetEulerAnglesZXZr(z_angle, x_angle, 0.0);
    EXPECT_NEAR(rotor.GetX(), q.GetX(), 1e-12);
    EXPECT_NEAR(rotor.GetY(), q.GetY(), 1e-12);
    EXPECT_NEAR(rotor.GetZ(), q.GetZ(), 1e-12);
    EXPECT_NEAR(rotor.GetW(), q.GetW(), 1e-12);

    res = rotor.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);

    res = q.rotate(vec, false);
    EXPECT_NEAR(0.0, res.GetX(), 1e-12);
    EXPECT_NEAR(0.0, res.GetY(), 1e-12);
    EXPECT_NEAR(1.0, res.GetZ(), 1e-12);

    double alpha, beta, gamma;
    q.GetEulerAnglesZXZr(alpha, beta, gamma);
    EXPECT_NEAR(z_angle, alpha, 1e-12);
    EXPECT_NEAR(x_angle, beta, 1e-12);
    EXPECT_NEAR(0.0, gamma, 1e-12);

    unsigned int n_rand = 100;

    for(unsigned int i=0; i<n_rand; ++i) {
        double z0 = (RandomDouble() * 2 - 1) * M_PI;
        double x0 = (RandomDouble() * M_PI);
        double z1 = (RandomDouble() * 2 - 1) * M_PI;
        q.SetEulerAnglesZXZr(z0, x0, z1);
        q.GetEulerAnglesZXZr(alpha, beta, gamma);
        EXPECT_NEAR(z0, alpha, 1e-12);
        EXPECT_NEAR(x0, beta, 1e-12);
        EXPECT_NEAR(z1, gamma, 1e-12);

        Quaternion z0_rotor;
        Quaternion x0_rotor;
        Quaternion z1_rotor;
        z0_rotor.SetAxisAngle(z_axis, z0);
        x0_rotor.SetAxisAngle(x_axis, x0);
        z1_rotor.SetAxisAngle(z_axis, z1);
        rotor = z0_rotor * x0_rotor * z1_rotor;
        EXPECT_NEAR(rotor.GetX(), q.GetX(), 1e-12);
        EXPECT_NEAR(rotor.GetY(), q.GetY(), 1e-12);
        EXPECT_NEAR(rotor.GetZ(), q.GetZ(), 1e-12);
        EXPECT_NEAR(rotor.GetW(), q.GetW(), 1e-12);

        q.SetEulerAngles(EulerAngles(EulerOrder::ZXZr, z0, x0, z1));
        EXPECT_NEAR(rotor.GetX(), q.GetX(), 1e-12);
        EXPECT_NEAR(rotor.GetY(), q.GetY(), 1e-12);
        EXPECT_NEAR(rotor.GetZ(), q.GetZ(), 1e-12);
        EXPECT_NEAR(rotor.GetW(), q.GetW(), 1e-12);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

