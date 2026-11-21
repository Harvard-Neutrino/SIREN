
#include <cmath>
#include <iostream>

#include <gtest/gtest.h>

#include "SIREN/math/Vector3D.h"

using namespace siren::math;

TEST(Comparison, Comparison_equal)
{
    Vector3D A;
    Vector3D B;
    EXPECT_TRUE(A == B);
    Vector3D* C = new Vector3D(1., 2., 3.);
    Vector3D* D = new Vector3D(1., 2., 3.);
    EXPECT_TRUE(*C == *D);
    D->SetCartesianCoordinates(0., 0., 0.);
    EXPECT_TRUE(A == *D);
    B.SetSphericalCoordinates(0.1, 0.2, 0.3);
    D->SetSphericalCoordinates(0.1, 0.2, 0.3);
    EXPECT_TRUE(B == *D);
}

TEST(Comparison, Comparison_not_equal)
{
    Vector3D A;
    Vector3D B;
    B.SetCartesianCoordinates(1., 2., 3.);
    EXPECT_TRUE(A != B);
    Vector3D C;
    C.SetSphericalCoordinates(0.1, 0.2, 0.3);
    EXPECT_TRUE(A != C);
}

TEST(Assignment, Copyconstructor)
{
    Vector3D A(1, 2, 3);
    Vector3D B(A);
    EXPECT_TRUE(A == B);
    Vector3D C(2, 3, 4);
    B = C;
    EXPECT_TRUE(C == B);
}

TEST(Assignment, Operator)
{
    Vector3D A;
    Vector3D B;
    A.SetCartesianCoordinates(1, 2, 3);
    EXPECT_TRUE(A != B);
    B = A;
    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    Vector3D A;
    Vector3D B;
    EXPECT_TRUE(A == B);
    Vector3D* C = new Vector3D(1, 2, 3);
    Vector3D* D = new Vector3D(1, 2, 3);
    C->SetSphericalCoordinates(0.1, 0.2, 0.3);
    D->SetSphericalCoordinates(0.1, 0.2, 0.3);
    EXPECT_TRUE(*C == *D);

    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
}

TEST(Addition, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C.SetCartesianCoordinates(3, 6, 9);
    D = A + B;
    EXPECT_TRUE(C == D);
    D.SetCartesianCoordinates(0, 0, 0);
    EXPECT_TRUE(D != C);
    D = B + A;
    EXPECT_TRUE(D == C);
}

TEST(AdditionAssignment, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C.SetCartesianCoordinates(3, 6, 9);
    A += B;
    EXPECT_TRUE(C == A);
    A.SetCartesianCoordinates(1, 2, 3);
    EXPECT_TRUE(A != C);
    B += A;
    EXPECT_TRUE(B == C);
}

TEST(Subtraction, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    Vector3D E;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C.SetCartesianCoordinates(3, 6, 9);
    E.SetCartesianCoordinates(-1, -2, -3);
    D = -A;
    EXPECT_TRUE(D == E);
    D = C - B;
    EXPECT_TRUE(D == A);
    D = B - C;
    EXPECT_TRUE(D == E);
}

TEST(Scaling, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    double factor1 = 2.;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    C = A * factor1;
    EXPECT_TRUE(C == B);
    C.SetCartesianCoordinates(0, 0, 0);
    EXPECT_TRUE(C != B);
    C = factor1 * A;
    EXPECT_TRUE(C == B);
}

TEST(ScalingAssignment, Operator)
{
    Vector3D A;
    Vector3D B;
    double factor1 = 2.;
    A.SetCartesianCoordinates(1, 2, 3);
    B.SetCartesianCoordinates(2, 4, 6);
    A *= factor1;
    EXPECT_TRUE(A == B);
    A.SetCartesianCoordinates(1, 2, 3);
}

TEST(ScalarProduct, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    double factor1 = 5;
    double result  = 0;
    A.SetCartesianCoordinates(1, 1, 1);
    B.SetCartesianCoordinates(2, 1, 2);
    EXPECT_TRUE(result != factor1);
    result = scalar_product(A, B);
    EXPECT_TRUE(result == factor1);
    result = 0.;
    result = scalar_product(B, A);
    EXPECT_TRUE(result == factor1);
}

TEST(VectorProduct, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Vector3D D;
    Vector3D E;
    A.SetCartesianCoordinates(1, 0, 0);
    B.SetCartesianCoordinates(0, 1, 0);
    C.SetCartesianCoordinates(0, 0, 1);
    D = vector_product(A, B);
    EXPECT_TRUE(D == C);
    D = -vector_product(B, A);
    EXPECT_TRUE(D == C);
    D = -vector_product(A, C);
    EXPECT_TRUE(D == B);
    D = vector_product(B, C);
    EXPECT_TRUE(D == A);
}

TEST(Magnitude, Operator)
{
    Vector3D A;
    Vector3D B;
    double sum = 3;
    double result;
    A.SetCartesianCoordinates(1, 2, 2);
    B.SetCartesianCoordinates(1, 1, 1);
    result = A.magnitude();
    EXPECT_TRUE(result == sum);
    result = B.magnitude();
    EXPECT_TRUE(result != sum);
}

TEST(Normalize, Operator)
{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    A.SetCartesianCoordinates(3, 0, 0);
    B.SetCartesianCoordinates(1, 1, 1);
    C.SetCartesianCoordinates(1, 0, 0);
    EXPECT_TRUE(A != C);
    EXPECT_TRUE(B != C);
    A.normalize();
    C.SetSphericalCoordinates(1,0,0);
    EXPECT_TRUE(A == C);
    B.normalize();
    EXPECT_TRUE(B != C);
}

// Normalizing the zero vector must not divide by zero: stay finite and leave the
// vector unchanged at magnitude 0 (the contract degenerate DetectorModel p1 - p0
// directions rely on).
TEST(Normalize, ZeroVectorDoesNotProduceNaN)
{
    Vector3D Z;
    Z.SetCartesianCoordinates(0.0, 0.0, 0.0);
    // Must not abort and must not produce NaN/Inf.
    EXPECT_NO_THROW(Z.normalize());
    std::array<double, 3> z = Z;
    EXPECT_FALSE(std::isnan(z[0]));
    EXPECT_FALSE(std::isnan(z[1]));
    EXPECT_FALSE(std::isnan(z[2]));
    EXPECT_TRUE(std::isfinite(z[0]));
    EXPECT_TRUE(std::isfinite(z[1]));
    EXPECT_TRUE(std::isfinite(z[2]));
    // The zero vector is left unchanged (magnitude stays exactly 0).
    EXPECT_DOUBLE_EQ(0.0, z[0]);
    EXPECT_DOUBLE_EQ(0.0, z[1]);
    EXPECT_DOUBLE_EQ(0.0, z[2]);
    EXPECT_DOUBLE_EQ(0.0, Z.magnitude());
}

// normalized() is the const sibling and must also be NaN-safe on a zero vector.
TEST(Normalize, ZeroVectorNormalizedIsFinite)
{
    Vector3D Z;
    Z.SetCartesianCoordinates(0.0, 0.0, 0.0);
    Vector3D N = Z.normalized();
    std::array<double, 3> n = N;
    EXPECT_FALSE(std::isnan(n[0]) || std::isnan(n[1]) || std::isnan(n[2]));
    EXPECT_DOUBLE_EQ(0.0, N.magnitude());
}

// Boundary just above the zero-length guard: a tiny but nonzero difference of two
// Earth-scale coordinates must normalize to a true unit vector, not be left as-is.
TEST(Normalize, TinyEarthScaleDifferenceIsFiniteUnit)
{
    // 1e-7 m is exactly representable relative to 6.371e6 m, so the subtraction is
    // exact, the magnitude nonzero, and normalize() must yield a unit vector.
    Vector3D p0;
    p0.SetCartesianCoordinates(6.371e6, 0.0, 0.0);
    Vector3D p1;
    p1.SetCartesianCoordinates(6.371e6 + 1e-7, 0.0, 0.0);
    Vector3D direction = p1 - p0;
    double mag = direction.magnitude();
    ASSERT_GT(mag, 0.0);
    direction.normalize();
    std::array<double, 3> d = direction;
    EXPECT_FALSE(std::isnan(d[0]) || std::isnan(d[1]) || std::isnan(d[2]));
    // Resulting vector is unit length.
    EXPECT_NEAR(1.0, direction.magnitude(), 1e-12);
    EXPECT_NEAR(1.0, d[0], 1e-12);
}

TEST(CalculateSphericalCoordinates, Conversion)
{
    Vector3D A;
    Vector3D B(1, 2, 2);
    A.SetCartesianCoordinates(1, 2, 2);
    EXPECT_TRUE(A == B);
    A.CalculateSphericalCoordinates();
    EXPECT_TRUE(A != B);
    B.SetSphericalCoordinates(3, std::atan2(2., 1.), std::acos(2. / 3.));
    double epsilon = std::numeric_limits<double>::epsilon();
    double error_factor = 2.;
    bool test_x = std::abs(A.GetX() - B.GetX()) < std::max(std::abs(A.GetX()), std::abs(B.GetX())) * epsilon * error_factor;
    bool test_y = std::abs(A.GetY() - B.GetY()) < std::max(std::abs(A.GetY()), std::abs(B.GetY())) * epsilon * error_factor;
    bool test_z = std::abs(A.GetZ() - B.GetZ()) < std::max(std::abs(A.GetZ()), std::abs(B.GetZ())) * epsilon * error_factor;
    bool test_theta = std::abs(A.GetTheta() - B.GetTheta()) < std::max(std::abs(A.GetTheta()), std::abs(B.GetTheta())) * epsilon * error_factor;
    bool test_phi = std::abs(A.GetPhi() - B.GetPhi()) < std::max(std::abs(A.GetPhi()), std::abs(B.GetPhi())) * epsilon * error_factor;
    EXPECT_TRUE(B == A || (test_x && test_y && test_z && test_theta && test_phi));
}

TEST(CalculateCartesianFromSpherical, Conversion)
{
    Vector3D A;
    Vector3D B;
    double epsilon = std::numeric_limits<double>::epsilon();
    A.SetCartesianCoordinates(1, 2, 2);
    B.SetSphericalCoordinates(3, std::atan2(2., 1.), std::acos(2. / 3.));
    B.CalculateCartesianFromSpherical();
    EXPECT_TRUE(A != B);
    B.SetSphericalCoordinates(0, 0, 0);
    double error_factor = 10.;
    bool test_x = std::abs(A.GetX() - B.GetX()) < std::max(std::abs(A.GetX()), std::abs(B.GetX())) * epsilon * error_factor;
    bool test_y = std::abs(A.GetY() - B.GetY()) < std::max(std::abs(A.GetY()), std::abs(B.GetY())) * epsilon * error_factor;
    bool test_z = std::abs(A.GetZ() - B.GetZ()) < std::max(std::abs(A.GetZ()), std::abs(B.GetZ())) * epsilon * error_factor;
    EXPECT_TRUE(A == B || (test_x && test_y && test_z)) << A << B;
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

