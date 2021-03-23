
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>
#include <gtest/gtest.h>

#include "earthmodel-service/DensityDist.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Polynomial.h"

using namespace earthmodel;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

Vector3D RandomVector() {
    return Vector3D(RandomDouble()*4-2, RandomDouble()*4-2, RandomDouble()*4-2);
}

Vector3D RandomDirection() {
    double az = RandomDouble()*2*M_PI;
    double zen = RandomDouble()*M_PI;
    Vector3D res;
    res.SetSphericalCoordinates(1.0, az, zen);
    res.CalculateCartesianFromSpherical();
    return res;
}

TEST(Constructor, AxisDistribution)
{
    CartesianAxis1D ax_A;
    RadialAxis1D ax_B;

    ConstantDistribution1D dist_A;
    PolynomialDistribution1D dist_B({});
    ExponentialDistribution1D dist_C(1.0);

    auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
    auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
    auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
    auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
    auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
    auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);
}

TEST(Constructor, Other)
{
    CartesianAxis1D ax_A;
    RadialAxis1D ax_B;

    ConstantDistribution1D dist_A;
    PolynomialDistribution1D dist_B({});
    ExponentialDistribution1D dist_C(1.0);

    auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
    auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
    auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
    auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
    auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
    auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> Ao(A);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> Bo(B);
    DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D> Co(C);
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> Do(D);
    DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D> Eo(E);
    DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D> Fo(F);
}

TEST(Constructor, Constant)
{
    CartesianAxis1D ax_A;
    RadialAxis1D ax_B;

    ConstantDistribution1D dist_A;

    auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
    auto B = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> Ao(A);
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> Bo(B);

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> Ad;
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> Bd;

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> Av(ax_A, 1.0);
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> Bv(ax_B, 1.0);
}

TEST(Constructor, Cartesian)
{
    CartesianAxis1D ax_A;

    PolynomialDistribution1D dist_A({});
    ExponentialDistribution1D dist_B(1.0);

    auto A = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_A);
    auto B = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_B);

    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> Ao(A);
    DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D> Bo(B);
}

TEST(Constructor, RadialPolynomial)
{
    using T = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>;

    RadialAxis1D ax_A;
    PolynomialDistribution1D dist_A({});
    Polynom poly({});
    std::vector<double> params;

    // T A; // Not implemented
    T B(ax_A, dist_A);
    T C(ax_A, poly);
    T D(ax_A, params);
    T E(B);
}

TEST(Comparison, Comparison_equal)
{
    CartesianAxis1D ax_A0;
    RadialAxis1D ax_B0;
    ConstantDistribution1D dist_A0;
    PolynomialDistribution1D dist_B0({});
    ExponentialDistribution1D dist_C0(1.0);

    CartesianAxis1D ax_A1;
    RadialAxis1D ax_B1;
    ConstantDistribution1D dist_A1;
    PolynomialDistribution1D dist_B1({});
    ExponentialDistribution1D dist_C1(1.0);

    auto A0 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A0, dist_A0);
    auto B0 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A0, dist_B0);
    auto C0 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A0, dist_C0);
    auto D0 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B0, dist_A0);
    auto E0 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B0, dist_B0);
    auto F0 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B0, dist_C0);

    auto A1 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A1, dist_A1);
    auto B1 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A1, dist_B1);
    auto C1 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A1, dist_C1);
    auto D1 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B1, dist_A1);
    auto E1 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B1, dist_B1);
    auto F1 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B1, dist_C1);

    EXPECT_TRUE(A0 == A1);
    EXPECT_TRUE(B0 == B1);
    EXPECT_TRUE(C0 == C1);
    EXPECT_TRUE(D0 == D1);
    EXPECT_TRUE(E0 == E1);
    EXPECT_TRUE(F0 == F1);

    CartesianAxis1D ax_A2(Vector3D(0,0,1),Vector3D(0,1,0));
    RadialAxis1D ax_B2(Vector3D(0,1,0));
    ConstantDistribution1D dist_A2(2.0);
    PolynomialDistribution1D dist_B2({1.0});
    ExponentialDistribution1D dist_C2(2.0);

    auto A2 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A2, dist_A0);
    auto B2 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A2, dist_B0);
    auto C2 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A2, dist_C0);
    auto D2 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B2, dist_A0);
    auto E2 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B2, dist_B0);
    auto F2 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B2, dist_C0);

    EXPECT_FALSE(A0 == A2);
    EXPECT_FALSE(B0 == B2);
    EXPECT_FALSE(C0 == C2);
    EXPECT_FALSE(D0 == D2);
    EXPECT_FALSE(E0 == E2);
    EXPECT_FALSE(F0 == F2);

    auto A3 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A0, dist_A2);
    auto B3 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A0, dist_B2);
    auto C3 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A0, dist_C2);
    auto D3 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B0, dist_A2);
    auto E3 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B0, dist_B2);
    auto F3 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B0, dist_C2);

    EXPECT_FALSE(A0 == A3);
    EXPECT_FALSE(B0 == B3);
    EXPECT_FALSE(C0 == C3);
    EXPECT_FALSE(D0 == D3);
    EXPECT_FALSE(E0 == E3);
    EXPECT_FALSE(F0 == F3);

    EXPECT_FALSE(A0 == B0);
    EXPECT_FALSE(A0 == C0);
    EXPECT_FALSE(A0 == D0);
    EXPECT_FALSE(A0 == E0);
    EXPECT_FALSE(A0 == F0);
    EXPECT_FALSE(B0 == C0);
    EXPECT_FALSE(B0 == D0);
    EXPECT_FALSE(B0 == E0);
    EXPECT_FALSE(B0 == F0);
    EXPECT_FALSE(C0 == D0);
    EXPECT_FALSE(C0 == E0);
    EXPECT_FALSE(C0 == F0);
    EXPECT_FALSE(D0 == E0);
    EXPECT_FALSE(D0 == F0);
    EXPECT_FALSE(E0 == F0);

    DensityDistribution* pA0 = A0.clone();
    DensityDistribution* pB0 = B0.clone();
    DensityDistribution* pC0 = C0.clone();
    DensityDistribution* pD0 = D0.clone();
    DensityDistribution* pE0 = E0.clone();
    DensityDistribution* pF0 = F0.clone();
    DensityDistribution* pA1 = A1.clone();
    DensityDistribution* pB1 = B1.clone();
    DensityDistribution* pC1 = C1.clone();
    DensityDistribution* pD1 = D1.clone();
    DensityDistribution* pE1 = E1.clone();
    DensityDistribution* pF1 = F1.clone();
    DensityDistribution* pA2 = A2.clone();
    DensityDistribution* pB2 = B2.clone();
    DensityDistribution* pC2 = C2.clone();
    DensityDistribution* pD2 = D2.clone();
    DensityDistribution* pE2 = E2.clone();
    DensityDistribution* pF2 = F2.clone();
    DensityDistribution* pA3 = A3.clone();
    DensityDistribution* pB3 = B3.clone();
    DensityDistribution* pC3 = C3.clone();
    DensityDistribution* pD3 = D3.clone();
    DensityDistribution* pE3 = E3.clone();
    DensityDistribution* pF3 = F3.clone();

    EXPECT_TRUE(*pA0 == *pA1);
    EXPECT_TRUE(*pB0 == *pB1);
    EXPECT_TRUE(*pC0 == *pC1);
    EXPECT_TRUE(*pD0 == *pD1);
    EXPECT_TRUE(*pE0 == *pE1);
    EXPECT_TRUE(*pF0 == *pF1);

    EXPECT_FALSE(*pA0 == *pA2);
    EXPECT_FALSE(*pB0 == *pB2);
    EXPECT_FALSE(*pC0 == *pC2);
    EXPECT_FALSE(*pD0 == *pD2);
    EXPECT_FALSE(*pE0 == *pE2);
    EXPECT_FALSE(*pF0 == *pF2);

    EXPECT_FALSE(*pA0 == *pA3);
    EXPECT_FALSE(*pB0 == *pB3);
    EXPECT_FALSE(*pC0 == *pC3);
    EXPECT_FALSE(*pD0 == *pD3);
    EXPECT_FALSE(*pE0 == *pE3);
    EXPECT_FALSE(*pF0 == *pF3);

    EXPECT_FALSE(*pA0 == *pB0);
    EXPECT_FALSE(*pA0 == *pC0);
    EXPECT_FALSE(*pA0 == *pD0);
    EXPECT_FALSE(*pA0 == *pE0);
    EXPECT_FALSE(*pA0 == *pF0);
    EXPECT_FALSE(*pB0 == *pC0);
    EXPECT_FALSE(*pB0 == *pD0);
    EXPECT_FALSE(*pB0 == *pE0);
    EXPECT_FALSE(*pB0 == *pF0);
    EXPECT_FALSE(*pC0 == *pD0);
    EXPECT_FALSE(*pC0 == *pE0);
    EXPECT_FALSE(*pC0 == *pF0);
    EXPECT_FALSE(*pD0 == *pE0);
    EXPECT_FALSE(*pD0 == *pF0);
    EXPECT_FALSE(*pE0 == *pF0);

    delete pA0;
    delete pB0;
    delete pC0;
    delete pD0;
    delete pE0;
    delete pF0;

    delete pA1;
    delete pB1;
    delete pC1;
    delete pD1;
    delete pE1;
    delete pF1;

    delete pA2;
    delete pB2;
    delete pC2;
    delete pD2;
    delete pE2;
    delete pF2;

    delete pA3;
    delete pB3;
    delete pC3;
    delete pD3;
    delete pE3;
    delete pF3;
}

TEST(Comparison, Comparison_not_equal)
{
    CartesianAxis1D ax_A0;
    RadialAxis1D ax_B0;
    ConstantDistribution1D dist_A0;
    PolynomialDistribution1D dist_B0({});
    ExponentialDistribution1D dist_C0(1.0);

    CartesianAxis1D ax_A1;
    RadialAxis1D ax_B1;
    ConstantDistribution1D dist_A1;
    PolynomialDistribution1D dist_B1({});
    ExponentialDistribution1D dist_C1(1.0);

    auto A0 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A0, dist_A0);
    auto B0 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A0, dist_B0);
    auto C0 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A0, dist_C0);
    auto D0 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B0, dist_A0);
    auto E0 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B0, dist_B0);
    auto F0 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B0, dist_C0);

    auto A1 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A1, dist_A1);
    auto B1 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A1, dist_B1);
    auto C1 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A1, dist_C1);
    auto D1 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B1, dist_A1);
    auto E1 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B1, dist_B1);
    auto F1 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B1, dist_C1);

    EXPECT_FALSE(A0 != A1);
    EXPECT_FALSE(B0 != B1);
    EXPECT_FALSE(C0 != C1);
    EXPECT_FALSE(D0 != D1);
    EXPECT_FALSE(E0 != E1);
    EXPECT_FALSE(F0 != F1);

    CartesianAxis1D ax_A2(Vector3D(0,0,1),Vector3D(0,1,0));
    RadialAxis1D ax_B2(Vector3D(0,1,0));
    ConstantDistribution1D dist_A2(2.0);
    PolynomialDistribution1D dist_B2({1.0});
    ExponentialDistribution1D dist_C2(2.0);

    auto A2 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A2, dist_A0);
    auto B2 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A2, dist_B0);
    auto C2 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A2, dist_C0);
    auto D2 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B2, dist_A0);
    auto E2 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B2, dist_B0);
    auto F2 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B2, dist_C0);

    EXPECT_TRUE(A0 != A2);
    EXPECT_TRUE(B0 != B2);
    EXPECT_TRUE(C0 != C2);
    EXPECT_TRUE(D0 != D2);
    EXPECT_TRUE(E0 != E2);
    EXPECT_TRUE(F0 != F2);

    auto A3 = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A0, dist_A2);
    auto B3 = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A0, dist_B2);
    auto C3 = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A0, dist_C2);
    auto D3 = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B0, dist_A2);
    auto E3 = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B0, dist_B2);
    auto F3 = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B0, dist_C2);

    EXPECT_TRUE(A0 != A3);
    EXPECT_TRUE(B0 != B3);
    EXPECT_TRUE(C0 != C3);
    EXPECT_TRUE(D0 != D3);
    EXPECT_TRUE(E0 != E3);
    EXPECT_TRUE(F0 != F3);

    EXPECT_TRUE(A0 != B0);
    EXPECT_TRUE(A0 != C0);
    EXPECT_TRUE(A0 != D0);
    EXPECT_TRUE(A0 != E0);
    EXPECT_TRUE(A0 != F0);
    EXPECT_TRUE(B0 != C0);
    EXPECT_TRUE(B0 != D0);
    EXPECT_TRUE(B0 != E0);
    EXPECT_TRUE(B0 != F0);
    EXPECT_TRUE(C0 != D0);
    EXPECT_TRUE(C0 != E0);
    EXPECT_TRUE(C0 != F0);
    EXPECT_TRUE(D0 != E0);
    EXPECT_TRUE(D0 != F0);
    EXPECT_TRUE(E0 != F0);

    DensityDistribution* pA0 = A0.clone();
    DensityDistribution* pB0 = B0.clone();
    DensityDistribution* pC0 = C0.clone();
    DensityDistribution* pD0 = D0.clone();
    DensityDistribution* pE0 = E0.clone();
    DensityDistribution* pF0 = F0.clone();
    DensityDistribution* pA1 = A1.clone();
    DensityDistribution* pB1 = B1.clone();
    DensityDistribution* pC1 = C1.clone();
    DensityDistribution* pD1 = D1.clone();
    DensityDistribution* pE1 = E1.clone();
    DensityDistribution* pF1 = F1.clone();
    DensityDistribution* pA2 = A2.clone();
    DensityDistribution* pB2 = B2.clone();
    DensityDistribution* pC2 = C2.clone();
    DensityDistribution* pD2 = D2.clone();
    DensityDistribution* pE2 = E2.clone();
    DensityDistribution* pF2 = F2.clone();
    DensityDistribution* pA3 = A3.clone();
    DensityDistribution* pB3 = B3.clone();
    DensityDistribution* pC3 = C3.clone();
    DensityDistribution* pD3 = D3.clone();
    DensityDistribution* pE3 = E3.clone();
    DensityDistribution* pF3 = F3.clone();

    EXPECT_FALSE(*pA0 != *pA1);
    EXPECT_FALSE(*pB0 != *pB1);
    EXPECT_FALSE(*pC0 != *pC1);
    EXPECT_FALSE(*pD0 != *pD1);
    EXPECT_FALSE(*pE0 != *pE1);
    EXPECT_FALSE(*pF0 != *pF1);

    EXPECT_TRUE(*pA0 != *pA2);
    EXPECT_TRUE(*pB0 != *pB2);
    EXPECT_TRUE(*pC0 != *pC2);
    EXPECT_TRUE(*pD0 != *pD2);
    EXPECT_TRUE(*pE0 != *pE2);
    EXPECT_TRUE(*pF0 != *pF2);

    EXPECT_TRUE(*pA0 != *pA3);
    EXPECT_TRUE(*pB0 != *pB3);
    EXPECT_TRUE(*pC0 != *pC3);
    EXPECT_TRUE(*pD0 != *pD3);
    EXPECT_TRUE(*pE0 != *pE3);
    EXPECT_TRUE(*pF0 != *pF3);

    EXPECT_TRUE(*pA0 != *pB0);
    EXPECT_TRUE(*pA0 != *pC0);
    EXPECT_TRUE(*pA0 != *pD0);
    EXPECT_TRUE(*pA0 != *pE0);
    EXPECT_TRUE(*pA0 != *pF0);
    EXPECT_TRUE(*pB0 != *pC0);
    EXPECT_TRUE(*pB0 != *pD0);
    EXPECT_TRUE(*pB0 != *pE0);
    EXPECT_TRUE(*pB0 != *pF0);
    EXPECT_TRUE(*pC0 != *pD0);
    EXPECT_TRUE(*pC0 != *pE0);
    EXPECT_TRUE(*pC0 != *pF0);
    EXPECT_TRUE(*pD0 != *pE0);
    EXPECT_TRUE(*pD0 != *pF0);
    EXPECT_TRUE(*pE0 != *pF0);

    delete pA0;
    delete pB0;
    delete pC0;
    delete pD0;
    delete pE0;
    delete pF0;

    delete pA1;
    delete pB1;
    delete pC1;
    delete pD1;
    delete pE1;
    delete pF1;

    delete pA2;
    delete pB2;
    delete pC2;
    delete pD2;
    delete pE2;
    delete pF2;

    delete pA3;
    delete pB3;
    delete pC3;
    delete pD3;
    delete pE3;
    delete pF3;
}

TEST(Assignment, Copyconstructor)
{
    CartesianAxis1D ax_A;
    RadialAxis1D ax_B;

    ConstantDistribution1D dist_A;
    PolynomialDistribution1D dist_B({});
    ExponentialDistribution1D dist_C(1.0);

    auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
    auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
    auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
    auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
    auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
    auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> Ao(A);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> Bo(B);
    DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D> Co(C);
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> Do(D);
    DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D> Eo(E);
    DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D> Fo(F);

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> Aa = A;
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> Ba = B;
    DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D> Ca = C;
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> Da = D;
    DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D> Ea = E;
    DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D> Fa = F;

    EXPECT_TRUE(A == Ao);
    EXPECT_TRUE(B == Bo);
    EXPECT_TRUE(C == Co);
    EXPECT_TRUE(D == Do);
    EXPECT_TRUE(E == Eo);
    EXPECT_TRUE(F == Fo);
    EXPECT_TRUE(A == Aa);
    EXPECT_TRUE(B == Ba);
    EXPECT_TRUE(C == Ca);
    EXPECT_TRUE(D == Da);
    EXPECT_TRUE(E == Ea);
    EXPECT_TRUE(F == Fa);
}

TEST(Evaluate, Axis_to_Distribution_connection)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }

        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*4-2);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
        auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
        auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            EXPECT_DOUBLE_EQ(A.Evaluate(position), dist_A.Evaluate(ax_A.GetX(position)));
            EXPECT_DOUBLE_EQ(B.Evaluate(position), dist_B.Evaluate(ax_A.GetX(position)));
            EXPECT_DOUBLE_EQ(C.Evaluate(position), dist_C.Evaluate(ax_A.GetX(position)));
            EXPECT_DOUBLE_EQ(D.Evaluate(position), dist_A.Evaluate(ax_B.GetX(position)));
            EXPECT_DOUBLE_EQ(E.Evaluate(position), dist_B.Evaluate(ax_B.GetX(position)));
            EXPECT_DOUBLE_EQ(F.Evaluate(position), dist_C.Evaluate(ax_B.GetX(position)));
        }
    }
}

TEST(Evaluate, ConstantDistribution)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        double val = RandomDouble()*10;
        ConstantDistribution1D dist_A(val);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            EXPECT_DOUBLE_EQ(A.Evaluate(position), val);
            EXPECT_DOUBLE_EQ(B.Evaluate(position), val);
        }
    }
}

TEST(Evaluate, CartesianDistribution)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }
        Vector3D axis = RandomDirection();
        Vector3D center = RandomVector();
        CartesianAxis1D ax_A(axis, center);

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*4-2);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            double x = (position - center)*axis;
            EXPECT_DOUBLE_EQ(A.Evaluate(position), dist_A.Evaluate(x));
            EXPECT_DOUBLE_EQ(B.Evaluate(position), dist_B.Evaluate(x));
            EXPECT_DOUBLE_EQ(C.Evaluate(position), dist_C.Evaluate(x));
        }
    }
}

TEST(Evaluate, RadialPolynomial)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }

        std::function<double(double)> eval = [&](double x)->double {
            double res = params[n-1];
            for(int i=n-2; i>=0; --i) {
                res = res*x + params[i];
            }
            return res;
        };

        Vector3D center = RandomVector();
        RadialAxis1D ax_A(center);

        PolynomialDistribution1D dist_A(params);

        auto A = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_A, dist_A);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            double x = (position - center).magnitude();
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.Evaluate(position), res);
        }
    }

}

TEST(Derivative, Axis_to_Distribution_connection)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }

        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*4-2);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
        auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
        auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            Vector3D direction = RandomDirection();
            double x_A = ax_A.GetX(position);
            double x_B = ax_B.GetX(position);
            double dxdt_A = ax_A.GetdX(position, direction);
            double dxdt_B = ax_B.GetdX(position, direction);
            EXPECT_DOUBLE_EQ(A.Derivative(position, direction), dist_A.Derivative(x_A)*dxdt_A);
            EXPECT_DOUBLE_EQ(B.Derivative(position, direction), dist_B.Derivative(x_A)*dxdt_A);
            EXPECT_DOUBLE_EQ(C.Derivative(position, direction), dist_C.Derivative(x_A)*dxdt_A);
            EXPECT_DOUBLE_EQ(D.Derivative(position, direction), dist_A.Derivative(x_B)*dxdt_B);
            EXPECT_DOUBLE_EQ(E.Derivative(position, direction), dist_B.Derivative(x_B)*dxdt_B);
            EXPECT_DOUBLE_EQ(F.Derivative(position, direction), dist_C.Derivative(x_B)*dxdt_B);
        }
    }
}

TEST(Derivative, ConstantDistribution)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        double val = RandomDouble()*10;
        ConstantDistribution1D dist_A(val);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            Vector3D direction = RandomDirection();
            EXPECT_DOUBLE_EQ(A.Derivative(position, direction), 0.0);
            EXPECT_DOUBLE_EQ(B.Derivative(position, direction), 0.0);
        }
    }
}

TEST(Derivative, CartesianDistribution)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }
        Vector3D axis = RandomDirection();
        Vector3D center = RandomVector();
        CartesianAxis1D ax_A(axis, center);

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*4-2);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            Vector3D direction = RandomDirection();
            double alpha = axis*direction;
            double x = (position - center)*axis;
            EXPECT_DOUBLE_EQ(A.Derivative(position, direction), dist_A.Derivative(x)*alpha);
            EXPECT_DOUBLE_EQ(B.Derivative(position, direction), dist_B.Derivative(x)*alpha);
            EXPECT_DOUBLE_EQ(C.Derivative(position, direction), dist_C.Derivative(x)*alpha);
        }
    }
}

TEST(Derivative, RadialPolynomial)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }

        std::function<double(double)> eval = [&](double x)->double {
            double res = params[n-1]*(n-1);
            for(int i=n-2; i>=1; --i) {
                res = res*x + params[i]*i;
            }
            return res;
        };

        Vector3D center = RandomVector();
        RadialAxis1D ax_A(center);

        PolynomialDistribution1D dist_A(params);

        auto A = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_A, dist_A);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D position = RandomVector();
            Vector3D direction = RandomDirection();
            double x = (position - center).magnitude();
            Vector3D r = (position - center);
            r.normalize();
            double dxdt = direction*r;
            EXPECT_DOUBLE_EQ(A.Derivative(position, direction), eval(x)*dxdt);
        }
    }
}

struct Quadrature
{
   //Abstract base class for elementary quadrature algorithms.
   int n; // Current level of refinement.

   virtual double next() = 0;
   //Returns the value of the integral at the nth stage of refinement.
   //The function next() must be defined in the derived class.
};

template<class T>
struct Trapzd: Quadrature
{
    double a, b, s; // Limits of integration and current value of integral.
    T &func;

    Trapzd() { };

    // func is function or functor to be integrated between limits: a and b
    Trapzd(T &funcc, const double aa, const double bb)
        : func(funcc), a(aa), b(bb)
    {
        n = 0;
    }

    // Returns the nth stage of refinement of the extended trapezoidal rule.
    // On the first call (n = 1), the routine returns the crudest estimate
    // of integral of f x / dx in [a,b]. Subsequent calls set n=2,3,... and
    // improve the accuracy by adding 2n - 2 additional interior points.
    double next()
    {
        double x, tnm, sum, del;
        int it, j;
        n++;

        if (n == 1)
        {
            return (s = 0.5 * (b-a) * (func(a) + func(b)));
        }
        else
        {
            for (it = 1, j = 1; j < n - 1; j++)
            {
                it <<= 1;
            }
            tnm = it;
            // This is the spacing of the points to be added.
            del = (b - a) / tnm;
            x = a + 0.5 * del;

            for (sum = 0.0,j = 0; j < it; j++, x += del)
            {
                sum += func(x);
            }
            // This replaces s by its refined value.
            s = 0.5 * (s + (b - a) * sum / tnm);
            return s;
        }
    }
};

template<class T>
double qtrap(T &func, const double a, const double b, const double eps = 1.0e-8)
{
    // Returns the integral of the function or functor func from a to b.
    // The constants EPS can be set to the desired fractional accuracy and
    // JMAX so that 2 to the power JMAX-1 is the maximum allowed number of
    // steps. integration is performed by the trapezoidal rule.

    const int JMAX = 25;
    double s, olds = 0.0; // Initial value of olds is arbitrary.

    Trapzd<T> t(func, a, b);

    for (int j = 0; j < JMAX; j++)
    {
        s = t.next();

        if (j > 5) // Avoid spurious early convergence.
        {
            if (abs(s - olds) < eps * abs(olds) || (s == 0.0 && olds == 0.0))
            {
                return s;
            }
        }
        olds = s;
    }
    EXPECT_TRUE(false) << "Too many steps in routine qtrap";
}

TEST(Integral, Axis_to_Distribution_connection)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*5+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*5);
        }

        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*1-0.5);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
        auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
        auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D p0 = RandomVector()*0.25;
            Vector3D p1 = RandomVector()*0.25;
            Vector3D direction = p1 - p0;
            double R = direction.magnitude();
            direction.normalize();
            double x0_A = ax_A.GetX(p0);
            double x0_B = ax_B.GetX(p0);
            double x1_A = ax_A.GetX(p1);
            double x1_B = ax_B.GetX(p1);
            std::function<double(double)> fA = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_A.Evaluate(ax_A.GetX(pos));
            };
            std::function<double(double)> fB = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_B.Evaluate(ax_A.GetX(pos));
            };
            std::function<double(double)> fC = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_C.Evaluate(ax_A.GetX(pos));
            };
            std::function<double(double)> fD = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_A.Evaluate(ax_B.GetX(pos));
            };
            std::function<double(double)> fE = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_B.Evaluate(ax_B.GetX(pos));
            };
            std::function<double(double)> fF = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_C.Evaluate(ax_B.GetX(pos));
            };
            double A_res = qtrap(fA, 0, R);
            double B_res = qtrap(fB, 0, R);
            double C_res = qtrap(fC, 0, R);
            double D_res = qtrap(fD, 0, R);
            double E_res = qtrap(fE, 0, R);
            double F_res = qtrap(fF, 0, R);
            EXPECT_NEAR(A.Integral(p0, p1), A_res, std::abs(A_res)*1e-8);
            EXPECT_NEAR(B.Integral(p0, p1), B_res, std::abs(B_res)*1e-8);
            EXPECT_NEAR(C.Integral(p0, p1), C_res, std::abs(C_res)*1e-8);
            EXPECT_NEAR(D.Integral(p0, p1), D_res, std::abs(D_res)*1e-8);
            EXPECT_NEAR(E.Integral(p0, p1), E_res, std::abs(E_res)*1e-4);
            EXPECT_NEAR(F.Integral(p0, p1), F_res, std::abs(F_res)*1e-4);

            EXPECT_NEAR(A.Integral(p0, direction, R), A_res, std::abs(A_res)*1e-8);
            EXPECT_NEAR(B.Integral(p0, direction, R), B_res, std::abs(B_res)*1e-8);
            EXPECT_NEAR(C.Integral(p0, direction, R), C_res, std::abs(C_res)*1e-8);
            EXPECT_NEAR(D.Integral(p0, direction, R), D_res, std::abs(D_res)*1e-8);
            EXPECT_NEAR(E.Integral(p0, direction, R), E_res, std::abs(E_res)*1e-4);
            EXPECT_NEAR(F.Integral(p0, direction, R), F_res, std::abs(F_res)*1e-4);

            EXPECT_NEAR(A.AntiDerivative(p1, direction)-A.AntiDerivative(p0, direction), A_res, std::abs(A_res)*1e-8);
            EXPECT_NEAR(B.AntiDerivative(p1, direction)-B.AntiDerivative(p0, direction), B_res, std::abs(B_res)*1e-8);
            EXPECT_NEAR(C.AntiDerivative(p1, direction)-C.AntiDerivative(p0, direction), C_res, std::abs(C_res)*1e-8);
            EXPECT_NEAR(D.AntiDerivative(p1, direction)-D.AntiDerivative(p0, direction), D_res, std::abs(D_res)*1e-8);
            EXPECT_NEAR(E.AntiDerivative(p1, direction)-E.AntiDerivative(p0, direction), E_res, std::abs(E_res)*1e-4);
            EXPECT_NEAR(F.AntiDerivative(p1, direction)-F.AntiDerivative(p0, direction), F_res, std::abs(F_res)*1e-4);
        }
    }
}

TEST(Integral, ConstantDistribution)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        double val = RandomDouble()*10;
        ConstantDistribution1D dist_A(val);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D p0 = RandomVector();
            Vector3D p1 = RandomVector();
            Vector3D direction = p1 - p0;
            double R = direction.magnitude();
            double res = R*val;
            double eps = std::abs(R*val)*1e-8;
            direction.normalize();
            EXPECT_NEAR(A.Integral(p0, p1), res, eps);
            EXPECT_NEAR(B.Integral(p0, p1), res, eps);

            EXPECT_NEAR(A.Integral(p0, direction, R), res, eps);
            EXPECT_NEAR(B.Integral(p0, direction, R), res, eps);

            EXPECT_NEAR(A.AntiDerivative(p1, direction)-A.AntiDerivative(p0, direction), res, eps);
            EXPECT_NEAR(B.AntiDerivative(p1, direction)-B.AntiDerivative(p0, direction), res, eps);
        }
    }
}

TEST(Integral, CartesianDistribution)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }
        Vector3D axis = RandomDirection();
        Vector3D center = RandomVector();
        CartesianAxis1D ax_A(axis, center);

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*4-2);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D p0 = RandomVector();
            Vector3D p1 = RandomVector();
            Vector3D direction = p1-p0;
            double R = direction.magnitude();
            direction.normalize();
            double alpha = axis*direction;
            double a = (p0 - center)*axis;
            double b = (p1 - center)*axis;

            double int_A = (dist_A.AntiDerivative(b) - dist_A.AntiDerivative(a))/alpha;
            double int_B = (dist_B.AntiDerivative(b) - dist_B.AntiDerivative(a))/alpha;
            double int_C = (dist_C.AntiDerivative(b) - dist_C.AntiDerivative(a))/alpha;
            EXPECT_NEAR(A.Integral(p0, p1), int_A, std::abs(int_A)*1e-8);
            EXPECT_NEAR(B.Integral(p0, p1), int_B, std::abs(int_B)*1e-8);
            EXPECT_NEAR(C.Integral(p0, p1), int_C, std::abs(int_C)*1e-8);

            EXPECT_NEAR(A.Integral(p0, direction, R), int_A, std::abs(int_A)*1e-8);
            EXPECT_NEAR(B.Integral(p0, direction, R), int_B, std::abs(int_B)*1e-8);
            EXPECT_NEAR(C.Integral(p0, direction, R), int_C, std::abs(int_C)*1e-8);

            EXPECT_NEAR(A.AntiDerivative(p1, direction)-A.AntiDerivative(p0, direction), int_A, std::abs(int_A)*1e-8);
            EXPECT_NEAR(B.AntiDerivative(p1, direction)-B.AntiDerivative(p0, direction), int_B, std::abs(int_B)*1e-8);
            EXPECT_NEAR(C.AntiDerivative(p1, direction)-C.AntiDerivative(p0, direction), int_C, std::abs(int_C)*1e-8);
        }
    }
}

double impact_parameter(Vector3D const & position, Vector3D const & direction, Vector3D const & center) {
    Vector3D cap = position - (position*direction)*direction;
    double b = (cap-center).magnitude();
    return b;
}

namespace {
double poly_0(double r, double b) {
    return sqrt(-std::pow(b,2) + std::pow(r,2));
}
double poly_1(double r, double b) {
    return (r*sqrt(-std::pow(b,2) + std::pow(r,2)))/2. + (std::pow(b,2)*atanh(r/sqrt(-std::pow(b,2) + std::pow(r,2))))/2.;
}
double poly_2(double r, double b) {
    return (sqrt(-std::pow(b,2) + std::pow(r,2))*(2*std::pow(b,2) + std::pow(r,2)))/3.;
}
double poly_3(double r, double b) {
    return (r*sqrt(-std::pow(b,2) + std::pow(r,2))*(3*std::pow(b,2) + 2*std::pow(r,2)) + 3*std::pow(b,4)*atanh(r/sqrt(-std::pow(b,2) + std::pow(r,2))))/8.;
}
double poly_4(double r, double b) {
    return (sqrt(-std::pow(b,2) + std::pow(r,2))*(8*std::pow(b,4) + 4*std::pow(b,2)*std::pow(r,2) + 3*std::pow(r,4)))/15.;
}
double poly_5(double r, double b) {
    return (r*sqrt(-std::pow(b,2) + std::pow(r,2))*(15*std::pow(b,4) + 10*std::pow(b,2)*std::pow(r,2) + 8*std::pow(r,4)) + 15*std::pow(b,6)*atanh(r/sqrt(-std::pow(b,2) + std::pow(r,2))))/48.;
}
double poly_6(double r, double b) {
    return (sqrt(-std::pow(b,2) + std::pow(r,2))*(16*std::pow(b,6) + 8*std::pow(b,4)*std::pow(r,2) + 6*std::pow(b,2)*std::pow(r,4) + 5*std::pow(r,6)))/35.;
}
double poly_7(double r, double b) {
    return (r*sqrt(-std::pow(b,2) + std::pow(r,2))*(105*std::pow(b,6) + 70*std::pow(b,4)*std::pow(r,2) + 56*std::pow(b,2)*std::pow(r,4) + 48*std::pow(r,6)) + 105*std::pow(b,8)*atanh(r/sqrt(-std::pow(b,2) + std::pow(r,2))))/384.;
}
double poly_8(double r, double b) {
    return (sqrt(-std::pow(b,2) + std::pow(r,2))*(128*std::pow(b,8) + 64*std::pow(b,6)*std::pow(r,2) + 48*std::pow(b,4)*std::pow(r,4) + 40*std::pow(b,2)*std::pow(r,6) + 35*std::pow(r,8)))/315.;
}
double poly_9(double r, double b) {
    return (r*sqrt(-std::pow(b,2) + std::pow(r,2))*(315*std::pow(b,8) + 210*std::pow(b,6)*std::pow(r,2) + 168*std::pow(b,4)*std::pow(r,4) + 144*std::pow(b,2)*std::pow(r,6) + 128*std::pow(r,8)) + 315*std::pow(b,10)*atanh(r/sqrt(-std::pow(b,2) + std::pow(r,2))))/1280.;
}
double poly_10(double r, double b) {
    return (sqrt(-std::pow(b,2) + std::pow(r,2))*(256*std::pow(b,10) + 128*std::pow(b,8)*std::pow(r,2) + 96*std::pow(b,6)*std::pow(r,4) + 80*std::pow(b,4)*std::pow(r,6) + 70*std::pow(b,2)*std::pow(r,8) + 63*std::pow(r,10)))/693.;
}
double poly_11(double r, double b) {
    return (r*sqrt(-std::pow(b,2) + std::pow(r,2))*(3465*std::pow(b,10) + 2310*std::pow(b,8)*std::pow(r,2) + 1848*std::pow(b,6)*std::pow(r,4) + 1584*std::pow(b,4)*std::pow(r,6) + 1408*std::pow(b,2)*std::pow(r,8) + 1280*std::pow(r,10)) + 3465*std::pow(b,12)*atanh(r/sqrt(-std::pow(b,2) + std::pow(r,2))))/15360.;
}
double poly_12(double r, double b) {
    return (sqrt(-std::pow(b,2) + std::pow(r,2))*(1024*std::pow(b,12) + 512*std::pow(b,10)*std::pow(r,2) + 384*std::pow(b,8)*std::pow(r,4) + 320*std::pow(b,6)*std::pow(r,6) + 280*std::pow(b,4)*std::pow(r,8) + 252*std::pow(b,2)*std::pow(r,10) + 231*std::pow(r,12)))/3003.;
}

template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}

double poly_int(std::vector<double> const & params, Vector3D const & p0, Vector3D const & p1, Vector3D const & center) {
    Vector3D direction = p1-p0;
    double R = direction.magnitude();
    direction.normalize();

    Vector3D cap = ((p0-center) - ((p0-center)*direction)*direction + (p1-center) - ((p1-center)*direction)*direction)*0.5 + center;
    double b = (cap-center).magnitude();

    EXPECT_NEAR((cap-center)*(p0-cap), 0.0, 1e-8);
    EXPECT_NEAR((cap-center)*(p1-cap), 0.0, 1e-8);
    EXPECT_NEAR((cap-center)*(p1-p0), 0.0, 1e-8);

    double t0 = (p0 - cap)*direction;
    double t1 = (p1 - cap)*direction;
    assert(t1 >= t0);

    double r0 = sqrt(t0*t0 + b*b);
    double r1 = sqrt(t1*t1 + b*b);

    EXPECT_NEAR(r0, (p0-center).magnitude(), (p0-center).magnitude()*1e-8);
    EXPECT_NEAR(r1, (p1-center).magnitude(), (p1-center).magnitude()*1e-8);
    EXPECT_DOUBLE_EQ(b, (cap-center).magnitude());
    EXPECT_TRUE(t1 >= t0);

    std::vector<std::pair<double,double>> r_points;

    if(signum(t0) == signum(t1)) {
        if(t0 > 0) {
            r_points.push_back({r1, 1});
            r_points.push_back({r0, -1});
        }
        else {
            r_points.push_back({r1, -1});
            r_points.push_back({r0, 1});
        }
    }
    else {
        if(t0 == 0) {
            r_points.push_back({r1, 1});
            r_points.push_back({r0, -1});
        }
        else if(t1 == 0) {
            r_points.push_back({r1, -1});
            r_points.push_back({r0, 1});
        }
        else {
            assert(t0 < 0 and t1 > 0);
            r_points.push_back({r0, 1});
            r_points.push_back({r1, 1});
            r_points.push_back({b, -2});
        }
    }

    std::vector<std::function<double(double,double)>> poly_ints{poly_0, poly_1, poly_2, poly_3, poly_4, poly_5, poly_6, poly_7, poly_8, poly_9, poly_10, poly_11, poly_12};

    assert(params.size() <= poly_ints.size());

    std::vector<double> terms;
    for(unsigned int i=0; i<params.size(); ++i) {
        for(unsigned int j=0; j<r_points.size(); ++j) {
            terms.push_back(params[i]*r_points[j].second*poly_ints[i](r_points[j].first, b));
        }
    }
    assert(terms.size() == 2*params.size() or terms.size() == 3*params.size());

    auto print = [&]() {
        std::stringstream s;
        s << "R " << R << "\n";
        s << "t0 " << t0 << "\n";
        s << "t1 " << t1 << "\n";
        s << "r0 " << r0 << "\n";
        s << "r1 " << r1 << "\n";
        s << "params ";
        for (auto a : params) {
            s << a << ' ';
        }
        s << "\n";
        s << "terms ";
        for (auto a : terms) {
            s << a << ' ';
        }
        return s.str();
    };

    if(r_points.size() == 2) {
        double res = terms[0] + terms[1];
        EXPECT_NEAR(res, R*params[0], std::abs(R*params[0])*1e-8) << print();
    }
    else {
        double res = terms[0] + terms[1] + terms[2];
        EXPECT_NEAR(res, R*params[0], std::abs(R*params[0])*1e-8) << print();
    }

    std::function<bool(double, double)> comp = [](double a, double b)->bool {
        return std::abs(a) < std::abs(b);
    };

    std::sort(terms.begin(), terms.end(), comp);

    double res = 0.0;
    for(unsigned int i=0; i<terms.size(); ++i) {
        res += terms[i];
    }

    return res;
}

}

TEST(Integral, RadialPolynomial)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        //unsigned int n = (int)(RandomDouble()*10+2);
        unsigned int n = 1;
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }

        std::function<double(double)> eval = [&](double x)->double {
            double res = params[n-1]*(n-1);
            for(int i=n-2; i>=1; --i) {
                res = res*x + params[i]*i;
            }
            return res;
        };

        Vector3D center = RandomVector();
        RadialAxis1D ax_A(center);

        PolynomialDistribution1D dist_A(params);

        auto A = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_A, dist_A);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D p0 = RandomVector();
            Vector3D p1 = RandomVector();
            Vector3D direction = p1-p0;
            double R = direction.magnitude();
            direction.normalize();
            double res = poly_int(params, p0, p1, center);
            EXPECT_NEAR(A.Integral(p0, p1), res, std::abs(res)*1e-8);
            EXPECT_NEAR(A.Integral(p0, direction, R), res, std::abs(res)*1e-8);
            EXPECT_NEAR(A.AntiDerivative(p1, direction) - A.AntiDerivative(p0, direction), res, std::abs(res)*1e-8);
        }
    }
}

TEST(InverseIntegral, Axis_to_Distribution_connection)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*5+2);
        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*5/i);
        }

        CartesianAxis1D ax_A(RandomDirection(), RandomVector());
        RadialAxis1D ax_B(RandomVector());

        ConstantDistribution1D dist_A(RandomDouble()*10);
        PolynomialDistribution1D dist_B(params);
        ExponentialDistribution1D dist_C(RandomDouble()*1-0.5);

        auto A = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(ax_A, dist_A);
        auto B = DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D>(ax_A, dist_B);
        auto C = DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, dist_C);
        auto D = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(ax_B, dist_A);
        auto E = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(ax_B, dist_B);
        auto F = DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D>(ax_B, dist_C);

        for(unsigned int i=0; i<N_RAND; ++i) {
            Vector3D p0 = RandomVector()*0.25;
            Vector3D p1 = RandomVector()*0.25;
            Vector3D direction = p1 - p0;
            double R = direction.magnitude();
            double int_fraction = RandomDouble()*(1.5-0.25);
            double A_int = A.Integral(p0, p1)*int_fraction;
            double B_int = B.Integral(p0, p1)*int_fraction;
            double C_int = C.Integral(p0, p1)*int_fraction;
            double D_int = D.Integral(p0, p1)*int_fraction;
            double E_int = E.Integral(p0, p1)*int_fraction;
            double F_int = F.Integral(p0, p1)*int_fraction;
            direction.normalize();
            double x0_A = ax_A.GetX(p0);
            double x0_B = ax_B.GetX(p0);
            double x1_A = ax_A.GetX(p1);
            double x1_B = ax_B.GetX(p1);
            std::function<double(double)> fA = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_A.Evaluate(ax_A.GetX(pos));
            };
            std::function<double(double)> fB = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_B.Evaluate(ax_A.GetX(pos));
            };
            std::function<double(double)> fC = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_C.Evaluate(ax_A.GetX(pos));
            };
            std::function<double(double)> fD = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_A.Evaluate(ax_B.GetX(pos));
            };
            std::function<double(double)> fE = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_B.Evaluate(ax_B.GetX(pos));
            };
            std::function<double(double)> fF = [&](double x)->double {
                Vector3D pos = p0 + direction*x;
                return dist_C.Evaluate(ax_B.GetX(pos));
            };

            std::function<double(double)> FA = [&](double x)->double {
                return qtrap(fA, 0, x) - A_int;
            };
            std::function<double(double)> FB = [&](double x)->double {
                return qtrap(fB, 0, x) - B_int;
            };
            std::function<double(double)> FC = [&](double x)->double {
                return qtrap(fC, 0, x) - C_int;
            };
            std::function<double(double)> FD = [&](double x)->double {
                return qtrap(fD, 0, x) - D_int;
            };
            std::function<double(double)> FE = [&](double x)->double {
                return qtrap(fE, 0, x) - E_int;
            };
            std::function<double(double)> FF = [&](double x)->double {
                return qtrap(fF, 0, x) - F_int;
            };

            double A_res;
            try {
                A_res = NewtonRaphson(FA, fA, 0, R, R/2);
            } catch(MathException& e) {
                A_res = -1;
            }
            double B_res;
            try {
                B_res = NewtonRaphson(FB, fB, 0, R, R/2);
            } catch(MathException& e) {
                B_res = -1;
            }
            double C_res;
            try {
                C_res = NewtonRaphson(FC, fC, 0, R, R/2);
            } catch(MathException& e) {
                C_res = -1;
            }
            double D_res;
            try {
                D_res = NewtonRaphson(FD, fD, 0, R, R/2);
            } catch(MathException& e) {
                D_res = -1;
            }
            double E_res;
            try {
                E_res = NewtonRaphson(FE, fE, 0, R, R/2);
            } catch(MathException& e) {
                E_res = -1;
            }
            double F_res;
            try {
                F_res = NewtonRaphson(FF, fF, 0, R, R/2);
            } catch(MathException& e) {
                F_res = -1;
            }

            EXPECT_NEAR(A.InverseIntegral(p0, direction, A_int, R), A_res, std::abs(A_res)*1e-8);
            EXPECT_NEAR(B.InverseIntegral(p0, direction, B_int, R), B_res, std::abs(B_res)*1e-7);
            EXPECT_NEAR(C.InverseIntegral(p0, direction, C_int, R), C_res, std::abs(C_res)*1e-8);
            EXPECT_NEAR(D.InverseIntegral(p0, direction, D_int, R), D_res, std::abs(D_res)*1e-8);
            EXPECT_NEAR(E.InverseIntegral(p0, direction, E_int, R), E_res, std::abs(E_res)*1e-8);
            EXPECT_NEAR(F.InverseIntegral(p0, direction, F_int, R), F_res, std::abs(F_res)*1e-8);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

