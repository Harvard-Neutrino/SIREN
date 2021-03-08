
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>
#include <gtest/gtest.h>

#include <earthmodel-service/DensityDist.h>
#include <earthmodel-service/Vector3D.h>
#include <earthmodel-service/Polynomial.h>

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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
