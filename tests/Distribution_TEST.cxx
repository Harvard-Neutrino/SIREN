
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

TEST(Comparison, Comparison_equal)
{
    Distribution1D* A = new ConstantDistribution1D();
    Distribution1D* B = new ConstantDistribution1D();
    EXPECT_TRUE(*A == *B);
    delete A;
    delete B;

    A = new PolynomialDistribution1D({});
    B = new PolynomialDistribution1D({});
    EXPECT_TRUE(*A == *B);
    delete A;
    delete B;

    A = new ExponentialDistribution1D(1.0);
    B = new ExponentialDistribution1D(1.0);
    EXPECT_TRUE(*A == *B);
    delete A;
    delete B;

    Distribution1D* C;
    A = new ConstantDistribution1D();
    B = new PolynomialDistribution1D({});
    C = new ExponentialDistribution1D(1.0);
    EXPECT_FALSE(*A == *B);
    EXPECT_FALSE(*B == *A);
    EXPECT_FALSE(*B == *C);
    EXPECT_FALSE(*C == *B);
    EXPECT_FALSE(*C == *A);
    EXPECT_FALSE(*A == *C);
    delete A;
    delete B;
    delete C;

    A = new ConstantDistribution1D(1.0);
    B = new ConstantDistribution1D(2.0);
    C = new ConstantDistribution1D(-1.0);
    EXPECT_FALSE(*A == *B);
    EXPECT_FALSE(*B == *A);
    EXPECT_FALSE(*B == *C);
    EXPECT_FALSE(*C == *B);
    EXPECT_FALSE(*C == *A);
    EXPECT_FALSE(*A == *C);
    delete A;
    delete B;
    delete C;

    A = new PolynomialDistribution1D({1.0});
    B = new PolynomialDistribution1D({2.0});
    C = new PolynomialDistribution1D({-1.0});
    EXPECT_FALSE(*A == *B);
    EXPECT_FALSE(*B == *A);
    EXPECT_FALSE(*B == *C);
    EXPECT_FALSE(*C == *B);
    EXPECT_FALSE(*C == *A);
    EXPECT_FALSE(*A == *C);
    delete A;
    delete B;
    delete C;

    A = new PolynomialDistribution1D({});
    B = new PolynomialDistribution1D({1.0});
    C = new PolynomialDistribution1D({1.0,1.0});
    EXPECT_FALSE(*A == *B);
    EXPECT_FALSE(*B == *A);
    EXPECT_FALSE(*B == *C);
    EXPECT_FALSE(*C == *B);
    EXPECT_FALSE(*C == *A);
    EXPECT_FALSE(*A == *C);
    delete A;
    delete B;
    delete C;

    A = new ExponentialDistribution1D(1.0);
    B = new ExponentialDistribution1D(2.0);
    C = new ExponentialDistribution1D(-1.0);
    EXPECT_FALSE(*A == *B);
    EXPECT_FALSE(*B == *A);
    EXPECT_FALSE(*B == *C);
    EXPECT_FALSE(*C == *B);
    EXPECT_FALSE(*C == *A);
    EXPECT_FALSE(*A == *C);
    delete A;
    delete B;
    delete C;
}

TEST(Comparison, Comparison_not_equal)
{
    Distribution1D* A = new ConstantDistribution1D();
    Distribution1D* B = new ConstantDistribution1D();
    EXPECT_FALSE(*A != *B);
    delete A;
    delete B;

    A = new PolynomialDistribution1D({});
    B = new PolynomialDistribution1D({});
    EXPECT_FALSE(*A != *B);
    delete A;
    delete B;

    A = new ExponentialDistribution1D(1.0);
    B = new ExponentialDistribution1D(1.0);
    EXPECT_FALSE(*A != *B);
    delete A;
    delete B;

    Distribution1D* C;
    A = new ConstantDistribution1D();
    B = new PolynomialDistribution1D({});
    C = new ExponentialDistribution1D(1.0);
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *A);
    EXPECT_TRUE(*B != *C);
    EXPECT_TRUE(*C != *B);
    EXPECT_TRUE(*C != *A);
    EXPECT_TRUE(*A != *C);
    delete A;
    delete B;
    delete C;

    A = new ConstantDistribution1D(1.0);
    B = new ConstantDistribution1D(2.0);
    C = new ConstantDistribution1D(-1.0);
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *A);
    EXPECT_TRUE(*B != *C);
    EXPECT_TRUE(*C != *B);
    EXPECT_TRUE(*C != *A);
    EXPECT_TRUE(*A != *C);
    delete A;
    delete B;
    delete C;

    A = new PolynomialDistribution1D({1.0});
    B = new PolynomialDistribution1D({2.0});
    C = new PolynomialDistribution1D({-1.0});
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *A);
    EXPECT_TRUE(*B != *C);
    EXPECT_TRUE(*C != *B);
    EXPECT_TRUE(*C != *A);
    EXPECT_TRUE(*A != *C);
    delete A;
    delete B;
    delete C;

    A = new PolynomialDistribution1D({});
    B = new PolynomialDistribution1D({1.0});
    C = new PolynomialDistribution1D({1.0,1.0});
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *A);
    EXPECT_TRUE(*B != *C);
    EXPECT_TRUE(*C != *B);
    EXPECT_TRUE(*C != *A);
    EXPECT_TRUE(*A != *C);
    delete A;
    delete B;
    delete C;

    A = new ExponentialDistribution1D(1.0);
    B = new ExponentialDistribution1D(2.0);
    C = new ExponentialDistribution1D(-1.0);
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *A);
    EXPECT_TRUE(*B != *C);
    EXPECT_TRUE(*C != *B);
    EXPECT_TRUE(*C != *A);
    EXPECT_TRUE(*A != *C);
    delete A;
    delete B;
    delete C;
}

TEST(Assignment, Copyconstructor)
{
    ConstantDistribution1D cA(1.0);
    ConstantDistribution1D cB = cA;
    ConstantDistribution1D cC(cA);
    EXPECT_TRUE(cA == cB);
    EXPECT_TRUE(cA == cC);

    PolynomialDistribution1D pA({1.0});
    PolynomialDistribution1D pB = pA;
    PolynomialDistribution1D pC(pA);
    EXPECT_TRUE(pA == pB);
    EXPECT_TRUE(pA == pC);

    ExponentialDistribution1D eA(1.0);
    ExponentialDistribution1D eB = eA;
    ExponentialDistribution1D eC(eA);
    EXPECT_TRUE(eA == eB);
    EXPECT_TRUE(eA == eC);
}

TEST(Copy, clone)
{
    ConstantDistribution1D A(1.0);
    PolynomialDistribution1D B({1.0});
    ExponentialDistribution1D C(1.0);

    Distribution1D* Ap = A.clone();
    Distribution1D* Bp = B.clone();
    Distribution1D* Cp = C.clone();

    EXPECT_TRUE(*Ap == A);
    EXPECT_TRUE(*Bp == B);
    EXPECT_TRUE(*Cp == C);

    EXPECT_TRUE((A == B) == (*Ap == *Bp));
    EXPECT_TRUE((A == C) == (*Ap == *Cp));
    EXPECT_TRUE((B == C) == (*Bp == *Cp));

    delete Ap;
    delete Bp;
    delete Cp;
}

TEST(Copy, create)
{
    ConstantDistribution1D A(1.0);
    PolynomialDistribution1D B({1.0});
    ExponentialDistribution1D C(1.0);

    std::shared_ptr<const Distribution1D> Ap = A.create();
    std::shared_ptr<const Distribution1D> Bp = B.create();
    std::shared_ptr<const Distribution1D> Cp = C.create();

    EXPECT_TRUE(A == *Ap);
    EXPECT_TRUE(B == *Bp);
    EXPECT_TRUE(C == *Cp);

    EXPECT_TRUE((A == B) == (*Ap == *Bp));
    EXPECT_TRUE((A == C) == (*Ap == *Cp));
    EXPECT_TRUE((B == C) == (*Bp == *Cp));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
