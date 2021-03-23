
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

TEST(Evaluate, Constant)
{
    unsigned int N_RAND = 100;
    double default_val = 1e-25;
    double val = 2.0;
    ConstantDistribution1D A;
    ConstantDistribution1D B(val);
    ConstantDistribution1D C(B);

    for(unsigned int i=0; i<N_RAND; ++i) {
        double x = RandomDouble()*40-20;
        EXPECT_DOUBLE_EQ(A.Evaluate(x), default_val);
        EXPECT_DOUBLE_EQ(B.Evaluate(x), val);
        EXPECT_DOUBLE_EQ(C.Evaluate(x), val);
    }
}

TEST(Evaluate, Polynomial)
{
    unsigned int N_RAND = 100;
    double default_val = 0;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);

        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }
        ASSERT_TRUE(params.size() == n);
        assert(params.size() == n);

        Polynom poly(params);

        PolynomialDistribution1D A({});
        PolynomialDistribution1D B(poly);
        PolynomialDistribution1D C(params);
        PolynomialDistribution1D D(A);
        PolynomialDistribution1D E(B);
        PolynomialDistribution1D F(C);

        std::function<double(double)> eval = [&](double x)->double {
            double res = params[n-1];
            for(int i=n-2; i>=0; --i) {
                res = res*x + params[i];
            }
            return res;
        };

        for(unsigned int i=0; i<N_RAND; ++i) {
            double x = RandomDouble()*40-20;
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.Evaluate(x), default_val);
            EXPECT_DOUBLE_EQ(B.Evaluate(x), res);
            EXPECT_DOUBLE_EQ(C.Evaluate(x), res);
            EXPECT_DOUBLE_EQ(D.Evaluate(x), default_val);
            EXPECT_DOUBLE_EQ(E.Evaluate(x), res);
            EXPECT_DOUBLE_EQ(F.Evaluate(x), res);
        }
    }
}

TEST(Evaluate, Exponential)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        double sigma = RandomDouble()*6-3;

        ExponentialDistribution1D A(sigma);
        ExponentialDistribution1D B(A);

        std::function<double(double)> eval = [&](double x)->double {
            return exp(sigma*x);
        };

        for(unsigned int i=0; i<N_RAND; ++i) {
            double x = RandomDouble()*10-5;
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.Evaluate(x), res);
            EXPECT_DOUBLE_EQ(B.Evaluate(x), res);
        }
    }
}

TEST(Derivative, Constant)
{
    unsigned int N_RAND = 100;
    double default_val = 1e-25;
    double val = 2.0;
    ConstantDistribution1D A;
    ConstantDistribution1D B(val);
    ConstantDistribution1D C(B);

    for(unsigned int i=0; i<N_RAND; ++i) {
        double x = RandomDouble()*40-20;
        EXPECT_DOUBLE_EQ(A.Derivative(x), 0.0);
        EXPECT_DOUBLE_EQ(B.Derivative(x), 0.0);
        EXPECT_DOUBLE_EQ(C.Derivative(x), 0.0);
    }
}

TEST(Derivative, Polynomial)
{
    unsigned int N_RAND = 100;
    double default_val = 0;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);

        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }
        ASSERT_TRUE(params.size() == n);
        assert(params.size() == n);

        Polynom poly(params);

        PolynomialDistribution1D A({});
        PolynomialDistribution1D B(poly);
        PolynomialDistribution1D C(params);
        PolynomialDistribution1D D(A);
        PolynomialDistribution1D E(B);
        PolynomialDistribution1D F(C);

        std::function<double(double)> eval = [&](double x)->double {
            double res = params[n-1]*(n-1);
            for(int i=n-2; i>=1; --i) {
                res = res*x + params[i]*i;
            }
            return res;
        };

        for(unsigned int i=0; i<N_RAND; ++i) {
            double x = RandomDouble()*40-20;
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.Derivative(x), default_val);
            EXPECT_DOUBLE_EQ(B.Derivative(x), res);
            EXPECT_DOUBLE_EQ(C.Derivative(x), res);
            EXPECT_DOUBLE_EQ(D.Derivative(x), default_val);
            EXPECT_DOUBLE_EQ(E.Derivative(x), res);
            EXPECT_DOUBLE_EQ(F.Derivative(x), res);
        }
    }
}

TEST(Derivative, Exponential)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        double sigma = RandomDouble()*6-3;

        ExponentialDistribution1D A(sigma);
        ExponentialDistribution1D B(A);

        std::function<double(double)> eval = [&](double x)->double {
            return sigma*exp(sigma*x);
        };

        for(unsigned int i=0; i<N_RAND; ++i) {
            double x = RandomDouble()*10-5;
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.Derivative(x), res);
            EXPECT_DOUBLE_EQ(B.Derivative(x), res);
        }
    }
}

TEST(AntiDerivative, Constant)
{
    unsigned int N_RAND = 100;
    double default_val = 1e-25;
    double val = 2.0;
    ConstantDistribution1D A;
    ConstantDistribution1D B(val);
    ConstantDistribution1D C(B);

    for(unsigned int i=0; i<N_RAND; ++i) {
        double x = RandomDouble()*40-20;
        EXPECT_DOUBLE_EQ(A.AntiDerivative(x), default_val*x);
        EXPECT_DOUBLE_EQ(B.AntiDerivative(x), val*x);
        EXPECT_DOUBLE_EQ(C.AntiDerivative(x), val*x);
    }
}

TEST(AntiDerivative, Polynomial)
{
    unsigned int N_RAND = 100;
    double default_val = 0;

    for(unsigned int p=0; p<N_RAND; ++p) {
        unsigned int n = (int)(RandomDouble()*10+2);

        std::vector<double> params;
        for(unsigned int i=1; i<(n+1); ++i) {
            params.push_back(RandomDouble()*20-10);
        }
        ASSERT_TRUE(params.size() == n);
        assert(params.size() == n);

        Polynom poly(params);

        PolynomialDistribution1D A({});
        PolynomialDistribution1D B(poly);
        PolynomialDistribution1D C(params);
        PolynomialDistribution1D D(A);
        PolynomialDistribution1D E(B);
        PolynomialDistribution1D F(C);

        std::function<double(double)> eval = [&](double x)->double {
            double res = params[n-1]/(n);
            for(int i=n-2; i>=0; --i) {
                res = res*x + params[i]/(i+1);
            }
            res *= x;
            return res;
        };

        for(unsigned int i=0; i<N_RAND; ++i) {
            double x = RandomDouble()*40-20;
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.AntiDerivative(x), default_val);
            EXPECT_DOUBLE_EQ(B.AntiDerivative(x), res);
            EXPECT_DOUBLE_EQ(C.AntiDerivative(x), res);
            EXPECT_DOUBLE_EQ(D.AntiDerivative(x), default_val);
            EXPECT_DOUBLE_EQ(E.AntiDerivative(x), res);
            EXPECT_DOUBLE_EQ(F.AntiDerivative(x), res);
        }
    }
}

TEST(AntiDerivative, Exponential)
{
    unsigned int N_RAND = 100;

    for(unsigned int p=0; p<N_RAND; ++p) {
        double sigma = RandomDouble()*6-3;

        ExponentialDistribution1D A(sigma);
        ExponentialDistribution1D B(A);

        std::function<double(double)> eval = [&](double x)->double {
            return exp(sigma*x)/sigma;
        };

        for(unsigned int i=0; i<N_RAND; ++i) {
            double x = RandomDouble()*10-5;
            double res = eval(x);
            EXPECT_DOUBLE_EQ(A.AntiDerivative(x), res);
            EXPECT_DOUBLE_EQ(B.AntiDerivative(x), res);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

