
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>

#include <gtest/gtest.h>

#include "earthmodel-service/Polynomial.h"

using namespace earthmodel;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Comparison, Comparison_equal)
{
    Polynom A({});
    Polynom B({});
    EXPECT_TRUE(A == B);

    std::vector<double> params_a{0,1,2,3,4,5};
    Polynom C(params_a);
    Polynom D(params_a);
    EXPECT_TRUE(C == D);

    Polynom E(A);
    Polynom F(C);
    EXPECT_TRUE(A == E);
    EXPECT_TRUE(C == F);

    std::vector<double> params_b{0,1,2,3,4,5,6};
    Polynom G(params_b);

    std::vector<double> params_c{0,1,2,3,4,-5};
    Polynom H(params_c);

    EXPECT_FALSE(A == C);
    EXPECT_FALSE(A == G);
    EXPECT_FALSE(A == H);
    EXPECT_FALSE(C == G);
    EXPECT_FALSE(C == H);
}

TEST(Comparison, Comparison_not_equal)
{
    Polynom A({});
    Polynom B({});
    EXPECT_FALSE(A != B);

    std::vector<double> params_a{0,1,2,3,4,5};
    Polynom C(params_a);
    Polynom D(params_a);
    EXPECT_FALSE(C != D);

    Polynom E(A);
    Polynom F(C);
    EXPECT_FALSE(A != E);
    EXPECT_FALSE(C != F);

    std::vector<double> params_b{0,1,2,3,4,5,6};
    Polynom G(params_b);

    std::vector<double> params_c{0,1,2,3,4,-5};
    Polynom H(params_c);

    EXPECT_TRUE(A != C);
    EXPECT_TRUE(A != G);
    EXPECT_TRUE(A != H);
    EXPECT_TRUE(C != G);
    EXPECT_TRUE(C != H);
}

TEST(Assignment, CopyConstructor)
{
    Polynom A({1, 2, 3});
    Polynom B(A);
    EXPECT_TRUE(A == B);
    Polynom C({2, 3, 4});
    B = C;
    EXPECT_TRUE(C == B);
}

TEST(Assignment, Operator)
{
    Polynom A({1, 2, 3});
    Polynom B({});
    EXPECT_TRUE(A != B);
    B = A;
    EXPECT_TRUE(A == B);
}

TEST(Get, GetCoefficient)
{
    std::vector<double> params{0,1,2,3,4,5,6,7,8,9};
    Polynom A(params);
    std::vector<double> coeff = A.GetCoefficient();
    EXPECT_TRUE(params.size() == coeff.size());
    for(unsigned int i=0; i<params.size(); ++i) {
        EXPECT_TRUE(params[i] == coeff[i]);
    }
}

TEST(Get, GetDerivative)
{
    std::vector<double> params{9,8,7,6,5,4,3,2,1,0};
    Polynom A(params);
    Polynom B = A.GetDerivative();
    std::vector<double> d_params = B.GetCoefficient();
    EXPECT_EQ(d_params.size(), params.size()-1);
    for(unsigned int i=1; i<params.size(); ++i) {
        double new_coeff = params[i]*i;
        EXPECT_DOUBLE_EQ(d_params[i-1], new_coeff);
    }
}

TEST(Get, GetAntiderivative)
{
    std::vector<double> params{9,8,7,6,5,4,3,2,1,0};
    Polynom A(params);
    double c = -20;
    Polynom B = A.GetAntiderivative(c);
    std::vector<double> d_params = B.GetCoefficient();
    EXPECT_EQ(d_params.size(), params.size()+1);
    EXPECT_DOUBLE_EQ(d_params[0], c);
    for(unsigned int i=0; i<params.size(); ++i) {
        double new_coeff = params[i] / (i+1);
        EXPECT_DOUBLE_EQ(d_params[i+1], new_coeff);
    }
}

TEST(PolynomialValue, Evaluate) {
	unsigned int N_RAND = 100;
    std::vector<double> zero_params{0,0,0,0,0,0,0,0,0,0};
    std::vector<double> selector;
    for(unsigned int i=0; i<zero_params.size(); ++i) {
        selector = zero_params;
        selector[i] = 1.0;
        Polynom poly(selector);
        for(unsigned int j=0; j<N_RAND; ++j) {
            double x = RandomDouble()*4-2;
            EXPECT_DOUBLE_EQ(poly.evaluate(x), std::pow(x, i)) << "Check single power";
		}
    }
    std::vector<double> params;
    for(unsigned int i=0; i<N_RAND; ++i) {
        params = {RandomDouble()*4-2, RandomDouble()*4-2, RandomDouble()*4-2, RandomDouble()*4-2, RandomDouble()*4-2};
        Polynom poly(params);
        for(unsigned int j=0; j<N_RAND; ++j) {
            double x = RandomDouble()*4-2;
            double actual = poly.evaluate(x);
            double expect = params[0] + x*(params[1] + x*(params[2] + x*(params[3] + x*(params[4]))));
            EXPECT_DOUBLE_EQ(actual, expect) << "Check 5 powers";
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

