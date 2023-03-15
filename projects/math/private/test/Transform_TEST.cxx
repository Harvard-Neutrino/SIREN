
#include <cmath>
#include <math.h>
#include <random>
#include <iostream>

#include <gtest/gtest.h>

#include "LeptonInjector/math/Interpolation.h"

using namespace LI::math;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

double f1(double x){
    return pow(x, 3) - 4;
}

double f2(double x){
    return std::sin(x);
}

TEST(IdentityTransform, Constructor) {
    ASSERT_NO_THROW(IdentityTransform<double>());
}

TEST(IdentityTransform, Function) {
    IdentityTransform<double> transform;

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        EXPECT_TRUE(x == transform.Function(x));
    }
}

TEST(IdentityTransform, Inverse) {
    IdentityTransform<double> transform;

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        EXPECT_TRUE(x == transform.Function(x));
    }
}

TEST(GenericTransform, Constructor) {
    std::function<double(double)> function = [](double)->double{return 0.0;};
    std::function<double(double)> inverse = [](double)->double{return 0.0;};
    ASSERT_NO_THROW(GenericTransform<double>(function, inverse));
}

TEST(GenericTransform, Function) {
    bool is_called = false;
    double expected = 0;
    double received = 0;

    std::function<double(double)> function = [&](double x)->double{
        is_called = true;
        received = x;
        return expected;
    };
    std::function<double(double)> inverse = [&](double x)->double{
        return 0.0;
    };

    GenericTransform<double> transform(function, inverse);

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        expected = RandomDouble();
        is_called = false;
        EXPECT_TRUE(expected == transform.Function(x));
        EXPECT_TRUE(is_called);
        EXPECT_TRUE(received == x);
    }
}

TEST(GenericTransform, Inverse) {
    bool is_called = false;
    double expected = 0;
    double received = 0;

    std::function<double(double)> function = [&](double x)->double{
        return 0.0;
    };
    std::function<double(double)> inverse = [&](double x)->double{
        is_called = true;
        received = x;
        return expected;
    };

    GenericTransform<double> transform(function, inverse);

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        expected = RandomDouble();
        is_called = false;
        EXPECT_TRUE(expected == transform.Inverse(x));
        EXPECT_TRUE(is_called);
        EXPECT_TRUE(received == x);
    }
}

TEST(LogTransform, Constructor) {
    ASSERT_NO_THROW(LogTransform<double>());
}

TEST(LogTransform, Function) {
    LogTransform<double> transform;

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        EXPECT_TRUE(log(x) == transform.Function(x));
        EXPECT_TRUE(log(x * 1e8) == transform.Function(x * 1e8));
        EXPECT_TRUE(std::isnan(transform.Function(-x)));
        EXPECT_TRUE(std::isnan(transform.Function(-x * 1e8)));
    }
}

TEST(LogTransform, Inverse) {
    LogTransform<double> transform;

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        EXPECT_TRUE(exp(x) == transform.Inverse(x));
        EXPECT_TRUE(exp(-x) == transform.Inverse(-x));
        EXPECT_TRUE(exp(x * 1e8) == transform.Inverse(x * 1e8));
        EXPECT_TRUE(exp(-x * 1e8) == transform.Inverse(-x * 1e8));
    }
}

TEST(SymLogTransform, Constructor) {
    ASSERT_NO_THROW(SymLogTransform<double>(1.0));
    ASSERT_NO_THROW(SymLogTransform<double>(-1.0));
    ASSERT_THROW(SymLogTransform<double>(0.0), std::runtime_error);
}

TEST(SymLogTransform, Function) {
    size_t M = 100;
    size_t N = 1000;
    for(size_t i=0; i<M; ++i) {
        double min_x = RandomDouble() * 2 + 1e-8;
        SymLogTransform<double> transform(min_x);
        SymLogTransform<double> n_transform(min_x);

        EXPECT_TRUE(transform.Function(0.0) == 0.0);

        for(size_t j=0; j<N; ++j) {
            double x0 = RandomDouble() * min_x;
            double x1 = RandomDouble() * min_x;
            double out0 = transform.Function(x0);
            double out1 = transform.Function(x1);
            EXPECT_TRUE(out0 > 0);
            EXPECT_TRUE(out1 > 0);
            EXPECT_NEAR(out0, n_transform.Function(x0), std::abs(out0) * 1e-8);
            EXPECT_NEAR(out1, n_transform.Function(x1), std::abs(out1) * 1e-8);
            EXPECT_TRUE(x0 == out0);
            EXPECT_TRUE(x1 == out1);
            EXPECT_TRUE((x0 - x1) == (out0 - out1));
            x0 = -x0;
            out0 = transform.Function(x0);
            EXPECT_TRUE(out0 < 0);
            EXPECT_TRUE(x0 == out0);
            EXPECT_TRUE((x0 - x1) == (out0 - out1));
            x1 = -x1;
            out1 = transform.Function(x1);
            EXPECT_TRUE(out1 < 0);
            EXPECT_TRUE(x1 == out1);
            EXPECT_TRUE((x0 - x1) == (out0 - out1));

            EXPECT_NEAR(out0, n_transform.Function(x0), std::abs(out0) * 1e-8);
            EXPECT_NEAR(out1, n_transform.Function(x1), std::abs(out1) * 1e-8);

            x0 = (RandomDouble() + 1) * min_x;
            out0 = transform.Function(x0);
            EXPECT_TRUE(out0 > 0);

            EXPECT_TRUE(x0 > x1);
            EXPECT_TRUE(out0 > out1);

            x1 = (RandomDouble() + 1) * min_x;
            out1 = transform.Function(x1);
            EXPECT_TRUE(out1 > 0);
            EXPECT_TRUE(out0 > min_x);
            EXPECT_TRUE(out1 > min_x);
            EXPECT_NEAR(out0, n_transform.Function(x0), std::abs(out0) * 1e-8);
            EXPECT_NEAR(out1, n_transform.Function(x1), std::abs(out1) * 1e-8);
            if(min_x > 1.0) {
                EXPECT_TRUE(out0 < x0);
                EXPECT_TRUE(out1 < x1);
            } else {
                if(log(x0) - x0 < log(min_x) - min_x) {
                    EXPECT_TRUE(out0 < x0);
                } else if(log(x0) - x0 > log(min_x) - min_x) {
                    EXPECT_TRUE(out0 > x0);
                }
                if(log(x1) - x1 < log(min_x) - min_x) {
                    EXPECT_TRUE(out1 < x1);
                } else if(log(x1) - x1 > log(min_x) - min_x) {
                    EXPECT_TRUE(out1 > x1);
                }
            }
            EXPECT_NEAR(log(x0) - log(x1), (out0 - out1), std::abs(log(x0) - log(x1)) * 1e-8);
            EXPECT_NEAR(log(x0) + min_x - log(min_x), out0, std::abs(log(x0) + min_x - log(min_x)) * 1e-8);
            EXPECT_NEAR(log(x1) + min_x - log(min_x), out1, std::abs(log(x1) + min_x - log(min_x)) * 1e-8);

            out0 = transform.Function(-x0);
            EXPECT_TRUE(out0 < 0);
            EXPECT_TRUE(out0 < -min_x);
            EXPECT_TRUE(out1 - out0 > 2.0 * min_x);
            EXPECT_NEAR(2.0 * min_x + log(x1) + log(x0) - 2.0 * log(min_x), out1 - out0, std::abs(2.0 * min_x + log(x1) - log(x0)) * 1e-8);

            out1 = transform.Function(-x1);
            EXPECT_TRUE(out1 < 0);
            EXPECT_TRUE(out1 < -min_x);
            EXPECT_NEAR(log(x1) - log(x0), (out0 - out1), std::abs(log(x0) - log(x1)) * 1e-8);
            EXPECT_NEAR(-log(x0) - min_x + log(min_x), out0, std::abs(-log(x0) - min_x + log(min_x)) * 1e-8);
            EXPECT_NEAR(-log(x1) - min_x + log(min_x), out1, std::abs(-log(x1) - min_x + log(min_x)) * 1e-8);

            EXPECT_NEAR(out0, n_transform.Function(-x0), std::abs(out0) * 1e-8);
            EXPECT_NEAR(out1, n_transform.Function(-x1), std::abs(out1) * 1e-8);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

