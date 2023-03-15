
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

TEST(LinearInterpolationOperator, Constructor) {
    ASSERT_NO_THROW(LinearInterpolationOperator<double>());
}

TEST(LinearInterpolationOperator, Operator) {
    LinearInterpolationOperator<double> op;
    size_t M = 100;
    size_t N = 1000;
    for(size_t i = 0; i<M; ++i) {
        double x0 = (RandomDouble() - 0.5) * 2;
        double dx = (RandomDouble() - 0.5) * 2;
        double x1 = x0 + dx;
        double y0 = (RandomDouble() - 0.5) * 2;
        double dy = (RandomDouble() - 0.5) * 2;
        double y1 = y0 + dy;
        EXPECT_NEAR(op(x0, y0, x1, y1, x0), y0, std::abs(y0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, x1), y1, std::abs(y1) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, (x0 + x1) / 2.0), (y0 + y1) / 2.0, std::abs((y0 + y1) / 2.0) * 1e-8);
        for(size_t j = 0; j<N; ++j) {
            double t0 = (RandomDouble() - 0.5) * 2 * dx + x0;
            double t1 = t0 + RandomDouble();
            double t2 = t0 + RandomDouble();
            double dt1 = t1 - t0;
            double dt2 = t2 - t0;
            EXPECT_TRUE(dt1 >= 0);
            EXPECT_TRUE(dt2 >= 0);
            double res0 = op(x0, y0, x1, y1, t0);
            double res1 = op(x0, y0, x1, y1, t1);
            double res2 = op(x0, y0, x1, y1, t2);
            double dres1 = res1 - res0;
            double dres2 = res2 - res0;
            double slope1 = dres1 / dt1;
            double slope2 = dres2 / dt2;
            EXPECT_NEAR(slope1, dy/dx, std::abs(dy/dx) * 1e8);
            EXPECT_NEAR(slope1, slope2, std::abs(slope1) * 1e8);

            double alpha = (RandomDouble() - 0.5) * 4.0;
            double t = x0 * alpha + x1 * (1.0 - alpha);
            double expect = y0 * alpha + y1 * (1.0 - alpha);
            double res = op(x0, y0, x1, y1, t);
            EXPECT_NEAR(res, expect, std::abs(expect) * 1e-8);
        }
    }
}

TEST(DropLinearInterpolationOperator, Constructor) {
    ASSERT_NO_THROW(DropLinearInterpolationOperator<double>());
}

TEST(DropLinearInterpolationOperator, Operator) {
    DropLinearInterpolationOperator<double> op;
    size_t M = 100;
    size_t N = 1000;
    for(size_t i = 0; i<M; ++i) {
        double x0 = (RandomDouble() - 0.5) * 2;
        double dx = (RandomDouble() - 0.5) * 2;
        double x1 = x0 + dx;
        double y0 = (RandomDouble() - 0.5) * 2;
        double dy = (RandomDouble() - 0.5) * 2;
        double y1 = y0 + dy;
        EXPECT_NEAR(op(x0, y0, x1, y1, x0), y0, std::abs(y0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, x1), y1, std::abs(y1) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, (x0 + x1) / 2.0), (y0 + y1) / 2.0, std::abs((y0 + y1) / 2.0) * 1e-8);
        for(size_t j = 0; j<N; ++j) {
            double t0 = (RandomDouble() - 0.5) * 2 * dx + x0;
            double t1 = t0 + RandomDouble();
            double t2 = t0 + RandomDouble();
            double dt1 = t1 - t0;
            double dt2 = t2 - t0;
            EXPECT_TRUE(dt1 >= 0);
            EXPECT_TRUE(dt2 >= 0);
            double res0 = op(x0, y0, x1, y1, t0);
            double res1 = op(x0, y0, x1, y1, t1);
            double res2 = op(x0, y0, x1, y1, t2);
            EXPECT_TRUE(op(x0, 0, x1, y1, t0) == 0);
            EXPECT_TRUE(op(x0, y0, x1, 0, t0) == 0);
            EXPECT_TRUE(op(x0, 0, x1, 0, t0) == 0);
            EXPECT_TRUE(op(x0, 0, x1, y1, t1) == 0);
            EXPECT_TRUE(op(x0, y0, x1, 0, t1) == 0);
            EXPECT_TRUE(op(x0, 0, x1, 0, t1) == 0);
            EXPECT_TRUE(op(x0, 0, x1, y1, t2) == 0);
            EXPECT_TRUE(op(x0, y0, x1, 0, t2) == 0);
            EXPECT_TRUE(op(x0, 0, x1, 0, t2) == 0);
            double dres1 = res1 - res0;
            double dres2 = res2 - res0;
            double slope1 = dres1 / dt1;
            double slope2 = dres2 / dt2;
            EXPECT_NEAR(slope1, dy/dx, std::abs(dy/dx) * 1e8);
            EXPECT_NEAR(slope1, slope2, std::abs(slope1) * 1e8);

            double alpha = (RandomDouble() - 0.5) * 4.0;
            double t = x0 * alpha + x1 * (1.0 - alpha);
            double expect = y0 * alpha + y1 * (1.0 - alpha);
            double res = op(x0, y0, x1, y1, t);
            EXPECT_NEAR(res, expect, std::abs(expect) * 1e-8);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

