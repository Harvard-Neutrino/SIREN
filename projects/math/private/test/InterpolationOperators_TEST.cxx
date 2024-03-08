
#include <cmath>
#include <math.h>
#include <random>
#include <iostream>

#include <gtest/gtest.h>

#include "SIREN/math/Interpolation.h"

using namespace SI::math;

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

TEST(BiLinearInterpolationOperator, Constructor) {
    ASSERT_NO_THROW(BiLinearInterpolationOperator<double>());
}

TEST(BiLinearInterpolationOperator, Operator) {
    BiLinearInterpolationOperator<double> op;
    size_t M = 100;
    size_t N = 1000;
    for(size_t i = 0; i<M; ++i) {
        double x0 = (RandomDouble() - 0.5) * 2;
        double dx = (RandomDouble() - 0.5) * 2;
        double x1 = x0 + dx;
        double y0 = (RandomDouble() - 0.5) * 2;
        double dy = (RandomDouble() - 0.5) * 2;
        double y1 = y0 + dy;
        double z00 = (RandomDouble() - 0.5) * 2;
        double z01 = (RandomDouble() - 0.5) * 2;
        double z10 = (RandomDouble() - 0.5) * 2;
        double z11 = (RandomDouble() - 0.5) * 2;
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x0, y0), z00, std::abs(z00) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x0, y1), z01, std::abs(z01) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x1, y0), z10, std::abs(z10) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x1, y1), z11, std::abs(z11) * 1e-8);
        double mid_x = (x0 + x1) / 2.0;
        double mid_y = (y0 + y1) / 2.0;
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, mid_x, mid_y), (z00 + z01 + z10 + z11)/4.0, std::abs((z00 + z01 + z10 + z11)/4.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x0, mid_y), (z00 + z01)/2.0, std::abs((z00 + z01)/2.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x1, mid_y), (z10 + z11)/2.0, std::abs((z10 + z11)/2.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, mid_x, y0), (z00 + z10)/2.0, std::abs((z00 + z10)/2.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, mid_x, y1), (z01 + z11)/2.0, std::abs((z01 + z11)/2.0) * 1e-8);
        double tx;
        double ty;
        double b0, b1, db;
        double z0, z1, dz;
        std::function<double()> evaluate = [&]()->double{
            return op(x0, y0, x1, y1, z00, z01, z10, z11, tx, ty);
        };
        std::function<void(double)> update_x = [&](double t)->void{
            tx = t;
        };
        std::function<void(double)> update_y = [&](double t)->void{
            ty = t;
        };
        std::function<void(double)> update;
        for(size_t l=0; l<4; ++l) {
            if(l == 0) {
                tx = x0;
                b0 = y0; b1 = y1;
                z0 = z00; z1 = z01;
                update = update_y;
            } else if(l == 1) {
                tx = x1;
                b0 = y0; b1 = y1;
                z0 = z10; z1 = z11;
                update = update_y;
            } else if(l == 2) {
                ty = y0;
                b0 = x0; b1 = x1;
                z0 = z00; z1 = z10;
                update = update_x;
            } else if(l == 3) {
                ty = y1;
                b0 = x0; b1 = x1;
                z0 = z01; z1 = z11;
                update = update_x;
            }
            db = b1 - b0;
            dz = z1 - z0;
            for(size_t j=0; j<N; ++j) {
                double t0 = (RandomDouble() - 0.25) * 2 * db + b0;
                double t1 = t0 + RandomDouble();
                double t2 = t0 + RandomDouble();
                update(t0);
                double res0 = evaluate();
                update(t1);
                double res1 = evaluate();
                update(t2);
                double res2 = evaluate();
                double dres1 = res1 - res0;
                double dres2 = res2 - res0;
                double slope1 = dres1 / (t1 - t0);
                double slope2 = dres2 / (t2 - t0);

                EXPECT_NEAR(slope1, dz/db, std::abs(dz/db) * 1e8);
                EXPECT_NEAR(slope1, slope2, std::abs(slope1) * 1e8);

                double alpha = (RandomDouble() - 0.5) * 4.0;
                double t = b0 * alpha + b1 * (1.0 - alpha);
                double expect = z0 * alpha + z1 * (1.0 - alpha);
                update(t);
                double res = evaluate();
                EXPECT_NEAR(res, expect, std::abs(expect) * 1e-8);
            }
        }
    }
}

TEST(DropBiLinearInterpolationOperator, Constructor) {
    ASSERT_NO_THROW(DropBiLinearInterpolationOperator<double>());
}

TEST(DropBiLinearInterpolationOperator, Operator) {
    DropBiLinearInterpolationOperator<double> op;
    size_t M = 100;
    size_t N = 1000;
    for(size_t i = 0; i<M; ++i) {
        double x0 = (RandomDouble() - 0.5) * 2;
        double dx = (RandomDouble() - 0.5) * 2;
        double x1 = x0 + dx;
        double y0 = (RandomDouble() - 0.5) * 2;
        double dy = (RandomDouble() - 0.5) * 2;
        double y1 = y0 + dy;
        double z00 = (RandomDouble() - 0.5) * 2;
        double z01 = (RandomDouble() - 0.5) * 2;
        double z10 = (RandomDouble() - 0.5) * 2;
        double z11 = (RandomDouble() - 0.5) * 2;
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x0, y0), z00, std::abs(z00) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x0, y1), z01, std::abs(z01) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x1, y0), z10, std::abs(z10) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x1, y1), z11, std::abs(z11) * 1e-8);
        double mid_x = (x0 + x1) / 2.0;
        double mid_y = (y0 + y1) / 2.0;
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, mid_x, mid_y), (z00 + z01 + z10 + z11)/4.0, std::abs((z00 + z01 + z10 + z11)/4.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x0, mid_y), (z00 + z01)/2.0, std::abs((z00 + z01)/2.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, x1, mid_y), (z10 + z11)/2.0, std::abs((z10 + z11)/2.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, mid_x, y0), (z00 + z10)/2.0, std::abs((z00 + z10)/2.0) * 1e-8);
        EXPECT_NEAR(op(x0, y0, x1, y1, z00, z01, z10, z11, mid_x, y1), (z01 + z11)/2.0, std::abs((z01 + z11)/2.0) * 1e-8);
        double tx;
        double ty;
        double b0, b1, db;
        double z0, z1, dz;
        std::function<double()> drop_eval_00 = [&]()->double{
            return op(x0, y0, x1, y1, 0, 0, z10, z11, tx, ty);
        };
        std::function<double()> drop_eval_01 = [&]()->double{
            return op(x0, y0, x1, y1, z00, 0, z10, 0, tx, ty);
        };
        std::function<double()> drop_eval_10 = [&]()->double{
            return op(x0, y0, x1, y1, z00, z01, 0, 0, tx, ty);
        };
        std::function<double()> drop_eval_11 = [&]()->double{
            return op(x0, y0, x1, y1, 0, z01, 0, z11, tx, ty);
        };
        std::function<double()> evaluate = [&]()->double{
            return op(x0, y0, x1, y1, z00, z01, z10, z11, tx, ty);
        };
        std::function<void(double)> update_x = [&](double t)->void{
            tx = t;
        };
        std::function<void(double)> update_y = [&](double t)->void{
            ty = t;
        };
        std::function<void(double)> update;
        for(size_t l=0; l<4; ++l) {
            if(l == 0) {
                tx = x0;
                b0 = y0; b1 = y1;
                z0 = z00; z1 = z01;
                update = update_y;
            } else if(l == 1) {
                tx = x1;
                b0 = y0; b1 = y1;
                z0 = z10; z1 = z11;
                update = update_y;
            } else if(l == 2) {
                ty = y0;
                b0 = x0; b1 = x1;
                z0 = z00; z1 = z10;
                update = update_x;
            } else if(l == 3) {
                ty = y1;
                b0 = x0; b1 = x1;
                z0 = z01; z1 = z11;
                update = update_x;
            }
            db = b1 - b0;
            dz = z1 - z0;
            for(size_t j=0; j<N; ++j) {
                double t0 = (RandomDouble() - 0.25) * 2 * db + b0;
                double t1 = t0 + RandomDouble();
                double t2 = t0 + RandomDouble();
                update(t0);
                EXPECT_TRUE(drop_eval_00() == 0);
                EXPECT_TRUE(drop_eval_01() == 0);
                EXPECT_TRUE(drop_eval_10() == 0);
                EXPECT_TRUE(drop_eval_11() == 0);
                double res0 = evaluate();
                update(t1);
                EXPECT_TRUE(drop_eval_00() == 0);
                EXPECT_TRUE(drop_eval_01() == 0);
                EXPECT_TRUE(drop_eval_10() == 0);
                EXPECT_TRUE(drop_eval_11() == 0);
                double res1 = evaluate();
                update(t2);
                EXPECT_TRUE(drop_eval_00() == 0);
                EXPECT_TRUE(drop_eval_01() == 0);
                EXPECT_TRUE(drop_eval_10() == 0);
                EXPECT_TRUE(drop_eval_11() == 0);
                double res2 = evaluate();
                double dres1 = res1 - res0;
                double dres2 = res2 - res0;
                double slope1 = dres1 / (t1 - t0);
                double slope2 = dres2 / (t2 - t0);

                EXPECT_NEAR(slope1, dz/db, std::abs(dz/db) * 1e8);
                EXPECT_NEAR(slope1, slope2, std::abs(slope1) * 1e8);

                double alpha = (RandomDouble() - 0.5) * 4.0;
                double t = b0 * alpha + b1 * (1.0 - alpha);
                double expect = z0 * alpha + z1 * (1.0 - alpha);
                update(t);
                double res = evaluate();
                EXPECT_NEAR(res, expect, std::abs(expect) * 1e-8);
            }
        }
    }
}

TEST(SimplexLinearInterpolationOperator, Constructor) {
    ASSERT_NO_THROW(SimplexLinearInterpolationOperator<double>());
}

TEST(SimplexLinearInterpolationOperator, Operator) {
    using Tri = typename DelaunayIndexer2D<double>::Tri;
    SimplexLinearInterpolationOperator<double> op;
    std::shared_ptr<Tri> tri(new Tri());
    size_t M = 100;
    for(size_t i = 0; i<M; ++i) {
        double x0 = (RandomDouble() - 0.5) * 2; tri->p0.x = x0;
        double x1 = (RandomDouble() - 0.5) * 2; tri->p1.x = x1;
        double x2 = (RandomDouble() - 0.5) * 2; tri->p2.x = x2;
        double y0 = (RandomDouble() - 0.5) * 2; tri->p0.y = y0;
        double y1 = (RandomDouble() - 0.5) * 2; tri->p1.y = y1;
        double y2 = (RandomDouble() - 0.5) * 2; tri->p2.y = y2;
        double z0 = (RandomDouble() - 0.5) * 2;
        double z1 = (RandomDouble() - 0.5) * 2;
        double z2 = (RandomDouble() - 0.5) * 2;
        EXPECT_NEAR(op(x0, y0, tri.get(), z0, z1, z2), z0, std::abs(z0) * 1e-8);
        EXPECT_NEAR(op(x1, y1, tri.get(), z0, z1, z2), z1, std::abs(z1) * 1e-8);
        EXPECT_NEAR(op(x2, y2, tri.get(), z0, z1, z2), z2, std::abs(z2) * 1e-8);
        double mid_x0 = (x1 + x2) / 2.0;
        double mid_x1 = (x0 + x2) / 2.0;
        double mid_x2 = (x0 + x1) / 2.0;
        double mid_y0 = (y1 + y2) / 2.0;
        double mid_y1 = (y0 + y2) / 2.0;
        double mid_y2 = (y0 + y1) / 2.0;
        EXPECT_NEAR(op(mid_x0, mid_y0, tri.get(), z0, z1, z2), (z1 + z2)/2.0, std::abs((z1 + z2)/2.0) * 1e-8);
        EXPECT_NEAR(op(mid_x1, mid_y1, tri.get(), z0, z1, z2), (z0 + z2)/2.0, std::abs((z0 + z2)/2.0) * 1e-8);
        EXPECT_NEAR(op(mid_x2, mid_y2, tri.get(), z0, z1, z2), (z1 + z0)/2.0, std::abs((z1 + z0)/2.0) * 1e-8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

