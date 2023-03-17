
#include <cmath>
#include <tuple>
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

TEST(RegularIndexer1D, Constructor) {
    std::vector<double> x = {0.0, 1.0};
    ASSERT_NO_THROW(RegularIndexer1D<double> A(x));
}

TEST(RegularIndexer1D, EmptyConstructor) {
    std::vector<double> x = {};
    EXPECT_THROW(RegularIndexer1D<double> A(x), std::runtime_error);
}

TEST(RegularIndexer1D, ZeroRangeConstructor) {
    std::vector<double> x = {0.0, 0.0};
    EXPECT_THROW(RegularIndexer1D<double> A(x), std::runtime_error);
}

TEST(RegularIndexer1D, SortedIndex) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double min = RandomDouble();
        double max = min + RandomDouble();
        size_t n_points = size_t(RandomDouble() * 100 + 2);
        double delta = (max - min) / (n_points - 1);
        std::vector<double> points(n_points);
        for(size_t j=0; j<n_points; ++j) {
            points[j] = (min + delta * j);
        }
        RegularIndexer1D<double> indexer(points);
        std::tuple<int, int> index = {0, 1};
        EXPECT_TRUE(index == indexer(min));
        EXPECT_TRUE(index == indexer(min - 1));
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
        EXPECT_TRUE(index == indexer(max + 1));
        for(size_t j=1; j<n_points-1; ++j) {
            index = {j, j+1};
            EXPECT_TRUE(index == indexer(min + delta * j + 1e-4 * delta));
        }
    }
}

TEST(RegularIndexer1D, UnsortedIndex) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double min = RandomDouble();
        double max = min + RandomDouble();
        size_t n_points = size_t(RandomDouble() * 100 + 2);
        double delta = (max - min) / (n_points - 1);
        std::vector<double> points(n_points);
        for(size_t j=0; j<n_points; ++j) {
            points[j] = (min + delta * j);
        }
        std::shuffle(std::begin(points), std::end(points), rng_);
        RegularIndexer1D<double> indexer(points);
        std::tuple<int, int> index = {0, 1};
        EXPECT_TRUE(index == indexer(min));
        EXPECT_TRUE(index == indexer(min - 1));
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
        EXPECT_TRUE(index == indexer(max + 1));
        for(size_t j=1; j<n_points-1; ++j) {
            index = {j, j+1};
            EXPECT_TRUE(index == indexer(min + delta * j + 1e-4 * delta));
        }
    }
}

TEST(IrregularIndexer1D, Constructor) {
    std::vector<double> x = {0.0, 1.0};
    ASSERT_NO_THROW(IrregularIndexer1D<double> A(x));
}

TEST(IrregularIndexer1D, EmptyConstructor) {
    std::vector<double> x = {};
    EXPECT_THROW(IrregularIndexer1D<double> A(x), std::runtime_error);
}

TEST(IrregularIndexer1D, ZeroRangeConstructor) {
    std::vector<double> x = {0.0, 0.0};
    EXPECT_THROW(IrregularIndexer1D<double> A(x), std::runtime_error);
}

TEST(IrregularIndexer1D, SortedIndex) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double min = RandomDouble();
        double range = RandomDouble();
        double max = min + range;
        size_t n_points = size_t(RandomDouble() * 100 + 2);
        std::vector<double> points(n_points);
        for(size_t j=0; j<n_points; ++j) {
            points[j] = RandomDouble() * range + min;
        }
        std::sort(points.begin(), points.end());
        IrregularIndexer1D<double> indexer(points);
        std::tuple<int, int> index = {0, 1};
        EXPECT_TRUE(index == indexer(min));
        EXPECT_TRUE(index == indexer(min - 1));
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
        EXPECT_TRUE(index == indexer(max + 1));
        for(size_t j=1; j<n_points-1; ++j) {
            index = {j, j+1};
            EXPECT_TRUE(index == indexer(points[j] + (points[j+1] - points[j]) * 1e-4));
        }
    }
}

TEST(IrregularIndexer1D, UnsortedIndex) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double min = RandomDouble();
        double range = RandomDouble();
        double max = min + range;
        size_t n_points = size_t(RandomDouble() * 100 + 2);
        std::vector<double> points(n_points);
        for(size_t j=0; j<n_points; ++j) {
            points[j] = RandomDouble() * range + min;
        }
        std::sort(points.begin(), points.end());
        IrregularIndexer1D<double> indexer(points);
        std::tuple<int, int> index = {0, 1};
        EXPECT_TRUE(index == indexer(min));
        EXPECT_TRUE(index == indexer(min - 1));
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
        EXPECT_TRUE(index == indexer(max + 1));
        for(size_t j=1; j<n_points-1; ++j) {
            index = {j, j+1};
            EXPECT_TRUE(index == indexer(points[j] + (points[j+1] - points[j]) * 1e-4));
        }
    }
}

TEST(IrregularIndexer1D, RegularGrid) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double min = RandomDouble();
        double max = min + RandomDouble();
        size_t n_points = size_t(RandomDouble() * 100 + 2);
        double delta = (max - min) / (n_points - 1);
        std::vector<double> points(n_points);
        for(size_t j=0; j<n_points; ++j) {
            points[j] = (min + delta * j);
        }
        IrregularIndexer1D<double> indexer(points);
        IrregularIndexer1D<double> reg_indexer(points);
        std::tuple<int, int> index = {0, 1};
        EXPECT_TRUE(index == indexer(min));
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
        for(size_t j=1; j<n_points-1; ++j) {
            index = {j, j+1};
            EXPECT_TRUE(index == indexer(min + delta * j + 1e-4 * delta));
            EXPECT_TRUE(indexer(min + delta * j + 1e-4 * delta) == reg_indexer(min + delta * j + 1e-4 * delta));
        }
    }
}

TEST(RegularGridIndexer2D, RandomIndex) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x_min = RandomDouble();
        double x_max = x_min + RandomDouble();
        size_t x_n_points = size_t(RandomDouble() * 100 + 2);
        double x_delta = (x_max - x_min) / (x_n_points - 1);
        std::vector<double> x_points(x_n_points);
        for(size_t j=0; j<x_n_points; ++j) {
            x_points[j] = (x_min + x_delta * j);
        }

        double y_min = RandomDouble();
        double y_max = y_min + RandomDouble();
        size_t y_n_points = size_t(RandomDouble() * 100 + 2);
        double y_delta = (y_max - y_min) / (y_n_points - 1);
        std::vector<double> y_points(y_n_points);
        for(size_t j=0; j<y_n_points; ++j) {
            y_points[j] = (y_min + y_delta * j);
        }

        std::vector<double> x(x_n_points*y_n_points);
        std::vector<double> y(x_n_points*y_n_points);
        std::vector<size_t> idxs(x_n_points*y_n_points);
        for(size_t i=0; i<x_n_points*y_n_points; ++i) {
            idxs[i] = i;
        }
        std::shuffle(idxs.begin(), idxs.end(), rng_);
        std::map<std::tuple<size_t, size_t>, size_t> z_map;
        size_t k = 0;
        for(size_t i=0; i<x_n_points; ++i) {
            for(size_t j=0; j<y_n_points; ++j) {
                x[idxs[k]] = x_points[i];
                y[idxs[k]] = y_points[j];
                z_map[{i, j}] = idxs[k];
                k += 1;
            }
        }

        RegularGridIndexer2D<double> indexer(x, y);

        for(size_t i=0; i<x_n_points-1; ++i) {
            double x0 = x_points[i];
            double x1 = x_points[i+1];
            double x_mid = (x0 + x1) / 2.0;
            for(size_t j=0; j<y_n_points-1; ++j) {
                double y0 = y_points[j];
                double y1 = y_points[j+1];
                double y_mid = (y0 + y1) / 2.0;
                size_t k00 = z_map[{i, j}];
                size_t k01 = z_map[{i, j+1}];
                size_t k10 = z_map[{i+1, j}];
                size_t k11 = z_map[{i+1, j+1}];
                Index2D index = indexer(x_mid, y_mid);
                EXPECT_TRUE(i == index.x0);
                EXPECT_TRUE(i+1 == index.x1);
                EXPECT_TRUE(j == index.y0);
                EXPECT_TRUE(j+1 == index.y1);
                EXPECT_TRUE(k00 == index.z00);
                EXPECT_TRUE(k01 == index.z01);
                EXPECT_TRUE(k10 == index.z10);
                EXPECT_TRUE(k11 == index.z11);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

