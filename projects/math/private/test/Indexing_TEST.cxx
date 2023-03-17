
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
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
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
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
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
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
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
        index = {n_points-2, n_points-1};
        EXPECT_TRUE(index == indexer(max));
        for(size_t j=1; j<n_points-1; ++j) {
            index = {j, j+1};
            EXPECT_TRUE(index == indexer(points[j] + (points[j+1] - points[j]) * 1e-4));
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

