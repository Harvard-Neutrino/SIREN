
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

TEST(LinearDelaunayInterpolator2D, Constructor) {
    using Simplex = typename IDelaBella2<double>::Simplex;
    using Vertex = typename IDelaBella2<double>::Vertex;
    size_t N = 100;
    size_t M = 1000;
    for(size_t m=0; m<N; ++m) {
        double x_min = RandomDouble();
        double x_range = RandomDouble();
        double x_max = x_min + x_range;
        double y_min = RandomDouble();
        double y_range = RandomDouble();
        double y_max = y_min + y_range;
        size_t n_points = size_t(RandomDouble() * 1000 + 4);
        std::vector<double> x_points(n_points);
        std::vector<double> y_points(n_points);
        std::vector<double> z_points(n_points);
        for(size_t j=0; j<n_points - 4; ++j) {
            x_points[j] = RandomDouble() * y_range + y_min;
            y_points[j] = RandomDouble() * y_range + y_min;
        }
        x_points[n_points - 4] = x_min;
        x_points[n_points - 3] = x_min;
        x_points[n_points - 2] = x_max;
        x_points[n_points - 1] = x_max;
        y_points[n_points - 4] = y_min;
        y_points[n_points - 3] = y_max;
        y_points[n_points - 2] = y_min;
        y_points[n_points - 1] = y_max;

        std::function<double(double, double)> f = [&](double x, double y)->double {
            return x * y + x + y + 2.0;
        };
        for(size_t j=0; j<n_points; ++j) {
            z_points[j] = f(x_points[j], y_points[j]);
        }

        std::shared_ptr<Transform<double>> x_transform(new IdentityTransform<double>());
        std::shared_ptr<Transform<double>> y_transform(new IdentityTransform<double>());
        std::shared_ptr<Transform<double>> z_transform(new IdentityTransform<double>());

        std::shared_ptr<SimplexLinearInterpolationOperator<double>> op(new SimplexLinearInterpolationOperator<double>);

        LinearDelaunayInterpolator2D<double> interp(x_points, y_points, z_points, x_transform, y_transform, z_transform, op);
    }
}

TEST(LinearDelaunayInterpolator2D, Interpolate) {
    using Simplex = typename IDelaBella2<double>::Simplex;
    using Vertex = typename IDelaBella2<double>::Vertex;
    size_t N = 100;
    size_t M = 1000;
    for(size_t m=0; m<N; ++m) {
        double x_min = RandomDouble();
        double x_range = RandomDouble();
        double x_max = x_min + x_range;
        double y_min = RandomDouble();
        double y_range = RandomDouble();
        double y_max = y_min + y_range;
        size_t n_points = size_t(RandomDouble() * 1000);
        std::vector<double> x_points(n_points + 4);
        std::vector<double> y_points(n_points + 4);
        std::vector<double> z_points(n_points + 4);
        for(size_t j=0; j<n_points; ++j) {
            x_points[j] = RandomDouble() * x_range + x_min;
            y_points[j] = RandomDouble() * y_range + y_min;
        }

        x_points[n_points + 0] = x_min;
        y_points[n_points + 0] = y_min;

        x_points[n_points + 1] = x_max;
        y_points[n_points + 1] = y_min;

        x_points[n_points + 2] = x_max;
        y_points[n_points + 2] = y_max;

        x_points[n_points + 3] = x_min;
        y_points[n_points + 3] = y_max;


        std::function<double(double, double)> f = [&](double x, double y)->double {
            return x * y + x + y + 2.0;
        };
        for(size_t j=0; j<n_points+4; ++j) {
            z_points[j] = f(x_points[j], y_points[j]);
        }

        std::shared_ptr<Transform<double>> x_transform(new IdentityTransform<double>());
        std::shared_ptr<Transform<double>> y_transform(new IdentityTransform<double>());
        std::shared_ptr<Transform<double>> z_transform(new IdentityTransform<double>());

        std::shared_ptr<SimplexLinearInterpolationOperator<double>> op(new SimplexLinearInterpolationOperator<double>);

        LinearDelaunayInterpolator2D<double> interp(x_points, y_points, z_points, x_transform, y_transform, z_transform, op);
        double x = x_range/2.0 + x_min;
        double y = y_range/2.0 + y_min;
        double z = f(x, y);
        EXPECT_NEAR(z, interp(x, y), std::abs(z) * 1e-2);
        for(size_t i=0; i<M; ++i) {
            double x = RandomDouble() * x_range + x_min;
            double y = RandomDouble() * y_range + y_min;
            double z = f(x, y);
            EXPECT_NEAR(z, interp(x, y), std::abs(z) * 1e-2);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

