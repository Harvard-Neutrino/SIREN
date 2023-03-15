
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

TEST(IdentityTransform, RoundTrip) {
    IdentityTransform<double> transform;

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        double y = transform.Function(x);
        double res = transform.Inverse(y);
        EXPECT_NEAR(x, res, std::abs(x) * 1e-8);
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

TEST(LogTransform, RoundTrip) {
    LogTransform<double> transform;

    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        double y = transform.Function(x);
        double res = transform.Inverse(y);
        EXPECT_NEAR(x, res, std::abs(x) * 1e-8);
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
            EXPECT_NEAR(2.0 * min_x + log(x1) + log(x0) - 2.0 * log(min_x), out1 - out0, std::abs(2.0 * min_x + log(x1) + log(x0) - 2.0 * log(min_x)) * 1e-8);

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

TEST(SymLogTransform, Inverse) {
    size_t M = 100;
    size_t N = 1000;
    for(size_t i=0; i<M; ++i) {
        double min_x = RandomDouble() * 2 + 1e-8;
        SymLogTransform<double> transform(min_x);
        SymLogTransform<double> n_transform(min_x);

        EXPECT_TRUE(transform.Inverse(0.0) == 0.0);

        for(size_t j=0; j<N; ++j) {
            double out0 = RandomDouble() * min_x;
            double out1 = RandomDouble() * min_x;
            double x0 = transform.Inverse(out0);
            double x1 = transform.Inverse(out1);
            EXPECT_TRUE(x0 > 0);
            EXPECT_TRUE(x1 > 0);
            EXPECT_NEAR(x0, n_transform.Inverse(out0), std::abs(x0) * 1e-8);
            EXPECT_NEAR(x1, n_transform.Inverse(out1), std::abs(x1) * 1e-8);
            EXPECT_TRUE(x0 == out0);
            EXPECT_TRUE(x1 == out1);
            EXPECT_TRUE((x0 - x1) == (out0 - out1));
            out0 = -out0;
            x0 = transform.Inverse(out0);
            EXPECT_TRUE(x0 < 0);
            EXPECT_TRUE(x0 == out0);
            EXPECT_TRUE((x0 - x1) == (out0 - out1));
            out1 = -x1;
            x1 = transform.Inverse(out1);
            EXPECT_TRUE(x1 < 0);
            EXPECT_TRUE(x1 == out1);
            EXPECT_TRUE((x0 - x1) == (out0 - out1));

            EXPECT_NEAR(x0, n_transform.Inverse(out0), std::abs(out0) * 1e-8);
            EXPECT_NEAR(x1, n_transform.Inverse(out1), std::abs(out1) * 1e-8);

            out0 = min_x + log(RandomDouble() * 2.0 + min_x) - log(min_x);
            EXPECT_TRUE(out0 > min_x);
            x0 = transform.Inverse(out0);
            EXPECT_TRUE(x0 > 0);
            EXPECT_TRUE(x0 > min_x);
            EXPECT_TRUE(x0 > x1);
            EXPECT_TRUE(out0 > out1);

            out1 = min_x + log(RandomDouble() * 2.0 + min_x) - log(min_x);
            x1 = transform.Inverse(out1);
            EXPECT_TRUE(x1 > 0);
            EXPECT_TRUE(x1 > min_x);
            EXPECT_NEAR(x0, n_transform.Inverse(out0), std::abs(out0) * 1e-8);
            EXPECT_NEAR(x1, n_transform.Inverse(out1), std::abs(out1) * 1e-8);
            if(min_x > 1.0) {
                EXPECT_TRUE(out0 < x0);
                EXPECT_TRUE(out1 < x1);
            } else {
                if(out0 < x0) {
                    EXPECT_TRUE(log(x0) - x0 < log(min_x) - min_x);
                } else if(out0 > x0) {
                    EXPECT_TRUE(log(x0) - x0 > log(min_x) - min_x);
                }
                if(out1 < x1) {
                    EXPECT_TRUE(log(x1) - x1 < log(min_x) - min_x);
                } else if(out1 > x1) {
                    EXPECT_TRUE(log(x1) - x1 > log(min_x) - min_x);
                }
            }
            EXPECT_NEAR(log(x0) - log(x1), (out0 - out1), std::abs(log(x0) - log(x1)) * 1e-8);
            EXPECT_NEAR(log(x0) + min_x - log(min_x), out0, std::abs(log(x0) + min_x - log(min_x)) * 1e-8);
            EXPECT_NEAR(log(x1) + min_x - log(min_x), out1, std::abs(log(x1) + min_x - log(min_x)) * 1e-8);

            x0 = transform.Inverse(-out0);
            EXPECT_TRUE(x0 < 0);
            EXPECT_TRUE(x0 < -min_x);
            EXPECT_TRUE(x1 - x0 > 2.0 * min_x);
            EXPECT_NEAR(2.0 * min_x + log(x1) + log(-x0) - 2.0 * log(min_x), out1 + out0, std::abs(2.0 * min_x + log(x1) + log(-x0) - 2.0 * log(min_x)) * 1e-8);
            EXPECT_NEAR((exp(out1) + exp(out0))*exp(log(min_x) - min_x), x1 - x0, std::abs((exp(out1) + exp(out0))*exp(log(min_x) - min_x)) * 1e-8);

            x1 = transform.Inverse(-out1);
            EXPECT_TRUE(x1 < 0);
            EXPECT_TRUE(x1 < -min_x);
            EXPECT_NEAR(log(-x1) - log(-x0), (out1 - out0), std::abs(log(-x1) - log(-x0)) * 1e-8);
            EXPECT_NEAR(-exp(out0 + log(min_x) - min_x), x0, std::abs(-exp(out0 - log(min_x) + min_x)) * 1e-8);
            EXPECT_NEAR(-exp(out1 + log(min_x) - min_x), x1, std::abs(-exp(out1 - log(min_x) + min_x)) * 1e-8);

            EXPECT_NEAR(x0, n_transform.Inverse(-out0), std::abs(x0) * 1e-8);
            EXPECT_NEAR(x1, n_transform.Inverse(-out1), std::abs(x1) * 1e-8);
        }
    }
}

TEST(SymLogTransform, RoundTrip) {
    size_t M = 100;
    size_t N = 1000;
    for(size_t i=0; i<M; ++i) {
        double min_x = RandomDouble() * 4 - 2;
        if(min_x == 0)
            min_x += int(2.0 * RandomDouble() - 1) * 1e-8;
        SymLogTransform<double> transform(min_x);
        for(size_t j=0; j<N; ++j) {
            double x = (RandomDouble() * 4 - 2) * std::abs(min_x);
            double y = transform.Function(x);
            double res = transform.Inverse(y);
            EXPECT_NEAR(x, res, std::abs(x) * 1e-8);
        }
    }
}

TEST(RangeTransform, Constructor) {
    size_t N = 1000;
    for(size_t i=0; i<N; ++i) {
        double x = RandomDouble();
        double y = RandomDouble();
        ASSERT_NO_THROW(RangeTransform<double>(x, x+y));
        ASSERT_NO_THROW(RangeTransform<double>(x+y, x));
        ASSERT_NO_THROW(RangeTransform<double>(-x, y));
        ASSERT_NO_THROW(RangeTransform<double>(x, -y));
        ASSERT_NO_THROW(RangeTransform<double>(-x, -x-y));
        ASSERT_NO_THROW(RangeTransform<double>(-x-y, -x));
        ASSERT_THROW(RangeTransform<double>(x, x), std::runtime_error);
        ASSERT_THROW(RangeTransform<double>(-x, -x), std::runtime_error);
    }
}

TEST(RangeTransform, Function) {
    size_t M = 100;
    size_t N = 1000;
    for(size_t i=0; i<M; ++i) {
        double x = RandomDouble() * 2 - 1;
        double range = RandomDouble() * 2 - 1;
        double y = range + x;
        RangeTransform<double> transform(x, y);
        for(size_t j=0; j<N; ++j) {
            double t0 = RandomDouble() * range + x;
            double res0 = transform.Function(t0);
            EXPECT_TRUE(res0 >= 0);
            EXPECT_TRUE(res0 <= 1);
            t0 = (RandomDouble() - 0.5) * range * 2 + x;
            res0 = transform.Function(t0);
            EXPECT_NEAR(res0, (t0 - x) / range, std::abs(res0) * 1e-8);
            double t1 = t0 + RandomDouble();
            double res1 = transform.Function(t1);
            if(range > 0) {
                EXPECT_TRUE(res1 - res0 > 0);
            } else {
                EXPECT_TRUE(res1 - res0 < 0);
            }
            if(std::abs(range) > 1) {
                EXPECT_TRUE(std::abs(res1 - res0) < std::abs(t1 - t0));
            } else {
                EXPECT_TRUE(std::abs(res1 - res0) >= std::abs(t1 - t0));
            }
        }
    }
}

TEST(RangeTransform, Inverse) {
    size_t M = 100;
    size_t N = 1000;
    for(size_t i=0; i<M; ++i) {
        double x = RandomDouble() * 2 - 1;
        double range = RandomDouble() * 2 - 1;
        double y = range + x;
        RangeTransform<double> transform(x, y);
        for(size_t j=0; j<N; ++j) {
            double t0 = RandomDouble();
            double res0 = transform.Inverse(t0);
            if(range > 0) {
                EXPECT_TRUE(res0 >= x);
                EXPECT_TRUE(res0 <= y);
            } else {
                EXPECT_TRUE(res0 <= x);
                EXPECT_TRUE(res0 >= y);
            }
            t0 = (RandomDouble() - 0.5) * 2;
            res0 = transform.Inverse(t0);
            EXPECT_NEAR(res0, t0 * range + x, std::abs(res0) * 1e-8);
            double t1 = t0 + RandomDouble();
            double res1 = transform.Inverse(t1);
            if(range > 0) {
                EXPECT_TRUE(res1 - res0 > 0);
            } else {
                EXPECT_TRUE(res1 - res0 < 0);
            }
            if(std::abs(range) > 1) {
                EXPECT_TRUE(std::abs(res1 - res0) > std::abs(t1 - t0));
            } else {
                EXPECT_TRUE(std::abs(res1 - res0) <= std::abs(t1 - t0));
            }
        }
    }
}

TEST(RangeTransform, RoundTrip) {
    size_t M = 100;
    size_t N = 1000;
    for(size_t i=0; i<M; ++i) {
        double x = RandomDouble() * 2 - 1;
        double range = RandomDouble() * 2 - 1;
        double y = range + x;
        RangeTransform<double> transform(x, y);
        for(size_t j=0; j<N; ++j) {
            double t0 = (RandomDouble() - 0.5) * range * 2 + x;
            double res0 = transform.Function(t0);
            double res_t = transform.Inverse(res0);
            EXPECT_NEAR(t0, res_t, std::abs(t0) * 1e-8);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

