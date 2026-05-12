#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

// Forward-declare the helper functions from Weighter.tcc (linked via SIREN lib)
// to avoid ODR violations from including the full header.
namespace siren { namespace injection {
    double one_minus_exp_of_negative(double x);
    double log_one_minus_exp_of_negative(double x);
}}

using namespace siren::injection;

// ---------------------------------------------------------------------------
// one_minus_exp_of_negative
// ---------------------------------------------------------------------------

TEST(WeighterHelpers, OneMinusExpSmallX) {
    double x = 1e-5;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-15);
}

TEST(WeighterHelpers, OneMinusExpMediumX) {
    double x = 0.5;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, OneMinusExpLargeX) {
    double x = 5.0;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, OneMinusExpVerySmallX) {
    // Exercises the Taylor expansion branch
    double x = 1e-8;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-15);
}

TEST(WeighterHelpers, OneMinusExpAtBranchPoint) {
    // Near the 0.1 branch boundary
    double x = 0.099;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);

    x = 0.101;
    result = one_minus_exp_of_negative(x);
    exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, OneMinusExpResultBounded) {
    // Result should always be in [0, 1) for non-negative x
    for(double x = 0.0; x <= 20.0; x += 0.1) {
        double result = one_minus_exp_of_negative(x);
        EXPECT_GE(result, 0.0) << "Failed at x=" << x;
        EXPECT_LT(result, 1.0) << "Failed at x=" << x;
    }
}

// ---------------------------------------------------------------------------
// log_one_minus_exp_of_negative
// ---------------------------------------------------------------------------

TEST(WeighterHelpers, LogOneMinusExpSmallX) {
    double x = 1e-5;
    double result = log_one_minus_exp_of_negative(x);
    double exact = std::log(1.0 - std::exp(-x));
    EXPECT_NEAR(result, exact, 1e-10);
}

TEST(WeighterHelpers, LogOneMinusExpMidX) {
    double x = 1.5;
    double result = log_one_minus_exp_of_negative(x);
    double exact = std::log(1.0 - std::exp(-x));
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, LogOneMinusExpLargeX) {
    // Exercises the exp-series branch (x > 3)
    double x = 5.0;
    double result = log_one_minus_exp_of_negative(x);
    double exact = std::log(1.0 - std::exp(-x));
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, LogOneMinusExpAtBranchPoints) {
    // Test near the 0.1 and 3.0 branch boundaries
    for(double x : {0.09, 0.11, 2.99, 3.01}) {
        double result = log_one_minus_exp_of_negative(x);
        double exact = std::log(1.0 - std::exp(-x));
        EXPECT_NEAR(result, exact, 1e-9) << "Failed at x=" << x;
    }
}

TEST(WeighterHelpers, LogOneMinusExpIsNegative) {
    // log(1 - exp(-x)) is always negative for x > 0
    for(double x = 0.001; x <= 20.0; x += 0.1) {
        double result = log_one_minus_exp_of_negative(x);
        EXPECT_LT(result, 0.0) << "Failed at x=" << x;
    }
}
