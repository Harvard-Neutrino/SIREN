
#include <cmath>
#include <random>
#include <iostream>

#include <gtest/gtest.h>

#include "LeptonInjector/math/Matrix3D.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/math/EulerAngles.h"

using namespace LI::math;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(EulerOrder, RoundTrip) {
    EulerFrame frames[] = {EulerFrame::Static, EulerFrame::Rotating};
    unsigned int int_frames[] = {(unsigned int)(EulerFrame::Static), (unsigned int)(EulerFrame::Rotating)};
    EulerRepetition reps[] = {EulerRepetition::No, EulerRepetition::Yes};
    unsigned int int_reps[] = {(unsigned int)(EulerRepetition::No), (unsigned int)(EulerRepetition::Yes)};
    EulerParity parities[] = {EulerParity::Even, EulerParity::Odd};
    unsigned int int_parities[] = {(unsigned int)(EulerParity::Even), (unsigned int)(EulerParity::Odd)};
    EulerAxis axes[] = {EulerAxis::X, EulerAxis::Y, EulerAxis::Z};
    unsigned int int_axes[] = {(unsigned int)(EulerAxis::X), (unsigned int)(EulerAxis::Y), (unsigned int)(EulerAxis::Z)};
    for(size_t i_frame=0; i_frame<2; ++i_frame) {
        for(size_t i_rep=0; i_rep<2; ++i_rep) {
            for(size_t i_parity=0; i_parity<2; ++i_parity) {
                for(size_t i_axis=0; i_axis<3; ++i_axis) {
                    unsigned int expected_int_order = 0;
                    EulerFrame frame = frames[i_frame];
                    unsigned int int_frame = int_frames[i_frame];
                    EulerRepetition rep = reps[i_rep];
                    unsigned int int_rep = int_reps[i_rep];
                    EulerParity parity = parities[i_parity];
                    unsigned int int_parity = int_parities[i_parity];
                    EulerAxis axis = axes[i_axis];
                    unsigned int int_axis = int_axes[i_axis];
                    expected_int_order |= (int_frame & 1) | ((int_rep & 1) << 1) | ((int_parity & 1) << 2) | ((int_axis & 3) << 3);
                    unsigned int int_order = GetEulerOrder(axis, parity, rep, frame);
                    EXPECT_TRUE(expected_int_order == int_order);
                    EXPECT_TRUE(frame == GetEulerFrame(int_order));
                    EXPECT_TRUE(rep == GetEulerRepetition(int_order));
                    EXPECT_TRUE(parity == GetEulerParity(int_order));
                    EXPECT_TRUE(axis == GetEulerAxis(int_order));
                    ASSERT_NO_THROW(EulerOrder order = (EulerOrder)(int_order););
                }
            }
        }
    }
}

TEST(EulerOrder, AxisI) {
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XYZs) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XYXs) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XZYs) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XZXs) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YZXs) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YZYs) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YXZs) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YXYs) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZXYs) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZXZs) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZYXs) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZYZs) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZYXr) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XYXr) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YZXr) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XZXr) == X);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XZYr) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YZYr) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZXYr) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YXYr) == Y);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::YXZr) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZXZr) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::XYZr) == Z);
    EXPECT_TRUE(GetEulerAxisI(EulerOrder::ZYZr) == Z);
}

TEST(EulerOrder, AxisJ) {
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XYZs) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XYXs) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XZYs) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XZXs) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YZXs) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YZYs) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YXZs) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YXYs) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZXYs) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZXZs) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZYXs) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZYZs) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZYXr) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XYXr) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YZXr) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XZXr) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XZYr) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YZYr) == Z);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZXYr) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YXYr) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::YXZr) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZXZr) == X);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::XYZr) == Y);
    EXPECT_TRUE(GetEulerAxisJ(EulerOrder::ZYZr) == Y);
}

TEST(EulerOrder, AxisK) {
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XYZs) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XYXs) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XZYs) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XZXs) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YZXs) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YZYs) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YXZs) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YXYs) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZXYs) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZXZs) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZYXs) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZYZs) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZYXr) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XYXr) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YZXr) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XZXr) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XZYr) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YZYr) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZXYr) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YXYr) == Z);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::YXZr) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZXZr) == Y);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::XYZr) == X);
    EXPECT_TRUE(GetEulerAxisK(EulerOrder::ZYZr) == X);
}

TEST(EulerOrder, AxisH) {
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XYZs) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XYXs) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XZYs) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XZXs) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YZXs) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YZYs) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YXZs) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YXYs) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZXYs) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZXZs) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZYXs) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZYZs) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZYXr) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XYXr) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YZXr) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XZXr) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XZYr) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YZYr) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZXYr) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YXYr) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::YXZr) == Y);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZXZr) == Z);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::XYZr) == X);
    EXPECT_TRUE(GetEulerAxisH(EulerOrder::ZYZr) == Z);
}

TEST(Constructor, Default) {
    ASSERT_NO_THROW(EulerAngles A;);
    EulerAngles angles;
    EXPECT_TRUE(angles.GetOrder() == EulerOrder::ZXZr);
    EXPECT_TRUE(angles.GetAlpha() == 0);
    EXPECT_TRUE(angles.GetBeta() == 0);
    EXPECT_TRUE(angles.GetGamma() == 0);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

