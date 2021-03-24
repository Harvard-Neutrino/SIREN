
#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include "earthmodel-service/EarthModel.h"

#include "FakeMaterialModel.h"

using namespace earthmodel;

TEST(Constructor, Default)
{
    EXPECT_NO_THROW(EarthModel A);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyPathEmptyModel)
{
    EXPECT_THROW(EarthModel A("", "", materials_file), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyPathEmptyModelEmptyMaterials)
{
    EXPECT_THROW(EarthModel A("", "", ""), char const *);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

