
#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include "earthmodel-service/EarthModel.h"

#include "MaterialFake.h"

using namespace earthmodel;

TEST(Constructor, Default)
{
    EarthModel A;
}

TEST_F(MaterialTest, EarthModelConstructor)
{
    EarthModel A("", "", materials_file);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

