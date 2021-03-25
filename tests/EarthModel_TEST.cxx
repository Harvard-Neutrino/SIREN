
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
#include "FakeEarthModel.h"

using namespace earthmodel;

TEST(Constructor, Default)
{
    EXPECT_NO_THROW(EarthModel A);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyModel)
{
    EXPECT_THROW(EarthModel A("", materials_file), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyPathEmptyModel)
{
    EXPECT_THROW(EarthModel A("", "", materials_file), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyModelEmptyMaterials)
{
    EXPECT_THROW(EarthModel A("", ""), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyPathEmptyModelEmptyMaterials)
{
    EXPECT_THROW(EarthModel A("", "", ""), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyPathBadModel)
{
    EXPECT_THROW(EarthModel("", std::tmpnam(nullptr), materials_file), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorBadModel)
{
    EXPECT_THROW(EarthModel(std::tmpnam(nullptr), materials_file), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorBadModelEmptyMaterial)
{
    EXPECT_THROW(EarthModel(std::tmpnam(nullptr), ""), char const *);
}

TEST_F(FakeLegacyEarthModelTest, LoadConcentricShellsFromLegacyFile)
{
    EarthModel A;
    ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
    EXPECT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, 1000, -1));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

