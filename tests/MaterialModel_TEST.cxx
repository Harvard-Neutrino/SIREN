
#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include "earthmodel-service/MaterialModel.h"

#include "FakeMaterialModel.h"

using namespace earthmodel;

TEST(Constructor, Default)
{
    MaterialModel A;
}

TEST_F(FakeMaterialModelTest, FileConstructor)
{
    ASSERT_NO_THROW(MaterialModel(materials_file));
    MaterialModel A(materials_file);
    ASSERT_NO_THROW(MaterialModel("", materials_file));
    MaterialModel B("", materials_file);
}

TEST(Constructor, EmptyFile)
{
    EXPECT_THROW(MaterialModel A(""), char const *);
}

TEST(Constructor, BadFile)
{
    EXPECT_THROW(MaterialModel A(std::tmpnam(nullptr)), char const *);
}

TEST_F(FakeMaterialModelTest, DuplicateFile)
{
    ASSERT_NO_THROW(MaterialModel(materials_file));
    MaterialModel A(materials_file);
    ASSERT_NO_THROW(A.AddModelFile(materials_file));
    ASSERT_NO_THROW(MaterialModel(std::vector<std::string>{materials_file, materials_file}));
    ASSERT_NO_THROW(A.HasMaterial(0));

    ASSERT_NO_THROW(MaterialModel("", materials_file));
    MaterialModel B("", materials_file);
    ASSERT_NO_THROW(B.AddModelFile(materials_file));
    ASSERT_NO_THROW(MaterialModel("", {materials_file, materials_file}));
    ASSERT_NO_THROW(B.HasMaterial(0));
}

std::string print_v(std::vector<std::string> v) {
    std::stringstream s;
    for(auto ss : v) {
        s << ss << " ";
    }
    return s.str();
}

template <typename T, typename U>
std::string print_m(std::map<T, U> v) {
    std::stringstream s;
    for(auto ss : v) {
        s << ss.first << ", " << ss.second << std::endl;;
    }
    return s.str();
}

template <typename T, typename U>
std::string print_vp(std::vector<std::pair<T, U>> v) {
    std::stringstream s;
    for(auto ss : v) {
        s << ss.first << ", " << ss.second << std::endl;;
    }
    return s.str();
}

TEST_F(FakeMaterialModelTest, MaterialCount)
{
    ASSERT_NO_THROW(MaterialModel("", materials_file));

    unsigned int N_RAND = 1000;
    for(unsigned int i=0; i<N_RAND; ++i) {
        clear();
        create_file({}, 12);
        MaterialModel A(materials_file);

        int material_count = 0;
        std::vector<std::string> names;
        while(A.HasMaterial(material_count)) {
            names.push_back(A.GetMaterialName(material_count));
            ASSERT_LT(material_count, material_names_.size());
            EXPECT_EQ(std::string(names.back()), std::string(material_names_[material_count])) << file_contents;
            EXPECT_EQ(A.GetMaterialMap(material_count).size(), material_components_[material_names_[material_count]].size()) << A.GetMaterialName(material_count) << " " << material_names_[material_count] << std::endl << print_m(A.GetMaterialMap(material_count)) << std::endl << std::endl << print_vp(material_components_[material_names_[material_count]]);
            material_count += 1;
        }

        ASSERT_EQ(material_count, material_names_.size()) << file_contents;
    }
}

TEST_F(FakeMaterialModelTest, DuplicateFileMaterialCount)
{
    ASSERT_NO_THROW(MaterialModel(materials_file));

    unsigned int N_RAND = 1000;
    for(unsigned int i=0; i<N_RAND; ++i) {
        clear();
        create_file({}, 12);
        MaterialModel A(materials_file);

        int material_count = 0;
        std::vector<std::string> names;
        while(A.HasMaterial(material_count)) {
            names.push_back(A.GetMaterialName(material_count));
            ASSERT_LT(material_count, material_names_.size());
            EXPECT_EQ(std::string(names.back()), std::string(material_names_[material_count])) << file_contents;
            EXPECT_EQ(A.GetMaterialMap(material_count).size(), material_components_[material_names_[material_count]].size()) << file_contents;
            material_count += 1;
        }

        int first_material_count = material_count;
        A.AddModelFile(materials_file);

        material_count = 0;
        names.clear();
        while(A.HasMaterial(material_count)) {
            names.push_back(A.GetMaterialName(material_count));
            ASSERT_LT(material_count, material_names_.size());
            EXPECT_EQ(std::string(names.back()), std::string(material_names_[material_count]));
            EXPECT_EQ(A.GetMaterialMap(material_count).size(), material_components_[material_names_[material_count]].size());
            material_count += 1;
        }

        EXPECT_EQ(material_count, material_names_.size()) << file_contents;
        EXPECT_EQ(material_count, first_material_count) << file_contents;
    }
}

TEST_F(FakeMaterialModelTest, MaterialId)
{
    ASSERT_NO_THROW(MaterialModel(materials_file));

    unsigned int N_RAND = 100;
    for(unsigned int i=0; i<N_RAND; ++i) {
        clear();
        create_file({}, 12);
        MaterialModel A(materials_file);

        for(unsigned int i=0; i<material_names_.size(); ++i) {
            std::string name = material_names_[i];
            ASSERT_TRUE(A.HasMaterial(name));
            ASSERT_TRUE(A.HasMaterial(i));
            EXPECT_EQ(A.GetMaterialName(i), name);
            EXPECT_EQ(A.GetMaterialId(name), i);
        }
    }
}

TEST_F(FakeMaterialModelTest, MaterialPNE)
{
    ASSERT_NO_THROW(MaterialModel(materials_file));

    unsigned int N_RAND = 100;
    for(unsigned int i=0; i<N_RAND; ++i) {
        clear();
        create_file({}, 12);
        MaterialModel A(materials_file);

        for(unsigned int i=0; i<material_names_.size(); ++i) {
            std::string name = material_names_[i];
            ASSERT_TRUE(A.HasMaterial(name));
            ASSERT_TRUE(A.HasMaterial(i));
            EXPECT_DOUBLE_EQ(A.GetPNERatio(i), pne_ratios_[i]);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

