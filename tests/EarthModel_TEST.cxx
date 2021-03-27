
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

TEST(DefaultMaterials, VacuumOnly)
{
    EarthModel A;
    MaterialModel materials = A.GetMaterials();

    int material_count = 0;
    while(materials.HasMaterial(material_count)) {
        material_count += 1;
    }
    EXPECT_EQ(material_count, 1);

    std::string name = "VACUUM";

    ASSERT_TRUE(materials.HasMaterial(name));

    int id = materials.GetMaterialId(name);

    ASSERT_EQ(0, id);
    EXPECT_DOUBLE_EQ(7.0/8.0, materials.GetPNERatio(id));

    int component_id = 1000070080;
    std::map<int, double> material_map = materials.GetMaterialMap(id);
    ASSERT_EQ(1, material_map.size());
    ASSERT_EQ(1, material_map.count(component_id));
    EXPECT_DOUBLE_EQ(1.0, material_map[component_id]);
}

/*
    EarthSector sector;
    sector.material_id = materials_.GetMaterialId("VACUUM");
    sector.level = sectors_.size();
    sector.geo = Sphere(Vector3D(0,0,0), std::numeric_limits<double>::infinity(), 0).create();
    sector.density = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>().create(); // Use the universe_mean_density from GEANT4
    sectors_.push_back(sector);
*/

TEST(DefaultSectors, VacuumOnly)
{
    EarthModel A;
    std::vector<EarthSector> sectors = A.GetSectors();
    ASSERT_EQ(1, sectors.size());
    EarthSector sector = sectors[0];
    EXPECT_EQ(0, sector.material_id);
    EXPECT_EQ(0, sector.level);
    Sphere geo(Vector3D(0,0,0), std::numeric_limits<double>::infinity(), 0);
    DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> density;
    EXPECT_EQ(geo, *sector.geo);
    EXPECT_EQ(density, *sector.density);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyModel)
{
    EXPECT_THROW(EarthModel A("", materials_file), char const *);
}

TEST_F(FakeMaterialModelTest, EarthModelConstructorEmptyPathEmptyModel)
{
    EXPECT_THROW(EarthModel A("", "", materials_file), char const *);
}

TEST(Constructor, EarthModelConstructorEmptyModelEmptyMaterials)
{
    EXPECT_THROW(EarthModel A("", ""), char const *);
}

TEST(Constructor, EarthModelConstructorEmptyPathEmptyModelEmptyMaterials)
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

TEST_F(FakeLegacyEarthModelTest, LegacyFileLoadWithNoIce)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        EarthModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        EXPECT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
    }
}

TEST_F(FakeLegacyEarthModelTest, LegacyFileLoadWithIce)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        EarthModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyEarthModelFile::RandomDouble()*180;
        EXPECT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
    }
}

TEST_F(FakeLegacyEarthModelTest, LegacyFileLayerNames)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        EarthModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyEarthModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A.GetSectors();
        EXPECT_EQ(layer_names.size(), sectors.size() - 1);
        unsigned int n_layers = std::min(layer_names.size(), sectors.size() - 1);
        for(unsigned int j=0; j<n_layers; ++j) {
            EXPECT_EQ(layer_names[j], sectors[j+1].name);
        }
    }
}

/*
    std::vector<std::string> layer_names;
    std::vector<double> layer_thicknesses;
    std::vector<double> layer_radii;
    std::vector<Polynom> layer_densities;
    std::vector<std::string> layer_materials;
    std::vector<std::string> layer_types;
*/

TEST_F(FakeLegacyEarthModelTest, LegacyFileLayerThickness)
{

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

