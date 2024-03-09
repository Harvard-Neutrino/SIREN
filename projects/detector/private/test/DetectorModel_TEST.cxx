
#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <gtest/gtest.h>

#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/detector/DensityDistribution.h"
#include "SIREN/detector/Distribution1D.h"
#include "SIREN/detector/Axis1D.h"
#include "SIREN/detector/CartesianAxis1D.h"
#include "SIREN/detector/ConstantDistribution1D.h"
#include "SIREN/detector/ConstantDensityDistribution.h"
#include "SIREN/detector/PolynomialDistribution1D.h"
#include "SIREN/detector/ExponentialDistribution1D.h"
#include "SIREN/detector/DensityDistribution1D.h"

#include "FakeDetectorModel.h"
#include "FakeMaterialModel.h"

using namespace siren::math;
using namespace siren::geometry;
using namespace siren::detector;
using namespace siren::utilities;
using namespace siren::dataclasses;

TEST(Constructor, Default)
{
    EXPECT_NO_THROW(DetectorModel A);
}

TEST(DefaultMaterials, VacuumOnly)
{
    DetectorModel A;
    MaterialModel materials = A.GetMaterials();

    int material_count = 0;
    while(materials.HasMaterial(material_count)) {
        material_count += 1;
    }
    EXPECT_EQ(material_count, 1);

    std::string name = "VACUUM";

    ASSERT_TRUE(materials.HasMaterial(name));

    int id = materials.GetMaterialId(name);

    double material_nucleons_per_gram = materials.GetTargetParticleFraction(id, ParticleType::Nucleon);
    double material_neutrons_per_gram = materials.GetTargetParticleFraction(id, ParticleType::Neutron);
    double material_protons_per_gram = materials.GetTargetParticleFraction(id, ParticleType::PPlus);
    double material_electrons_per_gram = materials.GetTargetParticleFraction(id, ParticleType::EMinus);

    const double nucleons_per_amu = 0.9943511899082073921408545013734389386622278134067745520545013436;
    const double nucleons_per_gram = nucleons_per_amu * Constants::avogadro;

    const double protons_per_amu = 0.8464425697382936033875383253155218234388231142436674312768944729;
    const double protons_per_gram = protons_per_amu * Constants::avogadro;

    const double neutrons_per_amu = 0.1479086201699137887533161760579171152234046991631071207776068707;
    const double neutrons_per_gram = neutrons_per_amu * Constants::avogadro;

    const double electrons_per_gram = protons_per_gram;

    ASSERT_EQ(0, id);
    EXPECT_DOUBLE_EQ(material_nucleons_per_gram, nucleons_per_gram);
    EXPECT_DOUBLE_EQ(material_neutrons_per_gram, neutrons_per_gram);
    EXPECT_DOUBLE_EQ(material_protons_per_gram, protons_per_gram);
    EXPECT_DOUBLE_EQ(material_electrons_per_gram, electrons_per_gram);

    // int component_id = 1000070080;
    // std::map<int, double> material_map = materials.GetMaterialMassFracs(id);
    // ASSERT_EQ(1, material_map.size());
    // ASSERT_EQ(1, material_map.count(component_id));
    // EXPECT_DOUBLE_EQ(1.0, material_map[component_id]);
}

/*
    DetectorSector sector;
    sector.material_id = materials_.GetMaterialId("VACUUM");
    sector.level = sectors_.size();
    sector.geo = Sphere(Vector3D(0,0,0), std::numeric_limits<double>::infinity(), 0).create();
    sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>().create(); // Use the universe_mean_density from GEANT4
    sectors_.push_back(sector);
*/

TEST(DefaultSectors, VacuumOnly)
{
    DetectorModel A;
    std::vector<DetectorSector> sectors = A.GetSectors();
    ASSERT_EQ(1, sectors.size());
    DetectorSector sector = sectors[0];
    EXPECT_EQ(0, sector.material_id);
    EXPECT_EQ(std::numeric_limits<int>::min(), sector.level);
    Sphere geo(Vector3D(0,0,0),std::numeric_limits<double>::infinity(), 0);
    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> density;
    EXPECT_EQ(geo, *sector.geo);
    EXPECT_EQ(density, *sector.density);
}

TEST_F(FakeMaterialModelTest, DetectorModelConstructorEmptyModel)
{
    EXPECT_THROW(DetectorModel A("", materials_file), std::runtime_error);
}

TEST_F(FakeMaterialModelTest, DetectorModelConstructorEmptyPathEmptyModel)
{
    EXPECT_THROW(DetectorModel A("", "", materials_file), std::runtime_error);
}

TEST(Constructor, DetectorModelConstructorEmptyModelEmptyMaterials)
{
    EXPECT_THROW(DetectorModel A("", ""), std::runtime_error);
}

TEST(Constructor, DetectorModelConstructorEmptyPathEmptyModelEmptyMaterials)
{
    EXPECT_THROW(DetectorModel A("", "", ""), std::runtime_error);
}

TEST_F(FakeMaterialModelTest, DetectorModelConstructorEmptyPathBadModel)
{
    EXPECT_THROW(DetectorModel("", std::tmpnam(nullptr), materials_file), std::runtime_error);
}

TEST_F(FakeMaterialModelTest, DetectorModelConstructorBadModel)
{
    EXPECT_THROW(DetectorModel(std::tmpnam(nullptr), materials_file), std::runtime_error);
}

TEST_F(FakeMaterialModelTest, DetectorModelConstructorBadModelEmptyMaterial)
{
    EXPECT_THROW(DetectorModel(std::tmpnam(nullptr), ""), std::runtime_error);
}

TEST(Path, SetGet)
{
    DetectorModel A;
    std::string result;
    std::string expect;
    A.SetPath("a");
    result = A.GetPath();
    expect = "a";
    EXPECT_EQ(expect, result);
    result = A.GetPath();
    EXPECT_EQ(expect, result);
    A.SetPath("b");
    result = A.GetPath();
    expect = "b";
    EXPECT_EQ(expect, result);
    A.SetPath("d");
    A.SetPath("e");
    result = A.GetPath();
    expect = "e";
    EXPECT_EQ(expect, result);
}

TEST(EdgeCases, ColumnDepthWithEqualPoints)
{
    DetectorModel A;
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        Vector3D v;
        EXPECT_EQ(0.0, A.GetColumnDepthInCGS(DetectorPosition(v), DetectorPosition(v)));
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileLoadWithNoIce)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        EXPECT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileLoadWithIce)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        EXPECT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileLayerNames)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        EXPECT_EQ(layer_names.size(), sectors.size() - 1);
        unsigned int n_layers = std::min(layer_names.size(), sectors.size() - 1);
        for(unsigned int j=0; j<n_layers; ++j) {
            EXPECT_EQ(layer_names[j], sectors[j+1].name);
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileLayerMaterials)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        EXPECT_EQ(layer_materials.size(), sectors.size() - 1);
        unsigned int n_layers = std::min(layer_materials.size(), sectors.size() - 1);
        MaterialModel materials = A.GetMaterials();
        for(unsigned int j=0; j<n_layers; ++j) {
            EXPECT_EQ(layer_materials[j], materials.GetMaterialName(sectors[j+1].material_id));
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileSectorTypes)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        unsigned int n_layers = sectors.size()-1;
        ASSERT_EQ(n_layers, layer_thicknesses.size());
        double max_radius = 0.0;
        for(unsigned int j=0; j<n_layers; ++j) {
            unsigned int sector_index = j+1;
            unsigned int layer_index = j;
            DetectorSector & sector = sectors[sector_index];
            std::shared_ptr<const Geometry> geo = sector.geo;
            Placement placement = geo->GetPlacement();
            Sphere const* sphere = dynamic_cast<Sphere const*>(geo.get());
            ASSERT_TRUE(sphere);
            double ice_offset = placement.GetPosition().magnitude();

            if(ice_offset == 0.0) {
                EXPECT_GT(sphere->GetRadius(), max_radius);
                EXPECT_LE(sphere->GetInnerRadius(), max_radius);
            }
            else {
                EXPECT_GT(sphere->GetRadius()+ice_offset, max_radius);
                if(sphere->GetInnerRadius() > 0) {
                    EXPECT_LE(sphere->GetInnerRadius()+ice_offset, max_radius);
                }
            }
            double layer_thickness = layer_thicknesses[layer_index];
            double layer_radius = layer_radii[layer_index];
            if(ice_offset == 0.0) {
                ASSERT_DOUBLE_EQ(layer_radius, sphere->GetRadius());
            }
            else {
                ASSERT_DOUBLE_EQ(layer_radius, sphere->GetRadius() + ice_offset);
            }
            if(ice_offset == 0.0) {
                EXPECT_NEAR(layer_thickness, sphere->GetRadius() - max_radius, std::max(std::abs(layer_thickness), std::abs(sphere->GetRadius() - max_radius))*1e-6);
            }
            else {
                EXPECT_NEAR(layer_thickness, sphere->GetRadius() - max_radius + ice_offset, std::max(std::abs(layer_thickness), std::abs(sphere->GetRadius() - max_radius + ice_offset))*1e-6);
            }
            if(ice_offset == 0.0) {
                max_radius = sphere->GetRadius();
            }
            else {
                max_radius = sphere->GetRadius() + ice_offset;
            }
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileGetMassDensityCachedIntersections)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        double max_radius = *std::max_element(layer_radii.begin(), layer_radii.end());
        max_depth = std::min(max_depth, max_radius);
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        Vector3D p0 = RandomVector(max_radius);
        Vector3D direction = RandomVector(1.0);
        direction.normalize();
        Geometry::IntersectionList intersections = A.GetIntersections(DetectorPosition(p0), DetectorDirection(direction));
        for(unsigned int j=0; j<10; ++j) {
            Vector3D p1 = p0 + direction * (FakeLegacyDetectorModelFile::RandomDouble()*max_radius*2 - max_radius);
            double expect = A.GetMassDensity(DetectorPosition(p1));
            double density = A.GetMassDensity(intersections, DetectorPosition(p1));
            EXPECT_DOUBLE_EQ(density, expect) << i << " " << j;
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileIntegralCachedIntersections)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        double max_radius = *std::max_element(layer_radii.begin(), layer_radii.end());
        max_depth = std::min(max_depth, max_radius);
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        Vector3D p0 = RandomVector(max_radius);
        Vector3D direction = RandomVector(1.0);
        direction.normalize();
        Geometry::IntersectionList intersections = A.GetIntersections(DetectorPosition(p0), DetectorDirection(direction));
        for(unsigned int j=0; j<10; ++j) {
            Vector3D p1 = p0 + direction * (FakeLegacyDetectorModelFile::RandomDouble()*max_radius*2 - max_radius);
            Vector3D p2 = p0 + direction * (FakeLegacyDetectorModelFile::RandomDouble()*max_radius*2 - max_radius);
            double expect = A.GetColumnDepthInCGS(DetectorPosition(p1), DetectorPosition(p2));
            double integral = A.GetColumnDepthInCGS(intersections, DetectorPosition(p1), DetectorPosition(p2));
            EXPECT_NEAR(integral, expect, std::max(std::abs(integral), std::abs(expect))*1e-8);
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileInverseIntegralCachedIntersections)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        double max_radius = *std::max_element(layer_radii.begin(), layer_radii.end());
        max_depth = std::min(max_depth, max_radius);
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = FakeLegacyDetectorModelFile::RandomDouble()*180;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        Vector3D p0 = RandomVector(max_radius);
        Vector3D p1 = RandomVector(max_radius);
        Vector3D direction = p1 - p0;
        direction.normalize();
        double integral = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));
        Geometry::IntersectionList intersections = A.GetIntersections(DetectorPosition(p0), DetectorDirection(direction));
        for(unsigned int j=0; j<10; ++j) {
            double expect = A.DistanceForColumnDepthFromPoint(DetectorPosition(p0), DetectorDirection(direction), integral);
            double distance = A.DistanceForColumnDepthFromPoint(intersections, DetectorPosition(p0), DetectorDirection(direction), integral);
            EXPECT_NEAR(distance, expect, std::max(std::abs(distance), std::abs(expect))*1e-8);
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantIntegralInternal)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        p0 = A.ToDet(GeometryPosition(p0));
        p1 = A.ToDet(GeometryPosition(p1));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));
        EXPECT_DOUBLE_EQ((p1 - p0).magnitude() * 100 * rho, sum);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantGetMassDensityInternal)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        p0 = A.ToDet(GeometryPosition(p0));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_dist = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density_dist);
        double rho = density_dist->Evaluate(Vector3D());
        double density = A.GetMassDensity(DetectorPosition(p0));
        EXPECT_DOUBLE_EQ(density, rho);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantIntegralNested)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        ASSERT_NO_THROW(reset(2, 1));
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(3, sectors.size());
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);
        EXPECT_GE(sphere_1->GetRadius(), sphere_0->GetRadius());
        double max_radius = sphere_1->GetRadius();
        double min_radius = sphere_0->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        GeometryPosition p0_geo(p0);
        GeometryPosition p1_geo(p1);
        DetectorPosition p0_det(A.ToDet(p0_geo));
        DetectorPosition p1_det(A.ToDet(p1_geo));
        double distance = (p1-p0).magnitude();
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        double rho_0 = density_0->Evaluate(Vector3D());
        double rho_1 = density_1->Evaluate(Vector3D());
        ASSERT_LE(p0.magnitude(), max_radius);
        ASSERT_LE(p1.magnitude(), max_radius);
        bool in_0 = p0.magnitude() <= sphere_0->GetRadius();
        bool in_1 = p1.magnitude() <= sphere_0->GetRadius();
        double integral = 0.0;
        double sum = A.GetColumnDepthInCGS(p0_det, p1_det);
        if(in_0 and in_1) {
            integral = rho_0 * (p1-p0).magnitude();
            ASSERT_NEAR(integral * 100, sum, std::min(std::abs(sum), std::abs(integral * 100)) * 1e-8);
        }
        else {
            Vector3D direction = p1_geo - p0_geo;
            direction.normalize();
            std::vector<Geometry::Intersection> intersections = sphere_0->Intersections(p0_geo, direction);
            if(intersections.size() > 1) {
                double dist_0 = intersections[0].distance;
                double dist_1 = intersections[1].distance;
                if(dist_0 > 0 and dist_1 > 0) {
                    double near = std::min(dist_0, dist_1);
                    if((not in_0) and in_1) {
                        integral += rho_1*near;
                        integral += rho_0*((p1-p0).magnitude() - near);
                        ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-12);
                    } else if(in_0 and not in_1) {
                        integral += rho_0*near;
                        integral += rho_1*((p1-p0).magnitude() - near);
                        ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-12);
                    } else if((not in_0) and (not in_1)) {
                        double far = std::max(std::max(dist_0, 0.0), std::max(dist_1, 0.0));
                        if(near < distance) {
                            integral += rho_1*near;
                            if(far < distance) {
                                integral += rho_0*(far-near);
                                integral += rho_1*(distance-far);
                            } else {
                                integral += rho_0*(distance - near);
                            }
                        } else {
                            integral += rho_1*distance;
                        }
                        ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-12);
                    }
                }
                else if (dist_0 <= 0 and dist_1 <= 0) {
                    integral = rho_1 * (p1-p0).magnitude();
                    ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-12);
                }
                else {
                    double dist = std::max(dist_0, dist_1);
                    if((not in_0) and in_1) {
                        assert(false);
                    } else if(in_0 and not in_1) {
                        integral += rho_0*dist;
                        integral += rho_1*((p1-p0).magnitude() - dist);
                        ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-12);
                    } else if((not in_0) and (not in_1)) {
                        assert(false);
                    }
                }
            }
            else {
                integral = rho_1 * (p1-p0).magnitude();
                ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-8);
            }
        }
        ASSERT_NEAR(integral * 100, sum, std::min(std::abs(integral * 100), std::abs(sum)) * 1e-12);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantGetMassDensityNested)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        ASSERT_NO_THROW(reset(2, 1));
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(3, sectors.size());
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);
        EXPECT_GE(sphere_1->GetRadius(), sphere_0->GetRadius());
        double max_radius = sphere_1->GetRadius();
        double min_radius = sphere_0->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        GeometryPosition p0_geo(p0);
        DetectorPosition p0_det(A.ToDet(p0_geo));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        double rho_0 = density_0->Evaluate(Vector3D());
        double rho_1 = density_1->Evaluate(Vector3D());
        ASSERT_LE(min_radius, max_radius);
        ASSERT_LE(p0_geo->magnitude(), max_radius);
        bool in_0 = p0_geo->magnitude() <= sphere_0->GetRadius();
        double density = A.GetMassDensity(p0_det);
        if(in_0) {
            ASSERT_DOUBLE_EQ(density, rho_0);
        }
        else {
            ASSERT_DOUBLE_EQ(density, rho_1);
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantIntegralIntersecting)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }
        double radius = FakeLegacyDetectorModelFile::RandomDouble()*1000;

        DetectorSector upper_sector;
        Vector3D upper_center(-radius/4.0,0,0);
        upper_sector.name = "upper";
        upper_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        upper_sector.level = -1;
        upper_sector.geo = Sphere(upper_center, radius, 0).create();
        upper_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector lower_sector;
        Vector3D lower_center(radius/4.0,0,0);
        lower_sector.name = "lower";
        lower_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        lower_sector.level = -2;
        lower_sector.geo = Sphere(lower_center, radius, 0).create();
        lower_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(lower_sector);
        A.AddSector(upper_sector);

        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(3, sectors.size());
        DetectorSector sector_vacuum = sectors[0];
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        ASSERT_EQ(lower_sector.name, sector_0.name);
        ASSERT_EQ(upper_sector.name, sector_1.name);
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);

        Vector3D p0 = RandomVector(0, radius*2.0);
        Vector3D p1 = RandomVector(0, radius*2.0);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_vacuum = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_vacuum.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_vacuum);
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        double rho_vacuum = density_vacuum->Evaluate(Vector3D());
        double rho_lower = density_0->Evaluate(lower_center);
        double rho_upper = density_1->Evaluate(upper_center);

        double integral = 0.0;
        double sum = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));

        std::vector<Geometry::Intersection> lower_intersections = sphere_0->Intersections(p0, direction);
        std::vector<Geometry::Intersection> upper_intersections = sphere_1->Intersections(p0, direction);

        bool p0_in_lower = (p0-lower_center).magnitude() < radius;
        bool p0_in_upper = (p0-upper_center).magnitude() < radius;
        bool p1_in_lower = (p1-lower_center).magnitude() < radius;
        bool p1_in_upper = (p1-upper_center).magnitude() < radius;

        if(p0_in_upper) { // Start in the full sphere
            if(p1_in_upper) { // End in full sphere
                integral += rho_upper*distance;
            } else if(p1_in_lower) { // End in the partial sphere
                int upper_exit_index = 0;
                if(upper_intersections[0].distance < upper_intersections[1].distance) {
                    upper_exit_index = 1;
                }
                double dist_in_upper = upper_intersections[upper_exit_index].distance;
                integral += rho_upper*dist_in_upper;
                if((upper_intersections[upper_exit_index].position - lower_center).magnitude() < radius) { // Transition directly between sectors
                    double dist_in_lower = distance - dist_in_upper;
                    integral += dist_in_lower*rho_lower;
                } else { // Go through vacuum to transition between sectors
                    int lower_entry_index = 0;
                    if(lower_intersections[0].distance > lower_intersections[1].distance) {
                        lower_entry_index = 1;
                    }
                    double dist_in_lower = (p1 - lower_intersections[lower_entry_index].position).magnitude();
                    integral += rho_lower*dist_in_lower;
                    double dist_in_vacuum = (distance - (dist_in_upper + dist_in_lower));
                    integral += rho_vacuum*dist_in_vacuum;
                }

            } else { // End in vacuum
                int upper_exit_index = 0;
                if(upper_intersections[0].distance < upper_intersections[1].distance) {
                    upper_exit_index = 1;
                }
                double dist_in_upper = upper_intersections[upper_exit_index].distance;
                integral += rho_upper*dist_in_upper;
                if((upper_intersections[upper_exit_index].position - lower_center).magnitude() < radius) { // Transition first to partial sphere
                    int lower_exit_index = 0;
                    if(lower_intersections[0].distance < lower_intersections[1].distance) {
                        lower_exit_index = 1;
                    }
                    double dist_in_lower = (lower_intersections[lower_exit_index].distance - dist_in_upper);
                    integral += rho_lower*dist_in_lower;
                    double dist_in_vacuum = distance - (dist_in_upper - dist_in_lower);
                    integral += rho_vacuum*dist_in_vacuum;
                } else { // Go straight to vacuum
                    double dist_in_vacuum = distance - dist_in_upper;
                    integral += rho_vacuum*dist_in_vacuum;
                }
            }
        } else if(p0_in_lower) { // Start in the partial sphere
            if(p1_in_lower and not p1_in_upper) { // End in the partial sphere
                double u_dist_0 = upper_intersections.size() > 1 ? upper_intersections[0].distance : 0.0;
                double u_dist_1 = upper_intersections.size() > 1 ? upper_intersections[1].distance : 0.0;
                if((u_dist_0 > 0 and u_dist_0 < distance) or (u_dist_1 > 0 and u_dist_1 < distance)) { // Passes through full sphere
                    double distance_in_upper = std::abs(u_dist_0 - u_dist_1);
                    double distance_in_lower = distance - distance_in_upper;
                    integral += rho_lower*distance_in_lower;
                    integral += rho_upper*distance_in_upper;
                } else { // Entire path in partial sphere
                    integral += rho_lower*distance;
                }

            } else if(p1_in_upper) { // End in the full sphere
                int upper_entry_index = 0;
                if(upper_intersections[0].distance > upper_intersections[1].distance) {
                    upper_entry_index = 1;
                }
                if((upper_intersections[upper_entry_index].position - lower_center).magnitude() < radius) { // Transition straight to full sphere
                    double dist_in_lower = upper_intersections[upper_entry_index].distance;
                    double dist_in_upper = distance - dist_in_lower;
                    integral += rho_lower*dist_in_lower;
                    integral += rho_upper*dist_in_upper;
                } else { // Go through vacuum to transition between sectors
                    double dist_in_upper = distance - upper_intersections[upper_entry_index].distance;
                    integral += rho_upper*dist_in_upper;
                    int lower_exit_index = 0;
                    if(lower_intersections[0].distance < lower_intersections[1].distance) {
                        lower_exit_index = 1;
                    }
                    double dist_in_lower = lower_intersections[lower_exit_index].distance;
                    integral += rho_lower*dist_in_lower;
                    double dist_in_vacuum = (distance - (dist_in_upper + dist_in_lower));
                    integral += rho_vacuum*dist_in_vacuum;
                }
            } else { // End in vacuum
                double u_dist_0 = upper_intersections.size() > 1 ? upper_intersections[0].distance : 0.0;
                double u_dist_1 = upper_intersections.size() > 1 ? upper_intersections[1].distance : 0.0;
                if((u_dist_0 > 0 and u_dist_0 < distance) or (u_dist_1 > 0 and u_dist_1 < distance)) { // Passes through full sphere
                    double distance_in_upper = std::abs(u_dist_0 - u_dist_1);
                    double lower_exit_distance = std::max(lower_intersections[0].distance, lower_intersections[1].distance);
                    double distance_in_lower = std::min(u_dist_0, u_dist_1);
                    if(lower_exit_distance > u_dist_0 and lower_exit_distance > u_dist_1) {
                        distance_in_lower = lower_exit_distance - distance_in_upper;
                    }
                    double distance_in_vacuum = distance - (distance_in_lower + distance_in_upper);
                    integral += rho_lower*distance_in_lower;
                    integral += rho_upper*distance_in_upper;
                    integral += rho_vacuum*distance_in_vacuum;
                } else { // Goes straight to vacuum
                    double dist_in_lower = std::max(lower_intersections[0].distance, lower_intersections[1].distance);
                    double dist_in_vacuum = distance - dist_in_lower;
                    integral += rho_lower*dist_in_lower;
                    integral += rho_vacuum*dist_in_vacuum;
                }

            }
        } else { // Start in vacuum
            double u_dist_0 = upper_intersections.size() > 1 ? upper_intersections[0].distance : 0.0;
            double u_dist_1 = upper_intersections.size() > 1 ? upper_intersections[1].distance : 0.0;
            double l_dist_0 = lower_intersections.size() > 1 ? lower_intersections[0].distance : 0.0;
            double l_dist_1 = lower_intersections.size() > 1 ? lower_intersections[1].distance : 0.0;
            bool hits_upper = (u_dist_0 > 0 and u_dist_0 < distance) or (u_dist_1 > 0 and u_dist_1 < distance);
            bool hits_lower = (l_dist_0 > 0 and l_dist_0 < distance) or (l_dist_1 > 0 and l_dist_1 < distance);
            if((not p1_in_upper) and (not p1_in_lower)) { // End in vacuum
                if((not hits_upper) and (not hits_lower)) { // Full path is in vacuum
                    integral += rho_vacuum*distance;
                } else if(hits_upper and (not hits_lower)) {
                    double dist_in_upper = std::abs(u_dist_1 - u_dist_0);
                    double dist_in_vacuum = distance - dist_in_upper;
                    integral += rho_vacuum*dist_in_vacuum;
                    integral += rho_upper*dist_in_upper;
                } else if(hits_lower and (not hits_upper)) {
                    double dist_in_lower = std::abs(l_dist_1 - l_dist_0);
                    double dist_in_vacuum = distance - dist_in_lower;
                    integral += rho_vacuum*dist_in_vacuum;
                    integral += rho_lower*dist_in_lower;
                } else if(hits_upper and hits_lower) {
                    double lower_min = std::min(l_dist_0, l_dist_1);
                    double lower_max = std::max(l_dist_0, l_dist_1);
                    double upper_min = std::min(u_dist_0, u_dist_1);
                    double upper_max = std::max(u_dist_0, u_dist_1);
                    if(lower_min < upper_min and upper_min < lower_max and lower_max < upper_max) {
                        // l (L) u (U) \l u
                        integral += rho_vacuum*lower_min;
                        integral += rho_lower*(upper_min - lower_min);
                        integral += rho_upper*(upper_max - upper_min);
                        integral += rho_vacuum*(distance - upper_max);
                    } else if(upper_min < lower_min and lower_min < upper_max and upper_max < lower_max) {
                        // u (U) \l u (L) l
                        integral += rho_vacuum*upper_min;
                        integral += rho_upper*(upper_max - upper_min);
                        integral += rho_lower*(lower_max - upper_max);
                        integral += rho_vacuum*(distance - lower_max);

                    } else if (upper_min < lower_min and lower_max < upper_max) {
                        // u \l (U) \l u
                        integral += rho_vacuum*upper_min;
                        integral += rho_upper*(upper_max - upper_min);
                        integral += rho_vacuum*(distance - upper_max);
                    } else if (lower_min < upper_min and upper_max < lower_max) {
                        // l u u l
                        integral += rho_vacuum*lower_min;
                        integral += rho_lower*(upper_min - lower_min);
                        integral += rho_upper*(upper_max - upper_min);
                        integral += rho_lower*(lower_max - upper_max);
                        integral += rho_vacuum*(distance - lower_max);
                    }
                }
            } else if(p1_in_upper) { // Ends in full sphere
                if((not hits_upper) and (not hits_lower)) { // Full path is in vacuum
                    assert(false);
                } else if(hits_upper and (not hits_lower)) {
                    double dist_in_vacuum = std::min(u_dist_0, u_dist_1);
                    double dist_in_upper = distance - dist_in_vacuum;
                    integral += rho_vacuum*dist_in_vacuum;
                    integral += rho_upper*dist_in_upper;
                } else if(hits_lower and (not hits_upper)) {
                    assert(false);
                } else if(hits_upper and hits_lower) {
                    double lower_min = std::min(l_dist_0, l_dist_1);
                    double lower_max = std::max(l_dist_0, l_dist_1);
                    double upper_min = std::min(u_dist_0, u_dist_1);
                    if(lower_min < upper_min and upper_min < lower_max) {
                        // l (L) u (U) \l u
                        integral += rho_vacuum*lower_min;
                        integral += rho_lower*(upper_min - lower_min);
                        integral += rho_upper*(distance - upper_min);
                    } else if(upper_min < lower_min) {
                        integral += rho_vacuum*upper_min;
                        integral += rho_upper*(distance - upper_min);
                    } else if (lower_min < upper_min) {
                        // l u u l
                        integral += rho_vacuum*lower_min;
                        integral += rho_lower*(upper_min - lower_min);
                        integral += rho_upper*(distance - upper_min);
                    }
                }
            } else if(p1_in_lower) { // Ends in partial sphere
                if((not hits_upper) and (not hits_lower)) { // Full path is in vacuum
                    assert(false);
                } else if(hits_upper and (not hits_lower)) {
                    assert(false);
                } else if(hits_lower and (not hits_upper)) {
                    double dist_in_vacuum = std::min(l_dist_0, l_dist_1);
                    double dist_in_lower = distance - dist_in_vacuum;
                    integral += rho_vacuum*dist_in_vacuum;
                    integral += rho_lower*dist_in_lower;
                } else if(hits_upper and hits_lower) {
                    double lower_min = std::min(l_dist_0, l_dist_1);
                    double lower_max = std::max(l_dist_0, l_dist_1);
                    double upper_min = std::min(u_dist_0, u_dist_1);
                    double upper_max = std::max(u_dist_0, u_dist_1);
                    if(lower_min < upper_min and upper_min < lower_max and lower_max < upper_max) {
                        assert(false);
                    } else if(upper_min < lower_min and lower_min < upper_max and upper_max < lower_max) {
                        // u (U) \l u (L) l
                        integral += rho_vacuum*upper_min;
                        integral += rho_upper*(upper_max - upper_min);
                        integral += rho_lower*(distance - upper_max);

                    } else if (upper_min < lower_min and lower_max < upper_max) {
                        // u \l (U) \l u
                        assert(false);
                    } else if (lower_min < upper_min and upper_max < lower_max) {
                        // l u u l
                        integral += rho_vacuum*lower_min;
                        integral += rho_lower*(upper_min - lower_min);
                        integral += rho_upper*(upper_max - upper_min);
                        integral += rho_lower*(distance - upper_max);
                    }
                }
            }
        }
        EXPECT_DOUBLE_EQ(integral * 100, sum);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantGetMassDensityIntersecting)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }
        double radius = FakeLegacyDetectorModelFile::RandomDouble()*1000;

        DetectorSector upper_sector;
        Vector3D upper_center(-radius/4.0,0,0);
        upper_sector.name = "upper";
        upper_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        upper_sector.level = -1;
        upper_sector.geo = Sphere(upper_center, radius, 0).create();
        upper_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector lower_sector;
        Vector3D lower_center(radius/4.0,0,0);
        lower_sector.name = "lower";
        lower_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        lower_sector.level = -2;
        lower_sector.geo = Sphere(lower_center, radius, 0).create();
        lower_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(lower_sector);
        A.AddSector(upper_sector);

        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(3, sectors.size());
        DetectorSector sector_vacuum = sectors[0];
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        ASSERT_EQ(lower_sector.name, sector_0.name);
        ASSERT_EQ(upper_sector.name, sector_1.name);
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);

        Vector3D p0 = RandomVector(0, radius*2.0);
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_vacuum = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_vacuum.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_vacuum);
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        double rho_vacuum = density_vacuum->Evaluate(Vector3D());
        double rho_lower = density_0->Evaluate(lower_center);
        double rho_upper = density_1->Evaluate(upper_center);

        double density = A.GetMassDensity(DetectorPosition(p0));

        bool p0_in_lower = (p0-lower_center).magnitude() < radius;
        bool p0_in_upper = (p0-upper_center).magnitude() < radius;

        if(p0_in_upper) { // Start in the full sphere
            EXPECT_DOUBLE_EQ(density, rho_upper);
        } else if(p0_in_lower) { // Start in the partial sphere
            EXPECT_DOUBLE_EQ(density, rho_lower);
        } else { // Start in vacuum
            EXPECT_DOUBLE_EQ(density, rho_vacuum);
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantIntegralHidden)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }
        double radius = FakeLegacyDetectorModelFile::RandomDouble()*1000;

        DetectorSector upper_sector;
        Vector3D upper_center(0,0,0);
        upper_sector.name = "upper";
        upper_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        upper_sector.level = -1;
        upper_sector.geo = Sphere(upper_center, radius, 0).create();
        upper_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector lower_sector;
        Vector3D lower_center(0,0,0);
        lower_sector.name = "lower";
        lower_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        lower_sector.level = -2;
        lower_sector.geo = Sphere(lower_center, radius/2.0, 0).create();
        lower_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(lower_sector);
        A.AddSector(upper_sector);

        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(3, sectors.size());
        DetectorSector sector_vacuum = sectors[0];
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        ASSERT_EQ(lower_sector.name, sector_0.name);
        ASSERT_EQ(upper_sector.name, sector_1.name);
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);

        Vector3D p0 = RandomVector(0, radius*2.0);
        Vector3D p1 = RandomVector(0, radius*2.0);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_vacuum = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_vacuum.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_vacuum);
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        double rho_vacuum = density_vacuum->Evaluate(Vector3D());
        double rho_lower = density_0->Evaluate(lower_center);
        double rho_upper = density_1->Evaluate(upper_center);

        std::vector<Geometry::Intersection> intersections = sphere_1->Intersections(p0, direction);

        bool p0_in_upper = (p0 - upper_center).magnitude() < radius;
        bool p1_in_upper = (p1 - upper_center).magnitude() < radius;

        double sum = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));
        double integral = 0.0;
        if(p0_in_upper) {
            if(p1_in_upper) {
                integral += rho_upper*distance;
                ASSERT_DOUBLE_EQ(integral * 100, sum) << sum << " " << rho_upper*distance << " " << rho_lower*distance;
            } else {
                double dist_in_upper = std::max(intersections[0].distance, intersections[1].distance);
                integral += rho_upper*dist_in_upper;
                integral += rho_vacuum*(distance - dist_in_upper);
                ASSERT_DOUBLE_EQ(integral * 100, sum);
            }
        } else {
            if(p1_in_upper) {
                double dist_in_vacuum = std::min(intersections[0].distance, intersections[1].distance);
                integral += rho_vacuum*dist_in_vacuum;
                integral += rho_upper*(distance - dist_in_vacuum);
                    ASSERT_DOUBLE_EQ(integral * 100, sum);
            } else {
                if(intersections.size() > 1) {
                    if(intersections[0].distance > 0 and intersections[1].distance > 0 and distance > std::max(intersections[0].distance, intersections[1].distance)) {
                        double dist_in_upper = std::abs(intersections[1].distance - intersections[0].distance);
                        integral += rho_upper*dist_in_upper;
                        integral += rho_vacuum*(distance - dist_in_upper);
                        ASSERT_DOUBLE_EQ(integral * 100, sum);
                    } else {
                        integral += rho_vacuum*distance;
                        ASSERT_DOUBLE_EQ(integral * 100, sum);
                    }
                } else {
                    integral += rho_vacuum*distance;
                    ASSERT_DOUBLE_EQ(integral * 100, sum);
                }
            }
        }
        EXPECT_DOUBLE_EQ(integral * 100, sum);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantGetMassDensityHidden)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }
        double radius = FakeLegacyDetectorModelFile::RandomDouble()*1000;

        DetectorSector upper_sector;
        Vector3D upper_center(0,0,0);
        upper_sector.name = "upper";
        upper_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        upper_sector.level = -1;
        upper_sector.geo = Sphere(upper_center, radius, 0).create();
        upper_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector lower_sector;
        Vector3D lower_center(0,0,0);
        lower_sector.name = "lower";
        lower_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        lower_sector.level = -2;
        lower_sector.geo = Sphere(lower_center, radius/2.0, 0).create();
        lower_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(lower_sector);
        A.AddSector(upper_sector);

        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(3, sectors.size());
        DetectorSector sector_vacuum = sectors[0];
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        ASSERT_EQ(lower_sector.name, sector_0.name);
        ASSERT_EQ(upper_sector.name, sector_1.name);
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);

        Vector3D p0 = RandomVector(0, radius*2.0);
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_vacuum = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_vacuum.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_vacuum);
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        double rho_vacuum = density_vacuum->Evaluate(Vector3D());
        double rho_upper = density_1->Evaluate(upper_center);

        bool p0_in_upper = (p0 - upper_center).magnitude() < radius;

        double density = A.GetMassDensity(DetectorPosition(p0));
        if(p0_in_upper) {
            EXPECT_DOUBLE_EQ(density, rho_upper);
        } else {
            EXPECT_DOUBLE_EQ(density, rho_vacuum);
        }
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantInverseIntegralInternal)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A.GetSectors();
        unsigned int n_layers = std::min(layer_names.size(), sectors.size() - 1);
        double max_radius = 0.0;
        for(unsigned int j=0; j<n_layers; ++j) {
            unsigned int sector_index = j+1;
            std::shared_ptr<const Geometry> geo = sectors[sector_index].geo;
            Sphere const* sphere = dynamic_cast<Sphere const*>(geo.get());
            max_radius = std::max(max_radius, sphere->GetRadius());
        }
        double min_radius = 0.0;

        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        p0 = A.ToDet(GeometryPosition(p0));
        p1 = A.ToDet(GeometryPosition(p1));

        double sum = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();

        double found_distance = A.DistanceForColumnDepthFromPoint(DetectorPosition(p0), DetectorDirection(direction), sum);

        EXPECT_NEAR(distance, found_distance, std::max(std::abs(distance), std::abs(found_distance))*1e-12);
    }
}
TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantInverseIntegralHidden)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }
        double radius = FakeLegacyDetectorModelFile::RandomDouble()*1000;

        DetectorSector upper_sector;
        Vector3D upper_center(0,0,0);
        upper_sector.name = "upper";
        upper_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        upper_sector.level = -1;
        upper_sector.geo = Sphere(upper_center, radius, 0).create();
        upper_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector lower_sector;
        Vector3D lower_center(0,0,0);
        lower_sector.name = "lower";
        lower_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        lower_sector.level = -2;
        lower_sector.geo = Sphere(lower_center, radius/2.0, 0).create();
        lower_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector bg_sector;
        Vector3D bg_center(0,0,0);
        bg_sector.name = "bg";
        bg_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        bg_sector.level = -100;
        bg_sector.geo = Sphere(bg_center, radius*100, 0).create();
        bg_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(lower_sector);
        A.AddSector(upper_sector);
        A.AddSector(bg_sector);

        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(4, sectors.size());
        DetectorSector sector_vacuum = sectors[0];
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        ASSERT_EQ(lower_sector.name, sector_0.name);
        ASSERT_EQ(upper_sector.name, sector_1.name);
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);

        Vector3D p0 = RandomVector(0, radius*2.0);
        Vector3D p1 = RandomVector(0, radius*2.0);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_vacuum = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_vacuum.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_vacuum);
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);

        double integral = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));

        double found_distance = A.DistanceForColumnDepthFromPoint(DetectorPosition(p0), DetectorDirection(direction), integral);
        EXPECT_NEAR(distance, found_distance, std::max(std::abs(distance), std::abs(found_distance))*1e-12);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantInverseIntegralIntersecting)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }
        double radius = FakeLegacyDetectorModelFile::RandomDouble()*1000;

        DetectorSector upper_sector;
        Vector3D upper_center(-radius/4.0,0,0);
        upper_sector.name = "upper";
        upper_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        upper_sector.level = -1;
        upper_sector.geo = Sphere(upper_center, radius, 0).create();
        upper_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector lower_sector;
        Vector3D lower_center(radius/4.0,0,0);
        lower_sector.name = "lower";
        lower_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        lower_sector.level = -2;
        lower_sector.geo = Sphere(lower_center, radius, 0).create();
        lower_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        DetectorSector bg_sector;
        Vector3D bg_center(0,0,0);
        bg_sector.name = "bg";
        bg_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        bg_sector.level = -100;
        bg_sector.geo = Sphere(bg_center, radius*100, 0).create();
        bg_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(lower_sector);
        A.AddSector(upper_sector);
        A.AddSector(bg_sector);

        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(4, sectors.size());
        DetectorSector sector_vacuum = sectors[0];
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        ASSERT_EQ(lower_sector.name, sector_0.name);
        ASSERT_EQ(upper_sector.name, sector_1.name);
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);

        Vector3D p0 = RandomVector(0, radius*2.0);
        Vector3D p1 = RandomVector(0, radius*2.0);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_vacuum = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_vacuum.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_vacuum);
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);

        double integral = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));

        double found_distance = A.DistanceForColumnDepthFromPoint(DetectorPosition(p0), DetectorDirection(direction), integral);
        EXPECT_NEAR(distance, found_distance, std::max(std::abs(distance), std::abs(found_distance))*1e-12);
    }
}

TEST_F(FakeLegacyDetectorModelTest, LegacyFileConstantInverseIntegralNested)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        reset();
        ASSERT_NO_THROW(reset(2, 1));
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));

        MaterialModel materials = A.GetMaterials();
        int material_count = 0;
        while(materials.HasMaterial(material_count)) {
            material_count += 1;
        }

        double max_depth = 5000;
        double max_radius = *std::max_element(layer_radii.begin(), layer_radii.end());
        max_depth = std::min(max_depth, max_radius);

        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A.LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));

        DetectorSector bg_sector;
        Vector3D bg_center(0,0,0);
        bg_sector.name = "bg";
        bg_sector.material_id = FakeLegacyDetectorModelFile::RandomDouble()*material_count;
        bg_sector.level = -100;
        bg_sector.geo = Sphere(bg_center, max_radius*100, 0).create();
        bg_sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(FakeLegacyDetectorModelFile::RandomDouble()*15).create();

        A.AddSector(bg_sector);
        std::vector<DetectorSector> sectors = A.GetSectors();
        ASSERT_EQ(4, sectors.size());
        DetectorSector sector_0 = sectors[1];
        DetectorSector sector_1 = sectors[2];
        Sphere const * sphere_0 = dynamic_cast<Sphere const *>(sector_0.geo.get());
        Sphere const * sphere_1 = dynamic_cast<Sphere const *>(sector_1.geo.get());
        ASSERT_TRUE(sphere_0);
        ASSERT_TRUE(sphere_1);
        EXPECT_GE(sphere_1->GetRadius(), sphere_0->GetRadius());
        max_radius = sphere_1->GetRadius();
        double min_radius = sphere_0->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_0 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_0.density.get());
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density_1 = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector_1.density.get());
        ASSERT_TRUE(density_0);
        ASSERT_TRUE(density_1);
        ASSERT_LE(p0.magnitude(), max_radius);
        ASSERT_LE(p1.magnitude(), max_radius);
        double integral = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1));

        double found_distance = A.DistanceForColumnDepthFromPoint(DetectorPosition(p0), DetectorDirection(direction), integral);
        EXPECT_NEAR(distance, found_distance, std::max(std::abs(distance), std::abs(found_distance))*1e-12);
    }
}

TEST_F(FakeDetectorModelTest, FileLoad)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset());
        DetectorModel A;
        ASSERT_NO_THROW(A.LoadMaterialModel(materials_file));
        EXPECT_NO_THROW(A.LoadDetectorModel(model_file));
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

