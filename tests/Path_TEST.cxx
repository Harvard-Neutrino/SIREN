
#include <math.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <gtest/gtest.h>

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"

#include "LeptonInjector/Particle.h"

#include "FakeMaterialModel.h"
#include "FakeEarthModel.h"

using namespace earthmodel;

TEST(DefaultConstructor, NoThrow)
{
    EXPECT_NO_THROW(Path());
}

TEST(DefaultConstructor, HasNone)
{
    Path A;
    EXPECT_FALSE(A.HasEarthModel());
    EXPECT_FALSE(A.HasPoints());
    EXPECT_FALSE(A.HasIntersections());
}

TEST(DefaultConstructor, MembersAreDefault)
{
    Path A;
    EXPECT_EQ(std::shared_ptr<const EarthModel>(), A.GetEarthModel());
    EXPECT_EQ(Vector3D(), A.GetFirstPoint());
    EXPECT_EQ(Vector3D(), A.GetLastPoint());
    EXPECT_EQ(Vector3D(), A.GetDirection());
    EXPECT_EQ(0, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(EarthModelConstructor, NoThrow)
{
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    EXPECT_NO_THROW(Path A(EMp));
}

TEST(EarthModelConstructor, MemberValues)
{
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Path A(EMp);
    EXPECT_EQ(EMp, A.GetEarthModel());
    EXPECT_EQ(Vector3D(), A.GetFirstPoint());
    EXPECT_EQ(Vector3D(), A.GetLastPoint());
    EXPECT_EQ(Vector3D(), A.GetDirection());
    EXPECT_EQ(0, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(PointsConstructor, NoThrow)
{
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    EXPECT_NO_THROW(Path A(EMp, B, C));
}

TEST(PointsConstructor, MemberValues)
{
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    Path A(EMp, B, C);
    EXPECT_EQ(EMp, A.GetEarthModel());
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(RayConstructor, NoThrow)
{
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    EXPECT_NO_THROW(Path A(EMp, B, direction, distance));
}

TEST(RayConstructor, MemberValues)
{
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    Path A(EMp, B, direction, distance);
    EXPECT_EQ(EMp, A.GetEarthModel());
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(B + direction * distance, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(EarthModel, SetGet) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Path A;
    A.SetEarthModel(EMp);
    EXPECT_EQ(EMp, A.GetEarthModel());
}

TEST(EnsureEarthModel, Throw) {
    Path A;
    EXPECT_THROW(A.EnsureEarthModel(), std::runtime_error);
}

TEST(EnsureEarthModel, NoThrow) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Path A(EMp);
    EXPECT_NO_THROW(A.EnsureEarthModel());
}

TEST(Points, SetGet) {
    Vector3D A(1,2,3);
    Vector3D B(4,6,8);
    Vector3D direction = B - A;
    double distance = direction.magnitude();
    direction.normalize();
    Path C;
    C.SetPoints(A, B);
    EXPECT_EQ(A, C.GetFirstPoint());
    EXPECT_EQ(B, C.GetLastPoint());
    EXPECT_EQ(direction, C.GetDirection());
    EXPECT_EQ(distance, C.GetDistance());
}

TEST(Points, SetGetRay) {
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    Path A;
    A.SetPointsWithRay(B, direction, distance);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(B + direction * distance, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST(EnsurePoints, Throw) {
    Path A;
    EXPECT_THROW(A.EnsurePoints(), std::runtime_error);
}

TEST(EnsurePoints, NoThrow) {
    Vector3D A(1,2,3);
    Vector3D B(4,6,8);
    Vector3D direction = B - A;
    double distance = direction.magnitude();
    direction.normalize();
    Path C;
    C.SetPoints(A, B);
    EXPECT_NO_THROW(C.EnsurePoints());
}

TEST(EnsurePoints, NoThrowRay) {
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    direction.normalize();
    Path A;
    A.SetPointsWithRay(B, direction, distance);
    EXPECT_NO_THROW(A.EnsurePoints());
}

TEST(EnsureIntersections, Throw) {
    Path A;
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    A = Path(EMp);
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    A = Path();
    A.SetPoints(B, C);
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    A = Path();
    A.SetPointsWithRay(B, direction, distance);
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);
}

TEST(EnsureIntersections, NoThrow) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A;

    A = Path();
    A.SetEarthModel(EMp);
    A.SetPoints(B, C);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path();
    A.SetEarthModel(EMp);
    A.SetPointsWithRay(B, direction, distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp);
    A.SetPoints(B, C);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp);
    A.SetPointsWithRay(B, direction, distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp, B, C);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp, B, direction, distance);
    EXPECT_NO_THROW(A.EnsureIntersections());
}

TEST(PointManipulation, Flip) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    Vector3D reverse = -direction;
    double distance = direction.magnitude();
    direction.normalize();
    reverse.normalize();

    Path A(EMp, B, C);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());

    A.Flip();
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(reverse, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());

    A.Flip();
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST(PointManipulation, ExtendFromEndByDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, B, C);

    double extra_distance = 101;
    A.ExtendFromEndByDistance(extra_distance);

    Vector3D end = C + direction * extra_distance;
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = -1;
    end = C + direction * extra_distance;
    A.ExtendFromEndByDistance(extra_distance);
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = -101;
    end = C + direction * extra_distance;
    A.ExtendFromEndByDistance(extra_distance);
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ExtendFromEndByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance/3.0, extra_distance);
        Vector3D end = inner_p1 + direction * extra_distance;
        EXPECT_EQ(end, P.GetLastPoint());
        EXPECT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance/3.0, extra_distance);
        end = inner_p1 + direction * extra_distance;
        EXPECT_EQ(end, P.GetLastPoint());
        EXPECT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance*1.5, extra_distance);
        end = inner_p1 + direction * extra_distance;
        EXPECT_EQ(inner_p0, P.GetLastPoint());
        EXPECT_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ExtendFromStartByDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, B, C);

    double extra_distance = 101;
    A.ExtendFromStartByDistance(extra_distance);
    Vector3D end = B + reverse * extra_distance;
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = -1;
    end = B + reverse * extra_distance;
    A.ExtendFromStartByDistance(extra_distance);
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = -101;
    end = B + reverse * extra_distance;
    A.ExtendFromStartByDistance(extra_distance);
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ExtendFromStartByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, -direction, extra_column_depth);
        ASSERT_DOUBLE_EQ(distance/3.0, extra_distance);
        ASSERT_DOUBLE_EQ(A->DistanceForColumnDepthFromPoint(P.GetFirstPoint(), -direction, extra_column_depth), A->DistanceForColumnDepthFromPoint(P.GetIntersections(), P.GetFirstPoint(), -direction, extra_column_depth));
        ASSERT_DOUBLE_EQ(extra_distance, P.GetDistanceFromStartInReverse(extra_column_depth));
        P.ExtendFromStartByColumnDepth(extra_column_depth);
        Vector3D end = inner_p0 - direction * extra_distance;
        ASSERT_EQ(end, P.GetFirstPoint());
        ASSERT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance/3.0, extra_distance);
        end = inner_p0 - direction * extra_distance;
        EXPECT_EQ(end, P.GetFirstPoint());
        EXPECT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance*1.5, extra_distance);
        end = inner_p0 - direction * extra_distance;
        EXPECT_EQ(inner_p1, P.GetFirstPoint());
        EXPECT_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ShrinkFromEndByDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, B, C);

    double extra_distance = -101;
    A.ShrinkFromEndByDistance(extra_distance);

    Vector3D end = C - direction * extra_distance;
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = 1;
    end = C - direction * extra_distance;
    A.ShrinkFromEndByDistance(extra_distance);
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = 101;
    end = C - direction * extra_distance;
    A.ShrinkFromEndByDistance(extra_distance);
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ShrinkFromEndByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance/3.0, extra_distance);
        Vector3D end = inner_p1 - direction * extra_distance;
        EXPECT_EQ(end, P.GetLastPoint());
        EXPECT_EQ(distance - extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance/3.0, extra_distance);
        end = inner_p1 - direction * extra_distance;
        EXPECT_EQ(end, P.GetLastPoint());
        EXPECT_EQ(distance - extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance*1.5, extra_distance);
        end = inner_p1 - direction * extra_distance;
        EXPECT_EQ(inner_p0, P.GetLastPoint());
        EXPECT_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ShrinkFromStartByDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, B, C);

    double extra_distance = -101;
    A.ShrinkFromStartByDistance(extra_distance);
    Vector3D end = B - reverse * extra_distance;
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = 1;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartByDistance(extra_distance);
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    extra_distance = 101;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartByDistance(extra_distance);
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ShrinkFromStartByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, extra_column_depth);
        ASSERT_DOUBLE_EQ(distance/3.0, extra_distance);
        ASSERT_DOUBLE_EQ(extra_distance, P.GetDistanceFromStartInReverse(extra_column_depth));
        P.ShrinkFromStartByColumnDepth(extra_column_depth);
        Vector3D end = inner_p0 + direction * extra_distance;
        ASSERT_EQ(end, P.GetFirstPoint());
        ASSERT_EQ(distance - extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance/3.0, extra_distance);
        end = inner_p0 + direction * extra_distance;
        EXPECT_EQ(end, P.GetFirstPoint());
        EXPECT_EQ(distance - extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance*1.5, extra_distance);
        end = inner_p0 + direction * extra_distance;
        EXPECT_EQ(inner_p1, P.GetFirstPoint());
        EXPECT_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ExtendFromEndToDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, B, C);

    double target_distance = 101;
    double extra_distance = target_distance - distance;
    A.ExtendFromEndToDistance(target_distance);

    Vector3D end = C + direction * extra_distance;
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = 1;
    extra_distance = target_distance - distance;
    end = C + direction * extra_distance;
    A.ExtendFromEndToDistance(target_distance);
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = -101;
    extra_distance = target_distance - distance;
    end = C + direction * extra_distance;
    A.ExtendFromEndToDistance(target_distance);
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ExtendFromEndToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        double target_distance = distance + extra_distance;
        double target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, target_column_depth);
        ASSERT_DOUBLE_EQ(distance + distance/3.0, target_distance);
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistanceFromStartAlongPath(target_column_depth));
        P.ExtendFromEndToColumnDepth(target_column_depth);
        Vector3D end = inner_p0 + direction * target_distance;
        ASSERT_TRUE((end - P.GetLastPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetLastPoint().magnitude()));
        ASSERT_DOUBLE_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, target_column_depth);
        P.ExtendFromEndToColumnDepth(target_column_depth);
        ASSERT_DOUBLE_EQ(distance - distance/3.0, target_distance);
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistanceFromStartAlongPath(target_column_depth));
        end = inner_p0 + direction * target_distance;
        ASSERT_EQ(inner_p1, P.GetLastPoint());
        ASSERT_DOUBLE_EQ(distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, target_column_depth);
        P.ExtendFromEndToColumnDepth(target_column_depth);
        ASSERT_DOUBLE_EQ(distance - distance*1.5, target_distance);
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistanceFromStartAlongPath(target_column_depth));
        end = inner_p0 + direction * target_distance;
        ASSERT_EQ(inner_p1, P.GetLastPoint());
        ASSERT_DOUBLE_EQ(distance, P.GetDistance());
    }
}


TEST(PointManipulation, ExtendFromStartToDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, B, C);

    double target_distance = 101;
    double extra_distance = target_distance - distance;
    A.ExtendFromStartToDistance(target_distance);
    Vector3D end = B + reverse * extra_distance;
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = 1;
    extra_distance = target_distance - distance;
    end = B + reverse * extra_distance;
    A.ExtendFromStartToDistance(target_distance);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(distance , A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = -101;
    extra_distance = target_distance - distance;
    end = B + reverse * extra_distance;
    A.ExtendFromStartToDistance(target_distance);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ExtendFromStartToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        double target_distance = distance + extra_distance;
        double target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, target_column_depth);
        ASSERT_DOUBLE_EQ(distance + distance/3.0, target_distance);
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistanceFromEndInReverse(target_column_depth));
        P.ExtendFromStartToColumnDepth(target_column_depth);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ(target_column_depth, sum);
        Vector3D end = inner_p1 - direction * target_distance;
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistance());
        ASSERT_TRUE((end - P.GetFirstPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetFirstPoint().magnitude()));

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, target_column_depth);
        P.ExtendFromStartToColumnDepth(target_column_depth);
        ASSERT_DOUBLE_EQ(distance - distance/3.0, target_distance);
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistanceFromEndInReverse(target_column_depth));
        end = inner_p1 - direction * target_distance;
        ASSERT_EQ(inner_p0, P.GetFirstPoint());
        ASSERT_DOUBLE_EQ(distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, target_column_depth);
        P.ExtendFromStartToColumnDepth(target_column_depth);
        ASSERT_DOUBLE_EQ(distance - distance*1.5, target_distance);
        ASSERT_DOUBLE_EQ(target_distance, P.GetDistanceFromEndAlongPath(target_column_depth));
        end = inner_p1 - direction * target_distance;
        ASSERT_EQ(inner_p0, P.GetFirstPoint());
        ASSERT_DOUBLE_EQ(distance, P.GetDistance());
    }
}

TEST(PointManipulation, ShrinkFromEndToDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, B, C);

    double target_distance = 101;
    double extra_distance = distance - target_distance;
    A.ShrinkFromEndToDistance(target_distance);

    Vector3D end = C - direction * extra_distance;
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = 1;
    extra_distance = distance - target_distance;
    end = C - direction * extra_distance;
    A.ShrinkFromEndToDistance(target_distance);
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = -101;
    extra_distance = distance - target_distance;
    end = C - direction * extra_distance;
    A.ShrinkFromEndToDistance(target_distance);
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ShrinkFromEndToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        double target_distance = distance - extra_distance;
        double target_column_depth = target_distance * rho * 100;
        P.ShrinkFromEndToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance/3.0, extra_distance);
        Vector3D end = inner_p1 - direction * extra_distance;
        ASSERT_TRUE((end - P.GetLastPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetLastPoint().magnitude()));
        EXPECT_DOUBLE_EQ(distance - extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance - extra_distance;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromEndToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance/3.0, extra_distance);
        end = inner_p1 - direction * extra_distance;
        EXPECT_EQ(inner_p1, P.GetLastPoint());
        EXPECT_DOUBLE_EQ(distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance - extra_distance;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromEndToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance*1.5, extra_distance);
        end = inner_p1 - direction * extra_distance;
        EXPECT_EQ(inner_p0, P.GetLastPoint());
        EXPECT_DOUBLE_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ShrinkFromStartToDistance) {
    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, B, C);

    double target_distance = 101;
    double extra_distance = distance - target_distance;
    A.ShrinkFromStartToDistance(target_distance);
    Vector3D end = B - reverse * extra_distance;
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = 1;
    extra_distance = distance - target_distance;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartToDistance(target_distance);
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, B, C);
    target_distance = -101;
    extra_distance = distance - target_distance;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartToDistance(target_distance);
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyEarthModelTest, ShrinkFromStartToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
        P = Path(A, p0, p1);
        sum = P.GetColumnDepthInBounds();
        ASSERT_DOUBLE_EQ((p1 - p0).magnitude() * rho * 100, sum);

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double target_distance = distance + extra_distance;
        double extra_column_depth = extra_distance * rho * 100;
        double target_column_depth = target_distance * rho * 100;
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, -direction, extra_column_depth);
        ASSERT_DOUBLE_EQ(distance/3.0, extra_distance);
        ASSERT_DOUBLE_EQ(extra_distance, P.GetDistanceFromStartInReverse(extra_column_depth));
        P.ShrinkFromStartToColumnDepth(target_column_depth);
        Vector3D end = inner_p0 + direction * extra_distance;
        ASSERT_EQ(inner_p0, P.GetFirstPoint());
        ASSERT_EQ(distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        target_distance = distance + extra_distance;
        extra_column_depth = extra_distance * rho * 100;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromStartToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, -direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(-distance/3.0, extra_distance);
        end = inner_p0 - direction * extra_distance;
        ASSERT_TRUE((end - P.GetFirstPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetFirstPoint().magnitude()));
        EXPECT_DOUBLE_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, inner_p0, inner_p1);
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        target_distance = distance + extra_distance;
        extra_column_depth = extra_distance * rho * 100;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromStartToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p0, direction, extra_column_depth);
        EXPECT_DOUBLE_EQ(distance*1.5, extra_distance);
        end = inner_p0 + direction * extra_distance;
        ASSERT_TRUE((inner_p0 - P.GetFirstPoint()).magnitude() < 1e-6 * std::max(inner_p0.magnitude(), P.GetFirstPoint().magnitude()));
        EXPECT_EQ(distance, P.GetDistance());
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude() * rho * 100, sum);
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthFromStartInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromStartInBounds(target_distance);
        EXPECT_NEAR(target_column_depth, sum, 1e-6 * std::max(target_column_depth, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartInBounds(target_distance);
        EXPECT_DOUBLE_EQ(0, sum);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartInBounds(target_distance);
        EXPECT_DOUBLE_EQ(distance * rho * 100, sum);
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthFromEndInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromEndInBounds(target_distance);
        EXPECT_NEAR(target_column_depth, sum, 1e-6 * std::max(target_column_depth, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndInBounds(target_distance);
        EXPECT_DOUBLE_EQ(0, sum);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndInBounds(target_distance);
        EXPECT_DOUBLE_EQ(distance * rho * 100, sum);
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthFromStartAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromStartAlongPath(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartAlongPath(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartAlongPath(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthFromEndAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromEndAlongPath(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndAlongPath(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = -distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndAlongPath(target_distance);
        EXPECT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthFromStartInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromStartInReverse(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartInReverse(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = -distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartInReverse(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));
    }
}

TEST_F(FakeLegacyEarthModelTest, GetColumnDepthFromEndInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromEndInReverse(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndInReverse(target_distance);
        ASSERT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndInReverse(target_distance);
        EXPECT_NEAR(target_column_depth, sum, 1e-6 * std::max(std::abs(target_column_depth), std::abs(sum)));
    }
}

/////////////////////////////////////////////////////////

TEST_F(FakeLegacyEarthModelTest, GetDistanceFromStartInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromStartInBounds(target_column_depth);
        EXPECT_NEAR(target_distance, sum, 1e-6 * std::max(target_distance, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartInBounds(target_column_depth);
        EXPECT_DOUBLE_EQ(0, sum);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartInBounds(target_column_depth);
        EXPECT_DOUBLE_EQ(distance, sum);
    }
}

TEST_F(FakeLegacyEarthModelTest, GetDistanceFromEndInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromEndInBounds(target_column_depth);
        EXPECT_NEAR(target_distance, sum, 1e-6 * std::max(target_distance, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndInBounds(target_column_depth);
        EXPECT_DOUBLE_EQ(0, sum);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndInBounds(target_column_depth);
        EXPECT_DOUBLE_EQ(distance, sum);
    }
}

TEST_F(FakeLegacyEarthModelTest, GetDistanceFromStartAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromStartAlongPath(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartAlongPath(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartAlongPath(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));
    }
}

TEST_F(FakeLegacyEarthModelTest, GetDistanceFromEndAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromEndAlongPath(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndAlongPath(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = -distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndAlongPath(target_column_depth);
        EXPECT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));
    }
}

TEST_F(FakeLegacyEarthModelTest, GetDistanceFromStartInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromStartInReverse(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartInReverse(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = -distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartInReverse(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));
    }
}

TEST_F(FakeLegacyEarthModelTest, GetDistanceFromEndInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<EarthModel> A(new EarthModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyEarthModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<EarthSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        EarthSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, inner_p0, inner_p1);
        DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<RadialAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromEndInReverse(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndInReverse(target_column_depth);
        ASSERT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndInReverse(target_column_depth);
        EXPECT_NEAR(target_distance, sum, 1e-6 * std::max(std::abs(target_distance), std::abs(sum)));
    }
}


// TEST()

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

