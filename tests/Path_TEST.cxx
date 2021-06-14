
#include <math.h>
#include <cmath>
#include <iostream>

#include <gtest/gtest.h>

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"

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
    EXPECT_FALSE(A.HasNucleonColumnDepth());
    EXPECT_FALSE(A.HasElectronColumnDepth());
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
    EXPECT_THROW(A.EnsureEarthModel(), const char*);
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
    EXPECT_THROW(A.EnsurePoints(), const char*);
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
    EXPECT_THROW(A.EnsureIntersections(), const char*);

    std::shared_ptr<const EarthModel> EMp(new EarthModel());
    A = Path(EMp);
    EXPECT_THROW(A.EnsureIntersections(), const char*);

    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    A = Path();
    A.SetPoints(B, C);
    EXPECT_THROW(A.EnsureIntersections(), const char*);

    A = Path();
    A.SetPointsWithRay(B, direction, distance);
    EXPECT_THROW(A.EnsureIntersections(), const char*);
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
        EXPECT_DOUBLE_EQ((inner_p1 - inner_p0).magnitude()*rho, sum);
    }
}

// TEST()

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

