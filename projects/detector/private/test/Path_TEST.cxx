
#include <math.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include <gtest/gtest.h>

#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/Sphere.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/DetectorModel.h"
#include "LeptonInjector/detector/Coordinates.h"
#include "LeptonInjector/detector/DensityDistribution.h"
#include "LeptonInjector/detector/DensityDistribution1D.h"
#include "LeptonInjector/detector/Distribution1D.h"
#include "LeptonInjector/detector/Axis1D.h"
#include "LeptonInjector/detector/RadialAxis1D.h"
#include "LeptonInjector/detector/CartesianAxis1D.h"
#include "LeptonInjector/detector/ConstantDistribution1D.h"
#include "LeptonInjector/detector/PolynomialDistribution1D.h"
#include "LeptonInjector/detector/ExponentialDistribution1D.h"

#include "LeptonInjector/dataclasses/Particle.h"

#include "FakeMaterialModel.h"
#include "FakeDetectorModel.h"

using namespace LI::detector;
using namespace LI::geometry;

#define EXPECT_NEAR_REL(A, B, C) \
    EXPECT_NEAR(A, B, C * std::max(std::abs(A), std::abs(B)));

#define EXPECT_VECTOR3D_NEAR_REL(A, B, C) \
    EXPECT_NEAR(A.GetX(), B.GetX(), C * std::max(std::abs(A.GetX()), std::abs(B.GetX()))); \
    EXPECT_NEAR(A.GetY(), B.GetY(), C * std::max(std::abs(A.GetY()), std::abs(B.GetY()))); \
    EXPECT_NEAR(A.GetZ(), B.GetZ(), C * std::max(std::abs(A.GetZ()), std::abs(B.GetZ())));

TEST(DefaultConstructor, NoThrow)
{
    EXPECT_NO_THROW(Path());
}

TEST(DefaultConstructor, HasNone)
{
    Path A;
    EXPECT_FALSE(A.HasDetectorModel());
    EXPECT_FALSE(A.HasPoints());
    EXPECT_FALSE(A.HasIntersections());
}

TEST(DefaultConstructor, MembersAreDefault)
{
    Path A;
    EXPECT_EQ(std::shared_ptr<const DetectorModel>(), A.GetDetectorModel());
    EXPECT_EQ(Vector3D(), A.GetFirstPoint());
    EXPECT_EQ(Vector3D(), A.GetLastPoint());
    EXPECT_EQ(Vector3D(), A.GetDirection());
    EXPECT_EQ(0, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(DetectorModelConstructor, NoThrow)
{
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    EXPECT_NO_THROW(Path A(EMp));
}

TEST(DetectorModelConstructor, MemberValues)
{
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Path A(EMp);
    EXPECT_EQ(EMp, A.GetDetectorModel());
    EXPECT_EQ(Vector3D(), A.GetFirstPoint());
    EXPECT_EQ(Vector3D(), A.GetLastPoint());
    EXPECT_EQ(Vector3D(), A.GetDirection());
    EXPECT_EQ(0, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(PointsConstructor, NoThrow)
{
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    EXPECT_NO_THROW(Path A(EMp, DetectorPosition(B), DetectorPosition(C)));
}

TEST(PointsConstructor, MemberValues)
{
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    Path A(EMp, DetectorPosition(B), DetectorPosition(C));
    EXPECT_EQ(EMp, A.GetDetectorModel());
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(RayConstructor, NoThrow)
{
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1, 2, 3);
    Vector3D C(4, 6, 8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    EXPECT_NO_THROW(Path A(EMp, DetectorPosition(B), DetectorDirection(direction), distance));
}

TEST(RayConstructor, MemberValues)
{
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
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
    Path A(EMp, DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_EQ(EMp, A.GetDetectorModel());
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(B + direction * distance, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
    EXPECT_EQ(Geometry::IntersectionList(), A.GetIntersections());
}

TEST(DetectorModel, SetGet) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Path A;
    A.SetDetectorModel(EMp);
    EXPECT_EQ(EMp, A.GetDetectorModel());
}

TEST(EnsureDetectorModel, Throw) {
    Path A;
    EXPECT_THROW(A.EnsureDetectorModel(), std::runtime_error);
}

TEST(EnsureDetectorModel, NoThrow) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Path A(EMp);
    EXPECT_NO_THROW(A.EnsureDetectorModel());
}

TEST(Points, SetGet) {
    Vector3D A(1,2,3);
    Vector3D B(4,6,8);
    Vector3D direction = B - A;
    double distance = direction.magnitude();
    direction.normalize();
    Path C;
    C.SetPoints(DetectorPosition(A), DetectorPosition(B));
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
    A.SetPointsWithRay(DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(B + direction * distance, A.GetLastPoint());
    EXPECT_EQ(direction, A.GetDirection());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST(EnsurePoints, ThrowDefault) {
    Path A;
    EXPECT_THROW(A.EnsurePoints(), std::runtime_error);
}

TEST(EnsurePoints, ThrowDetectorPositions) {
    Vector3D A(1,2,3);
    Vector3D B(4,6,8);
    Vector3D direction = B - A;
    direction.normalize();
    Path C;
    C.SetPoints(DetectorPosition(A), DetectorPosition(B));
    EXPECT_THROW(C.EnsurePoints(), std::runtime_error);
}

TEST(EnsurePoints, NoThrowGeometryPositions) {
    Vector3D A(1,2,3);
    Vector3D B(4,6,8);
    Vector3D direction = B - A;
    direction.normalize();
    Path C;
    C.SetPoints(GeometryPosition(A), GeometryPosition(B));
    EXPECT_NO_THROW(C.EnsurePoints());
}

TEST(EnsurePoints, ThrowRayDetectorPositions) {
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
    A.SetPointsWithRay(DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_THROW(A.EnsurePoints(), std::runtime_error);
}

TEST(EnsurePoints, NoThrowRayGeometryPositions) {
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
    A.SetPointsWithRay(GeometryPosition(B), GeometryDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsurePoints());
}

TEST(EnsureIntersections, Throw) {
    Path A;
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    A = Path(EMp);
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();
    A = Path();
    A.SetPoints(DetectorPosition(B), DetectorPosition(C));
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    A = Path();
    A.SetPointsWithRay(DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    A = Path();
    A.SetPoints(GeometryPosition(B), GeometryPosition(C));
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);

    A = Path();
    A.SetPointsWithRay(GeometryPosition(B), GeometryDirection(direction), distance);
    EXPECT_THROW(A.EnsureIntersections(), std::runtime_error);
}

TEST(EnsureIntersections, NoThrow) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A;

    A = Path();
    A.SetDetectorModel(EMp);
    A.SetPoints(DetectorPosition(B), DetectorPosition(C));
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path();
    A.SetDetectorModel(EMp);
    A.SetPointsWithRay(DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path();
    A.SetDetectorModel(EMp);
    A.SetPoints(GeometryPosition(B), GeometryPosition(C));
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path();
    A.SetDetectorModel(EMp);
    A.SetPointsWithRay(GeometryPosition(B), GeometryDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp);
    A.SetPoints(DetectorPosition(B), DetectorPosition(C));
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp);
    A.SetPointsWithRay(DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp);
    A.SetPoints(GeometryPosition(B), GeometryPosition(C));
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp);
    A.SetPointsWithRay(GeometryPosition(B), GeometryDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp, DetectorPosition(B), DetectorDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp, GeometryPosition(B), GeometryPosition(C));
    EXPECT_NO_THROW(A.EnsureIntersections());

    A = Path(EMp, GeometryPosition(B), GeometryDirection(direction), distance);
    EXPECT_NO_THROW(A.EnsureIntersections());
}

TEST(PointManipulation, Flip) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    Vector3D reverse = -direction;
    double distance = direction.magnitude();
    direction.normalize();
    reverse.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));
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
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double extra_distance = 101;
    A.ExtendFromEndByDistance(extra_distance);

    Vector3D end = C + direction * extra_distance;
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = -1;
    end = C + direction * extra_distance;
    A.ExtendFromEndByDistance(extra_distance);
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = -101;
    end = C + direction * extra_distance;
    A.ExtendFromEndByDistance(extra_distance);
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ExtendFromEndByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D direction = p1 - p0;
        DetectorPosition p0_det = A->ToDet(GeometryPosition(p0));
        DetectorPosition p1_det = A->ToDet(GeometryPosition(p1));
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        DetectorPosition inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        DetectorPosition inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        DetectorDirection direction_det = A->ToDet(GeometryDirection(direction));
        Path P(A, GeometryPosition(inner_p0), GeometryPosition(inner_p1));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, GeometryPosition(p0), GeometryPosition(p1));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, GeometryPosition(inner_p0), GeometryPosition(inner_p1));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1_det, direction_det, extra_column_depth);
        EXPECT_NEAR_REL(distance/3.0, extra_distance, 1e-8);
        Vector3D end = inner_p1_det + direction_det.get() * extra_distance;
        EXPECT_NEAR((end - P.GetLastPoint()).magnitude(), 0, 1e-8);
        EXPECT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, GeometryPosition(inner_p0), GeometryPosition(inner_p1));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1_det, direction_det, extra_column_depth);
        EXPECT_NEAR_REL(-distance/3.0, extra_distance, 1e-8);
        end = inner_p1_det + direction_det.get() * extra_distance;
        EXPECT_NEAR((end - P.GetLastPoint()).magnitude(), 0, 1e-8);
        EXPECT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, GeometryPosition(inner_p0), GeometryPosition(inner_p1));
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(inner_p1_det, direction_det, extra_column_depth);
        EXPECT_NEAR_REL(-distance*1.5, extra_distance, 1e-8);
        end = inner_p1_det + direction_det.get() * extra_distance;
        EXPECT_NEAR((inner_p0_det - P.GetLastPoint())->magnitude(), 0, 1e-8);
        EXPECT_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ExtendFromStartByDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double extra_distance = 101;
    A.ExtendFromStartByDistance(extra_distance);
    Vector3D end = B + reverse * extra_distance;
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = -1;
    end = B + reverse * extra_distance;
    A.ExtendFromStartByDistance(extra_distance);
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = -101;
    end = B + reverse * extra_distance;
    A.ExtendFromStartByDistance(extra_distance);
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ExtendFromStartByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1_det - p0_det;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1_det - inner_p0_det;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(distance/3.0, extra_distance, 1e-8);
        EXPECT_NEAR_REL(A->DistanceForColumnDepthFromPoint(P.GetFirstPoint(), DetectorDirection(-direction), extra_column_depth), A->DistanceForColumnDepthFromPoint(P.GetIntersections(), P.GetFirstPoint(), DetectorDirection(-direction), extra_column_depth), 1e-8);
        EXPECT_NEAR_REL(extra_distance, P.GetDistanceFromStartInReverse(extra_column_depth), 1e-8);
        P.ExtendFromStartByColumnDepth(extra_column_depth);
        Vector3D end = inner_p0_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(end, Vector3D(P.GetFirstPoint()), 1e-8);
        ASSERT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(-distance/3.0, extra_distance, 1e-8);
        end = inner_p0_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(end, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_EQ(distance + extra_distance, P.GetDistance());

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ExtendFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(-distance*1.5, extra_distance, 1e-8);
        end = inner_p0_det - direction * extra_distance;
        EXPECT_EQ(inner_p1_det, P.GetFirstPoint());
        EXPECT_EQ(0, P.GetDistance());
    }
}

TEST(PointManipulation, ShrinkFromEndByDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double extra_distance = -101;
    A.ShrinkFromEndByDistance(extra_distance);

    Vector3D end = C - direction * extra_distance;
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = 1;
    end = C - direction * extra_distance;
    A.ShrinkFromEndByDistance(extra_distance);
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = 101;
    end = C - direction * extra_distance;
    A.ShrinkFromEndByDistance(extra_distance);
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ShrinkFromEndByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(distance/3.0, extra_distance, 1e-8);
        Vector3D end = inner_p1_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(end, Vector3D(P.GetLastPoint()), 1e-8);
        EXPECT_NEAR_REL(distance - extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(-distance/3.0, extra_distance, 1e-8);
        end = inner_p1_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(end, Vector3D(P.GetLastPoint()), 1e-8);
        EXPECT_NEAR_REL(distance - extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromEndByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(distance*1.5, extra_distance, 1e-8);
        end = inner_p1_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p0_det, Vector3D(P.GetLastPoint()), 1e-8);
        EXPECT_NEAR_REL(0.0, P.GetDistance(), 1e-8);
    }
}

TEST(PointManipulation, ShrinkFromStartByDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double extra_distance = -101;
    A.ShrinkFromStartByDistance(extra_distance);
    Vector3D end = B - reverse * extra_distance;
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = 1;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartByDistance(extra_distance);
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    extra_distance = 101;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartByDistance(extra_distance);
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ShrinkFromStartByColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1_det - p0_det;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), extra_column_depth);
        EXPECT_NEAR_REL(distance/3.0, extra_distance, 1e-8);
        EXPECT_NEAR_REL(extra_distance, P.GetDistanceFromStartInReverse(extra_column_depth), 1e-8);
        P.ShrinkFromStartByColumnDepth(extra_column_depth);
        Vector3D end = inner_p0_det + direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(end, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_NEAR_REL(distance - extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), extra_column_depth);
        EXPECT_NEAR_REL(-distance/3.0, extra_distance, 1e-8);
        end = inner_p0_det + direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(end, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_NEAR_REL(distance - extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        P.ShrinkFromStartByColumnDepth(extra_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), extra_column_depth);
        EXPECT_NEAR_REL(distance*1.5, extra_distance, 1e-8);
        end = inner_p0_det + direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p1_det, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_NEAR_REL(0.0, P.GetDistance(), 1e-8);
    }
}

TEST(PointManipulation, ExtendFromEndToDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double target_distance = 101;
    double extra_distance = target_distance - distance;
    A.ExtendFromEndToDistance(target_distance);

    Vector3D end = C + direction * extra_distance;
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = 1;
    extra_distance = target_distance - distance;
    end = C + direction * extra_distance;
    A.ExtendFromEndToDistance(target_distance);
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = -101;
    extra_distance = target_distance - distance;
    end = C + direction * extra_distance;
    A.ExtendFromEndToDistance(target_distance);
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ExtendFromEndToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double target_distance = distance + extra_distance;
        double target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), target_column_depth);
        EXPECT_NEAR_REL(distance + distance/3.0, target_distance, 1e-8);
        EXPECT_NEAR_REL(target_distance, P.GetDistanceFromStartAlongPath(target_column_depth), 1e-8);
        P.ExtendFromEndToColumnDepth(target_column_depth);
        Vector3D end = inner_p0_det + direction * target_distance;
        ASSERT_TRUE((end - P.GetLastPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetLastPoint()->magnitude()));
        EXPECT_NEAR_REL(distance + extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), target_column_depth);
        P.ExtendFromEndToColumnDepth(target_column_depth);
        EXPECT_NEAR_REL(distance - distance/3.0, target_distance, 1e-8);
        EXPECT_NEAR_REL(target_distance, P.GetDistanceFromStartAlongPath(target_column_depth), 1e-8);
        end = inner_p0_det + direction * target_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p1_det, Vector3D(P.GetLastPoint()), 1e-8);
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), target_column_depth);
        P.ExtendFromEndToColumnDepth(target_column_depth);
        EXPECT_NEAR_REL(distance - distance*1.5, target_distance, 1e-8);
        EXPECT_NEAR_REL(target_distance, P.GetDistanceFromStartAlongPath(target_column_depth), 1e-8);
        end = inner_p0_det + direction * target_distance;
        ASSERT_EQ(inner_p1_det, P.GetLastPoint());
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);
    }
}


TEST(PointManipulation, ExtendFromStartToDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double target_distance = 101;
    double extra_distance = target_distance - distance;
    A.ExtendFromStartToDistance(target_distance);
    Vector3D end = B + reverse * extra_distance;
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance + extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = 1;
    extra_distance = target_distance - distance;
    end = B + reverse * extra_distance;
    A.ExtendFromStartToDistance(target_distance);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(distance , A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = -101;
    extra_distance = target_distance - distance;
    end = B + reverse * extra_distance;
    A.ExtendFromStartToDistance(target_distance);
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(distance, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ExtendFromStartToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double target_distance = distance + extra_distance;
        double target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), target_column_depth);
        EXPECT_NEAR_REL(distance + distance/3.0, target_distance, 1e-8);
        EXPECT_NEAR_REL(target_distance, P.GetDistanceFromEndInReverse(target_column_depth), 1e-8);
        P.ExtendFromStartToColumnDepth(target_column_depth);
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL(target_column_depth, sum, 1e-8);
        Vector3D end = inner_p1_det - direction * target_distance;
        EXPECT_NEAR_REL(target_distance, P.GetDistance(), 1e-8);
        ASSERT_TRUE((end - P.GetFirstPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetFirstPoint()->magnitude()));

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), target_column_depth);
        P.ExtendFromStartToColumnDepth(target_column_depth);
        EXPECT_NEAR_REL(distance - distance/3.0, target_distance, 1e-8);
        EXPECT_NEAR_REL(target_distance, P.GetDistanceFromEndInReverse(target_column_depth), 1e-8);
        end = inner_p1_det - direction * target_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p0_det, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance*1.5;
        target_distance = distance + extra_distance;
        target_column_depth = target_distance * rho * 100;
        target_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), target_column_depth);
        P.ExtendFromStartToColumnDepth(target_column_depth);
        EXPECT_NEAR_REL(distance - distance*1.5, target_distance, 1e-8);
        EXPECT_NEAR_REL(target_distance, P.GetDistanceFromEndAlongPath(target_column_depth), 1e-8);
        end = inner_p1_det - direction * target_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p0_det, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);
    }
}

TEST(PointManipulation, ShrinkFromEndToDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    direction.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double target_distance = 101;
    double extra_distance = distance - target_distance;
    A.ShrinkFromEndToDistance(target_distance);

    Vector3D end = C - direction * extra_distance;
    EXPECT_EQ(C, A.GetLastPoint());
    EXPECT_EQ(distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = 1;
    extra_distance = distance - target_distance;
    end = C - direction * extra_distance;
    A.ShrinkFromEndToDistance(target_distance);
    EXPECT_EQ(end, A.GetLastPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = -101;
    extra_distance = distance - target_distance;
    end = C - direction * extra_distance;
    A.ShrinkFromEndToDistance(target_distance);
    EXPECT_EQ(B, A.GetLastPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ShrinkFromEndToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double extra_column_depth = extra_distance * rho * 100;
        double target_distance = distance - extra_distance;
        double target_column_depth = target_distance * rho * 100;
        P.ShrinkFromEndToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(distance/3.0, extra_distance, 1e-8);
        Vector3D end = inner_p1_det - direction * extra_distance;
        ASSERT_TRUE((end - P.GetLastPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetLastPoint()->magnitude()));
        EXPECT_NEAR_REL(distance - extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance - extra_distance;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromEndToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(-distance/3.0, extra_distance, 1e-8);
        end = inner_p1_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p1_det, Vector3D(P.GetLastPoint()), 1e-8);
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        extra_column_depth = extra_distance * rho * 100;
        target_distance = distance - extra_distance;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromEndToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p1_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(distance*1.5, extra_distance, 1e-8);
        end = inner_p1_det - direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p0_det, Vector3D(P.GetLastPoint()), 1e-8);
        EXPECT_NEAR_REL(0.0, P.GetDistance(), 1e-8);
    }
}

TEST(PointManipulation, ShrinkFromStartToDistance) {
    std::shared_ptr<const DetectorModel> EMp(new DetectorModel());
    Vector3D B(1,2,3);
    Vector3D C(4,6,8);
    Vector3D direction = C - B;
    double distance = direction.magnitude();
    Vector3D reverse = -direction;
    direction.normalize();
    reverse.normalize();

    Path A(EMp, DetectorPosition(B), DetectorPosition(C));

    double target_distance = 101;
    double extra_distance = distance - target_distance;
    A.ShrinkFromStartToDistance(target_distance);
    Vector3D end = B - reverse * extra_distance;
    EXPECT_EQ(B, A.GetFirstPoint());
    EXPECT_EQ(distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = 1;
    extra_distance = distance - target_distance;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartToDistance(target_distance);
    EXPECT_EQ(end, A.GetFirstPoint());
    EXPECT_EQ(distance - extra_distance, A.GetDistance());

    A = Path(EMp, DetectorPosition(B), DetectorPosition(C));
    target_distance = -101;
    extra_distance = distance - target_distance;
    end = B - reverse * extra_distance;
    A.ShrinkFromStartToDistance(target_distance);
    EXPECT_EQ(C, A.GetFirstPoint());
    EXPECT_EQ(0, A.GetDistance());
}

TEST_F(FakeLegacyDetectorModelTest, ShrinkFromStartToColumnDepth) {
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        ASSERT_TRUE(p0.magnitude() < max_radius);
        ASSERT_TRUE(p1.magnitude() < max_radius);
        ASSERT_TRUE(inner_p0.magnitude() < max_radius);
        ASSERT_TRUE(inner_p1.magnitude() < max_radius);
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
        P = Path(A, DetectorPosition(p0_det), DetectorPosition(p1_det));
        sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((p1 - p0).magnitude() * rho * 100, sum, 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        double extra_distance = distance/3.0;
        double target_distance = distance + extra_distance;
        double extra_column_depth = extra_distance * rho * 100;
        double target_column_depth = target_distance * rho * 100;
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(distance/3.0, extra_distance, 1e-8);
        EXPECT_NEAR_REL(extra_distance, P.GetDistanceFromStartInReverse(extra_column_depth), 1e-8);
        P.ShrinkFromStartToColumnDepth(target_column_depth);
        Vector3D end = inner_p0_det + direction * extra_distance;
        EXPECT_VECTOR3D_NEAR_REL(inner_p0_det, Vector3D(P.GetFirstPoint()), 1e-8);
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = -distance/3.0;
        target_distance = distance + extra_distance;
        extra_column_depth = extra_distance * rho * 100;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromStartToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(-direction), extra_column_depth);
        EXPECT_NEAR_REL(-distance/3.0, extra_distance, 1e-8);
        end = inner_p0_det - direction * extra_distance;
        ASSERT_TRUE((end - P.GetFirstPoint()).magnitude() < 1e-6 * std::max(end.magnitude(), P.GetFirstPoint()->magnitude()));
        EXPECT_NEAR_REL(distance + extra_distance, P.GetDistance(), 1e-8);

        P = Path(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        P.EnsureIntersections();
        extra_distance = distance*1.5;
        target_distance = distance + extra_distance;
        extra_column_depth = extra_distance * rho * 100;
        target_column_depth = target_distance * rho * 100;
        P.ShrinkFromStartToColumnDepth(target_column_depth);
        extra_distance = A->DistanceForColumnDepthFromPoint(DetectorPosition(inner_p0_det), DetectorDirection(direction), extra_column_depth);
        EXPECT_NEAR_REL(distance*1.5, extra_distance, 1e-8);
        end = inner_p0_det + direction * extra_distance;
        ASSERT_TRUE((inner_p0_det - P.GetFirstPoint()).magnitude() < 1e-6 * std::max(inner_p0.magnitude(), P.GetFirstPoint()->magnitude()));
        EXPECT_NEAR_REL(distance, P.GetDistance(), 1e-8);
    }
}

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());
        double sum = P.GetColumnDepthInBounds();
        EXPECT_NEAR_REL((inner_p1 - inner_p0).magnitude() * rho * 100, sum, 1e-8);
    }
}

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthFromStartInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromStartInBounds(target_distance);
        EXPECT_NEAR(target_column_depth, sum, 1e-6 * std::max(target_column_depth, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartInBounds(target_distance);
        EXPECT_NEAR_REL(0.0, sum, 1e-8);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromStartInBounds(target_distance);
        EXPECT_NEAR_REL(distance * rho * 100, sum, 1e-8);
    }
}

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthFromEndInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetColumnDepthFromEndInBounds(target_distance);
        EXPECT_NEAR(target_column_depth, sum, 1e-6 * std::max(target_column_depth, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndInBounds(target_distance);
        EXPECT_NEAR_REL(0.0, sum, 1e-8);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetColumnDepthFromEndInBounds(target_distance);
        EXPECT_NEAR_REL(distance * rho * 100, sum, 1e-8);
    }
}

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthFromStartAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthFromEndAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthFromStartInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetColumnDepthFromEndInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetDistanceFromStartInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromStartInBounds(target_column_depth);
        EXPECT_NEAR(target_distance, sum, 1e-6 * std::max(target_distance, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartInBounds(target_column_depth);
        EXPECT_NEAR_REL(0.0, sum, 1e-8);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromStartInBounds(target_column_depth);
        EXPECT_NEAR_REL(distance, sum, 1e-8);
    }
}

TEST_F(FakeLegacyDetectorModelTest, GetDistanceFromEndInBounds)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
        ASSERT_TRUE(density);
        double rho = density->Evaluate(Vector3D());

        double target_distance = distance / 3.0;
        double target_column_depth = target_distance * rho * 100;
        double sum = P.GetDistanceFromEndInBounds(target_column_depth);
        EXPECT_NEAR(target_distance, sum, 1e-6 * std::max(target_distance, sum));

        target_distance = -distance / 3.0;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndInBounds(target_column_depth);
        EXPECT_NEAR_REL(0.0, sum, 1e-8);

        target_distance = distance * 1.5;
        target_column_depth = target_distance * rho * 100;
        sum = P.GetDistanceFromEndInBounds(target_column_depth);
        EXPECT_NEAR_REL(distance, sum, 1e-8);
    }
}

TEST_F(FakeLegacyDetectorModelTest, GetDistanceFromStartAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetDistanceFromEndAlongPath)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetDistanceFromStartInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

TEST_F(FakeLegacyDetectorModelTest, GetDistanceFromEndInReverse)
{
    unsigned int N_rand = 1000;
    for(unsigned int i=0; i<N_rand; ++i) {
        ASSERT_NO_THROW(reset(1, 1));
        std::shared_ptr<DetectorModel> A(new DetectorModel());
        ASSERT_NO_THROW(A->LoadMaterialModel(materials_file));
        double max_depth = 5000;
        max_depth = std::min(max_depth, *std::max_element(layer_radii.begin(), layer_radii.end()));
        double depth = FakeLegacyDetectorModelFile::RandomDouble()*max_depth;
        double ice_angle = -1;
        ASSERT_NO_THROW(A->LoadConcentricShellsFromLegacyFile(model_file, depth, ice_angle));
        std::vector<DetectorSector> sectors = A->GetSectors();
        ASSERT_EQ(2, sectors.size());
        DetectorSector sector = sectors[1];
        Sphere const * sphere = dynamic_cast<Sphere const *>(sector.geo.get());
        ASSERT_TRUE(sphere);
        double max_radius = sphere->GetRadius();
        double min_radius = sphere->GetInnerRadius();
        Vector3D p0 = RandomVector(max_radius, min_radius);
        Vector3D p1 = RandomVector(max_radius, min_radius);
        Vector3D p0_det = A->ToDet(GeometryPosition(p0));
        Vector3D p1_det = A->ToDet(GeometryPosition(p1));
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        Vector3D inner_p0 = p0 + direction * distance / 4.0;
        Vector3D inner_p1 = p1 - direction * distance / 4.0;
        Vector3D inner_p0_det = A->ToDet(GeometryPosition(inner_p0));
        Vector3D inner_p1_det = A->ToDet(GeometryPosition(inner_p1));
        direction = inner_p1 - inner_p0;
        distance = direction.magnitude();
        direction.normalize();
        Path P(A, DetectorPosition(inner_p0_det), DetectorPosition(inner_p1_det));
        DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const * density = dynamic_cast<DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> const *>(sector.density.get());
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

