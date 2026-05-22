// Tests for geometry operator semantics (==, !=, <, copy, swap) and
// validation of degenerate construction parameters. Also includes
// explicit multi-intersection position verification for hollow shapes
// and DistanceToClosestApproach.

#include <cmath>
#include <memory>
#include <vector>
#include <utility>

#include <gtest/gtest.h>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/ExtrPoly.h"
#include "SIREN/geometry/CutTube.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/Placement.h"

using namespace siren::geometry;
using namespace siren::math;

// =========================================================================
// Equality / inequality operators
// =========================================================================

TEST(GeometryOps, EqualitySameType) {
    Sphere a;
    Sphere b;
    EXPECT_TRUE(a == b);

    Geometry* c = new Sphere();
    Geometry* d = new Sphere();
    EXPECT_TRUE(*c == *d);
    delete c;
    delete d;
}

TEST(GeometryOps, EqualityDifferentType) {
    Sphere a;
    Box b;
    EXPECT_TRUE(a != b);

    Geometry* c = new Sphere();
    Geometry* d = new Box();
    EXPECT_TRUE(*c != *d);
    delete c;
    delete d;
}

TEST(GeometryOps, EqualityExtrPoly) {
    ExtrPoly a;
    ExtrPoly b;
    EXPECT_TRUE(a == b);

    ExtrPoly c;
    Box d;
    EXPECT_TRUE(c != d);
}

// =========================================================================
// Copy / assignment / swap
// =========================================================================

TEST(GeometryOps, CopyConstructor) {
    Sphere a(Vector3D(1, 0, 0), 20, 10);
    Sphere b(a);
    EXPECT_TRUE(a == b);
}

TEST(GeometryOps, AssignmentOperator) {
    Sphere a;
    Sphere b(Vector3D(0, 0, 0), 2, 1);
    EXPECT_TRUE(a != b);
    b = a;
    EXPECT_TRUE(a == b);
}

TEST(GeometryOps, Swap) {
    Sphere a;
    Sphere b;
    Sphere c(Vector3D(1, 2, 3), 4, 3);
    Sphere d(Vector3D(1, 2, 3), 4, 3);
    EXPECT_TRUE(c == d);
    a.swap(c);
    EXPECT_TRUE(a == d);
    EXPECT_TRUE(b == c);
}

TEST(GeometryOps, CrossTypeAssignNoEffect) {
    Geometry* c = new Sphere();
    Geometry* d = new Box();
    *d = *c;
    EXPECT_FALSE(*c == *d);
    delete c;
    delete d;
}

TEST(GeometryOps, SameTypeAssignThroughBase) {
    Geometry* c = new Sphere();
    Geometry* e = new Sphere(Vector3D(1, 0, 0), 20, 10);
    *c = *e;
    EXPECT_TRUE(*c == *e);
    delete c;
    delete e;
}

TEST(GeometryOps, ExtrPolyCopyAndAssign) {
    std::vector<std::vector<double>> poly = {{0,0},{1,0},{1,1},{0,1}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {ExtrPoly::ZSection(0, off, 1)};
    ExtrPoly a(poly, zsecs);
    ExtrPoly b(a);
    EXPECT_TRUE(a == b);

    ExtrPoly c;
    EXPECT_TRUE(a != c);
    c = a;
    EXPECT_TRUE(a == c);
}

TEST(GeometryOps, ExtrPolySwap) {
    std::vector<std::vector<double>> poly = {{0,0},{1,0},{1,1},{0,1}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {ExtrPoly::ZSection(0, off, 1)};
    ExtrPoly a(poly, zsecs);
    ExtrPoly b(poly, zsecs);
    ExtrPoly c;
    ExtrPoly d;
    EXPECT_TRUE(a == b);
    c.swap(a);
    EXPECT_TRUE(c == b);
    EXPECT_TRUE(d == a);
}

// =========================================================================
// DistanceToClosestApproach
// =========================================================================

TEST(GeometryOps, DistanceToClosestApproachSphere) {
    Sphere a;
    Vector3D position(1, -1, 0);
    Vector3D direction(0, 1, 0);
    direction.CalculateSphericalCoordinates();
    EXPECT_NEAR(a.DistanceToClosestApproach(position, direction), 1.0, 1e-9);
}

TEST(GeometryOps, DistanceToClosestApproachExtrPoly) {
    ExtrPoly a;
    Vector3D position(1, -1, 0);
    Vector3D direction(0, 1, 0);
    direction.CalculateSphericalCoordinates();
    EXPECT_NEAR(a.DistanceToClosestApproach(position, direction), 1.0, 1e-9);
}

// =========================================================================
// Validation: degenerate parameters should throw
// =========================================================================

TEST(GeometryOps, ConeZeroOuterRadiiThrows) {
    EXPECT_THROW(Cone(0, 0, 0, 0, 10), std::invalid_argument);
    EXPECT_NO_THROW(Cone(0, 5, 0, 0, 8));
    EXPECT_NO_THROW(Cone(0, 0, 0, 5, 8));
}

TEST(GeometryOps, TrdZeroHeightThrows) {
    EXPECT_THROW(Trd(5, 3, 4, 2, 0), std::invalid_argument);
    EXPECT_THROW(Trd(5, 3, 4, 2, -1), std::invalid_argument);
    EXPECT_NO_THROW(Trd(5, 3, 4, 2, 6));
}

TEST(GeometryOps, CutTubeNegativeRminThrows) {
    Vector3D low_norm(0, 0, -1);
    Vector3D high_norm(0, 0, 1);
    EXPECT_THROW(CutTube(-5.0, 10.0, 5.0, low_norm, high_norm), std::invalid_argument);
    EXPECT_THROW(CutTube(0.0, -3.0, 5.0, low_norm, high_norm), std::invalid_argument);
    EXPECT_THROW(CutTube(-2.0, -1.0, 5.0, low_norm, high_norm), std::invalid_argument);
    EXPECT_NO_THROW(CutTube(8.0, 5.0, 3.0, low_norm, high_norm));
    EXPECT_NO_THROW(CutTube(0.0, 5.0, 3.0, low_norm, high_norm));

    Placement p;
    EXPECT_THROW(CutTube(p, -5.0, 10.0, 5.0, low_norm, high_norm), std::invalid_argument);
    EXPECT_NO_THROW(CutTube(p, 8.0, 5.0, 3.0, low_norm, high_norm));
}

// =========================================================================
// Multi-intersection explicit position verification for hollow shapes
// =========================================================================

TEST(GeometryOps, HollowSphereTraversal) {
    Sphere hollow(10, 5);

    // Axial ray through center: expect 4 hits at z = 10, 5, -5, -10
    Vector3D origin(0, 0, 20);
    Vector3D dir(0, 0, -1);
    auto hits = hollow.Intersections(origin, dir);
    ASSERT_EQ(hits.size(), 4u);
    EXPECT_TRUE(hits[0].entering);
    EXPECT_NEAR(hits[0].position.GetZ(), 10.0, 1e-6);
    EXPECT_FALSE(hits[1].entering);
    EXPECT_NEAR(hits[1].position.GetZ(), 5.0, 1e-6);
    EXPECT_TRUE(hits[2].entering);
    EXPECT_NEAR(hits[2].position.GetZ(), -5.0, 1e-6);
    EXPECT_FALSE(hits[3].entering);
    EXPECT_NEAR(hits[3].position.GetZ(), -10.0, 1e-6);

    // Offset ray at y=8 misses inner shell: expect 2 hits
    Vector3D origin2(0, 8, 20);
    auto hits2 = hollow.Intersections(origin2, dir);
    ASSERT_EQ(hits2.size(), 2u);
    EXPECT_TRUE(hits2[0].entering);
    EXPECT_NEAR(hits2[0].position.GetZ(), 6.0, 1e-6);
    EXPECT_FALSE(hits2[1].entering);
    EXPECT_NEAR(hits2[1].position.GetZ(), -6.0, 1e-6);
}

TEST(GeometryOps, HollowCylinderTraversal) {
    Cylinder hollow(10, 5, 20);

    // Radial ray through center: expect 4 hits at x = 10, 5, -5, -10
    Vector3D origin(20, 0, 0);
    Vector3D dir(-1, 0, 0);
    auto hits = hollow.Intersections(origin, dir);
    ASSERT_EQ(hits.size(), 4u);
    EXPECT_TRUE(hits[0].entering);
    EXPECT_NEAR(hits[0].position.GetX(), 10.0, 1e-6);
    EXPECT_FALSE(hits[1].entering);
    EXPECT_NEAR(hits[1].position.GetX(), 5.0, 1e-6);
    EXPECT_TRUE(hits[2].entering);
    EXPECT_NEAR(hits[2].position.GetX(), -5.0, 1e-6);
    EXPECT_FALSE(hits[3].entering);
    EXPECT_NEAR(hits[3].position.GetX(), -10.0, 1e-6);

    // Axial ray through bore: expect 0 hits
    Vector3D bore_origin(0, 3, 20);
    Vector3D bore_dir(0, 0, -1);
    auto bore_hits = hollow.Intersections(bore_origin, bore_dir);
    EXPECT_EQ(bore_hits.size(), 0u);
}

TEST(GeometryOps, HollowConeTraversal) {
    // rmin1=1, rmax1=5, rmin2=1, rmax2=3, z=8 -> z from -4 to +4
    // At midplane: outer=4, inner=1
    Cone hollow(1, 5, 1, 3, 8);

    Vector3D origin(10, 0, 0);
    Vector3D dir(-1, 0, 0);
    auto hits = hollow.Intersections(origin, dir);
    ASSERT_EQ(hits.size(), 4u);
    EXPECT_TRUE(hits[0].entering);
    EXPECT_NEAR(hits[0].position.GetX(), 4.0, 1e-6);
    EXPECT_FALSE(hits[1].entering);
    EXPECT_NEAR(hits[1].position.GetX(), 1.0, 1e-6);
    EXPECT_TRUE(hits[2].entering);
    EXPECT_NEAR(hits[2].position.GetX(), -1.0, 1e-6);
    EXPECT_FALSE(hits[3].entering);
    EXPECT_NEAR(hits[3].position.GetX(), -4.0, 1e-6);
}

// =========================================================================
// Z-plane ordering: constructor acceptance and rejection
// =========================================================================

TEST(GeometryOps, PolyconeAcceptsAscendingZ) {
    EXPECT_NO_THROW(Polycone({-5, 0, 5}, {0, 0, 0}, {3, 5, 3}));
}

TEST(GeometryOps, PolyconeAcceptsDescendingZ) {
    EXPECT_NO_THROW(Polycone({5, 0, -5}, {0, 0, 0}, {3, 5, 3}));
}

TEST(GeometryOps, PolyconeAcceptsDuplicateZ) {
    // Step profile: narrow-wide-narrow
    EXPECT_NO_THROW(Polycone({-5, -2, -2, 2, 2, 5},
                             {0, 0, 0, 0, 0, 0},
                             {3, 3, 5, 5, 3, 3}));
}

TEST(GeometryOps, PolyconeRejectsNonMonotonicZ) {
    EXPECT_THROW(Polycone({-5, 5, 0}, {0, 0, 0}, {3, 3, 5}), std::runtime_error);
    EXPECT_THROW(Polycone({0, 5, -5, 10}, {0, 0, 0, 0}, {3, 3, 3, 3}), std::runtime_error);
}

TEST(GeometryOps, PolyhedraAcceptsAscendingZ) {
    EXPECT_NO_THROW(Polyhedra(6, 0, {-5, 5}, {0, 0}, {4, 4}));
}

TEST(GeometryOps, PolyhedraAcceptsDescendingZ) {
    EXPECT_NO_THROW(Polyhedra(6, 0, {5, -5}, {0, 0}, {4, 4}));
}

TEST(GeometryOps, PolyhedraAcceptsDuplicateZ) {
    EXPECT_NO_THROW(Polyhedra(6, 0, {-5, -1, -1, 1, 1, 5},
                              {0, 0, 0, 0, 0, 0},
                              {3, 3, 5, 5, 3, 3}));
}

TEST(GeometryOps, PolyhedraRejectsNonMonotonicZ) {
    EXPECT_THROW(Polyhedra(6, 0, {-5, 5, 0}, {0, 0, 0}, {3, 3, 5}), std::runtime_error);
}

// =========================================================================
// Entering flag consistency for hollow shapes
// =========================================================================

TEST(GeometryOps, EnteringFlagConsistency) {
    Sphere hollow_sphere(10, 5);
    Cylinder hollow_cyl(10, 5, 20);
    Cone hollow_cone(1, 5, 1, 3, 8);

    double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    double inv_sqrt3 = 1.0 / std::sqrt(3.0);
    struct RayDef { double ox, oy, oz, dx, dy, dz; };
    RayDef rays[10] = {
        { 30, 0, 0,  -1, 0, 0}, {-30, 0, 0,  1, 0, 0},
        { 0, 30, 0,   0,-1, 0}, { 0,-30, 0,  0, 1, 0},
        { 0, 0, 30,   0, 0,-1}, { 0, 0,-30,  0, 0, 1},
        { 30, 30, 0,  -inv_sqrt2, -inv_sqrt2, 0},
        { 0, 30, 30,  0, -inv_sqrt2, -inv_sqrt2},
        { 30, 0, 30,  -inv_sqrt2, 0, -inv_sqrt2},
        { 30, 30, 30, -inv_sqrt3, -inv_sqrt3, -inv_sqrt3},
    };

    Geometry const * shapes[3] = {&hollow_sphere, &hollow_cyl, &hollow_cone};
    const char * names[3] = {"Sphere", "Cylinder", "Cone"};

    for(int s = 0; s < 3; ++s) {
        for(int r = 0; r < 10; ++r) {
            Vector3D origin(rays[r].ox, rays[r].oy, rays[r].oz);
            Vector3D dir(rays[r].dx, rays[r].dy, rays[r].dz);
            auto hits = shapes[s]->Intersections(origin, dir);
            if(hits.empty()) continue;

            EXPECT_TRUE(hits[0].entering) << names[s] << " ray " << r;
            EXPECT_FALSE(hits.back().entering) << names[s] << " ray " << r;
            EXPECT_EQ(hits.size() % 2, 0u) << names[s] << " ray " << r;

            for(size_t i = 1; i < hits.size(); ++i) {
                EXPECT_GE(hits[i].distance, hits[i-1].distance) << names[s] << " ray " << r;
                EXPECT_NE(hits[i].entering, hits[i-1].entering) << names[s] << " ray " << r;
            }
        }
    }
}
