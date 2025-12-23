#include <limits>
#include <cmath>

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

#include "SIREN/math/EulerQuaternionConversions.h" // for QFromZXZr

using namespace siren::math;
using namespace siren::geometry;
using namespace siren::detector;
using namespace siren::utilities;
using namespace siren::dataclasses;

namespace {

siren::math::Vector3D RandomUnit(std::mt19937 & rng) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    double z = 2.0 * u01(rng) - 1.0;
    double phi = 2.0 * M_PI * u01(rng);
    double r_xy = std::sqrt(std::max(0.0, 1.0 - z*z));
    return siren::math::Vector3D(r_xy * std::cos(phi), r_xy * std::sin(phi), z);
}

siren::math::Vector3D RandomPointInBall(std::mt19937 & rng, double rmax) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    // radius ~ r^3 uniform
    double rr = std::cbrt(u01(rng)) * rmax;
    siren::math::Vector3D dir = RandomUnit(rng);
    return dir * rr;
}

// A deterministic toy model with distinct constant densities per sector.
// Levels: 0 = world, 1 = outer sphere, 2 = inner sphere.
DetectorModel MakeToyThreeSphereModel() {
    DetectorModel m;

    // Force identity transform (avoid “surprise quaternion” issues).
    m.SetDetectorOrigin(GeometryPosition(siren::math::Vector3D(0,0,0)));
    m.SetDetectorRotation(QFromZXZr(0.0, 0.0, 0.0));

    // Replace default sectors with an explicit finite world so intersections are well-behaved.
    m.ClearSectors();

    int mat = m.GetMaterials().GetMaterialId("VACUUM");

    DetectorSector world;
    world.name = "world";
    world.material_id = mat;
    world.level = 0;
    world.geo = Sphere(siren::math::Vector3D(0,0,0), 1000.0, 0.0).create();
    world.density = ConstantDensityDistribution(0.123).create();
    m.AddSector(world);

    DetectorSector outer;
    outer.name = "outer";
    outer.material_id = mat;
    outer.level = 1;
    outer.geo = Sphere(siren::math::Vector3D(0,0,0), 10.0, 0.0).create();
    outer.density = ConstantDensityDistribution(2.0).create();
    m.AddSector(outer);

    DetectorSector inner;
    inner.name = "inner";
    inner.material_id = mat;
    inner.level = 2;
    inner.geo = Sphere(siren::math::Vector3D(0,0,0), 5.0, 0.0).create();
    inner.density = ConstantDensityDistribution(7.0).create();
    m.AddSector(inner);

    return m;
}

int ExpectedToyLevel(double r, double eps = 1e-6) {
    if(r < 5.0 - eps) return 2;
    if(r < 10.0 - eps) return 1;
    return 0;
}

} // namespace

TEST(ContainingSector, DirectionInvarianceAwayFromBoundaries_ToyModel)
{
    DetectorModel m = MakeToyThreeSphereModel();
    std::mt19937 rng(0x12345678);

    // Directions should not matter away from boundaries.
    for(int i=0; i<400; ++i) {
        Vector3D p = RandomPointInBall(rng, 20.0);
        double r = p.magnitude();

        // Avoid boundary ambiguity.
        if(std::abs(r - 5.0) < 1e-3 || std::abs(r - 10.0) < 1e-3) { --i; continue; }

        DetectorSector expect = m.GetContainingSector(DetectorPosition(p));
        int expect_level = expect.level;

        for(int j=0; j<16; ++j) {
            Vector3D dir = RandomUnit(rng);
            auto isects = m.GetIntersections(DetectorPosition(p), DetectorDirection(dir));
            DetectorSector got = m.GetContainingSector(isects, DetectorPosition(p));
            EXPECT_EQ(got.level, expect_level)
                << "p=(" << p.GetX() << "," << p.GetY() << "," << p.GetZ() << ")"
                << " dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ() << ")";
        }
    }
}

TEST(ContainingSector, CachedIntersectionsAgreeWithMassDensity_ToyModel)
{
    DetectorModel m = MakeToyThreeSphereModel();
    std::mt19937 rng(0x12345678);

    // Use GetMassDensity(intersections, p) as an oracle for which sector contains p,
    // because GetMassDensity already implements the robust (<=, >=) interval check.
    for(int i=0; i<200; ++i) {
        Vector3D pref = RandomPointInBall(rng, 8.0);  // keep it well inside outer sphere often
        Vector3D dir  = RandomUnit(rng);

        auto isects = m.GetIntersections(DetectorPosition(pref), DetectorDirection(dir));

        std::uniform_real_distribution<double> ut(-30.0, 30.0);
        for(int k=0; k<64; ++k) {
            double t = ut(rng);
            Vector3D p = pref + dir * t;

            double r = p.magnitude();
            if(r > 999.0) continue; // stay inside world

            double rho = m.GetMassDensity(isects, DetectorPosition(p));
            DetectorSector s = m.GetContainingSector(isects, DetectorPosition(p));

            ASSERT_TRUE((bool)s.density) << "GetContainingSector returned an uninitialized sector";

            // Constant density distributions -> should match to tight tolerance.
            double rho_s = s.density->Evaluate(m.ToGeo(DetectorPosition(p)));
            EXPECT_NEAR(rho_s, rho, 1e-12)
                << "pref=(" << pref.GetX() << "," << pref.GetY() << "," << pref.GetZ() << ")"
                << " p=("    << p.GetX()    << "," << p.GetY()    << "," << p.GetZ()    << ")"
                << " dir=("  << dir.GetX()  << "," << dir.GetY()  << "," << dir.GetZ()  << ")"
                << " t=" << t
                << " sector=" << s.name << " level=" << s.level;
        }
    }
}

TEST(ContainingSector, CachedIntersectionsMatchesLocalForPointsOnLine_ToyModel)
{
    DetectorModel m = MakeToyThreeSphereModel();

    Vector3D dir(0.3, -0.4, 0.5);
    dir.normalize();

    // Reference point outside the outer sphere, along +dir
    Vector3D pref = dir * 12.0;

    // Cached intersections from pref along dir
    auto cached = m.GetIntersections(DetectorPosition(pref), DetectorDirection(dir));

    const double eps = 1e-6;

    // ds < 0 means "behind" pref w.r.t. +dir, which is where the old bug fails (dot was always +1).
    std::vector<double> ds = {
        +1.0,
        -1.0,

        // cross outer boundary at |12+ds| = 10  -> ds = -2 (and ds = -22 on the far side)
        -2.0 - eps, -2.0 + eps,
        -22.0 - eps, -22.0 + eps,

        // cross inner boundary at |12+ds| = 5   -> ds = -7 (and ds = -17 on the far side)
        -7.0 - eps, -7.0 + eps,
        -17.0 - eps, -17.0 + eps,

        // deep points in each region
        -8.0,   // inside inner
        -3.0,   // outer shell
        -25.0   // vacuum
    };

    for(double d : ds) {
        Vector3D p = pref + dir * d;

        // "Local" truth: compute intersections at p itself
        auto local = m.GetIntersections(DetectorPosition(p), DetectorDirection(dir));
        DetectorSector expect = m.GetContainingSector(local, DetectorPosition(p));

        // What we actually care about: reusing cached intersections
        DetectorSector got = m.GetContainingSector(cached, DetectorPosition(p));

        EXPECT_EQ(got.level, expect.level) << "ds=" << d;
    }
}


TEST(ContainingSector, ExitBoundaryPointFromIntersectionIsHandled_ToyModel)
{
    DetectorModel m = MakeToyThreeSphereModel();
    std::mt19937 rng(0x12345678);

    // This test is designed to catch the strict `end_point > 0` bug.
    // We construct p exactly at an *exit* intersection of the inner sphere along a random direction.
    // For that segment, end_point should be ~0 and GetContainingSector must still agree with GetMassDensity.
    Vector3D pref(0,0,0); // inside inner sphere

    for(int i=0; i<256; ++i) {
        Vector3D dir = RandomUnit(rng);
        auto isects = m.GetIntersections(DetectorPosition(pref), DetectorDirection(dir));

        // Find the inner-sphere exit intersection in the forward direction (distance > 0, entering==false).
        double d_exit = std::numeric_limits<double>::quiet_NaN();
        for(auto const & it : isects.intersections) {
            if(it.hierarchy == 2 && !it.entering && std::isfinite(it.distance) && it.distance > 0.0) {
                d_exit = it.distance;
                break;
            }
        }
        ASSERT_TRUE(std::isfinite(d_exit));

        Vector3D p_exit = pref + dir * d_exit;

        // Sanity: p_exit should be on the r=5 sphere up to floating error.
        // We don't REQUIRE exactness; we WANT tiny floating error here.
        ASSERT_NEAR(p_exit.magnitude(), 5.0, 1e-8);

        double rho = m.GetMassDensity(isects, DetectorPosition(p_exit));
        DetectorSector s = m.GetContainingSector(isects, DetectorPosition(p_exit));

        ASSERT_TRUE((bool)s.density) << "GetContainingSector returned an uninitialized sector at exit boundary";

        double rho_s = s.density->Evaluate(m.ToGeo(DetectorPosition(p_exit)));

        EXPECT_NEAR(rho_s, rho, 1e-12)
            << "dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ() << ")"
            << " d_exit=" << d_exit
            << " p_exit=(" << p_exit.GetX() << "," << p_exit.GetY() << "," << p_exit.GetZ() << ")"
            << " sector=" << s.name << " level=" << s.level;
    }
}

TEST(ContainingSector, EpsilonAroundBoundariesIsStableAcrossDirections_ToyModel)
{
    DetectorModel m = MakeToyThreeSphereModel();
    std::mt19937 rng(0x12345678);

    // Two points: just inside and just outside the inner boundary at r=5, along x-axis.
    // We check that direction choice doesn't flip the classification.
    const double eps = 1e-6;
    std::vector<Vector3D> pts = {
        Vector3D(5.0 - eps, 0, 0), // inside inner
        Vector3D(5.0 + eps, 0, 0), // outside inner (but inside outer)
        Vector3D(10.0 - eps, 0, 0), // inside outer
        Vector3D(10.0 + eps, 0, 0), // outside outer (world)
    };

    for(auto const & p : pts) {
        int expect = ExpectedToyLevel(p.magnitude(), 1e-7);

        for(int j=0; j<64; ++j) {
            Vector3D dir = RandomUnit(rng);
            auto isects = m.GetIntersections(DetectorPosition(p), DetectorDirection(dir));
            DetectorSector s = m.GetContainingSector(isects, DetectorPosition(p));
            EXPECT_EQ(s.level, expect)
                << "p=(" << p.GetX() << "," << p.GetY() << "," << p.GetZ() << ")"
                << " r=" << p.magnitude()
                << " dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ() << ")"
                << " got=" << s.level << " (" << s.name << ")"
                << " expect=" << expect;
        }
    }
}

