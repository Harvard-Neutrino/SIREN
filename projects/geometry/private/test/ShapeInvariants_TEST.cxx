// Shape invariant tests
//
// Validates geometric properties that must hold for all shapes:
// 1. Intersection count (convex: 0 or 2; hollow: 0,2,4; CSG: even)
// 2. Enter/exit ordering (strictly alternating)
// 3. Reciprocity (reversed ray hits the same surface points)
// 4. AABB consistency (intersections lie within the world bounding box)
// 5. AABB tightness (surface samples lie within the bounding box)
// 6. Distance-position consistency (position == origin + t*direction)
// 7. Monte Carlo cross-section (hit fraction matches expected area ratio)

#include <cmath>
#include <random>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <functional>

#include <gtest/gtest.h>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/ExtrPoly.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Torus.h"
#include "SIREN/geometry/EllipticalTube.h"
#include "SIREN/geometry/CutTube.h"
#include "SIREN/geometry/Trap.h"
#include "SIREN/geometry/Ellipsoid.h"
#include "SIREN/geometry/Para.h"
#include "SIREN/geometry/GenericPolycone.h"
#include "SIREN/geometry/GeometryMesh.h"

using namespace siren::geometry;
using namespace siren::math;

namespace {

std::mt19937 rng_(54321);
std::uniform_real_distribution<double> uniform(-1.0, 1.0);

Vector3D RandomPoint(double scale) {
    return Vector3D(uniform(rng_) * scale, uniform(rng_) * scale, uniform(rng_) * scale);
}

Vector3D RandomDirection() {
    double dx, dy, dz, mag;
    do {
        dx = uniform(rng_); dy = uniform(rng_); dz = uniform(rng_);
        mag = std::sqrt(dx*dx + dy*dy + dz*dz);
    } while(mag < 1e-12);
    return Vector3D(dx/mag, dy/mag, dz/mag);
}

// Registry of test shapes with metadata
struct ShapeEntry {
    std::string name;
    std::shared_ptr<Geometry> geo;
    int max_intersections; // max allowed per ray (2 for convex, 4 for hollow, -1 for any even)
    double known_cross_section; // for Monte Carlo test, -1 if unknown
    double mc_scale; // sampling box half-width for MC test
};

std::vector<ShapeEntry> MakeShapes() {
    std::vector<ShapeEntry> shapes;

    // Convex solids (must return 0 or 2 intersections)
    shapes.push_back({"Box(10,8,6)", Box(10, 8, 6).create(), 2, -1, 0});
    shapes.push_back({"Sphere(5,0)", Sphere(5, 0).create(), 2, -1, 0});
    shapes.push_back({"Cylinder(4,0,10)", Cylinder(4, 0, 10).create(), 2, -1, 0});
    shapes.push_back({"Cone(0,5,0,3,8)", Cone(0, 5, 0, 3, 8).create(), 2, -1, 0});
    shapes.push_back({"Trd(5,3,4,2,6)", Trd(5, 3, 4, 2, 6).create(), 2, -1, 0});

    // Pointed cones (apex at one end, #21)
    shapes.push_back({"ConePointed(0,5,0,0,8)", Cone(0, 5, 0, 0, 8).create(), 2, -1, 0});
    shapes.push_back({"ConePointedReverse(0,0,0,5,8)", Cone(0, 0, 0, 5, 8).create(), 2, -1, 0});

    // Extruded polygon (#20)
    {
        std::vector<std::vector<double>> polygon = {{-3,-3},{3,-3},{3,3},{-3,3}};
        double off[2] = {0, 0};
        std::vector<ExtrPoly::ZSection> zsecs = {
            ExtrPoly::ZSection(-4, off, 1.0),
            ExtrPoly::ZSection(4, off, 1.0)
        };
        shapes.push_back({"ExtrPoly(square)", ExtrPoly(polygon, zsecs).create(), 2, -1, 0});
    }
    {
        std::vector<std::vector<double>> polygon = {{-3,0},{0,-3},{3,0},{0,3}};
        double off[2] = {0, 0};
        std::vector<ExtrPoly::ZSection> zsecs = {
            ExtrPoly::ZSection(-2, off, 1.5),
            ExtrPoly::ZSection(2, off, 0.5)
        };
        shapes.push_back({"ExtrPoly(diamond,tapered)", ExtrPoly(polygon, zsecs).create(), 2, -1, 0});
    }

    // Hollow shapes (0, 2, or 4 intersections)
    shapes.push_back({"SphereHollow(5,2)", Sphere(5, 2).create(), 4, -1, 0});
    shapes.push_back({"CylinderHollow(5,2,10)", Cylinder(5, 2, 10).create(), 4, -1, 0});
    shapes.push_back({"ConeHollow(1,5,1,3,8)", Cone(1, 5, 1, 3, 8).create(), 4, -1, 0});

    // Polyhedra variants (#20 triangle, #21 hollow, zero-radius apex)
    shapes.push_back({"Polyhedra3(triangle)", Polyhedra(3, 0, {-4,4}, {0,0}, {5,5}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra6Hollow", Polyhedra(6, 0, {-5,5}, {2,2}, {5,5}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra4Pyramid", Polyhedra(4, M_PI/4, {-3,5}, {0,0}, {5,0}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra6PyramidRev", Polyhedra(6, 0, {-4,4}, {0,0}, {0,6}).create(), -1, -1, 0});

    // Multi-segment shapes
    shapes.push_back({"Polycone", Polycone({-5,-2,0,3,5}, {0,0,0,0,0}, {3,5,4,6,2}).create(), -1, -1, 0});
    shapes.push_back({"PolyconeHollow", Polycone({-4,0,4}, {1,2,1}, {5,6,5}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra6", Polyhedra(6, 0, {-5,0,5}, {0,0,0}, {4,6,4}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra4", Polyhedra(4, M_PI/4, {-3,3}, {0,0}, {5,5}).create(), -1, -1, 0});

    // Rotated-placement shapes (#22)
    {
        Quaternion q(std::cos(0.5), std::sin(0.5)*0.577, std::sin(0.5)*0.577, std::sin(0.5)*0.577);
        shapes.push_back({"ConeRotated", Cone(Placement(Vector3D(1,2,3), q), 0, 5, 0, 3, 8).create(), 2, -1, 0});
        shapes.push_back({"TrdRotated", Trd(Placement(Vector3D(-1,1,-2), q), 5, 3, 4, 2, 6).create(), 2, -1, 0});
        shapes.push_back({"PolyconeRotated", Polycone(Placement(Vector3D(2,0,-1), q), {-4,0,4}, {0,0,0}, {3,5,3}).create(), -1, -1, 0});
        shapes.push_back({"PolyhedraRotated", Polyhedra(Placement(Vector3D(0,-2,1), q), 6, 0, {-3,3}, {0,0}, {4,4}).create(), -1, -1, 0});
    }

    // Boolean CSG (any even number of intersections)
    {
        auto s1 = Sphere(Placement(Vector3D(-2,0,0)), 5, 0).create();
        auto s2 = Sphere(Placement(Vector3D(2,0,0)), 5, 0).create();
        shapes.push_back({"BoolUnion", BooleanGeometry(BooleanOperation::UNION, s1, s2).create(), -1, -1, 0});
    }
    {
        auto sp = Sphere(5, 0).create();
        auto bx = Box(8, 8, 8).create();
        shapes.push_back({"BoolIntersect", BooleanGeometry(BooleanOperation::INTERSECTION, sp, bx).create(), -1, -1, 0});
    }
    {
        auto sp = Sphere(5, 0).create();
        auto cy = Cylinder(2, 0, 12).create();
        shapes.push_back({"BoolSubtract", BooleanGeometry(BooleanOperation::SUBTRACTION, sp, cy).create(), -1, -1, 0});
    }

    // MC test shapes: sphere has known cross-section (pi*r^2 circle from any direction through center)
    shapes.push_back({"Sphere(r=5)_MC", Sphere(5, 0).create(), 2, M_PI * 25, 8});
    // Box has known cross-section along z-axis: x*y = 10*8 = 80
    shapes.push_back({"Box(10,8,6)_MC", Box(10, 8, 6).create(), 2, 80, 8});
    // Cylinder cross-section along z-axis: pi*r^2 = pi*16
    shapes.push_back({"Cylinder(r=4,z=10)_MC", Cylinder(4, 0, 10).create(), 2, M_PI * 16, 8});

    // Torus shapes (max 4 intersections for solid, -1 for hollow)
    shapes.push_back({"Torus(R=10,r=3,0)", Torus(10, 3, 0).create(), -1, -1, 0});
    shapes.push_back({"TorusHollow(R=10,3,1)", Torus(10, 3, 1).create(), -1, -1, 0});
    shapes.push_back({"TorusThin(R=8,r=0.5,0)", Torus(8, 0.5, 0).create(), -1, -1, 0});
    // Self-intersecting torus (rtor < rmax): tube passes through center
    shapes.push_back({"TorusSelfIntersecting(R=3,r=5,0)", Torus(3, 5, 0).create(), -1, -1, 0});
    shapes.push_back({"TorusSelfIntersectingHollow(R=3,5,2)", Torus(3, 5, 2).create(), -1, -1, 0});

    // Partial spheres
    shapes.push_back({"SphereThetaHemi", Sphere(5, 0, 0, 2*M_PI, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SpherePhiHalf", Sphere(5, 0, 0, M_PI, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"SpherePhiTheta", Sphere(5, 0, 0, M_PI, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SphereHollowThetaHemi", Sphere(5, 2, 0, 2*M_PI, 0, M_PI/2).create(), -1, -1, 0});

    // Partial torus
    shapes.push_back({"TorusQuarter", Torus(10, 3, 0, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"TorusHalf", Torus(10, 3, 0, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"TorusHollowQuarter", Torus(10, 3, 1, 0, M_PI/2).create(), -1, -1, 0});

    // Polycone with step-change (#8)
    shapes.push_back({"PolyconeStep", Polycone({-0.05, 0.05, 0.05, 0.15}, {0, 0, 0, 0}, {0.03, 0.03, 0.05, 0.05}).create(), -1, -1, 0});

    // GenericPolycone variants
    shapes.push_back({"GenericPolyconeSolid", GenericPolycone({0, 5, 0}, {-5, 0, 5}).create(), -1, -1, 0});
    shapes.push_back({"GenericPolyconeHollow", GenericPolycone({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5}).create(), -1, -1, 0});
    shapes.push_back({"GenericPolyconeUShape", GenericPolycone({5, 5, 3, 3}, {-4, 4, 4, -4}).create(), -1, -1, 0});
    shapes.push_back({"GenericPolyconeCrystal", GenericPolycone({5,38,40,40,30,30,15,15,0,0,5}, {0,0,70,80,80,77,77,80,80,60,60}).create(), -1, -1, 0});
    {
        Quaternion q(std::cos(0.5), std::sin(0.5)*0.577, std::sin(0.5)*0.577, std::sin(0.5)*0.577);
        shapes.push_back({"GenericPolyconeRotated", GenericPolycone(Placement(Vector3D(1,-1,2), q), {3,4.5,5,3.5,3,2}, {-5,0,5,5,0,-5}).create(), -1, -1, 0});
    }

    // ExtrPoly with offset
    {
        std::vector<std::vector<double>> polygon = {{-2,-2},{2,-2},{2,2},{-2,2}};
        double off_bot[2] = {0, 0};
        double off_top[2] = {1, 1};
        std::vector<ExtrPoly::ZSection> zsecs = {
            ExtrPoly::ZSection(-3, off_bot, 1.0),
            ExtrPoly::ZSection(3, off_top, 1.0)
        };
        shapes.push_back({"ExtrPoly(offset)", ExtrPoly(polygon, zsecs).create(), 2, -1, 0});
    }

    // ExtrPoly with multiple z-sections
    {
        std::vector<std::vector<double>> polygon = {{-2,-2},{2,-2},{2,2},{-2,2}};
        double off[2] = {0, 0};
        std::vector<ExtrPoly::ZSection> zsecs = {
            ExtrPoly::ZSection(-4, off, 1.0),
            ExtrPoly::ZSection(0, off, 1.5),
            ExtrPoly::ZSection(4, off, 0.8)
        };
        shapes.push_back({"ExtrPoly(multiZ)", ExtrPoly(polygon, zsecs).create(), -1, -1, 0});
    }

    // New solid types
    shapes.push_back({"EllipticalTube(3,2,5)", EllipticalTube(3, 2, 5).create(), 2, -1, 0});
    shapes.push_back({"CutTubeFlat(0,5,8,(0,0,-1),(0,0,1))", CutTube(0, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1)).create(), 2, -1, 0});
    shapes.push_back({"CutTubeHollow(2,5,8,(0,0,-1),(0,0,1))", CutTube(2, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1)).create(), 4, -1, 0});
    shapes.push_back({"CutTubeTilted(0,5,8)", CutTube(0, 5, 8, Vector3D(0,-0.2,-0.98), Vector3D(0.1,0,0.995)).create(), 2, -1, 0});
    shapes.push_back({"Trap(symmetric)", Trap(5, 0, 0, 3, 4, 4, 0, 2, 3, 3, 0).create(), 2, -1, 0});
    shapes.push_back({"Trap(general)", Trap(5, 0.1, 0.2, 3, 4, 5, 0.15, 2, 3, 4, 0.1).create(), 2, -1, 0});
    shapes.push_back({"Ellipsoid(5,3,4)", Ellipsoid(5, 3, 4).create(), 2, -1, 0});
    shapes.push_back({"Ellipsoid(5,3,4,-2,3)", Ellipsoid(5, 3, 4, -2, 3).create(), 4, -1, 0});
    shapes.push_back({"Para(4,3,5,0,0,0)", Para(4, 3, 5, 0, 0, 0).create(), 2, -1, 0});
    shapes.push_back({"Para(4,3,5,0.3,0.2,0.5)", Para(4, 3, 5, 0.3, 0.2, 0.5).create(), 2, -1, 0});

    return shapes;
}

} // anonymous namespace

// =========================================================================
// Test 1: Intersection count invariants
// =========================================================================
TEST(ShapeInvariants, IntersectionCount) {
    auto shapes = MakeShapes();
    int N = 10000;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue; // skip MC shapes
        int violations = 0;
        for(int i = 0; i < N; ++i) {
            Vector3D pos = RandomPoint(20);
            Vector3D dir = RandomDirection();
            auto isects = entry.geo->Intersections(pos, dir);
            int n = (int)isects.size();

            // Must be even
            if(n % 2 != 0) {
                violations++;
                if(violations <= 3) {
                    std::cerr << entry.name << ": odd intersection count " << n
                              << " at (" << pos.GetX() << "," << pos.GetY() << "," << pos.GetZ() << ")" << std::endl;
                }
                continue;
            }

            // For shapes with known max, check upper bound
            if(entry.max_intersections > 0 && n > entry.max_intersections) {
                violations++;
                if(violations <= 3) {
                    std::cerr << entry.name << ": " << n << " intersections (max " << entry.max_intersections << ")" << std::endl;
                }
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " intersection count violations";
    }
}

// =========================================================================
// Test 2: Enter/exit ordering
// =========================================================================
TEST(ShapeInvariants, EnterExitOrdering) {
    auto shapes = MakeShapes();
    int N = 10000;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;
        int violations = 0;
        for(int i = 0; i < N; ++i) {
            Vector3D pos = RandomPoint(20);
            Vector3D dir = RandomDirection();
            auto isects = entry.geo->Intersections(pos, dir);

            // Filter to just the forward intersections (positive distance) and check alternation
            bool last_entering = false;
            bool first = true;
            for(auto const & isect : isects) {
                if(!first) {
                    // Consecutive same-type is a violation
                    if(isect.entering == last_entering) {
                        violations++;
                        if(violations <= 3) {
                            std::cerr << entry.name << ": consecutive "
                                      << (isect.entering ? "entering" : "exiting")
                                      << " at dist=" << isect.distance << std::endl;
                        }
                        break;
                    }
                }
                first = false;
                last_entering = isect.entering;
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " enter/exit ordering violations";
    }
}

// =========================================================================
// Test 3: Ray reversal reciprocity
// =========================================================================
TEST(ShapeInvariants, Reciprocity) {
    auto shapes = MakeShapes();
    int N = 5000;
    double tol = 1e-6;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;
        int violations = 0;
        for(int i = 0; i < N; ++i) {
            Vector3D pos = RandomPoint(20);
            Vector3D dir = RandomDirection();
            Vector3D neg_dir(-dir.GetX(), -dir.GetY(), -dir.GetZ());

            auto fwd = entry.geo->Intersections(pos, dir);
            auto rev = entry.geo->Intersections(pos, neg_dir);

            // Both must find the same number of intersections
            if(fwd.size() != rev.size()) {
                violations++;
                if(violations <= 3) {
                    std::cerr << entry.name << ": forward=" << fwd.size()
                              << " reverse=" << rev.size() << std::endl;
                }
                continue;
            }

            // The intersection positions should match (reversed order, same points)
            // Forward sorted by distance ascending; reverse has negated distances
            for(size_t j = 0; j < fwd.size(); ++j) {
                size_t k = fwd.size() - 1 - j;
                double dx = fwd[j].position.GetX() - rev[k].position.GetX();
                double dy = fwd[j].position.GetY() - rev[k].position.GetY();
                double dz = fwd[j].position.GetZ() - rev[k].position.GetZ();
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                if(dist > tol) {
                    violations++;
                    if(violations <= 3) {
                        std::cerr << entry.name << ": position mismatch at pair " << j
                                  << " dist=" << dist << std::endl;
                    }
                    break;
                }
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " reciprocity violations";
    }
}

// =========================================================================
// Test 4: AABB consistency (intersections lie within bounding box)
// =========================================================================
TEST(ShapeInvariants, AABBConsistency) {
    auto shapes = MakeShapes();
    int N = 10000;
    double tol = 1e-6; // allow slight overshoot due to floating point

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;
        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid()) continue;
        // Skip shapes with infinite bounds
        if(!std::isfinite(box.min_corner.GetX())) continue;

        int violations = 0;
        for(int i = 0; i < N; ++i) {
            Vector3D pos = RandomPoint(20);
            Vector3D dir = RandomDirection();
            auto isects = entry.geo->Intersections(pos, dir);

            for(auto const & isect : isects) {
                double x = isect.position.GetX();
                double y = isect.position.GetY();
                double z = isect.position.GetZ();
                if(x < box.min_corner.GetX() - tol || x > box.max_corner.GetX() + tol ||
                   y < box.min_corner.GetY() - tol || y > box.max_corner.GetY() + tol ||
                   z < box.min_corner.GetZ() - tol || z > box.max_corner.GetZ() + tol) {
                    violations++;
                    if(violations <= 3) {
                        std::cerr << entry.name << ": intersection at ("
                                  << x << "," << y << "," << z << ") outside AABB ["
                                  << box.min_corner.GetX() << ".." << box.max_corner.GetX() << "]["
                                  << box.min_corner.GetY() << ".." << box.max_corner.GetY() << "]["
                                  << box.min_corner.GetZ() << ".." << box.max_corner.GetZ() << "]" << std::endl;
                    }
                }
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " AABB consistency violations";
    }
}

// =========================================================================
// Test 5: AABB tightness (surface samples lie within bounding box)
// =========================================================================
TEST(ShapeInvariants, AABBTightness) {
    auto shapes = MakeShapes();
    int N = 10000;
    double tol = 1e-6;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;
        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid()) continue;
        if(!std::isfinite(box.min_corner.GetX())) continue;

        // Shoot rays and collect surface points; verify they're inside the AABB
        int violations = 0;
        for(int i = 0; i < N; ++i) {
            // Shoot from inside the AABB to ensure we hit the surface
            Vector3D center = box.Centroid();
            Vector3D dir = RandomDirection();
            auto isects = entry.geo->Intersections(center, dir);

            for(auto const & isect : isects) {
                double x = isect.position.GetX();
                double y = isect.position.GetY();
                double z = isect.position.GetZ();
                if(x < box.min_corner.GetX() - tol || x > box.max_corner.GetX() + tol ||
                   y < box.min_corner.GetY() - tol || y > box.max_corner.GetY() + tol ||
                   z < box.min_corner.GetZ() - tol || z > box.max_corner.GetZ() + tol) {
                    violations++;
                    if(violations <= 3) {
                        std::cerr << entry.name << ": surface point outside AABB" << std::endl;
                    }
                }
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " AABB tightness violations";
    }
}

// =========================================================================
// Test 6: Distance-position consistency (pos == origin + t*dir)
// =========================================================================
TEST(ShapeInvariants, DistancePositionConsistency) {
    auto shapes = MakeShapes();
    int N = 10000;
    double tol = 1e-6;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;
        int violations = 0;
        for(int i = 0; i < N; ++i) {
            Vector3D pos = RandomPoint(20);
            Vector3D dir = RandomDirection();
            auto isects = entry.geo->Intersections(pos, dir);

            for(auto const & isect : isects) {
                // Expected position: pos + t * dir
                double ex = pos.GetX() + isect.distance * dir.GetX();
                double ey = pos.GetY() + isect.distance * dir.GetY();
                double ez = pos.GetZ() + isect.distance * dir.GetZ();

                double dx = isect.position.GetX() - ex;
                double dy = isect.position.GetY() - ey;
                double dz = isect.position.GetZ() - ez;
                double err = std::sqrt(dx*dx + dy*dy + dz*dz);

                if(err > tol) {
                    violations++;
                    if(violations <= 3) {
                        std::cerr << entry.name << ": position error " << err
                                  << " at dist=" << isect.distance << std::endl;
                    }
                }
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " distance-position violations";
    }
}

// =========================================================================
// Test 7: Monte Carlo cross-section
// =========================================================================
TEST(ShapeInvariants, MonteCarloCrossSection) {
    auto shapes = MakeShapes();
    int N = 100000;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section <= 0) continue;

        double half = entry.mc_scale;
        double sample_area = (2 * half) * (2 * half);
        int hits = 0;

        // Shoot parallel rays along z-axis from a square grid above the shape
        std::uniform_real_distribution<double> pos_dist(-half, half);
        for(int i = 0; i < N; ++i) {
            double x = pos_dist(rng_);
            double y = pos_dist(rng_);
            Vector3D pos(x, y, 20); // above the shape
            Vector3D dir(0, 0, -1); // shooting downward
            auto isects = entry.geo->Intersections(pos, dir);
            if(!isects.empty()) hits++;
        }

        double measured_area = sample_area * (double)hits / N;
        double expected_area = entry.known_cross_section;
        double relative_error = std::fabs(measured_area - expected_area) / expected_area;

        // Allow 5% tolerance for MC statistical error at 100k samples
        EXPECT_LT(relative_error, 0.05)
            << entry.name << ": measured=" << measured_area
            << " expected=" << expected_area
            << " rel_error=" << relative_error;
    }
}

// =========================================================================
// Test 8: Surface-boundary points (#23)
// Rays originating exactly on a surface with various directions
// =========================================================================
TEST(ShapeInvariants, SurfaceBoundaryPoints) {
    // Test shapes with known surface points
    struct SurfaceTest {
        std::string name;
        std::shared_ptr<Geometry> geo;
        Vector3D surface_point;
    };

    std::vector<SurfaceTest> tests = {
        {"Sphere(5) at +x", Sphere(5, 0).create(), Vector3D(5, 0, 0)},
        {"Sphere(5) at +y", Sphere(5, 0).create(), Vector3D(0, 5, 0)},
        {"Sphere(5) at +z", Sphere(5, 0).create(), Vector3D(0, 0, 5)},
        {"Box(10,8,6) at +x face", Box(10, 8, 6).create(), Vector3D(5, 0, 0)},
        {"Box(10,8,6) at edge", Box(10, 8, 6).create(), Vector3D(5, 4, 0)},
        {"Cylinder(4,0,10) barrel", Cylinder(4, 0, 10).create(), Vector3D(4, 0, 0)},
        {"Cylinder(4,0,10) top", Cylinder(4, 0, 10).create(), Vector3D(0, 0, 5)},
        {"Cone(0,5,0,3,8) barrel mid", Cone(0, 5, 0, 3, 8).create(), Vector3D(4, 0, 0)},
        {"Trd(5,3,4,2,6) top face", Trd(5, 3, 4, 2, 6).create(), Vector3D(0, 0, 6)},
    };

    for(auto const & t : tests) {
        // Shoot rays in several directions from the surface point
        Vector3D dirs[] = {
            Vector3D(1, 0, 0), Vector3D(-1, 0, 0),
            Vector3D(0, 1, 0), Vector3D(0, 0, 1),
            Vector3D(0, 0, -1),
        };
        for(auto const & dir : dirs) {
            auto isects = t.geo->Intersections(t.surface_point, dir);
            // Intersection count must be even
            EXPECT_EQ(isects.size() % 2, 0u)
                << t.name << " dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ()
                << "): odd count " << isects.size();
            // Enter/exit must alternate
            for(size_t i = 1; i < isects.size(); ++i) {
                EXPECT_NE(isects[i].entering, isects[i-1].entering)
                    << t.name << ": non-alternating at index " << i;
            }
        }
    }
}

// =========================================================================
// Test 9: Axis-aligned rays (#23)
// Exercise dz==0, dx==0, dy==0 branches in cap/face tests
// =========================================================================
TEST(ShapeInvariants, AxisAlignedRays) {
    auto shapes = MakeShapes();

    Vector3D dirs[] = {
        Vector3D(1, 0, 0), Vector3D(0, 1, 0), Vector3D(0, 0, 1),
        Vector3D(-1, 0, 0), Vector3D(0, -1, 0), Vector3D(0, 0, -1),
    };

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;
        for(auto const & dir : dirs) {
            // Shoot from several positions
            for(int i = 0; i < 1000; ++i) {
                Vector3D pos = RandomPoint(20);
                auto isects = entry.geo->Intersections(pos, dir);
                int n = (int)isects.size();

                EXPECT_EQ(n % 2, 0)
                    << entry.name << " axis-aligned: odd count " << n;

                // Check alternation
                for(int j = 1; j < n; ++j) {
                    EXPECT_NE(isects[j].entering, isects[j-1].entering)
                        << entry.name << " axis-aligned: non-alternating at " << j;
                }
            }
        }
    }
}

// =========================================================================
// Test 10: Tangent rays on curved surfaces (identified gap)
// =========================================================================
TEST(ShapeInvariants, TangentRays) {
    double eps = 1e-9;

    // --- Sphere of radius 5 ---
    auto sphere = Sphere(5, 0).create();
    {
        // Tangent ray at y=5 going in x-direction
        Vector3D tangent_pos(0, 5, 0);
        Vector3D tangent_dir(1, 0, 0);
        auto isects = sphere->Intersections(tangent_pos, tangent_dir);
        // Tangent: 0 or 2 intersections, both valid
        EXPECT_TRUE(isects.size() == 0 || isects.size() == 2)
            << "Sphere tangent ray: expected 0 or 2, got " << isects.size();
        EXPECT_EQ(isects.size() % 2, 0u) << "Sphere tangent: odd count";

        // Ray just outside (y = 5 + eps): should miss
        Vector3D outside_pos(0, 5 + eps, 0);
        auto outside_isects = sphere->Intersections(outside_pos, tangent_dir);
        EXPECT_EQ(outside_isects.size(), 0u)
            << "Sphere outside tangent: expected 0, got " << outside_isects.size();

        // Ray just inside (y = 5 - eps): should hit
        Vector3D inside_pos(0, 5 - eps, 0);
        auto inside_isects = sphere->Intersections(inside_pos, tangent_dir);
        EXPECT_EQ(inside_isects.size(), 2u)
            << "Sphere inside tangent: expected 2, got " << inside_isects.size();
        EXPECT_EQ(inside_isects.size() % 2, 0u) << "Sphere inside tangent: odd count";
    }

    // --- Cylinder of radius 4 ---
    auto cyl = Cylinder(4, 0, 10).create();
    {
        // Tangent ray at y=4, z=0 going in x-direction
        Vector3D tangent_pos(0, 4, 0);
        Vector3D tangent_dir(1, 0, 0);
        auto isects = cyl->Intersections(tangent_pos, tangent_dir);
        EXPECT_TRUE(isects.size() == 0 || isects.size() == 2)
            << "Cylinder tangent ray: expected 0 or 2, got " << isects.size();
        EXPECT_EQ(isects.size() % 2, 0u) << "Cylinder tangent: odd count";

        // Ray just outside (y = 4 + eps): should miss
        Vector3D outside_pos(0, 4 + eps, 0);
        auto outside_isects = cyl->Intersections(outside_pos, tangent_dir);
        EXPECT_EQ(outside_isects.size(), 0u)
            << "Cylinder outside tangent: expected 0, got " << outside_isects.size();

        // Ray just inside (y = 4 - eps): should hit
        Vector3D inside_pos(0, 4 - eps, 0);
        auto inside_isects = cyl->Intersections(inside_pos, tangent_dir);
        EXPECT_EQ(inside_isects.size(), 2u)
            << "Cylinder inside tangent: expected 2, got " << inside_isects.size();
        EXPECT_EQ(inside_isects.size() % 2, 0u) << "Cylinder inside tangent: odd count";
    }
}

// =========================================================================
// Fix 3: Horizontal rays at Polycone/Polyhedra section boundaries
// A horizontal ray exactly at an internal z-boundary must be claimed by
// exactly one section (half-open interval [z_lo, z_hi)). Such a ray DOES
// pierce the barrel surface (it crosses the cross-section disc/polygon at
// that z); the prior strict (z_lo, z_hi) lateral test wrongly dropped both
// boundary hits, returning 0 and reporting interior points as outside.
// =========================================================================
TEST(ShapeInvariants, HorizontalRayAtSectionBoundary) {
    // Polycone with internal boundary at z=0 where rmax=5. The z=0 cross
    // section is a disc of radius 5, so a horizontal ray along x at z=0
    // through the axis enters at x=-5 and exits at x=+5: exactly 2 hits.
    Polycone pc({-5, 0, 5}, {0, 0, 0}, {3, 5, 3});
    auto pc_hits = pc.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
    ASSERT_EQ(pc_hits.size(), 2u)
        << "Horizontal ray at internal z-boundary must produce 2 barrel hits";
    double pc_xmin = std::min(pc_hits[0].position.GetX(), pc_hits[1].position.GetX());
    double pc_xmax = std::max(pc_hits[0].position.GetX(), pc_hits[1].position.GetX());
    EXPECT_NEAR(pc_xmin, -5.0, 1e-6) << "Entry must lie on barrel at x=-5, z=0";
    EXPECT_NEAR(pc_xmax,  5.0, 1e-6) << "Exit must lie on barrel at x=+5, z=0";
    EXPECT_NEAR(pc_hits[0].position.GetZ(), 0.0, 1e-9);
    EXPECT_NEAR(pc_hits[1].position.GetZ(), 0.0, 1e-9);

    // Containment via a true horizontal ray at the z-boundary.
    EXPECT_TRUE(pc.IsInside(Vector3D(1, 0, 0), Vector3D(1, 0, 0)))
        << "Polycone interior point at z=0 boundary, horizontal ray";
    EXPECT_FALSE(pc.IsInside(Vector3D(6, 0, 0), Vector3D(1, 0, 0)))
        << "Polycone exterior point at z=0 boundary, horizontal ray";

    // Polyhedra (hexagon at z=0): a horizontal ray offset in y so it
    // crosses two lateral faces (not a vertex) of the convex cross
    // section exactly twice at the z=0 boundary.
    Polyhedra ph(6, 0, {-5, 0, 5}, {0, 0, 0}, {3, 5, 3});
    auto ph_hits = ph.Intersections(Vector3D(-20, 1, 0), Vector3D(1, 0, 0));
    EXPECT_EQ(ph_hits.size(), 2u)
        << "Polyhedra horizontal ray at z-boundary must produce 2 hits";
    EXPECT_TRUE(ph.IsInside(Vector3D(1, 1, 0), Vector3D(1, 0, 0)))
        << "Polyhedra interior point at z=0 boundary, horizontal ray";

    // Multi-section polycone: internal boundary at z=-2 where rmax=5.
    Polycone pc2({-5, -2, 0, 3, 5}, {0, 0, 0, 0, 0}, {3, 5, 4, 6, 2});
    auto pc2_hits = pc2.Intersections(Vector3D(-20, 0, -2), Vector3D(1, 0, 0));
    ASSERT_EQ(pc2_hits.size(), 2u)
        << "Horizontal ray at z=-2 boundary must produce 2 barrel hits";
    EXPECT_NEAR(std::min(pc2_hits[0].position.GetX(), pc2_hits[1].position.GetX()),
                -5.0, 1e-6);
    EXPECT_NEAR(std::max(pc2_hits[0].position.GetX(), pc2_hits[1].position.GetX()),
                 5.0, 1e-6);
    EXPECT_TRUE(pc2.IsInside(Vector3D(1, 0, -2), Vector3D(1, 0, 0)))
        << "Polycone interior point at z=-2 boundary, horizontal ray";
}

// =========================================================================
// Fix 4: Tangent rays on curved surfaces
// Tangent rays (determinant == 0) are treated as misses (no surface crossing).
// Near-tangent rays just inside the surface must still produce 2 hits.
// =========================================================================
TEST(ShapeInvariants, TangentRaysCone) {
    // Cone with outer radius 5 at z=0: tangent at y=5
    Cone cone(0, 5, 0, 5, 10); // cylinder-like cone
    Vector3D tangent_pos(0, 5, 0);
    Vector3D tangent_dir(1, 0, 0);
    auto isects = cone.Intersections(tangent_pos, tangent_dir);
    // Tangent: must be even (0 or 2)
    EXPECT_EQ(isects.size() % 2, 0u)
        << "Cone tangent: odd count " << isects.size();

    // Just outside: no hits
    Vector3D outside(0, 5 + 1e-6, 0);
    auto miss = cone.Intersections(outside, tangent_dir);
    EXPECT_EQ(miss.size(), 0u) << "Cone outside tangent should miss";

    // Just inside: 2 hits
    Vector3D inside(0, 5 - 1e-6, 0);
    auto hit = cone.Intersections(inside, tangent_dir);
    EXPECT_EQ(hit.size(), 2u) << "Cone inside tangent should hit twice";
}

TEST(ShapeInvariants, TangentRaysPolycone) {
    // Polycone section with outer radius 5: tangent at y=5
    Polycone pc({-5, 5}, {0, 0}, {5, 5});
    Vector3D tangent_pos(0, 5, 0);
    Vector3D tangent_dir(1, 0, 0);
    auto isects = pc.Intersections(tangent_pos, tangent_dir);
    EXPECT_EQ(isects.size() % 2, 0u)
        << "Polycone tangent: odd count " << isects.size();

    // Just inside: should hit
    Vector3D inside(0, 5 - 1e-6, 0);
    auto hit = pc.Intersections(inside, tangent_dir);
    EXPECT_EQ(hit.size(), 2u) << "Polycone inside tangent should hit twice";
}

// =========================================================================
// Torus quartic solver stability: various problematic ray configurations
// =========================================================================

TEST(ShapeInvariants, TorusTangentRays) {
    Torus torus(10, 3, 0);
    double R = 10, r = 3;
    int odd_count = 0;
    int total = 0;

    // Rays tangent to the tube-center circle in the z=0 plane at various angles
    for(int i = 0; i < 100; ++i) {
        double theta = 2.0 * M_PI * i / 100.0;
        double cx = R * std::cos(theta);
        double cy = R * std::sin(theta);
        // Tangent direction to the circle at this point
        double tx = -std::sin(theta);
        double ty = std::cos(theta);
        // Ray from outside the tube at the tube-center height (z=0)
        Vector3D pos(cx, cy, r + 0.5);
        Vector3D dir(tx, ty, 0);
        dir.normalize();
        auto hits = torus.Intersections(pos, dir);
        if(hits.size() % 2 != 0) odd_count++;
        total++;
    }
    EXPECT_EQ(odd_count, 0)
        << "Torus tangent rays: " << odd_count << "/" << total << " produced odd intersection count";
}

TEST(ShapeInvariants, TorusFarFieldRays) {
    Torus torus(10, 3, 0);
    double R = 10, r = 3;
    int wrong = 0;
    int total = 0;

    // Rays from far away that clearly pass through the tube center
    double distances[] = {100, 1000, 10000, 50000};
    for(double D : distances) {
        for(int i = 0; i < 50; ++i) {
            double theta = 2.0 * M_PI * i / 50.0;
            double cx = R * std::cos(theta);
            double cy = R * std::sin(theta);
            // Start far back along x, aimed at tube center
            Vector3D pos(-D, cy, 0);
            Vector3D dir(1, 0, 0);
            auto hits = torus.Intersections(pos, dir);
            total++;
            if(hits.size() % 2 != 0) {
                wrong++;
            } else if(hits.size() == 0 && std::fabs(cy) < R + r - 0.1 && std::fabs(cy) > R - r + 0.1) {
                // Should have hit the tube
                wrong++;
            }
        }
    }
    EXPECT_EQ(wrong, 0)
        << "Torus far-field rays: " << wrong << "/" << total << " incorrect";
}

TEST(ShapeInvariants, TorusSymmetryPlaneRays) {
    Torus torus(10, 3, 0);
    int odd_count = 0;
    int total = 0;

    // Rays in the z=0 plane through various points
    for(int i = 0; i < 200; ++i) {
        double angle = 2.0 * M_PI * i / 200.0;
        Vector3D pos(0, 0, 0);
        Vector3D dir(std::cos(angle), std::sin(angle), 0);
        auto hits = torus.Intersections(pos, dir);
        if(hits.size() % 2 != 0) odd_count++;
        total++;
    }
    EXPECT_EQ(odd_count, 0)
        << "Torus z=0 plane rays: " << odd_count << "/" << total << " produced odd count";
}

TEST(ShapeInvariants, TorusZAxisRay) {
    // Ray along z-axis through the center hole of the torus (rxy=0).
    // Solid torus R=10, r=3 occupies rxy in [7,13]; the z-axis (rxy=0)
    // never reaches the body, so the count is exactly zero (knowable).
    Torus torus(10, 3, 0);
    Vector3D pos(0, 0, -20);
    Vector3D dir(0, 0, 1);
    auto hits = torus.Intersections(pos, dir);
    EXPECT_EQ(hits.size(), 0u) << "Z-axis ray must produce no intersections";
}

// Accepted hits must lie on the implicit torus surface
// F = (sqrt(x^2+y^2) - R)^2 + z^2 - r^2 = 0. Before root polishing, the
// Ferrari solver's positional error scaled as R^2/r and no test ever
// checked surface-landing.
TEST(ShapeInvariants, TorusHitOnSurface) {
    auto surfaceF = [](double R, double r, Vector3D const & p) {
        double x = p.GetX(), y = p.GetY(), z = p.GetZ();
        double rxy = std::sqrt(x*x + y*y);
        return (rxy - R) * (rxy - R) + z*z - r*r;
    };

    // Case 1: axis-aligned ray in the z=0 plane crosses both tube lobes;
    // surface crossings are at |x| = R-r = 7 and R+r = 13 (4 hits).
    {
        double R = 10, r = 3;
        Torus torus(R, r, 0);
        auto hits = torus.Intersections(Vector3D(-50, 0, 0), Vector3D(1, 0, 0));
        ASSERT_EQ(hits.size(), 4u) << "z=0 axis ray must cross 4 surfaces";
        for(size_t i = 0; i < hits.size(); ++i) {
            EXPECT_NEAR(surfaceF(R, r, hits[i].position), 0.0, 1e-7)
                << "hit " << i << " is off the torus surface";
            double ax = std::fabs(hits[i].position.GetX());
            EXPECT_TRUE(std::fabs(ax - 7.0) < 1e-6 || std::fabs(ax - 13.0) < 1e-6)
                << "hit " << i << " x=" << hits[i].position.GetX();
        }
    }

    // Case 2: vertical ray through the +x tube center: 2 hits at z = +/-r.
    {
        double R = 10, r = 3;
        Torus torus(R, r, 0);
        auto hits = torus.Intersections(Vector3D(10, 0, -20), Vector3D(0, 0, 1));
        ASSERT_EQ(hits.size(), 2u);
        for(size_t i = 0; i < hits.size(); ++i)
            EXPECT_NEAR(surfaceF(R, r, hits[i].position), 0.0, 1e-7)
                << "lobe hit " << i << " off surface";
    }

    // Case 3: oblique general ray (all direction components nonzero) to
    // exercise the non-biquadratic Ferrari path. Must hit and land on
    // the surface; assert it is not vacuous.
    {
        double R = 10, r = 3;
        Torus torus(R, r, 0);
        Vector3D dir(1.0, 0.1, 0.05);
        dir.normalize();
        auto hits = torus.Intersections(Vector3D(-40, 2, -1), dir);
        EXPECT_EQ(hits.size() % 2, 0u) << "oblique ray odd count";
        EXPECT_GE(hits.size(), 2u) << "oblique ray must pierce the torus";
        for(size_t i = 0; i < hits.size(); ++i)
            EXPECT_NEAR(surfaceF(R, r, hits[i].position), 0.0, 1e-6)
                << "oblique hit " << i << " off surface";
    }

    // Case 4: thin large-R torus. The old tolerance 1e-4*(r^2+R^2)
    // accepted points ~0.5 units off here; the fix must keep hits on
    // the surface to ~1e-6.
    {
        double R = 100, r = 1;
        Torus torus(R, r, 0);
        auto hits = torus.Intersections(Vector3D(-200, 0, 0), Vector3D(1, 0, 0));
        ASSERT_EQ(hits.size(), 4u);
        for(size_t i = 0; i < hits.size(); ++i)
            EXPECT_NEAR(surfaceF(R, r, hits[i].position), 0.0, 1e-6)
                << "thin-torus hit " << i << " off surface";
    }
}

// =========================================================================
// Geometry::less() null safety for all shapes
// =========================================================================

TEST(ShapeInvariants, OperatorLessCrossType) {
    // Cross-type comparisons via operator< should not crash.
    // operator< dispatches to less() only for same-type pairs,
    // but the less() null guard is the last line of defense.
    Box box(10, 8, 6);
    Sphere sphere(5, 0);
    Cylinder cyl(4, 0, 10);
    Cone cone(0, 5, 0, 3, 8);
    Trd trd(5, 3, 4, 2, 6);
    Polycone pc({-5, 5}, {0, 0}, {5, 5});
    Polyhedra ph(6, 0, {-3, 3}, {0, 0}, {4, 4});
    Torus tor(10, 3, 0);
    auto boolean = BooleanGeometry(BooleanOperation::UNION,
        std::const_pointer_cast<const Geometry>(sphere.create()),
        std::const_pointer_cast<const Geometry>(box.create()));

    std::vector<Geometry*> shapes = {&box, &sphere, &cyl, &cone, &trd, &pc, &ph, &tor, &boolean};

    // All cross-type pairs should produce a consistent ordering without crashing
    for(size_t i = 0; i < shapes.size(); ++i) {
        for(size_t j = 0; j < shapes.size(); ++j) {
            bool a_lt_b = *shapes[i] < *shapes[j];
            bool b_lt_a = *shapes[j] < *shapes[i];
            if(i != j) {
                // Antisymmetry: can't have both a<b and b<a
                EXPECT_FALSE(a_lt_b && b_lt_a)
                    << "Ordering violation between " << i << " and " << j;
            }
        }
    }

    // Same-type comparison via operator< exercises less() with correct dynamic_cast
    Box box2(8, 8, 8);
    EXPECT_NO_THROW({
        bool result = box < box2;
        (void)result;
    });
    Sphere sphere2(3, 0);
    EXPECT_NO_THROW({
        bool result = sphere < sphere2;
        (void)result;
    });
}

// =========================================================================
// Cone apex rays: ray through the tip where the surface normal is undefined
// and the quadratic has coincident roots.
// =========================================================================
// Apex enter/exit semantics, not just parity. A ray crossing a cone tip
// has a knowable count AND a knowable entering flag at the apex: the apex
// hit is "entering" when the ray moves toward the body side of the tip.
// (Parity-only mod-2 here masked a hardcoded entering=false at the apex.)
TEST(ShapeInvariants, ConeApexRays) {
    // Pointed cone: base radius 5 at z=-5, apex at z=+5 (body below apex).
    Cone cone(0, 5, 0, 0, 10);

    // Up the z-axis: enter base cap (z=-5, entering), exit at apex (z=+5).
    {
        auto h = cone.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        ASSERT_EQ(h.size(), 2u) << "axis ray through apex";
        EXPECT_NEAR(h[0].position.GetZ(), -5.0, 1e-6);
        EXPECT_TRUE(h[0].entering) << "base-cap hit must be entering";
        EXPECT_NEAR(h[1].position.GetZ(),  5.0, 1e-6);
        EXPECT_FALSE(h[1].entering) << "apex hit (ray exiting upward) must be exiting";
    }

    // Down the z-axis from above: enter AT the apex (z=+5), exit base cap.
    // This is the case the parity-only test could not see: the apex hit
    // must be entering=true here.
    {
        auto h = cone.Intersections(Vector3D(0, 0, 20), Vector3D(0, 0, -1));
        ASSERT_EQ(h.size(), 2u) << "reverse axis ray through apex";
        EXPECT_NEAR(h[0].position.GetZ(), 5.0, 1e-6);
        EXPECT_TRUE(h[0].entering) << "apex hit (ray entering downward) must be entering";
        EXPECT_NEAR(h[1].position.GetZ(), -5.0, 1e-6);
        EXPECT_FALSE(h[1].entering) << "base-cap exit must be exiting";
    }

    // Transverse ray exactly through the apex point: the tip is a single
    // point (radius 0), so a perpendicular ray is a measure-zero grazing
    // touch -- exactly zero crossings, not merely an even count.
    {
        auto h = cone.Intersections(Vector3D(-10, 0, 5), Vector3D(1, 0, 0));
        EXPECT_EQ(h.size(), 0u) << "transverse apex ray must produce no crossings";
    }

    // Diagonal ray through the apex: enters the body, exits at the apex.
    {
        Vector3D dir(0.5, 0, 1);
        dir.normalize();
        auto h = cone.Intersections(Vector3D(-5, 0, -5), dir);
        ASSERT_EQ(h.size(), 2u) << "diagonal ray through apex";
        EXPECT_NE(h[0].entering, h[1].entering) << "one entering, one exiting";
        EXPECT_NEAR(h[1].position.GetZ(), 5.0, 1e-6) << "second hit is the apex";
        EXPECT_FALSE(h[1].entering) << "apex exit (ray leaving upward)";
    }

    // Reverse-pointed cone: apex at z=-5, base radius 5 at z=+5 (body
    // above apex). Up the z-axis: enter AT the apex, exit base cap.
    Cone cone_rev(0, 0, 0, 5, 10);
    {
        auto h = cone_rev.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        ASSERT_EQ(h.size(), 2u) << "reverse cone axis ray through apex";
        EXPECT_NEAR(h[0].position.GetZ(), -5.0, 1e-6);
        EXPECT_TRUE(h[0].entering) << "bottom-apex hit must be entering";
        EXPECT_NEAR(h[1].position.GetZ(), 5.0, 1e-6);
        EXPECT_FALSE(h[1].entering) << "top base-cap must be exiting";
    }
}

// =========================================================================
// GenericPolycone vertex-aimed rays: rays aimed directly at polygon
// vertices where two conical/disk surface segments meet.
// =========================================================================
TEST(ShapeInvariants, GenericPolyconeVertexRays) {
    // Solid triangle: vertices at (0,-5), (5,0), (0,5)
    GenericPolycone gpc({0, 5, 0}, {-5, 0, 5});

    // Ray aimed at the widest vertex (5, 0) from outside
    {
        auto hits = gpc.Intersections(Vector3D(10, 0, 0), Vector3D(-1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone vertex-aimed (widest): odd count " << hits.size();
    }

    // Ray aimed at an r=0 tip vertex from outside
    {
        auto hits = gpc.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone vertex-aimed (r=0 tip): odd count " << hits.size();
    }

    // Hollow hexagon: vertices at z=-5, 0, 5 with r changing
    GenericPolycone gpc_hex({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5});

    // Ray aimed at the junction between outer and inner surface at z=5
    {
        auto hits = gpc_hex.Intersections(Vector3D(10, 0, 5), Vector3D(-1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone hex vertex at z=5: odd count " << hits.size();
    }

    // Ray aimed at the junction between outer and inner at z=-5 (bottom closure)
    {
        auto hits = gpc_hex.Intersections(Vector3D(10, 0, -5), Vector3D(-1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone hex vertex at z=-5: odd count " << hits.size();
    }
}

// =========================================================================
// GenericPolycone z-axis rays through r=0 vertices: the conical surface
// degenerates to a point on the z-axis.
// =========================================================================
TEST(ShapeInvariants, GenericPolyconeZAxisRays) {
    // Solid triangle with r=0 tips at z=-5 and z=5
    GenericPolycone gpc({0, 5, 0}, {-5, 0, 5});

    // Z-axis ray through both tips
    {
        auto hits = gpc.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone z-axis through tips: odd count " << hits.size();
    }

    // Crystal profile with r=0 bottom section
    GenericPolycone crystal({5,38,40,40,30,30,15,15,0,0,5}, {0,0,70,80,80,77,77,80,80,60,60});
    {
        auto hits = crystal.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Crystal z-axis ray: odd count " << hits.size();
    }

    // Near-axis ray (epsilon off z-axis)
    {
        auto hits = crystal.Intersections(Vector3D(1e-8, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Crystal near-axis ray: odd count " << hits.size();
    }
}

// =========================================================================
// GenericPolycone tangent rays
// =========================================================================
TEST(ShapeInvariants, TangentRaysGenericPolycone) {
    double eps = 1e-6;

    // Hollow hexagon with max outer radius 5 at z=5
    GenericPolycone gpc({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5});

    // Tangent at y=5 (max outer radius)
    {
        auto hits = gpc.Intersections(Vector3D(0, 5, 3), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone tangent: odd count " << hits.size();
    }

    // Just outside
    {
        auto hits = gpc.Intersections(Vector3D(0, 5 + eps, 3), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 0u) << "GenericPolycone outside tangent should miss";
    }

    // Just inside
    {
        auto hits = gpc.Intersections(Vector3D(0, 5 - eps, 3), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone inside tangent: odd count " << hits.size();
    }
}

// =========================================================================
// GenericPolycone horizontal rays at polygon vertex z-values
// =========================================================================
TEST(ShapeInvariants, GenericPolyconeHorizontalBoundaryRays) {
    // Crystal profile with many z-transitions
    GenericPolycone crystal({5,38,40,40,30,30,15,15,0,0,5}, {0,0,70,80,80,77,77,80,80,60,60});

    // Horizontal rays at each unique z-value in the profile
    double z_vals[] = {0, 60, 70, 77, 80};
    for(double z : z_vals) {
        Vector3D pos(20, 0, z);
        Vector3D dir(-1, 0, 0);
        auto hits = crystal.Intersections(pos, dir);
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Crystal horizontal ray at z=" << z << ": odd count " << hits.size();
    }

    // Also test containment at vertex z-values with z-direction ray
    EXPECT_TRUE(crystal.IsInside(Vector3D(20, 0, 40), Vector3D(0, 0, 1)))
        << "Crystal interior at z=40 should be inside";
    EXPECT_FALSE(crystal.IsInside(Vector3D(0, 0, 40), Vector3D(0, 0, 1)))
        << "Crystal z-axis at z=40 should be outside (inside borehole)";
}

// =========================================================================
// Far-field rays: numerical stability at large distances.
// Tests that quadratic discriminant computation doesn't lose precision.
// =========================================================================
TEST(ShapeInvariants, FarFieldRays) {
    double distances[] = {1e2, 1e4, 1e6};
    int wrong = 0;
    int total = 0;

    // --- Cone ---
    Cone cone(0, 5, 0, 3, 10);
    for(double D : distances) {
        auto hits = cone.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
        // Should hit (ray through center)
        if(hits.size() == 0) wrong++;
        total++;
    }

    // --- Polycone ---
    Polycone pc({-5, 0, 5}, {1, 2, 1}, {3, 5, 3});
    for(double D : distances) {
        auto hits = pc.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
    }

    // --- GenericPolycone ---
    GenericPolycone gpc({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5});
    for(double D : distances) {
        auto hits = gpc.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
    }

    // --- Cylinder ---
    Cylinder cyl(5, 0, 10);
    for(double D : distances) {
        auto hits = cyl.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
    }

    // --- Sphere ---
    Sphere sphere(5, 0);
    for(double D : distances) {
        auto hits = sphere.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
    }

    // --- EllipticalTube ---
    EllipticalTube etube(3, 5, 10);
    for(double D : distances) {
        auto hits = etube.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
    }

    // --- Box ---
    Box box(10, 8, 6);
    for(double D : distances) {
        auto hits = box.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        if(hits.size() % 2 != 0) wrong++;
        total++;
    }

    EXPECT_EQ(wrong, 0)
        << "Far-field rays: " << wrong << "/" << total << " incorrect";
}

// =========================================================================
// EllipticalTube tangent rays at both semi-axes
// =========================================================================
TEST(ShapeInvariants, TangentRaysEllipticalTube) {
    double eps = 1e-6;
    EllipticalTube etube(3, 5, 10);

    // Tangent at y=5 (major semi-axis), ray in x-direction
    {
        auto hits = etube.Intersections(Vector3D(0, 5, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "EllipticalTube tangent y=5: odd count";
    }
    {
        auto hits = etube.Intersections(Vector3D(0, 5 + eps, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 0u) << "EllipticalTube outside y-tangent should miss";
    }
    {
        auto hits = etube.Intersections(Vector3D(0, 5 - eps, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u) << "EllipticalTube inside y-tangent should hit twice";
    }

    // Tangent at x=3 (minor semi-axis), ray in y-direction
    {
        auto hits = etube.Intersections(Vector3D(3, 0, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "EllipticalTube tangent x=3: odd count";
    }
    {
        auto hits = etube.Intersections(Vector3D(3 + eps, 0, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size(), 0u) << "EllipticalTube outside x-tangent should miss";
    }
    {
        auto hits = etube.Intersections(Vector3D(3 - eps, 0, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size(), 2u) << "EllipticalTube inside x-tangent should hit twice";
    }
}

// =========================================================================
// Trd edge-skimming and vertex-aimed rays
// =========================================================================
TEST(ShapeInvariants, TrdEdgeAndVertexRays) {
    // Trd: half-widths dx1=5,dx2=3,dy1=4,dy2=2, dz=6
    // Centered at origin. Vertices at (+-dx, +-dy, +-dz)
    Trd trd(5, 3, 4, 2, 6);

    // Ray aimed at a vertex: bottom corner (+5, +4, -6)
    {
        Vector3D dir = Vector3D(5, 4, -6);
        dir.normalize();
        auto hits = trd.Intersections(Vector3D(5*2, 4*2, -6*2), Vector3D(-dir.GetX(), -dir.GetY(), -dir.GetZ()));
        EXPECT_EQ(hits.size() % 2, 0u) << "Trd vertex-aimed ray: odd count " << hits.size();
    }

    // Ray skimming along a bottom edge (y=4 face at z=-6, from x=-10 to x=+10)
    {
        auto hits = trd.Intersections(Vector3D(-10, 4, -6), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "Trd edge-skimming ray: odd count " << hits.size();
    }

    // Ray along the z-axis through center
    {
        auto hits = trd.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "Trd z-axis should hit twice";
    }
}

// =========================================================================
// ExtrPoly edge and face-parallel rays
// =========================================================================
TEST(ShapeInvariants, ExtrPolyEdgeRays) {
    // Square extrusion: polygon (+-3, +-3), z from -4 to 4
    std::vector<std::vector<double>> polygon = {{-3,-3},{3,-3},{3,3},{-3,3}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off, 1.0),
        ExtrPoly::ZSection(4, off, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);

    // Ray along a lateral edge (x=3, y=3, varying z)
    {
        auto hits = ep.Intersections(Vector3D(3, 3, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u) << "ExtrPoly lateral edge ray: odd count " << hits.size();
    }

    // Ray along a lateral face (x=3, y=0, varying z) -- face-grazing
    {
        auto hits = ep.Intersections(Vector3D(3, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u) << "ExtrPoly face-grazing ray: odd count " << hits.size();
    }

    // Ray aimed at a bottom polygon vertex (3, 3, -4) from outside
    {
        Vector3D dir(3, 3, -4);
        dir.normalize();
        auto hits = ep.Intersections(Vector3D(6, 6, -8), Vector3D(-dir.GetX(), -dir.GetY(), -dir.GetZ()));
        EXPECT_EQ(hits.size() % 2, 0u) << "ExtrPoly vertex-aimed ray: odd count " << hits.size();
    }

    // Ray through center
    {
        auto hits = ep.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "ExtrPoly center z-axis should hit twice";
    }
}

// =========================================================================
// CutTube rays parallel to and near cut planes
// =========================================================================
TEST(ShapeInvariants, CutTubeCutPlaneRays) {
    // CutTube with tilted cuts: 30-degree tilt on both ends
    double tilt = 30.0 * M_PI / 180.0;
    Vector3D low_norm(std::sin(tilt), 0, -std::cos(tilt));
    Vector3D high_norm(-std::sin(tilt), 0, std::cos(tilt));
    CutTube ct(2, 5, 10, low_norm, high_norm);

    // Ray parallel to the low cut plane, just inside
    {
        Vector3D dir(std::cos(tilt), 0, std::sin(tilt));
        dir.normalize();
        auto hits = ct.Intersections(Vector3D(-20, 0, -9), dir);
        EXPECT_EQ(hits.size() % 2, 0u) << "CutTube low-plane parallel ray: odd count " << hits.size();
    }

    // Ray parallel to the high cut plane, just inside
    {
        Vector3D dir(std::cos(tilt), 0, std::sin(tilt));
        dir.normalize();
        auto hits = ct.Intersections(Vector3D(-20, 0, 9), dir);
        EXPECT_EQ(hits.size() % 2, 0u) << "CutTube high-plane parallel ray: odd count " << hits.size();
    }

    // Ray through the intersection of cut plane and barrel surface
    {
        auto hits = ct.Intersections(Vector3D(5, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u) << "CutTube barrel/plane edge ray: odd count " << hits.size();
    }

    // Ray through center
    {
        auto hits = ct.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u) << "CutTube z-axis ray: odd count " << hits.size();
    }
}

// =========================================================================
// BooleanGeometry coincident surface rays: ray through the touching
// point of two shapes in a union/intersection/subtraction.
// =========================================================================
TEST(ShapeInvariants, BooleanCoincidentSurfaceRays) {
    // Two spheres touching at origin: centers at (-5,0,0) and (5,0,0), radius 5
    auto s1 = Sphere(Placement(Vector3D(-5, 0, 0)), 5, 0).create();
    auto s2 = Sphere(Placement(Vector3D(5, 0, 0)), 5, 0).create();

    // Union: ray through the touching point along x
    {
        auto bg = BooleanGeometry(BooleanOperation::UNION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "Union touching spheres x-ray: odd count " << hits.size();
    }

    // Intersection: ray through the touching point along x
    {
        auto bg = BooleanGeometry(BooleanOperation::INTERSECTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "Intersection touching spheres x-ray: odd count " << hits.size();
    }

    // Subtraction: ray through the touching point along x
    {
        auto bg = BooleanGeometry(BooleanOperation::SUBTRACTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "Subtraction touching spheres x-ray: odd count " << hits.size();
    }

    // Overlapping spheres: centers at (-2,0,0) and (2,0,0), radius 5
    auto s3 = Sphere(Placement(Vector3D(-2, 0, 0)), 5, 0).create();
    auto s4 = Sphere(Placement(Vector3D(2, 0, 0)), 5, 0).create();

    // Union: ray through the overlap region
    {
        auto bg = BooleanGeometry(BooleanOperation::UNION, s3, s4);
        auto hits = bg.Intersections(Vector3D(0, -10, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "Union overlapping spheres y-ray: odd count " << hits.size();
    }
}

// =========================================================================
// Para edge-skimming and vertex-aimed rays
// =========================================================================
TEST(ShapeInvariants, ParaEdgeRays) {
    // Para: dx=5, dy=4, dz=6, alpha=0.3, theta=0.2, phi=0.1
    Para para(5, 4, 6, 0.3, 0.2, 0.1);

    // Ray through center along z
    {
        auto hits = para.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "Para z-axis should hit twice";
    }

    // Ray along x at the z-extent boundary
    {
        auto hits = para.Intersections(Vector3D(-20, 0, 6), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u) << "Para z-boundary x-ray: odd count " << hits.size();
    }

    // Axis-aligned rays from various directions
    Vector3D dirs[] = {
        Vector3D(1,0,0), Vector3D(0,1,0), Vector3D(0,0,1),
        Vector3D(-1,0,0), Vector3D(0,-1,0), Vector3D(0,0,-1)
    };
    for(auto const & d : dirs) {
        auto hits = para.Intersections(Vector3D(0, 0, 0), d);
        EXPECT_EQ(hits.size() % 2, 0u) << "Para origin ray: odd count " << hits.size();
    }
}

// =========================================================================
// TriangularMesh vertex-aimed and edge-aligned rays.
// Uses a cube mesh (12 triangles, 8 vertices, 18 edges).
// =========================================================================
TEST(ShapeInvariants, TriangularMeshVertexAndEdgeRays) {
    double s = 3;
    Vector3D v[8] = {
        Vector3D(-s,-s,-s), Vector3D(s,-s,-s), Vector3D(s,s,-s), Vector3D(-s,s,-s),
        Vector3D(-s,-s, s), Vector3D(s,-s, s), Vector3D(s,s, s), Vector3D(-s,s, s)
    };
    std::vector<std::array<Vector3D, 3>> triangles;
    auto addQuad = [&](Vector3D a, Vector3D b, Vector3D c, Vector3D d) {
        triangles.push_back({{a, b, c}});
        triangles.push_back({{a, c, d}});
    };
    addQuad(v[3], v[2], v[1], v[0]); // -z
    addQuad(v[4], v[5], v[6], v[7]); // +z
    addQuad(v[0], v[1], v[5], v[4]); // -y
    addQuad(v[1], v[2], v[6], v[5]); // +x
    addQuad(v[2], v[3], v[7], v[6]); // +y
    addQuad(v[3], v[0], v[4], v[7]); // -x
    TriangularMesh mesh(triangles);

    // Ray aimed at a vertex (3,3,3) from outside
    {
        Vector3D dir(1, 1, 1);
        dir.normalize();
        auto hits = mesh.Intersections(Vector3D(10, 10, 10), Vector3D(-dir.GetX(), -dir.GetY(), -dir.GetZ()));
        EXPECT_EQ(hits.size() % 2, 0u) << "Mesh vertex-aimed ray: odd count " << hits.size();
    }

    // Ray along an edge (x=3, y=3, varying z)
    {
        auto hits = mesh.Intersections(Vector3D(s, s, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u) << "Mesh edge-aligned ray: odd count " << hits.size();
    }

    // Ray along a face edge (x=3, y=0, z=-3 -> z=3 diagonal of +x face)
    {
        auto hits = mesh.Intersections(Vector3D(s, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u) << "Mesh face-diagonal ray: odd count " << hits.size();
    }

    // Ray through center
    {
        auto hits = mesh.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "Mesh center z-axis should hit twice";
    }
}

// =========================================================================
// Ellipsoid vertex/axis rays: rays along each semi-axis and at cut boundaries
// =========================================================================
TEST(ShapeInvariants, EllipsoidAxisRays) {
    Ellipsoid ell(3, 5, 7, -4, 6);

    // Ray along each axis through center
    {
        auto hits = ell.Intersections(Vector3D(-10, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u) << "Ellipsoid x-axis should hit twice";
    }
    {
        auto hits = ell.Intersections(Vector3D(0, -10, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size(), 2u) << "Ellipsoid y-axis should hit twice";
    }
    {
        auto hits = ell.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "Ellipsoid z-axis should hit twice";
    }

    // Ray exactly at z-cut boundary
    {
        auto hits = ell.Intersections(Vector3D(-10, 0, 6), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Ellipsoid at upper z-cut boundary: odd count " << hits.size();
    }
    {
        auto hits = ell.Intersections(Vector3D(-10, 0, -4), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Ellipsoid at lower z-cut boundary: odd count " << hits.size();
    }
}
