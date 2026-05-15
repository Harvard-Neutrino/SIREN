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
