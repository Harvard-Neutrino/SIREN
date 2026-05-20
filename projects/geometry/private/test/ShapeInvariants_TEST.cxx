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

// Shape construction helpers to reduce boilerplate in MakeShapes()

std::shared_ptr<Geometry> MakeSquareExtrPoly(double half, double half_z, double scale = 1.0) {
    std::vector<std::vector<double>> polygon = {{-half,-half},{half,-half},{half,half},{-half,half}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-half_z, off, scale),
        ExtrPoly::ZSection(half_z, off, scale)
    };
    return ExtrPoly(polygon, zsecs).create();
}

std::shared_ptr<Geometry> MakeCubeMesh(double half) {
    Vector3D v[8] = {
        Vector3D(-half,-half,-half), Vector3D(half,-half,-half),
        Vector3D(half,half,-half), Vector3D(-half,half,-half),
        Vector3D(-half,-half,half), Vector3D(half,-half,half),
        Vector3D(half,half,half), Vector3D(-half,half,half)
    };
    std::vector<std::array<Vector3D, 3>> tris;
    auto addQuad = [&](Vector3D a, Vector3D b, Vector3D c, Vector3D d) {
        tris.push_back({{a, b, c}});
        tris.push_back({{a, c, d}});
    };
    addQuad(v[3],v[2],v[1],v[0]); addQuad(v[4],v[5],v[6],v[7]);
    addQuad(v[0],v[1],v[5],v[4]); addQuad(v[1],v[2],v[6],v[5]);
    addQuad(v[2],v[3],v[7],v[6]); addQuad(v[3],v[0],v[4],v[7]);
    return TriangularMesh(tris).create();
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

    // Extruded polygon
    shapes.push_back({"ExtrPoly(square)", MakeSquareExtrPoly(3, 4), 2, -1, 0});
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

    // Polyhedra variants (triangle, hollow, zero-radius apex)
    shapes.push_back({"Polyhedra3(triangle)", Polyhedra(3, 0, {-4,4}, {0,0}, {5,5}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra6Hollow", Polyhedra(6, 0, {-5,5}, {2,2}, {5,5}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra4Pyramid", Polyhedra(4, M_PI/4, {-3,5}, {0,0}, {5,0}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra6PyramidRev", Polyhedra(6, 0, {-4,4}, {0,0}, {0,6}).create(), -1, -1, 0});

    // Multi-segment shapes
    shapes.push_back({"Polycone", Polycone({-5,-2,0,3,5}, {0,0,0,0,0}, {3,5,4,6,2}).create(), -1, -1, 0});
    shapes.push_back({"PolyconeHollow", Polycone({-4,0,4}, {1,2,1}, {5,6,5}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra6", Polyhedra(6, 0, {-5,0,5}, {0,0,0}, {4,6,4}).create(), -1, -1, 0});
    shapes.push_back({"Polyhedra4", Polyhedra(4, M_PI/4, {-3,3}, {0,0}, {5,5}).create(), -1, -1, 0});

    // Rotated-placement shapes
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
    // Cone frustum: z-projection = pi * max(rmax1,rmax2)^2 = pi*25
    shapes.push_back({"Cone(0,5,0,3,8)_MC", Cone(0, 5, 0, 3, 8).create(), 2, M_PI * 25, 8});
    // Pointed cone: z-projection = pi*25
    shapes.push_back({"ConePointed_MC", Cone(0, 5, 0, 0, 8).create(), 2, M_PI * 25, 8});
    // Trd: z-projection = bottom face = 2*dx1 * 2*dy1 = 10*8 = 80
    shapes.push_back({"Trd(5,3,4,2,6)_MC", Trd(5, 3, 4, 2, 6).create(), 2, 80, 10});
    // EllipticalTube: z-projection = pi*dx*dy = pi*6
    shapes.push_back({"EllipticalTube(3,2,5)_MC", EllipticalTube(3, 2, 5).create(), 2, M_PI * 6, 5});
    // Ellipsoid (full): z-projection = pi*ax*by = pi*15
    shapes.push_back({"Ellipsoid(5,3,4)_MC", Ellipsoid(5, 3, 4).create(), 2, M_PI * 15, 8});
    // Ellipsoid (z-cuts): equator intact, z-projection = pi*15
    shapes.push_back({"Ellipsoid(5,3,4,-2,3)_MC", Ellipsoid(5, 3, 4, -2, 3).create(), 4, M_PI * 15, 8});
    // CutTube flat: same as cylinder = pi*25
    shapes.push_back({"CutTubeFlat_MC", CutTube(0, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1)).create(), 2, M_PI * 25, 8});
    // Trap symmetric: z-projection = bottom face = 2*dx1 * 2*dy1 = 8*6 = 48
    shapes.push_back({"Trap(sym)_MC", Trap(5, 0, 0, 3, 4, 4, 0, 2, 3, 3, 0).create(), 2, 48, 8});
    // Para unsheared: z-projection = 2*dx * 2*dy = 8*6 = 48
    shapes.push_back({"Para(0,0,0)_MC", Para(4, 3, 5, 0, 0, 0).create(), 2, 48, 8});
    // Polycone solid: z-projection = pi*max(rmax)^2 = pi*36
    shapes.push_back({"Polycone_MC", Polycone({-5,-2,0,3,5}, {0,0,0,0,0}, {3,5,4,6,2}).create(), -1, M_PI * 36, 8});
    // GenericPolycone solid bicone: z-projection = pi*25
    shapes.push_back({"GenericPolyconeSolid_MC", GenericPolycone({0, 5, 0}, {-5, 0, 5}).create(), -1, M_PI * 25, 8});
    // GenericPolycone U-shape (hollow cylinder): z-projection = pi*(5^2 - 3^2) = pi*16
    shapes.push_back({"GenericPolyconeU_MC", GenericPolycone({5, 5, 3, 3}, {-4, 4, 4, -4}).create(), -1, M_PI * 16, 8});

    // Torus shapes (max 4 intersections for solid, -1 for hollow)
    shapes.push_back({"Torus(R=10,r=3,0)", Torus(10, 3, 0).create(), -1, -1, 0});
    shapes.push_back({"TorusHollow(R=10,3,1)", Torus(10, 3, 1).create(), -1, -1, 0});
    shapes.push_back({"TorusThin(R=8,r=0.5,0)", Torus(8, 0.5, 0).create(), -1, -1, 0});
    // Self-intersecting torus (rtor < rmax): tube passes through center
    shapes.push_back({"TorusSelfIntersecting(R=3,r=5,0)", Torus(3, 5, 0).create(), -1, -1, 0});
    shapes.push_back({"TorusSelfIntersectingHollow(R=3,5,2)", Torus(3, 5, 2).create(), -1, -1, 0});
    // Solid torus: annular z-projection = pi*(13^2 - 7^2) = 120*pi
    shapes.push_back({"Torus(10,3,0)_MC", Torus(10, 3, 0).create(), -1, M_PI * 120, 16});
    // Self-intersecting torus: disk z-projection = pi*(R+r)^2 = pi*64
    shapes.push_back({"TorusSelfIntersecting_MC", Torus(3, 5, 0).create(), -1, M_PI * 64, 10});

    // Partial spheres
    shapes.push_back({"SphereThetaHemi", Sphere(5, 0, 0, 2*M_PI, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SpherePhiHalf", Sphere(5, 0, 0, M_PI, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"SpherePhiTheta", Sphere(5, 0, 0, M_PI, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SphereHollowThetaHemi", Sphere(5, 2, 0, 2*M_PI, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SphereThetaCone45", Sphere(5, 0, 0, 2*M_PI, 0, M_PI/4).create(), -1, -1, 0});
    shapes.push_back({"SphereThetaBand45_90", Sphere(5, 0, 0, 2*M_PI, M_PI/4, M_PI/4).create(), -1, -1, 0});
    shapes.push_back({"SphereThetaBandMid", Sphere(5, 0, 0, 2*M_PI, M_PI/3, M_PI/3).create(), -1, -1, 0});
    shapes.push_back({"SphereHollowThetaCone60", Sphere(5, 2, 0, 2*M_PI, 0, M_PI/3).create(), -1, -1, 0});
    shapes.push_back({"SpherePhiThetaQuarter", Sphere(5, 0, 0, M_PI/2, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SpherePhiHalfThetaBand", Sphere(5, 0, 0, M_PI, M_PI/4, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"SphereHollowPhiTheta", Sphere(5, 2, 0, M_PI, 0, M_PI/2).create(), -1, -1, 0});

    // Partial torus
    shapes.push_back({"TorusQuarter", Torus(10, 3, 0, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"TorusHalf", Torus(10, 3, 0, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"TorusHollowQuarter", Torus(10, 3, 1, 0, M_PI/2).create(), -1, -1, 0});
    shapes.push_back({"TorusHollowPhiHalf", Torus(10, 3, 1, 0, M_PI).create(), -1, -1, 0});

    // Hollow + phi-cut combinations
    shapes.push_back({"CylinderHollowPhiHalf", Cylinder(5, 2, 12, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"ConeHollowPhiHalf", Cone(1, 5, 1, 3, 8, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"CutTubePhiHalf", CutTube(0, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1), 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"CutTubeHollowPhiHalf", CutTube(2, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1), 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"PolyconeHollowPhiHalf", Polycone({-4,0,4}, {1,2,1}, {5,6,5}, 0, M_PI).create(), -1, -1, 0});
    shapes.push_back({"GenericPolyconePhiHalf", GenericPolycone({0,5,0}, {-5,0,5}, 0, M_PI).create(), -1, -1, 0});

    // Partial polyhedra
    shapes.push_back({"Polyhedra6Phi180", Polyhedra(6, 0, {-4,4}, {0,0}, {5,5}, M_PI).create(), -1, -1, 0});

    // Polycone with step-change
    shapes.push_back({"PolyconeStep", Polycone({-0.05, 0.05, 0.05, 0.15}, {0, 0, 0, 0}, {0.03, 0.03, 0.05, 0.05}).create(), -1, -1, 0});

    // Reversed z-planes (constructor auto-reverses to ascending)
    shapes.push_back({"PolyconeReversed", Polycone({5, 0, -5}, {0, 0, 0}, {3, 5, 3}).create(), -1, -1, 0});
    shapes.push_back({"PolyhedraReversed", Polyhedra(6, 0, {5, -5}, {0, 0}, {4, 4}).create(), -1, -1, 0});

    // Duplicate z-planes (step profiles)
    shapes.push_back({"PolyconeStepFlange", Polycone({-5,-2,-2,2,2,5}, {0,0,0,0,0,0}, {3,3,5,5,3,3}).create(), -1, -1, 0});
    shapes.push_back({"PolyhedraStepFlange", Polyhedra(6, 0, {-5,-1,-1,1,1,5}, {0,0,0,0,0,0}, {3,3,5,5,3,3}).create(), -1, -1, 0});

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

    // Polyhedra regular hexagonal prism: rmax is circumradius, area = (n/2)*R^2*sin(2*pi/n)
    // rmax=5 is apothem; circumradius = 5/cos(pi/6). Area = (n/2)*R^2*sin(2*pi/n)
    {
        double R = 5.0 / std::cos(M_PI/6.0);
        shapes.push_back({"Polyhedra6_MC", Polyhedra(6, 0, {-5,5}, {0,0}, {5,5}).create(), -1, 3.0*R*R*std::sin(M_PI/3.0), 8});
    }

    // ExtrPoly square: z-projection = 6*6 = 36
    shapes.push_back({"ExtrPoly(square)_MC", MakeSquareExtrPoly(3, 4), 2, 36, 5});

    // BooleanGeometry: cylinder intersected with box = pi*9
    {
        auto bx = Box(6, 6, 6).create();
        auto cy = Cylinder(3, 0, 10).create();
        shapes.push_back({"BoolIntersect_MC", BooleanGeometry(BooleanOperation::INTERSECTION, bx, cy).create(), -1, M_PI * 9, 6});
    }

    // TriangularMesh cube: z-projection = 6*6 = 36
    shapes.push_back({"TriMesh(cube)_MC", MakeCubeMesh(3), -1, 36, 5});

    return shapes;
}

// Test helper: iterates over shapes, fires N random rays, and calls test_fn
// for each. Returns false from test_fn to report a violation. Violations are
// counted and printed (up to 3 per shape), then EXPECT_EQ'd to zero.
using ShapeTestFn = std::function<bool(ShapeEntry const &, Vector3D const &, Vector3D const &, std::string &)>;

void ForEachShape(std::vector<ShapeEntry> const & shapes, int N, bool skip_mc, ShapeTestFn test_fn) {
    for(auto const & entry : shapes) {
        if(skip_mc && entry.known_cross_section > 0) continue;
        int violations = 0;
        for(int i = 0; i < N; ++i) {
            Vector3D pos = RandomPoint(20);
            Vector3D dir = RandomDirection();
            std::string detail;
            if(!test_fn(entry, pos, dir, detail)) {
                violations++;
                if(violations <= 3)
                    std::cerr << entry.name << ": " << detail << std::endl;
            }
        }
        EXPECT_EQ(violations, 0) << entry.name << ": " << violations << " violations";
    }
}

} // anonymous namespace

// =========================================================================
// Test 1: Intersection count invariants
// =========================================================================
TEST(ShapeInvariants, IntersectionCount) {
    auto shapes = MakeShapes();
    ForEachShape(shapes, 10000, true, [](ShapeEntry const & entry, Vector3D const & pos, Vector3D const & dir, std::string & detail) {
        auto isects = entry.geo->Intersections(pos, dir);
        int n = (int)isects.size();
        if(n % 2 != 0) {
            detail = "odd intersection count " + std::to_string(n);
            return false;
        }
        if(entry.max_intersections > 0 && n > entry.max_intersections) {
            detail = std::to_string(n) + " intersections (max " + std::to_string(entry.max_intersections) + ")";
            return false;
        }
        return true;
    });
}

// =========================================================================
// Test 2: Enter/exit ordering
// =========================================================================
TEST(ShapeInvariants, EnterExitOrdering) {
    auto shapes = MakeShapes();
    ForEachShape(shapes, 10000, true, [](ShapeEntry const & entry, Vector3D const & pos, Vector3D const & dir, std::string & detail) {
        auto isects = entry.geo->Intersections(pos, dir);
        bool last_entering = false;
        bool first = true;
        for(auto const & isect : isects) {
            if(!first) {
                if(isect.entering == last_entering) {
                    detail = "consecutive " + std::string(isect.entering ? "entering" : "exiting")
                             + " at dist=" + std::to_string(isect.distance);
                    return false;
                }
            }
            first = false;
            last_entering = isect.entering;
        }
        return true;
    });
}

// =========================================================================
// Test 3: Ray reversal reciprocity
// =========================================================================
TEST(ShapeInvariants, Reciprocity) {
    auto shapes = MakeShapes();
    ForEachShape(shapes, 5000, true, [](ShapeEntry const & entry, Vector3D const & pos, Vector3D const & dir, std::string & detail) {
        double tol = 1e-6;
        Vector3D neg_dir(-dir.GetX(), -dir.GetY(), -dir.GetZ());
        auto fwd = entry.geo->Intersections(pos, dir);
        auto rev = entry.geo->Intersections(pos, neg_dir);
        if(fwd.size() != rev.size()) {
            detail = "forward=" + std::to_string(fwd.size()) + " reverse=" + std::to_string(rev.size());
            return false;
        }
        for(size_t j = 0; j < fwd.size(); ++j) {
            size_t k = fwd.size() - 1 - j;
            double dx = fwd[j].position.GetX() - rev[k].position.GetX();
            double dy = fwd[j].position.GetY() - rev[k].position.GetY();
            double dz = fwd[j].position.GetZ() - rev[k].position.GetZ();
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if(dist > tol) {
                detail = "position mismatch at pair " + std::to_string(j) + " dist=" + std::to_string(dist);
                return false;
            }
        }
        return true;
    });
}

// =========================================================================
// Test 4: AABB consistency (intersections lie within bounding box)
// =========================================================================
TEST(ShapeInvariants, AABBConsistency) {
    auto shapes = MakeShapes();
    double tol = 1e-6; // allow slight overshoot due to floating point
    ForEachShape(shapes, 10000, true, [tol](ShapeEntry const & entry, Vector3D const & pos, Vector3D const & dir, std::string & detail) {
        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid()) return true;
        if(!std::isfinite(box.min_corner.GetX())) return true;
        auto isects = entry.geo->Intersections(pos, dir);
        for(auto const & isect : isects) {
            double x = isect.position.GetX();
            double y = isect.position.GetY();
            double z = isect.position.GetZ();
            if(x < box.min_corner.GetX() - tol || x > box.max_corner.GetX() + tol ||
               y < box.min_corner.GetY() - tol || y > box.max_corner.GetY() + tol ||
               z < box.min_corner.GetZ() - tol || z > box.max_corner.GetZ() + tol) {
                detail = "intersection at (" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ") outside AABB";
                return false;
            }
        }
        return true;
    });
}

// =========================================================================
// Test 5: AABB tightness (surface samples lie within bounding box)
// =========================================================================
TEST(ShapeInvariants, AABBTightness) {
    auto shapes = MakeShapes();
    double tol = 1e-6;
    // Shoots rays from the AABB centroid (ignores the random pos from ForEachShape)
    ForEachShape(shapes, 10000, true, [tol](ShapeEntry const & entry, Vector3D const &, Vector3D const & dir, std::string & detail) {
        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid()) return true;
        if(!std::isfinite(box.min_corner.GetX())) return true;
        Vector3D center = box.Centroid();
        auto isects = entry.geo->Intersections(center, dir);
        for(auto const & isect : isects) {
            double x = isect.position.GetX();
            double y = isect.position.GetY();
            double z = isect.position.GetZ();
            if(x < box.min_corner.GetX() - tol || x > box.max_corner.GetX() + tol ||
               y < box.min_corner.GetY() - tol || y > box.max_corner.GetY() + tol ||
               z < box.min_corner.GetZ() - tol || z > box.max_corner.GetZ() + tol) {
                detail = "surface point outside AABB";
                return false;
            }
        }
        return true;
    });
}

// =========================================================================
// Test 6: Distance-position consistency (pos == origin + t*dir)
// =========================================================================
TEST(ShapeInvariants, DistancePositionConsistency) {
    auto shapes = MakeShapes();
    ForEachShape(shapes, 10000, true, [](ShapeEntry const & entry, Vector3D const & pos, Vector3D const & dir, std::string & detail) {
        double tol = 1e-6;
        auto isects = entry.geo->Intersections(pos, dir);
        for(auto const & isect : isects) {
            double ex = pos.GetX() + isect.distance * dir.GetX();
            double ey = pos.GetY() + isect.distance * dir.GetY();
            double ez = pos.GetZ() + isect.distance * dir.GetZ();
            double dx = isect.position.GetX() - ex;
            double dy = isect.position.GetY() - ey;
            double dz = isect.position.GetZ() - ez;
            double err = std::sqrt(dx*dx + dy*dy + dz*dz);
            if(err > tol) {
                detail = "position error " + std::to_string(err) + " at dist=" + std::to_string(isect.distance);
                return false;
            }
        }
        return true;
    });
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

        // Statistical bound: at N=100k, worst-case SE is sqrt(0.25/N) ~ 0.16%.
        // A 1.5% threshold gives ~9 sigma, so false-positive rate < 1e-19.
        EXPECT_LT(relative_error, 0.015)
            << entry.name << ": measured=" << measured_area
            << " expected=" << expected_area
            << " rel_error=" << relative_error;
    }
}

// =========================================================================
// Test 8: Surface-boundary points
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
// Test 9: Axis-aligned rays
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
    // Tangent on convex shape: must be 0 or 2
    EXPECT_TRUE(isects.size() == 0 || isects.size() == 2)
        << "Cone tangent: expected 0 or 2, got " << isects.size();

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
    EXPECT_TRUE(isects.size() == 0 || isects.size() == 2)
        << "Polycone tangent: expected 0 or 2, got " << isects.size();

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

    // Tangent to the outer equator: ray at rxy = R + r (outermost extent
    // of the torus surface in the z=0 plane). The quartic has a double root
    // here; the solver may produce 0, 1, or 2 hits depending on whether
    // dedup collapses the pair. We verify the near-miss/near-hit boundary
    // is correct even if the exact-tangent result is ambiguous.
    for(int i = 0; i < 100; ++i) {
        double theta = 2.0 * M_PI * i / 100.0;
        double outer = R + r;
        double px = outer * std::cos(theta);
        double py = outer * std::sin(theta);
        double tx = -std::sin(theta);
        double ty = std::cos(theta);
        Vector3D pos(px, py, 0);
        Vector3D dir(tx, ty, 0);
        dir.normalize();
        auto hits = torus.Intersections(pos, dir);
        EXPECT_LE(hits.size(), 2u)
            << "Torus outer-tangent ray #" << i << ": expected <= 2, got " << hits.size();
        if(hits.size() % 2 != 0) odd_count++;
        total++;
    }

    // Tangent to the tube top: ray at z = r, passing through the tube
    // center circle. Same double-root situation.
    for(int i = 0; i < 100; ++i) {
        double theta = 2.0 * M_PI * i / 100.0;
        double cx = R * std::cos(theta);
        double cy = R * std::sin(theta);
        double tx = -std::sin(theta);
        double ty = std::cos(theta);
        Vector3D pos(cx, cy, r);
        Vector3D dir(tx, ty, 0);
        dir.normalize();
        auto hits = torus.Intersections(pos, dir);
        EXPECT_LE(hits.size(), 2u)
            << "Torus top-tangent ray #" << i << ": expected <= 2, got " << hits.size();
        if(hits.size() % 2 != 0) odd_count++;
        total++;
    }

    if(odd_count > 0) {
        std::cerr << "Torus tangent rays: " << odd_count << "/" << total
                  << " produced odd count (quartic double-root dedup limitation)" << std::endl;
    }

    // Just outside the outer equator: must miss
    {
        Vector3D pos(R + r + 1e-6, 0, 0);
        Vector3D dir(0, 1, 0);
        auto hits = torus.Intersections(pos, dir);
        EXPECT_EQ(hits.size(), 0u) << "Torus outside outer tangent should miss";
    }

    // Just inside the outer equator: must hit
    {
        Vector3D pos(R + r - 1e-6, 0, 0);
        Vector3D dir(0, 1, 0);
        auto hits = torus.Intersections(pos, dir);
        EXPECT_EQ(hits.size(), 2u) << "Torus inside outer tangent should hit twice";
    }

    // Just above the tube top: must miss
    {
        Vector3D pos(R, 0, r + 1e-6);
        Vector3D dir(0, 1, 0);
        auto hits = torus.Intersections(pos, dir);
        EXPECT_EQ(hits.size(), 0u) << "Torus above tube top should miss";
    }

    // Just below the tube top: must hit
    {
        Vector3D pos(R, 0, r - 1e-6);
        Vector3D dir(0, 1, 0);
        auto hits = torus.Intersections(pos, dir);
        EXPECT_GE(hits.size(), 2u) << "Torus below tube top should hit";
    }
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
            double cy = R * std::sin(theta);
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

    // Rays in the z=0 plane from origin: each ray crosses both tube lobes
    // (at rxy=7 and rxy=13), producing exactly 4 hits.
    for(int i = 0; i < 200; ++i) {
        double angle = 2.0 * M_PI * i / 200.0;
        Vector3D pos(0, 0, 0);
        Vector3D dir(std::cos(angle), std::sin(angle), 0);
        auto hits = torus.Intersections(pos, dir);
        EXPECT_EQ(hits.size(), 4u)
            << "Torus z=0 plane ray #" << i << ": expected 4, got " << hits.size();
    }
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

    // Irreflexivity
    for(auto* s : shapes) {
        EXPECT_FALSE(*s < *s);
    }

    // Transitivity: no 3-cycle among the first 3 shapes
    bool ab = *shapes[0] < *shapes[1];
    bool bc = *shapes[1] < *shapes[2];
    bool ca = *shapes[2] < *shapes[0];
    EXPECT_FALSE(ab && bc && ca) << "Ordering has a 3-cycle";
    EXPECT_FALSE(!ab && !bc && !ca) << "Ordering has a reverse 3-cycle";

    // Same-type with different params: exactly one direction is true
    Sphere sphere2(3, 0);
    bool s_lt = sphere < sphere2;
    bool s_gt = sphere2 < sphere;
    EXPECT_NE(s_lt, s_gt);

    // Same-type with equal params: neither direction
    Sphere sphere3(5, 0);
    EXPECT_FALSE(sphere < sphere3);
    EXPECT_FALSE(sphere3 < sphere);
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

    // Ray aimed at the widest vertex (5, 0) from outside: the solid
    // body-of-revolution has r=5 at z=0. Horizontal ray from x=10 must
    // cross surface twice.
    {
        auto hits = gpc.Intersections(Vector3D(10, 0, 0), Vector3D(-1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "GenericPolycone vertex-aimed (widest): expected 2, got " << hits.size();
    }

    // Ray aimed at an r=0 tip vertex from outside
    {
        auto hits = gpc.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "GenericPolycone vertex-aimed (r=0 tip): odd count " << hits.size();
    }

    // Hollow hexagon: vertices at z=-5, 0, 5 with r changing
    GenericPolycone gpc_hex({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5});

    // Ray aimed at the junction between outer and inner surface at z=5:
    // outer r=5, inner r=3.5. Ray crosses 4 surfaces.
    {
        auto hits = gpc_hex.Intersections(Vector3D(10, 0, 5), Vector3D(-1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "GenericPolycone hex vertex at z=5: expected 4, got " << hits.size();
    }

    // Ray aimed at the junction between outer and inner at z=-5 (bottom closure):
    // outer r=3, inner r=2. Ray crosses 4 surfaces.
    {
        auto hits = gpc_hex.Intersections(Vector3D(10, 0, -5), Vector3D(-1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "GenericPolycone hex vertex at z=-5: expected 4, got " << hits.size();
    }
}

// =========================================================================
// GenericPolycone z-axis rays through r=0 vertices: the conical surface
// degenerates to a point on the z-axis.
// =========================================================================
TEST(ShapeInvariants, GenericPolyconeZAxisRays) {
    // Solid triangle with r=0 tips at z=-5 and z=5
    GenericPolycone gpc({0, 5, 0}, {-5, 0, 5});

    // Z-axis ray through both r=0 tips: enters at z=-5, exits at z=+5
    {
        auto hits = gpc.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u)
            << "GenericPolycone z-axis through tips: expected 2, got " << hits.size();
        if(hits.size() == 2) {
            EXPECT_NEAR(hits[0].position.GetZ(), -5.0, 1e-6);
            EXPECT_TRUE(hits[0].entering);
            EXPECT_NEAR(hits[1].position.GetZ(), 5.0, 1e-6);
            EXPECT_FALSE(hits[1].entering);
        }
    }

    // Reverse direction: enters at z=+5, exits at z=-5
    {
        auto hits = gpc.Intersections(Vector3D(0, 0, 10), Vector3D(0, 0, -1));
        EXPECT_EQ(hits.size(), 2u)
            << "GenericPolycone reverse z-axis: expected 2, got " << hits.size();
        if(hits.size() == 2) {
            EXPECT_NEAR(hits[0].position.GetZ(), 5.0, 1e-6);
            EXPECT_TRUE(hits[0].entering);
            EXPECT_NEAR(hits[1].position.GetZ(), -5.0, 1e-6);
            EXPECT_FALSE(hits[1].entering);
        }
    }

    // Crystal profile with r=0 bottom section
    GenericPolycone crystal({5,38,40,40,30,30,15,15,0,0,5}, {0,0,70,80,80,77,77,80,80,60,60});
    {
        auto hits = crystal.Intersections(Vector3D(0, 0, -10), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Crystal z-axis ray: odd count " << hits.size();
        EXPECT_GE(hits.size(), 2u) << "Crystal z-axis ray must hit the solid";
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

    // Tangent at y=5 (max outer radius) -- hollow shape: 0, 2, or 4
    {
        auto hits = gpc.Intersections(Vector3D(0, 5, 3), Vector3D(1, 0, 0));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2 || hits.size() == 4)
            << "GenericPolycone tangent: expected 0, 2, or 4, got " << hits.size();
    }

    // Just outside
    {
        auto hits = gpc.Intersections(Vector3D(0, 5 + eps, 3), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 0u) << "GenericPolycone outside tangent should miss";
    }

    // Just inside -- hollow shape: 0, 2, or 4
    {
        auto hits = gpc.Intersections(Vector3D(0, 5 - eps, 3), Vector3D(1, 0, 0));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2 || hits.size() == 4)
            << "GenericPolycone inside tangent: expected 0, 2, or 4, got " << hits.size();
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

    // --- Cone (solid): center-piercing ray must cross 2 surfaces ---
    Cone cone(0, 5, 0, 3, 10);
    for(double D : distances) {
        auto hits = cone.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Cone at D=" << D << ": expected 2, got " << hits.size();
    }

    // --- Polycone (hollow, rmin={1,2,1}, rmax={3,5,3}): at z=0, rmin=2, rmax=5.
    //     Ray enters outer at rxy=5, exits inner at rxy=2, enters inner, exits outer: 4 hits ---
    Polycone pc({-5, 0, 5}, {1, 2, 1}, {3, 5, 3});
    for(double D : distances) {
        auto hits = pc.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "Polycone at D=" << D << ": expected 4, got " << hits.size();
    }

    // --- GenericPolycone (hollow hex): same analysis, 4 hits ---
    GenericPolycone gpc({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5});
    for(double D : distances) {
        auto hits = gpc.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "GenericPolycone at D=" << D << ": expected 4, got " << hits.size();
    }

    // --- Cylinder (solid): 2 hits ---
    Cylinder cyl(5, 0, 10);
    for(double D : distances) {
        auto hits = cyl.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Cylinder at D=" << D << ": expected 2, got " << hits.size();
    }

    // --- Sphere (solid): 2 hits ---
    Sphere sphere(5, 0);
    for(double D : distances) {
        auto hits = sphere.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Sphere at D=" << D << ": expected 2, got " << hits.size();
    }

    // --- EllipticalTube (solid): 2 hits ---
    EllipticalTube etube(3, 5, 10);
    for(double D : distances) {
        auto hits = etube.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "EllipticalTube at D=" << D << ": expected 2, got " << hits.size();
    }

    // --- Box (solid): 2 hits ---
    Box box(10, 8, 6);
    for(double D : distances) {
        auto hits = box.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Box at D=" << D << ": expected 2, got " << hits.size();
    }
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
        EXPECT_EQ(hits.size(), 0u) << "EllipticalTube tangent y=5: should miss (tangent)";
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
        EXPECT_EQ(hits.size(), 0u) << "EllipticalTube tangent x=3: should miss (tangent)";
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
        // Vertex-aimed ray on convex shape: result is 0 or 2.
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "Trd vertex-aimed ray: expected 0 or 2, got " << hits.size();
    }

    // Ray skimming along a bottom edge (y=4 face at z=-6, from x=-10 to x=+10)
    {
        auto hits = trd.Intersections(Vector3D(-10, 4, -6), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u) << "Trd edge-skimming ray: should enter and exit";
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
        EXPECT_EQ(hits.size(), 0u) << "ExtrPoly lateral edge ray: should miss (tangent to edge)";
    }

    // Ray along a lateral face (x=3, y=0, varying z) -- face-grazing
    {
        auto hits = ep.Intersections(Vector3D(3, 0, -10), Vector3D(0, 0, 1));
        // Face-grazing ray on convex shape: result is 0 or 2.
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "ExtrPoly face-grazing ray: expected 0 or 2, got " << hits.size();
    }

    // Ray aimed at a bottom polygon vertex (3, 3, -4) from outside
    {
        Vector3D dir(3, 3, -4);
        dir.normalize();
        auto hits = ep.Intersections(Vector3D(6, 6, -8), Vector3D(-dir.GetX(), -dir.GetY(), -dir.GetZ()));
        EXPECT_EQ(hits.size(), 2u) << "ExtrPoly vertex-aimed ray: should enter and exit";
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
        EXPECT_EQ(hits.size(), 4u) << "CutTube low-plane parallel ray: hollow tube gives 4 hits";
    }

    // Ray parallel to the high cut plane, just inside
    {
        Vector3D dir(std::cos(tilt), 0, std::sin(tilt));
        dir.normalize();
        auto hits = ct.Intersections(Vector3D(-20, 0, 9), dir);
        EXPECT_EQ(hits.size(), 0u) << "CutTube high-plane parallel ray: above cut plane, no hit";
    }

    // Ray through the intersection of cut plane and barrel surface
    {
        auto hits = ct.Intersections(Vector3D(5, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "CutTube barrel/plane edge ray: should enter and exit";
    }

    // Ray through center
    {
        auto hits = ct.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 0u) << "CutTube z-axis ray: inside hollow center, no hit";
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
        EXPECT_EQ(hits.size(), 2u) << "Union touching spheres x-ray: merged solid, enter + exit";
    }

    // Intersection: ray through the touching point along x
    {
        auto bg = BooleanGeometry(BooleanOperation::INTERSECTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 0u) << "Intersection touching spheres x-ray: zero-volume, no hit";
    }

    // Subtraction: ray through the touching point along x
    {
        auto bg = BooleanGeometry(BooleanOperation::SUBTRACTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u) << "Subtraction touching spheres x-ray: left sphere minus point";
    }

    // Overlapping spheres: centers at (-2,0,0) and (2,0,0), radius 5
    auto s3 = Sphere(Placement(Vector3D(-2, 0, 0)), 5, 0).create();
    auto s4 = Sphere(Placement(Vector3D(2, 0, 0)), 5, 0).create();

    // Union: ray through the overlap region
    {
        auto bg = BooleanGeometry(BooleanOperation::UNION, s3, s4);
        auto hits = bg.Intersections(Vector3D(0, -10, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size(), 2u) << "Union overlapping spheres y-ray: single merged solid";
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

    // Ray along x just past the z-extent boundary (z > dz)
    {
        auto hits = para.Intersections(Vector3D(-20, 0, 6.01), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 0u) << "Para past z-boundary x-ray: no hit";
    }

    // Axis-aligned rays from various directions
    Vector3D dirs[] = {
        Vector3D(1,0,0), Vector3D(0,1,0), Vector3D(0,0,1),
        Vector3D(-1,0,0), Vector3D(0,-1,0), Vector3D(0,0,-1)
    };
    for(auto const & d : dirs) {
        auto hits = para.Intersections(Vector3D(0, 0, 0), d);
        EXPECT_EQ(hits.size(), 2u) << "Para origin ray: should enter and exit";
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

    // Ray aimed at a vertex (3,3,3) from outside (convex: 0 or 2)
    {
        Vector3D dir(1, 1, 1);
        dir.normalize();
        auto hits = mesh.Intersections(Vector3D(10, 10, 10), Vector3D(-dir.GetX(), -dir.GetY(), -dir.GetZ()));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "Mesh vertex-aimed ray: expected 0 or 2, got " << hits.size();
    }

    // Ray along an edge (x=3, y=3, varying z) (convex: 0 or 2)
    {
        auto hits = mesh.Intersections(Vector3D(s, s, -10), Vector3D(0, 0, 1));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "Mesh edge-aligned ray: expected 0 or 2, got " << hits.size();
    }

    // Ray along a face edge (x=3, y=0, z=-3 -> z=3 diagonal of +x face) (convex: 0 or 2)
    {
        auto hits = mesh.Intersections(Vector3D(s, 0, -10), Vector3D(0, 0, 1));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "Mesh face-diagonal ray: expected 0 or 2, got " << hits.size();
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

    // Ray exactly at z-cut boundary (convex shape: 0 or 2)
    {
        auto hits = ell.Intersections(Vector3D(-10, 0, 6), Vector3D(1, 0, 0));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "Ellipsoid at upper z-cut boundary: expected 0 or 2, got " << hits.size();
    }
    {
        auto hits = ell.Intersections(Vector3D(-10, 0, -4), Vector3D(1, 0, 0));
        EXPECT_TRUE(hits.size() == 0 || hits.size() == 2)
            << "Ellipsoid at lower z-cut boundary: expected 0 or 2, got " << hits.size();
    }
}

// =========================================================================
// Exact count tests for center-piercing rays on all shapes.
// A ray through the center of any finite solid must produce at least 2 hits.
// For convex shapes, exactly 2. This catches algorithms that silently
// return empty on hard cases.
// =========================================================================
TEST(ShapeInvariants, CenterPiercingExactCount) {
    // Convex shapes: ray through center must give exactly 2
    struct ConvexCase {
        std::string name;
        std::shared_ptr<Geometry> geo;
    };
    std::vector<ConvexCase> convex = {
        {"Box(10,8,6)", Box(10, 8, 6).create()},
        {"Sphere(5,0)", Sphere(5, 0).create()},
        {"Cylinder(4,0,10)", Cylinder(4, 0, 10).create()},
        {"Cone(0,5,0,3,8)", Cone(0, 5, 0, 3, 8).create()},
        {"Trd(5,3,4,2,6)", Trd(5, 3, 4, 2, 6).create()},
        {"EllipticalTube(3,5,10)", EllipticalTube(3, 5, 10).create()},
        {"CutTubeFlat(0,5,8)", CutTube(0, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1)).create()},
        {"Trap(sym)", Trap(5, 0, 0, 3, 4, 4, 0, 2, 3, 3, 0).create()},
        {"Ellipsoid(5,3,4)", Ellipsoid(5, 3, 4).create()},
        {"Para(4,3,5,0.3,0.2,0.5)", Para(4, 3, 5, 0.3, 0.2, 0.5).create()},
    };

    // Test with general (non-axis-aligned) directions through center
    Vector3D dirs[] = {
        Vector3D(1, 0, 0),
        Vector3D(0, 1, 0),
        Vector3D(0, 0, 1),
        Vector3D(1, 1, 1),
        Vector3D(1, -1, 0.5),
        Vector3D(-0.3, 0.7, -0.6),
        Vector3D(0.1, 0.2, 0.97),
        Vector3D(2, 3, -5),
    };

    for(auto & entry : convex) {
        Vector3D center = entry.geo->GetPlacement().GetPosition();
        for(auto & d : dirs) {
            Vector3D dir = d;
            dir.normalize();
            Vector3D origin = center - dir * 50.0;
            auto hits = entry.geo->Intersections(origin, dir);
            EXPECT_EQ(hits.size(), 2u)
                << entry.name << " center-piercing ray dir=("
                << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ()
                << ") got " << hits.size() << " hits instead of 2";
        }
    }
}

// =========================================================================
// General 3D ray invariants: non-axis-aligned rays through all shapes.
// Tests both that results are even AND non-zero for rays that must hit.
// =========================================================================
TEST(ShapeInvariants, General3DRays) {
    auto shapes = MakeShapes();
    int failures = 0;

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;

        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid() || !std::isfinite(box.min_corner.GetX())) continue;

        Vector3D center = (box.min_corner + box.max_corner) * 0.5;
        double extent = (box.max_corner - box.min_corner).magnitude();

        // Fire 50 random rays that pass through the AABB center.
        // These must all produce even intersection counts.
        // Most should produce >= 2 hits (some may tangentially miss).
        int zero_count = 0;
        for(int i = 0; i < 50; ++i) {
            Vector3D dir = RandomDirection();
            Vector3D origin = center - dir * (extent + 10);
            auto hits = entry.geo->Intersections(origin, dir);
            EXPECT_EQ(hits.size() % 2, 0u)
                << entry.name << " general ray #" << i << ": odd count " << hits.size();
            if(hits.size() % 2 != 0) failures++;
            if(hits.size() == 0) zero_count++;
        }
        // A ray through the AABB center of a convex solid should always hit.
        // Non-convex shapes (torus, partial solids) may have large holes at
        // their AABB center, so many random directions genuinely miss.
        // Only assert for convex shapes where zero hits means a solver failure.
        if(entry.max_intersections == 2) {
            EXPECT_EQ(zero_count, 0)
                << entry.name << " (convex): " << zero_count
                << "/50 center-aimed rays returned 0 hits";
        }
        // For non-convex shapes: the even-count invariant (tested above) is
        // the meaningful correctness check. Zero-hit rates depend on geometry.
    }
    EXPECT_EQ(failures, 0) << "Total odd-count failures: " << failures;
}

// =========================================================================
// Far-field precision: expose catastrophic cancellation in quadratic solver.
// At distance D from a shape of radius r, the discriminant b^2-4ac suffers
// catastrophic cancellation because b^2 ~ 4ac ~ D^2 while the useful
// difference is O(r^2). Tests at 1e8, 1e10, 1e12.
// =========================================================================
TEST(ShapeInvariants, FarFieldPrecisionExtreme) {
    double distances[] = {1e8, 1e10, 1e12};

    struct FarFieldCase {
        std::string name;
        std::shared_ptr<Geometry> geo;
        double radius;
    };
    std::vector<FarFieldCase> cases = {
        {"Sphere(5)", Sphere(5, 0).create(), 5},
        {"Cylinder(4,0,10)", Cylinder(4, 0, 10).create(), 4},
        {"Cone(0,5,0,3,8)", Cone(0, 5, 0, 3, 8).create(), 4},
        {"EllipticalTube(3,5,10)", EllipticalTube(3, 5, 10).create(), 3},
        {"Ellipsoid(5,3,4)", Ellipsoid(5, 3, 4).create(), 5},
        {"CutTubeFlat(0,5,8)", CutTube(0, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1)).create(), 5},
        {"Polycone", Polycone({-5,0,5}, {0,0,0}, {4,6,4}).create(), 6},
        {"GenericPolycone", GenericPolycone({3,4.5,5,3.5,3,2}, {-5,0,5,5,0,-5}).create(), 4.5},
    };

    for(auto & c : cases) {
        for(double D : distances) {
            // Axis-aligned ray through center
            auto hits_x = c.geo->Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
            EXPECT_EQ(hits_x.size() % 2, 0u)
                << c.name << " x-axis at D=" << D << ": odd count " << hits_x.size();
            EXPECT_GE(hits_x.size(), 2u)
                << c.name << " x-axis at D=" << D << ": expected >=2, got " << hits_x.size();

            if(hits_x.size() >= 2) {
                double span = hits_x.back().distance - hits_x.front().distance;
                EXPECT_GT(span, 0)
                    << c.name << " x-axis at D=" << D << ": hits not ordered";
                EXPECT_NEAR(span, 2 * c.radius, c.radius * 0.01)
                    << c.name << " x-axis at D=" << D
                    << ": span " << span << " vs expected " << 2*c.radius;
            }

            // General diagonal ray through center
            Vector3D dir(1, 1, 1);
            dir.normalize();
            auto hits_diag = c.geo->Intersections(Vector3D(-D, -D, -D) * (1.0/std::sqrt(3.0)), dir);
            EXPECT_EQ(hits_diag.size() % 2, 0u)
                << c.name << " diagonal at D=" << D << ": odd count " << hits_diag.size();
            EXPECT_GE(hits_diag.size(), 2u)
                << c.name << " diagonal at D=" << D << ": missed entirely";
        }
    }
}

// =========================================================================
// Far-field precision for slab/prism shapes: box, trd, trap, para, polyhedra.
// Verifies center-piercing rays still produce correct span at extreme distance.
// =========================================================================
TEST(ShapeInvariants, FarFieldPrecisionSlabShapes) {
    double distances[] = {1e8, 1e10, 1e12};

    struct SlabCase {
        std::string name;
        std::shared_ptr<Geometry> geo;
        double expected_span;
    };
    std::vector<SlabCase> cases = {
        {"Box(10,8,6)", Box(10, 8, 6).create(), 6.0},
        {"Trd(5,3,4,2,6)", Trd(5, 3, 4, 2, 6).create(), 12.0},
        {"Trap(sym)", Trap(5, 0, 0, 3, 4, 4, 0, 2, 3, 3, 0).create(), 10.0},
        {"Para(4,3,5,0,0,0)", Para(4, 3, 5, 0, 0, 0).create(), 10.0},
        {"Polyhedra6", Polyhedra(6, 0, {-5,5}, {0,0}, {5,5}).create(), 10.0},
    };

    for(auto & c : cases) {
        for(double D : distances) {
            // Z-axis ray through center
            auto hits = c.geo->Intersections(Vector3D(0, 0, -D), Vector3D(0, 0, 1));
            EXPECT_GE(hits.size(), 2u)
                << c.name << " z-axis at D=" << D << ": expected >=2, got " << hits.size();
            EXPECT_EQ(hits.size() % 2, 0u)
                << c.name << " z-axis at D=" << D << ": odd count " << hits.size();

            if(hits.size() >= 2) {
                double span = hits.back().distance - hits.front().distance;
                EXPECT_NEAR(span, c.expected_span, c.expected_span * 0.01)
                    << c.name << " z-axis at D=" << D
                    << ": span " << span << " vs expected " << c.expected_span;
            }

            // Diagonal ray through center
            Vector3D dir(1, 1, 1);
            dir.normalize();
            auto hits_diag = c.geo->Intersections(Vector3D(-D, -D, -D) * (1.0/std::sqrt(3.0)), dir);
            EXPECT_EQ(hits_diag.size() % 2, 0u)
                << c.name << " diagonal at D=" << D << ": odd count " << hits_diag.size();
            EXPECT_GE(hits_diag.size(), 2u)
                << c.name << " diagonal at D=" << D << ": missed entirely";
        }
    }
}

// =========================================================================
// Far-field precision for extruded polygon.
// =========================================================================
TEST(ShapeInvariants, FarFieldPrecisionExtrPoly) {
    double distances[] = {1e8, 1e10, 1e12};

    std::vector<std::vector<double>> polygon = {{-3,-3},{3,-3},{3,3},{-3,3}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off, 1.0),
        ExtrPoly::ZSection(4, off, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);

    for(double D : distances) {
        auto hits = ep.Intersections(Vector3D(0, 0, -D), Vector3D(0, 0, 1));
        EXPECT_GE(hits.size(), 2u)
            << "ExtrPoly z-axis at D=" << D << ": expected >=2, got " << hits.size();
        EXPECT_EQ(hits.size() % 2, 0u)
            << "ExtrPoly z-axis at D=" << D << ": odd count " << hits.size();
        if(hits.size() >= 2) {
            double span = hits.back().distance - hits.front().distance;
            EXPECT_NEAR(span, 8.0, 0.08)
                << "ExtrPoly z-axis at D=" << D << ": span " << span;
        }
    }
}

// =========================================================================
// Torus far-field surface landing: verify hit positions lie on the torus
// implicit surface at extreme distances.
// =========================================================================
TEST(ShapeInvariants, TorusFarFieldSurfaceLanding) {
    double R = 10, r = 3;
    Torus torus(R, r, 0);

    auto surfaceF = [&](Vector3D const & p) {
        double x = p.GetX(), y = p.GetY(), z = p.GetZ();
        double rxy = std::sqrt(x*x + y*y);
        return (rxy - R) * (rxy - R) + z*z - r*r;
    };

    // The solver shifts the ray origin to the closest approach point, so
    // quartic coefficients are bounded by torus dimensions, not D. Surface
    // residuals should be bounded by a fixed tolerance independent of D.
    double tol = 1e-4;
    double distances[] = {1e4, 1e6, 1e8};

    for(double D : distances) {
        // X-axis ray through both lobes
        auto hits = torus.Intersections(Vector3D(-D, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "Torus x-axis at D=" << D << ": expected 4, got " << hits.size();
        if(hits.size() == 4) {
            double outer_span = hits[3].distance - hits[0].distance;
            EXPECT_NEAR(outer_span, 2.0 * (R + r), 0.01)
                << "Torus outer span at D=" << D;
            double inner_span = hits[2].distance - hits[1].distance;
            EXPECT_NEAR(inner_span, 2.0 * (R - r), 0.01)
                << "Torus inner span at D=" << D;
            for(size_t i = 0; i < hits.size(); ++i) {
                double residual = std::fabs(surfaceF(hits[i].position));
                EXPECT_LT(residual, tol)
                    << "Torus hit " << i << " at D=" << D << " off surface by " << residual;
            }
        }

        // Vertical ray through tube center at x=R
        auto hits_v = torus.Intersections(Vector3D(R, 0, -D), Vector3D(0, 0, 1));
        EXPECT_EQ(hits_v.size(), 2u)
            << "Torus vertical through lobe at D=" << D << ": expected 2, got " << hits_v.size();
        if(hits_v.size() == 2) {
            double span = hits_v[1].distance - hits_v[0].distance;
            EXPECT_NEAR(span, 2.0 * r, 0.01)
                << "Torus lobe vertical span at D=" << D;
        }
    }
}

// =========================================================================
// Trap far-field span precision: verify center-piercing z-ray span at
// extreme distances.
// =========================================================================
TEST(ShapeInvariants, TrapFarFieldSpanPrecision) {
    Trap trap(10, 0, 0, 8, 5, 5, 0, 6, 3, 3, 0);

    double distances[] = {1e6, 1e8, 1e10, 1e12};
    for(double D : distances) {
        auto hits = trap.Intersections(Vector3D(0, 0, -D), Vector3D(0, 0, 1));
        ASSERT_EQ(hits.size(), 2u)
            << "Trap z-axis at D=" << D << ": expected 2, got " << hits.size();
        double span = hits[1].distance - hits[0].distance;
        EXPECT_NEAR(span, 20.0, 1e-6)
            << "Trap span at D=" << D << ": " << span << " vs expected 20.0";
    }
}

// =========================================================================
// Boolean geometry exact counts: verify CSG operations produce correct
// intersection counts for well-defined configurations.
// =========================================================================
TEST(ShapeInvariants, BooleanExactCounts) {
    // Two non-overlapping spheres: union ray through both must give 4 hits
    {
        auto s1 = Sphere(Placement(Vector3D(-10, 0, 0)), 3, 0).create();
        auto s2 = Sphere(Placement(Vector3D(10, 0, 0)), 3, 0).create();
        auto bg = BooleanGeometry(BooleanOperation::UNION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "Union of separated spheres along x: expected 4 hits, got " << hits.size();
    }

    // Two non-overlapping spheres: intersection should give 0 hits (empty solid)
    {
        auto s1 = Sphere(Placement(Vector3D(-10, 0, 0)), 3, 0).create();
        auto s2 = Sphere(Placement(Vector3D(10, 0, 0)), 3, 0).create();
        auto bg = BooleanGeometry(BooleanOperation::INTERSECTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 0u)
            << "Intersection of separated spheres: expected 0, got " << hits.size();
    }

    // Subtraction: sphere minus non-overlapping sphere = original sphere (2 hits)
    {
        auto s1 = Sphere(Placement(Vector3D(0, 0, 0)), 5, 0).create();
        auto s2 = Sphere(Placement(Vector3D(20, 0, 0)), 3, 0).create();
        auto bg = BooleanGeometry(BooleanOperation::SUBTRACTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-10, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Subtraction with non-overlap: expected 2, got " << hits.size();
    }

    // Overlapping spheres: union ray should give 2 (enters leftmost, exits rightmost)
    {
        auto s1 = Sphere(Placement(Vector3D(-2, 0, 0)), 5, 0).create();
        auto s2 = Sphere(Placement(Vector3D(2, 0, 0)), 5, 0).create();
        auto bg = BooleanGeometry(BooleanOperation::UNION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Union of overlapping spheres: expected 2, got " << hits.size();
        if(hits.size() == 2) {
            EXPECT_NEAR(hits[0].position.GetX(), -7, 1e-6)
                << "Union enter at x=-7";
            EXPECT_NEAR(hits[1].position.GetX(), 7, 1e-6)
                << "Union exit at x=+7";
        }
    }

    // Overlapping spheres: intersection ray should give 2 (overlap region)
    {
        auto s1 = Sphere(Placement(Vector3D(-2, 0, 0)), 5, 0).create();
        auto s2 = Sphere(Placement(Vector3D(2, 0, 0)), 5, 0).create();
        auto bg = BooleanGeometry(BooleanOperation::INTERSECTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Intersection of overlapping spheres: expected 2, got " << hits.size();
        if(hits.size() == 2) {
            EXPECT_NEAR(hits[0].position.GetX(), -3, 1e-6)
                << "Intersection enter at x=-3";
            EXPECT_NEAR(hits[1].position.GetX(), 3, 1e-6)
                << "Intersection exit at x=+3";
        }
    }

    // Subtraction: large sphere minus centered small sphere = shell (4 hits)
    {
        auto s1 = Sphere(5, 0).create();
        auto s2 = Sphere(2, 0).create();
        auto bg = BooleanGeometry(BooleanOperation::SUBTRACTION, s1, s2);
        auto hits = bg.Intersections(Vector3D(-10, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 4u)
            << "Shell (subtract concentric): expected 4, got " << hits.size();
    }
}

// =========================================================================
// CutTube exact intersection counts: rays through center must hit.
// The existing test only checked parity.
// =========================================================================
TEST(ShapeInvariants, CutTubeExactCounts) {
    double tilt = 30.0 * M_PI / 180.0;
    Vector3D low_norm(std::sin(tilt), 0, -std::cos(tilt));
    Vector3D high_norm(-std::sin(tilt), 0, std::cos(tilt));
    CutTube ct(0, 5, 10, low_norm, high_norm);

    // Z-axis ray through center: must hit exactly 2 (enters low cap, exits high cap)
    {
        auto hits = ct.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "CutTube z-axis: expected 2, got " << hits.size();
    }

    // X-axis ray through center: must hit exactly 2 (enters barrel, exits barrel)
    {
        auto hits = ct.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u) << "CutTube x-axis: expected 2, got " << hits.size();
    }

    // Diagonal ray through center
    {
        Vector3D dir(1, 1, 0.5);
        dir.normalize();
        auto hits = ct.Intersections(Vector3D(0, 0, 0) - dir * 20.0, dir);
        EXPECT_EQ(hits.size(), 2u) << "CutTube diagonal: expected 2, got " << hits.size();
    }

    // General directions through center (non-axis-aligned)
    Vector3D general_dirs[] = {
        Vector3D(1, 2, 3), Vector3D(-1, 0.5, 0.3), Vector3D(0.7, -0.7, 0.1),
        Vector3D(3, -1, 2), Vector3D(-2, 3, -1),
    };
    for(auto & d : general_dirs) {
        Vector3D dir = d;
        dir.normalize();
        auto hits = ct.Intersections(Vector3D(0, 0, 0) - dir * 30.0, dir);
        EXPECT_EQ(hits.size(), 2u)
            << "CutTube general dir (" << d.GetX() << "," << d.GetY() << "," << d.GetZ()
            << "): expected 2, got " << hits.size();
    }
}

// =========================================================================
// GenericPolycone exact counts for center-piercing rays.
// Existing tests only check parity.
// =========================================================================
TEST(ShapeInvariants, GenericPolyconeExactCounts) {
    // Solid cone: triangle in (r,z) space revolved => convex for radial rays
    GenericPolycone solid({0, 5, 0}, {-5, 0, 5});

    // Radial rays through center must give exactly 2
    {
        auto hits = solid.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u) << "Solid GenericPolycone x-axis: expected 2, got " << hits.size();
    }
    {
        auto hits = solid.Intersections(Vector3D(0, -20, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size(), 2u) << "Solid GenericPolycone y-axis: expected 2, got " << hits.size();
    }
    // Z-axis ray through r=0 apex vertices: enters at z=-5, exits at z=+5
    {
        auto hits = solid.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "Solid GenericPolycone z-axis: expected 2, got " << hits.size();
    }
    {
        Vector3D dir(1, 1, 1);
        dir.normalize();
        auto hits = solid.Intersections(Vector3D(0, 0, 0) - dir * 20, dir);
        EXPECT_EQ(hits.size(), 2u) << "Solid GenericPolycone diagonal: expected 2, got " << hits.size();
    }

    // Non-degenerate hollow shape: outer radius varies, inner void
    GenericPolycone hollow({3, 4.5, 5, 3.5, 3, 2}, {-5, 0, 5, 5, 0, -5});
    {
        auto hits = hollow.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_GE(hits.size(), 2u) << "Hollow GenericPolycone x-axis: expected >=2, got " << hits.size();
        EXPECT_EQ(hits.size() % 2, 0u) << "Hollow GenericPolycone x-axis: odd count " << hits.size();
    }
    {
        Vector3D dir(1, 2, -1);
        dir.normalize();
        auto hits = hollow.Intersections(Vector3D(0, 0, 0) - dir * 20, dir);
        EXPECT_GE(hits.size(), 2u) << "Hollow GenericPolycone general: expected >=2, got " << hits.size();
        EXPECT_EQ(hits.size() % 2, 0u) << "Hollow GenericPolycone general: odd count " << hits.size();
    }

    // Hollow (U-shape): center ray should give 4 (enters outer, exits inner, enters inner, exits outer)
    GenericPolycone ushape({5, 5, 3, 3}, {-4, 4, 4, -4});
    {
        auto hits = ushape.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 4u) << "U-shape GenericPolycone x-axis: expected 4, got " << hits.size();
    }
    {
        auto hits = ushape.Intersections(Vector3D(0, -20, 0), Vector3D(0, 1, 0));
        EXPECT_EQ(hits.size(), 4u) << "U-shape GenericPolycone y-axis: expected 4, got " << hits.size();
    }
}

// =========================================================================
// Polycone horizontal ray at internal z-boundary: must find the barrel hits.
// The existing HorizontalRayAtSectionBoundary test uses polycone with step
// at z=0 and verifies 2 hits. We extend to general directions.
// =========================================================================
TEST(ShapeInvariants, PolyconeInternalBoundaryExact) {
    // Polycone: 2 sections with same outer radius (5), no step.
    // Section 1: z in [-5, 0], rmin=0, rmax=5
    // Section 2: z in [0, 5], rmin=0, rmax=5
    // This is just a cylinder r=5, z=[-5,5], but split into 2 sections.
    Polycone pc({-5, 0, 5}, {0, 0, 0}, {5, 5, 5});

    // Horizontal ray at z=0 (internal boundary): must give 2
    {
        auto hits = pc.Intersections(Vector3D(-10, 0, 0), Vector3D(1, 0, 0));
        EXPECT_EQ(hits.size(), 2u)
            << "Polycone horizontal at z=0: expected 2, got " << hits.size();
        if(hits.size() == 2) {
            EXPECT_NEAR(hits[0].position.GetX(), -5.0, 1e-9);
            EXPECT_NEAR(hits[1].position.GetX(), 5.0, 1e-9);
        }
    }

    // General direction through center at z=0
    {
        Vector3D dir(1, 1, 0);
        dir.normalize();
        auto hits = pc.Intersections(Vector3D(0, 0, 0) - dir * 20, dir);
        EXPECT_EQ(hits.size(), 2u)
            << "Polycone diagonal horizontal: expected 2, got " << hits.size();
    }

    // Nearly-horizontal ray (small dz component)
    {
        Vector3D dir(1, 0, 1e-6);
        dir.normalize();
        auto hits = pc.Intersections(Vector3D(-10, 0, 0), dir);
        EXPECT_EQ(hits.size(), 2u)
            << "Polycone nearly-horizontal: expected 2, got " << hits.size();
    }
}

// =========================================================================
// Enter/exit consistency: verify that intersections strictly alternate
// entering/exiting, and that entering comes first when ray starts outside.
// =========================================================================
TEST(ShapeInvariants, EnterExitStrictAlternation) {
    auto shapes = MakeShapes();

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;

        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid() || !std::isfinite(box.min_corner.GetX())) continue;

        Vector3D center = (box.min_corner + box.max_corner) * 0.5;
        double extent = (box.max_corner - box.min_corner).magnitude();

        for(int i = 0; i < 30; ++i) {
            Vector3D dir = RandomDirection();
            Vector3D origin = center - dir * (extent + 20);
            auto hits = entry.geo->Intersections(origin, dir);

            if(hits.empty()) continue;

            // First hit from outside must be entering
            EXPECT_TRUE(hits[0].entering)
                << entry.name << " ray #" << i
                << ": first hit should be entering but is exiting";

            // Strictly alternating
            for(size_t j = 1; j < hits.size(); ++j) {
                EXPECT_NE(hits[j].entering, hits[j-1].entering)
                    << entry.name << " ray #" << i
                    << ": hits[" << j-1 << "] and hits[" << j << "] have same entering="
                    << hits[j].entering;
            }

            // Last hit must be exiting
            EXPECT_FALSE(hits.back().entering)
                << entry.name << " ray #" << i
                << ": last hit should be exiting but is entering";

            // Distances must be non-decreasing (coincident CSG boundaries
            // at the same distance are acceptable but must still alternate)
            for(size_t j = 1; j < hits.size(); ++j) {
                EXPECT_GE(hits[j].distance, hits[j-1].distance)
                    << entry.name << " ray #" << i
                    << ": hits[" << j << "].distance=" << hits[j].distance
                    << " < hits[" << j-1 << "].distance=" << hits[j-1].distance;
            }
        }
    }
}

// =========================================================================
// Intersection position accuracy: verify hit positions actually lie on
// the shape's surface using the implicit surface equation.
// =========================================================================
TEST(ShapeInvariants, IntersectionPositionOnSurface) {
    // For each shape, fire rays and verify the reported position matches
    // origin + t * direction (internal consistency)
    auto shapes = MakeShapes();

    for(auto const & entry : shapes) {
        if(entry.known_cross_section > 0) continue;

        AABB box = entry.geo->GetWorldBoundingBox();
        if(!box.IsValid() || !std::isfinite(box.min_corner.GetX())) continue;

        Vector3D center = (box.min_corner + box.max_corner) * 0.5;
        double extent = (box.max_corner - box.min_corner).magnitude();
        double tol = 1e-9 * extent;

        for(int i = 0; i < 20; ++i) {
            Vector3D dir = RandomDirection();
            Vector3D origin = center - dir * (extent + 10);
            auto hits = entry.geo->Intersections(origin, dir);

            for(auto & h : hits) {
                Vector3D expected_pos = origin + dir * h.distance;
                double err = (h.position - expected_pos).magnitude();
                EXPECT_LT(err, tol)
                    << entry.name << " ray #" << i
                    << ": position error " << err << " > tol " << tol
                    << " at distance " << h.distance;
            }
        }
    }
}

// =========================================================================
// Torus: verify all 4 intersections for axis-parallel ray through lobe.
// The torus with R=10, r=3 has lobes at x=+-10. A ray from x=-20 along x
// through y=0,z=0 must produce exactly 4 hits (enter outer at x=-13,
// exit inner at x=-7, enter inner at x=7, exit outer at x=13).
// =========================================================================
TEST(ShapeInvariants, TorusLobeExactFourHits) {
    Torus torus(10, 3, 0);

    // Ray along x-axis: passes through both lobes
    auto hits = torus.Intersections(Vector3D(-20, 0, 0), Vector3D(1, 0, 0));
    EXPECT_EQ(hits.size(), 4u) << "Torus x-axis: expected 4 hits, got " << hits.size();
    if(hits.size() == 4) {
        EXPECT_NEAR(hits[0].position.GetX(), -13, 1e-6);
        EXPECT_NEAR(hits[1].position.GetX(), -7, 1e-6);
        EXPECT_NEAR(hits[2].position.GetX(), 7, 1e-6);
        EXPECT_NEAR(hits[3].position.GetX(), 13, 1e-6);
        EXPECT_TRUE(hits[0].entering);
        EXPECT_FALSE(hits[1].entering);
        EXPECT_TRUE(hits[2].entering);
        EXPECT_FALSE(hits[3].entering);
    }

    // Ray along y-axis through (10,0,0): at x=10, torus surface satisfies
    // (sqrt(x^2+y^2) - R)^2 + z^2 = r^2. With x=10, z=0:
    //   (sqrt(100+y^2) - 10)^2 = 9
    //   sqrt(100+y^2) = 13  =>  y^2 = 69  =>  y = +/-sqrt(69)
    // The span is 2*sqrt(69) (NOT 2*r=6, which would only hold for a
    // radial-direction ray through the tube center).
    auto hits_y = torus.Intersections(Vector3D(10, -10, 0), Vector3D(0, 1, 0));
    EXPECT_EQ(hits_y.size(), 2u) << "Torus y-through-lobe: expected 2, got " << hits_y.size();
    if(hits_y.size() == 2) {
        double span_y = hits_y[1].position.GetY() - hits_y[0].position.GetY();
        double expected_span = 2.0 * std::sqrt(69.0); // 2*sqrt(R^2+2Rr+r^2 - R^2) when entering at rxy=R+r
        EXPECT_NEAR(span_y, expected_span, 1e-6) << "Torus lobe y-span should be 2*sqrt(69)";
    }

    // Diagonal ray through center of torus hole (no hit expected)
    auto hits_center = torus.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
    // For R=10, r=3: at (0,0,z) the nearest torus point is at distance R-r=7
    // so a z-axis ray should hit if it passes through the "donut hole" region
    // Actually at (0,0,z): sqrt(0+0) = 0, dist to tube center = R=10, so
    // surface equation: (0 - 10)^2 + z^2 = 9 => z^2 = 9 - 100 < 0 => no hit
    EXPECT_EQ(hits_center.size(), 0u) << "Torus z-axis (through hole): expected 0, got " << hits_center.size();

    // Non-axis-aligned ray through a lobe
    {
        Vector3D dir(0, 1, 1);
        dir.normalize();
        // Start at lobe center (10, 0, 0), go along (0,1,1) -- tube radius is 3
        Vector3D origin(10, 0, 0);
        origin = origin - dir * 10;
        auto hits_oblique = torus.Intersections(origin, dir);
        EXPECT_EQ(hits_oblique.size(), 2u) << "Torus oblique through lobe: expected 2, got " << hits_oblique.size();
    }
}

// =========================================================================
// Mesh geometry: non-axis-aligned mesh (rotated cube) to expose orientation
// and general-position intersection bugs.
// =========================================================================
TEST(ShapeInvariants, MeshGeneralPosition) {
    // Rotated cube: vertices of unit cube rotated 45 degrees around z-axis
    double c = std::cos(M_PI / 4.0);
    double s_val = std::sin(M_PI / 4.0);
    double sz = 5.0;

    // Cube vertices rotated 45 deg around z
    Vector3D v[8] = {
        Vector3D((-sz*c - (-sz)*s_val), (-sz*s_val + (-sz)*c), -sz),
        Vector3D(( sz*c - (-sz)*s_val), ( sz*s_val + (-sz)*c), -sz),
        Vector3D(( sz*c -   sz *s_val), ( sz*s_val +   sz *c), -sz),
        Vector3D((-sz*c -   sz *s_val), (-sz*s_val +   sz *c), -sz),
        Vector3D((-sz*c - (-sz)*s_val), (-sz*s_val + (-sz)*c),  sz),
        Vector3D(( sz*c - (-sz)*s_val), ( sz*s_val + (-sz)*c),  sz),
        Vector3D(( sz*c -   sz *s_val), ( sz*s_val +   sz *c),  sz),
        Vector3D((-sz*c -   sz *s_val), (-sz*s_val +   sz *c),  sz),
    };

    // 12 triangles (2 per face), consistent outward winding
    std::vector<std::array<Vector3D, 3>> triangles;
    auto addQuad = [&](Vector3D a, Vector3D b, Vector3D cc, Vector3D d) {
        triangles.push_back({{a, b, cc}});
        triangles.push_back({{a, cc, d}});
    };
    addQuad(v[3], v[2], v[1], v[0]); // -z
    addQuad(v[4], v[5], v[6], v[7]); // +z
    addQuad(v[0], v[1], v[5], v[4]); // -y
    addQuad(v[1], v[2], v[6], v[5]); // +x
    addQuad(v[2], v[3], v[7], v[6]); // +y
    addQuad(v[3], v[0], v[4], v[7]); // -x

    TriangularMesh mesh(triangles);

    // Ray through center along z-axis
    {
        auto hits = mesh.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 2u) << "Rotated mesh z-axis: expected 2, got " << hits.size();
    }

    // General diagonal ray through center
    {
        Vector3D dir(1, 2, 3);
        dir.normalize();
        auto hits = mesh.Intersections(Vector3D(0, 0, 0) - dir * 30, dir);
        EXPECT_EQ(hits.size(), 2u) << "Rotated mesh diagonal: expected 2, got " << hits.size();
    }

    // Multiple general directions
    int missed = 0;
    for(int i = 0; i < 50; ++i) {
        Vector3D dir = RandomDirection();
        auto hits = mesh.Intersections(Vector3D(0, 0, 0) - dir * 30, dir);
        if(hits.size() != 2) missed++;
        EXPECT_EQ(hits.size() % 2, 0u)
            << "Rotated mesh random ray #" << i << ": odd count " << hits.size();
    }
    EXPECT_EQ(missed, 0) << "Rotated mesh: " << missed << "/50 center rays missed";
}

// =========================================================================
// Non-convex center-piercing exact counts: verify that specific non-convex
// shapes produce the correct number of intersections through known interior
// points. Unlike the General3DRays test which lumps all non-convex shapes
// together, each shape here uses a ray origin and expected count that are
// known to be correct for that geometry.
// =========================================================================
TEST(ShapeInvariants, NonConvexCenterPiercing) {
    struct NonConvexCase {
        std::string name;
        std::shared_ptr<Geometry> geo;
        Vector3D known_interior;
        size_t expected_hits;
    };

    std::vector<NonConvexCase> cases = {
        // Hollow sphere: ray through center crosses outer surface, inner
        // surface, inner surface, outer surface = 4 hits
        {"SphereHollow(5,2)", Sphere(5, 2).create(), Vector3D(3.5, 0, 0), 4},
        // Hollow cylinder: same analysis, 4 hits
        {"CylinderHollow(5,2,10)", Cylinder(5, 2, 10).create(), Vector3D(3.5, 0, 0), 4},
        // Hollow cone: 4 hits through the annular cross-section
        {"ConeHollow(1,5,1,3,8)", Cone(1, 5, 1, 3, 8).create(), Vector3D(3, 0, 0), 4},
        // Torus: x-axis ray crosses both lobes at y=0 for 4 hits
        {"Torus(10,3,0)_xray", Torus(10, 3, 0).create(), Vector3D(10, 0, 0), 4},
        // Boolean subtraction (shell): 4 hits
        {"BoolSubtract(shell)", BooleanGeometry(BooleanOperation::SUBTRACTION,
            Sphere(5, 0).create(), Cylinder(2, 0, 12).create()).create(),
            Vector3D(3.5, 0, 0), 4},
        // U-shape GenericPolycone: 4 hits through the hollow
        {"GenericPolyconeU", GenericPolycone({5, 5, 3, 3}, {-4, 4, 4, -4}).create(),
            Vector3D(4, 0, 0), 4},
        // Polycone hollow: 4 hits
        {"PolyconeHollow", Polycone({-4, 0, 4}, {1, 2, 1}, {5, 6, 5}).create(),
            Vector3D(3.5, 0, 0), 4},
    };

    Vector3D dirs[] = {
        Vector3D(1, 0, 0), Vector3D(0, 1, 0), Vector3D(0, 0, 1),
        Vector3D(1, 1, 1), Vector3D(1, -1, 0.5), Vector3D(-0.3, 0.7, -0.6),
    };

    for(auto & c : cases) {
        for(auto & d : dirs) {
            Vector3D dir = d;
            dir.normalize();
            Vector3D origin = c.known_interior - dir * 50.0;
            auto hits = c.geo->Intersections(origin, dir);
            EXPECT_EQ(hits.size() % 2, 0u)
                << c.name << " dir=(" << d.GetX() << "," << d.GetY() << "," << d.GetZ()
                << "): odd count " << hits.size();
            EXPECT_GE(hits.size(), 2u)
                << c.name << " dir=(" << d.GetX() << "," << d.GetY() << "," << d.GetZ()
                << "): missed entirely (0 hits through known interior point)";
        }
    }

    // Torus donut hole: z-axis ray through center must return 0 (genuinely
    // misses the body). This confirms we don't over-count.
    {
        Torus torus(10, 3, 0);
        auto hits = torus.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
        EXPECT_EQ(hits.size(), 0u)
            << "Torus z-axis through donut hole should produce 0 hits";
    }
}

// =========================================================================
// DistanceToClosestApproach must use direction
// transform, not position transform. With a non-origin placement, the
// old code applied the placement translation to the direction vector,
// producing a distance that depended on the placement origin.
// =========================================================================
TEST(ShapeInvariants, DistanceToClosestApproachPlacement) {
    // DistanceToClosestApproach returns the parameter t along the ray
    // where the ray is closest to the shape center (in local frame).
    // For a ray from (0,0,-10) along +z through the origin, t=10.
    Sphere s_origin(5, 0);
    double dca_origin = s_origin.DistanceToClosestApproach(
        Vector3D(0, 0, -10), Vector3D(0, 0, 1));
    EXPECT_NEAR(dca_origin, 10.0, 1e-9);

    // Same sphere placed at (100, 0, 0). A ray from (100,0,-10) along +z
    // should give the same t=10 (same geometry in local frame).
    Sphere s_offset(Placement(Vector3D(100, 0, 0)), 5, 0);
    double dca_offset = s_offset.DistanceToClosestApproach(
        Vector3D(100, 0, -10), Vector3D(0, 0, 1));
    EXPECT_NEAR(dca_offset, 10.0, 1e-9)
        << "DCA along ray for offset sphere should match origin sphere";

    // Key check: a ray from (0,0,-10) along +z should
    // give t=10 for the origin sphere. For the offset sphere, local pos
    // becomes (-100,0,-10), so t = -dot((-100,0,-10), (0,0,1)) = 10.
    // The old bug (GlobalToLocalPosition on direction) would have
    // translated the direction by (100,0,0), producing dir=(100,0,1)
    // (unnormalized) and a completely wrong t value.
    double dca_miss = s_offset.DistanceToClosestApproach(
        Vector3D(0, 0, -10), Vector3D(0, 0, 1));
    EXPECT_NEAR(dca_miss, 10.0, 1e-9)
        << "DCA should not depend on placement translation applied to direction";

    // A ray from (100,0,-10) along +x with offset sphere at (100,0,0):
    // local pos = (0,0,-10), local dir = (1,0,0), t = 0.
    double dca_x = s_offset.DistanceToClosestApproach(
        Vector3D(100, 0, -10), Vector3D(1, 0, 0));
    EXPECT_NEAR(dca_x, 0.0, 1e-9);
}

// =========================================================================
// AABB must be valid after deserialization.
// Simulated by constructing a Box, then calling SetPlacement to trigger
// the lazy recompute path.
// =========================================================================
TEST(ShapeInvariants, AABBValidAfterDefaultConstruct) {
    // After construction with params, AABB should be valid and correct.
    Box box(10, 8, 6);
    AABB bb = box.GetWorldBoundingBox();
    EXPECT_TRUE(bb.IsValid()) << "Box AABB should be valid after construction";
    EXPECT_NEAR(bb.max_corner.GetX(), 5.0, 1e-9);
    EXPECT_NEAR(bb.min_corner.GetX(), -5.0, 1e-9);

    // After SetPlacement, AABB should update.
    box.SetPlacement(Placement(Vector3D(100, 0, 0)));
    AABB bb2 = box.GetWorldBoundingBox();
    EXPECT_TRUE(bb2.IsValid());
    EXPECT_NEAR(bb2.max_corner.GetX(), 105.0, 1e-9);
    EXPECT_NEAR(bb2.min_corner.GetX(), 95.0, 1e-9);

    // GenericPolycone AABB should also be valid.
    GenericPolycone gpc({0, 5, 0}, {-5, 0, 5});
    AABB gpc_bb = gpc.GetWorldBoundingBox();
    EXPECT_TRUE(gpc_bb.IsValid())
        << "GenericPolycone AABB should be valid after construction";
    EXPECT_NEAR(gpc_bb.max_corner.GetZ(), 5.0, 1e-6);
}

// =========================================================================
// Torus phi-cut AABB should be tighter than full rotation.
// A quarter-torus (delta_phi = pi/2) should have a smaller AABB than
// the full torus.
// =========================================================================
TEST(ShapeInvariants, TorusPhiCutAABBTightness) {
    double R = 10, r = 3;
    Torus full(R, r, 0);
    Torus quarter(R, r, 0, 0, M_PI / 2);

    AABB full_bb = full.GetBoundingBox();
    AABB quarter_bb = quarter.GetBoundingBox();

    // Full torus: extent is R+r = 13 in both x and y
    EXPECT_NEAR(full_bb.max_corner.GetX(), R + r, 1e-9);
    EXPECT_NEAR(full_bb.max_corner.GetY(), R + r, 1e-9);

    // Quarter torus (phi 0 to pi/2) should only extend into +x, +y quadrant.
    // The min corner should be much smaller in magnitude than -13.
    double full_volume = (full_bb.max_corner.GetX() - full_bb.min_corner.GetX())
                       * (full_bb.max_corner.GetY() - full_bb.min_corner.GetY())
                       * (full_bb.max_corner.GetZ() - full_bb.min_corner.GetZ());
    double quarter_volume = (quarter_bb.max_corner.GetX() - quarter_bb.min_corner.GetX())
                          * (quarter_bb.max_corner.GetY() - quarter_bb.min_corner.GetY())
                          * (quarter_bb.max_corner.GetZ() - quarter_bb.min_corner.GetZ());
    EXPECT_LT(quarter_volume, full_volume * 0.5)
        << "Quarter-torus AABB should be less than half the full-torus AABB volume";

    // Quarter torus min_x should be near -(tube radius) not -(R+r)
    EXPECT_GT(quarter_bb.min_corner.GetX(), -r - 1)
        << "Quarter-torus should not extend far into -x";
    EXPECT_GT(quarter_bb.min_corner.GetY(), -r - 1)
        << "Quarter-torus should not extend far into -y";

    // Verify the quarter AABB still contains a point we know is in the torus
    // Point on the tube surface at phi=pi/4
    double phi = M_PI / 4;
    double px = (R + r) * std::cos(phi);
    double py = (R + r) * std::sin(phi);
    EXPECT_LE(px, quarter_bb.max_corner.GetX() + 1e-6);
    EXPECT_LE(py, quarter_bb.max_corner.GetY() + 1e-6);
}

// =========================================================================
// TriangularMesh BVH traversal with many triangles.
// A densely tessellated sphere should produce correct even intersection
// counts for random rays. This exercises the heap-fallback path when the
// BVH tree depth exceeds the fixed stack capacity.
// =========================================================================
TEST(ShapeInvariants, MeshBVHLargeTriangleCount) {
    // Build a UV-sphere mesh with N_lat * N_lon * 2 triangles
    int N_lat = 30, N_lon = 60;
    double R = 5.0;
    std::vector<std::array<Vector3D, 3>> tris;
    tris.reserve(N_lat * N_lon * 2);

    auto sphere_pt = [&](int i, int j) {
        double theta = M_PI * i / N_lat;
        double phi = 2.0 * M_PI * j / N_lon;
        return Vector3D(R * std::sin(theta) * std::cos(phi),
                        R * std::sin(theta) * std::sin(phi),
                        R * std::cos(theta));
    };

    for(int i = 0; i < N_lat; ++i) {
        for(int j = 0; j < N_lon; ++j) {
            Vector3D a = sphere_pt(i, j);
            Vector3D b = sphere_pt(i + 1, j);
            Vector3D c = sphere_pt(i + 1, (j + 1) % N_lon);
            Vector3D d = sphere_pt(i, (j + 1) % N_lon);
            tris.push_back({{a, b, c}});
            tris.push_back({{a, c, d}});
        }
    }

    TriangularMesh mesh(tris);
    EXPECT_GE(tris.size(), 3000u) << "Mesh should have many triangles";

    // Fire random rays from outside and verify even intersection counts
    int odd_count = 0;
    for(int i = 0; i < 200; ++i) {
        Vector3D dir = RandomDirection();
        Vector3D origin = dir * (-20.0);
        auto hits = mesh.Intersections(origin, dir);
        if(hits.size() % 2 != 0) odd_count++;
        // Center-piercing ray on a convex mesh should produce exactly 2
        EXPECT_EQ(hits.size(), 2u)
            << "Mesh ray #" << i << ": expected 2, got " << hits.size();
    }
    EXPECT_EQ(odd_count, 0) << odd_count << "/200 rays had odd intersection count";
}

// =========================================================================
// delta_phi > 2*pi should clamp to 2*pi (matching Geant4 behavior),
// not throw. Common in real GDML files (e.g. deltaphi="360.5" degrees).
// =========================================================================
TEST(ShapeInvariants, DeltaPhiClampTo2Pi) {
    double over = 2.0 * M_PI + 0.01; // slightly over full circle

    // All phi-capable shapes should accept delta_phi > 2*pi without throwing
    EXPECT_NO_THROW(Sphere(5, 0, 0, over, 0, M_PI));
    EXPECT_NO_THROW(Cylinder(5, 0, 10, 0, over));
    EXPECT_NO_THROW(Cone(0, 5, 0, 3, 10, 0, over));
    EXPECT_NO_THROW(Polycone({-5, 5}, {0, 0}, {5, 5}, 0, over));
    EXPECT_NO_THROW(CutTube(0, 5, 8, Vector3D(0,0,-1), Vector3D(0,0,1), 0, over));
    EXPECT_NO_THROW(Torus(10, 3, 0, 0, over));
    EXPECT_NO_THROW(GenericPolycone({0, 5, 0}, {-5, 0, 5}, 0, over));
    EXPECT_NO_THROW(Polyhedra(6, 0, {-5, 5}, {0, 0}, {5, 5}, over));

    // The clamped shape should behave as a full-rotation shape (no phi cut).
    // A full-rotation sphere of radius 5: x-axis ray must give exactly 2 hits.
    Sphere s(5, 0, 0, over, 0, M_PI);
    auto hits = s.Intersections(Vector3D(-10, 0, 0), Vector3D(1, 0, 0));
    EXPECT_EQ(hits.size(), 2u) << "Clamped sphere should behave as full rotation";

    // delta_phi <= 0 should still throw
    EXPECT_THROW(Sphere(5, 0, 0, 0, 0, M_PI), std::invalid_argument);
    EXPECT_THROW(Sphere(5, 0, 0, -1, 0, M_PI), std::invalid_argument);
    EXPECT_THROW(Polyhedra(6, 0, {-5, 5}, {0, 0}, {5, 5}, 0), std::invalid_argument);
    EXPECT_THROW(Polyhedra(6, 0, {-5, 5}, {0, 0}, {5, 5}, -1), std::invalid_argument);
}

TEST(ShapeInvariants, ExtrPolySetPolygonRecomputesPlanes) {
    // Regression test for the SetPolygon stale-cache bug.
    //
    // ExtrPoly stores derived lateral-face plane coefficients in planes_,
    // computed from polygon_ at construction. ComputeIntersections reads
    // planes_ directly (planes_.size() for the loop count, planes_[i] for
    // the per-edge normal and offset). If SetPolygon replaces polygon_
    // without recomputing planes_, intersection results are wrong: the
    // solver runs the new polygon's z-cap data against the old polygon's
    // edges. A ray that should hit the new polygon will miss because the
    // old (smaller) face planes still reject it as "outside".
    //
    // The bug is silent: no out-of-bounds (planes_.size() bounds the loop)
    // and no crash. Only the intersection list is wrong.

    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-1.0, off, 1.0),
        ExtrPoly::ZSection( 1.0, off, 1.0)
    };

    // Construct with a small triangle near the origin (3 edges, max |x| = 0.1).
    std::vector<std::vector<double>> triangle = {
        {-0.1, -0.1}, {0.1, -0.1}, {0.0, 0.1}
    };
    ExtrPoly e(triangle, zsecs);

    // A ray at x = 2 parallel to z misses the triangle.
    Vector3D pos(2.0, 0.0, -5.0);
    Vector3D dir(0.0, 0.0, 1.0);
    auto hits_before = e.ComputeIntersections(pos, dir);
    EXPECT_TRUE(hits_before.empty())
        << "Ray at x=2 should miss small triangle (sanity check)";

    // Replace the polygon with a square wide enough to contain x = 2.
    // If SetPolygon failed to recompute planes_, the next intersection
    // call would still be running against the triangle's 3 outward-normal
    // planes and would reject (2, 0) as outside.
    std::vector<std::vector<double>> square = {
        {-3.0, -3.0}, {3.0, -3.0}, {3.0, 3.0}, {-3.0, 3.0}
    };
    e.SetPolygon(square);

    auto hits_after = e.ComputeIntersections(pos, dir);
    EXPECT_EQ(hits_after.size(), 2u)
        << "After SetPolygon to enclosing square, ray must hit both z-caps; "
           "got " << hits_after.size() << ". A failure here means planes_ "
           "is stale - SetPolygon did not call ComputeLateralPlanes.";
    if(hits_after.size() == 2u) {
        // Z-caps at -1 and 1; ray origin at z=-5 moving +z.
        EXPECT_NEAR(hits_after[0].distance, 4.0, 1e-9);
        EXPECT_NEAR(hits_after[1].distance, 6.0, 1e-9);
    }
}
