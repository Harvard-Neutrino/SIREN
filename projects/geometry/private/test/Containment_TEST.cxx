// Containment validation tests
//
// For each geometry shape, compare the intersection-based IsInside() result
// against an independent analytical containment check. This validates that
// the ray intersection calculations produce correct enter/exit pairs.

#include <cmath>
#include <random>
#include <iostream>

#include <gtest/gtest.h>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
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

std::mt19937 rng_(12345);
std::uniform_real_distribution<double> uniform(-1.0, 1.0);
std::uniform_real_distribution<double> uniform01(0.0, 1.0);

Vector3D RandomPoint(double scale) {
    return Vector3D(uniform(rng_) * scale, uniform(rng_) * scale, uniform(rng_) * scale);
}

Vector3D RandomDirection() {
    double dx = uniform(rng_);
    double dy = uniform(rng_);
    double dz = uniform(rng_);
    double mag = std::sqrt(dx*dx + dy*dy + dz*dz);
    if(mag < 1e-12) return Vector3D(1, 0, 0);
    return Vector3D(dx/mag, dy/mag, dz/mag);
}

// Independent containment checks -- no intersection math involved

bool InsideBox(Vector3D const & p, double x, double y, double z) {
    return std::fabs(p.GetX()) < 0.5*x
        && std::fabs(p.GetY()) < 0.5*y
        && std::fabs(p.GetZ()) < 0.5*z;
}

bool InsideSphere(Vector3D const & p, double r_outer, double r_inner) {
    double r2 = p.GetX()*p.GetX() + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
    return r2 < r_outer*r_outer && r2 > r_inner*r_inner;
}

bool InsideCylinder(Vector3D const & p, double r_outer, double r_inner, double z) {
    double rxy2 = p.GetX()*p.GetX() + p.GetY()*p.GetY();
    return rxy2 < r_outer*r_outer
        && rxy2 > r_inner*r_inner
        && std::fabs(p.GetZ()) < 0.5*z;
}

bool InsideCone(Vector3D const & p, double rmin1, double rmax1, double rmin2, double rmax2, double z) {
    double hz = 0.5 * z;
    double pz = p.GetZ();
    if(pz <= -hz || pz >= hz) return false;
    // Interpolation fraction: 0 at -hz, 1 at +hz
    double frac = (pz + hz) / z;
    double rmax_at_z = rmax1 + (rmax2 - rmax1) * frac;
    double rmin_at_z = rmin1 + (rmin2 - rmin1) * frac;
    double rxy2 = p.GetX()*p.GetX() + p.GetY()*p.GetY();
    return rxy2 < rmax_at_z*rmax_at_z && rxy2 > rmin_at_z*rmin_at_z;
}

bool InsideTrd(Vector3D const & p, double dx1, double dx2, double dy1, double dy2, double dz) {
    double pz = p.GetZ();
    if(std::fabs(pz) >= dz) return false;
    double frac = (pz + dz) / (2.0 * dz);
    double xh = dx1 + (dx2 - dx1) * frac;
    double yh = dy1 + (dy2 - dy1) * frac;
    return std::fabs(p.GetX()) < xh && std::fabs(p.GetY()) < yh;
}

bool InsidePolycone(Vector3D const & p, std::vector<double> const & zp,
                    std::vector<double> const & rmin, std::vector<double> const & rmax) {
    double pz = p.GetZ();
    double rxy2 = p.GetX()*p.GetX() + p.GetY()*p.GetY();
    // Find the z-segment
    for(size_t i = 0; i + 1 < zp.size(); ++i) {
        if(pz > zp[i] && pz < zp[i+1]) {
            double frac = (pz - zp[i]) / (zp[i+1] - zp[i]);
            double rmax_z = rmax[i] + (rmax[i+1] - rmax[i]) * frac;
            double rmin_z = rmin[i] + (rmin[i+1] - rmin[i]) * frac;
            return rxy2 < rmax_z*rmax_z && rxy2 > rmin_z*rmin_z;
        }
    }
    return false;
}

// For Polyhedra: check if (x,y) is inside a regular polygon.
// The Polyhedra shape places vertices at (r*cos(angle), r*sin(angle)), so r is
// the circumradius (vertex distance from center). Each face edge connects two
// adjacent vertices; test that the point is on the interior side of every edge.
bool InsideRegularPolygon(double px, double py, int n, double r, double start_phi) {
    double dphi = 2.0 * M_PI / n;
    for(int k = 0; k < n; ++k) {
        double a0 = start_phi + k * dphi;
        double a1 = start_phi + (k + 1) * dphi;
        // Edge from vertex k to vertex k+1
        double vx0 = r * std::cos(a0), vy0 = r * std::sin(a0);
        double vx1 = r * std::cos(a1), vy1 = r * std::sin(a1);
        // Inward-pointing normal: rotate edge vector 90 degrees clockwise
        double ex = vx1 - vx0, ey = vy1 - vy0;
        double nx = ey, ny = -ex; // outward normal (pointing away from center)
        // Point is outside if on the outward side of any edge
        if(nx * (px - vx0) + ny * (py - vy0) > 0) return false;
    }
    return true;
}

bool InsidePolyhedra(Vector3D const & p, int n, double start_phi,
                     std::vector<double> const & zp,
                     std::vector<double> const & rmin, std::vector<double> const & rmax) {
    double pz = p.GetZ();
    for(size_t i = 0; i + 1 < zp.size(); ++i) {
        if(pz > zp[i] && pz < zp[i+1]) {
            double frac = (pz - zp[i]) / (zp[i+1] - zp[i]);
            double rmax_z = rmax[i] + (rmax[i+1] - rmax[i]) * frac;
            double rmin_z = rmin[i] + (rmin[i+1] - rmin[i]) * frac;
            bool in_outer = InsideRegularPolygon(p.GetX(), p.GetY(), n, rmax_z, start_phi);
            bool in_inner = (rmin_z > 0) && InsideRegularPolygon(p.GetX(), p.GetY(), n, rmin_z, start_phi);
            return in_outer && !in_inner;
        }
    }
    return false;
}

// Helper: test a shape with N random points, comparing IsInside against the
// independent check. Uses multiple random directions per point since IsInside
// result should be direction-independent for interior points.
void ValidateContainment(Geometry const & geo, Placement const & placement,
                         std::function<bool(Vector3D const &)> analytic_inside,
                         double scale, int N, std::string const & label) {
    int mismatches = 0;
    for(int i = 0; i < N; ++i) {
        Vector3D global_pos = RandomPoint(scale);
        // Transform to local coords for the analytic check
        Vector3D local_pos = placement.GlobalToLocalPosition(global_pos);
        bool analytic = analytic_inside(local_pos);

        // Test with several random directions -- containment shouldn't depend on direction
        // (except on the surface boundary, which we avoid with strict inequalities)
        for(int d = 0; d < 3; ++d) {
            Vector3D dir = RandomDirection();
            bool intersection_based = geo.IsInside(global_pos, dir);
            if(analytic != intersection_based) {
                mismatches++;
                if(mismatches <= 3) {
                    std::cerr << label << " mismatch #" << mismatches
                              << ": local=(" << local_pos.GetX() << "," << local_pos.GetY() << "," << local_pos.GetZ() << ")"
                              << " dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ() << ")"
                              << " analytic=" << analytic << " intersection=" << intersection_based
                              << std::endl;
                }
            }
        }
    }
    EXPECT_EQ(mismatches, 0) << label << ": " << mismatches << " mismatches out of " << N*3 << " checks";
}

} // anonymous namespace

// =========================================================================
// Box
// =========================================================================
TEST(Containment, Box) {
    double x = 10, y = 8, z = 6;
    Placement pl(Vector3D(5, -3, 2));
    Box box(pl, x, y, z);
    ValidateContainment(box, pl, [=](Vector3D const & p) {
        return InsideBox(p, x, y, z);
    }, 20, 10000, "Box");
}

TEST(Containment, BoxRotated) {
    double x = 10, y = 8, z = 6;
    Quaternion q(std::cos(0.3), std::sin(0.3)*0.577, std::sin(0.3)*0.577, std::sin(0.3)*0.577);
    Placement pl(Vector3D(1, 2, 3), q);
    Box box(pl, x, y, z);
    ValidateContainment(box, pl, [=](Vector3D const & p) {
        return InsideBox(p, x, y, z);
    }, 20, 10000, "BoxRotated");
}

// =========================================================================
// Sphere
// =========================================================================
TEST(Containment, Sphere) {
    double r = 5;
    Placement pl(Vector3D(-2, 1, 3));
    Sphere sphere(pl, r, 0);
    ValidateContainment(sphere, pl, [=](Vector3D const & p) {
        return InsideSphere(p, r, 0);
    }, 15, 10000, "Sphere");
}

TEST(Containment, SphereHollow) {
    double r_out = 8, r_in = 3;
    Placement pl(Vector3D(0, 0, 0));
    Sphere sphere(pl, r_out, r_in);
    ValidateContainment(sphere, pl, [=](Vector3D const & p) {
        return InsideSphere(p, r_out, r_in);
    }, 15, 10000, "SphereHollow");
}

// =========================================================================
// Cylinder
// =========================================================================
TEST(Containment, Cylinder) {
    double r = 5, z = 12;
    Placement pl(Vector3D(3, 0, -1));
    Cylinder cyl(pl, r, 0, z);
    ValidateContainment(cyl, pl, [=](Vector3D const & p) {
        return InsideCylinder(p, r, 0, z);
    }, 15, 10000, "Cylinder");
}

TEST(Containment, CylinderHollow) {
    double r_out = 6, r_in = 2, z = 10;
    Placement pl(Vector3D(0, 0, 0));
    Cylinder cyl(pl, r_out, r_in, z);
    ValidateContainment(cyl, pl, [=](Vector3D const & p) {
        return InsideCylinder(p, r_out, r_in, z);
    }, 15, 10000, "CylinderHollow");
}

// =========================================================================
// Cone
// =========================================================================
TEST(Containment, ConeSolid) {
    double rmin1 = 0, rmax1 = 8, rmin2 = 0, rmax2 = 3, z = 10;
    Placement pl(Vector3D(0, 0, 0));
    Cone cone(pl, rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, pl, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 15, 10000, "ConeSolid");
}

TEST(Containment, ConeHollow) {
    double rmin1 = 2, rmax1 = 6, rmin2 = 1, rmax2 = 4, z = 8;
    Placement pl(Vector3D(1, -1, 2));
    Cone cone(pl, rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, pl, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 15, 10000, "ConeHollow");
}

TEST(Containment, ConeCylinder) {
    // Degenerate: cone with equal radii = cylinder
    double r = 5, z = 10;
    Placement pl(Vector3D(0, 0, 0));
    Cone cone(pl, 0, r, 0, r, z);
    ValidateContainment(cone, pl, [=](Vector3D const & p) {
        return InsideCone(p, 0, r, 0, r, z);
    }, 15, 10000, "ConeCylinder");
}

// =========================================================================
// Trd
// =========================================================================
TEST(Containment, Trd) {
    double dx1 = 5, dx2 = 3, dy1 = 4, dy2 = 2, dz = 6;
    Placement pl(Vector3D(0, 0, 0));
    Trd trd(pl, dx1, dx2, dy1, dy2, dz);
    ValidateContainment(trd, pl, [=](Vector3D const & p) {
        return InsideTrd(p, dx1, dx2, dy1, dy2, dz);
    }, 15, 10000, "Trd");
}

TEST(Containment, TrdBox) {
    // Degenerate: trd with equal half-widths = box
    double dx = 4, dy = 3, dz = 5;
    Placement pl(Vector3D(2, -1, 0));
    Trd trd(pl, dx, dx, dy, dy, dz);
    ValidateContainment(trd, pl, [=](Vector3D const & p) {
        return InsideTrd(p, dx, dx, dy, dy, dz);
    }, 15, 10000, "TrdBox");
}

// =========================================================================
// Polycone
// =========================================================================
TEST(Containment, Polycone) {
    std::vector<double> zp = {-5, -2, 0, 3, 5};
    std::vector<double> rmin = {0, 0, 0, 0, 0};
    std::vector<double> rmax = {3, 5, 4, 6, 2};
    Placement pl(Vector3D(0, 0, 0));
    Polycone pc(pl, zp, rmin, rmax);
    ValidateContainment(pc, pl, [&](Vector3D const & p) {
        return InsidePolycone(p, zp, rmin, rmax);
    }, 10, 10000, "Polycone");
}

TEST(Containment, PolyconeHollow) {
    std::vector<double> zp = {-4, 0, 4};
    std::vector<double> rmin = {1, 2, 1};
    std::vector<double> rmax = {5, 6, 5};
    Placement pl(Vector3D(0, 0, 0));
    Polycone pc(pl, zp, rmin, rmax);
    ValidateContainment(pc, pl, [&](Vector3D const & p) {
        return InsidePolycone(p, zp, rmin, rmax);
    }, 10, 10000, "PolyconeHollow");
}

// =========================================================================
// Polyhedra
// =========================================================================
TEST(Containment, Polyhedra) {
    int n = 6;
    double start_phi = 0;
    std::vector<double> zp = {-5, 0, 5};
    std::vector<double> rmin = {0, 0, 0};
    std::vector<double> rmax = {4, 6, 4};
    Placement pl(Vector3D(0, 0, 0));
    Polyhedra ph(pl, n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, pl, [&](Vector3D const & p) {
        return InsidePolyhedra(p, n, start_phi, zp, rmin, rmax);
    }, 10, 10000, "Polyhedra");
}

TEST(Containment, PolyhedraSquare) {
    // 4-sided polyhedra = tapered square prism
    int n = 4;
    double start_phi = M_PI / 4.0; // rotated 45 degrees so faces are axis-aligned
    std::vector<double> zp = {-3, 3};
    std::vector<double> rmin = {0, 0};
    std::vector<double> rmax = {5, 5};
    Placement pl(Vector3D(0, 0, 0));
    Polyhedra ph(pl, n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, pl, [&](Vector3D const & p) {
        return InsidePolyhedra(p, n, start_phi, zp, rmin, rmax);
    }, 10, 10000, "PolyhedraSquare");
}

// =========================================================================
// BooleanGeometry
// =========================================================================
TEST(Containment, BooleanUnion) {
    // Two overlapping spheres
    auto left = Sphere(Placement(Vector3D(-2, 0, 0)), 5, 0).create();
    auto right = Sphere(Placement(Vector3D(2, 0, 0)), 5, 0).create();
    BooleanGeometry bg(BooleanOperation::UNION, left, right);
    Placement pl; // identity
    ValidateContainment(bg, pl, [](Vector3D const & p) {
        double r2_left = (p.GetX()+2)*(p.GetX()+2) + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
        double r2_right = (p.GetX()-2)*(p.GetX()-2) + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
        return r2_left < 25 || r2_right < 25;
    }, 10, 10000, "BooleanUnion");
}

TEST(Containment, BooleanSubtraction) {
    // Sphere with a cylindrical hole
    auto sphere = Sphere(5, 0).create();
    auto hole = Cylinder(2, 0, 12).create();
    BooleanGeometry bg(BooleanOperation::SUBTRACTION, sphere, hole);
    Placement pl;
    ValidateContainment(bg, pl, [](Vector3D const & p) {
        double r2 = p.GetX()*p.GetX() + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
        double rxy2 = p.GetX()*p.GetX() + p.GetY()*p.GetY();
        bool in_sphere = r2 < 25;
        bool in_cyl = rxy2 < 4 && std::fabs(p.GetZ()) < 6;
        return in_sphere && !in_cyl;
    }, 10, 10000, "BooleanSubtraction");
}

TEST(Containment, BooleanIntersection) {
    // Intersection of a sphere and a box
    auto sphere = Sphere(5, 0).create();
    auto box = Box(8, 8, 8).create();
    BooleanGeometry bg(BooleanOperation::INTERSECTION, sphere, box);
    Placement pl;
    ValidateContainment(bg, pl, [](Vector3D const & p) {
        double r2 = p.GetX()*p.GetX() + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
        bool in_sphere = r2 < 25;
        bool in_box = std::fabs(p.GetX()) < 4 && std::fabs(p.GetY()) < 4 && std::fabs(p.GetZ()) < 4;
        return in_sphere && in_box;
    }, 10, 10000, "BooleanIntersection");
}

// =========================================================================
// Pointed Cone (#21)
// =========================================================================
TEST(Containment, ConePointed) {
    double rmin1 = 0, rmax1 = 5, rmin2 = 0, rmax2 = 0, z = 10;
    Placement pl(Vector3D(0, 0, 0));
    Cone cone(pl, rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, pl, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 10, 10000, "ConePointed");
}

TEST(Containment, ConePointedReverse) {
    double rmin1 = 0, rmax1 = 0, rmin2 = 0, rmax2 = 5, z = 10;
    Placement pl(Vector3D(0, 0, 0));
    Cone cone(pl, rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, pl, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 10, 10000, "ConePointedReverse");
}

// =========================================================================
// ExtrPoly (#20)
// =========================================================================
TEST(Containment, ExtrPolySquare) {
    std::vector<std::vector<double>> polygon = {{-3,-3},{3,-3},{3,3},{-3,3}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off, 1.0),
        ExtrPoly::ZSection(4, off, 1.0)
    };
    Placement pl(Vector3D(0, 0, 0));
    ExtrPoly ep(pl, polygon, zsecs);
    ValidateContainment(ep, pl, [](Vector3D const & p) {
        return std::fabs(p.GetX()) < 3 && std::fabs(p.GetY()) < 3 && std::fabs(p.GetZ()) < 4;
    }, 8, 10000, "ExtrPolySquare");
}

TEST(Containment, ExtrPolyTapered) {
    // Tapered: scale 1.5 at bottom, 0.5 at top
    std::vector<std::vector<double>> polygon = {{-2,-2},{2,-2},{2,2},{-2,2}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-3, off, 1.5),
        ExtrPoly::ZSection(3, off, 0.5)
    };
    Placement pl(Vector3D(0, 0, 0));
    ExtrPoly ep(pl, polygon, zsecs);
    ValidateContainment(ep, pl, [](Vector3D const & p) {
        double frac = (p.GetZ() + 3.0) / 6.0;
        double s = 1.5 + (0.5 - 1.5) * frac;
        return std::fabs(p.GetX()) < 2.0 * s
            && std::fabs(p.GetY()) < 2.0 * s
            && p.GetZ() > -3.0 && p.GetZ() < 3.0;
    }, 8, 10000, "ExtrPolyTapered");
}

// =========================================================================
// Rotated shapes (#22)
// =========================================================================
TEST(Containment, ConeRotated) {
    double rmin1 = 0, rmax1 = 6, rmin2 = 0, rmax2 = 3, z = 10;
    Quaternion q(std::cos(0.4), std::sin(0.4)*0.577, std::sin(0.4)*0.577, std::sin(0.4)*0.577);
    Placement pl(Vector3D(2, -1, 3), q);
    Cone cone(pl, rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, pl, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 20, 10000, "ConeRotated");
}

TEST(Containment, TrdRotated) {
    double dx1 = 5, dx2 = 3, dy1 = 4, dy2 = 2, dz = 6;
    Quaternion q(std::cos(0.3), std::sin(0.3)*0.577, std::sin(0.3)*0.577, std::sin(0.3)*0.577);
    Placement pl(Vector3D(-1, 2, 0), q);
    Trd trd(pl, dx1, dx2, dy1, dy2, dz);
    ValidateContainment(trd, pl, [=](Vector3D const & p) {
        return InsideTrd(p, dx1, dx2, dy1, dy2, dz);
    }, 20, 10000, "TrdRotated");
}

TEST(Containment, PolyconeRotated) {
    std::vector<double> zp = {-4, 0, 4};
    std::vector<double> rmin = {0, 0, 0};
    std::vector<double> rmax = {3, 5, 3};
    Quaternion q(std::cos(0.5), std::sin(0.5)*0.577, std::sin(0.5)*0.577, std::sin(0.5)*0.577);
    Placement pl(Vector3D(1, 0, -2), q);
    Polycone pc(pl, zp, rmin, rmax);
    ValidateContainment(pc, pl, [&](Vector3D const & p) {
        return InsidePolycone(p, zp, rmin, rmax);
    }, 12, 10000, "PolyconeRotated");
}

TEST(Containment, PolyhedraRotated) {
    int n = 6;
    double start_phi = 0;
    std::vector<double> zp = {-3, 3};
    std::vector<double> rmin = {0, 0};
    std::vector<double> rmax = {4, 4};
    Quaternion q(std::cos(0.6), std::sin(0.6)*0.577, std::sin(0.6)*0.577, std::sin(0.6)*0.577);
    Placement pl(Vector3D(0, -2, 1), q);
    Polyhedra ph(pl, n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, pl, [&](Vector3D const & p) {
        return InsidePolyhedra(p, n, start_phi, zp, rmin, rmax);
    }, 12, 10000, "PolyhedraRotated");
}

// =========================================================================
// Polyhedra with zero-radius apex (pyramid shapes)
// =========================================================================
TEST(Containment, PolyhedraPyramid) {
    // 4-sided polyhedra that tapers to a point at the top: rmax=0 at z=5
    int n = 4;
    double start_phi = M_PI / 4.0;
    std::vector<double> zp = {-3, 5};
    std::vector<double> rmin = {0, 0};
    std::vector<double> rmax = {5, 0};
    Placement pl(Vector3D(0, 0, 0));
    Polyhedra ph(pl, n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, pl, [&](Vector3D const & p) {
        return InsidePolyhedra(p, n, start_phi, zp, rmin, rmax);
    }, 10, 10000, "PolyhedraPyramid");
}

TEST(Containment, PolyhedraPyramidReverse) {
    // 6-sided polyhedra that tapers to a point at the bottom: rmax=0 at z=-4
    int n = 6;
    double start_phi = 0;
    std::vector<double> zp = {-4, 4};
    std::vector<double> rmin = {0, 0};
    std::vector<double> rmax = {0, 6};
    Placement pl(Vector3D(0, 0, 0));
    Polyhedra ph(pl, n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, pl, [&](Vector3D const & p) {
        return InsidePolyhedra(p, n, start_phi, zp, rmin, rmax);
    }, 10, 10000, "PolyhedraPyramidReverse");
}

// =========================================================================
// AABB cache invalidation after SetPlacement
// =========================================================================
TEST(Containment, AABBCacheInvalidation) {
    Box box(10, 8, 6);
    // Initial placement at origin
    AABB bb1 = box.GetWorldBoundingBox();
    EXPECT_NEAR(bb1.max_corner.GetX(), 5.0, 1e-9);

    // Move placement to (100, 0, 0)
    box.SetPlacement(Placement(Vector3D(100, 0, 0)));
    AABB bb2 = box.GetWorldBoundingBox();
    // The AABB should reflect the new placement, not the cached old one
    EXPECT_NEAR(bb2.max_corner.GetX(), 105.0, 1e-9)
        << "AABB cache should be invalidated after SetPlacement";
    EXPECT_NEAR(bb2.min_corner.GetX(), 95.0, 1e-9);

    // Call again to verify cache hit returns correct value
    AABB bb3 = box.GetWorldBoundingBox();
    EXPECT_NEAR(bb3.max_corner.GetX(), 105.0, 1e-9);
}

// =========================================================================
// BooleanGeometry with ray from inside (#3 verification)
// =========================================================================
TEST(Containment, BooleanSubtractionFromInside) {
    // Sphere with a box cut out, test that rays from inside work
    auto sphere = Sphere(5, 0).create();
    auto box = Box(4, 4, 4).create();
    BooleanGeometry bg(BooleanOperation::SUBTRACTION, sphere, box);
    Placement pl;
    ValidateContainment(bg, pl, [](Vector3D const & p) {
        double r2 = p.GetX()*p.GetX() + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
        bool in_sphere = r2 < 25;
        bool in_box = std::fabs(p.GetX()) < 2 && std::fabs(p.GetY()) < 2 && std::fabs(p.GetZ()) < 2;
        return in_sphere && !in_box;
    }, 10, 10000, "BooleanSubtractionFromInside");
}
