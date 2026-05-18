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

bool InsideSpherePartial(Vector3D const & p, double r_outer, double r_inner,
                         double start_phi, double delta_phi,
                         double start_theta, double delta_theta) {
    double r2 = p.GetX()*p.GetX() + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
    if(r2 >= r_outer*r_outer || r2 <= r_inner*r_inner) return false;
    // Phi check
    double phi = std::atan2(p.GetY(), p.GetX());
    if(phi < 0) phi += 2.0 * M_PI;
    double sp = start_phi;
    while(sp < 0) sp += 2.0 * M_PI;
    while(sp >= 2.0 * M_PI) sp -= 2.0 * M_PI;
    double ep = sp + delta_phi;
    bool phi_ok;
    if(ep <= 2.0 * M_PI + 1e-9) {
        phi_ok = phi >= sp && phi <= ep;
    } else {
        phi_ok = phi >= sp || phi <= (ep - 2.0 * M_PI);
    }
    if(!phi_ok) return false;
    // Theta check
    double r = std::sqrt(r2);
    double ct = p.GetZ() / r;
    if(ct > 1.0) ct = 1.0;
    if(ct < -1.0) ct = -1.0;
    double theta = std::acos(ct);
    return theta >= start_theta && theta <= start_theta + delta_theta;
}

bool InsidePhiRange(double x, double y, double start_phi, double delta_phi) {
    double phi = std::fmod(std::atan2(y, x), 2.0 * M_PI);
    if(phi < 0) phi += 2.0 * M_PI;
    double sp = std::fmod(start_phi, 2.0 * M_PI);
    if(sp < 0) sp += 2.0 * M_PI;
    double ep = sp + delta_phi;
    if(ep <= 2.0 * M_PI + 1e-9) {
        return phi >= sp - 1e-9 && phi <= ep + 1e-9;
    } else {
        return phi >= sp - 1e-9 || phi <= std::fmod(ep, 2.0*M_PI) + 1e-9;
    }
}

bool InsideTorusPartial(Vector3D const & p, double R, double rmax, double rmin,
                        double start_phi, double delta_phi) {
    double rxy = std::sqrt(p.GetX()*p.GetX() + p.GetY()*p.GetY());
    double d2 = (rxy - R)*(rxy - R) + p.GetZ()*p.GetZ();
    if(d2 >= rmax*rmax) return false;
    if(rmin > 0 && d2 <= rmin*rmin) return false;
    // Phi check
    double phi = std::atan2(p.GetY(), p.GetX());
    if(phi < 0) phi += 2.0 * M_PI;
    double sp = start_phi;
    while(sp < 0) sp += 2.0 * M_PI;
    while(sp >= 2.0 * M_PI) sp -= 2.0 * M_PI;
    double ep = sp + delta_phi;
    if(ep <= 2.0 * M_PI + 1e-9) {
        return phi >= sp && phi <= ep;
    } else {
        return phi >= sp || phi <= (ep - 2.0 * M_PI);
    }
}

bool InsideTorus(Vector3D const & p, double R, double rmax, double rmin) {
    double rxy = std::sqrt(p.GetX()*p.GetX() + p.GetY()*p.GetY());
    double d2 = (rxy - R)*(rxy - R) + p.GetZ()*p.GetZ();
    return d2 < rmax*rmax && (rmin <= 0 || d2 > rmin*rmin);
}

// Helper: test a shape with N random points, comparing IsInside against the
// independent check. Uses multiple random directions per point since IsInside
// result should be direction-independent for interior points.
void ValidateContainment(Geometry const & geo,
                         std::function<bool(Vector3D const &)> analytic_inside,
                         double scale, int N, std::string const & label) {
    Placement placement = geo.GetPlacement();
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

bool InsideEllipticalTube(Vector3D const & p, double dx, double dy, double dz) {
    double ex = p.GetX() / dx;
    double ey = p.GetY() / dy;
    return ex*ex + ey*ey < 1.0 && std::fabs(p.GetZ()) < dz;
}

bool InsideCutTube(Vector3D const & p, double rmin, double rmax, double dz,
                   Vector3D const & low_norm, Vector3D const & high_norm) {
    double rxy2 = p.GetX()*p.GetX() + p.GetY()*p.GetY();
    if(rxy2 >= rmax*rmax || rxy2 <= rmin*rmin) return false;
    // Check between cut planes (normals point outward)
    double low_dot = low_norm.GetX()*p.GetX() + low_norm.GetY()*p.GetY() + low_norm.GetZ()*(p.GetZ() + dz);
    double high_dot = high_norm.GetX()*p.GetX() + high_norm.GetY()*p.GetY() + high_norm.GetZ()*(p.GetZ() - dz);
    return low_dot < 0 && high_dot < 0;
}

bool InsideTrap(Vector3D const & p, double dz, double theta, double phi,
                double dy1, double dx1, double dx2, double alpha1,
                double dy2, double dx3, double dx4, double alpha2) {
    // Compute the 8 vertices and 6 plane equations, then test the point
    // against all planes. This matches the geometry's plane-equation approach.
    double ta1 = std::tan(alpha1), ta2 = std::tan(alpha2);
    double tt = std::tan(theta);
    double cp = std::cos(phi), sp = std::sin(phi);
    double xc_lo = -dz * tt * cp, yc_lo = -dz * tt * sp;
    double xc_hi = +dz * tt * cp, yc_hi = +dz * tt * sp;

    double v[8][3] = {
        {xc_lo - dy1*ta1 - dx1, yc_lo - dy1, -dz},
        {xc_lo - dy1*ta1 + dx1, yc_lo - dy1, -dz},
        {xc_lo + dy1*ta1 + dx2, yc_lo + dy1, -dz},
        {xc_lo + dy1*ta1 - dx2, yc_lo + dy1, -dz},
        {xc_hi - dy2*ta2 - dx3, yc_hi - dy2, +dz},
        {xc_hi - dy2*ta2 + dx3, yc_hi - dy2, +dz},
        {xc_hi + dy2*ta2 + dx4, yc_hi + dy2, +dz},
        {xc_hi + dy2*ta2 - dx4, yc_hi + dy2, +dz}
    };

    // 6 faces defined by vertex indices (outward winding).
    // The first 3 vertices of each face determine the plane equation,
    // so their triangulation must match Trap::ComputePlanes exactly.
    int faces[6][4] = {
        {0, 1, 2, 3}, // bottom (z=-dz)
        {7, 6, 5, 4}, // top (z=+dz)
        {0, 1, 5, 4}, // front (-y)
        {3, 2, 6, 7}, // back (+y)
        {0, 3, 7, 4}, // left (-x)
        {1, 2, 6, 5}  // right (+x)
    };

    // Centroid for outward normal orientation
    double cx = 0, cy = 0, cz = 0;
    for(int i = 0; i < 8; ++i) { cx += v[i][0]; cy += v[i][1]; cz += v[i][2]; }
    cx /= 8; cy /= 8; cz /= 8;

    double px = p.GetX(), py = p.GetY(), pz = p.GetZ();
    for(int f = 0; f < 6; ++f) {
        int i0 = faces[f][0], i1 = faces[f][1], i2 = faces[f][2];
        double ex1 = v[i1][0]-v[i0][0], ey1 = v[i1][1]-v[i0][1], ez1 = v[i1][2]-v[i0][2];
        double ex2 = v[i2][0]-v[i0][0], ey2 = v[i2][1]-v[i0][1], ez2 = v[i2][2]-v[i0][2];
        double nx = ey1*ez2 - ez1*ey2;
        double ny = ez1*ex2 - ex1*ez2;
        double nz = ex1*ey2 - ey1*ex2;
        // Ensure outward (away from centroid)
        double cdot = nx*(v[i0][0]-cx) + ny*(v[i0][1]-cy) + nz*(v[i0][2]-cz);
        if(cdot < 0) { nx = -nx; ny = -ny; nz = -nz; }
        double nmag = std::sqrt(nx*nx + ny*ny + nz*nz);
        if(nmag > 0) { nx /= nmag; ny /= nmag; nz /= nmag; }
        double d = -(nx*v[i0][0] + ny*v[i0][1] + nz*v[i0][2]);
        double val = nx*px + ny*py + nz*pz + d;
        // Strictly inside means on the inner side of every face. (The
        // normal is unit, so val is a true signed distance; the prior
        // val > -0.01 skin eroded the oracle inward by 1cm and would
        // have hidden an undersized production Trap at its faces.)
        if(val > 0.0) return false;
    }
    return true;
}

bool InsideEllipsoid(Vector3D const & p, double ax, double by, double cz,
                     double zcut1, double zcut2) {
    double ex = p.GetX() / ax;
    double ey = p.GetY() / by;
    double ez = p.GetZ() / cz;
    return ex*ex + ey*ey + ez*ez < 1.0
        && p.GetZ() > zcut1 && p.GetZ() < zcut2;
}

bool InsidePara(Vector3D const & p, double dx, double dy, double dz,
                double alpha, double theta, double phi) {
    // Reverse the shear to get "box" coordinates
    double e3x = std::sin(theta) * std::cos(phi);
    double e3y = std::sin(theta) * std::sin(phi);
    double e3z = std::cos(theta);
    // Solve for w: p.z = w * dz * cos(theta)
    double w = p.GetZ() / (dz * e3z);
    if(std::fabs(w) >= 1.0) return false;
    // Solve for v: p.y - w*dz*e3y = v * dy
    double v = (p.GetY() - w * dz * e3y) / dy;
    if(std::fabs(v) >= 1.0) return false;
    // Solve for u: p.x - w*dz*e3x - v*dy*tan(alpha) = u * dx
    double u = (p.GetX() - w * dz * e3x - v * dy * std::tan(alpha)) / dx;
    if(std::fabs(u) >= 1.0) return false;
    return true;
}

} // anonymous namespace

// =========================================================================
// Box
// =========================================================================
TEST(Containment, Box) {
    double x = 10, y = 8, z = 6;
    Placement pl(Vector3D(5, -3, 2));
    Box box(pl, x, y, z);
    ValidateContainment(box, [=](Vector3D const & p) {
        return InsideBox(p, x, y, z);
    }, 20, 10000, "Box");
}

TEST(Containment, BoxRotated) {
    double x = 10, y = 8, z = 6;
    Quaternion q(std::cos(0.3), std::sin(0.3)*0.577, std::sin(0.3)*0.577, std::sin(0.3)*0.577);
    Placement pl(Vector3D(1, 2, 3), q);
    Box box(pl, x, y, z);
    ValidateContainment(box, [=](Vector3D const & p) {
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
    ValidateContainment(sphere, [=](Vector3D const & p) {
        return InsideSphere(p, r, 0);
    }, 15, 10000, "Sphere");
}

TEST(Containment, SphereHollow) {
    double r_out = 8, r_in = 3;
    Sphere sphere(r_out, r_in);
    ValidateContainment(sphere, [=](Vector3D const & p) {
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
    ValidateContainment(cyl, [=](Vector3D const & p) {
        return InsideCylinder(p, r, 0, z);
    }, 15, 10000, "Cylinder");
}

TEST(Containment, CylinderHollow) {
    double r_out = 6, r_in = 2, z = 10;
    Cylinder cyl(r_out, r_in, z);
    ValidateContainment(cyl, [=](Vector3D const & p) {
        return InsideCylinder(p, r_out, r_in, z);
    }, 15, 10000, "CylinderHollow");
}

// =========================================================================
// Cone
// =========================================================================
TEST(Containment, ConeSolid) {
    double rmin1 = 0, rmax1 = 8, rmin2 = 0, rmax2 = 3, z = 10;
    Cone cone(rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 15, 10000, "ConeSolid");
}

TEST(Containment, ConeHollow) {
    double rmin1 = 2, rmax1 = 6, rmin2 = 1, rmax2 = 4, z = 8;
    Placement pl(Vector3D(1, -1, 2));
    Cone cone(pl, rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 15, 10000, "ConeHollow");
}

TEST(Containment, ConeCylinder) {
    // Degenerate: cone with equal radii = cylinder
    double r = 5, z = 10;
    Cone cone(0, r, 0, r, z);
    ValidateContainment(cone, [=](Vector3D const & p) {
        return InsideCone(p, 0, r, 0, r, z);
    }, 15, 10000, "ConeCylinder");
}

// =========================================================================
// Trd
// =========================================================================
TEST(Containment, Trd) {
    double dx1 = 5, dx2 = 3, dy1 = 4, dy2 = 2, dz = 6;
    Trd trd(dx1, dx2, dy1, dy2, dz);
    ValidateContainment(trd, [=](Vector3D const & p) {
        return InsideTrd(p, dx1, dx2, dy1, dy2, dz);
    }, 15, 10000, "Trd");
}

TEST(Containment, TrdBox) {
    // Degenerate: trd with equal half-widths = box
    double dx = 4, dy = 3, dz = 5;
    Placement pl(Vector3D(2, -1, 0));
    Trd trd(pl, dx, dx, dy, dy, dz);
    ValidateContainment(trd, [=](Vector3D const & p) {
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
    Polycone pc(zp, rmin, rmax);
    ValidateContainment(pc, [&](Vector3D const & p) {
        return InsidePolycone(p, zp, rmin, rmax);
    }, 10, 10000, "Polycone");
}

TEST(Containment, PolyconeHollow) {
    std::vector<double> zp = {-4, 0, 4};
    std::vector<double> rmin = {1, 2, 1};
    std::vector<double> rmax = {5, 6, 5};
    Polycone pc(zp, rmin, rmax);
    ValidateContainment(pc, [&](Vector3D const & p) {
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
    Polyhedra ph(n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, [&](Vector3D const & p) {
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
    Polyhedra ph(n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, [&](Vector3D const & p) {
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
    ValidateContainment(bg, [](Vector3D const & p) {
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
    ValidateContainment(bg, [](Vector3D const & p) {
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
    ValidateContainment(bg, [](Vector3D const & p) {
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
    Cone cone(rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 10, 10000, "ConePointed");
}

TEST(Containment, ConePointedReverse) {
    double rmin1 = 0, rmax1 = 0, rmin2 = 0, rmax2 = 5, z = 10;
    Cone cone(rmin1, rmax1, rmin2, rmax2, z);
    ValidateContainment(cone, [=](Vector3D const & p) {
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
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
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
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
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
    ValidateContainment(cone, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z);
    }, 20, 10000, "ConeRotated");
}

TEST(Containment, TrdRotated) {
    double dx1 = 5, dx2 = 3, dy1 = 4, dy2 = 2, dz = 6;
    Quaternion q(std::cos(0.3), std::sin(0.3)*0.577, std::sin(0.3)*0.577, std::sin(0.3)*0.577);
    Placement pl(Vector3D(-1, 2, 0), q);
    Trd trd(pl, dx1, dx2, dy1, dy2, dz);
    ValidateContainment(trd, [=](Vector3D const & p) {
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
    ValidateContainment(pc, [&](Vector3D const & p) {
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
    ValidateContainment(ph, [&](Vector3D const & p) {
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
    Polyhedra ph(n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, [&](Vector3D const & p) {
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
    Polyhedra ph(n, start_phi, zp, rmin, rmax);
    ValidateContainment(ph, [&](Vector3D const & p) {
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
    ValidateContainment(bg, [](Vector3D const & p) {
        double r2 = p.GetX()*p.GetX() + p.GetY()*p.GetY() + p.GetZ()*p.GetZ();
        bool in_sphere = r2 < 25;
        bool in_box = std::fabs(p.GetX()) < 2 && std::fabs(p.GetY()) < 2 && std::fabs(p.GetZ()) < 2;
        return in_sphere && !in_box;
    }, 10, 10000, "BooleanSubtractionFromInside");
}

// =========================================================================
// Polycone with step-change in radius (fix #8)
// =========================================================================
TEST(Containment, PolyconeStepChange) {
    // z-planes with consecutive equal z to represent a radius step:
    // cylinder r=0.03 from z=-0.05 to z=0.05, then steps to r=0.05 from z=0.05 to z=0.15
    std::vector<double> zp = {-0.05, 0.05, 0.05, 0.15};
    std::vector<double> rmin = {0, 0, 0, 0};
    std::vector<double> rmax = {0.03, 0.03, 0.05, 0.05};
    Polycone pc(zp, rmin, rmax);
    ValidateContainment(pc, [&](Vector3D const & p) {
        double pz = p.GetZ();
        double rxy2 = p.GetX()*p.GetX() + p.GetY()*p.GetY();
        // First segment: z in (-0.05, 0.05), r < 0.03
        if(pz > -0.05 && pz < 0.05) {
            return rxy2 < 0.03*0.03;
        }
        // Second segment: z in (0.05, 0.15), r < 0.05
        if(pz > 0.05 && pz < 0.15) {
            return rxy2 < 0.05*0.05;
        }
        return false;
    }, 0.2, 10000, "PolyconeStepChange");
}

// =========================================================================
// Hollow cylinder with rotation (identified gap)
// =========================================================================
TEST(Containment, CylinderHollowRotated) {
    double r_out = 6, r_in = 2, z = 10;
    Quaternion q(std::cos(0.4), std::sin(0.4)*0.577, std::sin(0.4)*0.577, std::sin(0.4)*0.577);
    Placement pl(Vector3D(1, -1, 2), q);
    Cylinder cyl(pl, r_out, r_in, z);
    ValidateContainment(cyl, [=](Vector3D const & p) {
        return InsideCylinder(p, r_out, r_in, z);
    }, 20, 10000, "CylinderHollowRotated");
}

// =========================================================================
// ExtrPoly with non-zero offset in z-sections (identified gap)
// =========================================================================
TEST(Containment, ExtrPolyWithOffset) {
    std::vector<std::vector<double>> polygon = {{-2,-2},{2,-2},{2,2},{-2,2}};
    double off_bot[2] = {0, 0};
    double off_top[2] = {1, 1};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-3, off_bot, 1.0),
        ExtrPoly::ZSection(3, off_top, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
        // At height z, the polygon center is offset by interpolated (ox, oy)
        double frac = (p.GetZ() + 3.0) / 6.0;
        double ox = 1.0 * frac;
        double oy = 1.0 * frac;
        return std::fabs(p.GetX() - ox) < 2.0
            && std::fabs(p.GetY() - oy) < 2.0
            && p.GetZ() > -3.0 && p.GetZ() < 3.0;
    }, 8, 10000, "ExtrPolyWithOffset");
}

// =========================================================================
// ExtrPoly with 3+ z-sections (identified gap)
// =========================================================================
TEST(Containment, ExtrPolyMultiSection) {
    std::vector<std::vector<double>> polygon = {{-2,-2},{2,-2},{2,2},{-2,2}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off, 1.0),
        ExtrPoly::ZSection(0, off, 1.5),
        ExtrPoly::ZSection(4, off, 0.8)
    };
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
        double pz = p.GetZ();
        if(pz <= -4.0 || pz >= 4.0) return false;
        double s;
        if(pz < 0.0) {
            // Between z=-4 (scale=1.0) and z=0 (scale=1.5)
            double frac = (pz + 4.0) / 4.0;
            s = 1.0 + (1.5 - 1.0) * frac;
        } else {
            // Between z=0 (scale=1.5) and z=4 (scale=0.8)
            double frac = pz / 4.0;
            s = 1.5 + (0.8 - 1.5) * frac;
        }
        return std::fabs(p.GetX()) < 2.0 * s
            && std::fabs(p.GetY()) < 2.0 * s;
    }, 8, 10000, "ExtrPolyMultiSection");
}

// =========================================================================
// Cone with zero height throws (fix #6)
// =========================================================================
TEST(Containment, ConeZeroHeightThrows) {
    EXPECT_THROW(Cone(0, 5, 0, 3, 0), std::invalid_argument);
}

// =========================================================================
// BooleanGeometry less-than null safety (fix #2)
// =========================================================================
TEST(Containment, BooleanLessNullSafety) {
    auto bg = BooleanGeometry(BooleanOperation::UNION, Sphere(5, 0).create(), Box(4, 4, 4).create());
    Box box(10, 10, 10);
    EXPECT_NO_THROW(bool result = bg < box; (void)result);
    EXPECT_NO_THROW(bool result = box < bg; (void)result);
}

// =========================================================================
// Fix 3: Containment at Polycone/Polyhedra internal z-boundaries
// Points exactly at internal z-boundaries must still have correct
// containment determination via horizontal rays.
// =========================================================================
TEST(Containment, PolyconeAtZBoundary) {
    // Polycone with internal boundary at z=0 where radius changes
    std::vector<double> zp = {-5, 0, 5};
    std::vector<double> rmin = {0, 0, 0};
    std::vector<double> rmax = {3, 5, 3};
    Polycone pc(zp, rmin, rmax);

    // Point inside at z=0 boundary: (1, 0, 0), rmax=5, r=1 is inside.
    // A horizontal ray at the exact internal z-boundary DOES pierce the
    // barrel (it crosses the z=0 disc), so containment must be correct for
    // a horizontal probe direction, not only off-axis ones.
    Vector3D inside_at_boundary(1, 0, 0);
    EXPECT_TRUE(pc.IsInside(inside_at_boundary, Vector3D(1, 0, 0)))
        << "Point (1,0,0) inside polycone at z=0, horizontal ray";
    EXPECT_TRUE(pc.IsInside(inside_at_boundary, Vector3D(1, 1, 1)))
        << "Point (1,0,0) inside polycone at z=0, diagonal ray";
    EXPECT_TRUE(pc.IsInside(inside_at_boundary, Vector3D(0, 0, 1)))
        << "Point (1,0,0) inside polycone at z=0, z-direction ray";

    // Point outside at z=0 boundary: (6, 0, 0), r=6 > rmax=5
    Vector3D outside_at_boundary(6, 0, 0);
    EXPECT_FALSE(pc.IsInside(outside_at_boundary, Vector3D(1, 0, 0)))
        << "Point (6,0,0) outside polycone at z=0, horizontal ray";
    EXPECT_FALSE(pc.IsInside(outside_at_boundary, Vector3D(0, 0, 1)))
        << "Point (6,0,0) outside polycone at z=0, z-direction ray";
}

TEST(Containment, PolyhedraAtZBoundary) {
    int n = 6;
    double start_phi = 0;
    std::vector<double> zp = {-5, 0, 5};
    std::vector<double> rmin = {0, 0, 0};
    std::vector<double> rmax = {3, 5, 3};
    Polyhedra ph(n, start_phi, zp, rmin, rmax);

    // Point inside at z=0 boundary with horizontal direction
    Vector3D inside_at_boundary(1, 0, 0);
    EXPECT_TRUE(ph.IsInside(inside_at_boundary, Vector3D(1, 0, 0)))
        << "Point (1,0,0) should be inside polyhedra at z=0 boundary";
    EXPECT_TRUE(ph.IsInside(inside_at_boundary, Vector3D(0, 1, 0)))
        << "Point (1,0,0) should be inside polyhedra at z=0 with y-direction";
}

// =========================================================================
// Torus
// =========================================================================

TEST(Containment, Torus) {
    Torus torus(10, 3, 0);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorus(p, 10, 3, 0);
    }, 15, 10000, "Torus");
}

TEST(Containment, TorusHollow) {
    Torus torus(10, 3, 1);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorus(p, 10, 3, 1);
    }, 15, 10000, "TorusHollow");
}

TEST(Containment, TorusRotated) {
    Quaternion q(std::cos(0.4), std::sin(0.4)*0.577, std::sin(0.4)*0.577, std::sin(0.4)*0.577);
    Placement pl(Vector3D(1, -1, 2), q);
    Torus torus(pl, 8, 2, 0);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorus(p, 8, 2, 0);
    }, 15, 10000, "TorusRotated");
}

TEST(Containment, TorusValidation) {
    // Major radius must be positive
    EXPECT_THROW(Torus(0, 3, 0), std::invalid_argument);
    // Outer tube radius must be positive
    EXPECT_THROW(Torus(10, 0, 0), std::invalid_argument);
    // Inner must be less than outer
    EXPECT_THROW(Torus(10, 3, 3), std::invalid_argument);
    EXPECT_THROW(Torus(10, 3, 4), std::invalid_argument);
    // Valid configurations
    EXPECT_NO_THROW(Torus(10, 3, 0));
    EXPECT_NO_THROW(Torus(10, 3, 2));
}

TEST(Containment, TorusSelfIntersecting) {
    // Self-intersecting torus: rtor < rmax (tube overlaps through center)
    Torus torus(3, 5, 0);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorus(p, 3, 5, 0);
    }, 3 + 5 + 1, 10000, "TorusSelfIntersecting");
}

TEST(Containment, TorusSelfIntersectingHollow) {
    Torus torus(3, 5, 2);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorus(p, 3, 5, 2);
    }, 3 + 5 + 1, 10000, "TorusSelfIntersectingHollow");
}

// =========================================================================
// Partial angular extent: Sphere
// =========================================================================

TEST(Containment, SphereThetaCut) {
    // Upper hemisphere: theta from 0 to pi/2
    double r = 5;
    double st = 0, dt = M_PI / 2.0;
    Sphere sphere(r, 0, 0, 2*M_PI, st, dt);
    ValidateContainment(sphere, [=](Vector3D const & p) {
        return InsideSpherePartial(p, r, 0, 0, 2*M_PI, st, dt);
    }, 8, 10000, "SphereThetaCut");
}

TEST(Containment, SpherePhiCut) {
    // Half sphere in phi: phi from 0 to pi
    double r = 5;
    double sp = 0, dp = M_PI;
    Sphere sphere(r, 0, sp, dp, 0, M_PI);
    ValidateContainment(sphere, [=](Vector3D const & p) {
        return InsideSpherePartial(p, r, 0, sp, dp, 0, M_PI);
    }, 8, 10000, "SpherePhiCut");
}

TEST(Containment, SpherePhiThetaCut) {
    // Quarter sphere: phi [0,pi], theta [0,pi/2]
    double r = 5;
    double sp = 0, dp = M_PI, st = 0, dt = M_PI / 2.0;
    Sphere sphere(r, 0, sp, dp, st, dt);
    ValidateContainment(sphere, [=](Vector3D const & p) {
        return InsideSpherePartial(p, r, 0, sp, dp, st, dt);
    }, 8, 10000, "SpherePhiThetaCut");
}

TEST(Containment, SphereHollowThetaCut) {
    // Hollow upper hemisphere
    double r_out = 5, r_in = 2;
    double st = 0, dt = M_PI / 2.0;
    Sphere sphere(r_out, r_in, 0, 2*M_PI, st, dt);
    ValidateContainment(sphere, [=](Vector3D const & p) {
        return InsideSpherePartial(p, r_out, r_in, 0, 2*M_PI, st, dt);
    }, 8, 10000, "SphereHollowThetaCut");
}

// =========================================================================
// Partial angular extent: Torus
// =========================================================================

TEST(Containment, TorusPhiCut) {
    // Quarter torus: phi from 0 to pi/2
    double R = 10, r = 3;
    double sp = 0, dp = M_PI / 2.0;
    Torus torus(R, r, 0, sp, dp);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorusPartial(p, R, r, 0, sp, dp);
    }, 15, 10000, "TorusPhiCut");
}

TEST(Containment, TorusPhiCutHalf) {
    // Half torus: phi from 0 to pi
    double R = 10, r = 3;
    double sp = 0, dp = M_PI;
    Torus torus(R, r, 0, sp, dp);
    ValidateContainment(torus, [=](Vector3D const & p) {
        return InsideTorusPartial(p, R, r, 0, sp, dp);
    }, 15, 10000, "TorusPhiCutHalf");
}

// =========================================================================
// Partial angular extent: Cylinder
// =========================================================================

TEST(Containment, CylinderPhiCut) {
    // Quarter cylinder: phi from 0 to pi/2
    double rmax = 5, rmin = 2, z = 10;
    double sp = 0, dp = M_PI / 2.0;
    Cylinder cyl(rmax, rmin, z, sp, dp);
    ValidateContainment(cyl, [=](Vector3D const & p) {
        return InsideCylinder(p, rmax, rmin, z) && InsidePhiRange(p.GetX(), p.GetY(), sp, dp);
    }, 10, 10000, "CylinderPhiCut");
}

// =========================================================================
// Partial angular extent: Cone
// =========================================================================

TEST(Containment, ConePhiCut) {
    // Quarter cone: phi from 0 to pi/2
    double rmin1 = 1, rmax1 = 6, rmin2 = 0, rmax2 = 3, z = 8;
    double sp = 0, dp = M_PI / 2.0;
    Cone cone(rmin1, rmax1, rmin2, rmax2, z, sp, dp);
    ValidateContainment(cone, [=](Vector3D const & p) {
        return InsideCone(p, rmin1, rmax1, rmin2, rmax2, z) && InsidePhiRange(p.GetX(), p.GetY(), sp, dp);
    }, 12, 10000, "ConePhiCut");
}

// =========================================================================
// Partial angular extent: Polycone
// =========================================================================

TEST(Containment, PolyconePhiCut) {
    // Quarter polycone: phi from 0 to pi/2
    std::vector<double> zp = {-5, 0, 5};
    std::vector<double> rmin_v = {1, 2, 1};
    std::vector<double> rmax_v = {4, 5, 3};
    double sp = 0, dp = M_PI / 2.0;
    Polycone pc(zp, rmin_v, rmax_v, sp, dp);
    ValidateContainment(pc, [&](Vector3D const & p) {
        return InsidePolycone(p, zp, rmin_v, rmax_v) && InsidePhiRange(p.GetX(), p.GetY(), sp, dp);
    }, 10, 10000, "PolyconePhiCut");
}

// =========================================================================
// Partial angular extent: CutTube
// =========================================================================

TEST(Containment, CutTubePhiCut) {
    // Quarter cut-tube: phi from 0 to pi/2
    double rmin = 1, rmax = 5, dz = 4;
    Vector3D low_norm(0, 0, -1);
    Vector3D high_norm(0, 0, 1);
    double sp = 0, dp = M_PI / 2.0;
    CutTube ct(rmin, rmax, dz, low_norm, high_norm, sp, dp);
    ValidateContainment(ct, [=](Vector3D const & p) {
        return InsideCutTube(p, rmin, rmax, dz, low_norm, high_norm) && InsidePhiRange(p.GetX(), p.GetY(), sp, dp);
    }, 10, 10000, "CutTubePhiCut");
}

// =========================================================================
// EllipticalTube
// =========================================================================
TEST(Containment, EllipticalTube) {
    double dx = 3, dy = 2, dz = 5;
    Placement pl(Vector3D(1, -2, 3));
    EllipticalTube et(pl, dx, dy, dz);
    ValidateContainment(et, [=](Vector3D const & p) {
        return InsideEllipticalTube(p, dx, dy, dz);
    }, 10, 10000, "EllipticalTube");
}

TEST(Containment, EllipticalTubeCircular) {
    double d = 4, dz = 6;
    EllipticalTube et(d, d, dz);
    ValidateContainment(et, [=](Vector3D const & p) {
        return InsideEllipticalTube(p, d, d, dz);
    }, 10, 10000, "EllipticalTubeCircular");
}

// =========================================================================
// CutTube
// =========================================================================
TEST(Containment, CutTubeFlat) {
    double rmin = 0, rmax = 5, dz = 4;
    Vector3D low_norm(0, 0, -1);
    Vector3D high_norm(0, 0, 1);
    CutTube ct(rmin, rmax, dz, low_norm, high_norm);
    ValidateContainment(ct, [=](Vector3D const & p) {
        return InsideCutTube(p, rmin, rmax, dz, low_norm, high_norm);
    }, 10, 10000, "CutTubeFlat");
}

TEST(Containment, CutTubeTilted) {
    double rmin = 0, rmax = 5, dz = 4;
    Vector3D low_norm(0.0, -0.2, -0.98);
    Vector3D high_norm(0.1, 0.0, 0.995);
    double lmag = std::sqrt(0.04 + 0.9604);
    low_norm = Vector3D(low_norm.GetX()/lmag, low_norm.GetY()/lmag, low_norm.GetZ()/lmag);
    double hmag = std::sqrt(0.01 + 0.990025);
    high_norm = Vector3D(high_norm.GetX()/hmag, high_norm.GetY()/hmag, high_norm.GetZ()/hmag);
    Placement pl(Vector3D(2, -1, 0));
    CutTube ct(pl, rmin, rmax, dz, low_norm, high_norm);
    ValidateContainment(ct, [=](Vector3D const & p) {
        return InsideCutTube(p, rmin, rmax, dz, low_norm, high_norm);
    }, 12, 10000, "CutTubeTilted");
}

TEST(Containment, CutTubeHollow) {
    double rmin = 2, rmax = 5, dz = 4;
    Vector3D low_norm(0, 0, -1);
    Vector3D high_norm(0, 0, 1);
    CutTube ct(rmin, rmax, dz, low_norm, high_norm);
    ValidateContainment(ct, [=](Vector3D const & p) {
        return InsideCutTube(p, rmin, rmax, dz, low_norm, high_norm);
    }, 10, 10000, "CutTubeHollow");
}

// =========================================================================
// Trap (General Trapezoid)
// =========================================================================
TEST(Containment, TrapAsTrd) {
    // theta=phi=alpha=0, symmetric => reduces to Trd
    double dz = 5, dy1 = 3, dx1 = 4, dx2 = 4, dy2 = 2, dx3 = 3, dx4 = 3;
    Trap trap(dz, 0, 0, dy1, dx1, dx2, 0, dy2, dx3, dx4, 0);
    ValidateContainment(trap, [=](Vector3D const & p) {
        return InsideTrap(p, dz, 0, 0, dy1, dx1, dx2, 0, dy2, dx3, dx4, 0);
    }, 12, 10000, "TrapAsTrd");
}

TEST(Containment, TrapGeneral) {
    // Exercise the full shear path (theta, phi, alpha all nonzero) with a
    // VALID trap: an affine shear of a box (constant dx, constant dy,
    // single alpha, uniform z-tilt) is a true parallelepiped, so its six
    // faces are planar by construction and containment is well defined.
    //
    // The earlier fixture used asymmetric dx1!=dx2, alpha1!=alpha2 etc.,
    // which yields non-planar side quads (an ill-defined G4Trap); the
    // analytic oracle and the solid then triangulate the quad
    // differently and disagreement is inherent, not a real bug. The old
    // InsideTrap "val > -0.01" skin existed only to paper over that.
    double dz = 5, theta = 0.15, phi = 0.3;
    double dy = 3, dx = 4, alpha = 0.2;
    Placement pl(Vector3D(1, -1, 2));
    Trap trap(pl, dz, theta, phi, dy, dx, dx, alpha, dy, dx, dx, alpha);
    ValidateContainment(trap, [=](Vector3D const & p) {
        return InsideTrap(p, dz, theta, phi, dy, dx, dx, alpha,
                          dy, dx, dx, alpha);
    }, 15, 10000, "TrapGeneral");
}

TEST(Containment, TrapAsymmetricAlpha0) {
    // Strong dx asymmetry with alpha=0 (planar faces, algorithm should be exact)
    double dz = 5, theta = 0, phi = 0;
    double dy1 = 3, dx1 = 4, dx2 = 1, alpha1 = 0;
    double dy2 = 3, dx3 = 4, dx4 = 1, alpha2 = 0;
    Trap trap(dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    ValidateContainment(trap, [=](Vector3D const & p) {
        return InsideTrap(p, dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    }, 12, 10000, "TrapAsymmetricAlpha0");
}

TEST(Containment, TrapNearTriangularAlpha0) {
    // Near-triangular (one face nearly collapsed)
    double dz = 5, theta = 0, phi = 0;
    double dy1 = 3, dx1 = 4, dx2 = 0.1, alpha1 = 0;
    double dy2 = 3, dx3 = 4, dx4 = 0.1, alpha2 = 0;
    Trap trap(dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    ValidateContainment(trap, [=](Vector3D const & p) {
        return InsideTrap(p, dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    }, 12, 10000, "TrapNearTriangularAlpha0");
}

TEST(Containment, TrapAsymWithTheta) {
    // Asymmetric dx with theta tilt (still alpha=0, so faces are planar)
    double dz = 5, theta = 0.2, phi = 0.5;
    double dy1 = 3, dx1 = 4, dx2 = 1, alpha1 = 0;
    double dy2 = 3, dx3 = 4, dx4 = 1, alpha2 = 0;
    Trap trap(dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    ValidateContainment(trap, [=](Vector3D const & p) {
        return InsideTrap(p, dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    }, 15, 10000, "TrapAsymWithTheta");
}

TEST(Containment, TrapTaperedAlpha0) {
    // Different top/bottom faces (full taper)
    double dz = 5, theta = 0, phi = 0;
    double dy1 = 3, dx1 = 4, dx2 = 2, alpha1 = 0;
    double dy2 = 2, dx3 = 3, dx4 = 1, alpha2 = 0;
    Trap trap(dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    ValidateContainment(trap, [=](Vector3D const & p) {
        return InsideTrap(p, dz, theta, phi, dy1, dx1, dx2, alpha1, dy2, dx3, dx4, alpha2);
    }, 12, 10000, "TrapTaperedAlpha0");
}

TEST(Containment, TrapAsymIntersectionCount) {
    // Verify that asymmetric alpha=0 Trap produces exactly 2 intersections
    // for center-piercing rays
    Trap trap(5, 0, 0, 3, 4, 1, 0, 3, 4, 1, 0);
    Vector3D dirs[] = {
        Vector3D(1, 0, 0), Vector3D(0, 1, 0), Vector3D(0, 0, 1),
        Vector3D(1, 1, 1), Vector3D(1, -1, 0.5), Vector3D(-0.3, 0.7, -0.6),
    };
    for(auto & d : dirs) {
        Vector3D dir = d;
        dir.normalize();
        Vector3D origin = Vector3D(0, 0, 0) - dir * 50.0;
        auto hits = trap.Intersections(origin, dir);
        EXPECT_EQ(hits.size(), 2u)
            << "TrapAsym center ray dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ()
            << ") got " << hits.size();
    }
}

TEST(Containment, TrapAsymEnterExitFlags) {
    // For an asymmetric Trap, verify entering/exiting flags are correct
    Trap trap(5, 0, 0, 3, 4, 1, 0, 3, 4, 1, 0);
    // Z-axis ray
    auto hits = trap.Intersections(Vector3D(0, 0, -20), Vector3D(0, 0, 1));
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_TRUE(hits[0].entering);
    EXPECT_FALSE(hits[1].entering);
    EXPECT_NEAR(hits[0].position.GetZ(), -5.0, 1e-6);
    EXPECT_NEAR(hits[1].position.GetZ(), 5.0, 1e-6);
}

// =========================================================================
// Ellipsoid
// =========================================================================
TEST(Containment, Ellipsoid) {
    double ax = 5, by = 3, cz = 4;
    Placement pl(Vector3D(2, -1, 0));
    Ellipsoid ell(pl, ax, by, cz);
    ValidateContainment(ell, [=](Vector3D const & p) {
        return InsideEllipsoid(p, ax, by, cz, -cz, cz);
    }, 10, 10000, "Ellipsoid");
}

TEST(Containment, EllipsoidWithZCuts) {
    double ax = 5, by = 3, cz = 4;
    double zcut1 = -2, zcut2 = 3;
    Ellipsoid ell(ax, by, cz, zcut1, zcut2);
    ValidateContainment(ell, [=](Vector3D const & p) {
        return InsideEllipsoid(p, ax, by, cz, zcut1, zcut2);
    }, 10, 10000, "EllipsoidWithZCuts");
}

TEST(Containment, EllipsoidSphere) {
    // Equal semi-axes = sphere
    double r = 5;
    Ellipsoid ell(r, r, r);
    ValidateContainment(ell, [=](Vector3D const & p) {
        return InsideEllipsoid(p, r, r, r, -r, r);
    }, 10, 10000, "EllipsoidSphere");
}

// =========================================================================
// Para (Parallelepiped)
// =========================================================================
TEST(Containment, ParaBox) {
    // alpha=theta=phi=0 => box
    double dx = 4, dy = 3, dz = 5;
    Para para(dx, dy, dz, 0, 0, 0);
    ValidateContainment(para, [=](Vector3D const & p) {
        return InsidePara(p, dx, dy, dz, 0, 0, 0);
    }, 10, 10000, "ParaBox");
}

TEST(Containment, ParaSheared) {
    double dx = 4, dy = 3, dz = 5;
    double alpha = 0.3, theta = 0.2, phi = 0.5;
    Placement pl(Vector3D(1, -2, 3));
    Para para(pl, dx, dy, dz, alpha, theta, phi);
    ValidateContainment(para, [=](Vector3D const & p) {
        return InsidePara(p, dx, dy, dz, alpha, theta, phi);
    }, 15, 10000, "ParaSheared");
}

TEST(Containment, TriangularMeshCube) {
    // Build a 2x2x2 cube centered at origin as a triangulated mesh
    double s = 3.0;
    Vector3D v[8] = {
        Vector3D(-s,-s,-s), Vector3D(s,-s,-s), Vector3D(s,s,-s), Vector3D(-s,s,-s),
        Vector3D(-s,-s, s), Vector3D(s,-s, s), Vector3D(s,s, s), Vector3D(-s,s, s)
    };
    std::vector<std::array<Vector3D, 3>> triangles;
    auto addQuad = [&](Vector3D a, Vector3D b, Vector3D c, Vector3D d) {
        triangles.push_back({{a, b, c}});
        triangles.push_back({{a, c, d}});
    };
    addQuad(v[3], v[2], v[1], v[0]); // -z face
    addQuad(v[4], v[5], v[6], v[7]); // +z face
    addQuad(v[0], v[1], v[5], v[4]); // -y face
    addQuad(v[1], v[2], v[6], v[5]); // +x face
    addQuad(v[2], v[3], v[7], v[6]); // +y face
    addQuad(v[3], v[0], v[4], v[7]); // -x face

    Placement pl(Vector3D(2, -1, 4));
    TriangularMesh mesh(pl, triangles);
    ValidateContainment(mesh, [=](Vector3D const & p) {
        return std::fabs(p.GetX()) < s && std::fabs(p.GetY()) < s && std::fabs(p.GetZ()) < s;
    }, s * 2, 10000, "TriangularMeshCube");
}

// =========================================================================
// GenericPolycone
// =========================================================================

bool InsideGenericPolycone(Vector3D const & p,
                           std::vector<double> const & rv,
                           std::vector<double> const & zv) {
    double r = std::sqrt(p.GetX()*p.GetX() + p.GetY()*p.GetY());
    double z = p.GetZ();
    size_t n = rv.size();
    int crossings = 0;
    for(size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        double ri = rv[i], zi = zv[i];
        double rj = rv[j], zj = zv[j];
        if((zi <= z && zj > z) || (zj <= z && zi > z)) {
            double r_cross = ri + (z - zi) / (zj - zi) * (rj - ri);
            if(r < r_cross) {
                crossings++;
            }
        }
    }
    return (crossings % 2) == 1;
}

TEST(Containment, GenericPolyconeSolid) {
    // Simple solid triangle outline (equivalent to a solid cone)
    std::vector<double> rv = {0, 5, 0};
    std::vector<double> zv = {-5, 0, 5};
    GenericPolycone gpc(rv, zv);
    ValidateContainment(gpc, [&](Vector3D const & p) {
        return InsideGenericPolycone(p, rv, zv);
    }, 6, 10000, "GenericPolyconeSolid");
}

TEST(Containment, GenericPolyconeHollow) {
    // Hollow shape: hexagonal outline like the celeritas test case
    std::vector<double> rv = {3, 4.5, 5, 3.5, 3, 2};
    std::vector<double> zv = {-5, 0, 5, 5, 0, -5};
    GenericPolycone gpc(rv, zv);
    ValidateContainment(gpc, [&](Vector3D const & p) {
        return InsideGenericPolycone(p, rv, zv);
    }, 6, 10000, "GenericPolyconeHollow");
}

TEST(Containment, GenericPolyconeUShape) {
    // U-shaped cross-section: creates a hollow annular shell
    std::vector<double> rv = {5, 5, 3, 3, 5};
    std::vector<double> zv = {-4, 4, 4, -4, -4};
    // Remove the degenerate closing edge (first == last)
    rv.pop_back(); zv.pop_back();
    GenericPolycone gpc(rv, zv);
    ValidateContainment(gpc, [&](Vector3D const & p) {
        return InsideGenericPolycone(p, rv, zv);
    }, 6, 10000, "GenericPolyconeUShape");
}

TEST(Containment, GenericPolyconePhiCut) {
    std::vector<double> rv = {0, 5, 0};
    std::vector<double> zv = {-5, 0, 5};
    double sp = M_PI / 4;
    double dp = M_PI;
    GenericPolycone gpc(rv, zv, sp, dp);
    ValidateContainment(gpc, [&](Vector3D const & p) {
        return InsideGenericPolycone(p, rv, zv) && InsidePhiRange(p.GetX(), p.GetY(), sp, dp);
    }, 6, 10000, "GenericPolyconePhiCut");
}

TEST(Containment, GenericPolyconeRotated) {
    std::vector<double> rv = {3, 4.5, 5, 3.5, 3, 2};
    std::vector<double> zv = {-5, 0, 5, 5, 0, -5};
    Quaternion q(std::cos(0.5), std::sin(0.5)*0.577, std::sin(0.5)*0.577, std::sin(0.5)*0.577);
    Placement pl(Vector3D(1, -1, 2), q);
    GenericPolycone gpc(pl, rv, zv);
    ValidateContainment(gpc, [&](Vector3D const & p) {
        return InsideGenericPolycone(p, rv, zv);
    }, 8, 10000, "GenericPolyconeRotated");
}

TEST(Containment, GenericPolyconeEquivalentToPolycone) {
    // A genericPolycone that traces the same shape as a regular polycone:
    // outer boundary going up, then inner boundary going back down.
    // Polycone: z={-5,0,5}, rmin={2,3,3.5}, rmax={3,4.5,5}
    std::vector<double> zp = {-5, 0, 5};
    std::vector<double> rmin_v = {2, 3, 3.5};
    std::vector<double> rmax_v = {3, 4.5, 5};
    Polycone pc(zp, rmin_v, rmax_v);

    // Equivalent genericPolycone outline: outer going up, inner going down
    std::vector<double> rv = {3, 4.5, 5, 3.5, 3, 2};
    std::vector<double> zv = {-5, 0, 5, 5, 0, -5};
    GenericPolycone gpc(rv, zv);

    int mismatches = 0;
    for(int i = 0; i < 30000; ++i) {
        Vector3D pos = RandomPoint(6);
        Vector3D dir = RandomDirection();
        bool pc_inside = pc.IsInside(pos, dir);
        bool gpc_inside = gpc.IsInside(pos, dir);
        if(pc_inside != gpc_inside) {
            mismatches++;
            if(mismatches <= 3) {
                std::cerr << "Equivalence mismatch #" << mismatches
                          << ": pos=(" << pos.GetX() << "," << pos.GetY() << "," << pos.GetZ() << ")"
                          << " polycone=" << pc_inside << " generic=" << gpc_inside
                          << std::endl;
            }
        }
    }
    EXPECT_EQ(mismatches, 0) << mismatches << " mismatches between Polycone and equivalent GenericPolycone";
}

TEST(Containment, GenericPolyconeCrystal) {
    // HPGe detector crystal profile modeled after LEGEND IC-type crystals.
    // Features: tapered outer surface, passivation groove at bottom,
    // central borehole cavity, multiple z-direction reversals.
    //
    // Cross-section (z=0 at top, z increases downward):
    //   - Cavity opening at top (r=5, z=0)
    //   - Tapered outer wall (r=38 at z=0 to r=40 at z=70)
    //   - Straight outer cylinder (r=40, z=70 to z=80)
    //   - Groove at bottom (outer r=30, inner r=15, z=77 to 80)
    //   - Flat bottom (r=0 to r=15 at z=80)
    //   - Inner borehole (r=0 at z=80 to z=60, then r=5 back up to z=0)
    std::vector<double> rv = {5, 38, 40, 40, 30, 30, 15, 15, 0, 0, 5};
    std::vector<double> zv = {0, 0, 70, 80, 80, 77, 77, 80, 80, 60, 60};
    GenericPolycone gpc(rv, zv);
    ValidateContainment(gpc, [&](Vector3D const & p) {
        return InsideGenericPolycone(p, rv, zv);
    }, 45, 20000, "GenericPolyconeCrystal");
}

// =========================================================================
// ExtrPoly extended coverage and cross-validation
// =========================================================================

TEST(Containment, ExtrPolyTriangle) {
    // Triangle prism: V0=(5,-5), V1=(-5,-5), V2=(0,5)
    std::vector<std::vector<double>> polygon = {{5,-5},{-5,-5},{0,5}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-5, off, 1.0),
        ExtrPoly::ZSection(5, off, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
        double px = p.GetX(), py = p.GetY();
        bool in_triangle = (py > -5.0 &&
                            py < 2.0*px + 5.0 &&
                            py < -2.0*px + 5.0);
        return in_triangle && p.GetZ() > -5.0 && p.GetZ() < 5.0;
    }, 8, 10000, "ExtrPolyTriangle");
}

TEST(Containment, ExtrPolyHexagon) {
    // Regular hexagon with circumradius 5, extruded along z
    double r = 5.0;
    std::vector<std::vector<double>> polygon;
    for(int k = 0; k < 6; ++k) {
        double a = k * M_PI / 3.0;
        polygon.push_back({r * std::cos(a), r * std::sin(a)});
    }
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off, 1.0),
        ExtrPoly::ZSection(4, off, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [=](Vector3D const & p) {
        bool in_hex = InsideRegularPolygon(p.GetX(), p.GetY(), 6, r, 0);
        return in_hex && p.GetZ() > -4.0 && p.GetZ() < 4.0;
    }, 8, 10000, "ExtrPolyHexagon");
}

TEST(Containment, ExtrPolyExtremeAspectRatio) {
    // Very elongated rectangle (100:1 aspect ratio)
    std::vector<std::vector<double>> polygon = {{-50,-0.5},{50,-0.5},{50,0.5},{-50,0.5}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-1, off, 1.0),
        ExtrPoly::ZSection(1, off, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
        return std::fabs(p.GetX()) < 50.0
            && std::fabs(p.GetY()) < 0.5
            && std::fabs(p.GetZ()) < 1.0;
    }, 55, 20000, "ExtrPolyExtremeAspectRatio");
}

TEST(Containment, ExtrPolyCombinedOffsetScale) {
    // Multi-section with both offset and scale varying
    std::vector<std::vector<double>> polygon = {{-2,-2},{2,-2},{2,2},{-2,2}};
    double off_bot[2] = {0, 0};
    double off_mid[2] = {0.5, 0.5};
    double off_top[2] = {1, 1};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off_bot, 1.0),
        ExtrPoly::ZSection(0, off_mid, 1.5),
        ExtrPoly::ZSection(4, off_top, 0.8)
    };
    ExtrPoly ep(polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
        double pz = p.GetZ();
        if(pz <= -4.0 || pz >= 4.0) return false;
        double s, ox, oy;
        if(pz < 0.0) {
            double frac = (pz + 4.0) / 4.0;
            s = 1.0 + (1.5 - 1.0) * frac;
            ox = 0.0 + (0.5 - 0.0) * frac;
            oy = 0.0 + (0.5 - 0.0) * frac;
        } else {
            double frac = pz / 4.0;
            s = 1.5 + (0.8 - 1.5) * frac;
            ox = 0.5 + (1.0 - 0.5) * frac;
            oy = 0.5 + (1.0 - 0.5) * frac;
        }
        return std::fabs(p.GetX() - ox) < 2.0 * s
            && std::fabs(p.GetY() - oy) < 2.0 * s;
    }, 10, 10000, "ExtrPolyCombinedOffsetScale");
}

TEST(Containment, ExtrPolyRotated) {
    // Square extrusion with rotated placement
    std::vector<std::vector<double>> polygon = {{-3,-3},{3,-3},{3,3},{-3,3}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-4, off, 1.0),
        ExtrPoly::ZSection(4, off, 1.0)
    };
    Quaternion q(std::cos(0.4), std::sin(0.4)*0.577, std::sin(0.4)*0.577, std::sin(0.4)*0.577);
    Placement pl(Vector3D(2, -1, 3), q);
    ExtrPoly ep(pl, polygon, zsecs);
    ValidateContainment(ep, [](Vector3D const & p) {
        return std::fabs(p.GetX()) < 3 && std::fabs(p.GetY()) < 3 && std::fabs(p.GetZ()) < 4;
    }, 12, 10000, "ExtrPolyRotated");
}

TEST(Containment, ExtrPolyVsBoxCrossValidation) {
    // Cross-validate that an ExtrPoly square prism agrees with a Box
    std::vector<std::vector<double>> polygon = {{-4,-3},{4,-3},{4,3},{-4,3}};
    double off[2] = {0, 0};
    std::vector<ExtrPoly::ZSection> zsecs = {
        ExtrPoly::ZSection(-5, off, 1.0),
        ExtrPoly::ZSection(5, off, 1.0)
    };
    ExtrPoly ep(polygon, zsecs);
    Box box(8, 6, 10); // full widths = 2*4, 2*3, 2*5

    int mismatches = 0;
    std::mt19937 rng(99999);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for(int i = 0; i < 30000; ++i) {
        Vector3D pos(u(rng) * 8, u(rng) * 6, u(rng) * 8);
        Vector3D dir(u(rng), u(rng), u(rng));
        double mag = std::sqrt(dir.GetX()*dir.GetX() + dir.GetY()*dir.GetY() + dir.GetZ()*dir.GetZ());
        if(mag < 1e-12) continue;
        dir = Vector3D(dir.GetX()/mag, dir.GetY()/mag, dir.GetZ()/mag);
        bool ep_inside = ep.IsInside(pos, dir);
        bool box_inside = box.IsInside(pos, dir);
        if(ep_inside != box_inside) {
            mismatches++;
            if(mismatches <= 3) {
                std::cerr << "ExtrPoly vs Box mismatch at ("
                    << pos.GetX() << "," << pos.GetY() << "," << pos.GetZ()
                    << ") ep=" << ep_inside << " box=" << box_inside << std::endl;
            }
        }
    }
    EXPECT_EQ(mismatches, 0) << "ExtrPoly vs Box: " << mismatches << " mismatches";
}
