// Tests for BVH spatial acceleration:
//   - Traversal finds all intersections (no missed sectors)
//   - Results match brute-force linear scan
//   - Edge cases: axis-aligned rays, single volume, empty model

#include <cmath>
#include <random>
#include <string>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/Torus.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/Path.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/detector/ConstantDensityDistribution.h"

using namespace siren::math;
using namespace siren::geometry;
using namespace siren::detector;

namespace {

std::mt19937 bvh_rng(42);

void BuildTestModel(DetectorModel & dm, int n_volumes) {
    dm.ClearSectors();

    DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = Sphere(100, 0).create();
    world.density = ConstantDensityDistribution(1.0).create();
    dm.AddSector(world);

    std::uniform_real_distribution<double> pos_dist(-30.0, 30.0);
    std::uniform_real_distribution<double> size_dist(2.0, 10.0);

    for(int i = 0; i < n_volumes; ++i) {
        double x = pos_dist(bvh_rng);
        double y = pos_dist(bvh_rng);
        double z = pos_dist(bvh_rng);

        DetectorSector sector;
        sector.name = "vol_" + std::to_string(i);
        sector.material_id = 1;
        sector.level = i + 1;
        sector.density = ConstantDensityDistribution(2.0).create();

        int shape_type = i % 8;
        Placement pl(Vector3D(x, y, z));
        if(shape_type == 0) {
            double r = size_dist(bvh_rng);
            sector.geo = Sphere(pl, r, 0).create();
        } else if(shape_type == 1) {
            double sx = size_dist(bvh_rng);
            double sy = size_dist(bvh_rng);
            double sz = size_dist(bvh_rng);
            sector.geo = Box(pl, sx, sy, sz).create();
        } else if(shape_type == 2) {
            double r = size_dist(bvh_rng);
            double h = size_dist(bvh_rng);
            sector.geo = Cylinder(pl, r, 0, h).create();
        } else if(shape_type == 3) {
            double rmin1 = size_dist(bvh_rng) * 0.2;
            double rmax1 = rmin1 + size_dist(bvh_rng);
            double rmin2 = size_dist(bvh_rng) * 0.2;
            double rmax2 = rmin2 + size_dist(bvh_rng);
            double h = size_dist(bvh_rng);
            sector.geo = Cone(pl, rmin1, rmax1, rmin2, rmax2, h).create();
        } else if(shape_type == 4) {
            double dx1 = size_dist(bvh_rng);
            double dx2 = size_dist(bvh_rng);
            double dy1 = size_dist(bvh_rng);
            double dy2 = size_dist(bvh_rng);
            double dz = size_dist(bvh_rng);
            sector.geo = Trd(pl, dx1, dx2, dy1, dy2, dz).create();
        } else if(shape_type == 5) {
            double h = size_dist(bvh_rng);
            double rmax = size_dist(bvh_rng);
            std::vector<double> zp = {-h, 0.0, h};
            std::vector<double> rn = {0.0, 0.0, 0.0};
            std::vector<double> rx = {rmax * 0.5, rmax, rmax * 0.5};
            sector.geo = Polycone(pl, zp, rn, rx).create();
        } else if(shape_type == 6) {
            double h = size_dist(bvh_rng);
            double rmax = size_dist(bvh_rng);
            std::vector<double> zp = {-h, 0.0, h};
            std::vector<double> rn = {0.0, 0.0, 0.0};
            std::vector<double> rx = {rmax * 0.5, rmax, rmax * 0.5};
            sector.geo = Polyhedra(pl, 6, 0.0, zp, rn, rx).create();
        } else {
            double rtor = size_dist(bvh_rng) + 3.0;
            double rmax = size_dist(bvh_rng) * 0.3;
            sector.geo = Torus(pl, rtor, rmax, 0.0).create();
        }
        dm.AddSector(sector);
    }
}

// Brute-force: compute intersections by testing every sector
Geometry::IntersectionList BruteForceIntersections(
        DetectorModel const & dm,
        Vector3D const & position,
        Vector3D const & direction) {
    Geometry::IntersectionList result;
    result.position = position;
    result.direction = direction;
    for(auto const & sector : dm.GetSectors()) {
        std::vector<Geometry::Intersection> hits = sector.geo->Intersections(position, direction);
        size_t prev = result.intersections.size();
        result.intersections.insert(result.intersections.end(), hits.begin(), hits.end());
        for(size_t j = prev; j < result.intersections.size(); ++j) {
            result.intersections[j].hierarchy = sector.level;
            result.intersections[j].matID = sector.material_id;
        }
    }
    DetectorModel::SortIntersections(result);
    return result;
}

bool CompareIntersectionLists(
        Geometry::IntersectionList const & a,
        Geometry::IntersectionList const & b) {
    if(a.intersections.size() != b.intersections.size()) return false;
    for(size_t j = 0; j < a.intersections.size(); ++j) {
        if(std::fabs(a.intersections[j].distance - b.intersections[j].distance) > 1e-6)
            return false;
        if(a.intersections[j].hierarchy != b.intersections[j].hierarchy)
            return false;
    }
    return true;
}

} // namespace

TEST(BVH, IntersectionsMatchBruteForce) {
    DetectorModel dm;
    BuildTestModel(dm, 30);

    std::uniform_real_distribution<double> pos_dist(-40.0, 40.0);
    std::uniform_real_distribution<double> dir_dist(-1.0, 1.0);

    int mismatches = 0;
    int total_rays = 200;

    for(int i = 0; i < total_rays; ++i) {
        Vector3D pos(pos_dist(bvh_rng), pos_dist(bvh_rng), pos_dist(bvh_rng));
        Vector3D dir(dir_dist(bvh_rng), dir_dist(bvh_rng), dir_dist(bvh_rng));
        if(dir.magnitude() < 1e-10) continue;
        dir.normalize();

        auto bvh_result = dm.GetIntersections(DetectorPosition(pos), DetectorDirection(dir));
        auto brute_result = BruteForceIntersections(dm, pos, dir);

        if(!CompareIntersectionLists(bvh_result, brute_result)) {
            mismatches++;
        }
    }
    EXPECT_EQ(mismatches, 0) << "BVH results diverged from brute-force on " << mismatches << "/" << total_rays << " rays";
}

TEST(BVH, AxisAlignedRays) {
    DetectorModel dm;
    BuildTestModel(dm, 15);

    Vector3D origins[] = {
        Vector3D(0, 0, 0),
        Vector3D(20, 0, 0),
        Vector3D(0, 20, 0),
        Vector3D(0, 0, 20),
    };
    Vector3D directions[] = {
        Vector3D(1, 0, 0),
        Vector3D(0, 1, 0),
        Vector3D(0, 0, 1),
        Vector3D(-1, 0, 0),
        Vector3D(0, -1, 0),
        Vector3D(0, 0, -1),
    };

    for(auto const & orig : origins) {
        for(auto const & dir : directions) {
            auto bvh_result = dm.GetIntersections(DetectorPosition(orig), DetectorDirection(dir));
            auto brute_result = BruteForceIntersections(dm, orig, dir);

            EXPECT_TRUE(CompareIntersectionLists(bvh_result, brute_result))
                << "Axis-aligned ray mismatch at origin=(" << orig.GetX() << "," << orig.GetY() << "," << orig.GetZ()
                << ") dir=(" << dir.GetX() << "," << dir.GetY() << "," << dir.GetZ() << ")"
                << " BVH=" << bvh_result.intersections.size()
                << " brute=" << brute_result.intersections.size();
        }
    }
}

TEST(BVH, LargeModelNoMisses) {
    DetectorModel dm;
    BuildTestModel(dm, 100);

    std::uniform_real_distribution<double> pos_dist(-50.0, 50.0);
    std::uniform_real_distribution<double> dir_dist(-1.0, 1.0);

    int mismatches = 0;
    int total_rays = 100;

    for(int i = 0; i < total_rays; ++i) {
        Vector3D pos(pos_dist(bvh_rng), pos_dist(bvh_rng), pos_dist(bvh_rng));
        Vector3D dir(dir_dist(bvh_rng), dir_dist(bvh_rng), dir_dist(bvh_rng));
        if(dir.magnitude() < 1e-10) continue;
        dir.normalize();

        auto bvh_result = dm.GetIntersections(DetectorPosition(pos), DetectorDirection(dir));
        auto brute_result = BruteForceIntersections(dm, pos, dir);

        if(!CompareIntersectionLists(bvh_result, brute_result)) {
            mismatches++;
        }
    }
    EXPECT_EQ(mismatches, 0) << "BVH missed intersections on " << mismatches << "/" << total_rays << " rays with 100 volumes";
}

TEST(BVH, SingleVolume) {
    DetectorModel dm;
    dm.ClearSectors();

    DetectorSector sector;
    sector.name = "only";
    sector.material_id = 0;
    sector.level = 0;
    sector.geo = Sphere(10, 0).create();
    sector.density = ConstantDensityDistribution(1.0).create();
    dm.AddSector(sector);

    auto result = dm.GetIntersections(
        DetectorPosition(Vector3D(0, 0, -20)),
        DetectorDirection(Vector3D(0, 0, 1)));

    EXPECT_EQ(result.intersections.size(), 2u);
}

TEST(BVH, EmptyModel) {
    DetectorModel dm;
    dm.ClearSectors();

    auto result = dm.GetIntersections(
        DetectorPosition(Vector3D(0, 0, 0)),
        DetectorDirection(Vector3D(0, 0, 1)));

    EXPECT_EQ(result.intersections.size(), 0u);
}

// =========================================================================
// BVH with degenerate BooleanGeometry AABB (inverted min>max)
// Ensures invalid AABBs don't bloat the BVH or cause false negatives.
// =========================================================================

TEST(BVH, DegenrateBooleanAABBDoesNotCorruptBVH) {
    DetectorModel dm;
    dm.ClearSectors();

    // World sphere
    DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = Sphere(100, 0).create();
    world.density = ConstantDensityDistribution(1.0).create();
    dm.AddSector(world);

    // Non-overlapping INTERSECTION boolean => produces invalid (inverted) AABB
    auto left = Box(Placement(Vector3D(-50, 0, 0)), 5, 5, 5).create();
    auto right = Box(Placement(Vector3D(50, 0, 0)), 5, 5, 5).create();
    auto degenerate = std::make_shared<BooleanGeometry>(
        BooleanOperation::INTERSECTION,
        std::const_pointer_cast<const Geometry>(left),
        std::const_pointer_cast<const Geometry>(right));

    DetectorSector degen_sector;
    degen_sector.name = "degenerate";
    degen_sector.material_id = 1;
    degen_sector.level = 1;
    degen_sector.geo = degenerate;
    degen_sector.density = ConstantDensityDistribution(2.0).create();
    dm.AddSector(degen_sector);

    // Add a normal box that should definitely be found by the BVH
    DetectorSector normal;
    normal.name = "normal_box";
    normal.material_id = 2;
    normal.level = 2;
    normal.geo = Box(5, 5, 5).create();
    normal.density = ConstantDensityDistribution(3.0).create();
    dm.AddSector(normal);

    // Ray through the origin should find the normal box
    auto result = dm.GetIntersections(
        DetectorPosition(Vector3D(0, 0, -20)),
        DetectorDirection(Vector3D(0, 0, 1)));

    // Should find at least 2 hits from the normal box + 2 from world
    bool found_normal = false;
    for(auto const & isect : result.intersections) {
        if(isect.matID == 2) found_normal = true;
    }
    EXPECT_TRUE(found_normal)
        << "Normal box should be found despite degenerate boolean in BVH";
}

// =========================================================================
// Path::SetPointsWithRay(DetectorPosition) collinearity reuse with rotated
// detector frame. Ensures the intersection cache is correctly invalidated
// when the detector has a non-identity rotation.
// =========================================================================

TEST(BVH, PathDetectorFrameCollinearity) {
    auto dm_ptr = std::make_shared<DetectorModel>();
    dm_ptr->ClearSectors();

    // Simple world + inner box
    DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = Sphere(100, 0).create();
    world.density = ConstantDensityDistribution(1.0).create();
    dm_ptr->AddSector(world);

    DetectorSector inner;
    inner.name = "target";
    inner.material_id = 1;
    inner.level = 1;
    inner.geo = Box(10, 10, 10).create();
    inner.density = ConstantDensityDistribution(2.0).create();
    dm_ptr->AddSector(inner);

    // 90-degree rotation around z-axis:
    // Detector +x maps to Geometry +y, Detector +y maps to Geometry -x
    double angle = M_PI / 2.0;
    Quaternion rot_z(0, 0, std::sin(angle / 2.0), std::cos(angle / 2.0));
    dm_ptr->SetDetectorRotation(rot_z);

    siren::detector::Path path(dm_ptr);

    // First call: ray along detector +x direction (= geometry +y direction)
    path.SetPointsWithRay(
        DetectorPosition(Vector3D(0, 0, 0)),
        DetectorDirection(Vector3D(1, 0, 0)),
        50.0);
    path.EnsureIntersections();
    auto isects1 = path.GetIntersections();

    // Second call: same ray in detector coords, but different geometry-frame ray.
    // With the rotation, detector +y maps to geometry -x.
    // If the frame mismatch bug were present, this might incorrectly reuse
    // the cached intersections from the first call.
    path.SetPointsWithRay(
        DetectorPosition(Vector3D(0, 0, 0)),
        DetectorDirection(Vector3D(0, 1, 0)),
        50.0);
    path.EnsureIntersections();
    auto isects2 = path.GetIntersections();

    // Both should have intersections (the box is at origin, both rays pass through it)
    EXPECT_GE(isects1.intersections.size(), 2u)
        << "First ray should intersect the box";
    EXPECT_GE(isects2.intersections.size(), 2u)
        << "Second ray should intersect the box";

    // The two rays are NOT collinear in geometry frame
    // (they are perpendicular: one along +y, one along -x)
    // so their intersections should be different
    bool directions_differ = false;
    if(isects1.intersections.size() == isects2.intersections.size()) {
        for(size_t i = 0; i < isects1.intersections.size(); ++i) {
            if(std::fabs(isects1.intersections[i].distance - isects2.intersections[i].distance) > 1e-6) {
                directions_differ = true;
                break;
            }
        }
    } else {
        directions_differ = true;
    }
    // If we have a symmetric box at origin, the distances might be identical by symmetry.
    // But the intersection positions should definitely differ.
    if(!directions_differ && !isects1.intersections.empty()) {
        auto p1 = isects1.intersections[0].position;
        auto p2 = isects2.intersections[0].position;
        double pos_diff = std::sqrt(
            std::pow(p1.GetX() - p2.GetX(), 2) +
            std::pow(p1.GetY() - p2.GetY(), 2) +
            std::pow(p1.GetZ() - p2.GetZ(), 2));
        EXPECT_GT(pos_diff, 0.1)
            << "Perpendicular rays should produce different intersection positions";
    }
}

// =========================================================================
// Path collinearity reuse with translated detector origin
// =========================================================================

TEST(BVH, PathDetectorOriginCollinearity) {
    auto dm_ptr = std::make_shared<DetectorModel>();
    dm_ptr->ClearSectors();

    DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = Sphere(200, 0).create();
    world.density = ConstantDensityDistribution(1.0).create();
    dm_ptr->AddSector(world);

    // Box at geometry origin
    DetectorSector target;
    target.name = "target";
    target.material_id = 1;
    target.level = 1;
    target.geo = Box(10, 10, 10).create();
    target.density = ConstantDensityDistribution(2.0).create();
    dm_ptr->AddSector(target);

    // Shift detector origin: detector (0,0,0) maps to geometry (50,0,0)
    dm_ptr->SetDetectorOrigin(GeometryPosition(Vector3D(50, 0, 0)));

    siren::detector::Path path(dm_ptr);

    // Ray from detector origin along +z. In geometry frame this starts at (50,0,0)
    // and goes along +z, missing the box at origin entirely.
    path.SetPointsWithRay(
        DetectorPosition(Vector3D(0, 0, 0)),
        DetectorDirection(Vector3D(0, 0, 1)),
        100.0);
    path.EnsureIntersections();
    auto isects1 = path.GetIntersections();

    // Count intersections with the target box (material_id = 1)
    int box_hits_1 = 0;
    for(auto const & isect : isects1.intersections) {
        if(isect.matID == 1) box_hits_1++;
    }

    // Now ray from detector (-50, 0, 0) along +z.
    // In geometry frame this is (0,0,0) along +z, which hits the box.
    path.SetPointsWithRay(
        DetectorPosition(Vector3D(-50, 0, 0)),
        DetectorDirection(Vector3D(0, 0, 1)),
        100.0);
    path.EnsureIntersections();
    auto isects2 = path.GetIntersections();

    int box_hits_2 = 0;
    for(auto const & isect : isects2.intersections) {
        if(isect.matID == 1) box_hits_2++;
    }

    // First ray should miss the box, second should hit it
    EXPECT_EQ(box_hits_1, 0) << "Ray from (50,0,0) in geo frame should miss box at origin";
    EXPECT_GE(box_hits_2, 2) << "Ray from (0,0,0) in geo frame should hit box at origin";
}
