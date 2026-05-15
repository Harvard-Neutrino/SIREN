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
#include "SIREN/geometry/Placement.h"
#include "SIREN/detector/DetectorModel.h"
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

        int shape_type = i % 3;
        Placement pl(Vector3D(x, y, z));
        if(shape_type == 0) {
            double r = size_dist(bvh_rng);
            sector.geo = Sphere(pl, r, 0).create();
        } else if(shape_type == 1) {
            double sx = size_dist(bvh_rng);
            double sy = size_dist(bvh_rng);
            double sz = size_dist(bvh_rng);
            sector.geo = Box(pl, sx, sy, sz).create();
        } else {
            double r = size_dist(bvh_rng);
            double h = size_dist(bvh_rng);
            sector.geo = Cylinder(pl, r, 0, h).create();
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
