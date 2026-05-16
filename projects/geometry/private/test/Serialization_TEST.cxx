// Serialization roundtrip tests for all geometry shapes.
// Verifies that serialize -> deserialize through cereal produces
// identical objects, including derived boolean flags that are
// recomputed during deserialization (Sphere phi/theta cuts,
// Torus phi cut).

#include <gtest/gtest.h>
#include <sstream>
#include <memory>
#include <cmath>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>

#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/Torus.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/EllipticalTube.h"
#include "SIREN/geometry/CutTube.h"
#include "SIREN/geometry/Trap.h"
#include "SIREN/geometry/Ellipsoid.h"
#include "SIREN/geometry/Para.h"
#include "SIREN/geometry/serializable.h"
#include "SIREN/math/Vector3D.h"

using namespace siren::geometry;
using namespace siren::math;

namespace {

// Serialize a shared_ptr<Geometry> to binary and back.
// Returns the deserialized pointer.
std::shared_ptr<Geometry> RoundtripGeometry(std::shared_ptr<Geometry> original) {
    std::stringstream ss;
    {
        cereal::BinaryOutputArchive oarchive(ss);
        oarchive(original);
    }
    std::shared_ptr<Geometry> loaded;
    {
        cereal::BinaryInputArchive iarchive(ss);
        iarchive(loaded);
    }
    return loaded;
}

// Count the number of forward intersections a geometry reports
// for a given ray.
int CountIntersections(Geometry const & geo,
                       Vector3D const & pos,
                       Vector3D const & dir) {
    auto hits = geo.Intersections(pos, dir);
    int count = 0;
    for(auto const & h : hits) {
        if(h.distance > Geometry::GEOMETRY_PRECISION) {
            ++count;
        }
    }
    return count;
}

} // namespace

// ---- Individual shape roundtrip tests ----

TEST(SerializationRoundtrip, Box) {
    auto box = std::make_shared<Box>(1.0, 2.0, 3.0);
    auto loaded = RoundtripGeometry(box);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*box, *loaded);
}

TEST(SerializationRoundtrip, Cylinder) {
    auto cyl = std::make_shared<Cylinder>(1.0, 0.5, 3.0);
    auto loaded = RoundtripGeometry(cyl);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*cyl, *loaded);
}

TEST(SerializationRoundtrip, SphereFull) {
    auto sphere_full = std::make_shared<Sphere>(5.0, 2.0);
    auto loaded = RoundtripGeometry(sphere_full);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*sphere_full, *loaded);
}

TEST(SerializationRoundtrip, SphereWithCuts) {
    // Sphere with phi and theta cuts -- tests that derived flags
    // (has_phi_cut_, has_theta_cut_) are recomputed on deserialization
    auto sphere_cut = std::make_shared<Sphere>(
        5.0, 2.0, 0.0, M_PI, 0.2, 1.0);
    auto loaded = RoundtripGeometry(sphere_cut);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*sphere_cut, *loaded);

    // Verify that intersection behavior matches, proving
    // that the derived flags were correctly recomputed.
    // Use a variety of rays: some that hit the cut region,
    // some that miss entirely.
    Vector3D origins[] = {
        Vector3D(3, 0, 0),    // from +x, passes through mid-theta
        Vector3D(0, 3, 0),    // from +y, passes through mid-theta
        Vector3D(-10, 0, 0),  // from -x, outside phi cut [0,pi]
        Vector3D(0, 0, 10),   // from +z, near theta=0
        Vector3D(0, 0, -10),  // from -z, near theta=pi
        Vector3D(2, 2, 4),    // diagonal, likely in cut region
    };
    Vector3D directions[] = {
        Vector3D(0, 0, 1),
        Vector3D(0, 0, 1),
        Vector3D(1, 0, 0),
        Vector3D(0, 0, -1),
        Vector3D(0, 0, 1),
        Vector3D(-1, -1, -1),
    };
    for(int i = 0; i < 6; ++i) {
        Vector3D dir = directions[i];
        double mag = std::sqrt(dir.GetX()*dir.GetX() + dir.GetY()*dir.GetY() + dir.GetZ()*dir.GetZ());
        dir.SetCartesianCoordinates(dir.GetX()/mag, dir.GetY()/mag, dir.GetZ()/mag);
        auto hits_orig = sphere_cut->Intersections(origins[i], dir);
        auto hits_load = loaded->Intersections(origins[i], dir);
        EXPECT_EQ(hits_orig.size(), hits_load.size())
            << "Intersection count mismatch for sphere with cuts, ray " << i
            << " (orig=" << hits_orig.size() << " loaded=" << hits_load.size() << ")";
    }
}

TEST(SerializationRoundtrip, Cone) {
    auto cone = std::make_shared<Cone>(1.0, 3.0, 0.5, 2.0, 5.0);
    auto loaded = RoundtripGeometry(cone);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*cone, *loaded);
}

TEST(SerializationRoundtrip, Trd) {
    auto trd = std::make_shared<Trd>(1.0, 2.0, 1.5, 2.5, 3.0);
    auto loaded = RoundtripGeometry(trd);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*trd, *loaded);
}

TEST(SerializationRoundtrip, Polycone) {
    auto polycone = std::make_shared<Polycone>(
        std::vector<double>{-5.0, 0.0, 5.0},
        std::vector<double>{0.0, 0.0, 0.0},
        std::vector<double>{3.0, 5.0, 3.0}
    );
    auto loaded = RoundtripGeometry(polycone);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*polycone, *loaded);
}

TEST(SerializationRoundtrip, Polyhedra) {
    auto polyhedra = std::make_shared<Polyhedra>(
        6, 0.0,
        std::vector<double>{-5.0, 0.0, 5.0},
        std::vector<double>{0.0, 0.0, 0.0},
        std::vector<double>{3.0, 5.0, 3.0}
    );
    auto loaded = RoundtripGeometry(polyhedra);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*polyhedra, *loaded);
}

TEST(SerializationRoundtrip, TorusFull) {
    auto torus_full = std::make_shared<Torus>(10.0, 3.0, 1.0);
    auto loaded = RoundtripGeometry(torus_full);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*torus_full, *loaded);
}

TEST(SerializationRoundtrip, TorusWithPhiCut) {
    // Torus with phi cut -- tests that has_phi_cut_ is
    // recomputed on deserialization
    auto torus_cut = std::make_shared<Torus>(10.0, 3.0, 1.0, 0.0, M_PI);
    auto loaded = RoundtripGeometry(torus_cut);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*torus_cut, *loaded);

    // Verify intersection behavior matches after roundtrip.
    // Shoot rays through the torus from multiple angles.
    Vector3D origins[] = {
        Vector3D(10, 0, -10),
        Vector3D(0, 10, -10),
        Vector3D(-10, 0, -10),
        Vector3D(10, 0, 10),
    };
    Vector3D dir(0, 0, 1);
    Vector3D dir_neg(0, 0, -1);
    for(int i = 0; i < 4; ++i) {
        int orig_count = CountIntersections(*torus_cut, origins[i], dir);
        int load_count = CountIntersections(*loaded, origins[i], dir);
        EXPECT_EQ(orig_count, load_count)
            << "Intersection count mismatch for torus with phi cut, ray " << i;
    }
    for(int i = 0; i < 4; ++i) {
        int orig_count = CountIntersections(*torus_cut, origins[i], dir_neg);
        int load_count = CountIntersections(*loaded, origins[i], dir_neg);
        EXPECT_EQ(orig_count, load_count)
            << "Intersection count mismatch (neg dir) for torus with phi cut, ray " << i;
    }
}

TEST(SerializationRoundtrip, BooleanGeometry) {
    // Subtraction: box minus cylinder
    auto box = std::make_shared<Box>(4.0, 4.0, 4.0);
    auto cyl = std::make_shared<Cylinder>(1.0, 0.0, 4.0);
    std::shared_ptr<const Geometry> box_base = box;
    std::shared_ptr<const Geometry> cyl_base = cyl;
    auto bool_geo = std::make_shared<BooleanGeometry>(
        BooleanOperation::SUBTRACTION,
        box_base,
        cyl_base
    );
    auto loaded = RoundtripGeometry(bool_geo);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*bool_geo, *loaded);
}

TEST(SerializationRoundtrip, EllipticalTube) {
    auto et = std::make_shared<EllipticalTube>(3.0, 2.0, 5.0);
    auto loaded = RoundtripGeometry(et);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*et, *loaded);
}

TEST(SerializationRoundtrip, CutTube) {
    auto ct = std::make_shared<CutTube>(1.0, 5.0, 8.0,
        siren::math::Vector3D(0, -0.2, -0.98), siren::math::Vector3D(0.1, 0, 0.995));
    auto loaded = RoundtripGeometry(ct);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*ct, *loaded);
}

TEST(SerializationRoundtrip, Trap) {
    auto trap = std::make_shared<Trap>(5.0, 0.1, 0.2, 3.0, 4.0, 5.0, 0.15, 2.0, 3.0, 4.0, 0.1);
    auto loaded = RoundtripGeometry(trap);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*trap, *loaded);
}

TEST(SerializationRoundtrip, Ellipsoid) {
    auto ell = std::make_shared<Ellipsoid>(5.0, 3.0, 4.0, -2.0, 3.0);
    auto loaded = RoundtripGeometry(ell);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*ell, *loaded);
}

TEST(SerializationRoundtrip, Para) {
    auto para = std::make_shared<Para>(4.0, 3.0, 5.0, 0.3, 0.2, 0.5);
    auto loaded = RoundtripGeometry(para);
    ASSERT_NE(loaded, nullptr);
    EXPECT_EQ(*para, *loaded);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
