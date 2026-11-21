
// Degenerate trajectory-direction behavior for DetectorModel depth queries.
// DetectorModel computes direction as (p1 - p0); at Earth-scale coordinates a
// sub-micrometer separation cancels catastrophically (coincident points give the
// zero vector). All depth queries must stay finite, non-negative, and not abort.
//
// NOTE: the unit-direction assert (std::abs(1.0 - std::abs(dot)) < 1e-6) only
// fires in an assertions-enabled build (Debug / RelWithDebInfo), not under NDEBUG;
// the finite/non-NaN return-value asserts below lock in behavior in either build.

#include <cmath>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/geometry/Geometry.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/dataclasses/Particle.h"

using namespace siren::math;
using namespace siren::geometry;
using namespace siren::detector;
using namespace siren::dataclasses;

namespace {

// distance_threshold in DetectorModel is 1e-5 m; pick an offset well below it
// but large enough to be exactly representable next to an Earth-scale x value.
constexpr double kEarthScaleX = 6.371e6;   // m, ~Earth radius
constexpr double kSubThresholdOffset = 1e-7; // m, << distance_threshold (1e-5)

bool IsFinite(double x) {
    return std::isfinite(x) && !std::isnan(x);
}

} // namespace

// Distinct but sub-threshold points with non-empty targets: must short-circuit
// to distance/total_decay_length instead of normalizing the corrupted direction.
TEST(DetectorModelDegenerateDirection, InteractionDepthSubThresholdWithTargets)
{
    DetectorModel A; // infinite vacuum sphere, no files

    Vector3D p0;
    p0.SetCartesianCoordinates(kEarthScaleX, 0.0, 0.0);
    Vector3D p1;
    p1.SetCartesianCoordinates(kEarthScaleX + kSubThresholdOffset, 0.0, 0.0);

    // p0 != p1 (early guard must not fire), yet separation is below threshold.
    ASSERT_TRUE(p0 != p1);
    double distance = (p1 - p0).magnitude();
    ASSERT_GT(distance, 0.0);
    ASSERT_LE(distance, 1e-5);

    std::vector<ParticleType> targets = {ParticleType::Nucleon};
    std::vector<double> total_cross_sections = {1.0e-27}; // cm^2, arbitrary
    double total_decay_length = 1.0e3; // m, arbitrary finite, nonzero

    double depth = std::numeric_limits<double>::quiet_NaN();
    // Must NOT abort on assert(std::abs(1.0 - std::abs(dot)) < 1e-6).
    ASSERT_NO_THROW(depth = A.GetInteractionDepthInCGS(
        DetectorPosition(p0), DetectorPosition(p1),
        targets, total_cross_sections, total_decay_length));

    EXPECT_TRUE(IsFinite(depth)) << "InteractionDepth must be finite, got " << depth;

    // Sub-threshold result is the decay-only limit (continuous with targets.empty()).
    double expected = distance / total_decay_length;
    EXPECT_NEAR(expected, depth, std::abs(expected) * 1e-9 + 1e-30);
    EXPECT_GE(depth, 0.0);
}

// With empty targets, GetInteractionDepthInCGS takes the decay-only branch.
// Confirm the sub-threshold result is continuous with that branch and finite.
TEST(DetectorModelDegenerateDirection, InteractionDepthSubThresholdDecayOnly)
{
    DetectorModel A;

    Vector3D p0;
    p0.SetCartesianCoordinates(kEarthScaleX, 0.0, 0.0);
    Vector3D p1;
    p1.SetCartesianCoordinates(kEarthScaleX + kSubThresholdOffset, 0.0, 0.0);

    double distance = (p1 - p0).magnitude();
    std::vector<ParticleType> targets; // empty: decay-only branch
    std::vector<double> total_cross_sections;
    double total_decay_length = 2.5e2;

    double depth = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NO_THROW(depth = A.GetInteractionDepthInCGS(
        DetectorPosition(p0), DetectorPosition(p1),
        targets, total_cross_sections, total_decay_length));

    EXPECT_TRUE(IsFinite(depth));
    EXPECT_NEAR(distance / total_decay_length, depth,
                std::abs(distance / total_decay_length) * 1e-9 + 1e-30);
}

// Column depth on a sub-threshold but distinct segment normalizes a tiny (~1e-7)
// direction; must complete without aborting and return finite, non-negative.
TEST(DetectorModelDegenerateDirection, ColumnDepthSubThresholdIsFinite)
{
    DetectorModel A;

    Vector3D p0;
    p0.SetCartesianCoordinates(kEarthScaleX, 0.0, 0.0);
    Vector3D p1;
    p1.SetCartesianCoordinates(kEarthScaleX + kSubThresholdOffset, 0.0, 0.0);

    double cd = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NO_THROW(cd = A.GetColumnDepthInCGS(DetectorPosition(p0), DetectorPosition(p1)));
    EXPECT_TRUE(IsFinite(cd)) << "ColumnDepth must be finite, got " << cd;
    EXPECT_GE(cd, 0.0);
}

// GetParticleColumnDepth on the same sub-threshold segment must also return
// finite per-target values without aborting.
TEST(DetectorModelDegenerateDirection, ParticleColumnDepthSubThresholdIsFinite)
{
    DetectorModel A;

    Vector3D p0;
    p0.SetCartesianCoordinates(kEarthScaleX, 0.0, 0.0);
    Vector3D p1;
    p1.SetCartesianCoordinates(kEarthScaleX + kSubThresholdOffset, 0.0, 0.0);

    std::vector<ParticleType> targets = {ParticleType::Nucleon, ParticleType::EMinus};

    std::vector<double> counts;
    Geometry::IntersectionList intersections = A.GetIntersections(
        DetectorPosition(p0), DetectorDirection(Vector3D(1.0, 0.0, 0.0)));
    ASSERT_NO_THROW(counts = A.GetParticleColumnDepth(
        intersections, DetectorPosition(p0), DetectorPosition(p1), targets));

    ASSERT_EQ(targets.size(), counts.size());
    for(double c : counts) {
        EXPECT_TRUE(IsFinite(c)) << "ParticleColumnDepth entry must be finite, got " << c;
        EXPECT_GE(c, 0.0);
    }
}

// Exact coincidence (p0 == p1) at Earth-scale: (p1 - p0) is the zero vector. The
// early guard plus Vector3D::normalize()'s zero-length guard must keep every depth
// query at 0, finite, no NaN.
TEST(DetectorModelDegenerateDirection, CoincidentPointsReturnZeroFinite)
{
    DetectorModel A;

    Vector3D v;
    v.SetCartesianCoordinates(kEarthScaleX, 1.2e6, -3.4e5);

    // Column depth between identical points is exactly zero.
    double cd = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NO_THROW(cd = A.GetColumnDepthInCGS(DetectorPosition(v), DetectorPosition(v)));
    EXPECT_TRUE(IsFinite(cd));
    EXPECT_DOUBLE_EQ(0.0, cd);

    // Interaction depth between identical points is exactly zero (early guard).
    std::vector<ParticleType> targets = {ParticleType::Nucleon};
    std::vector<double> total_cross_sections = {1.0e-27};
    double total_decay_length = 1.0e3;
    double depth = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NO_THROW(depth = A.GetInteractionDepthInCGS(
        DetectorPosition(v), DetectorPosition(v),
        targets, total_cross_sections, total_decay_length));
    EXPECT_TRUE(IsFinite(depth));
    EXPECT_DOUBLE_EQ(0.0, depth);

    // Particle column depth between identical points is all zeros.
    Geometry::IntersectionList intersections = A.GetIntersections(
        DetectorPosition(v), DetectorDirection(Vector3D(1.0, 0.0, 0.0)));
    std::vector<double> counts;
    ASSERT_NO_THROW(counts = A.GetParticleColumnDepth(
        intersections, DetectorPosition(v), DetectorPosition(v), targets));
    ASSERT_EQ(targets.size(), counts.size());
    for(double c : counts) {
        EXPECT_TRUE(IsFinite(c));
        EXPECT_DOUBLE_EQ(0.0, c);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
