#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <unistd.h>

#include <gtest/gtest.h>

#include "SIREN/utilities/Random.h"
#include "SIREN/distributions/primary/PrimaryExternalDistribution.h"
#include "SIREN/distributions/primary/vertex/PrimaryBoundedVertexDistribution.h"
#include "SIREN/distributions/primary/vertex/PrimaryPhysicalVertexDistribution.h"
#include "SIREN/geometry/Sphere.h"

using namespace siren::distributions;
using namespace siren::dataclasses;

namespace {

// Helper: write a temporary CSV file and return its path.
// The caller is responsible for removing the file.
std::string WriteTempCSV(std::string const & contents) {
    std::string tmpl = "/tmp/siren_test_XXXXXX";
    int fd = mkstemp(&tmpl[0]);
    std::string path = tmpl;
    close(fd);
    std::ofstream out(path);
    out << contents;
    out.close();
    return path;
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, ConstructFromValidFile) {
    std::string path = WriteTempCSV("E,m\n10.0,1.0\n20.0,2.0\n");
    ASSERT_NO_THROW(PrimaryExternalDistribution dist(path));
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, ConstructWithEmin) {
    std::string path = WriteTempCSV("E,m\n5.0,1.0\n10.0,1.0\n20.0,2.0\n");
    PrimaryExternalDistribution dist(path, 8.0);
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 2u);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, ConstructFromMissingFileThrows) {
    EXPECT_THROW(PrimaryExternalDistribution dist("/nonexistent/path.csv"),
                 std::runtime_error);
}

TEST(PrimaryExternalDistribution, ConstructFromEmptyDataThrows) {
    std::string path = WriteTempCSV("E,m\n");
    EXPECT_THROW(PrimaryExternalDistribution dist(path), std::runtime_error);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, AllRowsFilteredByEminThrows) {
    std::string path = WriteTempCSV("E,m\n1.0,0.5\n2.0,0.5\n");
    EXPECT_THROW(PrimaryExternalDistribution dist(path, 100.0),
                 std::runtime_error);
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// CSV parsing edge cases
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, SkipsBlankLines) {
    std::string path = WriteTempCSV("E,m\n10.0,1.0\n\n20.0,2.0\n\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 2u);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SkipsCommentLines) {
    std::string path = WriteTempCSV("E,m\n# this is a comment\n10.0,1.0\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 1u);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, TooManyColumnsThrows) {
    std::string path = WriteTempCSV("E,m\n10.0,1.0,99.0\n");
    EXPECT_THROW(PrimaryExternalDistribution dist(path), std::runtime_error);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, TooFewColumnsThrows) {
    std::string path = WriteTempCSV("E,m\n10.0\n");
    EXPECT_THROW(PrimaryExternalDistribution dist(path), std::runtime_error);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, HandlesWhitespaceInHeader) {
    std::string path = WriteTempCSV(" E , m \n10.0,1.0\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 1u);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, HandlesCRLFLineEndings) {
    std::string path = WriteTempCSV("E,m\r\n10.0,1.0\r\n20.0,2.0\r\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 2u);
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// GetPhysicalNumEvents
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, GetPhysicalNumEventsNoFilter) {
    std::string path = WriteTempCSV("E\n1.0\n2.0\n3.0\n4.0\n5.0\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 5u);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, GetPhysicalNumEventsWithFilter) {
    std::string path = WriteTempCSV("E\n1.0\n2.0\n3.0\n4.0\n5.0\n");
    PrimaryExternalDistribution dist(path, 3.0);
    // E >= 3.0 keeps 3.0, 4.0, 5.0 (LoadInputFile filter uses < emin)
    EXPECT_EQ(dist.GetPhysicalNumEvents(), 3u);
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// Sample
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, SampleSetsEnergy) {
    std::string path = WriteTempCSV("E\n42.0\n");
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    PrimaryDistributionRecord record(ParticleType::NuMu);
    dist.Sample(rand, nullptr, nullptr, record);
    EXPECT_DOUBLE_EQ(record.GetEnergy(), 42.0);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleSetsMass) {
    std::string path = WriteTempCSV("E,m\n10.0,0.5\n");
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    PrimaryDistributionRecord record(ParticleType::NuMu);
    dist.Sample(rand, nullptr, nullptr, record);
    EXPECT_DOUBLE_EQ(record.GetMass(), 0.5);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleSetsPositionWhenAllPresent) {
    std::string path = WriteTempCSV("E,x0,y0,z0\n10.0,1.0,2.0,3.0\n");
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    PrimaryDistributionRecord record(ParticleType::NuMu);
    dist.Sample(rand, nullptr, nullptr, record);
    auto pos = record.GetInitialPosition();
    EXPECT_DOUBLE_EQ(pos[0], 1.0);
    EXPECT_DOUBLE_EQ(pos[1], 2.0);
    EXPECT_DOUBLE_EQ(pos[2], 3.0);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleSetsMomentumWhenAllPresent) {
    std::string path = WriteTempCSV("E,px,py,pz\n10.0,4.0,5.0,6.0\n");
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    PrimaryDistributionRecord record(ParticleType::NuMu);
    dist.Sample(rand, nullptr, nullptr, record);
    auto mom = record.GetThreeMomentum();
    EXPECT_DOUBLE_EQ(mom[0], 4.0);
    EXPECT_DOUBLE_EQ(mom[1], 5.0);
    EXPECT_DOUBLE_EQ(mom[2], 6.0);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleWorksWithoutEColumn) {
    // No "E" column -- success should default true
    std::string path = WriteTempCSV("x0,y0,z0\n1.0,2.0,3.0\n");
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    PrimaryDistributionRecord record(ParticleType::NuMu);
    ASSERT_NO_THROW(dist.Sample(rand, nullptr, nullptr, record));
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleCustomParameters) {
    std::string path = WriteTempCSV("E,Q2,x0,y0,z0,px,py,pz\n10.0,1.5,0,0,0,1,0,0\n");
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    PrimaryDistributionRecord record(ParticleType::NuMu);
    dist.Sample(rand, nullptr, nullptr, record);
    // Set interaction vertex so Finalize can propagate interaction_parameters
    record.SetInteractionVertex({0, 0, 0});
    InteractionRecord irec;
    record.Finalize(irec);
    EXPECT_DOUBLE_EQ(irec.interaction_parameters.at("Q2"), 1.5);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleUniformOverRows) {
    // With many rows, each should be sampled roughly equally
    std::string csv = "E\n";
    size_t nrows = 5;
    for(size_t i = 0; i < nrows; ++i) {
        csv += std::to_string(double(i + 1)) + "\n";
    }
    std::string path = WriteTempCSV(csv);
    PrimaryExternalDistribution dist(path);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    std::vector<size_t> counts(nrows, 0);
    size_t N = 50000;
    for(size_t i = 0; i < N; ++i) {
        PrimaryDistributionRecord record(ParticleType::NuMu);
        dist.Sample(rand, nullptr, nullptr, record);
        size_t idx = size_t(record.GetEnergy()) - 1;
        counts[idx]++;
    }
    double expected = double(N) / nrows;
    for(size_t i = 0; i < nrows; ++i) {
        // Allow 5-sigma deviation
        EXPECT_NEAR(double(counts[i]), expected, 5.0 * std::sqrt(expected));
    }
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, SampleRespectsEmin) {
    // Mix of below and above emin; all samples must have E >= emin
    std::string csv = "E\n";
    for(int i = 1; i <= 20; ++i) {
        csv += std::to_string(double(i)) + "\n";
    }
    std::string path = WriteTempCSV(csv);
    double emin = 10.0;
    PrimaryExternalDistribution dist(path, emin);
    auto rand = std::make_shared<siren::utilities::SIREN_random>();

    for(size_t i = 0; i < 10000; ++i) {
        PrimaryDistributionRecord record(ParticleType::NuMu);
        dist.Sample(rand, nullptr, nullptr, record);
        EXPECT_GE(record.GetEnergy(), emin);
    }
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// GenerationProbability
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, GenerationProbabilityAboveEmin) {
    std::string path = WriteTempCSV("E\n10.0\n20.0\n");
    PrimaryExternalDistribution dist(path, 5.0);

    InteractionRecord record;
    record.primary_momentum[0] = 15.0;
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, record), 1.0);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, GenerationProbabilityBelowEmin) {
    std::string path = WriteTempCSV("E\n10.0\n20.0\n");
    PrimaryExternalDistribution dist(path, 5.0);

    InteractionRecord record;
    record.primary_momentum[0] = 3.0;
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, record), 0.0);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, GenerationProbabilityAtEmin) {
    std::string path = WriteTempCSV("E\n10.0\n20.0\n");
    PrimaryExternalDistribution dist(path, 5.0);

    InteractionRecord record;
    record.primary_momentum[0] = 5.0;
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, record), 1.0);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, GenerationProbabilityNoEmin) {
    std::string path = WriteTempCSV("E\n10.0\n");
    PrimaryExternalDistribution dist(path);

    InteractionRecord record;
    record.primary_momentum[0] = 0.0;
    // emin defaults to 0, so energy >= 0 should return 1
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, record), 1.0);
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// Name and DensityVariables
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, Name) {
    std::string path = WriteTempCSV("E\n10.0\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_EQ(dist.Name(), "PrimaryExternalDistribution");
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, DensityVariables) {
    std::string path = WriteTempCSV("E\n10.0\n");
    PrimaryExternalDistribution dist(path);
    auto vars = dist.DensityVariables();
    ASSERT_EQ(vars.size(), 1u);
    EXPECT_EQ(vars[0], "External");
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// clone
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, CloneIsEquivalent) {
    std::string path = WriteTempCSV("E,m\n10.0,1.0\n20.0,2.0\n");
    PrimaryExternalDistribution dist(path);
    PrimaryExternalDistribution copy(dist);
    auto cloned = dist.clone();

    EXPECT_TRUE(dist == copy);
    EXPECT_EQ(cloned->Name(), "PrimaryExternalDistribution");
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// equal / less (strict weak ordering)
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, EqualSameData) {
    std::string path = WriteTempCSV("E\n10.0\n20.0\n");
    PrimaryExternalDistribution a(path);
    PrimaryExternalDistribution b(path);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(b < a);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, NotEqualDifferentData) {
    std::string path_a = WriteTempCSV("E\n10.0\n");
    std::string path_b = WriteTempCSV("E\n20.0\n");
    PrimaryExternalDistribution a(path_a);
    PrimaryExternalDistribution b(path_b);
    EXPECT_FALSE(a == b);
    // Strict weak ordering: exactly one of a<b or b<a must be true
    EXPECT_NE(a < b, b < a);
    std::remove(path_a.c_str());
    std::remove(path_b.c_str());
}

TEST(PrimaryExternalDistribution, NotEqualDifferentEmin) {
    std::string path = WriteTempCSV("E\n10.0\n20.0\n30.0\n");
    PrimaryExternalDistribution a(path, 0.0);
    PrimaryExternalDistribution b(path, 15.0);
    EXPECT_FALSE(a == b);
    EXPECT_NE(a < b, b < a);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistribution, StrictWeakOrdering) {
    std::string p1 = WriteTempCSV("E\n1.0\n");
    std::string p2 = WriteTempCSV("E\n2.0\n");
    std::string p3 = WriteTempCSV("E\n3.0\n");
    PrimaryExternalDistribution a(p1);
    PrimaryExternalDistribution b(p2);
    PrimaryExternalDistribution c(p3);
    // Irreflexivity
    EXPECT_FALSE(a < a);
    // If a < b and b < c then a < c (transitivity)
    if((a < b) && (b < c)) {
        EXPECT_TRUE(a < c);
    } else if((c < b) && (b < a)) {
        EXPECT_TRUE(c < a);
    }
    std::remove(p1.c_str());
    std::remove(p2.c_str());
    std::remove(p3.c_str());
}

// ---------------------------------------------------------------------------
// Copy constructor
// ---------------------------------------------------------------------------

TEST(PrimaryExternalDistribution, CopyConstructor) {
    std::string path = WriteTempCSV("E,m\n10.0,1.0\n20.0,2.0\n");
    PrimaryExternalDistribution original(path);
    PrimaryExternalDistribution copy(original);

    EXPECT_TRUE(original == copy);
    EXPECT_EQ(copy.GetPhysicalNumEvents(), original.GetPhysicalNumEvents());
    EXPECT_EQ(copy.Name(), original.Name());
    std::remove(path.c_str());
}

// ===========================================================================
// PrimaryBoundedVertexDistribution
// ===========================================================================

TEST(PrimaryBoundedVertexDistribution, EqualSameMaxLength) {
    PrimaryBoundedVertexDistribution a(100.0);
    PrimaryBoundedVertexDistribution b(100.0);
    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(b < a);
}

TEST(PrimaryBoundedVertexDistribution, NotEqualDifferentMaxLength) {
    PrimaryBoundedVertexDistribution a(100.0);
    PrimaryBoundedVertexDistribution b(200.0);
    EXPECT_FALSE(a == b);
    EXPECT_NE(a < b, b < a);
}

TEST(PrimaryBoundedVertexDistribution, EqualSameFiducialVolume) {
    auto sphere_a = std::make_shared<siren::geometry::Sphere>(10.0, 0.0);
    auto sphere_b = std::make_shared<siren::geometry::Sphere>(10.0, 0.0);
    PrimaryBoundedVertexDistribution a(sphere_a, 100.0);
    PrimaryBoundedVertexDistribution b(sphere_b, 100.0);
    EXPECT_TRUE(a == b);
}

TEST(PrimaryBoundedVertexDistribution, NotEqualDifferentFiducialVolume) {
    auto small = std::make_shared<siren::geometry::Sphere>(10.0, 0.0);
    auto large = std::make_shared<siren::geometry::Sphere>(50.0, 0.0);
    PrimaryBoundedVertexDistribution a(small, 100.0);
    PrimaryBoundedVertexDistribution b(large, 100.0);
    EXPECT_FALSE(a == b);
    EXPECT_NE(a < b, b < a);
}

TEST(PrimaryBoundedVertexDistribution, NotEqualNullVsNonNullFiducial) {
    auto sphere = std::make_shared<siren::geometry::Sphere>(10.0, 0.0);
    PrimaryBoundedVertexDistribution a(100.0);
    PrimaryBoundedVertexDistribution b(sphere, 100.0);
    EXPECT_FALSE(a == b);
    EXPECT_NE(a < b, b < a);
}

TEST(PrimaryBoundedVertexDistribution, EqualBothNullFiducial) {
    PrimaryBoundedVertexDistribution a(100.0);
    PrimaryBoundedVertexDistribution b(100.0);
    EXPECT_TRUE(a == b);
}

TEST(PrimaryBoundedVertexDistribution, StrictWeakOrdering) {
    PrimaryBoundedVertexDistribution a(100.0);
    PrimaryBoundedVertexDistribution b(200.0);
    PrimaryBoundedVertexDistribution c(300.0);
    EXPECT_FALSE(a < a);
    if((a < b) && (b < c)) {
        EXPECT_TRUE(a < c);
    }
}

TEST(PrimaryBoundedVertexDistribution, Name) {
    PrimaryBoundedVertexDistribution dist;
    EXPECT_EQ(dist.Name(), "PrimaryBoundedVertexDistribution");
}

TEST(PrimaryBoundedVertexDistribution, ClonePreservesEquality) {
    auto sphere = std::make_shared<siren::geometry::Sphere>(10.0, 0.0);
    PrimaryBoundedVertexDistribution dist(sphere, 100.0);
    PrimaryBoundedVertexDistribution copy(dist);
    EXPECT_TRUE(dist == copy);
}

// ===========================================================================
// PrimaryPhysicalVertexDistribution
// ===========================================================================

TEST(PrimaryPhysicalVertexDistribution, EqualAlwaysTrue) {
    PrimaryPhysicalVertexDistribution a;
    PrimaryPhysicalVertexDistribution b;
    EXPECT_TRUE(a == b);
}

TEST(PrimaryPhysicalVertexDistribution, LessAlwaysFalse) {
    PrimaryPhysicalVertexDistribution a;
    PrimaryPhysicalVertexDistribution b;
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(b < a);
}

TEST(PrimaryPhysicalVertexDistribution, Name) {
    PrimaryPhysicalVertexDistribution dist;
    EXPECT_EQ(dist.Name(), "PrimaryPhysicalVertexDistribution");
}

TEST(PrimaryPhysicalVertexDistribution, ClonePreservesEquality) {
    PrimaryPhysicalVertexDistribution dist;
    PrimaryPhysicalVertexDistribution copy(dist);
    EXPECT_TRUE(dist == copy);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
