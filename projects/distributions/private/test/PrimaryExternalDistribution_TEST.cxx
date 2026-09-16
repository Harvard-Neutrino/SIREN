#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <unistd.h>

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <sstream>
#include <cereal/archives/json.hpp>

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

TEST(PrimaryExternalDistribution, ConstructFromInMemoryRows) {
    PrimaryExternalDistribution dist(
        {"E", "x0", "y0", "z0"},
        {{10.0, 1.0, 2.0, 3.0}, {20.0, 4.0, 5.0, 6.0}},
        15.0);

    EXPECT_EQ(dist.GetPhysicalNumEvents(), 1u);
    EXPECT_TRUE(dist.RequiredVariables().empty());
    EXPECT_TRUE(dist.SetVariables().count(DistributionVariable::PrimaryEnergy));
    EXPECT_TRUE(dist.SetVariables().count(DistributionVariable::InitialPosition));
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

// ===========================================================================
// Segment mode ("length" column)
// ===========================================================================

namespace {

std::string SegmentCSV(double length) {
    std::ostringstream out;
    out << "E,px,py,pz,x0,y0,z0,m,length\n";
    out << "0.05,0.0,0.0,0.05,1.0,-2.0,0.5,0.0," << length << "\n";
    return out.str();
}

siren::dataclasses::InteractionRecord SampleOne(PrimaryExternalDistribution const & dist, unsigned int seed = 3) {
    auto rand = std::make_shared<siren::utilities::SIREN_random>(seed);
    PrimaryDistributionRecord record(ParticleType::Gamma);
    dist.Sample(rand, nullptr, nullptr, record);
    siren::dataclasses::InteractionRecord ir;
    record.Finalize(ir);
    return ir;
}

} // namespace

TEST(PrimaryExternalDistributionSegments, DeclaresVertexAndLongitudinalDensity) {
    std::string path = WriteTempCSV(SegmentCSV(0.2));
    PrimaryExternalDistribution metadata(path);
    EXPECT_FALSE(metadata.ProvidesExternalBounds());
    EXPECT_FALSE(metadata.SetVariables().count(DistributionVariable::InteractionVertex));
    PrimaryExternalDistribution dist(path);
    dist.SetSegmentColumn("length");
    EXPECT_TRUE(dist.ProvidesExternalBounds());
    EXPECT_TRUE(dist.SetVariables().count(DistributionVariable::InteractionVertex));
    EXPECT_TRUE(dist.SetVariables().count(DistributionVariable::InitialPosition));
    std::vector<std::string> density = dist.DensityVariables();
    ASSERT_EQ(density.size(), 2u);
    EXPECT_EQ(density[1], "PrimaryPositionLongitudinal");
    EXPECT_EQ(dist.PhysicalDensityVariables(), std::vector<std::string>{"External"});
    EXPECT_TRUE(dist.PhysicalDensityDiffers());
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, PointTablesAreUnchanged) {
    std::string path = WriteTempCSV("E,px,py,pz,x0,y0,z0,m\n0.05,0.0,0.0,0.05,1.0,-2.0,0.5,0.0\n");
    PrimaryExternalDistribution dist(path);
    EXPECT_FALSE(dist.ProvidesExternalBounds());
    EXPECT_FALSE(dist.SetVariables().count(DistributionVariable::InteractionVertex));
    EXPECT_EQ(dist.DensityVariables(), std::vector<std::string>{"External"});
    EXPECT_FALSE(dist.PhysicalDensityDiffers());
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, VertexAlongSegmentAndBounds) {
    double const L = 0.3;
    std::string path = WriteTempCSV(SegmentCSV(L));
    PrimaryExternalDistribution dist(path);
    dist.SetSegmentColumn("length");
    for (unsigned int seed = 1; seed < 40; ++seed) {
        siren::dataclasses::InteractionRecord ir = SampleOne(dist, seed);
        double s = ir.interaction_vertex[2] - 0.5;
        EXPECT_DOUBLE_EQ(ir.interaction_vertex[0], 1.0);
        EXPECT_DOUBLE_EQ(ir.interaction_vertex[1], -2.0);
        EXPECT_GE(s, 0.0);
        EXPECT_LE(s, L);
        EXPECT_DOUBLE_EQ(ir.interaction_parameters.at("length"), L);
        EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, ir), 1.0 / L);
        EXPECT_DOUBLE_EQ(dist.PhysicalDensity(nullptr, nullptr, ir), 1.0);
        auto bounds = dist.InjectionBounds(nullptr, nullptr, ir);
        EXPECT_DOUBLE_EQ(std::get<0>(bounds).GetZ(), 0.5);
        EXPECT_DOUBLE_EQ(std::get<1>(bounds).GetZ(), 0.5 + L);
        EXPECT_DOUBLE_EQ(std::get<1>(bounds).GetX(), 1.0);
    }
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, RejectsInconsistentHeaders) {
    for (std::string header : {"E,px,py,pz,m,length\n0.05,0,0,0.05,0,0.1\n",
                               "E,x0,y0,z0,m,length\n0.05,0,0,0,0,0.1\n",
                               "E,px,py,pz,x0,y0,z0,x,y,z,m,length\n0.05,0,0,0.05,0,0,0,0,0,0,0,0.1\n",
                               "E,px,py,pz,x0,y0,z0,m,length\n0.05,0,0,0.05,0,0,0,0,0.0\n",
                               "E,px,py,pz,x0,y0,z0,m,length\n0.05,0,0,0.05,0,0,0,0,-1.0\n"}) {
        std::string path = WriteTempCSV(header);
        EXPECT_THROW({ PrimaryExternalDistribution dist(path); dist.SetSegmentColumn("length"); }, std::runtime_error) << header;
        std::remove(path.c_str());
    }
}

TEST(PrimaryExternalDistributionSegments, SegmentModeSurvivesSerialization) {
    std::string path = WriteTempCSV(SegmentCSV(0.25));
    auto segment = std::make_shared<PrimaryExternalDistribution>(path);
    segment->SetSegmentColumn("length");
    std::shared_ptr<PrimaryInjectionDistribution> original = segment;
    std::stringstream buffer;
    {
        cereal::JSONOutputArchive archive(buffer);
        archive(cereal::make_nvp("Distribution", original));
    }
    std::shared_ptr<PrimaryInjectionDistribution> restored;
    {
        cereal::JSONInputArchive archive(buffer);
        archive(cereal::make_nvp("Distribution", restored));
    }
    auto dist = std::dynamic_pointer_cast<PrimaryExternalDistribution>(restored);
    ASSERT_NE(dist, nullptr);
    EXPECT_TRUE(dist->ProvidesExternalBounds());
    EXPECT_TRUE(dist->SetVariables().count(DistributionVariable::InteractionVertex));
    siren::dataclasses::InteractionRecord ir = SampleOne(*dist, 5);
    EXPECT_DOUBLE_EQ(dist->GenerationProbability(nullptr, nullptr, ir), 4.0);
    auto bounds = dist->InjectionBounds(nullptr, nullptr, ir);
    EXPECT_DOUBLE_EQ(std::get<1>(bounds).GetZ(), 0.75);
    EXPECT_TRUE(*original == *restored);
    EXPECT_EQ(buffer.str().find("SegmentColumn") == std::string::npos, false);
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, OffSupportRecordsHaveZeroDensity) {
    std::string path = WriteTempCSV(SegmentCSV(0.12));
    PrimaryExternalDistribution dist(path);
    dist.SetSegmentColumn("length");
    siren::dataclasses::InteractionRecord ir = SampleOne(dist, 7);
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, ir), 1.0 / 0.12);
    siren::dataclasses::InteractionRecord moved = ir;
    moved.interaction_vertex[2] = 0.5 + 0.13;  // beyond the end
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, moved), 0.0);
    moved = ir; moved.interaction_vertex[2] = 0.5 - 0.01;  // upstream
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, moved), 0.0);
    moved = ir; moved.interaction_vertex[0] = 1.001;  // transverse
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, moved), 0.0);
    moved = ir; moved.primary_initial_position[2] = 0.49;  // different start
    EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, moved), 0.0);
    // A shorter table over the same source never claims the longer one's tail.
    std::string short_path = WriteTempCSV(SegmentCSV(0.06));
    PrimaryExternalDistribution shorter(short_path);
    shorter.SetSegmentColumn("length");
    siren::dataclasses::InteractionRecord tail = ir;
    tail.interaction_vertex[2] = 0.5 + 0.09;
    EXPECT_DOUBLE_EQ(shorter.GenerationProbability(nullptr, nullptr, tail), 0.0);
    tail.interaction_vertex[2] = 0.5 + 0.03;
    EXPECT_DOUBLE_EQ(shorter.GenerationProbability(nullptr, nullptr, tail), 1.0 / 0.06);
    std::remove(path.c_str());
    std::remove(short_path.c_str());
}

TEST(PrimaryExternalDistributionSegments, OldArchivesWithLengthColumnStayPointTables) {
    // A version-1 archive written before segment mode existed.
    std::string path = WriteTempCSV(SegmentCSV(0.12));
    std::shared_ptr<PrimaryInjectionDistribution> point = std::make_shared<PrimaryExternalDistribution>(path);
    std::stringstream buffer;
    {
        cereal::JSONOutputArchive archive(buffer);
        archive(cereal::make_nvp("Distribution", point));
    }
    std::string text = buffer.str();
    // Downgrade the stored class version and drop the new field.
    size_t sc = text.find("\"SegmentColumn\"");
    ASSERT_NE(sc, std::string::npos);
    size_t line_start = text.rfind('\n', sc);
    size_t line_end = text.find('\n', sc);
    text.erase(line_start, line_end - line_start);
    // remove trailing comma left on the previous field
    size_t comma = text.rfind(',', line_start);
    text.erase(comma, 1);
    size_t ver = text.find("\"cereal_class_version\": 2");
    ASSERT_NE(ver, std::string::npos);
    text.replace(ver, std::string("\"cereal_class_version\": 2").size(), "\"cereal_class_version\": 1");
    std::stringstream input(text);
    std::shared_ptr<PrimaryInjectionDistribution> loaded;
    {
        cereal::JSONInputArchive archive(input);
        archive(cereal::make_nvp("Distribution", loaded));
    }
    auto restored = std::dynamic_pointer_cast<PrimaryExternalDistribution>(loaded);
    ASSERT_NE(restored, nullptr);
    EXPECT_FALSE(restored->ProvidesExternalBounds());
    EXPECT_EQ(restored->GetSegmentColumn(), "");
    EXPECT_FALSE(restored->SetVariables().count(DistributionVariable::InteractionVertex));
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, NonintegralCachedRowIsRejected) {
    std::string path = WriteTempCSV(SegmentCSV(0.12));
    PrimaryExternalDistribution dist(path);
    dist.SetSegmentColumn("length");
    siren::dataclasses::InteractionRecord ir = SampleOne(dist, 9);
    for (double bad : {0.75, -0.5, 1.0, 17.0, 1e18,
                       std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
        siren::dataclasses::InteractionRecord probe = ir;
        probe.interaction_parameters["PrimaryExternalDistribution_row"] = bad;
        EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, probe), 0.0) << bad;
        auto bounds = dist.InjectionBounds(nullptr, nullptr, probe);
        EXPECT_DOUBLE_EQ((std::get<1>(bounds) - std::get<0>(bounds)).magnitude(), 0.0) << bad;
    }
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, RejectedSetSegmentColumnPreservesState) {
    std::string path = WriteTempCSV(SegmentCSV(0.12));
    PrimaryExternalDistribution dist(path);
    dist.SetSegmentColumn("length");
    PrimaryExternalDistribution reference(dist);
    for (std::string bad : {"absent", "m", "E", "x0"}) {
        EXPECT_THROW(dist.SetSegmentColumn(bad), std::runtime_error) << bad;
        EXPECT_TRUE(dist == reference) << bad;
        EXPECT_TRUE(dist.ProvidesExternalBounds()) << bad;
        EXPECT_EQ(dist.GetSegmentColumn(), "length") << bad;
        siren::dataclasses::InteractionRecord ir = SampleOne(dist, 11);
        EXPECT_DOUBLE_EQ(dist.GenerationProbability(nullptr, nullptr, ir), 1.0 / 0.12) << bad;
    }
    // A successful update still switches semantics.
    dist.SetSegmentColumn("");
    EXPECT_FALSE(dist.ProvidesExternalBounds());
    std::remove(path.c_str());
}

TEST(PrimaryExternalDistributionSegments, PooledTablesMayNameTheirColumnsDifferently) {
    std::string path = WriteTempCSV(SegmentCSV(0.12));
    std::string other_csv = SegmentCSV(0.12);
    other_csv.replace(other_csv.find("length"), 6, "step_length");
    std::string other_path = WriteTempCSV(other_csv);
    PrimaryExternalDistribution a(path);
    a.SetSegmentColumn("length");
    PrimaryExternalDistribution b(other_path);
    b.SetSegmentColumn("step_length");
    // Each table evaluates the other's records without demanding its own
    // column name: identical geometry means on-support, so both report 1/L.
    siren::dataclasses::InteractionRecord from_a = SampleOne(a, 13);
    siren::dataclasses::InteractionRecord from_b = SampleOne(b, 13);
    EXPECT_DOUBLE_EQ(a.GenerationProbability(nullptr, nullptr, from_b), 1.0 / 0.12);
    EXPECT_DOUBLE_EQ(b.GenerationProbability(nullptr, nullptr, from_a), 1.0 / 0.12);
    auto bounds = b.InjectionBounds(nullptr, nullptr, from_a);
    EXPECT_DOUBLE_EQ(std::get<1>(bounds).GetZ(), 0.5 + 0.12);
    std::remove(path.c_str());
    std::remove(other_path.c_str());
}

TEST(PrimaryExternalDistributionSegments, RowLayoutMismatchIgnoresLengthsWeightsAndColumnNames) {
    std::vector<std::string> keys = {"E", "px", "py", "pz", "x0", "y0", "z0", "m", "length", "weight"};
    std::vector<std::vector<double>> rows = {
        {0.05, 0.0, 0.0, 0.05, 1.0, -2.0, 0.5, 0.0, 0.12, 1.0},
        {0.05, 0.0, 0.0, 0.05, 1.0, -1.9, 0.5, 0.0, 0.12, 2.0}};
    PrimaryExternalDistribution reference(keys, rows);
    reference.SetSegmentColumn("length");
    // Identical content, a copy, different lengths, different physical
    // weights, explicit sampling weights, a renamed length column, and a
    // table without the weight column are all the same layout.
    EXPECT_EQ(reference.RowLayoutMismatch(reference), "");
    PrimaryExternalDistribution copy(reference);
    EXPECT_EQ(reference.RowLayoutMismatch(copy), "");
    std::vector<std::vector<double>> shorter = rows;
    shorter[0][8] = 0.06;
    shorter[1][9] = 7.0;
    PrimaryExternalDistribution other(keys, shorter, std::vector<double>{1.0, 9.0});
    other.SetSegmentColumn("length");
    EXPECT_EQ(reference.RowLayoutMismatch(other), "");
    EXPECT_EQ(other.RowLayoutMismatch(reference), "");
    std::vector<std::string> renamed = keys;
    renamed[8] = "step_length";
    PrimaryExternalDistribution named(renamed, rows);
    named.SetSegmentColumn("step_length");
    EXPECT_EQ(reference.RowLayoutMismatch(named), "");
    EXPECT_EQ(named.RowLayoutMismatch(reference), "");
    std::vector<std::string> unweighted(keys.begin(), keys.end() - 1);
    std::vector<std::vector<double>> unweighted_rows;
    for (auto row : rows) { row.pop_back(); unweighted_rows.push_back(row); }
    PrimaryExternalDistribution no_weight(unweighted, unweighted_rows);
    no_weight.SetSegmentColumn("length");
    EXPECT_EQ(reference.RowLayoutMismatch(no_weight), "");
}

TEST(PrimaryExternalDistributionSegments, RowLayoutMismatchReportsReorderingSubsetsAndColumns) {
    std::vector<std::string> keys = {"E", "px", "py", "pz", "x0", "y0", "z0", "m", "length"};
    std::vector<std::vector<double>> rows = {
        {0.05, 0.0, 0.0, 0.05, 1.0, -2.0, 0.5, 0.0, 0.12},
        {0.05, 0.0, 0.0, 0.05, 1.0, -1.9, 0.5, 0.0, 0.12}};
    PrimaryExternalDistribution reference(keys, rows);
    reference.SetSegmentColumn("length");
    PrimaryExternalDistribution reversed(keys, {rows[1], rows[0]});
    reversed.SetSegmentColumn("length");
    std::string mismatch = reference.RowLayoutMismatch(reversed);
    EXPECT_NE(mismatch.find("row 0"), std::string::npos) << mismatch;
    EXPECT_NE(mismatch.find("y0"), std::string::npos) << mismatch;
    PrimaryExternalDistribution subset(keys, {rows[0]});
    subset.SetSegmentColumn("length");
    EXPECT_NE(reference.RowLayoutMismatch(subset).find("row counts differ"), std::string::npos);
    EXPECT_NE(subset.RowLayoutMismatch(reference).find("row counts differ"), std::string::npos);
    std::vector<std::string> extra = keys;
    extra.push_back("material");
    std::vector<std::vector<double>> extra_rows = rows;
    for (auto & row : extra_rows) row.push_back(1.0);
    PrimaryExternalDistribution annotated(extra, extra_rows);
    annotated.SetSegmentColumn("length");
    EXPECT_NE(reference.RowLayoutMismatch(annotated).find("\"material\""), std::string::npos);
    // The layout compares row identity only: a point table over the same rows
    // whose "length" is metadata shares it (either table's segment column is
    // ignored), whichever side asks.
    PrimaryExternalDistribution point(keys, rows);
    EXPECT_EQ(reference.RowLayoutMismatch(point), "");
    EXPECT_EQ(point.RowLayoutMismatch(reference), "");
}

TEST(PrimaryExternalDistributionSegments, ForeignRowIsOffSupportBeforeItsSamplingDensity) {
    // A record from a two-row table evaluated by a one-row subset with
    // explicit sampling weights: support is decided first, so the answer is
    // zero density, not an out-of-range row error.
    std::vector<std::string> keys = {"E", "px", "py", "pz", "x0", "y0", "z0", "m", "length"};
    std::vector<std::vector<double>> rows = {
        {0.05, 0.0, 0.0, 0.05, 1.0, -2.0, 0.5, 0.0, 0.12},
        {0.05, 0.0, 0.0, 0.05, 1.0, -1.9, 0.5, 0.0, 0.12}};
    PrimaryExternalDistribution full(keys, rows, std::vector<double>{1.0, 1.0});
    full.SetSegmentColumn("length");
    PrimaryExternalDistribution subset(keys, {rows[0]}, std::vector<double>{1.0});
    subset.SetSegmentColumn("length");
    bool saw_row_one = false;
    for (unsigned int seed = 1; seed < 40 && !saw_row_one; ++seed) {
        siren::dataclasses::InteractionRecord ir = SampleOne(full, seed);
        if (ir.interaction_parameters.at("PrimaryExternalDistribution_row") != 1.0) continue;
        saw_row_one = true;
        EXPECT_DOUBLE_EQ(subset.GenerationProbability(nullptr, nullptr, ir), 0.0);
    }
    EXPECT_TRUE(saw_row_one);
}

TEST(PrimaryExternalDistributionSegments, SupportToleranceDoesNotExtendShortSegments) {
    // A 0.1 nm proposal must not claim the records of a 1 nm proposal beyond
    // its own end point: the tolerance is a few ulps of the coordinates, not
    // a fixed floor of 1e-9 m.
    std::vector<std::string> keys = {"E", "px", "py", "pz", "x0", "y0", "z0", "m", "length"};
    PrimaryExternalDistribution shorter(keys, {{0.05, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 1e-10}});
    shorter.SetSegmentColumn("length");
    PrimaryExternalDistribution longer(keys, {{0.05, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 1e-9}});
    longer.SetSegmentColumn("length");
    unsigned int inside = 0, outside = 0;
    for (unsigned int seed = 1; seed <= 200; ++seed) {
        siren::dataclasses::InteractionRecord ir = SampleOne(longer, seed);
        EXPECT_DOUBLE_EQ(longer.GenerationProbability(nullptr, nullptr, ir), 1e9);
        double density = shorter.GenerationProbability(nullptr, nullptr, ir);
        if (ir.interaction_vertex[2] <= 1e-10) {
            EXPECT_DOUBLE_EQ(density, 1e10);
            ++inside;
        } else {
            EXPECT_DOUBLE_EQ(density, 0.0);
            ++outside;
        }
    }
    EXPECT_GT(inside, 5u);
    EXPECT_GT(outside, 150u);
    // Far from the origin the tolerance follows the coordinate scale, so a
    // table's own records along an oblique short segment stay on support.
    PrimaryExternalDistribution far(keys, {{0.05, 0.01, 0.02, 0.02, -23.0, 4.0, 0.3, 0.0, 1e-7}});
    far.SetSegmentColumn("length");
    for (unsigned int seed = 1; seed <= 50; ++seed) {
        siren::dataclasses::InteractionRecord ir = SampleOne(far, seed);
        EXPECT_DOUBLE_EQ(far.GenerationProbability(nullptr, nullptr, ir), 1e7);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
