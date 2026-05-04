#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <gtest/gtest.h>

#include "SIREN/utilities/Integration.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/distributions/primary/energy/TabulatedFluxDistribution.h"
#include "SIREN/distributions/primary/energy_direction/Tabulated2DFluxDistribution.h"

using namespace siren::distributions;
using namespace siren::dataclasses;

// ---------------------------------------------------------------------------
// 1. Integration utilities
// ---------------------------------------------------------------------------

TEST(TrapezoidIntegrate, LinearFunction) {
    // f(x)=2x on [0,1]: exact integral = 1.0
    // Trapezoid rule is exact for linear functions.
    std::vector<double> x = {0.0, 0.25, 0.5, 0.75, 1.0};
    std::vector<double> y = {0.0, 0.5, 1.0, 1.5, 2.0};
    double result = siren::utilities::trapezoidIntegrate(x, y);
    EXPECT_DOUBLE_EQ(result, 1.0);
}

TEST(TrapezoidIntegrate, ConstantFunction) {
    // f(x)=3 on [0,5]: integral = 15.0
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y(6, 3.0);
    double result = siren::utilities::trapezoidIntegrate(x, y);
    EXPECT_DOUBLE_EQ(result, 15.0);
}

TEST(TrapezoidIntegrate, SingleInterval) {
    // Two points: f(0)=0, f(1)=4 => area = 0.5*(0+4)*1 = 2
    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {0.0, 4.0};
    double result = siren::utilities::trapezoidIntegrate(x, y);
    EXPECT_DOUBLE_EQ(result, 2.0);
}

TEST(TrapezoidIntegrate, MismatchedSizesThrow) {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0};
    EXPECT_THROW(siren::utilities::trapezoidIntegrate(x, y), std::runtime_error);
}

// ---------------------------------------------------------------------------
// 1b. Bounded trapezoidIntegrate
// ---------------------------------------------------------------------------

TEST(TrapezoidIntegrateBounded, FullRangeMatchesUnbounded) {
    // Bounds equal to the full node range should give the same result
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 1.0, 4.0, 9.0, 16.0};
    double full = siren::utilities::trapezoidIntegrate(x, y);
    double bounded = siren::utilities::trapezoidIntegrate(x, y, 0.0, 4.0);
    EXPECT_DOUBLE_EQ(full, bounded);
}

TEST(TrapezoidIntegrateBounded, SubrangeOnNodes) {
    // Bounds that land exactly on nodes: integrate [1,3] of f(x)=3
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y(6, 3.0);
    double result = siren::utilities::trapezoidIntegrate(x, y, 1.0, 3.0);
    EXPECT_DOUBLE_EQ(result, 6.0);
}

TEST(TrapezoidIntegrateBounded, SubrangeBetweenNodes) {
    // f(x)=2x on nodes [0, 1, 2, 3, 4].  Integrate over [0.5, 3.5].
    // Exact integral of 2x from 0.5 to 3.5 = [x^2] = 12.25 - 0.25 = 12.0
    // Trapezoid is exact for linear functions, so this should be exact.
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 2.0, 4.0, 6.0, 8.0};
    double result = siren::utilities::trapezoidIntegrate(x, y, 0.5, 3.5);
    EXPECT_DOUBLE_EQ(result, 12.0);
}

TEST(TrapezoidIntegrateBounded, BoundsWithinSingleInterval) {
    // Both xmin and xmax fall inside [1, 2].  f(x)=2x.
    // Integral of 2x from 1.25 to 1.75 = [x^2] = 3.0625 - 1.5625 = 1.5
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 2.0, 4.0, 6.0};
    double result = siren::utilities::trapezoidIntegrate(x, y, 1.25, 1.75);
    EXPECT_DOUBLE_EQ(result, 1.5);
}

TEST(TrapezoidIntegrateBounded, BoundsOutsideNodesClamp) {
    // Bounds wider than the data: should clamp to the node range.
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {2.0, 2.0, 2.0};  // constant
    double result = siren::utilities::trapezoidIntegrate(x, y, -10.0, 100.0);
    EXPECT_DOUBLE_EQ(result, 4.0);  // 2 * (3-1) = 4
}

TEST(TrapezoidIntegrateBounded, ZeroWidthRange) {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 2.0};
    EXPECT_DOUBLE_EQ(siren::utilities::trapezoidIntegrate(x, y, 1.5, 1.5), 0.0);
    // Also on a node
    EXPECT_DOUBLE_EQ(siren::utilities::trapezoidIntegrate(x, y, 1.0, 1.0), 0.0);
}

TEST(TrapezoidIntegrateBounded, NarrowingBoundsReducesIntegral) {
    // Motivating use case: a flux table spanning [1,10], integrated first
    // over the full range, then over a narrower sub-range.
    std::vector<double> x = {1.0, 2.0, 4.0, 7.0, 10.0};
    std::vector<double> y = {0.5, 1.0, 3.0, 2.0, 0.5};
    double full    = siren::utilities::trapezoidIntegrate(x, y);
    double bounded = siren::utilities::trapezoidIntegrate(x, y, 2.0, 7.0);
    EXPECT_GT(full, bounded);
    EXPECT_GT(bounded, 0.0);
}

TEST(TrapezoidIntegrateBounded, NonUniformSpacing) {
    // Nodes at 0, 1, 5, 6 with y = 0, 2, 10, 12 (f(x)=2x).
    // Integrate [0.5, 5.5].  Exact integral of 2x = [x^2] = 30.25 - 0.25 = 30.
    // Trapezoid is exact for linear, so the bounded version should agree.
    std::vector<double> x = {0.0, 1.0, 5.0, 6.0};
    std::vector<double> y = {0.0, 2.0, 10.0, 12.0};
    double result = siren::utilities::trapezoidIntegrate(x, y, 0.5, 5.5);
    EXPECT_DOUBLE_EQ(result, 30.0);
}

TEST(TrapezoidIntegrateBounded, BoundsCompletelyOutside) {
    // Range entirely below the nodes
    std::vector<double> x = {5.0, 6.0, 7.0};
    std::vector<double> y = {1.0, 1.0, 1.0};
    EXPECT_DOUBLE_EQ(siren::utilities::trapezoidIntegrate(x, y, 0.0, 4.0), 0.0);
    // Range entirely above the nodes
    EXPECT_DOUBLE_EQ(siren::utilities::trapezoidIntegrate(x, y, 8.0, 10.0), 0.0);
}

TEST(TrapezoidIntegrateBounded, BoundsOnNodeEdgesExact) {
    // xmin and xmax land exactly on node boundaries
    // f(x) = x on [0,1,2,3,4,5], integrate [2,4]
    // Exact integral of x from 2 to 4 = [x^2/2] = 8 - 2 = 6
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double result = siren::utilities::trapezoidIntegrate(x, y, 2.0, 4.0);
    EXPECT_DOUBLE_EQ(result, 6.0);
}

TEST(SimpsonIntegrate2D, ConstantFunction) {
    // f(x,y)=1 over [0,1]x[0,1] => integral = 1.0
    auto f = [](double, double) -> double { return 1.0; };
    double result = siren::utilities::simpsonIntegrate2D(f, 0.0, 1.0, 0.0, 1.0);
    EXPECT_NEAR(result, 1.0, 1e-6);
}

TEST(SimpsonIntegrate2D, SeparableFunction) {
    // f(x,y)=x*y over [0,1]x[0,1] => integral = 0.25
    auto f = [](double x, double y) -> double { return x * y; };
    double result = siren::utilities::simpsonIntegrate2D(f, 0.0, 1.0, 0.0, 1.0);
    EXPECT_NEAR(result, 0.25, 1e-6);
}

TEST(SimpsonIntegrate2D, OddFunctionZeroIntegral) {
    // f(x,y)=x over [-1,1]x[0,1] => integral = 0
    // Tests the divide-by-zero fix for near-zero estimates.
    auto f = [](double x, double) -> double { return x; };
    double result = siren::utilities::simpsonIntegrate2D(f, -1.0, 1.0, 0.0, 1.0);
    EXPECT_NEAR(result, 0.0, 1e-10);
}

TEST(SimpsonIntegrate2D, NegativeToleranceThrows) {
    auto f = [](double, double) -> double { return 1.0; };
    EXPECT_THROW(
        siren::utilities::simpsonIntegrate2D(f, 0.0, 1.0, 0.0, 1.0, -1e-3),
        std::runtime_error
    );
}

// ---------------------------------------------------------------------------
// 2. TabulatedFluxDistribution romberg flag
// ---------------------------------------------------------------------------

TEST(TabulatedFluxDistribution, RombergFlag) {
    // Simple triangular flux: f(E) rises linearly.
    // Use enough points so Romberg and trapezoid give different results.
    std::vector<double> energies = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> flux     = {1.0, 4.0, 9.0, 16.0, 25.0}; // E^2

    TabulatedFluxDistribution distR(energies, flux, false, true);   // romberg
    TabulatedFluxDistribution distT(energies, flux, false, false);  // trapezoid

    double integralR = distR.GetIntegral();
    double integralT = distT.GetIntegral();

    EXPECT_GT(integralR, 0.0);
    EXPECT_GT(integralT, 0.0);
    // Romberg interpolates between nodes; trapezoid uses nodes directly.
    // For a non-linear function they should differ.
    EXPECT_NE(integralR, integralT);
}

// ---------------------------------------------------------------------------
// 3. Tabulated2DFluxDistribution
// ---------------------------------------------------------------------------

// Helper: build a flat 3x3 grid with flux=1 everywhere.
// energies [1,2,3], cosZeniths [-1,0,1], flux all 1.0
// The grid is stored in row-major order: for each energy, iterate cosZeniths.
static void MakeFlatGrid(std::vector<double>& energies,
                         std::vector<double>& cosZeniths,
                         std::vector<double>& flux) {
    double e_vals[] = {1.0, 2.0, 3.0};
    double c_vals[] = {-1.0, 0.0, 1.0};
    energies.clear();
    cosZeniths.clear();
    flux.clear();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            energies.push_back(e_vals[i]);
            cosZeniths.push_back(c_vals[j]);
            flux.push_back(1.0);
        }
    }
}

// --- Construction ---

TEST(Tabulated2DFlux, ConstructFromVectors) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    ASSERT_NO_THROW(Tabulated2DFluxDistribution dist(e, cz, f));
}

TEST(Tabulated2DFlux, ConstructFromVectorsWithBounds) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    ASSERT_NO_THROW(
        Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f)
    );
}

TEST(Tabulated2DFlux, ConstructEmptyVectorsThrows) {
    std::vector<double> e, cz, f;
    EXPECT_THROW(
        Tabulated2DFluxDistribution dist(e, cz, f),
        std::runtime_error
    );
}

TEST(Tabulated2DFlux, ConstructFromFile) {
    // Write a unique temporary flux table file, then load it.
    char tmpfile[] = "/tmp/siren_test_flux2d_XXXXXX";
    int fd = mkstemp(tmpfile);
    ASSERT_NE(fd, -1);
    {
        std::ofstream out(tmpfile);
        double e_vals[] = {1.0, 2.0, 3.0};
        double c_vals[] = {-1.0, 0.0, 1.0};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                out << e_vals[i] << " " << c_vals[j] << " " << 1.0 << "\n";
            }
        }
    }
    close(fd);
    ASSERT_NO_THROW(Tabulated2DFluxDistribution dist(tmpfile));
    std::remove(tmpfile);
}

// --- PDF evaluation ---

TEST(Tabulated2DFlux, IntegralConstantFlux) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(e, cz, f);
    double integral = dist.GetIntegral();
    // For constant flux=1, integral = (eMax-eMin)*(czMax-czMin) = 2*2 = 4
    EXPECT_NEAR(integral, 4.0, 0.1);
}

TEST(Tabulated2DFlux, SamplePDFConstant) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(e, cz, f);
    double integral = dist.GetIntegral();
    // Normalized PDF for constant flux = 1/integral
    double pdf_val = dist.SamplePDF(2.0, 0.0);
    EXPECT_NEAR(pdf_val, 1.0 / integral, 1e-6);
}

TEST(Tabulated2DFlux, SampleUnnormedPDF) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(e, cz, f);
    double val = dist.SampleUnnormedPDF(2.0, 0.0);
    EXPECT_NEAR(val, 1.0, 1e-6);
}

TEST(Tabulated2DFlux, GetNodes) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(e, cz, f);

    std::vector<double> eNodes = dist.GetEnergyNodes();
    std::vector<double> czNodes = dist.GetCosZenithNodes();
    ASSERT_EQ(eNodes.size(), 9u);
    ASSERT_EQ(czNodes.size(), 9u);
    // Check first and last
    EXPECT_DOUBLE_EQ(eNodes.front(), 1.0);
    EXPECT_DOUBLE_EQ(eNodes.back(), 3.0);
    EXPECT_DOUBLE_EQ(czNodes.front(), -1.0);
    EXPECT_DOUBLE_EQ(czNodes.back(), 1.0);
}

// --- Sampling (SampleEnergyAndDirection) ---

TEST(Tabulated2DFlux, SamplingBounds) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    PrimaryDistributionRecord record(ParticleType::NuMu);
    record.SetMass(0.0);
    record.SetEnergy(1.0);
    record.SetDirection({0.0, 0.0, 1.0});

    const size_t N = 10000;
    for (size_t i = 0; i < N; ++i) {
        PrimaryDistributionRecord rec(ParticleType::NuMu);
        rec.SetMass(0.0);
        rec.SetEnergy(1.0);
        rec.SetDirection({0.0, 0.0, 1.0});

        auto result = dist.SampleEnergyAndDirection(rand, nullptr, nullptr, rec);
        double energy = result.first;
        siren::math::Vector3D dir = result.second;

        // Energy within bounds
        EXPECT_GE(energy, 1.0);
        EXPECT_LE(energy, 3.0);

        // Direction is a unit vector
        double mag = dir.magnitude();
        EXPECT_NEAR(mag, 1.0, 1e-10);

        // cos(zenith) = dir.z should be in [-1, 1]
        double cosZ = dir.GetZ();
        EXPECT_GE(cosZ, -1.0 - 1e-10);
        EXPECT_LE(cosZ, 1.0 + 1e-10);
    }
}

// --- GenerationProbability ---

TEST(Tabulated2DFlux, GenerationProbabilityInBounds) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    InteractionRecord record;
    record.primary_momentum[0] = 2.0;  // energy
    // Direction along z: px=0, py=0, pz=2
    record.primary_momentum[1] = 0.0;
    record.primary_momentum[2] = 0.0;
    record.primary_momentum[3] = 2.0;

    double prob = dist.GenerationProbability(nullptr, nullptr, record);
    EXPECT_GT(prob, 0.0);
}

TEST(Tabulated2DFlux, GenerationProbabilityBelowMin) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    InteractionRecord record;
    record.primary_momentum[0] = 0.5;  // below energyMin=1
    record.primary_momentum[1] = 0.0;
    record.primary_momentum[2] = 0.0;
    record.primary_momentum[3] = 0.5;

    double prob = dist.GenerationProbability(nullptr, nullptr, record);
    EXPECT_EQ(prob, 0.0);
}

TEST(Tabulated2DFlux, GenerationProbabilityAboveMax) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    InteractionRecord record;
    record.primary_momentum[0] = 5.0;  // above energyMax=3
    record.primary_momentum[1] = 0.0;
    record.primary_momentum[2] = 0.0;
    record.primary_momentum[3] = 5.0;

    double prob = dist.GenerationProbability(nullptr, nullptr, record);
    EXPECT_EQ(prob, 0.0);
}

// --- Equality and ordering ---

TEST(Tabulated2DFlux, EqualitySameData) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution a(1.0, 3.0, -1.0, 1.0, e, cz, f);
    Tabulated2DFluxDistribution b(1.0, 3.0, -1.0, 1.0, e, cz, f);
    EXPECT_TRUE(a == b);
}

TEST(Tabulated2DFlux, InequalityDifferentBounds) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution a(1.0, 3.0, -1.0, 1.0, e, cz, f);
    Tabulated2DFluxDistribution b(1.5, 3.0, -1.0, 1.0, e, cz, f);
    EXPECT_FALSE(a == b);
}

TEST(Tabulated2DFlux, LessWithNullDynamicCast) {
    // less() with a non-Tabulated2DFluxDistribution should return false.
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution a(1.0, 3.0, -1.0, 1.0, e, cz, f);

    // Use a NormalizationConstant as a different WeightableDistribution subclass.
    NormalizationConstant other(1.0);
    EXPECT_FALSE(a < other);
}

TEST(Tabulated2DFlux, Irreflexivity) {
    // Strict weak ordering: !(a < a)
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution a(1.0, 3.0, -1.0, 1.0, e, cz, f);
    EXPECT_FALSE(a < a);
}

// --- Name ---

TEST(Tabulated2DFlux, Name) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(e, cz, f);
    EXPECT_EQ(dist.Name(), "Tabulated2DFluxDistribution");
}

// --- Clone ---

TEST(Tabulated2DFlux, Clone) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    auto cloned = dist.clone();
    ASSERT_NE(cloned, nullptr);

    // The clone should be equal to the original via the WeightableDistribution
    // operator==.
    const WeightableDistribution& orig_ref = dist;
    const WeightableDistribution& clone_ref = *cloned;
    EXPECT_TRUE(orig_ref == clone_ref);
}

// --- Validation and edge cases ---

TEST(Tabulated2DFlux, SetCosZenithBoundsInvalidThrows) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(e, cz, f);
    // Outside [-1, 1]
    EXPECT_THROW(dist.SetCosZenithBounds(-2.0, 1.0), std::runtime_error);
    EXPECT_THROW(dist.SetCosZenithBounds(-1.0, 1.5), std::runtime_error);
    // min > max
    EXPECT_THROW(dist.SetCosZenithBounds(0.5, -0.5), std::runtime_error);
    // Valid should not throw
    EXPECT_NO_THROW(dist.SetCosZenithBounds(-0.5, 0.5));
}

TEST(Tabulated2DFlux, SamplingAfterBoundsChange) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    // Sample once to initialize MH state
    auto rand = std::make_shared<siren::utilities::SIREN_random>();
    PrimaryDistributionRecord rec(ParticleType::NuMu);
    rec.SetMass(0.0);
    rec.SetEnergy(1.0);
    rec.SetDirection({0.0, 0.0, 1.0});
    dist.SampleEnergyAndDirection(rand, nullptr, nullptr, rec);

    // Narrow the energy bounds and sample again -- results should be valid
    dist.SetEnergyBounds(1.5, 2.5);
    for (int i = 0; i < 1000; ++i) {
        PrimaryDistributionRecord r(ParticleType::NuMu);
        r.SetMass(0.0);
        r.SetEnergy(1.0);
        r.SetDirection({0.0, 0.0, 1.0});
        auto result = dist.SampleEnergyAndDirection(rand, nullptr, nullptr, r);
        EXPECT_GE(result.first, 1.5);
        EXPECT_LE(result.first, 2.5);
    }
}

TEST(Tabulated2DFlux, GenerationProbabilityZeroMomentum) {
    std::vector<double> e, cz, f;
    MakeFlatGrid(e, cz, f);
    Tabulated2DFluxDistribution dist(1.0, 3.0, -1.0, 1.0, e, cz, f);

    InteractionRecord record;
    record.primary_momentum[0] = 2.0;
    record.primary_momentum[1] = 0.0;
    record.primary_momentum[2] = 0.0;
    record.primary_momentum[3] = 0.0;  // zero 3-momentum

    double prob = dist.GenerationProbability(nullptr, nullptr, record);
    EXPECT_EQ(prob, 0.0);
}

TEST(Tabulated2DFlux, UnsortedInputBoundsCorrect) {
    // Input in non-sorted order: energy goes 3, 1, 2
    std::vector<double> e   = {3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    std::vector<double> cz  = {1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 0.0};
    std::vector<double> flux(9, 1.0);
    Tabulated2DFluxDistribution dist(e, cz, flux);

    // Bounds should be derived from min/max, not first/last
    double integral = dist.GetIntegral();
    // Integral of constant 1 over [1,3] x [-1,1] = 4
    EXPECT_NEAR(integral, 4.0, 0.5);
}

TEST(Tabulated2DFlux, MalformedFileThrows) {
    char tmpfile[] = "/tmp/siren_test_malformed_XXXXXX";
    int fd = mkstemp(tmpfile);
    ASSERT_NE(fd, -1);
    {
        std::ofstream out(tmpfile);
        out << "1.0 2.0 3.0\n";
        out << "bad_data\n";  // only one token, not three
        out << "4.0 5.0 6.0\n";
    }
    close(fd);
    EXPECT_THROW(Tabulated2DFluxDistribution dist(tmpfile), std::runtime_error);
    std::remove(tmpfile);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
