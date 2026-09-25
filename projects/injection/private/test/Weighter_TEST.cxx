#include <cmath>
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/VertexWeightingMode.h"
#include "SIREN/detector/ConstantDensityDistribution.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/primary/vertex/SphereVolumePositionDistribution.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/Weighter.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

using namespace siren::injection;

// ---------------------------------------------------------------------------
// one_minus_exp_of_negative
// ---------------------------------------------------------------------------

TEST(WeighterHelpers, OneMinusExpSmallX) {
    double x = 1e-5;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-15);
}

TEST(WeighterHelpers, OneMinusExpMediumX) {
    double x = 0.5;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, OneMinusExpLargeX) {
    double x = 5.0;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, OneMinusExpVerySmallX) {
    // Exercises the Taylor expansion branch
    double x = 1e-8;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-15);
}

TEST(WeighterHelpers, OneMinusExpAtBranchPoint) {
    // Near the 0.1 branch boundary
    double x = 0.099;
    double result = one_minus_exp_of_negative(x);
    double exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);

    x = 0.101;
    result = one_minus_exp_of_negative(x);
    exact = 1.0 - std::exp(-x);
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, OneMinusExpResultBounded) {
    // Result should always be in [0, 1) for non-negative x
    for(double x = 0.0; x <= 20.0; x += 0.1) {
        double result = one_minus_exp_of_negative(x);
        EXPECT_GE(result, 0.0) << "Failed at x=" << x;
        EXPECT_LT(result, 1.0) << "Failed at x=" << x;
    }
}

// ---------------------------------------------------------------------------
// log_one_minus_exp_of_negative
// ---------------------------------------------------------------------------

TEST(WeighterHelpers, LogOneMinusExpSmallX) {
    double x = 1e-5;
    double result = log_one_minus_exp_of_negative(x);
    double exact = std::log(1.0 - std::exp(-x));
    EXPECT_NEAR(result, exact, 1e-10);
}

TEST(WeighterHelpers, LogOneMinusExpMidX) {
    double x = 1.5;
    double result = log_one_minus_exp_of_negative(x);
    double exact = std::log(1.0 - std::exp(-x));
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, LogOneMinusExpLargeX) {
    // Exercises the exp-series branch (x > 3)
    double x = 5.0;
    double result = log_one_minus_exp_of_negative(x);
    double exact = std::log(1.0 - std::exp(-x));
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(WeighterHelpers, LogOneMinusExpAtBranchPoints) {
    // Test near the 0.1 and 3.0 branch boundaries
    for(double x : {0.09, 0.11, 2.99, 3.01}) {
        double result = log_one_minus_exp_of_negative(x);
        double exact = std::log(1.0 - std::exp(-x));
        EXPECT_NEAR(result, exact, 1e-9) << "Failed at x=" << x;
    }
}

TEST(WeighterHelpers, LogOneMinusExpIsNegative) {
    // log(1 - exp(-x)) is always negative for x > 0
    for(double x = 0.001; x <= 20.0; x += 0.1) {
        double result = log_one_minus_exp_of_negative(x);
        EXPECT_LT(result, 0.0) << "Failed at x=" << x;
    }
}

// ---------------------------------------------------------------------------
// Weighter::EventWeight guard behavior
//
// The fixtures below use VertexWeightingMode::Fixed() on the injection and
// physical process, so InteractionProbability and NormalizedPositionProbability
// (which integrate density along the path) are never evaluated. Only the
// channel-selection and final-state-density factors are, and both reduce to
// an exact 1.0 for this single-channel DummyCrossSection with a matching
// signature, so the per-vertex generation/physical probabilities are
// predictable. A test position distribution supplies a controlled density.
// ---------------------------------------------------------------------------

namespace {

class GuardPositionDistribution : public siren::distributions::SphereVolumePositionDistribution {
    double density;
public:
    explicit GuardPositionDistribution(double density)
        : SphereVolumePositionDistribution(siren::geometry::Sphere(50.0, 0.0)),
          density(density) {}

    double GenerationProbability(
            std::shared_ptr<siren::detector::DetectorModel const>,
            std::shared_ptr<siren::interactions::InteractionCollection const>,
            siren::dataclasses::InteractionRecord const &) const override {
        return density;
    }
};

struct WeighterGuardFixture {
    std::shared_ptr<siren::injection::Injector> injector;
    std::shared_ptr<siren::injection::Weighter> weighter;
    siren::dataclasses::InteractionTree tree;
};

// `events_to_inject` seeds the realized-count normalization: Weighter reads
// EventsToInject() while InjectedEvents() is still 0 (i.e. before any
// generation). Physical normalization and generation density are independent,
// so invalid inputs and arithmetic overflow can be tested separately.
WeighterGuardFixture BuildWeighterGuardFixture(unsigned int events_to_inject,
                                                double physical_normalization = 1.0,
                                                double generation_density = 1.0) {
    siren::dataclasses::ParticleType primary_type = siren::dataclasses::ParticleType::NuMu;
    siren::dataclasses::ParticleType target_type = siren::dataclasses::ParticleType::Nucleon;

    std::shared_ptr<siren::detector::DetectorModel> detector_model =
        std::make_shared<siren::detector::DetectorModel>();
    detector_model->ClearSectors();
    siren::detector::DetectorSector world;
    world.name = "world";
    world.material_id = 0;  // default "VACUUM" material (real nucleon content)
    world.level = 0;
    world.geo = siren::geometry::Sphere(100.0, 0.0).create();
    world.density = siren::detector::ConstantDensityDistribution(1.0).create();
    detector_model->AddSector(world);

    std::shared_ptr<siren::interactions::DummyCrossSection> xs =
        std::make_shared<siren::interactions::DummyCrossSection>();
    std::vector<std::shared_ptr<siren::interactions::CrossSection>> xs_vec = {xs};
    std::shared_ptr<siren::interactions::InteractionCollection> int_col =
        std::make_shared<siren::interactions::InteractionCollection>(primary_type, xs_vec);

    std::shared_ptr<siren::injection::PrimaryInjectionProcess> primary_inj =
        std::make_shared<siren::injection::PrimaryInjectionProcess>(primary_type, int_col);
    primary_inj->SetWeightingMode(siren::dataclasses::VertexWeightingMode::Fixed());
    primary_inj->AddPrimaryInjectionDistribution(
        std::make_shared<GuardPositionDistribution>(generation_density));

    std::shared_ptr<siren::injection::PhysicalProcess> primary_phys =
        std::make_shared<siren::injection::PhysicalProcess>(primary_type, int_col);
    primary_phys->SetWeightingMode(siren::dataclasses::VertexWeightingMode::Fixed());
    primary_phys->AddPhysicalDistribution(
        std::make_shared<siren::distributions::NormalizationConstant>(physical_normalization));

    std::shared_ptr<siren::utilities::SIREN_random> random =
        std::make_shared<siren::utilities::SIREN_random>(1234);
    std::shared_ptr<siren::injection::Injector> injector =
        std::make_shared<siren::injection::Injector>(
            events_to_inject, detector_model, primary_inj,
            std::vector<std::shared_ptr<siren::injection::SecondaryInjectionProcess>>{},
            random);

    std::shared_ptr<siren::injection::Weighter> weighter =
        std::make_shared<siren::injection::Weighter>(
            std::vector<std::shared_ptr<siren::injection::Injector>>{injector},
            detector_model, primary_phys,
            std::vector<std::shared_ptr<siren::injection::PhysicalProcess>>{});

    siren::dataclasses::InteractionRecord record;
    record.signature.primary_type = primary_type;
    record.signature.target_type = target_type;
    record.signature.secondary_types = {primary_type, target_type};
    record.primary_mass = 0.0;
    record.primary_momentum = {2.0, 0.0, 0.0, 2.0};
    record.interaction_vertex = {0.0, 0.0, 10.0};

    siren::dataclasses::InteractionTree tree;
    tree.add_entry(record);

    return WeighterGuardFixture{injector, weighter, tree};
}

} // namespace

// generation_probability <= 0 (here: EventsToInject() == 0, the realized-count
// seed used before any generation) must raise WeightCalculationError rather
// than let the reciprocal in EventWeight blow up silently.
TEST(WeighterGuards, GenerationProbabilityNonpositiveThrowsWeightCalculationError) {
    WeighterGuardFixture fixture = BuildWeighterGuardFixture(
        /*events_to_inject=*/0);
    EXPECT_EQ(fixture.injector->InjectedEvents(), 0u);
    EXPECT_EQ(fixture.injector->EventsToInject(), 0u);
    EXPECT_THROW(fixture.weighter->EventWeight(fixture.tree),
                 siren::utilities::WeightCalculationError);
}

// physical_probability == 0 (here: a NormalizationConstant(0.0) on the
// physical process) must yield an exact 0.0 weight without raising --
// distinct from the generation-side guard above.
TEST(WeighterGuards, PhysicalProbabilityZeroGivesZeroWeightWithoutThrowing) {
    WeighterGuardFixture fixture = BuildWeighterGuardFixture(
        /*events_to_inject=*/100, /*physical_normalization=*/0.0);
    double weight = 0.0;
    EXPECT_NO_THROW(weight = fixture.weighter->EventWeight(fixture.tree));
    EXPECT_EQ(weight, 0.0);
}

TEST(WeighterGuards, InvalidPhysicalProbabilityThrows) {
    for(double probability : {-1.0, std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
        SCOPED_TRACE(probability);
        auto fixture = BuildWeighterGuardFixture(1, probability);
        EXPECT_THROW(fixture.weighter->EventWeight(fixture.tree),
                     siren::utilities::WeightCalculationError);
    }
}

TEST(WeighterGuards, InvalidGenerationProbabilityThrowsEvenWithZeroPhysicalProbability) {
    for(double probability : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
        SCOPED_TRACE(probability);
        for(double physical : {0.0, 1.0}) {
            auto fixture = BuildWeighterGuardFixture(1, physical, probability);
            EXPECT_THROW(fixture.weighter->EventWeight(fixture.tree),
                         siren::utilities::WeightCalculationError);
        }
    }
}

TEST(WeighterGuards, FinitePositiveProbabilitiesGiveExpectedWeight) {
    auto fixture = BuildWeighterGuardFixture(100, 2.0, 0.25);
    EXPECT_DOUBLE_EQ(fixture.weighter->EventWeight(fixture.tree), 2.0 / 25.0);
}

TEST(WeighterGuards, InverseWeightOverflowThrowsInsteadOfReturningZero) {
    auto fixture = BuildWeighterGuardFixture(1, 1e-308, 1e100);
    EXPECT_THROW(fixture.weighter->EventWeight(fixture.tree),
                 siren::utilities::WeightCalculationError);
}

TEST(WeighterGuards, WeightOverflowThrows) {
    // The first ratio remains representable; the second underflows to zero.
    for(double generation : {1.0, 1e-100}) {
        auto fixture = BuildWeighterGuardFixture(1, 1e308, generation * 0.01);
        EXPECT_THROW(fixture.weighter->EventWeight(fixture.tree),
                     siren::utilities::WeightCalculationError);
    }
}

TEST(WeighterGuards, PooledInverseWeightOverflowThrows) {
    auto first = BuildWeighterGuardFixture(1, 1.0, 1e308);
    auto second = BuildWeighterGuardFixture(1, 1.0, 1e308);
    EXPECT_GT(first.weighter->EventWeight(first.tree), 0.0);
    Weighter pooled({first.injector, second.injector},
                    first.weighter->GetDetectorModel(),
                    first.weighter->GetPrimaryPhysicalProcess());
    EXPECT_THROW(pooled.EventWeight(first.tree),
                 siren::utilities::WeightCalculationError);
}

TEST(WeighterGuards, ZeroPhysicalProbabilityDoesNotHideInvalidPooledInjector) {
    auto first = BuildWeighterGuardFixture(1, 0.0);
    auto second = BuildWeighterGuardFixture(0, 0.0);
    Weighter pooled({first.injector, second.injector},
                    first.weighter->GetDetectorModel(),
                    first.weighter->GetPrimaryPhysicalProcess());
    EXPECT_THROW(pooled.EventWeight(first.tree),
                 siren::utilities::WeightCalculationError);
}

TEST(WeighterGuards, EmptyInjectorPoolThrowsInsteadOfReturningInfinity) {
    auto fixture = BuildWeighterGuardFixture(1);
    Weighter empty({}, fixture.weighter->GetDetectorModel(),
                   fixture.weighter->GetPrimaryPhysicalProcess());
    EXPECT_THROW(empty.EventWeight(fixture.tree),
                 siren::utilities::WeightCalculationError);
}

TEST(WeighterGuards, ProcessWeightRejectsInvalidProbabilitiesAndOverflow) {
    for(auto const & probabilities : std::vector<std::pair<double, double>>{
            {-1.0, 1.0}, {1.0, -1.0}, {0.0, 0.0}, {1.0, 0.0},
            {std::numeric_limits<double>::infinity(), 1.0},
            {1.0, std::numeric_limits<double>::quiet_NaN()}, {1e308, 0.01}}) {
        auto fixture = BuildWeighterGuardFixture(1, probabilities.first, probabilities.second);
        PrimaryProcessWeighter process(fixture.weighter->GetPrimaryPhysicalProcess(),
                                       fixture.injector->GetPrimaryProcess(),
                                       fixture.weighter->GetDetectorModel());
        auto const & datum = *fixture.tree.tree.front();
        auto bounds = fixture.injector->PrimaryInjectionBounds(datum.record);
        EXPECT_THROW(process.EventWeight(bounds, datum),
                     siren::utilities::WeightCalculationError);
    }
}

TEST(WeighterGuards, ProcessWeightPreservesZeroAndPositiveWeights) {
    for(double physical : {0.0, 2.0}) {
        auto fixture = BuildWeighterGuardFixture(1, physical, 0.25);
        PrimaryProcessWeighter process(fixture.weighter->GetPrimaryPhysicalProcess(),
                                       fixture.injector->GetPrimaryProcess(),
                                       fixture.weighter->GetDetectorModel());
        auto const & datum = *fixture.tree.tree.front();
        auto bounds = fixture.injector->PrimaryInjectionBounds(datum.record);
        EXPECT_DOUBLE_EQ(process.EventWeight(bounds, datum), physical / 0.25);
    }
}

TEST(WeighterBreakdown, InvalidPhysicalProbabilityIsFlaggedAndTotalIsNaN) {
    for(double probability : {-1.0, std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
        auto fixture = BuildWeighterGuardFixture(1, probability);
        auto breakdown = fixture.weighter->EventWeightWithBreakdown(fixture.tree);
        EXPECT_TRUE(std::isnan(breakdown.total));
        ASSERT_EQ(breakdown.vertices.size(), 1u);
        auto const & flags = breakdown.vertices.front().flags;
        std::string expected = probability < 0.0 ? "physical density negative"
                                               : "physical density non-finite";
        EXPECT_NE(std::find(flags.begin(), flags.end(), expected), flags.end());
    }
}

TEST(WeighterBreakdown, InvalidGenerationProbabilityIsFlaggedEvenWithZeroPhysicalSupport) {
    for(double probability : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
        auto fixture = BuildWeighterGuardFixture(1, 0.0, probability);
        auto breakdown = fixture.weighter->EventWeightWithBreakdown(fixture.tree);
        EXPECT_TRUE(std::isnan(breakdown.total));
        ASSERT_EQ(breakdown.vertices.size(), 1u);
        auto const & flags = breakdown.vertices.front().flags;
        std::string expected = probability == 0.0 ? "generation density zero"
                             : probability < 0.0 ? "generation density negative"
                                                 : "generation density non-finite";
        EXPECT_NE(std::find(flags.begin(), flags.end(), expected), flags.end());
    }
}

TEST(WeighterBreakdown, OverflowProducesNaNWithDiagnosticFlag) {
    for(auto const & probabilities : std::vector<std::pair<double, double>>{
            {1e-308, 1e100}, {1e308, 0.01}, {1e308, 1e-102}}) {
        auto fixture = BuildWeighterGuardFixture(1, probabilities.first, probabilities.second);
        auto breakdown = fixture.weighter->EventWeightWithBreakdown(fixture.tree);
        EXPECT_TRUE(std::isnan(breakdown.total));
        ASSERT_EQ(breakdown.vertices.size(), 1u);
        EXPECT_FALSE(breakdown.vertices.front().flags.empty());
    }
}

TEST(WeighterBreakdown, PooledInverseWeightOverflowProducesNaN) {
    auto first = BuildWeighterGuardFixture(1, 1.0, 1e308);
    auto second = BuildWeighterGuardFixture(1, 1.0, 1e308);
    Weighter pooled({first.injector, second.injector},
                    first.weighter->GetDetectorModel(),
                    first.weighter->GetPrimaryPhysicalProcess());
    auto breakdown = pooled.EventWeightWithBreakdown(first.tree);
    EXPECT_TRUE(std::isnan(breakdown.total));
    ASSERT_EQ(breakdown.vertices.size(), 2u);
    EXPECT_FALSE(breakdown.vertices.back().flags.empty());
}

TEST(WeighterBreakdown, ZeroSupportDoesNotHideInvalidPooledInjector) {
    auto first = BuildWeighterGuardFixture(1, 0.0);
    auto second = BuildWeighterGuardFixture(0, 0.0);
    Weighter pooled({first.injector, second.injector},
                    first.weighter->GetDetectorModel(),
                    first.weighter->GetPrimaryPhysicalProcess());
    EXPECT_TRUE(std::isnan(pooled.EventWeightWithBreakdown(first.tree).total));
}

TEST(WeighterBreakdown, ValidZeroAndPositiveTotalsMatchScalarWeight) {
    for(double physical : {0.0, 2.0}) {
        auto fixture = BuildWeighterGuardFixture(100, physical, 0.25);
        auto breakdown = fixture.weighter->EventWeightWithBreakdown(fixture.tree);
        EXPECT_DOUBLE_EQ(breakdown.total, fixture.weighter->EventWeight(fixture.tree));
        ASSERT_EQ(breakdown.vertices.size(), 1u);
        if(physical == 0.0) {
            EXPECT_DOUBLE_EQ(breakdown.total, 0.0);
            EXPECT_FALSE(breakdown.vertices.front().flags.empty());
        } else {
            EXPECT_TRUE(breakdown.vertices.front().flags.empty());
        }
    }
}

TEST(WeighterBreakdown, EmptyPoolHasNaNTotal) {
    auto fixture = BuildWeighterGuardFixture(1);
    Weighter empty({}, fixture.weighter->GetDetectorModel(),
                   fixture.weighter->GetPrimaryPhysicalProcess());
    EXPECT_TRUE(std::isnan(empty.EventWeightWithBreakdown(fixture.tree).total));
}
