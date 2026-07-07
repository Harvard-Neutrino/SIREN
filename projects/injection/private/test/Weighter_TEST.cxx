#include <cmath>
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
// Both fixtures below use VertexWeightingMode::Fixed() on the injection and
// physical process, so InteractionProbability and NormalizedPositionProbability
// (which integrate density along the path) are never evaluated. Only the
// channel-selection and final-state-density factors are, and both reduce to
// an exact 1.0 for this single-channel DummyCrossSection with a matching
// signature, so the per-vertex generation/physical probabilities are
// predictable and the two guard branches can be pinned exactly.
// ---------------------------------------------------------------------------

namespace {

struct WeighterGuardFixture {
    std::shared_ptr<siren::injection::Injector> injector;
    std::shared_ptr<siren::injection::Weighter> weighter;
    siren::dataclasses::InteractionTree tree;
};

// `events_to_inject` seeds the realized-count normalization: Weighter reads
// EventsToInject() while InjectedEvents() is still 0 (i.e. before any
// generation). `zero_physical_normalization` adds a NormalizationConstant(0.0)
// to the physical process; PhysicalProbability multiplies its result by that
// distribution's normalization, so it becomes exactly 0.0 without perturbing
// GenerationProbability, which never applies the physical normalization.
WeighterGuardFixture BuildWeighterGuardFixture(unsigned int events_to_inject,
                                                bool zero_physical_normalization) {
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
        std::make_shared<siren::distributions::SphereVolumePositionDistribution>(
            siren::geometry::Sphere(50.0, 0.0)));

    std::shared_ptr<siren::injection::PhysicalProcess> primary_phys =
        std::make_shared<siren::injection::PhysicalProcess>(primary_type, int_col);
    primary_phys->SetWeightingMode(siren::dataclasses::VertexWeightingMode::Fixed());
    if (zero_physical_normalization) {
        primary_phys->AddPhysicalDistribution(
            std::make_shared<siren::distributions::NormalizationConstant>(0.0));
    }

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
        /*events_to_inject=*/0, /*zero_physical_normalization=*/false);
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
        /*events_to_inject=*/100, /*zero_physical_normalization=*/true);
    double weight = 0.0;
    EXPECT_NO_THROW(weight = fixture.weighter->EventWeight(fixture.tree));
    EXPECT_EQ(weight, 0.0);
}
