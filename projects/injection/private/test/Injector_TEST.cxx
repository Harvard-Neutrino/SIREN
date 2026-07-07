#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/injection/FailureLedger.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"

using namespace siren::detector;
using namespace siren::injection;
using namespace siren::dataclasses;
using namespace siren::interactions;
using namespace siren::utilities;

// ---------------------------------------------------------------------------
// Fail-loud hardening: Injector process-setup + reset semantics
//
// These tests are self-contained (no external data files): they build a
// minimal PrimaryInjectionProcess over a DummyCrossSection and drive only the
// process-registration and counter-reset code paths.
// ---------------------------------------------------------------------------

namespace {

// A PrimaryInjectionProcess with an interaction collection but deliberately
// NO VertexPositionDistribution.  SetPrimaryProcess must reject it.
std::shared_ptr<PrimaryInjectionProcess> MakeVertexlessPrimaryProcess() {
    ParticleType primary = ParticleType::NuMu;
    std::shared_ptr<DummyCrossSection> xs = std::make_shared<DummyCrossSection>();
    std::vector<std::shared_ptr<CrossSection>> xs_vec = {xs};
    std::shared_ptr<InteractionCollection> interactions =
        std::make_shared<InteractionCollection>(primary, xs_vec);
    std::shared_ptr<PrimaryInjectionProcess> process =
        std::make_shared<PrimaryInjectionProcess>(primary, interactions);
    // Intentionally add no distributions at all -> no vertex distribution.
    return process;
}

// A VertexPositionDistribution that always throws a tagged InjectionFailure
// instead of sampling.  Every other override is an unreachable stub: Sample()
// throws before SamplePosition/GenerationProbability/etc. can be invoked.
class ThrowingVertexDistribution : public siren::distributions::VertexPositionDistribution {
public:
    std::tuple<siren::math::Vector3D, siren::math::Vector3D> SamplePosition(
        std::shared_ptr<siren::utilities::SIREN_random> /*rand*/,
        std::shared_ptr<siren::detector::DetectorModel const> /*detector_model*/,
        std::shared_ptr<siren::interactions::InteractionCollection const> /*interactions*/,
        siren::dataclasses::PrimaryDistributionRecord & /*record*/) const override {
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::NoTargetsOnPath, "forced no targets");
    }

    void Sample(
        std::shared_ptr<siren::utilities::SIREN_random> /*rand*/,
        std::shared_ptr<siren::detector::DetectorModel const> /*detector_model*/,
        std::shared_ptr<siren::interactions::InteractionCollection const> /*interactions*/,
        siren::dataclasses::PrimaryDistributionRecord & /*record*/) const override {
        // Throw directly rather than delegating to SamplePosition() so the
        // failure fires unconditionally, with no dependence on the base
        // class's Sample() wiring.
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::NoTargetsOnPath, "forced no targets");
    }

    double GenerationProbability(
        std::shared_ptr<siren::detector::DetectorModel const> /*detector_model*/,
        std::shared_ptr<siren::interactions::InteractionCollection const> /*interactions*/,
        siren::dataclasses::InteractionRecord const & /*record*/) const override {
        return 0.0;
    }

    std::string Name() const override {
        return "ThrowingVertexDistribution";
    }

    std::shared_ptr<siren::distributions::PrimaryInjectionDistribution> clone() const override {
        return std::make_shared<ThrowingVertexDistribution>(*this);
    }

    std::tuple<siren::math::Vector3D, siren::math::Vector3D> InjectionBounds(
        std::shared_ptr<siren::detector::DetectorModel const> /*detector_model*/,
        std::shared_ptr<siren::interactions::InteractionCollection const> /*interactions*/,
        siren::dataclasses::InteractionRecord const & /*interaction*/) const override {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>();
    }

protected:
    bool equal(siren::distributions::WeightableDistribution const & distribution) const override {
        return dynamic_cast<ThrowingVertexDistribution const *>(&distribution) != nullptr;
    }

    bool less(siren::distributions::WeightableDistribution const & /*distribution*/) const override {
        return false;
    }
};

// A PrimaryInjectionProcess whose only injection distribution is a
// ThrowingVertexDistribution, so SetPrimaryProcess accepts it (it satisfies
// the vertex-distribution requirement) but GenerateEvent's first Sample()
// call always fails with a tagged FailureReason.
std::shared_ptr<PrimaryInjectionProcess> MakeAlwaysFailingPrimaryProcess() {
    ParticleType primary = ParticleType::NuMu;
    std::shared_ptr<DummyCrossSection> xs = std::make_shared<DummyCrossSection>();
    std::vector<std::shared_ptr<CrossSection>> xs_vec = {xs};
    std::shared_ptr<InteractionCollection> interactions =
        std::make_shared<InteractionCollection>(primary, xs_vec);
    std::shared_ptr<PrimaryInjectionProcess> process =
        std::make_shared<PrimaryInjectionProcess>(primary, interactions);
    process->AddPrimaryInjectionDistribution(std::make_shared<ThrowingVertexDistribution>());
    return process;
}

} // namespace

// SetPrimaryProcess with a process that has no vertex/position distribution
// must reject it with a catchable AddProcessFailure.
TEST(InjectorHardening, SetPrimaryProcessWithoutVertexDistributionThrows) {
    std::shared_ptr<DetectorModel> detector_model = std::make_shared<DetectorModel>();
    std::shared_ptr<SIREN_random> random = std::make_shared<SIREN_random>(1234);

    Injector injector(10, detector_model, random);
    std::shared_ptr<PrimaryInjectionProcess> bad_process =
        MakeVertexlessPrimaryProcess();

    EXPECT_THROW(injector.SetPrimaryProcess(bad_process),
                 siren::utilities::AddProcessFailure);

    // The thrown type is also a std::runtime_error (backward-compat contract:
    // every SIREN typed exception derives from runtime_error).
    EXPECT_THROW(injector.SetPrimaryProcess(bad_process), std::runtime_error);
}

// The constructor that takes a primary process routes through SetPrimaryProcess,
// so a vertexless process must make construction throw rather than exit(0).
TEST(InjectorHardening, ConstructWithVertexlessPrimaryProcessThrows) {
    std::shared_ptr<DetectorModel> detector_model = std::make_shared<DetectorModel>();
    std::shared_ptr<SIREN_random> random = std::make_shared<SIREN_random>(1234);
    std::shared_ptr<PrimaryInjectionProcess> bad_process =
        MakeVertexlessPrimaryProcess();

    EXPECT_THROW(
        Injector(10, detector_model, bad_process, random),
        siren::utilities::AddProcessFailure);
}

// The zero-arg ResetInjectedEvents() overload keeps the
// injection quota (EventsToInject) while clearing the run counters.  The
// counted overload changes the quota.  This pins the two overloads' distinct
// contracts without needing a full generation chain: after the counted reset
// sets a known quota, the zero-arg reset must preserve it.
TEST(InjectorHardening, ResetInjectedEventsZeroArgPreservesQuota) {
    std::shared_ptr<DetectorModel> detector_model = std::make_shared<DetectorModel>();
    std::shared_ptr<SIREN_random> random = std::make_shared<SIREN_random>(1234);
    Injector injector(0, detector_model, random);

    // Counted overload sets a new quota and clears counters.
    injector.ResetInjectedEvents(777);
    EXPECT_EQ(injector.EventsToInject(), 777u);
    EXPECT_EQ(injector.InjectionAttempts(), 0u);
    EXPECT_EQ(injector.InjectedEvents(), 0u);

    // Zero-arg overload preserves the quota and clears the run counters.
    injector.ResetInjectedEvents();
    EXPECT_EQ(injector.EventsToInject(), 777u)
        << "zero-arg ResetInjectedEvents() must not change the quota";
    EXPECT_EQ(injector.InjectionAttempts(), 0u);
    EXPECT_EQ(injector.InjectedEvents(), 0u);

    // Counted overload can still change the quota afterwards.
    injector.ResetInjectedEvents(3);
    EXPECT_EQ(injector.EventsToInject(), 3u);
    injector.ResetInjectedEvents();
    EXPECT_EQ(injector.EventsToInject(), 3u);
}

// FailureLedger aggregates repeated (depth, parent_pdg, reason) keys into one
// entry whose count increments and whose exemplar is fixed at first-record.
TEST(FailureLedgerUnit, RecordAggregatesByKeyAndKeepsFirstExemplar) {
    FailureLedger ledger;

    ledger.Record(0, 14, FailureReason::NoTargetsOnPath, "first message");
    ledger.Record(0, 14, FailureReason::NoTargetsOnPath, "second message");
    ledger.Record(1, 14, FailureReason::KinematicallyForbidden, "other reason");

    EXPECT_EQ(ledger.entries.size(), 2u);

    FailureLedger::Key repeated_key{0, 14, FailureReason::NoTargetsOnPath};
    ASSERT_EQ(ledger.entries.count(repeated_key), 1u);
    EXPECT_EQ(ledger.entries.at(repeated_key).count, 2u);
    EXPECT_EQ(ledger.entries.at(repeated_key).exemplar, "first message");

    FailureLedger::Key distinct_key{1, 14, FailureReason::KinematicallyForbidden};
    ASSERT_EQ(ledger.entries.count(distinct_key), 1u);
    EXPECT_EQ(ledger.entries.at(distinct_key).count, 1u);
    EXPECT_EQ(ledger.entries.at(distinct_key).exemplar, "other reason");

    ledger.Clear();
    EXPECT_TRUE(ledger.entries.empty());
}

// A primary vertex distribution that always throws NoTargetsOnPath must
// surface as a single depth-0 ledger entry tagged with that reason, leave a
// non-empty partial tree behind, and be cleared by ResetInjectedEvents().
TEST(InjectorHardening, PrimaryFailureRecordsLedgerAndPartialTree) {
    std::shared_ptr<DetectorModel> detector_model = std::make_shared<DetectorModel>();
    std::shared_ptr<SIREN_random> random = std::make_shared<SIREN_random>(1234);
    std::shared_ptr<PrimaryInjectionProcess> process = MakeAlwaysFailingPrimaryProcess();

    Injector injector(1, detector_model, process, random);

    siren::dataclasses::InteractionTree tree = injector.GenerateEvent();
    EXPECT_TRUE(tree.tree.empty());

    FailureLedger const & ledger = injector.GetFailureLedger();
    EXPECT_EQ(ledger.entries.size(), 1u);

    int primary_pdg = static_cast<int>(ParticleType::NuMu);
    FailureLedger::Key key{0, primary_pdg, FailureReason::NoTargetsOnPath};
    ASSERT_EQ(ledger.entries.count(key), 1u);
    EXPECT_EQ(ledger.entries.at(key).count, 1u);
    EXPECT_NE(ledger.entries.at(key).exemplar.find("forced"), std::string::npos);

    EXPECT_GE(injector.GetLastFailedTree().tree.size(), 1u);

    injector.ResetInjectedEvents();
    EXPECT_TRUE(injector.GetFailureLedger().entries.empty());
}
