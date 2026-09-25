#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"

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
