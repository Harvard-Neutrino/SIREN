// Serialization round-trip tests for the charm-DIS decay classes and the Weighter.
//
// These guard two previously-broken behaviors:
//   1. CharmMesonDecay / CharmMesonDecay3Body are default-constructible but used
//      cereal save() + load_and_construct(). cereal does not invoke
//      load_and_construct for a default-constructible type, so the serialized
//      body was never read on load -- the object survived only because its
//      primary_types set is a fixed const member, while every byte of the body
//      was left in the stream, corrupting whatever followed. They now use a
//      member load(); the trailing-sentinel assertions below catch any
//      regression to a body-skipping load.
//   2. Weighter::LoadWeighter was a stub that printed "not yet supported" and
//      called exit(0). It now deserializes the weighter; the round-trip reaches
//      its assertions (it would have killed the test process before).

#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <cstdio>

#include <gtest/gtest.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>

#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/CharmMesonDecay.h"
#include "SIREN/interactions/CharmMesonDecay3Body.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Weighter.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/dataclasses/Particle.h"

using namespace siren::interactions;
using namespace siren::injection;
using siren::dataclasses::ParticleType;
using siren::detector::DetectorModel;

namespace {
// Round-trip a polymorphic Decay through a binary archive followed by a sentinel
// int. The sentinel reads back correctly only if the decay body was fully
// consumed on load -- i.e. it directly detects the body-skip serialization bug.
std::shared_ptr<Decay> roundtrip_decay_with_sentinel(std::shared_ptr<Decay> orig,
                                                     bool & sentinel_ok) {
    const int kSentinel = 0x5A5A5A;
    std::stringstream ss;
    {
        cereal::BinaryOutputArchive oarchive(ss);
        oarchive(orig);
        int s = kSentinel;
        oarchive(s);
    }
    std::shared_ptr<Decay> loaded;
    int s = 0;
    {
        cereal::BinaryInputArchive iarchive(ss);
        iarchive(loaded);
        iarchive(s);
    }
    sentinel_ok = (s == kSentinel);
    return loaded;
}

std::string read_file(std::string const & path) {
    std::ifstream is(path, std::ios::binary);
    std::stringstream buf;
    buf << is.rdbuf();
    return buf.str();
}
} // namespace

// --- Polymorphic Decay round-trips (load path consumes the full body) --------

TEST(CharmSerialization, CharmMesonDecayPolymorphicRoundTrip) {
    std::shared_ptr<Decay> orig = std::make_shared<CharmMesonDecay>(ParticleType::D0);
    bool sentinel_ok = false;
    std::shared_ptr<Decay> loaded = roundtrip_decay_with_sentinel(orig, sentinel_ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<CharmMesonDecay>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(loaded->equal(*orig));
    EXPECT_EQ(orig->GetPossibleSignatures().size(),
              loaded->GetPossibleSignatures().size());
    EXPECT_TRUE(sentinel_ok) << "decay body was not fully consumed on load";
}

TEST(CharmSerialization, CharmMesonDecay3BodyPolymorphicRoundTrip) {
    std::shared_ptr<Decay> orig = std::make_shared<CharmMesonDecay3Body>(ParticleType::D0);
    bool sentinel_ok = false;
    std::shared_ptr<Decay> loaded = roundtrip_decay_with_sentinel(orig, sentinel_ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<CharmMesonDecay3Body>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(loaded->equal(*orig));
    EXPECT_EQ(orig->GetPossibleSignatures().size(),
              loaded->GetPossibleSignatures().size());
    EXPECT_TRUE(sentinel_ok) << "decay body was not fully consumed on load";
}

// --- Weighter SaveWeighter/LoadWeighter with a charm decay process -----------

TEST(CharmSerialization, WeighterWithCharmDecaySaveLoad) {
    auto detector_model = std::make_shared<DetectorModel>();
    std::shared_ptr<Decay> decay = std::make_shared<CharmMesonDecay>(ParticleType::D0);
    auto ic = std::make_shared<InteractionCollection>(
        ParticleType::D0, std::vector<std::shared_ptr<Decay>>{decay});
    auto phys = std::make_shared<PhysicalProcess>(ParticleType::D0, ic);

    // Empty injectors keep the Weighter construction data-free (Initialize()'s
    // per-injector loop is then a no-op); the charm decay still rides through
    // the serialized PhysicalProcess graph.
    std::vector<std::shared_ptr<Injector>> no_injectors;
    std::vector<std::shared_ptr<PhysicalProcess>> no_secondaries;
    Weighter w(no_injectors, detector_model, phys, no_secondaries);

    std::string base = "/tmp/siren_charm_weighter_test";
    std::string base2 = base + "_rt";
    std::remove((base + ".siren_weighter").c_str());
    std::remove((base2 + ".siren_weighter").c_str());

    ASSERT_NO_THROW(w.SaveWeighter(base));   // writes <base>.siren_weighter

    // The filename constructor calls LoadWeighter (previously a stub that
    // exit(0)'d). Reaching the next statement at all proves loading is enabled.
    std::unique_ptr<Weighter> w2;
    ASSERT_NO_THROW(w2.reset(new Weighter(no_injectors, base)));
    ASSERT_NE(w2, nullptr);

    // Re-serialize the loaded weighter; a faithful round-trip is byte-identical.
    ASSERT_NO_THROW(w2->SaveWeighter(base2));
    std::string bytes1 = read_file(base + ".siren_weighter");
    std::string bytes2 = read_file(base2 + ".siren_weighter");
    EXPECT_FALSE(bytes1.empty());
    EXPECT_EQ(bytes1, bytes2);

    std::remove((base + ".siren_weighter").c_str());
    std::remove((base2 + ".siren_weighter").c_str());
}
