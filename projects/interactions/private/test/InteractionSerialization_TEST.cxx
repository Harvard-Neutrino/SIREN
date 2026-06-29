// Serialization round-trip tests for non-default-constructible interaction
// classes relying on load_and_construct. cereal only recognizes a STATIC
// load_and_construct; a non-static one is silently ignored, leaving the body
// unread and corrupting the rest of the archive. A trailing sentinel after each
// object reads back correctly only if the body was fully consumed, so it detects
// a load that skips the body. (DarkNewsDecay excluded: abstract/python-backed.)

#include <string>
#include <vector>
#include <set>
#include <memory>
#include <sstream>

#include <gtest/gtest.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/HNLDecay.h"
#include "SIREN/interactions/ElectroweakDecay.h"
#include "SIREN/interactions/HNLDipoleDecay.h"
#include "SIREN/interactions/ElasticScattering.h"
#include "SIREN/interactions/HNLDipoleFromTable.h"
#include "SIREN/dataclasses/Particle.h"

using namespace siren::interactions;
using siren::dataclasses::ParticleType;

namespace {
// Round-trip a polymorphic base pointer through a binary archive followed by a
// sentinel int (which reads back correctly only if the body was fully consumed).
template <typename Base>
std::shared_ptr<Base> roundtrip_with_sentinel(std::shared_ptr<Base> orig, bool & sentinel_ok) {
    const int kSentinel = 0x5A5A5A;
    std::stringstream ss;
    {
        cereal::BinaryOutputArchive oarchive(ss);
        oarchive(orig);
        int s = kSentinel;
        oarchive(s);
    }
    std::shared_ptr<Base> loaded;
    int s = 0;
    {
        cereal::BinaryInputArchive iarchive(ss);
        iarchive(loaded);
        iarchive(s);
    }
    sentinel_ok = (s == kSentinel);
    return loaded;
}
} // namespace

TEST(InteractionSerialization, HNLDecayRoundTrip) {
    std::shared_ptr<Decay> orig = std::make_shared<HNLDecay>(
        0.1, std::vector<double>{0.0, 0.0, 1e-3}, HNLDecay::ChiralNature::Majorana);
    bool ok = false;
    std::shared_ptr<Decay> loaded = roundtrip_with_sentinel<Decay>(orig, ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<HNLDecay>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(ok) << "HNLDecay body not consumed on load";
}

TEST(InteractionSerialization, ElectroweakDecayRoundTrip) {
    std::shared_ptr<Decay> orig = std::make_shared<ElectroweakDecay>(
        std::set<ParticleType>{ParticleType::TauMinus});
    bool ok = false;
    std::shared_ptr<Decay> loaded = roundtrip_with_sentinel<Decay>(orig, ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<ElectroweakDecay>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(ok) << "ElectroweakDecay body not consumed on load";
}

TEST(InteractionSerialization, HNLDipoleDecayRoundTrip) {
    std::shared_ptr<Decay> orig = std::make_shared<HNLDipoleDecay>(
        0.1, std::vector<double>{0.0, 1e-6, 0.0}, HNLDipoleDecay::ChiralNature::Majorana);
    bool ok = false;
    std::shared_ptr<Decay> loaded = roundtrip_with_sentinel<Decay>(orig, ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<HNLDipoleDecay>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(ok) << "HNLDipoleDecay body not consumed on load";
}

TEST(InteractionSerialization, ElasticScatteringRoundTrip) {
    std::shared_ptr<CrossSection> orig = std::make_shared<ElasticScattering>(
        std::set<ParticleType>{ParticleType::NuMu});
    bool ok = false;
    std::shared_ptr<CrossSection> loaded = roundtrip_with_sentinel<CrossSection>(orig, ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<ElasticScattering>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(ok) << "ElasticScattering body not consumed on load";
}

TEST(InteractionSerialization, HNLDipoleFromTableRoundTrip) {
    std::shared_ptr<CrossSection> orig = std::make_shared<HNLDipoleFromTable>(
        0.1, 1e-6, HNLDipoleFromTable::HelicityChannel::Conserving);
    bool ok = false;
    std::shared_ptr<CrossSection> loaded = roundtrip_with_sentinel<CrossSection>(orig, ok);
    ASSERT_NE(loaded, nullptr);
    EXPECT_NE(std::dynamic_pointer_cast<HNLDipoleFromTable>(loaded), nullptr);
    EXPECT_TRUE(orig->equal(*loaded));
    EXPECT_TRUE(ok) << "HNLDipoleFromTable body not consumed on load";
}
